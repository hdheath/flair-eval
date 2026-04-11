#!/usr/bin/env python3
"""
ted_component_diagnostic.py — Diagnose TED component score interactions with TP/FP status.

For each TED configuration and each end type (5′ TSS / 3′ TTS), produces:

Panel 1 — ROC curves: discriminative power of each component score (depth, model,
          annot, reality, depth+annot renormalized) for separating TP from FP isoforms.
Panel 2 — Component score violins: distribution of each component score split by TP/FP.
Panel 3 — Joint score heatmaps: 2-D log-ratio density (TP vs FP) for every score-pair.
Panel 4 — Signal vs score colored by TP/FP: scatter with boundary signal on x-axis.
Panel 5 — Weight sweep ternary surface: best F1 over (w_depth, w_model, w_annot) grid.
Panel 6 — Marginal model value: delta-AUC bar chart per end × config.
Summary — summary_stats.tsv with per-component AUCs, optimal weights, counts.

TP/FP classification: peak within ±50 bp (CAGE for 5′, QuantSeq for 3′).

Inputs match the SjcAltEndAnalysis pattern:
  --bed          label:path  (BED12+TED isoform files)
  --cage-peaks   BED6 CAGE peaks
  --qs-peaks     BED6 QuantSeq peaks
  --cage-plus/--cage-minus/--qs-plus/--qs-minus  bedGraph signal tracks
  --output       output directory
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from bisect import bisect_left
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
import numpy as np
from scipy.stats import spearmanr

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, PALETTE
    from signal_utils import load_signal_tracks, isoform_signal, SIG_WINDOW
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, PALETTE
    from evaluation.signal_utils import load_signal_tracks, isoform_signal, SIG_WINDOW

apply_rc()
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s")
log = logging.getLogger(__name__)

# ── Constants ───────────────────────────────────────────────────────────────

WINDOW = 50  # peak match distance

TED_SCORE_KEYS = [
    "TED_confidence", "TED_tss_reality", "TED_tts_reality",
    "TED_depth", "TED_tss_model", "TED_tts_model",
    "TED_tss_annot", "TED_tts_annot", "TED_tss_annot_dist", "TED_tts_annot_dist",
]

COMPONENTS = ["depth", "model", "annot", "reality"]
COMP_COLORS = {
    "depth":   PALETTE[4],   # blue
    "model":   PALETTE[5],   # vermillion
    "annot":   PALETTE[2],   # green
    "reality": PALETTE[0],   # orange
    "depth+annot": PALETTE[6],  # purple
}

COMP_PAIRS = [("depth", "model"), ("depth", "annot"), ("model", "annot")]


# ── Parsing ─────────────────────────────────────────────────────────────────

def parse_bed_with_ted_scores(path: str) -> List[dict]:
    """Parse BED12+ file with TED extra columns."""
    isoforms: list[dict] = []
    has_ted = None
    with open(path) as f:
        for line in f:
            if line.startswith(("#", "track")):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 12:
                continue
            if has_ted is None:
                has_ted = len(c) >= 20
                if not has_ted:
                    return []
            chrom, start, end = c[0], int(c[1]), int(c[2])
            name, strand = c[3], c[5]
            try:
                score = int(c[4])
            except ValueError:
                score = 0
            bc = int(c[9])
            bsz = [int(x) for x in c[10].rstrip(",").split(",") if x]
            bst = [int(x) for x in c[11].rstrip(",").split(",") if x]
            juncs = []
            for i in range(bc - 1):
                juncs.append((start + bst[i] + bsz[i], start + bst[i + 1]))
            tss = start if strand == "+" else end
            tts = end if strand == "+" else start
            iso = dict(
                chrom=chrom, start=start, end=end, name=name,
                score=score, strand=strand, tss=tss, tts=tts,
                junctions=tuple(juncs), n_exons=bc,
            )
            for j, key in enumerate(TED_SCORE_KEYS):
                idx = 12 + j
                if idx < len(c):
                    try:
                        iso[key] = float(c[idx])
                    except ValueError:
                        iso[key] = 0.0
                else:
                    iso[key] = 0.0
            isoforms.append(iso)
    return isoforms


def _parse_peaks_bed(path: str) -> Dict[Tuple[str, str], List[int]]:
    """Parse BED6 peaks → {(chrom, strand): sorted midpoints}."""
    peaks: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            chrom, start, end, strand = cols[0], int(cols[1]), int(cols[2]), cols[5]
            peaks[(chrom, strand)].append((start + end) // 2)
    for k in peaks:
        peaks[k].sort()
    return peaks


def _nearest_peak(pos: int, sorted_peaks: List[int]) -> Optional[int]:
    """Return nearest peak midpoint within WINDOW, or None."""
    if not sorted_peaks:
        return None
    idx = bisect_left(sorted_peaks, pos)
    best, best_d = None, WINDOW + 1
    for i in (idx - 1, idx):
        if 0 <= i < len(sorted_peaks):
            d = abs(pos - sorted_peaks[i])
            if d < best_d:
                best_d, best = d, sorted_peaks[i]
    return best if best_d <= WINDOW else None


# ── Per-config data extraction ──────────────────────────────────────────────

def _extract(
    isoforms: List[dict],
    peaks: Dict[Tuple[str, str], List[int]],
    cage_p, cage_m, qs_p, qs_m,
    end_type: str,
) -> dict:
    """Extract arrays of component scores, signal, and TP/FP labels."""
    prefix = "tss" if end_type == "tss" else "tts"
    depths, models, annots, realities, signals = [], [], [], [], []
    tp_labels = []  # 1=TP, 0=FP

    for iso in isoforms:
        pos = iso["tss"] if end_type == "tss" else iso["tts"]
        peak_list = peaks.get((iso["chrom"], iso["strand"]), [])
        is_tp = _nearest_peak(pos, peak_list) is not None

        sig_tss, sig_tts = isoform_signal(iso, cage_p, cage_m, qs_p, qs_m)
        sig = sig_tss if end_type == "tss" else sig_tts

        depths.append(iso.get("TED_depth", 0.0))
        models.append(iso.get(f"TED_{prefix}_model", 0.0))
        annots.append(iso.get(f"TED_{prefix}_annot", 0.0))
        realities.append(iso.get(f"TED_{prefix}_reality", 0.0))
        signals.append(sig)
        tp_labels.append(1 if is_tp else 0)

    return {
        "depth": np.asarray(depths),
        "model": np.asarray(models),
        "annot": np.asarray(annots),
        "reality": np.asarray(realities),
        "signal": np.asarray(signals),
        "tp": np.asarray(tp_labels, dtype=int),
    }


# ── ROC helpers ─────────────────────────────────────────────────────────────

def _roc_curve(labels: np.ndarray, scores: np.ndarray):
    """Return (fpr, tpr, auc) for TP=1 labels."""
    if len(np.unique(labels)) < 2 or np.std(scores) == 0:
        return np.array([0, 1]), np.array([0, 1]), 0.5
    order = np.argsort(-scores)
    labels_s = labels[order]
    n_pos = labels.sum()
    n_neg = len(labels) - n_pos
    tp_cum = np.cumsum(labels_s)
    fp_cum = np.cumsum(1 - labels_s)
    tpr = np.concatenate([[0], tp_cum / n_pos])
    fpr = np.concatenate([[0], fp_cum / n_neg])
    auc = np.trapz(tpr, fpr)
    return fpr, tpr, auc


# ── Panel 1: ROC curves ────────────────────────────────────────────────────

def plot_roc(
    data_by_label: dict[str, dict],
    end_label: str,
    outdir: Path,
) -> dict:
    """ROC curves per config.  Returns {label: {comp: auc}}."""
    labels = list(data_by_label.keys())
    n = len(labels)
    if n == 0:
        return {}

    fig, axes = plt.subplots(1, n, figsize=(min(W2, n * 3.0), 3.0),
                             squeeze=False, sharey=True)
    axes = axes[0]

    auc_results = {}

    for i, label in enumerate(labels):
        ax = axes[i]
        d = data_by_label[label]
        tp = d["tp"]
        n_tp = tp.sum()
        n_fp = len(tp) - n_tp
        auc_results[label] = {}

        # Component ROC curves
        score_sets = [
            ("depth",      d["depth"]),
            ("model",      d["model"]),
            ("annot",      d["annot"]),
            ("reality",    d["reality"]),
        ]
        # depth+annot renormalized
        da_raw = d["depth"] + d["annot"]
        da_max = da_raw.max() if da_raw.max() > 0 else 1.0
        score_sets.append(("depth+annot", da_raw / da_max))

        for comp, scores in score_sets:
            fpr, tpr, auc = _roc_curve(tp, scores)
            ax.plot(fpr, tpr, label=f"{comp} ({auc:.3f})",
                    color=COMP_COLORS[comp], linewidth=1.0)
            auc_results[label][comp] = auc

        ax.plot([0, 1], [0, 1], "--", color="grey", linewidth=0.5)
        ax.legend(fontsize=5, loc="lower right", frameon=False)
        ax.set_title(f"{label}\n(TP={n_tp}, FP={n_fp})", fontsize=6)

        if i == 0:
            style_ax(ax, xlabel="FPR", ylabel="TPR")
        else:
            style_ax(ax, xlabel="FPR")

    fig.suptitle(f"Component ROC — {end_label}", fontsize=8, y=1.02)
    fig.tight_layout()
    savefig(fig, outdir / f"roc_curves_{end_label}")
    return auc_results


# ── Panel 2: TP/FP violins ─────────────────────────────────────────────────

def plot_violins(
    data_by_label: dict[str, dict],
    end_label: str,
    outdir: Path,
) -> None:
    labels = list(data_by_label.keys())
    n = len(labels)
    if n == 0:
        return

    comps = COMPONENTS
    n_comp = len(comps)

    fig, axes = plt.subplots(n, n_comp, figsize=(W2, max(W1, n * 2.0)),
                             squeeze=False, sharex="col")

    for row, label in enumerate(labels):
        d = data_by_label[label]
        tp = d["tp"].astype(bool)

        for col, comp in enumerate(comps):
            ax = axes[row, col]
            vals_tp = d[comp][tp]
            vals_fp = d[comp][~tp]

            data = []
            positions = []
            colors = []
            tick_labels = []
            if len(vals_tp) > 0:
                data.append(vals_tp)
                positions.append(0)
                colors.append(PALETTE[2])  # green
                tick_labels.append(f"TP\n({len(vals_tp)})")
            if len(vals_fp) > 0:
                data.append(vals_fp)
                positions.append(1)
                colors.append(PALETTE[5])  # vermillion
                tick_labels.append(f"FP\n({len(vals_fp)})")

            if not data:
                ax.set_visible(False)
                continue

            parts = ax.violinplot(data, positions=positions, showmedians=True,
                                  showextrema=False)
            for j, pc in enumerate(parts["bodies"]):
                pc.set_facecolor(colors[j])
                pc.set_alpha(0.5)
            parts["cmedians"].set_color("black")

            ax.set_xticks(positions)
            ax.set_xticklabels(tick_labels, fontsize=5)

            if row == 0:
                ax.set_title(comp, fontsize=7)
            if col == 0:
                ax.set_ylabel(label, fontsize=5, rotation=0, ha="right",
                              va="center", labelpad=30)
            style_ax(ax)

    fig.suptitle(f"Component scores by TP/FP — {end_label}", fontsize=8, y=1.02)
    fig.tight_layout()
    savefig(fig, outdir / f"component_violins_{end_label}")


# ── Panel 3: Joint heatmaps ────────────────────────────────────────────────

def plot_joint_heatmaps(
    data_by_label: dict[str, dict],
    end_label: str,
    outdir: Path,
) -> None:
    labels = list(data_by_label.keys())
    n = len(labels)
    if n == 0:
        return

    n_pairs = len(COMP_PAIRS)
    # Layout: rows=configs, cols=3 pairs × 2 (TP, FP side by side)
    fig, axes = plt.subplots(n, n_pairs,
                             figsize=(n_pairs * 2.5, max(W1, n * 2.2)),
                             squeeze=False)

    bins = np.linspace(0, 1, 26)

    for row, label in enumerate(labels):
        d = data_by_label[label]
        tp = d["tp"].astype(bool)

        for col, (c1, c2) in enumerate(COMP_PAIRS):
            ax = axes[row, col]
            x_tp = d[c1][tp]
            y_tp = d[c2][tp]
            x_fp = d[c1][~tp]
            y_fp = d[c2][~tp]

            # 2D histograms
            h_tp, _, _ = np.histogram2d(x_tp, y_tp, bins=bins)
            h_fp, _, _ = np.histogram2d(x_fp, y_fp, bins=bins)

            # Normalize to density
            h_tp_n = h_tp / max(h_tp.sum(), 1)
            h_fp_n = h_fp / max(h_fp.sum(), 1)

            # Log-ratio: log2(TP / FP), with pseudocounts
            pseudo = 1e-4
            ratio = np.log2((h_tp_n + pseudo) / (h_fp_n + pseudo))

            im = ax.imshow(ratio.T, origin="lower", aspect="auto",
                           extent=[0, 1, 0, 1],
                           cmap="RdBu", vmin=-3, vmax=3)

            if row == 0:
                ax.set_title(f"{c1} × {c2}", fontsize=7)
            if col == 0:
                ax.set_ylabel(label, fontsize=5, rotation=0, ha="right",
                              va="center", labelpad=30)
            style_ax(ax, xlabel=c1 if row == n - 1 else None)
            if row < n - 1:
                ax.set_xticklabels([])
            if col > 0:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel(label, fontsize=5, rotation=0, ha="right",
                              va="center", labelpad=30)

    # Colorbar
    cbar = fig.colorbar(im, ax=axes, shrink=0.6, pad=0.02, label="log₂(TP/FP)")
    cbar.ax.tick_params(labelsize=5)

    fig.suptitle(f"Joint score log-ratio (TP vs FP) — {end_label}",
                 fontsize=8, y=1.02)
    fig.tight_layout()
    savefig(fig, outdir / f"joint_heatmap_{end_label}")


# ── Panel 4: Signal vs score colored by TP/FP ──────────────────────────────

def plot_signal_vs_score(
    data_by_label: dict[str, dict],
    end_label: str,
    outdir: Path,
) -> None:
    labels = list(data_by_label.keys())
    n = len(labels)
    if n == 0:
        return

    comps = COMPONENTS
    n_comp = len(comps)

    fig, axes = plt.subplots(n, n_comp, figsize=(W2, max(W1, n * 2.0)),
                             squeeze=False, sharex="col")

    sig_label = "CAGE signal" if end_label == "5prime" else "QuantSeq signal"

    for row, label in enumerate(labels):
        d = data_by_label[label]
        tp = d["tp"].astype(bool)
        x = np.log1p(d["signal"])

        for col, comp in enumerate(comps):
            ax = axes[row, col]
            y = d[comp]

            # Plot FP first (background), TP on top
            ax.scatter(x[~tp], y[~tp], s=3, alpha=0.3, color=PALETTE[5],
                       edgecolors="none", rasterized=True, label="FP", zorder=2)
            ax.scatter(x[tp], y[tp], s=3, alpha=0.3, color=PALETTE[2],
                       edgecolors="none", rasterized=True, label="TP", zorder=3)

            if row == 0:
                ax.set_title(comp, fontsize=7)
            if row == n - 1:
                style_ax(ax, xlabel=f"log₁₊({sig_label})")
            else:
                style_ax(ax)
                ax.set_xticklabels([])
            if col == 0:
                ax.set_ylabel(label, fontsize=5, rotation=0, ha="right",
                              va="center", labelpad=30)
            if row == 0 and col == n_comp - 1:
                ax.legend(fontsize=5, loc="upper right", frameon=False,
                          markerscale=3)

    fig.suptitle(f"Signal vs component score by TP/FP — {end_label}",
                 fontsize=8, y=1.02)
    fig.tight_layout()
    savefig(fig, outdir / f"signal_vs_score_{end_label}")


# ── Panel 5: Weight sweep ───────────────────────────────────────────────────

def _sweep_weights(d: dict, n_steps: int = 21) -> dict:
    """Sweep (w_depth, w_model, w_annot) simplex, compute best F1 at each point.

    Returns dict with arrays: w_depth, w_model, w_annot, best_f1, best_thresh.
    """
    tp = d["tp"]
    depth = d["depth"]
    model = d["model"]
    annot = d["annot"]

    grid = np.linspace(0, 1, n_steps)
    w_ds, w_ms, w_as, f1s, thresholds = [], [], [], [], []

    for wd in grid:
        for wm in grid:
            wa = 1.0 - wd - wm
            if wa < -0.01 or wa > 1.01:
                continue
            wa = max(0.0, wa)

            composite = wd * depth + wm * model + wa * annot
            best_f1, best_t = 0.0, 0.5

            for t in np.linspace(0.05, 0.95, 19):
                pred = (composite >= t).astype(int)
                tp_hit = ((pred == 1) & (tp == 1)).sum()
                fp_hit = ((pred == 1) & (tp == 0)).sum()
                fn_hit = ((pred == 0) & (tp == 1)).sum()
                prec = tp_hit / (tp_hit + fp_hit) if (tp_hit + fp_hit) > 0 else 0
                rec = tp_hit / (tp_hit + fn_hit) if (tp_hit + fn_hit) > 0 else 0
                f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0
                if f1 > best_f1:
                    best_f1, best_t = f1, t

            w_ds.append(wd)
            w_ms.append(wm)
            w_as.append(wa)
            f1s.append(best_f1)
            thresholds.append(best_t)

    return {
        "w_depth": np.asarray(w_ds),
        "w_model": np.asarray(w_ms),
        "w_annot": np.asarray(w_as),
        "best_f1": np.asarray(f1s),
        "best_thresh": np.asarray(thresholds),
    }


def plot_weight_sweep(
    data_by_label: dict[str, dict],
    end_label: str,
    outdir: Path,
) -> dict:
    """Ternary-like 2D heatmap: x=w_model, y=w_depth (w_annot=1-x-y).

    Returns {label: {w_depth, w_model, w_annot, f1, threshold}} for optimal.
    """
    labels = list(data_by_label.keys())
    n = len(labels)
    if n == 0:
        return {}

    fig, axes = plt.subplots(1, n, figsize=(min(W2, n * 3.0), 3.0),
                             squeeze=False, sharey=True)
    axes = axes[0]

    optimal = {}

    for i, label in enumerate(labels):
        ax = axes[i]
        d = data_by_label[label]
        sw = _sweep_weights(d, n_steps=21)

        best_idx = np.argmax(sw["best_f1"])
        opt = {
            "w_depth": sw["w_depth"][best_idx],
            "w_model": sw["w_model"][best_idx],
            "w_annot": sw["w_annot"][best_idx],
            "f1": sw["best_f1"][best_idx],
            "threshold": sw["best_thresh"][best_idx],
        }
        optimal[label] = opt

        # 2D scatter (w_model on x, w_depth on y, color = F1)
        sc = ax.scatter(sw["w_model"], sw["w_depth"], c=sw["best_f1"],
                        cmap="viridis", s=12, edgecolors="none",
                        vmin=0, vmax=max(sw["best_f1"].max(), 0.5))
        ax.scatter([opt["w_model"]], [opt["w_depth"]], marker="*",
                   s=100, color="red", edgecolors="black", linewidths=0.5,
                   zorder=5)
        ax.text(0.03, 0.03,
                f"Best F1={opt['f1']:.3f}\n"
                f"d={opt['w_depth']:.2f} m={opt['w_model']:.2f} a={opt['w_annot']:.2f}\n"
                f"t={opt['threshold']:.2f}",
                transform=ax.transAxes, fontsize=5, va="bottom",
                bbox=dict(facecolor="white", alpha=0.8, pad=1, edgecolor="none"))

        ax.set_title(label, fontsize=6)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        # Shade invalid region w_depth + w_model > 1
        ax.fill_between([0, 1], [1, 0], [1.05, 1.05], color="grey", alpha=0.15)

        if i == 0:
            style_ax(ax, xlabel="w_model", ylabel="w_depth")
        else:
            style_ax(ax, xlabel="w_model")

    fig.colorbar(sc, ax=axes, shrink=0.8, pad=0.02, label="Best F1")
    fig.suptitle(f"Weight sweep (w_annot = 1 − w_depth − w_model) — {end_label}",
                 fontsize=8, y=1.02)
    fig.tight_layout()
    savefig(fig, outdir / f"weight_sweep_{end_label}")
    return optimal


# ── Panel 6: Marginal model value ───────────────────────────────────────────

def plot_marginal_value(
    auc_5prime: dict,
    auc_3prime: dict,
    outdir: Path,
) -> None:
    """Bar chart: AUC(depth+annot+model) - AUC(depth+annot) per end × config."""
    labels_5 = list(auc_5prime.keys())
    labels_3 = list(auc_3prime.keys())
    all_labels = list(dict.fromkeys(labels_5 + labels_3))  # preserve order
    if not all_labels:
        return

    fig, ax = plt.subplots(figsize=(W2, W1))
    x = np.arange(len(all_labels))
    w = 0.35

    deltas_5 = []
    deltas_3 = []
    for label in all_labels:
        if label in auc_5prime:
            a = auc_5prime[label]
            deltas_5.append(a.get("reality", 0.5) - a.get("depth+annot", 0.5))
        else:
            deltas_5.append(0)
        if label in auc_3prime:
            a = auc_3prime[label]
            deltas_3.append(a.get("reality", 0.5) - a.get("depth+annot", 0.5))
        else:
            deltas_3.append(0)

    ax.bar(x - w / 2, deltas_5, w, label="5′ (TSS)", color=PALETTE[4])
    ax.bar(x + w / 2, deltas_3, w, label="3′ (TTS)", color=PALETTE[0])
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")

    ax.set_xticks(x)
    ax.set_xticklabels(all_labels, fontsize=5, rotation=45, ha="right")
    ax.legend(fontsize=6, frameon=False)
    style_ax(ax, ylabel="ΔAUC (reality − depth+annot)",
             title="Marginal discriminative value of model score")
    fig.tight_layout()
    savefig(fig, outdir / "marginal_value_bar")


# ── Summary TSV ─────────────────────────────────────────────────────────────

def write_summary(
    auc_5: dict, auc_3: dict,
    opt_5: dict, opt_3: dict,
    data_5: dict, data_3: dict,
    outdir: Path,
) -> None:
    all_labels = list(dict.fromkeys(list(auc_5.keys()) + list(auc_3.keys())))
    rows = []
    for label in all_labels:
        row = {"config": label}
        for end, aucs, opts, data in [
            ("5prime", auc_5, opt_5, data_5),
            ("3prime", auc_3, opt_3, data_3),
        ]:
            if label in aucs:
                for comp, val in aucs[label].items():
                    row[f"auc_{end}_{comp}"] = f"{val:.4f}"
            if label in opts:
                o = opts[label]
                row[f"opt_{end}_w_depth"] = f"{o['w_depth']:.2f}"
                row[f"opt_{end}_w_model"] = f"{o['w_model']:.2f}"
                row[f"opt_{end}_w_annot"] = f"{o['w_annot']:.2f}"
                row[f"opt_{end}_f1"] = f"{o['f1']:.4f}"
                row[f"opt_{end}_thresh"] = f"{o['threshold']:.2f}"
            if label in data:
                d = data[label]
                row[f"n_{end}_tp"] = str(d["tp"].sum())
                row[f"n_{end}_fp"] = str(len(d["tp"]) - d["tp"].sum())
        rows.append(row)

    if not rows:
        return
    fieldnames = list(rows[0].keys())
    for r in rows[1:]:
        for k in r:
            if k not in fieldnames:
                fieldnames.append(k)

    path = outdir / "summary_stats.tsv"
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    log.info("Summary → %s", path)


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bed", nargs="+", required=True,
                        help="label:path pairs for BED12+TED isoform files")
    parser.add_argument("--cage-peaks", required=True, help="CAGE peaks BED6")
    parser.add_argument("--qs-peaks", required=True, help="QuantSeq peaks BED6")
    parser.add_argument("--cage-plus", required=True, help="CAGE bedGraph (+ strand)")
    parser.add_argument("--cage-minus", required=True, help="CAGE bedGraph (- strand)")
    parser.add_argument("--qs-plus", required=True, help="QuantSeq bedGraph (+ strand)")
    parser.add_argument("--qs-minus", required=True, help="QuantSeq bedGraph (- strand)")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Parse inputs ────────────────────────────────────────────────────
    beds_by_method: dict[str, List[dict]] = {}
    for entry in args.bed:
        label, path = entry.split(":", 1)
        isos = parse_bed_with_ted_scores(path)
        if isos:
            beds_by_method[label] = isos
            log.info("Loaded %d TED isoforms for %s", len(isos), label)
        else:
            log.info("Skipping %s (no TED columns)", label)

    if not beds_by_method:
        log.warning("No TED-format BED files found; nothing to do")
        sys.exit(0)

    cage_peaks = _parse_peaks_bed(args.cage_peaks)
    qs_peaks = _parse_peaks_bed(args.qs_peaks)
    log.info("Loaded %d / %d CAGE / QuantSeq peak groups",
             len(cage_peaks), len(qs_peaks))

    cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
        args.cage_plus, args.cage_minus, args.qs_plus, args.qs_minus,
    )

    # ── Run per end type ────────────────────────────────────────────────
    auc_5prime: dict = {}
    auc_3prime: dict = {}
    opt_5prime: dict = {}
    opt_3prime: dict = {}
    data_5prime: dict = {}
    data_3prime: dict = {}

    for end_type, end_label, peaks, auc_store, opt_store, data_store in [
        ("tss", "5prime", cage_peaks, auc_5prime, opt_5prime, data_5prime),
        ("tts", "3prime", qs_peaks, auc_3prime, opt_3prime, data_3prime),
    ]:
        data_by_label: dict[str, dict] = {}
        for label, isoforms in beds_by_method.items():
            d = _extract(isoforms, peaks, cage_p, cage_m, qs_p, qs_m, end_type)
            n_tp = d["tp"].sum()
            n_fp = len(d["tp"]) - n_tp
            if n_tp == 0 or n_fp == 0:
                log.warning("%s %s: TP=%d FP=%d — skipping (need both classes)",
                            label, end_label, n_tp, n_fp)
                continue
            data_by_label[label] = d
            log.info("%s %s: %d isoforms (TP=%d, FP=%d)",
                     label, end_label, len(d["tp"]), n_tp, n_fp)

        if not data_by_label:
            log.warning("No valid data for %s, skipping all panels", end_label)
            continue

        data_store.update(data_by_label)

        # Panel 1 — ROC curves
        aucs = plot_roc(data_by_label, end_label, outdir)
        auc_store.update(aucs)

        # Panel 2 — Violins
        plot_violins(data_by_label, end_label, outdir)

        # Panel 3 — Joint heatmaps
        plot_joint_heatmaps(data_by_label, end_label, outdir)

        # Panel 4 — Signal vs score by TP/FP
        plot_signal_vs_score(data_by_label, end_label, outdir)

        # Panel 5 — Weight sweep
        opts = plot_weight_sweep(data_by_label, end_label, outdir)
        opt_store.update(opts)

    # Panel 6 — Marginal model value (cross-end)
    plot_marginal_value(auc_5prime, auc_3prime, outdir)

    # Summary TSV
    write_summary(auc_5prime, auc_3prime, opt_5prime, opt_3prime,
                  data_5prime, data_3prime, outdir)

    log.info("Done — output in %s", outdir)


if __name__ == "__main__":
    main()
