#!/usr/bin/env python3
"""
ted_score_vs_signal.py — Compare TED per-isoform scores to orthogonal boundary signal.

For each TED configuration, scatter-plots each score component (depth, model,
annotation proximity, composite reality) against orthogonal boundary signal
(CAGE for 5′ TSS, dRNA for 3′ TTS).  Non-TED tools (no extra BED columns)
are silently skipped.

Outputs per end type (5prime / 3prime):
  ted_score_scatter_{end}.{png,svg}   — 4-column scatter grid (rows=configs)
  ted_score_correlation_{end}.{png,svg} — grouped bar of Spearman ρ per component

Inputs:
  --bed       label:path  (BED12/22 isoform files, one per assembler/config)
  --cage-plus/--cage-minus/--qs-plus/--qs-minus  bedGraph signal tracks
  --output    output directory
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
from scipy.stats import spearmanr, gaussian_kde

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

# TED extra-column mapping (0-indexed from col 12)
# Actual order from flair/isoform_data.py TED_SCORE_COLUMNS:
# col12=TED_confidence, col13=TED_tss_reality, col14=TED_tts_reality,
# col15=TED_depth (shared), col16=TED_tss_model, col17=TED_tts_model,
# col18=TED_tss_annot, col19=TED_tts_annot, col20=TED_tss_annot_dist, col21=TED_tts_annot_dist
TED_SCORE_KEYS = [
    "TED_confidence", "TED_tss_reality", "TED_tts_reality",
    "TED_depth", "TED_tss_model", "TED_tts_model",
    "TED_tss_annot", "TED_tts_annot", "TED_tss_annot_dist", "TED_tts_annot_dist",
]

SCORE_COMPONENTS = ["depth", "model", "annot", "reality"]


# ── BED parser with TED scores ─────────────────────────────────────────────

def parse_bed_with_ted_scores(path: str) -> List[dict]:
    """Parse BED12+ file.  Returns list of isoform dicts with TED scores
    when extra columns are present (cols 13–22), otherwise empty list."""
    isoforms: list[dict] = []
    has_ted = None
    with open(path) as f:
        for line in f:
            if line.startswith(("#", "track")):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 12:
                continue
            # Detect TED output by presence of extra columns
            if has_ted is None:
                has_ted = len(c) >= 20  # at least 8 score cols
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
            # Parse TED scores
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


# ── Per-config data extraction ──────────────────────────────────────────────

def _extract_score_signal(
    isoforms: List[dict],
    cage_p, cage_m, qs_p, qs_m,
    end_type: str,
) -> dict:
    """For one assembler, extract arrays of each score component + boundary signal.

    Returns dict with keys: 'depth', 'model', 'annot', 'reality', 'signal' —
    each a numpy array.
    """
    prefix = "tss" if end_type == "tss" else "tts"
    depths, models, annots, realities, signals = [], [], [], [], []

    for iso in isoforms:
        sig_tss, sig_tts = isoform_signal(iso, cage_p, cage_m, qs_p, qs_m)
        sig = sig_tss if end_type == "tss" else sig_tts

        depths.append(iso.get("TED_depth", 0.0))  # depth is shared across ends
        models.append(iso.get(f"TED_{prefix}_model", 0.0))
        annots.append(iso.get(f"TED_{prefix}_annot", 0.0))
        realities.append(iso.get(f"TED_{prefix}_reality", 0.0))
        signals.append(sig)

    return {
        "depth": np.array(depths),
        "model": np.array(models),
        "annot": np.array(annots),
        "reality": np.array(realities),
        "signal": np.array(signals),
    }


# ── Scatter grid plot ───────────────────────────────────────────────────────

def plot_scatter_grid(
    data_by_label: dict[str, dict],
    end_label: str,
    outdir: Path,
) -> None:
    """4-column scatter grid: rows = TED configs, cols = score components."""
    labels = list(data_by_label.keys())
    n_rows = len(labels)
    if n_rows == 0:
        return

    fig, axes = plt.subplots(
        n_rows, 4,
        figsize=(W2, max(W1, n_rows * 1.5)),
        squeeze=False, sharex="col",
    )

    sig_label = "CAGE signal" if end_label == "5prime" else "dRNA signal"

    for row, label in enumerate(labels):
        d = data_by_label[label]
        for col, comp in enumerate(SCORE_COMPONENTS):
            ax = axes[row, col]
            x = d[comp]
            y = d["signal"]

            # Scatter
            ax.scatter(x, y, s=3, alpha=0.3, color=PALETTE[col],
                       edgecolors="none", rasterized=True)

            # Spearman ρ (skip if constant)
            if len(x) > 2 and np.std(x) > 0 and np.std(y) > 0:
                rho, pval = spearmanr(x, y)
                ax.text(0.97, 0.97, f"ρ={rho:.2f}",
                        transform=ax.transAxes, ha="right", va="top",
                        fontsize=5, fontstyle="italic")

            # Labels
            if row == 0:
                ax.set_title(comp, fontsize=7)
            if row == n_rows - 1:
                style_ax(ax, xlabel=f"TED {comp} score")
            else:
                style_ax(ax)
            if col == 0:
                ax.set_ylabel(label, fontsize=5, rotation=0, ha="right",
                              va="center", labelpad=30)

    fig.supylabel(sig_label, fontsize=7, x=0.02)
    fig.suptitle(f"TED scores vs boundary signal ({end_label})",
                 fontsize=8, y=1.01)
    fig.tight_layout()
    savefig(fig, outdir / f"ted_score_scatter_{end_label}")


# ── Boundary signal colored by score components ─────────────────────────────

COLOR_COLS = ["density", "depth", "model", "annot", "reality"]
COMPONENT_CMAPS = {
    "density": "magma",
    "depth":   "viridis",
    "model":   "viridis",
    "annot":   "viridis",
    "reality": "viridis",
}


def _kde_colors(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return log-density values for KDE coloring."""
    mask = np.isfinite(x) & np.isfinite(y)
    vals = np.full(len(x), np.nan)
    if mask.sum() < 10:
        return vals
    xy = np.vstack([x[mask], y[mask]])
    try:
        kde = gaussian_kde(xy)
        vals[mask] = np.log1p(kde(xy))
    except np.linalg.LinAlgError:
        vals[mask] = 0.0
    return vals


def plot_colored_scatter(
    data_by_label: dict[str, dict],
    end_label: str,
    outdir: Path,
) -> None:
    """5-column scatter: x=signal, y=reality, color=density/depth/model/annot/reality."""
    labels = list(data_by_label.keys())
    n_rows = len(labels)
    if n_rows == 0:
        return

    n_cols = len(COLOR_COLS)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(n_cols * 2.2, max(W1, n_rows * 1.8)),
        squeeze=False,
    )

    sig_label = "CAGE signal" if end_label == "5prime" else "dRNA signal"

    # Compute global signal range for shared x-axis
    all_sig = np.concatenate([d["signal"] for d in data_by_label.values()])
    sig_max = np.percentile(all_sig[all_sig > 0], 99) if np.any(all_sig > 0) else 1.0

    for row, label in enumerate(labels):
        d = data_by_label[label]
        x_raw = d["signal"]
        y_raw = d["reality"]

        # log1p for signal (heavily skewed)
        x = np.log1p(x_raw)
        x_max = np.log1p(sig_max)

        for col_i, ccol in enumerate(COLOR_COLS):
            ax = axes[row, col_i]

            # Determine color values
            if ccol == "density":
                c = _kde_colors(x, y_raw)
            else:
                c = d[ccol]

            cmap = COMPONENT_CMAPS[ccol]

            # Sort so high values on top
            order = np.argsort(c)
            xs, ys, cs = x[order], y_raw[order], c[order]

            # Normalize
            finite = cs[np.isfinite(cs)]
            if len(finite) > 0 and np.ptp(finite) > 0:
                norm = Normalize(vmin=np.nanmin(finite), vmax=np.nanmax(finite))
            else:
                norm = Normalize(0, 1)

            sc = ax.scatter(
                xs, ys, c=cs, cmap=cmap, norm=norm,
                s=2, alpha=0.5, edgecolors="none", rasterized=True,
            )

            # Spearman ρ annotation (reality vs signal)
            if col_i == 0 and len(x_raw) > 2 and np.std(x_raw) > 0 and np.std(y_raw) > 0:
                rho, _ = spearmanr(x_raw, y_raw)
                ax.text(0.97, 0.03, f"ρ={rho:.2f}",
                        transform=ax.transAxes, ha="right", va="bottom",
                        fontsize=5, fontstyle="italic",
                        bbox=dict(facecolor="white", alpha=0.7, pad=1, edgecolor="none"))

            # Axis limits
            ax.set_xlim(-0.05 * x_max, x_max * 1.05)
            ax.set_ylim(-0.02, 1.02)

            # Column titles
            if row == 0:
                ax.set_title(f"color = {ccol}", fontsize=7)
            # X-label on bottom row only
            if row == n_rows - 1:
                style_ax(ax, xlabel=f"log₁₊({sig_label})")
            else:
                style_ax(ax)
                ax.set_xticklabels([])
            # Row label
            if col_i == 0:
                ax.set_ylabel(label, fontsize=5, rotation=0, ha="right",
                              va="center", labelpad=35)
            else:
                ax.set_yticklabels([])

    fig.supylabel("TED reality score", fontsize=7, x=0.01)
    fig.suptitle(f"Boundary signal colored by TED components ({end_label})",
                 fontsize=8, y=1.01)
    fig.tight_layout()
    savefig(fig, outdir / f"boundary_signal_colored_{end_label}")


# ── Correlation summary bar chart ───────────────────────────────────────────

def plot_correlation_summary(
    data_by_label: dict[str, dict],
    end_label: str,
    outdir: Path,
) -> None:
    """Grouped bar: Spearman ρ (TED score component vs signal) per config."""
    labels = list(data_by_label.keys())
    if not labels:
        return

    n_groups = len(SCORE_COMPONENTS)
    n_bars = len(labels)
    bar_w = 0.8 / n_bars
    x = np.arange(n_groups)

    fig, ax = plt.subplots(figsize=(W2, W1))

    for i, label in enumerate(labels):
        d = data_by_label[label]
        rhos = []
        for comp in SCORE_COMPONENTS:
            xv, yv = d[comp], d["signal"]
            if len(xv) > 2 and np.std(xv) > 0 and np.std(yv) > 0:
                rho, _ = spearmanr(xv, yv)
            else:
                rho = 0.0
            rhos.append(rho)
        offset = (i - n_bars / 2 + 0.5) * bar_w
        ax.bar(x + offset, rhos, bar_w, label=label, alpha=0.85,
               color=PALETTE[i % len(PALETTE)])

    ax.set_xticks(x)
    ax.set_xticklabels(SCORE_COMPONENTS)
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.legend(fontsize=5, loc="best", frameon=False, ncol=max(1, n_bars // 4))
    sig_label = "CAGE signal" if end_label == "5prime" else "dRNA signal"
    style_ax(ax, ylabel=f"Spearman ρ (vs {sig_label})",
             title=f"TED score–signal correlation ({end_label})")
    fig.tight_layout()
    savefig(fig, outdir / f"ted_score_correlation_{end_label}")


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bed", nargs="+", required=True,
                        help="label:path pairs for BED12 isoform files")
    parser.add_argument("--cage-plus", required=True)
    parser.add_argument("--cage-minus", required=True)
    parser.add_argument("--qs-plus", required=True)
    parser.add_argument("--qs-minus", required=True)
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load signal tracks
    cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
        args.cage_plus, args.cage_minus, args.qs_plus, args.qs_minus,
    )

    # Parse BED files, skip non-TED
    ted_isoforms: dict[str, list[dict]] = {}
    for entry in args.bed:
        label, path = entry.split(":", 1)
        isos = parse_bed_with_ted_scores(path)
        if not isos:
            log.info("Skipping %s (no TED score columns)", label)
            continue
        ted_isoforms[label] = isos
        log.info("Loaded %d isoforms with TED scores for %s", len(isos), label)

    if not ted_isoforms:
        log.warning("No TED-scored isoform files found, nothing to plot")
        return

    # Generate plots per end type
    for end_type, end_label in [("tss", "5prime"), ("tts", "3prime")]:
        data_by_label: dict[str, dict] = {}
        for label, isos in ted_isoforms.items():
            data_by_label[label] = _extract_score_signal(
                isos, cage_p, cage_m, qs_p, qs_m, end_type,
            )
        plot_scatter_grid(data_by_label, end_label, outdir)
        plot_correlation_summary(data_by_label, end_label, outdir)
        plot_colored_scatter(data_by_label, end_label, outdir)

    log.info("Done — output in %s", outdir)


if __name__ == "__main__":
    main()
