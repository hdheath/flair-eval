#!/usr/bin/env python3
"""
isoform_end_roc.py — AUC-ROC curves for isoform end-calling accuracy.

Each isoform is treated as a binary classifier prediction:
  Score  = orthogonal signal at the called end (CAGE at TSS, dRNA at TTS).
  Label  = 1 (TP) if the called end is within --window bp of a reference peak,
           0 (FP) otherwise.

ROC curves are plotted with TPR (sensitivity) on Y and FPR (1-specificity) on X,
sweeping the signal threshold from high to low.  AUC is computed via the
trapezoidal rule and reported in the legend and in a summary TSV.

JC-deduplication is applied to the TP labels:
  Within each junction chain group, if multiple isoforms map to the SAME
  reference peak, only ONE counts as a TP.  This is consistent with
  compute_jc_deduplicated_precision_recall() used throughout the pipeline
  and prevents over-segmented assemblers from inflating their ROC curves.

Paired ROC:
  Score  = geometric mean of CAGE signal at TSS and dRNA signal at TTS.
           This is the joint signal discriminating full-isoform correctness.
  Label  = 1 (joint TP) only if BOTH ends independently hit their respective
           orthogonal peaks within --window bp.  JC-dedup applied: within each
           JC group, only ONE isoform may claim a given (tss_peak, tts_peak) pair
           as a TP.

Outputs (all in --output dir):
  roc_5prime.png          — TSS ROC curves, all methods overlaid
  roc_3prime.png          — TTS ROC curves, all methods overlaid
  roc_combined.png        — 5' and 3' side-by-side
  roc_paired.png          — Paired joint-TP ROC curves, all methods overlaid
  roc_auc_summary.tsv     — mode, end, AUC, n_isoforms, n_tp, n_fp

Usage:
    python isoform_end_roc.py \\
        --bed          label1:bed1.bed label2:bed2.bed ... \\
        --cage-peaks   cage_peaks.bed \\
        --drna-peaks   drna_peaks.bed \\
        --cage-plus    cage_plus.bg   --cage-minus  cage_minus.bg \\
        --qs-plus      qs_plus.bg     --qs-minus    qs_minus.bg \\
        --output       output_dir/ \\
        [--window 50]  [--verbose]
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler
    from signal_utils import parse_isoforms, load_signal_tracks, BedGraphTrack, isoform_signal
    from ted_end_precision import parse_peaks_bed, _nearest_annot
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler
    from evaluation.signal_utils import parse_isoforms, load_signal_tracks, BedGraphTrack, isoform_signal
    from evaluation.ted_end_precision import parse_peaks_bed, _nearest_annot

apply_rc()

# ── JC-dedup TP labelling ─────────────────────────────────────────────────────

def _jc_dedup_labels(
    isoforms: List[dict],
    peaks: Dict[Tuple[str, str], List[Tuple[int, int]]],
    end: str,   # "tss" or "tts"
    window: int,
) -> List[Tuple[float, int]]:
    """Return list of (signal_score, jc_dedup_tp_label) for every isoform.

    Labels are JC-deduplicated: within a junction chain group, if multiple
    isoforms map to the same peak, only ONE counts as a TP.  All others in that
    group that also match the same peak receive label=0 (they are redundant).

    signal_score is NaN when the bedGraph returns zero (excluded from ROC).
    Isoforms on chromosomes/strands absent from the peak file are labelled FP.
    """
    # Group by junction chain
    jc_groups: Dict[tuple, List[dict]] = defaultdict(list)
    for iso in isoforms:
        if iso["n_exons"] >= 2 and iso["junctions"]:
            key = (iso["chrom"], iso["strand"], iso["junctions"])
        else:
            # Single-exon: use rounded coordinates as group key (100 bp bins)
            rnd = lambda v: 100 * round(v / 100)
            key = (iso["chrom"], iso["strand"], "__se__",
                   rnd(iso["start"]), rnd(iso["end"]))
        jc_groups[key].append(iso)

    results: List[Tuple[float, int]] = []

    for key, members in jc_groups.items():
        # Which peaks has this JC group already claimed?
        claimed_peaks: set = set()

        for iso in members:
            ch, strand = iso["chrom"], iso["strand"]
            pos = iso["tss"] if end == "tss" else iso["tts"]
            sig = iso.get(f"_sig_{end}", np.nan)
            if sig == 0.0:
                sig = np.nan  # treat zero signal as missing (no bedGraph coverage)

            peak_iv = _nearest_annot(
                pos,
                peaks.get((ch, strand), []),
                window,
            )

            if peak_iv is None:
                # No matching peak → FP
                results.append((sig, 0))
            elif peak_iv in claimed_peaks:
                # Same peak already claimed by another isoform in this JC group → redundant FP
                results.append((sig, 0))
            else:
                claimed_peaks.add(peak_iv)
                results.append((sig, 1))

    return results


# ── Paired JC-dedup TP labelling ─────────────────────────────────────────────

def _jc_dedup_paired_labels(
    isoforms: List[dict],
    cage_peaks: Dict[Tuple[str, str], List[Tuple[int, int]]],
    drna_peaks: Dict[Tuple[str, str], List[Tuple[int, int]]],
    window: int,
) -> List[Tuple[float, int]]:
    """Return list of (joint_score, joint_tp_label) for every isoform.

    joint_score = geometric mean of CAGE signal at TSS and dRNA signal at TTS.
                  NaN if either signal is zero/missing.
    joint_tp    = 1 only if BOTH ends hit their respective orthogonal peaks AND
                  the (tss_peak, tts_peak) pair has not already been claimed by
                  another isoform in the same JC group (JC-dedup).
    """
    jc_groups: Dict[tuple, List[dict]] = defaultdict(list)
    for iso in isoforms:
        if iso["n_exons"] >= 2 and iso["junctions"]:
            key = (iso["chrom"], iso["strand"], iso["junctions"])
        else:
            rnd = lambda v: 100 * round(v / 100)
            key = (iso["chrom"], iso["strand"], "__se__",
                   rnd(iso["start"]), rnd(iso["end"]))
        jc_groups[key].append(iso)

    results: List[Tuple[float, int]] = []

    for key, members in jc_groups.items():
        claimed_pairs: set = set()

        for iso in members:
            ch, strand = iso["chrom"], iso["strand"]
            sig_tss = iso.get("_sig_tss", np.nan)
            sig_tts = iso.get("_sig_tts", np.nan)
            # Geometric mean — zero if either end is missing
            if sig_tss == 0.0 or sig_tts == 0.0 or np.isnan(sig_tss) or np.isnan(sig_tts):
                joint_score = np.nan
            else:
                joint_score = float(np.sqrt(sig_tss * sig_tts))

            cage_iv = _nearest_annot(iso["tss"], cage_peaks.get((ch, strand), []), window)
            drna_iv = _nearest_annot(iso["tts"], drna_peaks.get((ch, strand), []), window)

            if cage_iv is None or drna_iv is None:
                results.append((joint_score, 0))
            else:
                pair = (cage_iv, drna_iv)
                if pair in claimed_pairs:
                    results.append((joint_score, 0))
                else:
                    claimed_pairs.add(pair)
                    results.append((joint_score, 1))

    return results


# ── ROC computation ───────────────────────────────────────────────────────────

def _roc_curve(
    scores_labels: List[Tuple[float, int]],
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Compute ROC curve and AUC from (score, label) pairs.

    Pairs with NaN scores are excluded.  Returns (fpr, tpr, auc).
    """
    valid = [(s, l) for s, l in scores_labels if not np.isnan(s)]
    if not valid:
        return np.array([0.0, 1.0]), np.array([0.0, 1.0]), 0.5

    scores = np.array([s for s, _ in valid], dtype=float)
    labels = np.array([l for _, l in valid], dtype=int)

    n_pos = labels.sum()
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return np.array([0.0, 1.0]), np.array([0.0, 1.0]), 0.5

    # Sort by score descending (highest signal = most likely TP first)
    order = np.argsort(-scores)
    labels_sorted = labels[order]

    tps = np.cumsum(labels_sorted)
    fps = np.cumsum(1 - labels_sorted)

    tpr = np.concatenate([[0.0], tps / n_pos])
    fpr = np.concatenate([[0.0], fps / n_neg])

    # Trapezoidal AUC
    auc = float(np.trapz(tpr, fpr))
    return fpr, tpr, auc


# ── Per-method data collection ────────────────────────────────────────────────

def collect_roc_data(
    beds_by_method: Dict[str, List[dict]],
    cage_peaks: Dict[Tuple[str, str], List[Tuple[int, int]]],
    drna_peaks: Dict[Tuple[str, str], List[Tuple[int, int]]],
    window: int,
) -> Dict[str, dict]:
    """Build ROC data per method.

    Returns:
      {method: {
          "5prime": {"fpr": np.ndarray, "tpr": np.ndarray, "auc": float,
                     "n": int, "n_tp": int, "n_fp": int},
          "3prime": {...},
      }}
    """
    results: Dict[str, dict] = {}

    for method, isoforms in beds_by_method.items():
        entry: dict = {}

        for end, peaks in [("5prime", cage_peaks), ("3prime", drna_peaks)]:
            end_key = "tss" if end == "5prime" else "tts"
            pairs = _jc_dedup_labels(isoforms, peaks, end_key, window)
            n       = len(pairs)
            n_tp    = sum(l for s, l in pairs if not np.isnan(s))
            n_fp    = sum(1 - l for s, l in pairs if not np.isnan(s))
            fpr, tpr, auc = _roc_curve(pairs)
            entry[end] = {
                "fpr": fpr, "tpr": tpr, "auc": auc,
                "n": n, "n_tp": int(n_tp), "n_fp": int(n_fp),
            }

        # Paired: joint_tp label, geometric-mean signal score
        paired_pairs = _jc_dedup_paired_labels(isoforms, cage_peaks, drna_peaks, window)
        n_p    = len(paired_pairs)
        n_tp_p = sum(l for s, l in paired_pairs if not np.isnan(s))
        n_fp_p = sum(1 - l for s, l in paired_pairs if not np.isnan(s))
        fpr_p, tpr_p, auc_p = _roc_curve(paired_pairs)
        entry["paired"] = {
            "fpr": fpr_p, "tpr": tpr_p, "auc": auc_p,
            "n": n_p, "n_tp": int(n_tp_p), "n_fp": int(n_fp_p),
        }

        results[method] = entry

    return results


# ── Plotting ─────────────────────────────────────────────────────────────────

def _plot_roc_panel(
    ax: plt.Axes,
    roc_data: Dict[str, dict],
    end: str,           # "5prime" or "3prime"
    styler: ModeStyler,
    title: str,
) -> None:
    """Draw ROC curves for all methods onto ax."""
    # Random-classifier diagonal
    ax.plot([0, 1], [0, 1], color="#cccccc", linewidth=0.6,
            linestyle="--", zorder=0)

    # Sort methods by AUC descending so legend reads top-to-bottom
    methods_sorted = sorted(
        roc_data.keys(),
        key=lambda m: roc_data[m][end]["auc"],
        reverse=True,
    )

    handles = []
    for method in methods_sorted:
        d = roc_data[method][end]
        auc = d["auc"]
        fpr, tpr = d["fpr"], d["tpr"]
        dash = styler.dash(method)
        ls = "-" if dash == (1, 0) else (0, dash)
        line, = ax.plot(
            fpr, tpr,
            color=styler.color(method),
            linestyle=ls,
            linewidth=1.1,
            zorder=2,
        )
        label = f"{_short(method)}  AUC={auc:.3f}  (n={d['n']:,})"
        handles.append(styler.legend_handle(method, label=label, markersize=0))

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("FPR  (1 − specificity)", fontsize=7)
    ax.set_ylabel("TPR  (sensitivity)", fontsize=7)
    ax.set_title(title, fontsize=7)
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1))
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1))
    style_ax(ax)

    ax.legend(handles=handles, fontsize=5.5, frameon=False,
              loc="lower right", handlelength=1.8)


def plot_roc(
    roc_data: Dict[str, dict],
    styler: ModeStyler,
    out_5prime: Path,
    out_3prime: Path,
    out_combined: Path,
    out_paired: Path,
) -> None:
    """Produce four output figures."""
    # --- 5' only ---
    fig, ax = plt.subplots(1, 1, figsize=(W1 * 1.15, W1 * 1.15))
    _plot_roc_panel(ax, roc_data, "5prime", styler, "5′ TSS end accuracy  (CAGE signal)")
    fig.tight_layout(pad=0.4)
    savefig(fig, out_5prime)

    # --- 3' only ---
    fig, ax = plt.subplots(1, 1, figsize=(W1 * 1.15, W1 * 1.15))
    _plot_roc_panel(ax, roc_data, "3prime", styler, "3′ TTS end accuracy  (dRNA signal)")
    fig.tight_layout(pad=0.4)
    savefig(fig, out_3prime)

    # --- Combined side-by-side ---
    fig, axes = plt.subplots(1, 2, figsize=(W2, W1 * 1.15))
    _plot_roc_panel(axes[0], roc_data, "5prime", styler, "5′ TSS  (CAGE)")
    _plot_roc_panel(axes[1], roc_data, "3prime", styler, "3′ TTS  (dRNA)")
    fig.tight_layout(pad=0.4, w_pad=1.0)
    savefig(fig, out_combined)

    # --- Paired joint-TP ---
    fig, ax = plt.subplots(1, 1, figsize=(W1 * 1.15, W1 * 1.15))
    _plot_roc_panel(ax, roc_data, "paired", styler,
                   "Paired joint-TP  (√CAGE·dRNA signal,  both ends must hit peaks)")
    fig.tight_layout(pad=0.4)
    savefig(fig, out_paired)


# ── TSV summary ───────────────────────────────────────────────────────────────

def write_auc_tsv(roc_data: Dict[str, dict], path: Path) -> None:
    """Write per-method per-end AUC summary TSV."""
    with open(path, "w") as f:
        f.write("mode\tend\tauc\tn_isoforms\tn_tp\tn_fp\n")
        for method in roc_data:
            for end in ("5prime", "3prime", "paired"):
                d = roc_data[method].get(end)
                if d is None:
                    continue
                f.write(
                    f"{method}\t{end}\t{d['auc']:.6f}\t"
                    f"{d['n']}\t{d['n_tp']}\t{d['n_fp']}\n"
                )


# ── Helpers ───────────────────────────────────────────────────────────────────

def _short(mode: str) -> str:
    return (mode
            .replace("TED-", "")
            .replace("FLAIR-", "FL-")
            .replace("isoquant_", "IQ-"))


def _parse_label_path(entries: List[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for e in entries:
        if ":" not in e:
            continue
        label, path = e.split(":", 1)
        out[label] = path
    return out


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--bed",         nargs="+", required=True,
                   help="label:path pairs for BED12/GTF isoform files")
    p.add_argument("--cage-peaks",  default=None,
                   help="CAGE peak BED6 file (5′ reference)")
    p.add_argument("--drna-peaks",  default=None,
                   help="dRNA / QuantSeq peak BED6 file (3′ reference)")
    p.add_argument("--cage-plus",   required=True,
                   help="CAGE bedGraph, plus strand")
    p.add_argument("--cage-minus",  required=True,
                   help="CAGE bedGraph, minus strand")
    p.add_argument("--qs-plus",     required=True,
                   help="dRNA/QuantSeq bedGraph, plus strand")
    p.add_argument("--qs-minus",    required=True,
                   help="dRNA/QuantSeq bedGraph, minus strand")
    p.add_argument("--output",      required=True,
                   help="Output directory")
    p.add_argument("--window",      type=int, default=50,
                   help="bp tolerance for end-to-peak matching (default 50)")
    p.add_argument("--verbose",     action="store_true")
    args = p.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    # --- Load signal tracks ---
    if args.verbose:
        print("Loading signal tracks...", file=sys.stderr)
    cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
        args.cage_plus, args.cage_minus, args.qs_plus, args.qs_minus,
    )

    # --- Load peak files ---
    if not args.cage_peaks and not args.drna_peaks:
        print("WARNING: neither --cage-peaks nor --drna-peaks provided — "
              "cannot compute ROC. Exiting.", file=sys.stderr)
        sys.exit(0)

    if args.verbose:
        print("Loading peak files...", file=sys.stderr)
    cage_peaks = parse_peaks_bed(args.cage_peaks) if args.cage_peaks else {}
    drna_peaks = parse_peaks_bed(args.drna_peaks) if args.drna_peaks else {}

    if not cage_peaks and not drna_peaks:
        print("WARNING: peak files were provided but no peaks loaded — "
              "files may be empty. Exiting.", file=sys.stderr)
        sys.exit(0)

    # --- Parse isoform BED/GTF files and attach signal scores ---
    bed_paths = _parse_label_path(args.bed)
    beds_by_method: Dict[str, List[dict]] = {}

    for label, path in bed_paths.items():
        if not Path(path).exists():
            print(f"WARNING: {path} not found — skipping {label}", file=sys.stderr)
            continue
        isos = parse_isoforms(path)
        if not isos:
            continue

        # Attach CAGE signal at TSS and dRNA signal at TTS to each isoform
        for iso in isos:
            iso["_sig_tss"] = isoform_signal(iso, cage_p, cage_m, "tss")
            iso["_sig_tts"] = isoform_signal(iso, qs_p,   qs_m,   "tts")

        beds_by_method[label] = isos
        if args.verbose:
            print(f"  {label}: {len(isos)} isoforms", file=sys.stderr)

    if not beds_by_method:
        print("No isoform data loaded — exiting.", file=sys.stderr)
        sys.exit(1)

    # --- Compute ROC data ---
    if args.verbose:
        print("Computing JC-deduplicated ROC curves...", file=sys.stderr)
    roc_data = collect_roc_data(beds_by_method, cage_peaks, drna_peaks, args.window)

    if args.verbose:
        for method, d in roc_data.items():
            print(
                f"  {method}: 5′ AUC={d['5prime']['auc']:.3f} "
                f"(n={d['5prime']['n']})  "
                f"3′ AUC={d['3prime']['auc']:.3f} "
                f"(n={d['3prime']['n']})  "
                f"paired AUC={d['paired']['auc']:.3f} "
                f"(n_tp={d['paired']['n_tp']})",
                file=sys.stderr,
            )

    # --- Plot ---
    styler = ModeStyler(list(roc_data.keys()))
    plot_roc(
        roc_data, styler,
        out_5prime   = out / "roc_5prime.png",
        out_3prime   = out / "roc_3prime.png",
        out_combined = out / "roc_combined.png",
        out_paired   = out / "roc_paired.png",
    )

    # --- TSV summary ---
    write_auc_tsv(roc_data, out / "roc_auc_summary.tsv")

    print(f"Saved ROC outputs to {args.output}")


if __name__ == "__main__":
    main()
