#!/usr/bin/env python3
"""
True-positive overlap analysis: baseline vs. each other method.

For every non-baseline method, computes Jaccard-style peak overlap metrics
against the baseline method and produces per-end-type bar charts.

Definitions (for a given method M vs baseline B):
    shared_TPs  = peaks recovered by BOTH B and M
    novel_TPs   = peaks recovered by M but NOT B
    lost_TPs    = peaks recovered by B but NOT M
    total_TPs_M = shared_TPs + novel_TPs
    total_TPs_B = shared_TPs + lost_TPs
    total_peaks = all peaks in the peak set

Metrics:
    novelty_rate      = novel_TPs / total_TPs_M
    incremental_recall = novel_TPs / total_peaks
    swap_rate         = lost_TPs  / total_TPs_B

Input:
    Peak-reason TSV files (peak_id, score, read_count, reason) — one per
    method.  Pass as  label:path  pairs, grouped by end type.

Output (per end type):
    novelty_rate_{end}.png         — bar chart of novelty rate
    incremental_recall_{end}.png   — bar chart of incremental recall
    swap_rate_{end}.png            — bar chart of swap rate
    tp_overlap_summary_{end}.png   — grouped bar chart (shared/novel/lost)
    tp_overlap_metrics_{end}.tsv   — numeric table

Usage:
    python tp_overlap_plot.py \\
        --cage baseline:peaks_bl.tsv trust-end:peaks_te.tsv ... \\
        --drna baseline:peaks_bl.tsv trust-end:peaks_te.tsv ... \\
        --output output_dir/ \\
        [--baseline-label baseline] \\
        [--title-prefix "v3 — "] [--verbose]
"""

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Attempt to import shared pub_style; fall back gracefully.
try:
    from pub_style import (
        MODE_COLORS, PALETTE, apply_rc, style_ax, legend_outside, savefig,
        W1, W2,
    )
except ImportError:
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from pub_style import (
            MODE_COLORS, PALETTE, apply_rc, style_ax, legend_outside, savefig,
            W1, W2,
        )
    except ImportError:
        PALETTE = ['#0072B2', '#E69F00', '#009E73', '#D55E00',
                    '#CC79A7', '#56B4E9', '#F0E442', '#999999']
        MODE_COLORS = {}
        W1, W2 = 3.5, 7.2
        def apply_rc(): pass
        def style_ax(ax, **kw):
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        def legend_outside(fig_or_ax, **kw):
            return fig_or_ax.legend(**kw)
        def savefig(fig, path, **kw):
            path = Path(path)
            for fmt in ("png", "svg"):
                fig.savefig(path.with_suffix(f".{fmt}"), dpi=300,
                            bbox_inches="tight")
            plt.close(fig)

apply_rc()

logger = logging.getLogger(__name__)

# Overlap-specific colours
OVERLAP_COLORS = {
    "shared": "#009E73",   # bluish green — mutual recovery
    "novel":  "#0072B2",   # blue — gained by method
    "lost":   "#D55E00",   # vermillion — lost vs baseline
}


# ── Data loading ────────────────────────────────────────────────────────────

def _parse_label_path(spec: str) -> Tuple[str, Path]:
    """Parse 'label:path' or infer label from filename."""
    if ":" in spec:
        label, path_str = spec.split(":", 1)
        return label, Path(path_str)
    p = Path(spec)
    label = p.stem.split("_")[0]
    return label, p


def load_peak_reason_tsv(path: Path) -> Dict[str, str]:
    """Load peak reason TSV, return {peak_id: reason}."""
    result = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            result[row["peak_id"]] = row["reason"]
    return result


def get_tp_set(peak_reasons: Dict[str, str]) -> Set[str]:
    """Return set of peak_ids that are true positives (reason == 'recovered')."""
    return {pid for pid, reason in peak_reasons.items() if reason == "recovered"}


# ── Overlap computation ────────────────────────────────────────────────────

def compute_overlap(
    baseline_tps: Set[str],
    method_tps: Set[str],
    total_peaks: int,
) -> Dict[str, float]:
    """Compute overlap metrics between baseline and another method.

    Returns dict with keys: shared, novel, lost, total_method, total_baseline,
    total_peaks, novelty_rate, incremental_recall, swap_rate, jaccard.
    """
    shared = baseline_tps & method_tps
    novel = method_tps - baseline_tps
    lost = baseline_tps - method_tps

    total_m = len(method_tps)
    total_b = len(baseline_tps)

    novelty = len(novel) / total_m if total_m > 0 else 0.0
    inc_recall = len(novel) / total_peaks if total_peaks > 0 else 0.0
    swap = len(lost) / total_b if total_b > 0 else 0.0
    union = baseline_tps | method_tps
    jaccard = len(shared) / len(union) if len(union) > 0 else 0.0

    return {
        "shared": len(shared),
        "novel": len(novel),
        "lost": len(lost),
        "total_method": total_m,
        "total_baseline": total_b,
        "total_peaks": total_peaks,
        "novelty_rate": novelty,
        "incremental_recall": inc_recall,
        "swap_rate": swap,
        "jaccard": jaccard,
    }


# ── Plotting helpers ────────────────────────────────────────────────────────

def _mode_color(label: str) -> str:
    """Map mode label to colour via MODE_COLORS (fall back to PALETTE)."""
    # Try exact match
    if label in MODE_COLORS:
        return MODE_COLORS[label]
    # Try common aliases
    aliases = {
        "trust-end": "trust-ends",
        "ted-def": "ted-default",
        "ted-leaf": "ted-leaf-selection",
        "ted-sclip": "density-asymmetric-softclip",
        "ted-t060": "ted-strict-threshold",
        "ted-mcsa": "ted-minmax-norm",
        "ted-2d": "ted-2d-cluster",
    }
    alias = aliases.get(label, "")
    if alias in MODE_COLORS:
        return MODE_COLORS[alias]
    # Cycle palette
    idx = hash(label) % len(PALETTE)
    return PALETTE[idx]


def _plot_single_metric_bar(
    methods: List[str],
    values: List[float],
    end_label: str,
    metric_name: str,
    ylabel: str,
    title_prefix: str,
    output_dir: Path,
    filename: str,
    fmt_pct: bool = True,
    color_by_method: bool = True,
    fixed_color: Optional[str] = None,
) -> None:
    """Create a single bar chart for one metric across methods."""
    fig, ax = plt.subplots(figsize=(W1, W1 * 0.75))

    x = np.arange(len(methods))
    colors = [_mode_color(m) for m in methods] if color_by_method else \
             [fixed_color or PALETTE[0]] * len(methods)
    bars = ax.bar(x, values, width=0.65, color=colors, edgecolor="none")

    # Value labels atop each bar
    for bar_obj, v in zip(bars, values):
        label_text = f"{v:.1%}" if fmt_pct else f"{v:.0f}"
        ax.text(bar_obj.get_x() + bar_obj.get_width() / 2,
                bar_obj.get_height() + 0.005,
                label_text, ha="center", va="bottom", fontsize=6)

    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=45, ha="right")
    if fmt_pct:
        ax.set_ylim(0, max(max(values) * 1.15, 0.05) if values else 0.05)
    else:
        ax.set_ylim(0, max(max(values) * 1.15, 1) if values else 1)

    style_ax(ax, ylabel=ylabel)

    savefig(fig, output_dir / filename)


def _plot_tp_summary_bars(
    methods: List[str],
    overlaps: List[Dict[str, float]],
    end_label: str,
    title_prefix: str,
    output_dir: Path,
    filename: str,
) -> None:
    """Create grouped bar chart showing shared/novel/lost TP counts."""
    fig, ax = plt.subplots(figsize=(W2, W1 * 0.8))

    n = len(methods)
    x = np.arange(n)
    bar_w = 0.25

    shared_vals = [o["shared"] for o in overlaps]
    novel_vals = [o["novel"] for o in overlaps]
    lost_vals = [o["lost"] for o in overlaps]

    ax.bar(x - bar_w, shared_vals, bar_w, label="Shared TPs",
           color=OVERLAP_COLORS["shared"], edgecolor="none")
    ax.bar(x, novel_vals, bar_w, label="Novel TPs",
           color=OVERLAP_COLORS["novel"], edgecolor="none")
    ax.bar(x + bar_w, lost_vals, bar_w, label="Lost TPs",
           color=OVERLAP_COLORS["lost"], edgecolor="none")

    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=45, ha="right")
    style_ax(ax, ylabel="Peak count")
    legend_outside(ax, loc="upper right", bbox_to_anchor=(1.0, 1.0))

    savefig(fig, output_dir / filename)


def _plot_combined_rates(
    methods: List[str],
    overlaps: List[Dict[str, float]],
    end_label: str,
    title_prefix: str,
    output_dir: Path,
    filename: str,
) -> None:
    """Side-by-side novelty rate, incremental recall, swap rate."""
    fig, ax = plt.subplots(figsize=(W2, W1 * 0.8))

    n = len(methods)
    x = np.arange(n)
    bar_w = 0.25

    novelty_vals = [o["novelty_rate"] for o in overlaps]
    inc_recall_vals = [o["incremental_recall"] for o in overlaps]
    swap_vals = [o["swap_rate"] for o in overlaps]

    b1 = ax.bar(x - bar_w, novelty_vals, bar_w, label="Novelty rate",
                color=OVERLAP_COLORS["novel"], edgecolor="none")
    b2 = ax.bar(x, inc_recall_vals, bar_w, label="Incremental recall",
                color=OVERLAP_COLORS["shared"], edgecolor="none")
    b3 = ax.bar(x + bar_w, swap_vals, bar_w, label="Swap rate",
                color=OVERLAP_COLORS["lost"], edgecolor="none")

    # Value labels
    for bars in (b1, b2, b3):
        for bar_obj in bars:
            h = bar_obj.get_height()
            if h > 0.001:
                ax.text(bar_obj.get_x() + bar_obj.get_width() / 2,
                        h + 0.003, f"{h:.1%}", ha="center", va="bottom",
                        fontsize=5.5)

    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=45, ha="right")
    ax.set_ylim(0, max(max(novelty_vals + inc_recall_vals + swap_vals) * 1.2, 0.05))
    style_ax(ax, ylabel="Rate")
    legend_outside(ax, loc="upper right", bbox_to_anchor=(1.0, 1.0))

    savefig(fig, output_dir / filename)


# ── Main analysis ──────────────────────────────────────────────────────────

def analyse_end(
    label_paths: List[str],
    baseline_label: str,
    end_name: str,
    end_tag: str,
    output_dir: Path,
    title_prefix: str,
) -> bool:
    """Run overlap analysis for one end type (5prime/3prime).

    Returns True if plots were saved, False on error.
    """
    # Load data
    data: Dict[str, Dict[str, str]] = {}
    for spec in label_paths:
        label, path = _parse_label_path(spec)
        if not path.exists():
            logger.warning("File not found: %s", path)
            continue
        data[label] = load_peak_reason_tsv(path)

    if baseline_label not in data:
        logger.error("Baseline label '%s' not found in %s data", baseline_label, end_name)
        return False

    # All peaks (same peak set across methods — use baseline's key set)
    all_peaks = set(data[baseline_label].keys())
    total_peaks = len(all_peaks)
    baseline_tps = get_tp_set(data[baseline_label])

    logger.info("%s: %d total peaks, %d baseline TPs",
                end_name, total_peaks, len(baseline_tps))

    # Compute overlaps for each non-baseline method
    methods = []
    overlaps = []
    for label, peak_reasons in data.items():
        if label == baseline_label:
            continue
        method_tps = get_tp_set(peak_reasons)
        ov = compute_overlap(baseline_tps, method_tps, total_peaks)
        methods.append(label)
        overlaps.append(ov)
        logger.info("  %s: shared=%d  novel=%d  lost=%d  "
                     "novelty=%.1f%%  inc_recall=%.1f%%  swap=%.1f%%",
                     label, ov["shared"], ov["novel"], ov["lost"],
                     ov["novelty_rate"] * 100,
                     ov["incremental_recall"] * 100,
                     ov["swap_rate"] * 100)

    if not methods:
        logger.warning("No non-baseline methods for %s — skipping.", end_name)
        return False

    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Novelty rate bars
    _plot_single_metric_bar(
        methods, [o["novelty_rate"] for o in overlaps],
        end_name, "Novelty rate",
        "novel TPs / total TPs(method)", title_prefix,
        output_dir, f"novelty_rate_{end_tag}.png",
    )

    # 2. Incremental recall bars
    _plot_single_metric_bar(
        methods, [o["incremental_recall"] for o in overlaps],
        end_name, "Incremental recall",
        "novel TPs / total peaks", title_prefix,
        output_dir, f"incremental_recall_{end_tag}.png",
    )

    # 3. Swap rate bars
    _plot_single_metric_bar(
        methods, [o["swap_rate"] for o in overlaps],
        end_name, "Swap rate",
        "lost TPs / total TPs(baseline)", title_prefix,
        output_dir, f"swap_rate_{end_tag}.png",
    )

    # 4. Summary grouped bars (shared / novel / lost)
    _plot_tp_summary_bars(
        methods, overlaps, end_name, title_prefix,
        output_dir, f"tp_overlap_summary_{end_tag}.png",
    )

    # 5. Combined rates (novelty + incremental recall + swap side-by-side)
    _plot_combined_rates(
        methods, overlaps, end_name, title_prefix,
        output_dir, f"tp_overlap_rates_{end_tag}.png",
    )

    # 6. Write TSV
    tsv_path = output_dir / f"tp_overlap_metrics_{end_tag}.tsv"
    with open(tsv_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "method", "shared_TPs", "novel_TPs", "lost_TPs",
            "total_TPs_method", "total_TPs_baseline", "total_peaks",
            "novelty_rate", "incremental_recall", "swap_rate", "jaccard",
        ])
        for m, o in zip(methods, overlaps):
            writer.writerow([
                m, o["shared"], o["novel"], o["lost"],
                o["total_method"], o["total_baseline"], o["total_peaks"],
                f"{o['novelty_rate']:.4f}",
                f"{o['incremental_recall']:.4f}",
                f"{o['swap_rate']:.4f}",
                f"{o['jaccard']:.4f}",
            ])
    logger.info("  Wrote %s", tsv_path)

    return True


# ── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--cage", nargs="+", default=[],
        help="CAGE peak-reason TSVs as label:path pairs.")
    parser.add_argument(
        "--drna", nargs="+", default=[],
        help="dRNA peak-reason TSVs as label:path pairs.")
    parser.add_argument(
        "--output", "-o", required=True,
        help="Output directory for plots and tables.")
    parser.add_argument(
        "--baseline-label", default="baseline",
        help="Label of the baseline method (default: 'baseline').")
    parser.add_argument(
        "--title-prefix", default="",
        help="Optional prefix for plot titles.")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
    )

    out = Path(args.output)
    ok = True

    if args.cage:
        if not analyse_end(args.cage, args.baseline_label,
                           "5\u2032 (CAGE)", "5prime", out, args.title_prefix):
            ok = False
    else:
        logger.warning("No --cage inputs; skipping 5' analysis.")

    if args.drna:
        if not analyse_end(args.drna, args.baseline_label,
                           "3\u2032 (dRNA)", "3prime", out, args.title_prefix):
            ok = False
    else:
        logger.warning("No --drna inputs; skipping 3' analysis.")

    if ok:
        logger.info("TP overlap analysis complete → %s", out)
    else:
        logger.error("Some analyses failed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
