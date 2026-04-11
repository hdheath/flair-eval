#!/usr/bin/env python3
"""
Signal-vs-support scatter dashboard: tiled per-mode scatter plots.

Reads per-mode peak-reason TSVs (each containing per-peak records with
peak_id, score, read_count, reason) and tiles compact scatter plots into
an N×2 grid (one column per end type: CAGE / QuantSeq).  Each sub-plot
reproduces the per-run signal-vs-support diagnostic at thumbnail scale
for rapid cross-mode visual comparison.

The mode name is extracted from each filename using the standard naming
convention: ``{dataset}_{alignment}_{partition}_{mode}_transcriptome_{endtype}_peak_reasons.tsv``

Usage:
    python signal_support_dashboard.py \\
        --cage-tsvs  mode1_cage_peak_reasons.tsv mode2_cage_peak_reasons.tsv ... \\
        --quantseq-tsvs mode1_quantseq_peak_reasons.tsv mode2_quantseq_peak_reasons.tsv ... \\
        --output signal_support_dashboard.png \\
        [--title-prefix "test_name: "] [--verbose]
"""

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from pub_style import style_ax, savefig, REASON_COLORS as _PUB_REASON_COLORS, PALETTE


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

REASON_COLORS = _PUB_REASON_COLORS

REASON_LABELS = {
    "recovered":           "Recovered",
    "no_reads":            "No long reads",
    "single_exon_only":    "Single-exon only",
    "near_miss":           "Near miss",
    "reads_unassigned":    "Reads unassigned",
    "reads_redirected":    "Reads redirected",
    "proximal_apa":        "Proximal APA",
    "trailing_truncation": "Trailing / truncation",
    "other_missed":        "Other missed",
}

REASON_ORDER = [
    "recovered", "near_miss", "reads_redirected", "proximal_apa",
    "trailing_truncation", "single_exon_only", "reads_unassigned",
    "no_reads", "other_missed",
]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _extract_mode_from_filename(filename: str) -> str:
    """Extract the transcriptome_mode from a peak_reasons filename.

    Expected pattern:
        {dataset}_{align}_{partition}_{mode}_transcriptome_{endtype}_peak_reasons.tsv

    Falls back to the full stem if pattern doesn't match.
    """
    stem = Path(filename).stem  # remove .tsv
    # Try to match the standard naming convention
    m = re.search(r"_([^_]+)_transcriptome_(?:cage|quantseq)_peak_reasons$", stem)
    if m:
        return m.group(1)
    # Fallback: try broader pattern
    m = re.search(r"chr\d+_([^_]+(?:_[^_]+)*)_transcriptome_", stem)
    if m:
        return m.group(1)
    return stem


def load_peak_reason_tsvs(paths):
    """Load peak-reason TSVs, returning dict[mode] -> list of record dicts."""
    mode_records = {}
    for p in paths:
        p = Path(p)
        mode = _extract_mode_from_filename(p.name)
        records = []
        try:
            with open(p) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    try:
                        score = float(row.get("score", 0))
                    except (ValueError, TypeError):
                        score = 0.0
                    try:
                        read_count = int(float(row.get("read_count", 0)))
                    except (ValueError, TypeError):
                        read_count = 0
                    reason = row.get("reason", "other_missed")
                    records.append({
                        "peak_id": row.get("peak_id", ""),
                        "score": score,
                        "read_count": read_count,
                        "reason": reason,
                    })
        except Exception as e:
            print(f"Warning: could not read {p}: {e}", file=sys.stderr)
            continue
        if records:
            mode_records[mode] = records
    return mode_records


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _plot_scatter_panel(ax, records, title, compact=True):
    """Single signal-vs-support scatter plot for one mode × one end type.

    Args:
        ax: matplotlib Axes
        records: list of {peak_id, score, read_count, reason} dicts
        title: panel title (mode name)
        compact: if True, use smaller markers and minimal annotation
    """
    if not records:
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="#999999")
        ax.set_visible(True)
        return

    # Group by reason
    by_reason = {}
    for r in records:
        by_reason.setdefault(r["reason"], []).append(r)

    # Plot in reverse order so recovered is on top
    for reason in reversed(REASON_ORDER):
        if reason not in by_reason:
            continue
        pts = by_reason[reason]
        scores = [max(p["score"], 0.01) for p in pts]
        counts = [max(p["read_count"], 0.3) for p in pts]
        color = REASON_COLORS.get(reason, "#95a5a6")

        alpha = 0.80 if reason == "recovered" else 0.55
        size = 12 if compact else 25
        marker = "o" if reason == "recovered" else "s"
        zorder = 10 if reason == "recovered" else 5

        ax.scatter(
            scores, counts,
            c=color, s=size, alpha=alpha, marker=marker,
            edgecolors="none", linewidths=0, zorder=zorder,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")

    # Reference lines
    ax.axhline(y=3, color="gray", linestyle=":", linewidth=0.6, alpha=0.3)
    ax.axvline(x=1.0, color="gray", linestyle=":", linewidth=0.6, alpha=0.3)

    # Summary annotation
    n_total = len(records)
    n_recovered = sum(1 for r in records if r["reason"] == "recovered")
    pct = (n_recovered / n_total * 100) if n_total > 0 else 0

    ax.text(0.02, 0.98, f"{n_recovered}/{n_total} ({pct:.0f}%)",
            transform=ax.transAxes, fontsize=7, va="top", ha="left",
            color="#333333", fontweight="normal",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor="none"))

    ax.grid(True, alpha=0.15, which="both")
    ax.set_axisbelow(True)

    style_ax(ax, title=title)


def create_dashboard(cage_modes, quantseq_modes, output_path, title_prefix=""):
    """Create the tiled signal-vs-support scatter dashboard.

    Args:
        cage_modes: dict[mode] -> list of records (5' / CAGE)
        quantseq_modes: dict[mode] -> list of records (3' / QuantSeq)
        output_path: output PNG path
        title_prefix: optional prefix
    """
    # Determine all modes and sort
    all_modes = sorted(set(list(cage_modes.keys()) + list(quantseq_modes.keys())))
    if not all_modes:
        print("No peak-reason data to plot.", file=sys.stderr)
        return False

    has_cage = bool(cage_modes)
    has_quantseq = bool(quantseq_modes)
    n_cols = (1 if has_cage else 0) + (1 if has_quantseq else 0)
    if n_cols == 0:
        print("No data columns available.", file=sys.stderr)
        return False

    n_rows = len(all_modes)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(3.5 * n_cols, 2.5 * n_rows),
        squeeze=False,
    )

    col_idx = 0

    if has_cage:
        for row_idx, mode in enumerate(all_modes):
            ax = axes[row_idx, col_idx]
            records = cage_modes.get(mode, [])
            title = f"{mode}" if row_idx == 0 else mode
            _plot_scatter_panel(ax, records, title)
            if row_idx == 0:
                ax.set_title(f"CAGE (5\u2032) — {mode}", fontsize=7)
            else:
                ax.set_title(mode, fontsize=7)
            if row_idx == n_rows - 1:
                ax.set_xlabel("Peak Signal (TPM)", fontsize=7)
            ax.set_ylabel("Read Count", fontsize=8)
        col_idx += 1

    if has_quantseq:
        for row_idx, mode in enumerate(all_modes):
            ax = axes[row_idx, col_idx]
            records = quantseq_modes.get(mode, [])
            title = f"{mode}" if row_idx == 0 else mode
            _plot_scatter_panel(ax, records, title)
            if row_idx == 0:
                ax.set_title(f"QuantSeq (3\u2032) — {mode}", fontsize=7)
            else:
                ax.set_title(mode, fontsize=7)
            if row_idx == n_rows - 1:
                ax.set_xlabel("Peak Signal (TPM)", fontsize=7)
            if col_idx == 0:
                ax.set_ylabel("Read Count", fontsize=8)

    # Shared reason legend at bottom
    all_reasons_present = set()
    for records_list in list(cage_modes.values()) + list(quantseq_modes.values()):
        for r in records_list:
            all_reasons_present.add(r["reason"])
    legend_reasons = [r for r in REASON_ORDER if r in all_reasons_present]
    legend_patches = [
        Patch(facecolor=REASON_COLORS.get(r, "#95a5a6"), edgecolor="none",
              label=REASON_LABELS.get(r, r))
        for r in legend_reasons
    ]
    fig.legend(
        handles=legend_patches, loc="lower center",
        ncol=min(5, len(legend_patches)), fontsize=8,
        frameon=False,
    )

    fig.tight_layout(rect=(0, 0.05, 1, 1.0), h_pad=1.5)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    savefig(fig, output_path, dpi=300)
    print(f"Saved signal-vs-support dashboard to {output_path}")
    return True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Tiled signal-vs-support scatter dashboard across assembler modes",
    )
    parser.add_argument("--cage-tsvs", nargs="*", default=[],
                        help="Per-mode CAGE peak-reason TSVs")
    parser.add_argument("--quantseq-tsvs", nargs="*", default=[],
                        help="Per-mode QuantSeq peak-reason TSVs")
    parser.add_argument("--output", required=True, help="Output PNG path")
    parser.add_argument("--title-prefix", default="", help="Title prefix")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if not args.cage_tsvs and not args.quantseq_tsvs:
        print("No peak-reason TSVs provided — nothing to plot.", file=sys.stderr)
        sys.exit(0)

    cage_modes = {}
    quantseq_modes = {}

    if args.cage_tsvs:
        cage_modes = load_peak_reason_tsvs(args.cage_tsvs)
        if args.verbose:
            print(f"Loaded CAGE peak reasons for {len(cage_modes)} modes: {list(cage_modes.keys())}")

    if args.quantseq_tsvs:
        quantseq_modes = load_peak_reason_tsvs(args.quantseq_tsvs)
        if args.verbose:
            print(f"Loaded QuantSeq peak reasons for {len(quantseq_modes)} modes: {list(quantseq_modes.keys())}")

    if not cage_modes and not quantseq_modes:
        print("No valid peak-reason data found — skipping dashboard.", file=sys.stderr)
        sys.exit(0)

    success = create_dashboard(cage_modes, quantseq_modes, args.output, args.title_prefix)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
