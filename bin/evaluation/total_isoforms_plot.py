#!/usr/bin/env python3
"""
Total isoforms horizontal bar chart — one bar per method.

Reads evaluation TSV files (one per method) and plots a horizontal bar
chart of the total isoform count.

Usage:
    python total_isoforms_plot.py \\
        --eval eval1.tsv eval2.tsv ... \\
        --output output_dir/
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pub_style import style_ax, savefig, W1, MODE_COLORS, PALETTE


# ── Helpers ─────────────────────────────────────────────────────────────────

def _mode_color(mode: str) -> str:
    if mode in MODE_COLORS:
        return MODE_COLORS[mode]
    return PALETTE[hash(mode) % len(PALETTE)]


def _load_eval_row(path: str) -> dict | None:
    """Load the first data row from an evaluation TSV."""
    p = Path(path)
    if not p.exists():
        return None
    with open(p) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            return dict(row)
    return None


# ── Plotting ────────────────────────────────────────────────────────────────

def plot_total_isoforms(method_counts: dict, output_dir: str):
    """Horizontal bar chart of total isoform counts per method.

    Parameters
    ----------
    method_counts : dict[str, int]
        Mapping from method label to total isoform count.
    output_dir : str or Path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    methods = list(method_counts.keys())
    counts = [method_counts[m] for m in methods]
    colors = [_mode_color(m) for m in methods]

    fig, ax = plt.subplots(figsize=(W1, W1 * 0.78))
    y = np.arange(len(methods))
    bars = ax.barh(y, counts, 0.6, color=colors, edgecolor="none")
    xmax = max(counts) if counts else 1
    for bar, c in zip(bars, counts):
        ax.text(
            bar.get_width() + xmax * 0.02,
            bar.get_y() + bar.get_height() / 2,
            f"{c:,}",
            va="center", ha="left", fontsize=6, color="#333333",
        )

    ax.set_yticks(y)
    ax.set_yticklabels(methods)
    ax.invert_yaxis()
    style_ax(ax, xlabel="Total isoforms")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, xmax * 1.18)
    ax.set_axisbelow(True)
    fig.tight_layout(pad=0.3)
    savefig(fig, output_dir / "total_isoforms.png")


# ── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--eval", nargs="+", required=True,
        help="Evaluation TSV files (each must contain 'isoforms_observed' and "
             "'transcriptome_mode' columns)",
    )
    parser.add_argument(
        "--output", required=True,
        help="Output directory for the bar chart",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    method_counts: dict[str, int] = {}
    for tsv_path in args.eval:
        row = _load_eval_row(tsv_path)
        if row is None:
            print(f"WARNING: skipping missing file {tsv_path}", file=sys.stderr)
            continue
        mode = row.get("transcriptome_mode", Path(tsv_path).stem)
        count = int(float(row.get("isoforms_observed", 0)))
        method_counts[mode] = count
        if args.verbose:
            print(f"  {mode}: {count:,} isoforms")

    if not method_counts:
        print("No data loaded — skipping total isoforms plot", file=sys.stderr)
        sys.exit(1)

    plot_total_isoforms(method_counts, args.output)
    print(f"Saved total isoforms bar chart to {args.output}")


if __name__ == "__main__":
    main()
