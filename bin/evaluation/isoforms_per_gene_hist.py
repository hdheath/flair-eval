#!/usr/bin/env python3
"""
Isoforms-per-gene distribution plots — histogram and/or box plot.

Reads BED12 isoform files (one per method), extracts gene IDs from the name
field (ENST…_ENSG…), counts isoforms per gene, and plots:
  - Overlaid step histograms (default)
  - Side-by-side box plots (--box flag)

Usage:
    python isoforms_per_gene_hist.py \\
        --bed label1:bed1.bed label2:bed2.bed ... \\
        --output output_dir/

    python isoforms_per_gene_hist.py --box \\
        --bed label1:bed1.bed label2:bed2.bed ... \\
        --output output_dir/
"""

import argparse
import sys
from collections import Counter
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pub_style import style_ax, savefig, W1, PALETTE, MODE_COLORS
from signal_utils import parse_isoforms, gene_from_name


# ── Helpers ─────────────────────────────────────────────────────────────────

def _mode_color(mode: str) -> str:
    if mode in MODE_COLORS:
        return MODE_COLORS[mode]
    return PALETTE[hash(mode) % len(PALETTE)]


def _count_isoforms_per_gene(bed_path: str) -> Counter:
    """Parse BED12 or GTF, extract gene ID, return isoforms-per-gene counts."""
    gene_counts: Counter = Counter()
    isoforms = parse_isoforms(bed_path)
    for iso in isoforms:
        gid = gene_from_name(iso["name"])
        gene_counts[gid] += 1
    return gene_counts


# ── Plotting ────────────────────────────────────────────────────────────────

def plot_isoforms_per_gene_hist(
    method_gene_counts: dict,
    output_dir: str,
):
    """Create a multi-panel bar chart of isoforms-per-gene, one panel per method.

    Parameters
    ----------
    method_gene_counts : dict[str, Counter]
        Mapping from method label to Counter of {gene_id: n_isoforms}.
    output_dir : str or Path
        Directory to write the histogram PNG/SVG.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    methods = list(method_gene_counts.keys())
    n_methods = len(methods)
    if n_methods == 0:
        return

    # Determine common bin range: 1, 2, 3, ..., max capped at 10+
    max_val = 1
    for counts in method_gene_counts.values():
        if counts:
            max_val = max(max_val, max(counts.values()))
    bin_max = min(max_val + 1, 11)  # cap at 10+
    bin_vals = list(range(1, bin_max))

    # Compute proportions for each method
    method_props = {}
    for method, counts in method_gene_counts.items():
        if not counts:
            method_props[method] = [0.0] * len(bin_vals)
            continue
        total = len(counts)  # number of genes, not isoforms
        props = []
        for b in bin_vals:
            if b < bin_max - 1:
                n = sum(1 for v in counts.values() if v == b)
            else:
                n = sum(1 for v in counts.values() if v >= b)
            props.append(n / total if total > 0 else 0)
        method_props[method] = props

    # Layout: one row per method
    fig_h = max(1.2 * n_methods, 2.5)
    fig, axes = plt.subplots(n_methods, 1, figsize=(W1, fig_h), sharex=True,
                             squeeze=False)

    tick_labels = [str(v) if v < bin_max - 1 else f"{v}+" for v in bin_vals]
    x = np.arange(len(bin_vals))

    # Compute shared y-axis max across all methods
    y_max = 0.0
    for props in method_props.values():
        if props:
            y_max = max(y_max, max(props))
    y_max = y_max * 1.1  # add 10% headroom

    for idx, method in enumerate(methods):
        ax = axes[idx, 0]
        color = PALETTE[idx % len(PALETTE)]
        ax.bar(x, method_props[method], width=0.7, color=color, edgecolor="none",
               alpha=0.88)
        style_ax(ax, ylabel="Proportion", faint_y_grid=True)
        ax.set_title(method, fontsize=7, pad=3, loc="left", fontweight="semibold")
        ax.set_ylim(0, y_max)
        if idx == n_methods - 1:
            ax.set_xticks(x)
            ax.set_xticklabels(tick_labels, fontsize=7)
            ax.set_xlabel("Isoforms per gene", fontsize=7)
        else:
            ax.tick_params(labelbottom=False)

    fig.tight_layout(h_pad=0.6)
    savefig(fig, output_dir / "isoforms_per_gene_histogram.png", dpi=300)


def plot_total_genes_bar(
    method_gene_counts: dict,
    output_dir: str,
):
    """Horizontal bar chart of total genes detected per method."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    methods = list(method_gene_counts.keys())
    totals = [len(method_gene_counts[m]) for m in methods]

    fig, ax = plt.subplots(figsize=(W1, max(0.45 * len(methods), 1.5)))
    y = np.arange(len(methods))
    colors = [PALETTE[i % len(PALETTE)] for i in range(len(methods))]

    ax.barh(y, totals, height=0.6, color=colors, edgecolor="none", alpha=0.88)

    for i, v in enumerate(totals):
        ax.text(v + max(totals) * 0.02, i, str(v), va="center", fontsize=6.5)

    ax.set_yticks(y)
    ax.set_yticklabels(methods, fontsize=7)
    ax.invert_yaxis()
    style_ax(ax, xlabel="Total genes", faint_y_grid=False)
    ax.set_xlim(0, max(totals) * 1.15)

    fig.tight_layout()
    savefig(fig, output_dir / "total_genes_bar.png", dpi=300)


def plot_isoforms_per_gene_box(
    method_gene_counts: dict,
    output_dir: str,
):
    """Create a box plot of isoforms-per-gene across methods.

    Parameters
    ----------
    method_gene_counts : dict[str, Counter]
        Mapping from method label to Counter of {gene_id: n_isoforms}.
    output_dir : str or Path
        Directory to write the box plot PNG/SVG.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    methods = list(method_gene_counts.keys())
    all_counts = [list(method_gene_counts[m].values()) for m in methods]

    fig, ax = plt.subplots(figsize=(W1, W1 * 0.82))
    bp = ax.boxplot(
        all_counts,
        patch_artist=True,
        showfliers=False,
        widths=0.55,
        medianprops=dict(color="black", linewidth=0.6),
        whiskerprops=dict(linewidth=0.35),
        capprops=dict(linewidth=0.35),
    )
    for patch, m in zip(bp["boxes"], methods):
        patch.set_facecolor(_mode_color(m))
        patch.set_alpha(0.85)
        patch.set_edgecolor("none")

    # Diamond marker at the mean
    for i, counts in enumerate(all_counts):
        if counts:
            ax.scatter(
                i + 1, np.mean(counts),
                color="black", marker="D", s=10, zorder=5, linewidths=0,
            )

    ax.set_xticks(range(1, len(methods) + 1))
    ax.set_xticklabels(methods, rotation=40, ha="right")
    style_ax(ax, ylabel="Isoforms per gene")
    ax.set_axisbelow(True)
    fig.tight_layout(pad=0.3)
    savefig(fig, output_dir / "isoforms_per_gene_boxplot.png")


# ── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--bed", nargs="+", required=True,
        help="label:path pairs, e.g. baseline:path/to/baseline.isoforms.bed",
    )
    parser.add_argument(
        "--output", required=True,
        help="Output directory for the histogram plot",
    )
    parser.add_argument(
        "--box", action="store_true",
        help="Also generate a box plot (in addition to the histogram)",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    method_gene_counts = {}
    for entry in args.bed:
        if ":" not in entry:
            print(f"WARNING: skipping malformed entry '{entry}' (expected label:path)",
                  file=sys.stderr)
            continue
        label, path = entry.split(":", 1)
        if not Path(path).exists():
            print(f"WARNING: file not found: {path}", file=sys.stderr)
            continue
        counts = _count_isoforms_per_gene(path)
        method_gene_counts[label] = counts
        if args.verbose:
            print(f"  {label}: {sum(counts.values())} isoforms across "
                  f"{len(counts)} genes")

    if not method_gene_counts:
        print("No data loaded — skipping histogram", file=sys.stderr)
        sys.exit(1)

    plot_isoforms_per_gene_hist(method_gene_counts, args.output)
    plot_total_genes_bar(method_gene_counts, args.output)

    if args.box:
        plot_isoforms_per_gene_box(method_gene_counts, args.output)
        print(f"Saved isoforms-per-gene histogram + box plot to {args.output}")
    else:
        print(f"Saved isoforms-per-gene histogram to {args.output}")


if __name__ == "__main__":
    main()
