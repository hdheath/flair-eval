#!/usr/bin/env python3
"""Gene-level isoform variation proportions plot.

For each BED12 isoform file, classifies genes into:
  - Single isoform
  - Alt ends only (same SJC, different TSS/TTS)
  - Alt splicing only (different SJCs, each with one end combo)
  - Alt splicing + ends (different SJCs and alt ends within at least one)

Produces a stacked proportional bar chart comparing samples.
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from signal_utils import parse_isoforms
    from end_variation import classify_genes_by_variation
except ImportError:
    from evaluation.signal_utils import parse_isoforms
    from evaluation.end_variation import classify_genes_by_variation

GOLDEN_RATIO = 1.618
CATEGORIES = ["single_isoform", "alt_ends_only", "alt_splicing_only", "alt_splicing_and_ends"]
LABELS = ["Single isoform", "Alt ends only", "Alt splicing only", "Alt splicing + ends"]
COLORS = ["#999999", "#4C78A8", "#F58518", "#E45756"]


def plot(gene_counts_by_sample, output_path, title="Gene-Level Isoform Variation"):
    samples = list(gene_counts_by_sample.keys())
    n = len(samples)

    proportions = {cat: [] for cat in CATEGORIES}
    raw_counts = {cat: [] for cat in CATEGORIES}
    for s in samples:
        counts = gene_counts_by_sample[s]
        total = sum(counts.get(c, 0) for c in CATEGORIES)
        for cat in CATEGORIES:
            c = counts.get(cat, 0)
            raw_counts[cat].append(c)
            proportions[cat].append(c / total if total > 0 else 0)

    fig_width = max(5.0, 1.6 * n + 1.5)
    fig, ax = plt.subplots(figsize=(fig_width, fig_width / GOLDEN_RATIO))

    x = range(n)
    bottoms = [0.0] * n
    for cat, label, color in zip(CATEGORIES, LABELS, COLORS):
        vals = proportions[cat]
        bars = ax.bar(x, vals, bottom=bottoms, label=label, color=color,
                      edgecolor='white', linewidth=0.5, width=0.7)
        for i, (b, v, rc) in enumerate(zip(bars, vals, raw_counts[cat])):
            if v > 0.04:
                ax.text(b.get_x() + b.get_width() / 2,
                        bottoms[i] + v / 2,
                        str(rc), ha='center', va='center', fontsize=7,
                        color='white' if color in ("#999999", "#4C78A8", "#E45756") else 'black')
        bottoms = [b + v for b, v in zip(bottoms, vals)]

    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Proportion of genes", fontsize=8)
    ax.set_xticks(list(x))
    ax.set_xticklabels(samples, fontsize=7, rotation=30, ha='right')
    ax.legend(fontsize=7, loc='upper right', framealpha=0.9)
    ax.set_title(title, fontsize=8, fontweight='normal')
    ax.grid(True, alpha=0.25, linestyle='--', axis='y')
    ax.set_axisbelow(True)

    totals = [sum(raw_counts[cat][i] for cat in CATEGORIES) for i in range(n)]
    footer_parts = [f"{s}: {t} genes" for s, t in zip(samples, totals)]
    fig.text(0.5, 0.01, " | ".join(footer_parts), ha='center', fontsize=6)

    plt.tight_layout(rect=(0, 0.04, 1, 1.0))
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"[INFO] Saved plot to {output_path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bed", nargs="+", required=True,
        help="label:path pairs for BED12 isoform files",
    )
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--title", default="Gene-Level Isoform Variation")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    gene_counts = {}
    for item in args.bed:
        if ":" not in item:
            print(f"[WARN] Skipping malformed bed arg (no ':'): {item}", file=sys.stderr)
            continue
        label, bed_path = item.split(":", 1)
        p = Path(bed_path)
        if not p.exists():
            print(f"[WARN] BED file not found: {p}", file=sys.stderr)
            continue
        isos = parse_isoforms(p)
        counts = classify_genes_by_variation(isos)
        gene_counts[label] = counts
        if args.verbose:
            print(f"  {label}: {counts}", file=sys.stderr)

    if not gene_counts:
        print("[WARN] No valid samples, skipping plot.", file=sys.stderr)
        sys.exit(0)

    plot(gene_counts, out_dir / "gene_variation_proportions.png", title=args.title)


if __name__ == "__main__":
    main()
