#!/usr/bin/env python3
"""
jc_end_distance_plot.py — Within-junction-chain end distance distributions.

For each junction chain group with >1 isoform, computes all pairwise
distances between isoform 5' ends (TSS) and 3' ends (TTS).  Produces
two panel-grid figures (one for 5', one for 3') where:
  - Each row is a transcriptome mode (default, ted-2d, ted-1d2d)
  - Each column is a sample
  - Histograms show the distribution of pairwise distances

Usage:
  python jc_end_distance_plot.py \
      --beds mode1:sample1:/path/to/isoforms.bed mode1:sample2:/path/to/isoforms.bed ... \
      --outdir /path/to/output
"""

import argparse
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_bed12(bed_path):
    """Parse BED12 into list of dicts with junction chains and ends."""
    isoforms = []
    with open(bed_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            strand = parts[5]
            esizes = [int(x) for x in parts[10].rstrip(",").split(",") if x]
            estarts = [int(x) for x in parts[11].rstrip(",").split(",") if x]
            exons = [(start + estarts[i], start + estarts[i] + esizes[i])
                     for i in range(len(esizes))]
            introns = tuple((exons[x][1], exons[x + 1][0])
                            for x in range(len(exons) - 1))
            if not introns:
                continue  # skip single-exon
            if strand == "+":
                tss, tts = start, end
            else:
                tss, tts = end, start
            isoforms.append({
                "chrom": chrom, "strand": strand,
                "introns": introns, "tss": tss, "tts": tts,
            })
    return isoforms


def compute_pairwise_distances(isoforms):
    """Group by JC, compute pairwise |TSS-TSS| and |TTS-TTS| distances."""
    jc_groups = defaultdict(list)
    for iso in isoforms:
        key = (iso["chrom"], iso["strand"], iso["introns"])
        jc_groups[key].append(iso)

    tss_dists = []
    tts_dists = []
    for key, members in jc_groups.items():
        if len(members) < 2:
            continue
        for a, b in combinations(members, 2):
            tss_dists.append(abs(a["tss"] - b["tss"]))
            tts_dists.append(abs(a["tts"] - b["tts"]))
    return tss_dists, tts_dists


def main():
    parser = argparse.ArgumentParser(
        description="Plot within-JC end distance distributions across modes.")
    parser.add_argument("--beds", nargs="+", required=True,
                        help="mode:sample:/path/to/isoforms.bed triplets")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--max-dist", type=int, default=2000,
                        help="Max distance to show in histogram (default: 2000)")
    parser.add_argument("--bins", type=int, default=50)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Parse inputs: {(mode, sample): bed_path}
    entries = []
    for spec in args.beds:
        parts = spec.split(":", 2)
        if len(parts) != 3:
            print(f"ERROR: expected mode:sample:path, got {spec}", file=sys.stderr)
            sys.exit(1)
        entries.append({"mode": parts[0], "sample": parts[1], "path": parts[2]})

    # Determine grid: rows=modes, cols=samples (in order of appearance)
    mode_order = list(dict.fromkeys(e["mode"] for e in entries))
    sample_order = list(dict.fromkeys(e["sample"] for e in entries))

    # Compute distances
    data = {}
    for e in entries:
        isoforms = parse_bed12(e["path"])
        tss_d, tts_d = compute_pairwise_distances(isoforms)
        data[(e["mode"], e["sample"])] = {
            "tss": tss_d, "tts": tts_d, "n_isoforms": len(isoforms),
        }

    # Colors per mode
    mode_colors = {
        "default": "#636363",
        "ted-2d": "#3498db",
        "ted-1d2d": "#e74c3c",
    }
    # Fallback for unknown modes
    fallback_colors = ["#2ecc71", "#f39c12", "#9b59b6", "#1abc9c"]

    n_rows = len(mode_order)
    n_cols = len(sample_order)
    bin_edges = np.linspace(0, args.max_dist, args.bins + 1)

    for end_type, end_label, end_title in [("tss", "5prime", "5' (TSS)"),
                                            ("tts", "3prime", "3' (TTS)")]:
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 2.5 * n_rows),
                                 squeeze=False, sharex=True, sharey=False)
        fig.suptitle(f"Within-JC Pairwise End Distances — {end_title}", fontsize=14, y=1.02)

        for ri, mode in enumerate(mode_order):
            for ci, sample in enumerate(sample_order):
                ax = axes[ri][ci]
                key = (mode, sample)
                if key not in data:
                    ax.set_visible(False)
                    continue

                dists = data[key][end_type]
                clipped = [d for d in dists if d <= args.max_dist]
                color = mode_colors.get(mode, fallback_colors[ri % len(fallback_colors)])

                ax.hist(clipped, bins=bin_edges, color=color, alpha=0.8, edgecolor="white", linewidth=0.3)

                # Stats annotation
                if dists:
                    med = np.median(dists)
                    n_pairs = len(dists)
                    ax.axvline(med, color="black", linewidth=1, linestyle="--", alpha=0.7)
                    ax.text(0.97, 0.95, f"n={n_pairs}\nmed={med:.0f}bp",
                            transform=ax.transAxes, fontsize=7, ha="right", va="top",
                            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))
                else:
                    ax.text(0.5, 0.5, "No multi-iso JCs", transform=ax.transAxes,
                            ha="center", va="center", fontsize=8, color="gray")

                # Labels
                if ri == 0:
                    ax.set_title(sample, fontsize=10)
                if ci == 0:
                    ax.set_ylabel(mode, fontsize=10, fontweight="bold")
                if ri == n_rows - 1:
                    ax.set_xlabel("Pairwise distance (bp)")

                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)

        plt.tight_layout()
        out_path = outdir / f"jc_end_distances_{end_label}.png"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
