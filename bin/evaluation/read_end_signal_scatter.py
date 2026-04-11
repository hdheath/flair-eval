#!/usr/bin/env python3
"""
Read end-signal density scatter — per-sample KDE-coloured panels.

For each sample, plots TSS signal (CAGE) vs TTS signal (QuantSeq) for every
aligned *read* (not assembled isoform), coloured by KDE density.  Shows what
the raw data looks like against orthogonal signal before any assembly.

Usage:
    python read_end_signal_scatter.py \
        --bed label1:reads1.bed label2:reads2.bed ... \
        --cage-plus cage_plus.bg --cage-minus cage_minus.bg \
        --qs-plus qs_plus.bg --qs-minus qs_minus.bg \
        --output output_dir/
"""

import argparse
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from pub_style import savefig, style_ax, W2
from signal_utils import (
    parse_bed12,
    load_signal_tracks,
    isoform_signal,
    kde_density,
)

MAX_READS = 100_000
RNG_SEED = 42


def subsample(reads, max_n=MAX_READS):
    """Deterministic random subsample if len(reads) > max_n."""
    if len(reads) <= max_n:
        return reads
    rng = np.random.default_rng(RNG_SEED)
    idx = rng.choice(len(reads), size=max_n, replace=False)
    return [reads[i] for i in idx]


def plot_read_end_signal_scatter(
    reads_by_sample: dict,
    cage_p, cage_m, qs_p, qs_m,
    output_dir: str,
):
    """Multi-panel KDE-coloured density scatter of TSS vs TTS signal per read."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    samples = list(reads_by_sample.keys())
    n = len(samples)
    if n == 0:
        return

    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    pw = W2 / ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(W2, pw * nrows), squeeze=False)

    eps = 1e-3

    for idx, sample in enumerate(samples):
        ax = axes[idx // ncols][idx % ncols]
        reads = reads_by_sample[sample]
        total = len(reads)
        reads = subsample(reads)

        sigs = [isoform_signal(r, cage_p, cage_m, qs_p, qs_m) for r in reads]
        tss = np.array([s[0] for s in sigs]) + eps
        tts = np.array([s[1] for s in sigs]) + eps

        density = kde_density(tts, tss)
        order = np.argsort(density)

        ax.scatter(
            tts[order], tss[order],
            c=density[order], cmap="viridis_r", vmin=0, vmax=1,
            s=1.5, alpha=0.5, edgecolors="none", rasterized=True,
        )

        ax.set_xscale("log")
        ax.set_yscale("log")

        if idx >= n - ncols:
            ax.set_xlabel("TTS signal (QuantSeq)", fontsize=7)
        else:
            ax.set_xticklabels([])
        if idx % ncols == 0:
            ax.set_ylabel("TSS signal (CAGE)", fontsize=7)
        else:
            ax.set_yticklabels([])

        style_ax(ax)
        ax.text(0.04, 0.96, sample, transform=ax.transAxes,
                ha="left", va="top", fontsize=6, fontweight="bold")
        sub_note = f" ({MAX_READS:,} sampled)" if total > MAX_READS else ""
        ax.text(0.96, 0.04, f"n = {total:,}{sub_note}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=5, color="#666666")

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.subplots_adjust(left=0.07, right=0.90, bottom=0.10, top=0.97,
                        hspace=0.35, wspace=0.30)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.012, 0.7])
    norm = mcolors.Normalize(vmin=0, vmax=1)
    sm = plt.cm.ScalarMappable(cmap="viridis_r", norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("Relative density", fontsize=6)
    cbar.ax.tick_params(labelsize=5, length=1.5, width=0.3)

    savefig(fig, output_dir / "read_end_signal_scatter.png")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--bed", nargs="+", required=True,
        help="label:path pairs for reads BED12 files (one per sample)",
    )
    parser.add_argument("--cage-plus",  required=True, help="CAGE bedGraph (+ strand)")
    parser.add_argument("--cage-minus", required=True, help="CAGE bedGraph (- strand)")
    parser.add_argument("--qs-plus",    required=True, help="QuantSeq bedGraph (+ strand)")
    parser.add_argument("--qs-minus",   required=True, help="QuantSeq bedGraph (- strand)")
    parser.add_argument("--output",     required=True, help="Output directory")
    parser.add_argument("--verbose",    action="store_true")
    args = parser.parse_args()

    reads_by_sample = {}
    for entry in args.bed:
        if ":" not in entry:
            print(f"WARNING: skipping malformed entry '{entry}'", file=sys.stderr)
            continue
        label, path = entry.rsplit(":", 1)
        if not Path(path).exists():
            print(f"WARNING: file not found: {path}", file=sys.stderr)
            continue
        reads = parse_bed12(path)
        if reads:
            reads_by_sample[label] = reads
            if args.verbose:
                print(f"  {label}: {len(reads)} reads")

    if not reads_by_sample:
        print("No read data loaded — skipping", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print("  Loading signal tracks...")
    cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
        args.cage_plus, args.cage_minus, args.qs_plus, args.qs_minus,
    )

    plot_read_end_signal_scatter(reads_by_sample, cage_p, cage_m, qs_p, qs_m, args.output)
    print(f"Saved read end-signal scatter to {args.output}")


if __name__ == "__main__":
    main()
