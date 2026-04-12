#!/usr/bin/env python3
"""
Isoform end-signal density scatter — per-method KDE-coloured panels.

For each method, plots TSS signal (CAGE) vs TTS signal (dRNA) for every
isoform, coloured by KDE density.  Requires BED12 isoform files plus four
bedGraph signal tracks (CAGE +/- strand, dRNA +/- strand).

Usage:
    python end_signal_scatter_plot.py \\
        --bed label1:bed1.bed label2:bed2.bed ... \\
        --cage-plus cage_plus.bg --cage-minus cage_minus.bg \\
        --qs-plus qs_plus.bg --qs-minus qs_minus.bg \\
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
    parse_isoforms,
    load_signal_tracks,
    isoform_signal,
    kde_density,
)


# ── Plotting ────────────────────────────────────────────────────────────────

def plot_end_signal_scatter(
    beds_by_method: dict,
    cage_p, cage_m, qs_p, qs_m,
    output_dir: str,
):
    """Multi-panel KDE-coloured density scatter of TSS vs TTS signal.

    Parameters
    ----------
    beds_by_method : dict[str, list[dict]]
        Mapping from method label to list of BED12 isoform dicts.
    cage_p, cage_m, qs_p, qs_m : BedGraphTrack
        Signal tracks for CAGE +/- and dRNA +/-.
    output_dir : str or Path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    methods = list(beds_by_method.keys())
    n = len(methods)
    if n == 0:
        return

    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    pw = W2 / ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(W2, pw * nrows), squeeze=False)

    eps = 1e-3

    for idx, m in enumerate(methods):
        ax = axes[idx // ncols][idx % ncols]
        isoforms = beds_by_method[m]

        sigs = [isoform_signal(iso, cage_p, cage_m, qs_p, qs_m) for iso in isoforms]
        tss = np.array([s[0] for s in sigs]) + eps
        tts = np.array([s[1] for s in sigs]) + eps

        # KDE density colouring
        density = kde_density(tts, tss)
        order = np.argsort(density)

        ax.scatter(
            tts[order], tss[order],
            c=density[order], cmap="viridis_r", vmin=0, vmax=1,
            s=3, alpha=0.7, edgecolors="none", rasterized=True,
        )

        ax.set_xscale("log")
        ax.set_yscale("log")

        # Axis labels on edges only
        if idx >= n - ncols:
            ax.set_xlabel("TTS signal", fontsize=7)
        else:
            ax.set_xticklabels([])
        if idx % ncols == 0:
            ax.set_ylabel("TSS signal", fontsize=7)
        else:
            ax.set_yticklabels([])

        style_ax(ax)
        ax.text(0.04, 0.96, m, transform=ax.transAxes,
                ha="left", va="top", fontsize=6, fontweight="bold")
        ax.text(0.96, 0.04, f"n = {len(sigs):,}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=5, color="#666666")

    # Hide unused panels
    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    # Layout + shared colourbar
    fig.subplots_adjust(left=0.07, right=0.90, bottom=0.10, top=0.97,
                        hspace=0.35, wspace=0.30)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.012, 0.7])
    norm = mcolors.Normalize(vmin=0, vmax=1)
    sm = plt.cm.ScalarMappable(cmap="viridis_r", norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("Relative density", fontsize=6)
    cbar.ax.tick_params(labelsize=5, length=1.5, width=0.3)

    savefig(fig, output_dir / "end_signal_scatter.png")


# ── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--bed", nargs="+", required=True,
        help="label:path pairs for BED12 isoform files",
    )
    parser.add_argument("--cage-plus",  required=True, help="CAGE bedGraph (+ strand)")
    parser.add_argument("--cage-minus", required=True, help="CAGE bedGraph (- strand)")
    parser.add_argument("--qs-plus",    required=True, help="dRNA bedGraph (+ strand)")
    parser.add_argument("--qs-minus",   required=True, help="dRNA bedGraph (- strand)")
    parser.add_argument("--output",     required=True, help="Output directory")
    parser.add_argument("--verbose",    action="store_true")
    args = parser.parse_args()

    beds_by_method = {}
    for entry in args.bed:
        if ":" not in entry:
            print(f"WARNING: skipping malformed entry '{entry}'", file=sys.stderr)
            continue
        label, path = entry.split(":", 1)
        if not Path(path).exists():
            print(f"WARNING: file not found: {path}", file=sys.stderr)
            continue
        isoforms = parse_isoforms(path)
        if isoforms:
            beds_by_method[label] = isoforms
            if args.verbose:
                print(f"  {label}: {len(isoforms)} isoforms")

    if not beds_by_method:
        print("No BED data loaded — skipping", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print("  Loading signal tracks...")
    cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
        args.cage_plus, args.cage_minus, args.qs_plus, args.qs_minus,
    )

    plot_end_signal_scatter(beds_by_method, cage_p, cage_m, qs_p, qs_m, args.output)
    print(f"Saved end-signal scatter to {args.output}")


if __name__ == "__main__":
    main()
