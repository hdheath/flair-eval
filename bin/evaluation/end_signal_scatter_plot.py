#!/usr/bin/env python3
"""
Isoform end-signal density scatter — per-method KDE-coloured panels.

For each method, plots TSS signal (CAGE) vs TTS signal (dRNA) for every
isoform.  Produces two figures:

  end_signal_scatter.png        KDE-density coloured scatter (relative density)
  end_signal_scatter_counts.png Hexbin 2D histogram coloured by log10 isoform
                                count per bin — shows absolute concentration of
                                isoforms on real signal peaks.

Requires BED12 isoform files plus four bedGraph signal tracks
(CAGE +/- strand, dRNA +/- strand).

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


def plot_end_signal_hexbin(
    beds_by_method: dict,
    cage_p, cage_m, qs_p, qs_m,
    output_dir: Path,
    gridsize: int = 60,
):
    """Multi-panel 2D hexbin count scatter of TSS vs TTS signal.

    Same axes as the KDE scatter but colour encodes raw isoform counts
    per bin (log scale), making absolute density differences visible.
    """
    output_dir = Path(output_dir)
    methods = list(beds_by_method.keys())
    n = len(methods)
    if n == 0:
        return

    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    pw = W2 / ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(W2, pw * nrows), squeeze=False)

    eps = 1e-3
    # Collect all signal values first to set a shared colour scale
    all_tss, all_tts = [], []
    per_method: list[tuple[np.ndarray, np.ndarray]] = []
    for m in methods:
        sigs = [isoform_signal(iso, cage_p, cage_m, qs_p, qs_m)
                for iso in beds_by_method[m]]
        tss = np.log10(np.array([s[0] for s in sigs]) + eps)
        tts = np.log10(np.array([s[1] for s in sigs]) + eps)
        per_method.append((tss, tts))
        all_tss.extend(tss)
        all_tts.extend(tts)

    x_range = (min(all_tts), max(all_tts))
    y_range = (min(all_tss), max(all_tss))

    hb_ref = None  # store one hexbin for the shared colourbar

    for idx, m in enumerate(methods):
        ax = axes[idx // ncols][idx % ncols]
        tss, tts = per_method[idx]

        hb = ax.hexbin(
            tts, tss,
            gridsize=gridsize,
            bins="log",          # log-scale count colouring
            cmap="YlOrRd",
            extent=[x_range[0], x_range[1], y_range[0], y_range[1]],
            linewidths=0.0,
        )
        if hb_ref is None:
            hb_ref = hb

        # Restore readable tick labels (values are log10 of signal)
        for axis, rng in [(ax.xaxis, x_range), (ax.yaxis, y_range)]:
            tks = np.linspace(rng[0], rng[1], 5)
            axis.set_ticks(tks)
            axis.set_ticklabels([f"$10^{{{t:.1f}}}$" for t in tks], fontsize=4)

        if idx >= n - ncols:
            ax.set_xlabel("TTS signal (dRNA)", fontsize=7)
        if idx % ncols == 0:
            ax.set_ylabel("TSS signal (CAGE)", fontsize=7)

        style_ax(ax)
        ax.text(0.04, 0.96, m, transform=ax.transAxes,
                ha="left", va="top", fontsize=6, fontweight="bold")
        ax.text(0.96, 0.04, f"n = {len(tss):,}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=5, color="#666666")

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.subplots_adjust(left=0.07, right=0.90, bottom=0.10, top=0.97,
                        hspace=0.35, wspace=0.35)
    if hb_ref is not None:
        cbar_ax = fig.add_axes([0.92, 0.15, 0.012, 0.7])
        cbar = fig.colorbar(hb_ref, cax=cbar_ax)
        cbar.set_label("log₁₀ isoform count", fontsize=6)
        cbar.ax.tick_params(labelsize=5, length=1.5, width=0.3)

    savefig(fig, output_dir / "end_signal_scatter_counts.png")


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

    # Cap any single method's isoform set at 300,000.  KDE density colouring
    # is O(N²); tools like IsoSeq can emit >1M isoforms, which makes the KDE
    # step never finish.  Subsampling at this scale doesn't lose visual
    # density information (the scatter is already saturated past ~50K points).
    SCATTER_CAP = 300_000
    rng = np.random.default_rng(0)
    for label, isoforms in list(beds_by_method.items()):
        if len(isoforms) > SCATTER_CAP:
            idx = rng.choice(len(isoforms), size=SCATTER_CAP, replace=False)
            beds_by_method[label] = [isoforms[i] for i in idx]
            if args.verbose:
                print(f"  {label}: subsampled {len(isoforms):,} -> {SCATTER_CAP:,} "
                      f"(scatter density cap)", flush=True)

    if args.verbose:
        print("  Loading signal tracks...", flush=True)
    cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
        args.cage_plus, args.cage_minus, args.qs_plus, args.qs_minus,
    )

    output_dir = Path(args.output)
    plot_end_signal_scatter(beds_by_method, cage_p, cage_m, qs_p, qs_m, output_dir)
    plot_end_signal_hexbin(beds_by_method, cage_p, cage_m, qs_p, qs_m, output_dir)
    print(f"Saved end-signal scatter to {args.output}")


if __name__ == "__main__":
    main()
