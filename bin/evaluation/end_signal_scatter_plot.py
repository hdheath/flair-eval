#!/usr/bin/env python3
"""
Isoform end-signal density scatter - per-method zero-safe hexbin panels.

For each method, plots TSS signal (CAGE) vs TTS signal (dRNA) for every
isoform.  Produces a publication-oriented hexbin plot with marginal
histograms:

  end_signal_scatter.png

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

from pub_style import savefig, W2
from signal_hexbin import signal_hexbin_with_marginals, write_signal_summary
from signal_utils import (
    parse_isoforms,
    load_signal_tracks,
    isoform_signal,
)


# ── Plotting ────────────────────────────────────────────────────────────────

def plot_end_signal_scatter(
    beds_by_method: dict,
    cage_p, cage_m, qs_p, qs_m,
    output_dir: str,
    *,
    signal_max: float = 100.0,
    auto_range: bool = False,
    hex_gridsize: int = 46,
):
    """Multi-panel count-coloured hexbin scatter of TSS vs TTS signal.

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

    ncols = min(3, n)
    nrows = (n + ncols - 1) // ncols
    fig_h = max(2.6, (W2 / ncols) * 0.95 * nrows)
    fig = plt.figure(figsize=(W2, fig_h))
    outer = fig.add_gridspec(
        nrows, ncols, left=0.08, right=0.985, bottom=0.10, top=0.94,
        wspace=0.38, hspace=0.46,
    )

    rows = []
    for idx, m in enumerate(methods):
        ss = outer[idx // ncols, idx % ncols]
        isoforms = beds_by_method[m]

        sigs = [isoform_signal(iso, cage_p, cage_m, qs_p, qs_m) for iso in isoforms]
        tss = np.array([s[0] for s in sigs], dtype=float)
        tts = np.array([s[1] for s in sigs], dtype=float)

        fixed_range = None if auto_range else (signal_max, signal_max)
        _, info = signal_hexbin_with_marginals(
            tts,
            tss,
            fig=fig,
            gs=ss,
            xlabel="TTS signal (dRNA TPM)",
            ylabel="TSS signal (CAGE TPM)",
            title=m,
            fixed_range=fixed_range,
            auto_range=auto_range,
            hex_gridsize=hex_gridsize,
        )
        row = {
            k: v for k, v in info.items()
            if k not in {"ax_sc", "ax_top", "ax_right", "cax"}
        }
        row.update(method=m, n_isoforms=len(isoforms))
        rows.append(row)

    savefig(fig, output_dir / "end_signal_scatter.png")
    write_signal_summary(output_dir / "end_signal_scatter_summary.tsv", rows)


# ── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--bed", nargs="+", required=True,
        help="label:path pairs for BED12/GTF/GFF isoform files",
    )
    parser.add_argument("--cage-plus",  required=True, help="CAGE bedGraph (+ strand)")
    parser.add_argument("--cage-minus", required=True, help="CAGE bedGraph (- strand)")
    parser.add_argument("--qs-plus",    required=True, help="dRNA bedGraph (+ strand)")
    parser.add_argument("--qs-minus",   required=True, help="dRNA bedGraph (- strand)")
    parser.add_argument("--output",     required=True, help="Output directory")
    parser.add_argument(
        "--signal-max", type=float, default=100.0,
        help="Raw TPM upper axis limit when --auto-range is not set",
    )
    parser.add_argument(
        "--auto-range", action="store_true",
        help="Use per-panel data-driven signal axis limits instead of --signal-max",
    )
    parser.add_argument(
        "--hex-gridsize", type=int, default=46,
        help="Number of hexagons along the x-axis",
    )
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
        print("  Loading signal tracks...", flush=True)
    cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
        args.cage_plus, args.cage_minus, args.qs_plus, args.qs_minus,
    )

    output_dir = Path(args.output)
    plot_end_signal_scatter(
        beds_by_method,
        cage_p, cage_m, qs_p, qs_m,
        output_dir,
        signal_max=args.signal_max,
        auto_range=args.auto_range,
        hex_gridsize=args.hex_gridsize,
    )
    print(f"Saved end-signal scatter to {args.output}")


if __name__ == "__main__":
    main()
