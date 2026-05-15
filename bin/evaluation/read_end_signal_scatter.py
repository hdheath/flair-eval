#!/usr/bin/env python3
"""
Read end-signal density scatter - per-sample zero-safe hexbin panels.

For each sample, plots TSS signal (CAGE) vs TTS signal (dRNA) for every
aligned *read* (not assembled isoform), coloured by read count per hexbin.
Shows what the raw data looks like against orthogonal signal before assembly.

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

from pub_style import savefig, W2
from signal_hexbin import signal_hexbin_with_marginals, write_signal_summary
from signal_utils import (
    parse_bed12,
    load_signal_tracks,
    isoform_signal,
)

RNG_SEED = 42


def subsample(reads, max_n=0):
    """Deterministic random subsample when max_n is positive."""
    if not max_n or max_n <= 0 or len(reads) <= max_n:
        return reads
    rng = np.random.default_rng(RNG_SEED)
    idx = rng.choice(len(reads), size=max_n, replace=False)
    return [reads[i] for i in idx]


def plot_read_end_signal_scatter(
    reads_by_sample: dict,
    cage_p, cage_m, qs_p, qs_m,
    output_dir: str,
    *,
    max_reads: int = 0,
    signal_max: float = 100.0,
    auto_range: bool = False,
    signal_stat: str = "mean",
):
    """Multi-panel count-coloured hexbin scatter of TSS vs TTS signal per read."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    samples = list(reads_by_sample.keys())
    n = len(samples)
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
    for idx, sample in enumerate(samples):
        ss = outer[idx // ncols, idx % ncols]
        reads = reads_by_sample[sample]
        total = len(reads)
        reads = subsample(reads, max_n=max_reads)

        sigs = [
            isoform_signal(r, cage_p, cage_m, qs_p, qs_m, stat=signal_stat)
            for r in reads
        ]
        tss = np.array([s[0] for s in sigs], dtype=float)
        tts = np.array([s[1] for s in sigs], dtype=float)

        fixed_range = None if auto_range else (signal_max, signal_max)
        stat_label = {
            "mean": "mean TPM",
            "max": "max TPM",
            "sum": "sum TPM-bp",
        }[signal_stat]
        xlabel = f"TTS signal (dRNA {stat_label})" if idx // ncols == nrows - 1 else ""
        ylabel = f"TSS signal (CAGE {stat_label})" if idx % ncols == 0 else ""
        _, info = signal_hexbin_with_marginals(
            tts,
            tss,
            fig=fig,
            gs=ss,
            xlabel=xlabel,
            ylabel=ylabel,
            title=sample,
            fixed_range=fixed_range,
            auto_range=auto_range,
            hex_gridsize=46,
        )
        row = {
            k: v for k, v in info.items()
            if k not in {"ax_sc", "ax_top", "ax_right", "cax"}
        }
        row.update(
            sample=sample,
            total_reads=total,
            plotted_reads=len(reads),
            signal_stat=signal_stat,
        )
        rows.append(row)

    savefig(fig, output_dir / "read_end_signal_scatter.png")
    write_signal_summary(output_dir / "read_end_signal_scatter_summary.tsv", rows)


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
    parser.add_argument("--qs-plus",    required=True, help="dRNA bedGraph (+ strand)")
    parser.add_argument("--qs-minus",   required=True, help="dRNA bedGraph (- strand)")
    parser.add_argument("--output",     required=True, help="Output directory")
    parser.add_argument(
        "--max-reads", type=int, default=0,
        help="Deterministically subsample each read BED to this many reads; 0 uses all reads",
    )
    parser.add_argument(
        "--signal-max", type=float, default=100.0,
        help="Raw TPM upper axis limit when --auto-range is not set",
    )
    parser.add_argument(
        "--auto-range", action="store_true",
        help="Use per-panel data-driven signal axis limits instead of --signal-max",
    )
    parser.add_argument(
        "--signal-stat", choices=("mean", "max", "sum"), default="mean",
        help=(
            "Signal statistic over the +/-50 bp read-end window. Use max for "
            "any-support diagnostics; mean preserves the historical averaging."
        ),
    )
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

    plot_read_end_signal_scatter(
        reads_by_sample,
        cage_p, cage_m, qs_p, qs_m,
        args.output,
        max_reads=args.max_reads,
        signal_max=args.signal_max,
        auto_range=args.auto_range,
        signal_stat=args.signal_stat,
    )
    print(f"Saved read end-signal scatter to {args.output}")


if __name__ == "__main__":
    main()
