#!/usr/bin/env python3
"""Read vs isoform exon-length distribution plot.

Compares spliced read lengths from a BED12 reads file with spliced isoform
lengths from one or more BED12/GTF/GFF assembler outputs.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pub_style import ModeStyler, W1, W2, savefig, style_ax
from signal_utils import parse_isoforms


def read_exon_lengths(reads_bed: Path) -> np.ndarray:
    """Return one spliced length per read name from a BED12 reads file."""
    lengths = {}
    with open(reads_bed) as f:
        for line in f:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            name = cols[3]
            if name in lengths:
                continue
            try:
                block_sizes = [int(x) for x in cols[10].rstrip(",").split(",") if x]
            except ValueError:
                continue
            exon_len = sum(block_sizes)
            if exon_len > 0:
                lengths[name] = exon_len
    return np.asarray(list(lengths.values()), dtype=int)


def isoform_exon_lengths(path: Path) -> np.ndarray:
    """Return spliced isoform lengths from BED12/GTF/GFF."""
    vals = []
    for iso in parse_isoforms(path):
        length = int(iso.get("spliced_len") or (iso["end"] - iso["start"]))
        if length > 0:
            vals.append(length)
    return np.asarray(vals, dtype=int)


def load_supported_ids(counts_path: Path | None, min_support: int = 1) -> set[str] | None:
    """Return transcript IDs with count >= min_support from assembler count tables.

    Handles IsoQuant (`feature_id/count`), Bambu (`TXNAME/GENEID/sample_count`),
    and FLAIR count tables by using column 1 as the transcript ID and summing
    numeric fields to its right.
    """
    if counts_path is None or not counts_path.exists() or counts_path.stat().st_size == 0:
        return None
    supported: set[str] = set()
    n_parsed = 0
    with open(counts_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            numeric_values = []
            for value in parts[1:]:
                try:
                    numeric_values.append(float(value))
                except ValueError:
                    continue
            if not numeric_values:
                continue
            n_parsed += 1
            if sum(numeric_values) >= min_support:
                supported.add(parts[0])
    return supported if n_parsed else None


def is_supported_name(name: str, supported_ids: set[str]) -> bool:
    if name in supported_ids:
        return True
    if "_" in name and name.split("_", 1)[0] in supported_ids:
        return True
    return False


def isoform_exon_lengths_supported(
    path: Path,
    supported_ids: set[str] | None = None,
) -> np.ndarray:
    """Return spliced isoform lengths, optionally filtering by count support."""
    vals = []
    for iso in parse_isoforms(path):
        name = str(iso.get("name") or "")
        if supported_ids is not None and not is_supported_name(name, supported_ids):
            continue
        length = int(iso.get("spliced_len") or (iso["end"] - iso["start"]))
        if length > 0:
            vals.append(length)
    return np.asarray(vals, dtype=int)


def _parse_label_path_pairs(entries) -> dict[str, Path]:
    pairs = {}
    for entry in entries or []:
        if ":" not in entry:
            print(f"WARNING: skipping malformed entry '{entry}' (expected label:path)",
                  file=sys.stderr)
            continue
        label, path = entry.split(":", 1)
        pairs[label] = Path(path)
    return pairs


def _write_summary(path: Path, read_lens: np.ndarray, iso_lens_by_method: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as out:
        out.write("label\ttype\tn\tmedian_bp\tmean_bp\tp90_bp\n")
        rows = [("Reads", "reads", read_lens)]
        rows.extend((label, "isoforms", vals) for label, vals in iso_lens_by_method.items())
        for label, kind, vals in rows:
            vals = vals[np.isfinite(vals) & (vals > 0)]
            if vals.size:
                median = float(np.median(vals))
                mean = float(np.mean(vals))
                p90 = float(np.percentile(vals, 90))
            else:
                median = mean = p90 = 0.0
            out.write(f"{label}\t{kind}\t{vals.size}\t{median:.1f}\t{mean:.1f}\t{p90:.1f}\n")


def plot_exon_length_distributions(
    read_lens: np.ndarray,
    iso_lens_by_method: dict[str, np.ndarray],
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    non_empty = {m: vals for m, vals in iso_lens_by_method.items() if vals.size}
    if read_lens.size == 0 or not non_empty:
        print("No exon-length data to plot", file=sys.stderr)
        return

    all_vals = [read_lens] + list(non_empty.values())
    combined = np.concatenate(all_vals)
    lo = max(1, int(np.floor(combined.min())))
    hi = max(lo + 1, int(np.ceil(combined.max() * 1.05)))
    bins = np.logspace(np.log10(lo), np.log10(hi), 56)

    fig_width = W2 if len(non_empty) > 3 else W1
    fig, ax = plt.subplots(figsize=(fig_width, max(2.55, fig_width * 0.58)))

    read_counts, _ = np.histogram(read_lens, bins=bins)
    read_color = "#0072B2"
    ax.stairs(read_counts, bins, fill=True, color=read_color, alpha=0.16,
              linewidth=0.0, baseline=0)
    ax.stairs(read_counts, bins, color=read_color, linewidth=1.0, label="Reads")
    ax.axvline(np.median(read_lens), color=read_color, linewidth=0.7,
               linestyle=(0, (2, 2)), alpha=0.8)

    styler = ModeStyler(list(non_empty.keys()))
    for method, vals in non_empty.items():
        counts, _ = np.histogram(vals, bins=bins)
        color = styler.color(method)
        ax.stairs(counts, bins, color=color, linewidth=0.9, label=method)
        ax.axvline(np.median(vals), color=color, linewidth=0.55,
                   linestyle=(0, (2, 2)), alpha=0.75)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(bottom=0.8)
    style_ax(ax, xlabel="Spliced length (bp)", ylabel="Count", faint_y_grid=True)
    ax.tick_params(axis="both", which="both", direction="out", pad=2)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False,
              fontsize=6.5, handlelength=1.6, borderaxespad=0.2)
    fig.tight_layout(rect=(0, 0, 0.82, 1))
    savefig(fig, output_dir / "read_vs_isoform_exon_lengths.png")
    _write_summary(output_dir / "read_vs_isoform_exon_lengths.tsv", read_lens, non_empty)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads-bed", required=True, type=Path,
                        help="BED12 file with aligned reads")
    parser.add_argument("--bed", nargs="+", required=True,
                        help="label:path pairs for BED12/GTF/GFF isoform files")
    parser.add_argument("--counts", nargs="*", default=[],
                        help="Optional label:path pairs for transcript count support filtering")
    parser.add_argument("--min-support", type=int, default=1,
                        help="Minimum transcript count retained when --counts is supplied")
    parser.add_argument("--output", required=True, type=Path,
                        help="Output directory")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if not args.reads_bed.exists():
        print(f"Reads BED not found: {args.reads_bed}", file=sys.stderr)
        sys.exit(0)

    if args.verbose:
        print(f"Loading read exon lengths from {args.reads_bed}", file=sys.stderr)
    read_lens = read_exon_lengths(args.reads_bed)

    counts_by_method = _parse_label_path_pairs(args.counts)
    iso_lens_by_method = {}
    for label, path in _parse_label_path_pairs(args.bed).items():
        if not path.exists():
            print(f"WARNING: file not found: {path}", file=sys.stderr)
            continue
        supported_ids = load_supported_ids(counts_by_method.get(label), args.min_support)
        vals = isoform_exon_lengths_supported(path, supported_ids)
        iso_lens_by_method[label] = vals
        if args.verbose:
            median = int(np.median(vals)) if vals.size else 0
            support_msg = " count-filtered" if supported_ids is not None else ""
            print(f"  {label}: n={vals.size:,} median={median:,} bp{support_msg}",
                  file=sys.stderr)

    plot_exon_length_distributions(read_lens, iso_lens_by_method, args.output)


if __name__ == "__main__":
    main()
