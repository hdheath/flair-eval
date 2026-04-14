#!/usr/bin/env python3
"""
read_end_heatmap.py — Read-end distribution heatmap at called isoform ends.

The long-read analogue of end_signal_heatmap.py.  Instead of showing CAGE/dRNA
bedGraph signal at each isoform end, shows where the *reads assigned to each
isoform* actually land relative to the called end.

Each row = one isoform, sorted by read count descending.
Columns = bp offset from the called TSS (left panel) or TTS (right panel).
Colour = read density (reads per bp bin, normalised to reads-per-isoform).

A method with precise end-calling shows a tight band at offset 0 — reads cluster
around the called position.  Methods that over-segment or misplace ends show
diffuse distributions or off-centre peaks.

This plot is purely from the long-read data itself, with no orthogonal signal.

Outputs:
    read_end_heatmap_{method}.png     — per-method 2-panel (5′ / 3′) heatmap
    read_end_heatmap_summary_5prime.png
    read_end_heatmap_summary_3prime.png

Usage:
    python read_end_heatmap.py \\
        --bed          label1:isoforms1.bed label2:isoforms2.bed ... \\
        --read-map     label1:map1.txt      label2:map2.txt ...      \\
        --reads-bed    label1:reads1.bed    label2:reads2.bed ...    \\
        --output       output_dir/          \\
        [--flank 300]  [--bin-size 5]  [--max-isoforms 2000]  [--min-reads 2]
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.signal import savgol_filter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, str(Path(__file__).resolve().parent))

from pub_style import apply_rc, style_ax, savefig, W1, W2
from signal_utils import parse_isoforms

apply_rc()


# ── Parsers ───────────────────────────────────────────────────────────────────

def _parse_reads_bed(path: Path) -> Dict[str, dict]:
    """Parse a read_audit BED12 → {read_name: {chrom, strand, tss, tts}}.

    Strand-aware: tss = 5' end, tts = 3' end of the read.
    Skips track/browser header lines and keeps only the first alignment per
    read name (primary alignment).
    """
    reads: Dict[str, dict] = {}
    with open(path) as f:
        for line in f:
            if line.startswith(("track", "browser", "#")):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            name = cols[3]
            if name in reads:
                continue
            chrom  = cols[0]
            start  = int(cols[1])
            end    = int(cols[2])
            strand = cols[5]
            tss, tts = (start, end) if strand == "+" else (end, start)
            reads[name] = {"chrom": chrom, "strand": strand, "tss": tss, "tts": tts}
    return reads


# ── Matrix builder ────────────────────────────────────────────────────────────

def _parse_iso_to_reads(map_path: Path) -> Dict[str, List[str]]:
    """Parse isoform.read.map.txt → {isoform_id: [read_id, ...]}."""
    iso_to_reads: Dict[str, List[str]] = {}
    with open(map_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue
            iso_id = parts[0]
            reads = [r.strip() for r in parts[1].split(",") if r.strip()]
            if reads:
                iso_to_reads[iso_id] = reads
    return iso_to_reads


def _build_matrix(
    isos: List[dict],
    iso_to_reads: Dict[str, List[str]],
    read_ends: Dict[str, dict],
    flank: int,
    bin_size: int,
    max_isoforms: int,
    min_reads: int,
    end: str,  # "tss" or "tts"
) -> Tuple[np.ndarray, np.ndarray]:
    """Build (matrix, read_count_array) sorted by read count descending.

    matrix shape: (n_isoforms, n_bins)
    Each cell = number of reads whose end falls in that bp bin, normalised
    by the total reads assigned to that isoform (so each row sums to 1).
    """
    n_bins = (2 * flank) // bin_size
    rows = []

    for iso in isos:
        iso_name = iso["name"]
        read_ids = iso_to_reads.get(iso_name, [])
        if len(read_ids) < min_reads:
            continue

        iso_end_pos = iso["tss"] if end == "tss" else iso["tts"]
        strand = iso["strand"]

        counts = np.zeros(n_bins, dtype=np.float32)
        n_placed = 0

        for rid in read_ids:
            r = read_ends.get(rid)
            if r is None:
                continue
            if r["chrom"] != iso["chrom"]:
                continue

            # Genomic position of the read's end
            read_end_pos = r["tss"] if end == "tss" else r["tts"]

            # Offset relative to the isoform's called end
            # Strand-aware: upstream is negative, downstream is positive
            if strand == "+":
                offset = read_end_pos - iso_end_pos
            else:
                offset = iso_end_pos - read_end_pos

            # Bin index
            bin_idx = int((offset + flank) // bin_size)
            if 0 <= bin_idx < n_bins:
                counts[bin_idx] += 1
                n_placed += 1

        if n_placed == 0:
            continue

        # Binary: 1 where any read end falls, 0 elsewhere
        rows.append((len(read_ids), (counts > 0).astype(np.float32)))

    if not rows:
        return np.zeros((0, n_bins), dtype=np.float32), np.array([], dtype=int)

    rows.sort(key=lambda x: -x[0])
    if max_isoforms and len(rows) > max_isoforms:
        rows = rows[:max_isoforms]

    matrix = np.stack([r[1] for r in rows]).astype(np.float32)
    counts = np.array([r[0] for r in rows], dtype=int)
    return matrix, counts


# ── Plotting ──────────────────────────────────────────────────────────────────

def _smooth(arr: np.ndarray, window: int = 11) -> np.ndarray:
    if len(arr) < window:
        return arr
    return savgol_filter(arr, window_length=window, polyorder=3)


def _panel(
    fig, gs_inner,
    matrix: np.ndarray,
    flank: int,
    bin_size: int,
    color: str,
    end_label: str,
    show_ylabel: bool,
) -> None:
    """Draw metaplot + binary heatmap into a 2-row GridSpecFromSubplotSpec.

    Matrix is binary (1 = at least one read end in this bin, 0 = none).
    Metaplot shows fraction of isoforms with a read in each bin.
    No colorbar — presence/absence needs no scale.
    """
    n_bins = (2 * flank) // bin_size
    x = np.linspace(-flank, flank, n_bins, endpoint=False) + bin_size / 2
    n_iso = matrix.shape[0]

    ax_meta = fig.add_subplot(gs_inner[0])
    ax_heat = fig.add_subplot(gs_inner[1], sharex=ax_meta)

    # Metaplot: fraction of isoforms with a read end in each bin
    col_frac = np.mean(matrix, axis=0)
    ax_meta.fill_between(x, col_frac, alpha=0.25, color=color)
    ax_meta.plot(x, _smooth(col_frac), color=color, linewidth=1.2)
    ax_meta.axvline(0, color="0.3", linewidth=0.7, linestyle="--")
    ax_meta.set_xlim(-flank, flank)
    ax_meta.set_ylim(0, min(1.0, col_frac.max() * 1.2) if col_frac.max() > 0 else 1.0)
    ax_meta.tick_params(labelbottom=False, bottom=False)
    ax_meta.spines[["top", "right"]].set_visible(False)
    ax_meta.set_title(f"{end_label} (reads)", fontsize=7)
    ax_meta.set_ylabel("Fraction\nof isoforms", fontsize=5)

    # Heatmap — binary, no colorbar
    from matplotlib.colors import ListedColormap
    binary_cmap = ListedColormap(["white", color])
    ax_heat.imshow(
        matrix,
        aspect="auto",
        cmap=binary_cmap,
        vmin=0, vmax=1,
        interpolation="nearest",
        extent=[-flank, flank, n_iso, 0],
    )
    ax_heat.axvline(0, color="0.5", linewidth=0.6, linestyle="--", alpha=0.6)

    x_ticks = np.array([-flank, -flank // 2, 0, flank // 2, flank])
    ax_heat.set_xticks(x_ticks)
    ax_heat.set_xticklabels([str(t) for t in x_ticks], fontsize=6)
    ax_heat.set_xlabel(f"Offset from {end_label} (bp)", fontsize=6)

    if show_ylabel:
        ax_heat.set_ylabel(f"{n_iso} isoforms\n(↓ read support)", fontsize=6)


def plot_method_heatmap(
    method: str,
    tss_matrix: Optional[np.ndarray],
    tts_matrix: Optional[np.ndarray],
    flank: int,
    bin_size: int,
    output_path: Path,
) -> None:
    has_tss = tss_matrix is not None and tss_matrix.shape[0] > 0
    has_tts = tts_matrix is not None and tts_matrix.shape[0] > 0
    n_cols = (1 if has_tss else 0) + (1 if has_tts else 0)
    if n_cols == 0:
        return

    fig = plt.figure(figsize=(4.5 * n_cols, 6))
    fig.suptitle(method, fontsize=8, fontweight="bold")
    outer = gridspec.GridSpec(1, n_cols, figure=fig, wspace=0.4)

    col = 0
    if has_tss:
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[col], height_ratios=[1, 4], hspace=0.05)
        _panel(fig, inner, tss_matrix, flank, bin_size,
               "steelblue", "TSS", show_ylabel=True)
        col += 1
    if has_tts:
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[col], height_ratios=[1, 4], hspace=0.05)
        _panel(fig, inner, tts_matrix, flank, bin_size,
               "coral", "TTS", show_ylabel=(col == 0))

    fig.tight_layout()
    savefig(fig, output_path)


def plot_summary_heatmaps(
    matrices: Dict[str, np.ndarray],
    flank: int,
    bin_size: int,
    end_label: str,
    color: str,
    output_path: Path,
) -> None:
    methods = [m for m, mat in matrices.items()
               if mat is not None and mat.shape[0] > 0]
    if not methods:
        return

    n_cols = len(methods)
    fig = plt.figure(figsize=(3.5 * n_cols, 6))
    outer = gridspec.GridSpec(1, n_cols, figure=fig, wspace=0.35)

    for col, method in enumerate(methods):
        mat = matrices[method]
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[col], height_ratios=[1, 4], hspace=0.05)
        _panel(fig, inner, mat, flank, bin_size,
               color, end_label, show_ylabel=(col == 0))
        ax0 = fig.axes[col * 2]
        ax0.set_title(f"{method}\n{end_label} (reads)", fontsize=6)

    fig.tight_layout()
    savefig(fig, output_path)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--bed",        nargs="+", required=True,
                   help="label:path pairs for BED12 isoform files")
    p.add_argument("--read-map",   nargs="+", required=True,
                   help="label:path pairs for isoform.read.map.txt files")
    p.add_argument("--reads-bed",  required=True,
                   help="Single read_audit BED12 file shared across all methods")
    p.add_argument("--output",     required=True, help="Output directory")
    p.add_argument("--flank",      type=int, default=300,
                   help="Flank in bp each side of called end (default 300)")
    p.add_argument("--bin-size",   type=int, default=5,
                   help="Bin size in bp (default 5)")
    p.add_argument("--max-isoforms", type=int, default=2000,
                   help="Max isoforms per method (default 2000, 0=all)")
    p.add_argument("--min-reads",  type=int, default=2,
                   help="Min reads assigned to an isoform to include it (default 2)")
    p.add_argument("--verbose",    action="store_true")
    args = p.parse_args()

    def parse_label_paths(entries):
        result = {}
        for entry in entries:
            if ":" not in entry:
                print(f"WARNING: skipping malformed entry '{entry}' (no ':')", file=sys.stderr)
                continue
            label, path = entry.split(":", 1)
            if not Path(path).exists():
                print(f"WARNING: file not found: {path}", file=sys.stderr)
                continue
            result[label] = Path(path)
        return result

    isom_paths = parse_label_paths(args.bed)
    map_paths  = parse_label_paths(args.read_map)

    reads_bed_path = Path(args.reads_bed)
    if not reads_bed_path.exists():
        print(f"ERROR: reads BED not found: {reads_bed_path}", file=sys.stderr)
        sys.exit(1)

    labels = sorted(set(isom_paths) & set(map_paths))
    if not labels:
        print("No labels with both --bed and --read-map present", file=sys.stderr)
        sys.exit(1)

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    tss_matrices: Dict[str, np.ndarray] = {}
    tts_matrices: Dict[str, np.ndarray] = {}

    # Load the shared reads BED once — all methods share the same reads
    if args.verbose:
        print(f"Loading reads BED: {reads_bed_path}")
    read_ends = _parse_reads_bed(reads_bed_path)
    if args.verbose:
        print(f"  {len(read_ends)} reads loaded")

    for label in labels:
        if args.verbose:
            print(f"Processing {label}...")

        isos = parse_isoforms(isom_paths[label])
        iso_to_reads = _parse_iso_to_reads(map_paths[label])

        if args.verbose:
            print(f"  {len(isos)} isoforms, {len(iso_to_reads)} in map")

        tss_mat, _ = _build_matrix(
            isos, iso_to_reads, read_ends,
            args.flank, args.bin_size, args.max_isoforms, args.min_reads,
            end="tss",
        )
        tts_mat, _ = _build_matrix(
            isos, iso_to_reads, read_ends,
            args.flank, args.bin_size, args.max_isoforms, args.min_reads,
            end="tts",
        )
        tss_matrices[label] = tss_mat
        tts_matrices[label] = tts_mat

        if args.verbose:
            print(f"  TSS matrix: {tss_mat.shape[0]} isoforms, "
                  f"TTS matrix: {tts_mat.shape[0]} isoforms")

        plot_method_heatmap(
            method=label,
            tss_matrix=tss_mat if tss_mat.shape[0] > 0 else None,
            tts_matrix=tts_mat if tts_mat.shape[0] > 0 else None,
            flank=args.flank,
            bin_size=args.bin_size,
            output_path=out / f"read_end_heatmap_{label}.png",
        )

    plot_summary_heatmaps(
        matrices=tss_matrices,
        flank=args.flank, bin_size=args.bin_size,
        end_label="TSS", color="steelblue",
        output_path=out / "read_end_heatmap_summary_5prime.png",
    )
    plot_summary_heatmaps(
        matrices=tts_matrices,
        flank=args.flank, bin_size=args.bin_size,
        end_label="TTS", color="coral",
        output_path=out / "read_end_heatmap_summary_3prime.png",
    )

    print(f"Saved read-end heatmap outputs to {args.output}")


if __name__ == "__main__":
    main()
