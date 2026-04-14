#!/usr/bin/env python3
"""
end_signal_heatmap.py — Smarca4-style signal density heatmap at called isoform ends.

Each row = one isoform, sorted by read support (descending).
Columns = bp offset from the called TSS (left panel) or TTS (right panel).
Colour = CAGE bedGraph value (TSS) or dRNA bedGraph value (TTS).

A method with precise end-calling shows a bright band running vertically
at offset 0.  Methods that over-segment or misplace ends show diffuse or
off-centre signal.  One figure is produced per method; a summary figure
tiles all methods for direct comparison.

Outputs (all in --output dir):
    heatmap_{method}.png          — per-method 2-panel (5′ / 3′) heatmap
    heatmap_summary_5prime.png    — all methods tiled, TSS panel only
    heatmap_summary_3prime.png    — all methods tiled, TTS panel only

Usage:
    python end_signal_heatmap.py \\
        --bed          label1:bed1.bed label2:bed2.bed ... \\
        --read-map     label1:map1.txt label2:map2.txt ... \\
        --cage-plus    cage_plus.bg  --cage-minus  cage_minus.bg \\
        --qs-plus      qs_plus.bg    --qs-minus    qs_minus.bg \\
        --output       output_dir/   \\
        [--flank 300]  [--bin-size 5] [--max-isoforms 2000]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.signal import savgol_filter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2
    from signal_utils import (
        parse_isoforms, load_signal_tracks, BedGraphTrack,
        load_read_map, lookup_read_count,
    )
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2
    from evaluation.signal_utils import (
        parse_isoforms, load_signal_tracks, BedGraphTrack,
        load_read_map, lookup_read_count,
    )

apply_rc()


# ── Matrix builder ────────────────────────────────────────────────────────────

def _build_matrix(
    isos: List[dict],
    read_counts: Dict[str, int],
    track_plus: BedGraphTrack,
    track_minus: BedGraphTrack,
    flank: int,
    bin_size: int,
    max_isoforms: int,
    end: str,   # "tss" or "tts"
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (matrix, read_count_array) sorted by read support descending.

    matrix shape: (n_isoforms, n_bins)
    """
    n_bins = (2 * flank) // bin_size
    rows = []

    for iso in isos:
        ch, strand = iso["chrom"], iso["strand"]
        if end == "tss":
            pos = iso["start"] if strand == "+" else iso["end"]
            track = track_plus if strand == "+" else track_minus
        else:
            pos = iso["end"] if strand == "+" else iso["start"]
            track = track_plus if strand == "+" else track_minus

        rc = lookup_read_count(iso["name"], read_counts)
        raw = track.window_values(ch, pos - flank, pos + flank)

        # Reverse minus-strand so upstream is always left
        if strand == "-":
            raw = raw[::-1]

        # Bin-average
        usable = (len(raw) // bin_size) * bin_size
        if usable < n_bins * bin_size:
            continue  # window fell off chromosome edge
        binned = raw[:usable].reshape(-1, bin_size).mean(axis=1)[:n_bins]
        rows.append((rc, binned))

    if not rows:
        return np.zeros((0, n_bins), dtype=np.float32), np.array([], dtype=int)

    # Sort by read count descending
    rows.sort(key=lambda x: -x[0])
    if max_isoforms and len(rows) > max_isoforms:
        rows = rows[:max_isoforms]

    matrix = np.stack([r[1] for r in rows]).astype(np.float32)
    counts = np.array([r[0] for r in rows], dtype=int)
    return matrix, counts


# ── Plotting helpers ──────────────────────────────────────────────────────────

def _smooth(arr: np.ndarray, window: int = 11) -> np.ndarray:
    if len(arr) < window:
        return arr
    return savgol_filter(arr, window_length=window, polyorder=3)


def _panel(fig, gs_inner, matrix: np.ndarray, flank: int, bin_size: int,
           cmap: str, vmax: float, end_label: str, signal_label: str,
           show_ylabel: bool) -> None:
    """Draw metaplot + heatmap into a 2-row GridSpecFromSubplotSpec."""
    n_bins = (2 * flank) // bin_size
    x = np.linspace(-flank, flank, n_bins, endpoint=False) + bin_size / 2
    n_iso = matrix.shape[0]

    ax_meta = fig.add_subplot(gs_inner[0])
    ax_heat = fig.add_subplot(gs_inner[1], sharex=ax_meta)

    # Metaplot — clip to vmax so y-axis is on the same scale as the colorbar
    col_mean = np.nanmean(np.clip(matrix, 0, vmax), axis=0)
    ax_meta.fill_between(x, col_mean, alpha=0.25,
                         color="steelblue" if cmap == "Blues" else "coral")
    ax_meta.plot(x, _smooth(col_mean),
                 color="steelblue" if cmap == "Blues" else "coral",
                 linewidth=1.2)
    ax_meta.axvline(0, color="0.3", linewidth=0.7, linestyle="--")
    ax_meta.set_xlim(-flank, flank)
    ax_meta.set_ylim(bottom=0)
    ax_meta.tick_params(labelbottom=False, bottom=False)
    ax_meta.spines[["top", "right"]].set_visible(False)
    ax_meta.set_title(f"{end_label} ({signal_label})", fontsize=7)

    # Heatmap
    im = ax_heat.imshow(
        matrix,
        aspect="auto",
        cmap=cmap,
        vmin=0, vmax=vmax,
        interpolation="nearest",
        extent=[-flank, flank, n_iso, 0],
    )
    ax_heat.axvline(0, color="white", linewidth=0.6, linestyle="--", alpha=0.5)

    x_ticks = np.array([-flank, -flank // 2, 0, flank // 2, flank])
    ax_heat.set_xticks(x_ticks)
    ax_heat.set_xticklabels([str(t) for t in x_ticks], fontsize=6)
    ax_heat.set_xlabel(f"Offset from {end_label} (bp)", fontsize=6)

    if show_ylabel:
        ax_heat.set_ylabel(f"{n_iso} isoforms\n(↓ read support)", fontsize=6)

    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.03, pad=0.02)
    cbar.set_label(signal_label, fontsize=6)
    cbar.ax.tick_params(labelsize=5)


def plot_method_heatmap(
    method: str,
    tss_matrix: np.ndarray,
    tts_matrix: np.ndarray,
    flank: int,
    bin_size: int,
    output_path: Path,
) -> None:
    """Single per-method figure with 5′ and 3′ panels side by side."""
    has_tss = tss_matrix is not None and tss_matrix.shape[0] > 0
    has_tts = tts_matrix is not None and tts_matrix.shape[0] > 0
    n_cols = (1 if has_tss else 0) + (1 if has_tts else 0)
    if n_cols == 0:
        return

    vmax_tss = float(np.percentile(tss_matrix, 95)) if has_tss else 1.0
    vmax_tts = float(np.percentile(tts_matrix, 95)) if has_tts else 1.0
    # Avoid zero vmax
    vmax_tss = max(vmax_tss, 1e-6)
    vmax_tts = max(vmax_tts, 1e-6)

    fig = plt.figure(figsize=(4.5 * n_cols, 6))
    fig.suptitle(method, fontsize=8, fontweight="bold")
    outer = gridspec.GridSpec(1, n_cols, figure=fig, wspace=0.4)

    col = 0
    if has_tss:
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[col], height_ratios=[1, 4], hspace=0.05)
        _panel(fig, inner, tss_matrix, flank, bin_size,
               "Blues", vmax_tss, "TSS", "CAGE", show_ylabel=True)
        col += 1
    if has_tts:
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[col], height_ratios=[1, 4], hspace=0.05)
        _panel(fig, inner, tts_matrix, flank, bin_size,
               "Reds", vmax_tts, "TTS", "dRNA", show_ylabel=(col == 0))

    fig.tight_layout()
    savefig(fig, output_path)


def plot_summary_heatmaps(
    matrices: Dict[str, np.ndarray],
    flank: int,
    bin_size: int,
    end_label: str,
    signal_label: str,
    cmap: str,
    output_path: Path,
) -> None:
    """Tile all methods into one figure (one column per method)."""
    methods = [m for m, mat in matrices.items()
               if mat is not None and mat.shape[0] > 0]
    if not methods:
        return

    n_cols = len(methods)
    # Shared vmax across all methods for fair comparison
    all_vals = np.concatenate([matrices[m].ravel() for m in methods])
    vmax = float(np.percentile(all_vals, 95))
    vmax = max(vmax, 1e-6)

    fig = plt.figure(figsize=(3.5 * n_cols, 6))
    outer = gridspec.GridSpec(1, n_cols, figure=fig, wspace=0.35)

    for col, method in enumerate(methods):
        mat = matrices[method]
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[col], height_ratios=[1, 4], hspace=0.05)
        _panel(fig, inner, mat, flank, bin_size,
               cmap, vmax, end_label, signal_label,
               show_ylabel=(col == 0))
        # Add method name as column title via the metaplot axes
        ax0 = fig.axes[col * 2]   # metaplot ax for this column
        ax0.set_title(f"{method}\n{end_label} ({signal_label})", fontsize=6)

    fig.tight_layout()
    savefig(fig, output_path)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed",          nargs="+", required=True,
                   help="label:path pairs for BED12/GTF isoform files")
    p.add_argument("--read-map",     nargs="+", default=[],
                   help="label:path pairs for isoform read-map files (for row sorting)")
    p.add_argument("--cage-plus",    required=True)
    p.add_argument("--cage-minus",   required=True)
    p.add_argument("--qs-plus",      required=True)
    p.add_argument("--qs-minus",     required=True)
    p.add_argument("--output",       required=True, help="Output directory")
    p.add_argument("--flank",        type=int, default=300,
                   help="Flank in bp each side of called end (default 300)")
    p.add_argument("--bin-size",     type=int, default=5,
                   help="Bin size in bp (default 5)")
    p.add_argument("--max-isoforms", type=int, default=2000,
                   help="Max isoforms per method (default 2000, 0=all)")
    p.add_argument("--verbose",      action="store_true")
    args = p.parse_args()

    # Parse BED files
    beds_by_method: Dict[str, List[dict]] = {}
    for entry in args.bed:
        if ":" not in entry:
            continue
        label, path = entry.split(":", 1)
        if not Path(path).exists():
            print(f"WARNING: file not found: {path}", file=sys.stderr)
            continue
        isos = parse_isoforms(path)
        if isos:
            beds_by_method[label] = isos
            if args.verbose:
                print(f"  {label}: {len(isos)} isoforms")

    if not beds_by_method:
        print("No isoform data loaded — skipping", file=sys.stderr)
        sys.exit(1)

    # Parse read-maps
    read_maps: Dict[str, Dict[str, int]] = {}
    for entry in args.read_map:
        if ":" not in entry:
            continue
        label, path = entry.split(":", 1)
        if not Path(path).exists():
            continue
        rc = load_read_map(path)
        if rc:
            read_maps[label] = rc

    if args.verbose:
        print("Loading signal tracks...")
    cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
        args.cage_plus, args.cage_minus, args.qs_plus, args.qs_minus,
    )

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    tss_matrices: Dict[str, np.ndarray] = {}
    tts_matrices: Dict[str, np.ndarray] = {}

    for method, isos in beds_by_method.items():
        if args.verbose:
            print(f"  Building matrices for {method}...")
        rc = read_maps.get(method, {})

        tss_mat, _ = _build_matrix(isos, rc, cage_p, cage_m,
                                   args.flank, args.bin_size,
                                   args.max_isoforms, "tss")
        tts_mat, _ = _build_matrix(isos, rc, qs_p, qs_m,
                                   args.flank, args.bin_size,
                                   args.max_isoforms, "tts")
        tss_matrices[method] = tss_mat
        tts_matrices[method] = tts_mat

        # Per-method figure
        plot_method_heatmap(
            method=method,
            tss_matrix=tss_mat if tss_mat.shape[0] > 0 else None,
            tts_matrix=tts_mat if tts_mat.shape[0] > 0 else None,
            flank=args.flank,
            bin_size=args.bin_size,
            output_path=out / f"heatmap_{method}.png",
        )
        if args.verbose:
            print(f"    TSS: {tss_mat.shape[0]} rows, TTS: {tts_mat.shape[0]} rows")

    # Summary tiled figures
    plot_summary_heatmaps(
        matrices=tss_matrices,
        flank=args.flank, bin_size=args.bin_size,
        end_label="TSS", signal_label="CAGE", cmap="Blues",
        output_path=out / "heatmap_summary_5prime.png",
    )
    plot_summary_heatmaps(
        matrices=tts_matrices,
        flank=args.flank, bin_size=args.bin_size,
        end_label="TTS", signal_label="dRNA", cmap="Reds",
        output_path=out / "heatmap_summary_3prime.png",
    )

    print(f"Saved heatmap outputs to {args.output}")


if __name__ == "__main__":
    main()
