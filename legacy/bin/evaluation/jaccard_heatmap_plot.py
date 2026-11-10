#!/usr/bin/env python3
"""
Pairwise Jaccard heatmaps for splice junctions and transcript ends.

Produces three plots per run:
  1. Splice-junction Jaccard heatmap (single-column)
  2. TSS + TTS Jaccard heatmaps (double-column, side-by-side)

Reads BED12 isoform files (one per method).  Splice junctions are derived
from block starts/sizes; transcript ends are binned at 50 bp resolution.

Usage:
    python jaccard_heatmap_plot.py \\
        --bed label1:bed1.bed label2:bed2.bed ... \\
        --output output_dir/
"""

import argparse
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from pub_style import style_ax, savefig, W1, W2, MODE_COLORS, PALETTE
from signal_utils import parse_isoforms, END_BIN


# ── Colour helpers ──────────────────────────────────────────────────────────

def _mode_color(mode: str) -> str:
    if mode in MODE_COLORS:
        return MODE_COLORS[mode]
    return PALETTE[hash(mode) % len(PALETTE)]


# ── Heatmap rendering ──────────────────────────────────────────────────────

def _draw_heatmap_ax(ax, mat, labels, cbar_label="", vmin=None, vmax=1.0):
    """Cividis heatmap with masked diagonal on an Axes."""
    n = len(labels)
    masked = np.ma.array(mat, mask=np.eye(n, dtype=bool))
    if vmin is None:
        vmin = max(0, float(np.nanmin(masked)) - 0.02)
    im = ax.imshow(masked, cmap="cividis", vmin=vmin, vmax=vmax, aspect="equal")
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, rotation=40, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels)
    ax.tick_params(length=2, width=0.35)
    cbar = ax.figure.colorbar(im, ax=ax, shrink=0.78, pad=0.03, aspect=22)
    cbar.ax.tick_params(labelsize=6, length=1.5, width=0.3)
    if cbar_label:
        cbar.set_label(cbar_label, fontsize=7)
    # Annotate off-diagonal cells with the Jaccard value
    for i in range(n):
        for j in range(n):
            if i != j:
                val = mat[i, j]
                # Use white text on dark cells, black on light
                text_color = "white" if val < (vmin + vmax) / 2 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=5.5, color=text_color)
    # Grey diagonal cells
    for i in range(n):
        ax.add_patch(plt.Rectangle(
            (i - 0.5, i - 0.5), 1, 1,
            facecolor="#e8e8e8", edgecolor="none", zorder=2,
        ))


def _jaccard_matrix(sets_by_method, methods):
    """Compute pairwise Jaccard matrix from {method: set}."""
    n = len(methods)
    mat = np.ones((n, n))
    for i, mi in enumerate(methods):
        for j, mj in enumerate(methods):
            if i != j:
                inter = len(sets_by_method[mi] & sets_by_method[mj])
                union = len(sets_by_method[mi] | sets_by_method[mj])
                mat[i, j] = inter / union if union else 0.0
    return mat


def _jaccard_display_label(label: str) -> str:
    """Compact method labels for heatmap axes."""
    replacements = {
        "TED-low-3-high-5": "TED-low-3\nhigh-5",
    }
    return replacements.get(label, label)


def _offdiag_values(*mats):
    vals = []
    for mat in mats:
        n = mat.shape[0]
        vals.extend(mat[~np.eye(n, dtype=bool)])
    return np.asarray(vals, dtype=float)


def _draw_jaccard_triangle(ax, mat, labels, title, *, norm, show_ylabels=True):
    """Lower-triangle Jaccard heatmap with compact annotations."""
    n = len(labels)
    mask = np.triu(np.ones((n, n), dtype=bool), k=0)
    masked = np.ma.array(mat, mask=mask)
    im = ax.imshow(masked, cmap="cividis", norm=norm, aspect="equal")

    display_labels = [_jaccard_display_label(m) for m in labels]
    ax.set_xlim(-0.5, n - 1.5)
    ax.set_ylim(n - 0.5, 0.5)
    ax.set_xticks(np.arange(n - 1))
    ax.set_yticks(np.arange(1, n))
    ax.set_xticklabels(display_labels[:-1], rotation=45, ha="right",
                       rotation_mode="anchor")
    ax.set_yticklabels(display_labels[1:] if show_ylabels else [])
    ax.tick_params(axis="both", which="major", direction="out", length=2.5,
                   width=0.45, labelsize=6, pad=1)
    if not show_ylabels:
        ax.tick_params(axis="y", which="major", left=False)
    ax.set_title(title, fontsize=7, pad=3, fontweight="normal")

    ax.set_xticks(np.arange(-0.5, n - 0.5, 1), minor=True)
    ax.set_yticks(np.arange(0.5, n + 0.5, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.65)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    midpoint = norm.vmin + 0.55 * (norm.vmax - norm.vmin)
    for i in range(n):
        for j in range(i):
            val = float(mat[i, j])
            color = "white" if val < midpoint else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=5.4, color=color)
    ax.set_facecolor("white")
    return im


# ── Plot functions ──────────────────────────────────────────────────────────

def plot_sj_jaccard(beds_by_method, output_dir, vmin=None):
    """Splice-junction Jaccard heatmap (single-column width)."""
    junc_sets = {}
    for method, isoforms in beds_by_method.items():
        juncs = set()
        for iso in isoforms:
            for j in iso["junctions"]:
                juncs.add((iso["chrom"], iso["strand"]) + j)
        junc_sets[method] = juncs

    methods = [m for m in beds_by_method if m in junc_sets]
    if len(methods) < 2:
        print("  Skipping SJ Jaccard: need ≥2 methods", file=sys.stderr)
        return None

    mat = _jaccard_matrix(junc_sets, methods)
    fig, ax = plt.subplots(figsize=(W1, W1 * 0.92))
    _draw_heatmap_ax(ax, mat, methods, cbar_label="Jaccard index", vmin=vmin)
    fig.tight_layout(pad=0.3)
    savefig(fig, Path(output_dir) / "sj_jaccard_heatmap.png")
    return mat


def _compute_end_sets(beds_by_method):
    """Extract TSS and TTS sets from parsed isoforms."""
    tss_map, tts_map = {}, {}
    for method, isoforms in beds_by_method.items():
        tss, tts = set(), set()
        for iso in isoforms:
            if iso["strand"] == "+":
                tss.add((iso["chrom"], "+", iso["start"] // END_BIN))
                tts.add((iso["chrom"], "+", iso["end"]   // END_BIN))
            else:
                tss.add((iso["chrom"], "-", iso["end"]   // END_BIN))
                tts.add((iso["chrom"], "-", iso["start"] // END_BIN))
        tss_map[method] = tss
        tts_map[method] = tts
    return tss_map, tts_map


def plot_end_jaccard(beds_by_method, output_dir, vmin=None):
    """TSS + TTS Jaccard heatmaps side-by-side (double-column width)."""
    tss_map, tts_map = _compute_end_sets(beds_by_method)

    methods = [m for m in beds_by_method if m in tss_map]
    if len(methods) < 2:
        print("  Skipping end Jaccard: need ≥2 methods", file=sys.stderr)
        return None, None

    tss_mat = _jaccard_matrix(tss_map, methods)
    tts_mat = _jaccard_matrix(tts_map, methods)

    fig, axes = plt.subplots(1, 2, figsize=(W2, W2 * 0.44))
    _draw_heatmap_ax(axes[0], tss_mat, methods, cbar_label="TSS Jaccard", vmin=vmin)
    _draw_heatmap_ax(axes[1], tts_mat, methods, cbar_label="TTS Jaccard", vmin=vmin)
    fig.tight_layout(w_pad=2, pad=0.3)
    savefig(fig, Path(output_dir) / "end_jaccard_heatmap.png")
    return tss_mat, tts_mat


def plot_jaccard_summary(beds_by_method, output_dir):
    """Three-panel lower-triangle SJ/TSS/TTS Jaccard figure."""
    methods = list(beds_by_method.keys())
    if len(methods) < 2:
        print("  Skipping Jaccard summary: need >=2 methods", file=sys.stderr)
        return None

    junc_sets = {}
    for m, isos in beds_by_method.items():
        junc_sets[m] = {
            (iso["chrom"], iso["strand"]) + j
            for iso in isos
            for j in iso["junctions"]
        }
    sj_mat = _jaccard_matrix(junc_sets, methods)
    tss_map, tts_map = _compute_end_sets(beds_by_method)
    tss_mat = _jaccard_matrix(tss_map, methods)
    tts_mat = _jaccard_matrix(tts_map, methods)

    vals = _offdiag_values(sj_mat, tss_mat, tts_mat)
    vmin = max(0.0, np.floor((float(vals.min()) - 0.01) * 20.0) / 20.0) if vals.size else 0.0
    norm = Normalize(vmin=vmin, vmax=1.0)

    fig, axes = plt.subplots(1, 3, figsize=(W2, W2 * 0.34))
    fig.subplots_adjust(left=0.11, right=0.89, bottom=0.31, top=0.84, wspace=0.16)
    for idx, (ax, mat, title) in enumerate(zip(
        axes,
        (sj_mat, tss_mat, tts_mat),
        ("Splice junctions", "TSS", "TTS"),
    )):
        _draw_jaccard_triangle(ax, mat, methods, title, norm=norm,
                               show_ylabels=(idx == 0))

    sm = ScalarMappable(norm=norm, cmap="cividis")
    sm.set_array([])
    cax = fig.add_axes([0.915, 0.33, 0.014, 0.49])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("Jaccard index", fontsize=7)
    cbar.ax.tick_params(direction="out", labelsize=6, length=2, width=0.4, pad=1)
    savefig(fig, Path(output_dir) / "jaccard_summary_heatmap.png")
    return sj_mat, tss_mat, tts_mat


def compute_global_vmin(beds_by_method):
    """Compute all three Jaccard matrices and return the global vmin.

    This ensures SJ, TSS, and TTS heatmaps share the same color scale,
    making visual comparisons meaningful.
    """
    methods = list(beds_by_method.keys())
    if len(methods) < 2:
        return 0.0

    # Splice junctions
    junc_sets = {}
    for method, isoforms in beds_by_method.items():
        juncs = set()
        for iso in isoforms:
            for j in iso["junctions"]:
                juncs.add((iso["chrom"], iso["strand"]) + j)
        junc_sets[method] = juncs
    sj_mat = _jaccard_matrix(junc_sets, methods)

    # TSS / TTS
    tss_map, tts_map = _compute_end_sets(beds_by_method)
    tss_mat = _jaccard_matrix(tss_map, methods)
    tts_mat = _jaccard_matrix(tts_map, methods)

    # Find global min across all off-diagonal values
    all_vals = []
    n = len(methods)
    for mat in (sj_mat, tss_mat, tts_mat):
        for i in range(n):
            for j in range(n):
                if i != j:
                    all_vals.append(mat[i, j])

    if not all_vals:
        return 0.0
    return max(0, min(all_vals) - 0.02)


# ── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--bed", nargs="+", required=True,
        help="label:path pairs, e.g. baseline:path/to/baseline.isoforms.bed",
    )
    parser.add_argument(
        "--output", required=True,
        help="Output directory for plots",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    beds_by_method = {}
    for entry in args.bed:
        if ":" not in entry:
            print(f"WARNING: skipping malformed entry '{entry}' (expected label:path)",
                  file=sys.stderr)
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

    if len(beds_by_method) < 2:
        print("Need ≥2 BED files for Jaccard comparison — skipping", file=sys.stderr)
        sys.exit(0)

    # Compute shared color scale across SJ, TSS, and TTS heatmaps
    global_vmin = compute_global_vmin(beds_by_method)
    if args.verbose:
        print(f"  Global Jaccard vmin={global_vmin:.3f}")

    plot_jaccard_summary(beds_by_method, output_dir)
    plot_sj_jaccard(beds_by_method, output_dir, vmin=global_vmin)
    plot_end_jaccard(beds_by_method, output_dir, vmin=global_vmin)
    print(f"Saved Jaccard heatmaps to {output_dir}")


if __name__ == "__main__":
    main()
