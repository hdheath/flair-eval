#!/usr/bin/env python3
"""
internal_priming.py — Internal priming analysis at isoform TTS positions.

For each isoform, computes the A-content in a window downstream of the TTS
using the reference genome sequence as a proxy for internal priming likelihood.
Also computes distance to the nearest annotated polyA site from a GTF.

Produces three plots:
  A1. internal_priming_summary.png
      Fraction of isoforms with high A-content (>= threshold) at TTS,
      stratified by read-support bins. Per-method bars, all methods on same axes.

  A2. internal_priming_cross_tool.png
      Same A-content metric across ALL tools (FLAIR, FLAIR-permissive,
      IsoQuant, TED modes) — single grouped bar chart for direct comparison.

  A3. proximal_apa_distance_support.png
      Scatter: distance from TTS to nearest annotated polyA site (y-axis)
      vs read support (x-axis), coloured by within-50bp (likely real) vs
      distant (likely internal priming). Faceted by method.

Usage:
    python internal_priming.py \\
        --bed label1:bed1.bed label2:bed2.bed ... \\
        --read-map label1:map1.txt label2:map2.txt ... \\
        --genome genome.fa \\
        --gtf annotation.gtf \\
        --output output_dir/ \\
        [--window 30] \\
        [--a-threshold 0.60] \\
        [--apa-window 50]
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
    from signal_utils import parse_isoforms, tss_tts, load_read_map, lookup_read_count
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
    from evaluation.signal_utils import parse_isoforms, tss_tts, load_read_map, lookup_read_count

apply_rc()

# ── Constants ────────────────────────────────────────────────────────────────

DEFAULT_WINDOW    = 30    # bp downstream of TTS to scan for A-content
DEFAULT_A_THRESH  = 0.60  # A-fraction threshold for "likely internal priming"
DEFAULT_APA_WIN   = 50    # bp window for "proximal to polyA site"
SUPPORT_BINS      = [1, 2, 5, 10, 20, float("inf")]
SUPPORT_LABELS    = ["1", "2–4", "5–9", "10–19", "20+"]


# ── Genome FASTA ─────────────────────────────────────────────────────────────

def load_fasta(path: str | Path) -> Dict[str, str]:
    """Load a FASTA file into {chrom: sequence}."""
    seqs: Dict[str, str] = {}
    chrom = None
    buf: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if chrom is not None:
                    seqs[chrom] = "".join(buf).upper()
                chrom = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if chrom is not None:
        seqs[chrom] = "".join(buf).upper()
    return seqs


def a_content_downstream(chrom: str, tts: int, strand: str,
                          genome: Dict[str, str], window: int) -> Optional[float]:
    """Return fraction of A bases in a window downstream of the TTS.

    Downstream means: 3' of the TTS in the direction of transcription.
      + strand: sequence [tts, tts+window)
      - strand: sequence [tts-window, tts), then reverse-complement → A becomes T
    """
    seq = genome.get(chrom)
    if seq is None:
        return None
    if strand == "+":
        region = seq[max(0, tts):min(len(seq), tts + window)]
        return region.count("A") / len(region) if region else None
    else:
        region = seq[max(0, tts - window):min(len(seq), tts)]
        # downstream on minus strand → look for T (complement of A)
        return region.count("T") / len(region) if region else None


# ── GTF polyA site extraction ─────────────────────────────────────────────────

def extract_polya_sites(gtf_path: str | Path) -> Dict[str, List[int]]:
    """Extract annotated TTS positions per chrom from a GTF as polyA site proxies."""
    sites: Dict[str, List[int]] = defaultdict(list)
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip().split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            chrom  = cols[0]
            start  = int(cols[3]) - 1
            end    = int(cols[4])
            strand = cols[6]
            _, tts = tss_tts(start, end, strand)
            sites[chrom].append(tts)
    # sort each list
    for ch in sites:
        sites[ch] = sorted(set(sites[ch]))
    return sites


def nearest_apa_distance(chrom: str, tts: int,
                          polya_sites: Dict[str, List[int]]) -> Optional[int]:
    """Return distance (bp) from tts to nearest annotated polyA site."""
    import bisect
    positions = polya_sites.get(chrom)
    if not positions:
        return None
    idx = bisect.bisect_left(positions, tts)
    candidates = []
    if idx < len(positions):
        candidates.append(abs(positions[idx] - tts))
    if idx > 0:
        candidates.append(abs(positions[idx - 1] - tts))
    return min(candidates) if candidates else None


# ── Support bin assignment ───────────────────────────────────────────────────

def support_bin(n: int) -> int:
    for i, upper in enumerate(SUPPORT_BINS[1:]):
        if n < upper:
            return i
    return len(SUPPORT_BINS) - 2


# ── Per-isoform data building ─────────────────────────────────────────────────

def build_isoform_data(
    beds_by_method: Dict[str, List[dict]],
    read_maps:       Dict[str, Dict[str, int]],
    genome:          Dict[str, str],
    polya_sites:     Dict[str, List[int]],
    window:          int,
    a_thresh:        float,
) -> Dict[str, List[dict]]:
    """Return {method: [{reads, a_content, apa_dist, likely_ip}]}."""
    data: Dict[str, List[dict]] = {}
    for method, isoforms in beds_by_method.items():
        rc = read_maps.get(method, {})
        rows = []
        for iso in isoforms:
            reads  = _lookup_count(iso["name"], rc) if rc else 0
            a_frac = a_content_downstream(
                iso["chrom"], iso["tts"], iso["strand"], genome, window)
            apa_d  = nearest_apa_distance(iso["chrom"], iso["tts"], polya_sites)
            if a_frac is None:
                continue
            rows.append({
                "reads":     reads,
                "a_content": a_frac,
                "apa_dist":  apa_d,
                "likely_ip": a_frac >= a_thresh,
                "bin":       support_bin(reads),
            })
        data[method] = rows
    return data


# ── Plot A1: internal priming fraction by read-support bin ───────────────────

def plot_internal_priming_summary(
    data: Dict[str, List[dict]],
    output_dir: Path,
    a_thresh: float,
):
    methods = list(data.keys())
    if not methods:
        return
    styler = ModeStyler(methods)
    n_bins = len(SUPPORT_LABELS)
    x      = np.arange(n_bins)
    width  = 0.8 / max(len(methods), 1)

    fig, ax = plt.subplots(figsize=(W2, W1 * 0.75))
    handles = []
    for i, m in enumerate(methods):
        rows = data[m]
        fracs = []
        for b in range(n_bins):
            bin_rows = [r for r in rows if r["bin"] == b]
            if bin_rows:
                fracs.append(sum(r["likely_ip"] for r in bin_rows) / len(bin_rows))
            else:
                fracs.append(0.0)
        offset = (i - len(methods) / 2 + 0.5) * width
        bars = ax.bar(x + offset, fracs, width * 0.9,
                      color=styler.color(m), edgecolor="none", alpha=0.85)
        handles.append(styler.legend_handle(m, label=m, markersize=5))

    ax.set_xticks(x)
    ax.set_xticklabels(SUPPORT_LABELS, fontsize=7)
    ax.set_ylim(0, 1)
    style_ax(ax,
             xlabel="Read support (reads per isoform)",
             ylabel=f"Fraction with TTS A-content ≥ {a_thresh:.0%}")
    ax.axhline(a_thresh, color="0.5", lw=0.5, ls="--", zorder=0)
    legend_outside(fig, handles=handles, loc="upper right",
                   bbox_to_anchor=(1.0, 1.0), ncol=1, fontsize=6)
    fig.tight_layout(pad=0.3)
    savefig(fig, output_dir / "internal_priming_summary.png")


# ── Plot A2: cross-tool single grouped bar ───────────────────────────────────

def plot_cross_tool_comparison(
    data: Dict[str, List[dict]],
    output_dir: Path,
    a_thresh: float,
):
    methods = list(data.keys())
    if not methods:
        return
    styler  = ModeStyler(methods)
    fracs   = [
        (sum(r["likely_ip"] for r in data[m]) / len(data[m]))
        if data[m] else 0.0
        for m in methods
    ]
    x = np.arange(len(methods))

    fig, ax = plt.subplots(figsize=(W2, W1 * 0.65))
    colors  = [styler.color(m) for m in methods]
    ax.bar(x, fracs, 0.7, color=colors, edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=5, rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.axhline(a_thresh, color="0.5", lw=0.5, ls="--", zorder=0)
    style_ax(ax, ylabel=f"Fraction with TTS A-content ≥ {a_thresh:.0%}")
    fig.tight_layout(pad=0.3)
    savefig(fig, output_dir / "internal_priming_cross_tool.png")


# ── Plot A3: APA distance vs read support ────────────────────────────────────

def plot_proximal_apa(
    data: Dict[str, List[dict]],
    output_dir: Path,
    apa_window: int,
):
    methods = [m for m in data if any(r["apa_dist"] is not None for r in data[m])]
    if not methods:
        return

    styler = ModeStyler(methods)
    n = len(methods)
    ncols = min(n, 4)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(W2, W1 * 0.75 * nrows / max(nrows, 1)),
        squeeze=False,
    )

    for idx, m in enumerate(methods):
        ax   = axes[idx // ncols][idx % ncols]
        rows = [r for r in data[m] if r["apa_dist"] is not None and r["reads"] > 0]
        if not rows:
            ax.set_visible(False)
            continue

        xs       = np.array([r["reads"]   for r in rows], dtype=float)
        ys       = np.array([r["apa_dist"] for r in rows], dtype=float)
        proximal = ys <= apa_window

        ax.scatter(xs[~proximal], ys[~proximal],
                   s=1.5, alpha=0.25, color="#CC79A7",
                   edgecolors="none", rasterized=True, label="Distant")
        ax.scatter(xs[proximal], ys[proximal],
                   s=1.5, alpha=0.4, color="#009E73",
                   edgecolors="none", rasterized=True, label="Proximal")
        ax.axhline(apa_window, color="0.5", lw=0.5, ls="--")
        ax.set_xscale("log")
        style_ax(ax)
        ax.text(0.04, 0.97, m, transform=ax.transAxes,
                ha="left", va="top", fontsize=5, fontweight="bold")
        pct = proximal.sum() / len(proximal) * 100 if len(proximal) else 0
        ax.text(0.04, 0.88, f"{pct:.0f}% proximal",
                transform=ax.transAxes, ha="left", va="top", fontsize=5,
                color="#009E73")

        if idx >= n - ncols:
            ax.set_xlabel("Read support", fontsize=6)
        if idx % ncols == 0:
            ax.set_ylabel("Distance to nearest polyA site (bp)", fontsize=6)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    # Shared legend
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#009E73",
               markersize=5, label=f"≤ {apa_window} bp (proximal)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#CC79A7",
               markersize=5, label=f"> {apa_window} bp (likely IP)"),
    ]
    fig.legend(handles=handles, loc="lower right", fontsize=6, frameon=False,
               bbox_to_anchor=(1.0, 0.0))
    fig.tight_layout(pad=0.3)
    savefig(fig, output_dir / "proximal_apa_distance_support.png")


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bed",       nargs="+", required=True,
                        help="label:path pairs for BED12 isoform files")
    parser.add_argument("--read-map",  nargs="+", default=[],
                        help="label:path pairs for isoform read-map files")
    parser.add_argument("--genome",    required=True,
                        help="Reference genome FASTA")
    parser.add_argument("--gtf",       required=True,
                        help="GTF annotation (for polyA site proxies)")
    parser.add_argument("--output",    required=True,
                        help="Output directory")
    parser.add_argument("--window",    type=int, default=DEFAULT_WINDOW,
                        help=f"bp downstream of TTS for A-content scan (default: {DEFAULT_WINDOW})")
    parser.add_argument("--a-threshold", type=float, default=DEFAULT_A_THRESH,
                        help=f"A-fraction threshold for internal priming call (default: {DEFAULT_A_THRESH})")
    parser.add_argument("--apa-window", type=int, default=DEFAULT_APA_WIN,
                        help=f"bp for proximal polyA classification (default: {DEFAULT_APA_WIN})")
    parser.add_argument("--verbose",   action="store_true")
    args = parser.parse_args()

    # --- Load BED files ---
    beds_by_method: Dict[str, List[dict]] = {}
    for entry in args.bed:
        if ":" not in entry:
            print(f"WARNING: skipping malformed --bed entry '{entry}'", file=sys.stderr)
            continue
        label, path = entry.split(":", 1)
        if not Path(path).exists():
            print(f"WARNING: file not found: {path}", file=sys.stderr)
            continue
        isos = parse_isoforms(path)
        if isos:
            beds_by_method[label] = isos
            if args.verbose:
                print(f"  {label}: {len(isos)} isoforms", file=sys.stderr)

    if not beds_by_method:
        print("No isoform data loaded — skipping", file=sys.stderr)
        sys.exit(1)

    # --- Load read maps ---
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

    # --- Load genome ---
    if args.verbose:
        print("  Loading genome FASTA...", file=sys.stderr)
    genome = load_fasta(args.genome)

    # --- Load polyA sites from GTF ---
    if args.verbose:
        print("  Extracting polyA sites from GTF...", file=sys.stderr)
    polya_sites = extract_polya_sites(args.gtf)

    # --- Build per-isoform data ---
    if args.verbose:
        print("  Computing A-content and APA distances...", file=sys.stderr)
    data = build_isoform_data(
        beds_by_method, read_maps, genome, polya_sites,
        args.window, args.a_threshold,
    )

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_internal_priming_summary(data, output_dir, args.a_threshold)
    plot_cross_tool_comparison(data, output_dir, args.a_threshold)
    plot_proximal_apa(data, output_dir, args.apa_window)

    print(f"Saved internal priming plots to {args.output}")


if __name__ == "__main__":
    main()
