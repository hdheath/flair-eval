#!/usr/bin/env python3
"""
Plot histograms of pairwise end distances within splice junction chain (SJC) groups.

For each SJC that has >1 isoform (alternative promoters or polyadenylation sites),
compute pairwise distances between 5' ends (TSS) and 3' ends (TTS), then show
histograms comparing modes.

Usage:
    python sjc_end_distance_plot.py \
        --bed-dir results/.../transcriptome \
        --samples A549_cDNA A549_dRNA A549_ONT_cDNA WTC11_cDNA \
        --modes ted-2d ted-1d2d \
        --test-name 04_02_26_4samples_ted_comparison_v2 \
        --output sjc_end_distances.png
"""

import argparse
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    from pub_style import style_ax, savefig, MODE_COLORS
    from signal_utils import parse_bed12, group_by_junction_chain
except ImportError:
    from evaluation.pub_style import style_ax, savefig, MODE_COLORS
    from evaluation.signal_utils import parse_bed12, group_by_junction_chain


def compute_pairwise_distances(groups):
    """For SJC groups with >1 isoform, compute pairwise TSS and TTS distances.

    Args:
        groups: dict of sjc_key -> list of isoform dicts (from group_by_junction_chain).

    Returns (tss_distances, tts_distances) as lists of absolute distances.
    """
    tss_dists = []
    tts_dists = []
    for sjc_key, isoforms in groups.items():
        if len(isoforms) < 2:
            continue
        for a, b in combinations(isoforms, 2):
            tss_dists.append(abs(a["tss"] - b["tss"]))
            tts_dists.append(abs(a["tts"] - b["tts"]))
    return tss_dists, tts_dists


def plot_end_distance_histograms(mode_data, output_path, max_dist=2000, bin_width=50):
    """Create histogram figure with TSS and TTS distance distributions.

    mode_data: dict mode -> (tss_dists, tts_dists, n_groups, n_isoforms)
    """
    modes = sorted(mode_data.keys(),
                   key=lambda m: ['default', 'ted-2d', 'ted-1d2d'].index(m)
                   if m in ['default', 'ted-2d', 'ted-1d2d'] else 999)

    bins = np.arange(0, max_dist + bin_width, bin_width)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=False)

    for ax, end_type, end_idx, label in [
        (axes[0], 'tss', 0, "5\u2032 End (TSS) Distance"),
        (axes[1], 'tts', 1, "3\u2032 End (TTS) Distance"),
    ]:
        for mode in modes:
            dists = mode_data[mode][end_idx]
            if not dists:
                continue
            color = MODE_COLORS.get(mode, '#999999')
            n_groups = mode_data[mode][2]
            n_iso = mode_data[mode][3]
            clipped = [d for d in dists if d <= max_dist]
            n_beyond = len(dists) - len(clipped)
            lbl = f"{mode} (n={len(dists)}, {n_groups} SJCs)"
            ax.hist(clipped, bins=bins, alpha=0.55, color=color, edgecolor='white',
                    linewidth=0.3, label=lbl)

        ax.set_xlabel('Pairwise Distance (bp)', fontsize=9)
        ax.set_ylabel('Count', fontsize=9)
        ax.set_title(label, fontsize=10)
        style_ax(ax)
        ax.legend(fontsize=7, frameon=False)

    fig.suptitle('Distance Between Alternative Ends Sharing the Same SJC',
                 fontsize=11, y=1.02)
    fig.tight_layout()
    savefig(fig, Path(output_path), dpi=300)
    print(f"Saved SJC end distance plot to {output_path}")


def plot_end_distance_2d(mode_data, output_path, max_dist=2000):
    """2D scatter of TSS distance vs TTS distance per isoform pair.

    Points near x-axis → alternative promoter only (same TTS).
    Points near y-axis → alternative polyadenylation only (same TSS).
    Points off both axes → both ends differ.

    mode_data: dict mode -> (tss_dists, tts_dists, n_groups, n_isoforms)
    """
    modes = sorted(mode_data.keys(),
                   key=lambda m: ['default', 'ted-2d', 'ted-1d2d'].index(m)
                   if m in ['default', 'ted-2d', 'ted-1d2d'] else 999)
    # Skip modes with no data
    modes = [m for m in modes if len(mode_data[m][0]) > 0]
    if not modes:
        return

    n_modes = len(modes)
    fig, axes = plt.subplots(1, n_modes, figsize=(4.5 * n_modes, 4.2), squeeze=False)

    for col, mode in enumerate(modes):
        ax = axes[0][col]
        tss_d = np.array(mode_data[mode][0])
        tts_d = np.array(mode_data[mode][1])
        n_groups = mode_data[mode][2]
        n_pairs = len(tss_d)

        color = MODE_COLORS.get(mode, '#999999')

        # Classify pairs
        thresh = 50  # bp — "same end" threshold
        both_same = np.sum((tss_d <= thresh) & (tts_d <= thresh))
        tss_only = np.sum((tss_d > thresh) & (tts_d <= thresh))
        tts_only = np.sum((tss_d <= thresh) & (tts_d > thresh))
        both_diff = np.sum((tss_d > thresh) & (tts_d > thresh))

        # Clip for display
        tss_show = np.clip(tss_d, 0, max_dist)
        tts_show = np.clip(tts_d, 0, max_dist)

        ax.scatter(tss_show, tts_show, s=18, alpha=0.5, color=color,
                   edgecolors='white', linewidths=0.3, zorder=3)

        # Quadrant annotations
        txt_kw = dict(fontsize=6.5, color='#555555', ha='center',
                      bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))
        ax.text(max_dist * 0.75, max_dist * 0.06,
                f'Alt promoter only: {tss_only}', **txt_kw)
        ax.text(max_dist * 0.25, max_dist * 0.94,
                f'Alt polyA only: {tts_only}', **txt_kw)
        ax.text(max_dist * 0.75, max_dist * 0.94,
                f'Both differ: {both_diff}', **txt_kw)
        ax.text(max_dist * 0.25, max_dist * 0.06,
                f'Same ends: {both_same}', **txt_kw)

        # Reference lines at threshold
        ax.axhline(thresh, color='#aaaaaa', linewidth=0.6, linestyle='--', zorder=1)
        ax.axvline(thresh, color='#aaaaaa', linewidth=0.6, linestyle='--', zorder=1)

        ax.set_xlim(-max_dist * 0.03, max_dist * 1.05)
        ax.set_ylim(-max_dist * 0.03, max_dist * 1.05)
        ax.set_xlabel("5\u2032 (TSS) pairwise distance (bp)", fontsize=9)
        if col == 0:
            ax.set_ylabel("3\u2032 (TTS) pairwise distance (bp)", fontsize=9)
        ax.set_title(f"{mode}  ({n_pairs} pairs, {n_groups} SJCs)", fontsize=9)
        ax.set_aspect('equal')
        style_ax(ax)

    fig.suptitle('Alternative Ends Within Shared SJC Groups', fontsize=11, y=1.02)
    fig.tight_layout()
    savefig(fig, Path(output_path), dpi=300)
    print(f"Saved 2D end distance plot to {output_path}")


def plot_end_distance_2d_per_sample(sample_mode_data, output_path, max_dist=2000):
    """Per-sample faceted 2D scatter of TSS vs TTS distances."""
    samples = sorted(sample_mode_data.keys())
    # Collect modes that have data in at least one sample
    all_modes = []
    for s in samples:
        for m in sample_mode_data[s]:
            if m not in all_modes and len(sample_mode_data[s][m][0]) > 0:
                all_modes.append(m)
    mode_order = ['default', 'ted-2d', 'ted-1d2d']
    all_modes.sort(key=lambda m: mode_order.index(m) if m in mode_order else 999)
    if not all_modes:
        return

    n_modes = len(all_modes)
    n_samples = len(samples)

    fig, axes = plt.subplots(n_samples, n_modes,
                              figsize=(4.0 * n_modes, 3.8 * n_samples),
                              squeeze=False)

    thresh = 50
    for row, sample in enumerate(samples):
        for col, mode in enumerate(all_modes):
            ax = axes[row][col]
            data = sample_mode_data[sample].get(mode)
            if data is None or len(data[0]) == 0:
                ax.text(0.5, 0.5, 'no data', transform=ax.transAxes,
                        ha='center', va='center', fontsize=8, color='#999999')
                ax.set_title(f"{sample} — {mode}", fontsize=8)
                style_ax(ax)
                continue

            tss_d = np.array(data[0])
            tts_d = np.array(data[1])
            n_groups = data[2]
            color = MODE_COLORS.get(mode, '#999999')

            tss_only = np.sum((tss_d > thresh) & (tts_d <= thresh))
            tts_only = np.sum((tss_d <= thresh) & (tts_d > thresh))
            both_diff = np.sum((tss_d > thresh) & (tts_d > thresh))

            tss_show = np.clip(tss_d, 0, max_dist)
            tts_show = np.clip(tts_d, 0, max_dist)

            ax.scatter(tss_show, tts_show, s=14, alpha=0.5, color=color,
                       edgecolors='white', linewidths=0.2, zorder=3)
            ax.axhline(thresh, color='#aaaaaa', linewidth=0.5, linestyle='--', zorder=1)
            ax.axvline(thresh, color='#aaaaaa', linewidth=0.5, linestyle='--', zorder=1)

            txt_kw = dict(fontsize=5.5, color='#555555', ha='center',
                          bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0.5))
            ax.text(max_dist * 0.75, max_dist * 0.06,
                    f'Prom: {tss_only}', **txt_kw)
            ax.text(max_dist * 0.25, max_dist * 0.94,
                    f'PolyA: {tts_only}', **txt_kw)
            ax.text(max_dist * 0.75, max_dist * 0.94,
                    f'Both: {both_diff}', **txt_kw)

            ax.set_xlim(-max_dist * 0.03, max_dist * 1.05)
            ax.set_ylim(-max_dist * 0.03, max_dist * 1.05)
            ax.set_aspect('equal')
            if row == n_samples - 1:
                ax.set_xlabel("5\u2032 distance (bp)", fontsize=7)
            if col == 0:
                ax.set_ylabel("3\u2032 distance (bp)", fontsize=7)
            ax.set_title(f"{sample} — {mode} ({len(tss_d)}p, {n_groups} SJCs)",
                         fontsize=7)
            style_ax(ax)

    fig.suptitle('Alternative Ends Within Shared SJC Groups', fontsize=11, y=1.01)
    fig.tight_layout()
    savefig(fig, Path(output_path), dpi=300)
    print(f"Saved per-sample 2D end distance plot to {output_path}")


def plot_per_sample_end_distances(sample_mode_data, output_path, max_dist=2000, bin_width=50):
    """Create per-sample faceted histogram figure.

    sample_mode_data: dict sample -> mode -> (tss_dists, tts_dists, n_groups, n_isoforms)
    """
    samples = sorted(sample_mode_data.keys())
    n_samples = len(samples)
    if n_samples == 0:
        return

    fig, axes = plt.subplots(n_samples, 2, figsize=(10, 2.8 * n_samples),
                              squeeze=False, sharey=False)

    bins = np.arange(0, max_dist + bin_width, bin_width)
    mode_order = ['default', 'ted-2d', 'ted-1d2d']

    for row, sample in enumerate(samples):
        modes_data = sample_mode_data[sample]
        modes = sorted(modes_data.keys(),
                       key=lambda m: mode_order.index(m) if m in mode_order else 999)

        for col, (end_idx, label) in enumerate([(0, "5\u2032 (TSS)"), (1, "3\u2032 (TTS)")]):
            ax = axes[row][col]
            for mode in modes:
                dists = modes_data[mode][end_idx]
                if not dists:
                    continue
                color = MODE_COLORS.get(mode, '#999999')
                n_groups = modes_data[mode][2]
                clipped = [d for d in dists if d <= max_dist]
                lbl = f"{mode} (n={len(dists)}, {n_groups} SJCs)"
                ax.hist(clipped, bins=bins, alpha=0.55, color=color, edgecolor='white',
                        linewidth=0.3, label=lbl)

            ax.set_xlabel('Pairwise Distance (bp)', fontsize=8)
            ax.set_ylabel('Count', fontsize=8)
            title = f"{sample} — {label}"
            ax.set_title(title, fontsize=9)
            style_ax(ax)
            ax.legend(fontsize=6, frameon=False)

    fig.suptitle('Distance Between Alternative Ends Sharing the Same SJC',
                 fontsize=11, y=1.01)
    fig.tight_layout()
    savefig(fig, Path(output_path), dpi=300)
    print(f"Saved per-sample SJC end distance plot to {output_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--bed', nargs='+', required=True,
                        help='BED12 isoform files as label:path pairs. '
                             'Label format: sample::mode (e.g. A549_cDNA::TED-default)')
    parser.add_argument('--output', required=True,
                        help='Output directory for all plots')
    parser.add_argument('--max-dist', type=int, default=2000,
                        help='Maximum distance to show in histogram (bp)')
    parser.add_argument('--bin-width', type=int, default=50,
                        help='Histogram bin width (bp)')

    args = parser.parse_args()
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    # Parse label:path pairs → {sample: {mode: (tss, tts, n_groups, n_iso)}}
    combined_mode_data = defaultdict(lambda: [[], [], 0, 0])
    sample_mode_data = {}

    for entry in args.bed:
        if ':' not in entry:
            print(f"Warning: skipping malformed entry (no ':'): {entry}", file=sys.stderr)
            continue
        # Split on last colon so label can contain :: separator
        label, bed_path_str = entry.rsplit(':', 1)
        bed_path = Path(bed_path_str)

        # Label is "sample::mode" or just "mode"
        if '::' in label:
            sample, mode = label.split('::', 1)
        else:
            sample, mode = 'all', label

        if not bed_path.exists():
            print(f"Warning: {bed_path} not found, skipping", file=sys.stderr)
            continue

        groups = group_by_junction_chain(parse_bed12(bed_path))
        multi_groups = {k: v for k, v in groups.items() if len(v) > 1}
        tss_dists, tts_dists = compute_pairwise_distances(groups)
        n_groups = len(multi_groups)
        n_iso = sum(len(v) for v in multi_groups.values())

        print(f"  {sample} {mode}: {len(groups)} SJC groups, "
              f"{n_groups} with alt ends ({n_iso} isoforms), "
              f"{len(tss_dists)} TSS pairs, {len(tts_dists)} TTS pairs")

        if sample not in sample_mode_data:
            sample_mode_data[sample] = {}
        sample_mode_data[sample][mode] = (tss_dists, tts_dists, n_groups, n_iso)

        # Accumulate into combined
        combined_mode_data[mode][0].extend(tss_dists)
        combined_mode_data[mode][1].extend(tts_dists)
        combined_mode_data[mode][2] += n_groups
        combined_mode_data[mode][3] += n_iso

    # Convert combined to tuples
    combined = {mode: (data[0], data[1], data[2], data[3])
                for mode, data in combined_mode_data.items()}

    has_data = any(len(d[0]) > 0 for d in combined.values())

    if not has_data:
        print("No multi-isoform SJC groups found in any mode", file=sys.stderr)
        return

    # Per-sample 2D scatter (the primary output)
    plot_end_distance_2d_per_sample(sample_mode_data,
                                     outdir / 'sjc_alt_end_2d_per_sample.png',
                                     max_dist=args.max_dist)


if __name__ == '__main__':
    main()
