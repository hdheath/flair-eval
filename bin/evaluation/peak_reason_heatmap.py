#!/usr/bin/env python3
"""
Cross-mode peak reason heatmap.

Reads per-peak reason TSVs from multiple evaluation runs and produces a
categorical heatmap where:
  - Each row is a peak (sorted by signal strength)
  - Each column is a mode/version
  - Each cell is colored by the recovery reason

One figure per end type (TSS/CAGE, TTS/QuantSeq).

Usage:
    python peak_reason_heatmap.py \\
        --cage-tsvs  mode1_cage_peak_reasons.tsv mode2_cage_peak_reasons.tsv ... \\
        --quantseq-tsvs mode1_quantseq_peak_reasons.tsv mode2_quantseq_peak_reasons.tsv ... \\
        --output-prefix my_test_heatmap \\
        [--max-peaks 300] [--title-prefix "test: "] [--verbose]
"""

import argparse
import csv
import re
import sys
import warnings
from pathlib import Path
from typing import Dict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pub_style import REASON_COLORS as _PUB_REASON_COLORS, savefig


# ── Reason palette (kept in sync with plots.py) ─────────────────────────────

_REASON_COLORS = _PUB_REASON_COLORS

_REASON_LABELS = {
    'recovered':           'Recovered',
    'no_reads':            'No long reads',
    'single_exon_only':    'Single-exon only',
    'near_miss':           'Near miss',
    'reads_unassigned':    'Reads unassigned',
    'reads_redirected':    'Reads redirected',
    'proximal_apa':        'Proximal APA',
    'trailing_truncation': 'Trailing / truncation',
    'other_missed':        'Other missed',
}

_REASON_ORDER = [
    'recovered', 'near_miss', 'reads_redirected', 'proximal_apa',
    'trailing_truncation', 'single_exon_only', 'reads_unassigned',
    'no_reads', 'other_missed',
]

PLOT_DPI = 1000


# ── Helpers ──────────────────────────────────────────────────────────────────

def _save_figure(fig, output_path: Path) -> bool:
    """Save figure and close."""
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        savefig(fig, output_path, dpi=PLOT_DPI)
        return True
    except Exception as e:
        print(f"ERROR: Failed to save plot to {output_path}: {e}", file=sys.stderr)
        plt.close(fig)
        return False


def _reason_to_int(reason: str) -> int:
    """Map reason string to integer for colormap indexing."""
    _MAP = {r: i for i, r in enumerate(_REASON_ORDER)}
    return _MAP.get(reason, len(_REASON_ORDER) - 1)


# ── Main heatmap function ───────────────────────────────────────────────────

def _rank_array(values):
    """Return 0-based dense ranks (highest value = rank 0)."""
    indexed = sorted(enumerate(values), key=lambda x: -x[1])
    ranks = [0] * len(values)
    for rank, (orig_idx, _) in enumerate(indexed):
        ranks[orig_idx] = rank
    return ranks


def plot_peak_reason_heatmap(
    mode_peak_reasons: Dict[str, Dict[str, str]],
    mode_peak_scores: Dict[str, Dict[str, float]],
    output_path: Path,
    title: str = "Peak Recovery Heatmap",
    end_type: str = 'tss',
    max_peaks: int = 300,
    mode_peak_reads: Dict[str, Dict[str, float]] = None,
) -> bool:
    """Categorical heatmap: rows=peaks, columns=modes, cells colored by reason.

    Peaks are sorted by combined rank of orthogonal signal strength and
    long-read support (highest combined rank at top).  Two sidebars on the
    left show relative orthogonal signal and long-read support.

    Args:
        mode_peak_reasons: Dict[mode_name -> Dict[peak_id -> reason_string]]
        mode_peak_scores:  Dict[mode_name -> Dict[peak_id -> signal_score]]
        output_path: Where to save the PNG
        title: Plot suptitle
        end_type: 'tss' or 'tts'
        max_peaks: Maximum number of peaks to show (highest combined rank first)
        mode_peak_reads:   Dict[mode_name -> Dict[peak_id -> read_count]]
    """
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch

    if mode_peak_reads is None:
        mode_peak_reads = {}

    if not mode_peak_reasons:
        return False

    # Collect all peak IDs across all modes
    all_peak_ids = set()
    for reasons in mode_peak_reasons.values():
        all_peak_ids.update(reasons.keys())

    if not all_peak_ids:
        return False

    all_peak_list = list(all_peak_ids)

    # Compute mean signal per peak across modes
    peak_mean_signal = {}
    for pid in all_peak_list:
        scores = []
        for mode, score_map in mode_peak_scores.items():
            if pid in score_map:
                scores.append(score_map[pid])
        peak_mean_signal[pid] = (sum(scores) / len(scores)) if scores else 0.0

    # Compute mean long-read support per peak across modes
    peak_mean_reads = {}
    for pid in all_peak_list:
        counts = []
        for mode, read_map in mode_peak_reads.items():
            if pid in read_map:
                counts.append(read_map[pid])
        peak_mean_reads[pid] = (sum(counts) / len(counts)) if counts else 0.0

    # Sort by combined rank: rank by signal + rank by read support
    signal_vals_all = [peak_mean_signal[p] for p in all_peak_list]
    read_vals_all = [peak_mean_reads[p] for p in all_peak_list]
    signal_ranks = _rank_array(signal_vals_all)
    read_ranks = _rank_array(read_vals_all)
    combined = [sr + rr for sr, rr in zip(signal_ranks, read_ranks)]
    order = sorted(range(len(all_peak_list)), key=lambda i: combined[i])
    sorted_peaks = [all_peak_list[i] for i in order]
    if len(sorted_peaks) > max_peaks:
        sorted_peaks = sorted_peaks[:max_peaks]

    # Order modes consistently
    _MODE_ORDER = [
        "default", "density-asymmetric", "density-asymmetric-softclip",
        "density-plain", "density-strict", "k-means", "more-ends",
        "bambu_default", "isoquant_pacbio",
    ]
    modes = sorted(mode_peak_reasons.keys(),
                   key=lambda m: _MODE_ORDER.index(m) if m in _MODE_ORDER else 999)

    n_peaks = len(sorted_peaks)
    n_modes = len(modes)

    # Build the integer matrix (rows=peaks, cols=modes)
    n_reasons = len(_REASON_ORDER)
    ABSENT = n_reasons  # extra color for "peak not in this mode's evaluation"

    matrix = np.full((n_peaks, n_modes), ABSENT, dtype=int)
    for j, mode in enumerate(modes):
        reasons = mode_peak_reasons.get(mode, {})
        for i, pid in enumerate(sorted_peaks):
            if pid in reasons:
                matrix[i, j] = _reason_to_int(reasons[pid])

    # Build colormap: reason colors + absent color
    cmap_colors = [_REASON_COLORS[r] for r in _REASON_ORDER] + ['#f0f0f0']
    cmap = ListedColormap(cmap_colors)
    bounds = list(range(n_reasons + 2))
    norm = BoundaryNorm(bounds, cmap.N)

    # Figure sizing
    cell_w = 0.55
    cell_h = 0.12
    fig_w = max(6, 2.4 + cell_w * n_modes + 3.5)
    fig_h = max(4, 1.5 + cell_h * n_peaks)
    fig_h = min(fig_h, 50)

    has_reads = bool(mode_peak_reads)
    n_sidebars = 2 if has_reads else 1
    width_ratios = [0.08, 0.08, 1] if has_reads else [0.12, 1]
    fig, axes = plt.subplots(
        1, n_sidebars + 1, figsize=(fig_w, fig_h),
        gridspec_kw={'width_ratios': width_ratios, 'wspace': 0.02},
    )
    if n_sidebars == 1:
        axes = [axes[0], axes[1]]  # keep indexing consistent

    # ── Left sidebar 1: orthogonal signal strength ──
    ax_signal = axes[0]
    signal_vals = np.array([peak_mean_signal.get(pid, 0) for pid in sorted_peaks])
    max_sig = signal_vals.max() if signal_vals.max() > 0 else 1.0
    signal_norm = signal_vals / max_sig
    signal_img = signal_norm.reshape(-1, 1)
    ax_signal.imshow(signal_img, aspect='auto', cmap='Blues', vmin=0, vmax=1,
                     interpolation='nearest')
    ax_signal.set_xticks([0])
    sig_label = "5' Peaks" if end_type == 'tss' else "3' Peaks"
    ax_signal.set_xticklabels([sig_label], fontsize=7, rotation=90)
    ax_signal.set_yticks([])
    sort_desc = 'signal + read support' if has_reads else 'signal'
    ax_signal.set_ylabel(f'Peaks (n={n_peaks}, sorted by {sort_desc})', fontsize=7)

    # ── Left sidebar 2: long-read support ──
    if has_reads:
        ax_reads = axes[1]
        read_vals = np.array([peak_mean_reads.get(pid, 0) for pid in sorted_peaks])
        max_rd = read_vals.max() if read_vals.max() > 0 else 1.0
        read_norm = read_vals / max_rd
        read_img = read_norm.reshape(-1, 1)
        ax_reads.imshow(read_img, aspect='auto', cmap='Oranges', vmin=0, vmax=1,
                        interpolation='nearest')
        ax_reads.set_xticks([0])
        ax_reads.set_xticklabels(['Reads'], fontsize=7, rotation=90)
        ax_reads.set_yticks([])

    # ── Main heatmap ──
    ax_main = axes[n_sidebars]
    ax_main.imshow(matrix, aspect='auto', cmap=cmap, norm=norm,
                   interpolation='nearest')
    ax_main.set_xticks(range(n_modes))
    ax_main.set_xticklabels(modes, fontsize=8, rotation=45, ha='right')
    ax_main.set_yticks([])

    # Grid lines between cells
    for x in np.arange(-0.5, n_modes, 1):
        ax_main.axvline(x, color='white', linewidth=0.3)
    for y in np.arange(-0.5, n_peaks, max(1, n_peaks // 30)):
        ax_main.axhline(y, color='white', linewidth=0.15)

    # ── Legend ──
    present_reasons = set(matrix.flatten()) - {ABSENT}
    legend_patches = []
    for idx, reason in enumerate(_REASON_ORDER):
        if idx in present_reasons:
            legend_patches.append(
                Patch(facecolor=_REASON_COLORS[reason], edgecolor='gray',
                      linewidth=0.5, label=_REASON_LABELS[reason])
            )
    if ABSENT in set(matrix.flatten()):
        legend_patches.append(
            Patch(facecolor='#f0f0f0', edgecolor='gray',
                  linewidth=0.5, label='Not evaluated')
        )
    ax_main.legend(
        handles=legend_patches, loc='upper left',
        bbox_to_anchor=(1.02, 1.0), fontsize=7,
        frameon=True, framealpha=0.95, borderaxespad=0,
    )

    end_label = "5' (TSS)" if end_type == 'tss' else "3' (TTS)"
    # (no suptitle — pub-quality)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout(rect=(0, 0, 0.88, 0.97))
    return _save_figure(fig, output_path)


# Known assembler mode names (used for mode extraction from filenames).
# Longest patterns first so multi-word modes match before their suffixes.
_KNOWN_MODES = [
    'density-asymmetric-softclip', 'density-asymmetric', 'density-plain',
    'density-strict', 'isoquant_pacbio', 'bambu_default', 'more-ends',
    'k-means', 'default',
]


def extract_mode_from_filename(filename: str) -> str:
    """Extract the transcriptome_mode from a per-peak reason TSV filename.

    Expected pattern: ..._<mode>_transcriptome_{cage,quantseq}_peak_reasons.tsv

    We first try to match known mode names directly, then fall back to a regex.
    """
    stem = Path(filename).stem
    # Strategy 1: look for known mode names immediately before "_transcriptome_"
    for mode in _KNOWN_MODES:
        if f'_{mode}_transcriptome_' in stem:
            return mode
    # Strategy 2: regex — capture the segment right before '_transcriptome_'
    # Use a greedy prefix .* so we get the LAST match (mode, not alignment tag)
    m = re.search(r'_([a-zA-Z][a-zA-Z0-9_-]+)_transcriptome_', stem)
    if m:
        return m.group(1)
    return stem


def load_peak_reason_tsvs(tsv_paths):
    """Load per-peak reason TSVs and return per-mode dicts.

    Returns:
        mode_reasons: Dict[mode -> Dict[peak_id -> reason]]
        mode_scores:  Dict[mode -> Dict[peak_id -> score]]
        mode_reads:   Dict[mode -> Dict[peak_id -> read_count]]
    """
    mode_reasons = {}
    mode_scores = {}
    mode_reads = {}

    for tsv_path in tsv_paths:
        mode = extract_mode_from_filename(str(tsv_path))
        reasons = {}
        scores = {}
        reads = {}

        with open(tsv_path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                pid = row['peak_id']
                reasons[pid] = row['reason']
                try:
                    scores[pid] = float(row['score'])
                except (ValueError, TypeError):
                    scores[pid] = 0.0
                try:
                    reads[pid] = float(row.get('read_count', 0))
                except (ValueError, TypeError):
                    reads[pid] = 0.0

        mode_reasons[mode] = reasons
        mode_scores[mode] = scores
        mode_reads[mode] = reads

    return mode_reasons, mode_scores, mode_reads


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--cage-tsvs', nargs='*', default=[],
                        help='Per-peak reason TSVs for CAGE (5\') end')
    parser.add_argument('--quantseq-tsvs', nargs='*', default=[],
                        help='Per-peak reason TSVs for QuantSeq (3\') end')
    parser.add_argument('--output-prefix', required=True,
                        help='Output filename prefix (e.g. testname)')
    parser.add_argument('--max-peaks', type=int, default=300,
                        help='Maximum peaks to display (default: 300)')
    parser.add_argument('--title-prefix', default='',
                        help='Prefix for plot title')
    parser.add_argument('--verbose', action='store_true')

    args = parser.parse_args()

    if not args.cage_tsvs and not args.quantseq_tsvs:
        print("No per-peak reason TSVs provided, nothing to plot", file=sys.stderr)
        sys.exit(0)

    any_success = False

    if args.cage_tsvs:
        cage_reasons, cage_scores, cage_reads = load_peak_reason_tsvs(args.cage_tsvs)
        if args.verbose:
            for mode, reasons in cage_reasons.items():
                print(f"  CAGE mode '{mode}': {len(reasons)} peaks")
        out_path = Path(f"{args.output_prefix}_cage_peak_reason_heatmap.png")
        ok = plot_peak_reason_heatmap(
            mode_peak_reasons=cage_reasons,
            mode_peak_scores=cage_scores,
            output_path=out_path,
            title=f"{args.title_prefix}Peak Recovery by Mode",
            end_type='tss',
            max_peaks=args.max_peaks,
            mode_peak_reads=cage_reads,
        )
        if ok:
            print(f"Saved CAGE heatmap to {out_path}")
            any_success = True
        elif args.verbose:
            print("CAGE heatmap: skipped (no data or matplotlib unavailable)")

    if args.quantseq_tsvs:
        qs_reasons, qs_scores, qs_reads = load_peak_reason_tsvs(args.quantseq_tsvs)
        if args.verbose:
            for mode, reasons in qs_reasons.items():
                print(f"  QuantSeq mode '{mode}': {len(reasons)} peaks")
        out_path = Path(f"{args.output_prefix}_quantseq_peak_reason_heatmap.png")
        ok = plot_peak_reason_heatmap(
            mode_peak_reasons=qs_reasons,
            mode_peak_scores=qs_scores,
            output_path=out_path,
            title=f"{args.title_prefix}Peak Recovery by Mode",
            end_type='tts',
            max_peaks=args.max_peaks,
            mode_peak_reads=qs_reads,
        )
        if ok:
            print(f"Saved QuantSeq heatmap to {out_path}")
            any_success = True
        elif args.verbose:
            print("QuantSeq heatmap: skipped (no data or matplotlib unavailable)")

    sys.exit(0 if any_success else 1)


if __name__ == '__main__':
    main()
