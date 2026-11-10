#!/usr/bin/env python3
"""
Cross-run summary plots for peak recovery reasons (signal vs support analysis).

Reads merged evaluation TSVs containing per-peak reason counts and produces
a 2-row figure:
  Row 1: 5' (TSS / CAGE)
  Row 2: 3' (TTS / dRNA)

Each row has two panels:
  Left:  Stacked horizontal bar chart of reason proportions per mode
  Right: Recovery rate breakdown (recovered % with reason composition)

Usage:
    python signal_support_summary_plot.py \\
        --input eval1.tsv eval2.tsv ... \\
        --output signal_support_summary.png \\
        [--title-prefix "test_name: "] [--verbose]
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pub_style import style_ax, savefig, REASON_COLORS as _PUB_REASON_COLORS, MODE_COLORS as _PUB_MODE_COLORS, PALETTE

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Reason colors — colorblind-safe via pub_style
REASON_COLORS = _PUB_REASON_COLORS

REASON_LABELS = {
    'recovered':           'Recovered',
    'no_reads':            'No long reads',
    'single_exon_only':    'Single-exon only',
    'near_miss':           'Near miss',
    'reads_unassigned':    'Reads unassigned',
    'reads_redirected':    'Reads redirected',
    'proximal_apa':        'Proximal APA',
    'trailing_truncation': 'Trailing / truncation',
    'other_missed':        'Other missed',
    'no_isoform_model':    'No isoform model',
    'alignment_filtered':  'Alignment filtered',
    'end_absorbed':        'End absorbed',
    'end_spread':          'End spread',
    'low_signal':          'Low signal',
}

REASON_ORDER = [
    'recovered', 'near_miss', 'reads_redirected', 'proximal_apa',
    'trailing_truncation', 'single_exon_only', 'reads_unassigned',
    'no_reads', 'no_isoform_model', 'alignment_filtered',
    'end_absorbed', 'end_spread', 'low_signal', 'other_missed',
]

MODE_COLORS = _PUB_MODE_COLORS

DEFAULT_MODE_ORDER = [
    "default",
    "ted-2d",
    "ted-1d2d",
    "density-asymmetric",
    "density-asymmetric-softclip",
    "density-plain",
    "density-strict",
    "k-means",
    "more-ends",
    "bambu_default",
    "isoquant_pacbio",
]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_evaluation_tsvs(paths):
    """Load evaluation TSVs and return list of row dicts."""
    rows = []
    for p in paths:
        with open(p) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                rows.append(row)
    return rows


def extract_reason_counts(row, prefix):
    """Extract reason counts for a given end type prefix (5prime or 3prime)."""
    counts = {}
    for reason in REASON_ORDER:
        key = f"{prefix}_reason_{reason}"
        val = row.get(key, '')
        try:
            counts[reason] = int(float(val)) if val else 0
        except (ValueError, TypeError):
            counts[reason] = 0
    return counts


# ---------------------------------------------------------------------------
# Plotting — individual figures
# ---------------------------------------------------------------------------

def plot_reason_bar_single(modes, mode_reason_counts, output_path, title):
    """Standalone stacked horizontal bar chart of reason proportions per mode."""
    y_positions = np.arange(len(modes))
    bar_height = 0.65

    all_reasons = set()
    for counts in mode_reason_counts.values():
        all_reasons.update(k for k, v in counts.items() if v > 0)
    plot_reasons = [r for r in REASON_ORDER if r in all_reasons]

    fig, ax = plt.subplots(figsize=(max(4.5, 3.5), max(2.5, 0.5 * len(modes) + 1.0)))

    left_offsets = np.zeros(len(modes))
    for reason in plot_reasons:
        widths = []
        for mode in modes:
            total = sum(mode_reason_counts[mode].values())
            count = mode_reason_counts[mode].get(reason, 0)
            pct = (count / total * 100) if total > 0 else 0
            widths.append(pct)
        widths = np.array(widths)
        ax.barh(y_positions, widths, left=left_offsets, height=bar_height,
                color=REASON_COLORS.get(reason, '#95a5a6'),
                label=REASON_LABELS.get(reason, reason),
                edgecolor='white', linewidth=0.5)
        left_offsets += widths

    ax.set_yticks(y_positions)
    ax.set_yticklabels(modes, fontsize=7)
    ax.set_xlim(0, 100)
    ax.invert_yaxis()
    style_ax(ax, xlabel='Percentage of Peaks (%)')

    # Reason legend below
    from matplotlib.patches import Patch
    legend_patches = [Patch(facecolor=REASON_COLORS.get(r, '#95a5a6'),
                            edgecolor='white', label=REASON_LABELS.get(r, r))
                      for r in plot_reasons]
    ax.legend(handles=legend_patches, loc='upper center', bbox_to_anchor=(0.5, -0.15),
              ncol=min(4, len(legend_patches)), fontsize=6, frameon=False)

    fig.tight_layout()
    savefig(fig, Path(output_path), dpi=300)


def plot_recovery_rate_single(modes, mode_reason_counts, output_path, title):
    """Standalone bar chart comparing overall recovery rate per mode."""
    y_positions = np.arange(len(modes))
    bar_height = 0.65

    recovery_rates = []
    for mode in modes:
        total = sum(mode_reason_counts[mode].values())
        recovered = mode_reason_counts[mode].get('recovered', 0)
        rate = (recovered / total * 100) if total > 0 else 0
        recovery_rates.append(rate)

    fig, ax = plt.subplots(figsize=(max(4.0, 3.5), max(2.5, 0.5 * len(modes) + 1.0)))

    colors = [MODE_COLORS.get(mode, PALETTE[hash(mode) % len(PALETTE)]) for mode in modes]
    bars = ax.barh(y_positions, recovery_rates, height=bar_height,
                   color=colors, edgecolor='none', linewidth=0, alpha=0.85)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(modes, fontsize=7)
    ax.set_xlim(0, 105)
    ax.invert_yaxis()
    style_ax(ax, xlabel='Peak Recall (%)')

    for bar, rate, mode in zip(bars, recovery_rates, modes):
        total = sum(mode_reason_counts[mode].values())
        recovered = mode_reason_counts[mode].get('recovered', 0)
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2,
                f'{rate:.1f}% ({recovered}/{total})',
                va='center', fontsize=8, color='#333333')

    fig.tight_layout()
    savefig(fig, Path(output_path), dpi=300)


def create_summary_plots(rows, output_dir, title_prefix=""):
    """Create individual signal support summary figures in output_dir."""
    mode_5prime = {}
    mode_3prime = {}

    for row in rows:
        mode = row.get('pipeline_mode', row.get('transcriptome_mode', 'unknown'))
        counts_5 = extract_reason_counts(row, '5prime')
        counts_3 = extract_reason_counts(row, '3prime')

        if sum(counts_5.values()) > 0:
            if mode not in mode_5prime:
                mode_5prime[mode] = counts_5
            else:
                for reason, count in counts_5.items():
                    mode_5prime[mode][reason] = mode_5prime[mode].get(reason, 0) + count
        if sum(counts_3.values()) > 0:
            if mode not in mode_3prime:
                mode_3prime[mode] = counts_3
            else:
                for reason, count in counts_3.items():
                    mode_3prime[mode][reason] = mode_3prime[mode].get(reason, 0) + count

    if not mode_5prime and not mode_3prime:
        print("No peak reason data found in evaluation files", file=sys.stderr)
        return False

    all_modes = sorted(set(list(mode_5prime.keys()) + list(mode_3prime.keys())),
                       key=lambda m: DEFAULT_MODE_ORDER.index(m) if m in DEFAULT_MODE_ORDER else 999)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if mode_5prime:
        modes_5 = [m for m in all_modes if m in mode_5prime]
        plot_reason_bar_single(modes_5, mode_5prime,
                               output_dir / "reason_5prime.png",
                               f"{title_prefix}TSS (5\u2032 Peaks): Peak Recovery Reasons")
        plot_recovery_rate_single(modes_5, mode_5prime,
                                  output_dir / "recovery_rate_5prime.png",
                                  f"{title_prefix}TSS (5\u2032 Peaks): Peak Recall")

    if mode_3prime:
        modes_3 = [m for m in all_modes if m in mode_3prime]
        plot_reason_bar_single(modes_3, mode_3prime,
                               output_dir / "reason_3prime.png",
                               f"{title_prefix}TTS (3\u2032 Peaks): Peak Recovery Reasons")
        plot_recovery_rate_single(modes_3, mode_3prime,
                                  output_dir / "recovery_rate_3prime.png",
                                  f"{title_prefix}TTS (3\u2032 Peaks): Peak Recall")

    print(f"Saved signal support summary plots to {output_dir}")
    return True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', nargs='+', required=True,
                        help='Evaluation TSV files')
    parser.add_argument('--output', required=True,
                        help='Output directory for individual plots')
    parser.add_argument('--title-prefix', default='',
                        help='Prefix for plot title')
    parser.add_argument('--verbose', action='store_true')

    args = parser.parse_args()

    rows = load_evaluation_tsvs(args.input)
    if args.verbose:
        print(f"Loaded {len(rows)} evaluation rows from {len(args.input)} files")

    success = create_summary_plots(rows, args.output, args.title_prefix)
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
