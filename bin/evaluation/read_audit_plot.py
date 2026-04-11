#!/usr/bin/env python3
"""TED Read Audit Summary Plot — cross-mode comparison of read classifications.

Generates a multi-panel figure comparing how different TED (or FLAIR) modes
classify reads.

Panels
------
1. **Stacked bar**: proportion of reads in each classification category per mode.
2. **5' delta distribution**: violin/box of |delta_5p| for assigned reads per mode.
3. **3' delta distribution**: violin/box of |delta_3p| for assigned reads per mode.
4. **Assignment rate bar**: fraction of reads assigned vs unassigned per mode.

Input is one or more ``*.read_audit.tsv`` files.  Mode labels are inferred from
the filename or supplied explicitly via ``mode:path`` notation.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # noqa: E402

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Attempt to import shared pub_style; fall back gracefully.
try:
    from pub_style import (
        MODE_COLORS, PALETTE, apply_rc, style_ax, legend_outside, savefig,
    )
except ImportError:
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from pub_style import (
            MODE_COLORS, PALETTE, apply_rc, style_ax, legend_outside, savefig,
        )
    except ImportError:
        # Minimal fallback
        PALETTE = ['#0072B2', '#E69F00', '#009E73', '#D55E00',
                    '#CC79A7', '#56B4E9', '#F0E442', '#999999']
        MODE_COLORS = {}
        def apply_rc(): pass
        def style_ax(ax, **kw):
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        def legend_outside(fig_or_ax, **kw):
            return fig_or_ax.legend(**kw)
        def savefig(fig, path, dpi=300, **kw):
            fig.savefig(path, dpi=dpi, bbox_inches='tight', facecolor='white')
            plt.close(fig)


# ── Classification categories ───────────────────────────────────────────────

CLASS_ORDER = ['kept', 'reassigned_5p', 'reassigned_3p', 'reassigned_both', 'unassigned']

# Colours for stacked-bar classification buckets (distinct from mode colours)
CLASS_COLORS = {
    'kept':             '#009E73',  # bluish green
    'reassigned_5p':    '#0072B2',  # blue
    'reassigned_3p':    '#E69F00',  # orange
    'reassigned_both':  '#D55E00',  # vermillion
    'unassigned':       '#999999',  # grey
}

CLASS_LABELS = {
    'kept':             'Kept (both ends)',
    'reassigned_5p':    "Reassigned 5'",
    'reassigned_3p':    "Reassigned 3'",
    'reassigned_both':  'Reassigned both',
    'unassigned':       'Unassigned',
}


# ── Data loading ────────────────────────────────────────────────────────────

_MODE_PAT = re.compile(
    r'_([A-Za-z0-9][A-Za-z0-9_-]*?)_transcriptome\.read_audit\.tsv$'
)

def infer_mode_from_path(path: str) -> str:
    """Try to extract the transcriptome mode from a standard filename."""
    m = _MODE_PAT.search(Path(path).name)
    if m:
        return m.group(1)
    # Fallback: stem minus .read_audit
    stem = Path(path).stem
    if stem.endswith('.read_audit'):
        stem = stem[:-len('.read_audit')]
    return stem


def load_audit_tsvs(inputs: list[str]) -> pd.DataFrame:
    """Load one or more ``read_audit.tsv`` files, tagging each with a *mode*.

    Accepts either plain paths (mode inferred from filename) or
    ``mode:path`` pairs.
    """
    frames = []
    for entry in inputs:
        if ':' in entry and not entry.startswith('/'):
            mode, path = entry.split(':', 1)
        else:
            path = entry
            mode = infer_mode_from_path(path)
        try:
            df = pd.read_csv(path, sep='\t')
            df['mode'] = mode
            frames.append(df)
        except Exception as exc:
            logging.warning('Could not read %s: %s', path, exc)

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def _mode_color(mode: str) -> str:
    if mode in MODE_COLORS:
        return MODE_COLORS[mode]
    return PALETTE[hash(mode) % len(PALETTE)]


def _order_modes(modes: list[str]) -> list[str]:
    """Return modes in a visually meaningful order (baseline first, then sorted)."""
    priority = ['baseline', 'default', 'ted-default']
    ordered = [m for m in priority if m in modes]
    for m in sorted(modes):
        if m not in ordered:
            ordered.append(m)
    return ordered


# ── Plotting ────────────────────────────────────────────────────────────────

def create_read_audit_plots(
    df: pd.DataFrame,
    output_dir: Path,
    title_prefix: str = '',
) -> bool:
    """Create individual read-audit summary figures in output_dir.

    Returns True on success.
    """
    if df.empty:
        logging.error('No data to plot.')
        return False

    apply_rc()

    modes = _order_modes(list(df['mode'].unique()))
    n_modes = len(modes)
    fw = min(7.2, 4 + n_modes * 0.4)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Panel 1: Stacked classification bar
    fig, ax = plt.subplots(figsize=(fw, 3.5))
    _plot_stacked_bar(ax, df, modes, title_prefix)
    fig.tight_layout()
    savefig(fig, output_dir / 'read_classification.png')

    # Panel 2: Assignment rate
    fig, ax = plt.subplots(figsize=(fw, 3.0))
    _plot_assignment_rate(ax, df, modes, title_prefix)
    fig.tight_layout()
    savefig(fig, output_dir / 'assignment_rate.png')

    # Panel 3: 5' delta distribution
    fig, ax = plt.subplots(figsize=(fw, 3.0))
    _plot_delta_distribution(ax, df, modes, 'delta_5p', "5'", title_prefix)
    fig.tight_layout()
    savefig(fig, output_dir / 'delta_5prime.png')

    # Panel 4: 3' delta distribution
    fig, ax = plt.subplots(figsize=(fw, 3.0))
    _plot_delta_distribution(ax, df, modes, 'delta_3p', "3'", title_prefix)
    fig.tight_layout()
    savefig(fig, output_dir / 'delta_3prime.png')

    return True


def _plot_stacked_bar(ax, df: pd.DataFrame, modes: list[str], title_prefix: str):
    """Stacked horizontal bar — classification proportions."""
    # Compute proportions
    rows = []
    for mode in modes:
        sub = df[df['mode'] == mode]
        total = len(sub)
        if total == 0:
            continue
        row = {'mode': mode}
        for cls in CLASS_ORDER:
            row[cls] = (sub['classification'] == cls).sum() / total
        rows.append(row)

    if not rows:
        return

    prop_df = pd.DataFrame(rows)
    y_pos = np.arange(len(prop_df))

    left = np.zeros(len(prop_df))
    for cls in CLASS_ORDER:
        widths = prop_df[cls].values
        ax.barh(y_pos, widths, left=left, height=0.7,
                color=CLASS_COLORS[cls], label=CLASS_LABELS[cls],
                edgecolor='none')
        left += widths

    ax.set_yticks(y_pos)
    ax.set_yticklabels(prop_df['mode'])
    ax.set_xlim(0, 1)
    ax.set_xlabel('Proportion of reads')
    style_ax(ax)

    # Legend outside (below the plot)
    handles, labels = ax.get_legend_handles_labels()
    legend_outside(ax, handles=handles, labels=labels,
                   loc='upper left', bbox_to_anchor=(1.02, 1.0), ncol=1)


def _plot_assignment_rate(ax, df: pd.DataFrame, modes: list[str], title_prefix: str):
    """Bar chart — fraction of reads assigned (not unassigned)."""
    rates = []
    labels = []
    for mode in modes:
        sub = df[df['mode'] == mode]
        total = len(sub)
        if total == 0:
            continue
        assigned = (sub['classification'] != 'unassigned').sum()
        rates.append(assigned / total)
        labels.append(mode)

    if not rates:
        return

    x_pos = np.arange(len(rates))
    colors = [_mode_color(m) for m in labels]
    ax.bar(x_pos, rates, color=colors, edgecolor='none', width=0.7)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylim(0, 1.05)
    ax.set_ylabel('Assignment rate')
    style_ax(ax)


def _plot_delta_distribution(
    ax, df: pd.DataFrame, modes: list[str],
    col: str, label: str, title_prefix: str,
):
    """Box plot of |delta| (bp) for assigned reads."""
    assigned = df[df['classification'] != 'unassigned'].copy()
    if col not in assigned.columns:
        return

    assigned[col] = pd.to_numeric(assigned[col], errors='coerce')
    assigned = assigned.dropna(subset=[col])
    assigned['abs_delta'] = assigned[col].abs()

    # Collect per-mode data
    data_by_mode = []
    mode_labels = []
    for mode in modes:
        sub = assigned[assigned['mode'] == mode]['abs_delta']
        if len(sub) == 0:
            continue
        data_by_mode.append(sub.values)
        mode_labels.append(mode)

    if not data_by_mode:
        return

    bp = ax.boxplot(
        data_by_mode,
        positions=range(len(data_by_mode)),
        widths=0.6,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color='black', linewidth=1.2),
    )
    for i, (patch, mode) in enumerate(zip(bp['boxes'], mode_labels)):
        patch.set_facecolor(_mode_color(mode))
        patch.set_edgecolor('none')
        patch.set_alpha(0.8)

    ax.set_xticks(range(len(mode_labels)))
    ax.set_xticklabels(mode_labels, rotation=45, ha='right')
    ax.set_ylabel(f'|{label} delta| (bp)')
    style_ax(ax)


# ── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--input', '-i', nargs='+', required=True,
        help='read_audit.tsv files.  Use mode:path or plain path (mode inferred).')
    parser.add_argument(
        '--output', '-o', required=True,
        help='Output directory for individual plots.')
    parser.add_argument(
        '--title-prefix', default='',
        help='Optional title prefix for all panels.')
    parser.add_argument('--verbose', '-v', action='store_true')
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s  %(levelname)-8s  %(message)s',
    )

    df = load_audit_tsvs(args.input)
    if df.empty:
        logging.error('No data loaded — exiting.')
        sys.exit(1)

    modes = df['mode'].unique()
    logging.info('Loaded %d reads across %d modes: %s',
                 len(df), len(modes), list(modes))

    out = Path(args.output)

    if create_read_audit_plots(df, out, args.title_prefix):
        logging.info('Saved read-audit plots to %s', out)
    else:
        sys.exit(1)


if __name__ == '__main__':
    main()
