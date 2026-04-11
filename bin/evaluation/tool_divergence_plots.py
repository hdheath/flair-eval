#!/usr/bin/env python3
"""
Tool-divergence summary plots.

Creates a multi-panel figure:
  Panel 1: Heatmap of pairwise Jaccard index for SJ, TSS, TTS
  Panel 2: Grouped bar chart comparing SJ vs end Jaccard per tool pair
  Panel 3: Motif collapse summary — bar chart by motif category
  Panel 4: Distinct 3' ends histogram + collapse rate annotation

Usage:
    python tool_divergence_plots.py \\
        --eval-files *.evaluation.tsv \\
        --isoform-beds tool1:path1.bed tool2:path2.bed ... \\
        --reads-bed reads.bed \\
        --map-files tool1:map1.txt tool2:map2.txt ... \\
        --genome genome.fa \\
        --output divergence_summary.png \\
        [--window-5prime 50] [--window-3prime 5] \\
        [--title-prefix "Test: "] [--verbose]

Can also be imported and called programmatically.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import numpy as np
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# Publication-style helpers (style_ax, savefig, PALETTE, etc.)
try:
    from pub_style import style_ax, savefig, apply_rc, PALETTE, MODE_COLORS
    apply_rc()
except ImportError:
    # Fallback stubs so the script can still run without pub_style
    def style_ax(ax, **kwargs):
        if 'title' in kwargs:
            ax.set_title(kwargs['title'])
        if 'xlabel' in kwargs:
            ax.set_xlabel(kwargs['xlabel'])
        if 'ylabel' in kwargs:
            ax.set_ylabel(kwargs['ylabel'])
    def savefig(fig, path, **kwargs):
        fig.savefig(path, dpi=kwargs.get('dpi', 300), bbox_inches='tight')
        plt.close(fig)
    PALETTE = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
               '#D55E00', '#CC79A7', '#000000']
    MODE_COLORS = {}


def plot_jaccard_heatmaps(
    jaccard_result: dict,
    output_path: Path,
    title_prefix: str = "",
    dpi: int = 300,
) -> bool:
    """Create a 3-panel heatmap showing pairwise Jaccard for SJ, TSS, TTS.

    Args:
        jaccard_result: output from compute_pairwise_jaccard()
        output_path: where to save the PNG
        title_prefix: prepended to the figure title
        dpi: output resolution

    Returns True on success.
    """
    if not HAS_MATPLOTLIB:
        print("Warning: matplotlib not available, skipping plot", file=sys.stderr)
        return False

    tool_names = jaccard_result.get('tool_names', [])
    if len(tool_names) < 2:
        print("Warning: Need ≥2 tools for heatmap", file=sys.stderr)
        return False

    n = len(tool_names)
    matrices = {
        'Splice Junctions': jaccard_result.get('sj_jaccard_matrix', {}),
        'TSS (5\' ends)': jaccard_result.get('tss_jaccard_matrix', {}),
        'TTS (3\' ends)': jaccard_result.get('tts_jaccard_matrix', {}),
    }

    fig, axes = plt.subplots(1, 3, figsize=(5 * 3 + 1.5, 5))
    # (no suptitle — pub-quality)

    cmap = plt.cm.RdYlGn

    for ax_idx, (label, matrix) in enumerate(matrices.items()):
        ax = axes[ax_idx]

        # Build symmetric matrix
        data = np.ones((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                a, b = tool_names[i], tool_names[j]
                val = matrix.get((a, b), matrix.get((b, a), 0.0))
                data[i, j] = val
                data[j, i] = val

        im = ax.imshow(data, cmap=cmap, vmin=0, vmax=1, aspect='equal')

        # Annotate cells
        for i in range(n):
            for j in range(n):
                val = data[i, j]
                color = 'white' if val < 0.4 or val > 0.85 else 'black'
                ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                        fontsize=7, color=color, fontweight='normal')

        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(tool_names, rotation=45, ha='right', fontsize=7)
        ax.set_yticklabels(tool_names, fontsize=7)
        style_ax(ax, title=label)

    # Shared colorbar
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.65])
    fig.colorbar(im, cax=cbar_ax, label='Jaccard Index')

    fig.tight_layout(rect=[0, 0, 0.9, 1.0])
    savefig(fig, output_path, dpi=dpi)
    print(f"Saved Jaccard heatmap to {output_path}")
    return True


def plot_jaccard_comparison_bars(
    jaccard_result: dict,
    output_path: Path,
    title_prefix: str = "",
    dpi: int = 300,
) -> bool:
    """Grouped bar chart: SJ vs TSS vs TTS Jaccard per tool pair.

    Visually demonstrates the hypothesis that SJ agreement >> end agreement.
    """
    if not HAS_MATPLOTLIB:
        return False

    tool_names = jaccard_result.get('tool_names', [])
    sj_matrix = jaccard_result.get('sj_jaccard_matrix', {})
    tss_matrix = jaccard_result.get('tss_jaccard_matrix', {})
    tts_matrix = jaccard_result.get('tts_jaccard_matrix', {})

    from itertools import combinations
    pairs = list(combinations(tool_names, 2))
    if not pairs:
        return False

    pair_labels = [f"{a} vs {b}" for a, b in pairs]
    sj_vals = [sj_matrix.get(p, 0.0) for p in pairs]
    tss_vals = [tss_matrix.get(p, 0.0) for p in pairs]
    tts_vals = [tts_matrix.get(p, 0.0) for p in pairs]

    x = np.arange(len(pairs))
    width = 0.25

    fig, ax = plt.subplots(figsize=(max(6, len(pairs) * 2.5), 5))
    bars_sj = ax.bar(x - width, sj_vals, width, label='Splice Junctions',
                     color=PALETTE[2], edgecolor='none', alpha=0.85)
    bars_tss = ax.bar(x, tss_vals, width, label='TSS (5\' ends)',
                      color=PALETTE[4], edgecolor='none', alpha=0.85)
    bars_tts = ax.bar(x + width, tts_vals, width, label='TTS (3\' ends)',
                      color=PALETTE[0], edgecolor='none', alpha=0.85)

    # Value labels
    for bars in [bars_sj, bars_tss, bars_tts]:
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.01,
                    f'{h:.2f}', ha='center', va='bottom', fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(pair_labels, rotation=30, ha='right', fontsize=7)
    ax.set_ylim(0, 1.15)
    ax.axhline(0.9, color='gray', linestyle='--', alpha=0.3, linewidth=0.8)
    ax.text(len(pairs) - 0.5, 0.91, 'SJ threshold (90%)',
            fontsize=7, color='gray', ha='right')
    style_ax(ax, ylabel='Jaccard Index', xlabel='Tool Pair',
             title='SJ vs End Agreement Across Tools')
    ax.legend(fontsize=7, frameon=False)

    fig.tight_layout()
    savefig(fig, output_path, dpi=dpi)
    print(f"Saved Jaccard comparison bars to {output_path}")
    return True


def plot_motif_collapse_summary(
    collapse_result: dict,
    output_path: Path,
    title_prefix: str = "",
    dpi: int = 300,
) -> bool:
    """Two-panel motif collapse figure:
    Left: bar chart of collapse events by motif category
    Right: histogram of distinct ends per group + annotations
    """
    if not HAS_MATPLOTLIB:
        return False

    end_type = collapse_result.get('end_type', 'tts')
    end_label = "5\u2032" if end_type == 'tss' else "3\u2032"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # ── Left panel: collapse by category ──
    category_counts = collapse_result.get('_per_category_collapse_counts', {})
    cat_colors_3p = {'PAS': '#e74c3c', 'ARE': '#f39c12', 'miRNA': '#9b59b6'}
    cat_colors_5p = {'Inr': '#f39c12', 'TOP': '#9b59b6', 'PRTE': '#27ae60'}
    cat_colors = cat_colors_5p if end_type == 'tss' else cat_colors_3p

    if category_counts:
        categories = sorted(category_counts.keys())
        counts = [category_counts[c] for c in categories]
        bar_colors = [cat_colors.get(c, '#95a5a6') for c in categories]
        bars = ax1.bar(range(len(categories)), counts, color=bar_colors,
                       edgecolor='none', alpha=0.85)
        for bar, count in zip(bars, counts):
            ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                     str(count), ha='center', va='bottom', fontsize=7)
        ax1.set_xticks(range(len(categories)))
        ax1.set_xticklabels(categories, fontsize=7)
    style_ax(ax1, xlabel='Motif Category', ylabel='Collapse Events',
             title=f'{end_label} Collapse Events by Motif Category')

    # Summary text box
    total_groups = collapse_result.get('motif_collapse_total_groups', 0)
    affected = collapse_result.get('motif_collapse_affected_groups', 0)
    rate = collapse_result.get('motif_collapse_rate', 0)
    total_events = collapse_result.get('total_collapse_events', 0)
    motif_total = collapse_result.get('reads_with_motif_total', 0)
    motif_collapsed = collapse_result.get('reads_with_motif_collapsed', 0)
    motif_frac = collapse_result.get('motif_collapse_fraction', 0)

    summary_text = (
        f"Groups analyzed: {total_groups}\n"
        f"Groups affected: {affected} ({rate:.1%})\n"
        f"Total events: {total_events}\n"
        f"Motif reads: {motif_collapsed}/{motif_total} ({motif_frac:.1%})"
    )
    ax1.text(0.98, 0.98, summary_text, transform=ax1.transAxes,
             fontsize=8, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='wheat', alpha=0.7))

    # ── Right panel: histogram of distinct ends per group ──
    distinct_counts = collapse_result.get('_distinct_end_counts',
                                          collapse_result.get('_distinct_3p_counts', []))
    if distinct_counts:
        max_val = min(max(distinct_counts), 20)
        bins = range(2, max_val + 2)
        ax2.hist(distinct_counts, bins=bins, color=PALETTE[1], edgecolor='none',
                 alpha=0.85, rwidth=0.85)
        ax2.set_xticks(range(2, max_val + 1))
    from matplotlib.ticker import MaxNLocator
    ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
    style_ax(ax2, xlabel=f'Distinct {end_label} Ends per Junction Chain Group',
             ylabel='Number of Groups',
             title=f'Distribution of {end_label} End Diversity')

    mean_distinct = collapse_result.get('distinct_ends_per_group_mean',
                                        collapse_result.get('distinct_3p_ends_per_group_mean', 0))
    ax2.axvline(mean_distinct, color='red', linestyle='--', linewidth=1.5)
    ax2.text(mean_distinct + 0.2, ax2.get_ylim()[1] * 0.9,
             f'Mean: {mean_distinct:.1f}', fontsize=7, color='red')

    fig.tight_layout()
    savefig(fig, output_path, dpi=dpi)
    print(f"Saved {end_label} motif collapse summary to {output_path}")
    return True


def plot_motif_collapse_paired(
    collapse_5p: dict,
    collapse_3p: dict,
    output_path: Path,
    title_prefix: str = "",
    dpi: int = 300,
) -> bool:
    """Side-by-side 5\u2032 and 3\u2032 motif collapse for one assembler.

    Layout: 1 × 2
      Left:  5\u2032 RNA motif collapse by category (Inr / TOP)
      Right: 3\u2032 UTR motif collapse by category (PAS / ARE / miRNA)
    """
    if not HAS_MATPLOTLIB:
        return False

    fig, (ax5, ax3) = plt.subplots(1, 2, figsize=(14, 5.5))

    cat_colors_5p = {'Inr': '#f39c12', 'TOP': '#9b59b6', 'PRTE': '#27ae60'}
    cat_colors_3p = {'PAS': '#e74c3c', 'ARE': '#f39c12', 'miRNA': '#9b59b6'}

    def _draw_panel(ax, collapse_result, cat_colors, label):
        cat_counts = collapse_result.get('_per_category_collapse_counts', {})
        if cat_counts:
            categories = sorted(cat_counts.keys())
            counts = [cat_counts[c] for c in categories]
            bar_cols = [cat_colors.get(c, '#95a5a6') for c in categories]
            bars = ax.bar(range(len(categories)), counts, color=bar_cols,
                          edgecolor='none', alpha=0.85)
            for bar, cnt in zip(bars, counts):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                        str(cnt), ha='center', va='bottom', fontsize=7)
            ax.set_xticks(range(len(categories)))
            ax.set_xticklabels(categories, fontsize=7)
        style_ax(ax, xlabel='Motif Category', ylabel='Collapse Events',
                 title=f'{label} Motif Collapse')

        # Summary annotation
        total_groups = collapse_result.get('motif_collapse_total_groups', 0)
        affected = collapse_result.get('motif_collapse_affected_groups', 0)
        rate = collapse_result.get('motif_collapse_rate', 0)
        total_events = collapse_result.get('total_collapse_events', 0)
        motif_total = collapse_result.get('reads_with_motif_total', 0)
        motif_collapsed = collapse_result.get('reads_with_motif_collapsed', 0)
        motif_frac = collapse_result.get('motif_collapse_fraction', 0)
        txt = (
            f"Groups: {affected}/{total_groups} ({rate:.1%})\n"
            f"Events: {total_events}\n"
            f"Motif reads: {motif_collapsed}/{motif_total} ({motif_frac:.1%})"
        )
        ax.text(0.98, 0.98, txt, transform=ax.transAxes, fontsize=8,
                va='top', ha='right',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='wheat', alpha=0.7))

    _draw_panel(ax5, collapse_5p, cat_colors_5p, "5\u2032 Promoter")
    _draw_panel(ax3, collapse_3p, cat_colors_3p, "3\u2032 UTR")

    fig.tight_layout()
    savefig(fig, output_path, dpi=dpi)
    print(f"Saved paired motif collapse to {output_path}")
    return True


def plot_tool_divergence_dashboard(
    jaccard_result: dict,
    collapse_results_3p: Optional[Dict[str, dict]] = None,
    collapse_results_5p: Optional[Dict[str, dict]] = None,
    output_path: Path = None,
    title_prefix: str = "",
    dpi: int = 300,
) -> bool:
    """Cross-tool divergence dashboard.

    Top row: Jaccard comparison bars + SJ−End delta heatmap.
    Bottom row (if collapse data): per-tool motif collapse comparison
    with side-by-side 5\u2032 and 3\u2032 bars.
    """
    if not HAS_MATPLOTLIB:
        return False

    # Accept legacy single-dict form for backward compat
    if collapse_results_3p is not None and not isinstance(collapse_results_3p, dict):
        collapse_results_3p = None
    if collapse_results_5p is not None and not isinstance(collapse_results_5p, dict):
        collapse_results_5p = None

    def _has_data(d):
        return (d is not None and len(d) > 0
                and any(v.get('motif_collapse_total_groups', 0) > 0
                        for v in d.values()))

    has_3p = _has_data(collapse_results_3p)
    has_5p = _has_data(collapse_results_5p)
    has_collapse = has_3p or has_5p

    if has_collapse:
        fig = plt.figure(figsize=(16, 11))
        gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)
        ax_bars = fig.add_subplot(gs[0, 0])
        ax_heat = fig.add_subplot(gs[0, 1])
        ax_rate = fig.add_subplot(gs[1, 0])
        ax_frac = fig.add_subplot(gs[1, 1])
    else:
        fig, (ax_bars, ax_heat) = plt.subplots(1, 2, figsize=(14, 5.5))

    fig.suptitle('', fontsize=1)  # placeholder for tight_layout spacing

    # ── Panel 1: Jaccard comparison bars ──
    tool_names = jaccard_result.get('tool_names', [])
    sj_matrix = jaccard_result.get('sj_jaccard_matrix', {})
    tss_matrix = jaccard_result.get('tss_jaccard_matrix', {})
    tts_matrix = jaccard_result.get('tts_jaccard_matrix', {})

    from itertools import combinations
    pairs = list(combinations(tool_names, 2))

    if pairs:
        pair_labels = [f"{a}\nvs\n{b}" for a, b in pairs]
        sj_vals = [sj_matrix.get(p, 0.0) for p in pairs]
        tss_vals = [tss_matrix.get(p, 0.0) for p in pairs]
        tts_vals = [tts_matrix.get(p, 0.0) for p in pairs]

        x = np.arange(len(pairs))
        width = 0.25
        ax_bars.bar(x - width, sj_vals, width, label='Splice Junctions',
                    color=PALETTE[2], edgecolor='none', alpha=0.85)
        ax_bars.bar(x, tss_vals, width, label='TSS (5\')',
                    color=PALETTE[4], edgecolor='none', alpha=0.85)
        ax_bars.bar(x + width, tts_vals, width, label='TTS (3\')',
                    color=PALETTE[0], edgecolor='none', alpha=0.85)
        ax_bars.set_xticks(x)
        ax_bars.set_xticklabels(pair_labels, fontsize=8)
        ax_bars.set_ylim(0, 1.15)
        ax_bars.axhline(0.9, color='gray', linestyle='--', alpha=0.3, linewidth=0.8)
        ax_bars.legend(fontsize=8, loc='lower left', frameon=False)

    style_ax(ax_bars, ylabel='Jaccard Index', title='SJ vs End Agreement')

    # ── Panel 2: Compact heatmap (average across feature types) ──
    n = len(tool_names)
    if n >= 2:
        # Show difference: SJ Jaccard - End Jaccard (positive = ends less consistent)
        data = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                a, b = tool_names[i], tool_names[j]
                sj_j = sj_matrix.get((a, b), 0)
                end_j = (tss_matrix.get((a, b), 0) + tts_matrix.get((a, b), 0)) / 2
                delta = sj_j - end_j
                data[i, j] = delta
                data[j, i] = delta

        im = ax_heat.imshow(data, cmap='RdBu_r', vmin=-0.5, vmax=0.5, aspect='equal')
        for i in range(n):
            for j in range(n):
                if i != j:
                    val = data[i, j]
                    color = 'white' if abs(val) > 0.25 else 'black'
                    ax_heat.text(j, i, f'{val:+.2f}', ha='center', va='center',
                                fontsize=7, color=color, fontweight='normal')
                else:
                    ax_heat.text(j, i, '—', ha='center', va='center',
                                fontsize=7, color='gray')
        ax_heat.set_xticks(range(n))
        ax_heat.set_yticks(range(n))
        ax_heat.set_xticklabels(tool_names, rotation=45, ha='right', fontsize=7)
        ax_heat.set_yticklabels(tool_names, fontsize=7)
        plt.colorbar(im, ax=ax_heat, label='SJ Jaccard − End Jaccard\n(+ve = ends less consistent)',
                     shrink=0.8)

    style_ax(ax_heat, title='SJ \u2212 End Jaccard Delta')

    # ── Bottom row: per-tool motif collapse comparison (if data available) ──
    if has_collapse:
        # Gather per-tool metrics
        all_tools = sorted(set(
            list(collapse_results_3p.keys() if collapse_results_3p else []) +
            list(collapse_results_5p.keys() if collapse_results_5p else [])
        ))

        # Panel 3: collapse rate per tool (5' vs 3')
        x = np.arange(len(all_tools))
        width = 0.35
        rates_3p = [collapse_results_3p.get(t, {}).get('motif_collapse_rate', 0)
                     if collapse_results_3p else 0 for t in all_tools]
        rates_5p = [collapse_results_5p.get(t, {}).get('motif_collapse_rate', 0)
                     if collapse_results_5p else 0 for t in all_tools]

        if has_5p:
            ax_rate.bar(x - width / 2, rates_5p, width, label="5\u2032 (promoter)",
                        color=PALETTE[4], edgecolor='none', alpha=0.85)
        if has_3p:
            offset = width / 2 if has_5p else 0
            ax_rate.bar(x + offset, rates_3p, width, label="3\u2032 (UTR)",
                        color=PALETTE[0], edgecolor='none', alpha=0.85)
        ax_rate.set_xticks(x)
        ax_rate.set_xticklabels(all_tools, fontsize=7, rotation=30, ha='right')
        style_ax(ax_rate, ylabel='Collapse Rate', title='Motif Collapse Rate by Tool')
        ax_rate.legend(fontsize=8, frameon=False)
        ax_rate.set_ylim(0, max(max(rates_3p + rates_5p, default=0) * 1.2, 0.05))

        # Panel 4: motif-read collapse fraction per tool (5' vs 3')
        fracs_3p = [collapse_results_3p.get(t, {}).get('motif_collapse_fraction', 0)
                     if collapse_results_3p else 0 for t in all_tools]
        fracs_5p = [collapse_results_5p.get(t, {}).get('motif_collapse_fraction', 0)
                     if collapse_results_5p else 0 for t in all_tools]

        if has_5p:
            ax_frac.bar(x - width / 2, fracs_5p, width, label="5\u2032 (promoter)",
                        color=PALETTE[4], edgecolor='none', alpha=0.85)
        if has_3p:
            offset = width / 2 if has_5p else 0
            ax_frac.bar(x + offset, fracs_3p, width, label="3\u2032 (UTR)",
                        color=PALETTE[0], edgecolor='none', alpha=0.85)
        ax_frac.set_xticks(x)
        ax_frac.set_xticklabels(all_tools, fontsize=7, rotation=30, ha='right')
        style_ax(ax_frac, ylabel='Motif-Read Collapse Fraction',
                 title='Reads with Motif near Collapse Site')
        ax_frac.legend(fontsize=8, frameon=False)
        ax_frac.set_ylim(0, max(max(fracs_3p + fracs_5p, default=0) * 1.2, 0.05))

    savefig(fig, output_path, dpi=dpi)
    print(f"Saved divergence dashboard to {output_path}")
    return True


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Compute and plot cross-tool divergence metrics')

    parser.add_argument('--isoform-beds', nargs='+', required=True,
                        help='tool:path pairs, e.g. flair:isoforms.bed bambu:isoforms.bed')
    parser.add_argument('--reads-bed', type=Path, default=None,
                        help='Reads BED12 (for motif collapse)')
    parser.add_argument('--map-files', nargs='*', default=None,
                        help='tool:path pairs for read maps (for motif collapse)')
    parser.add_argument('--genome', type=Path, default=None,
                        help='Genome FASTA (for motif collapse)')
    parser.add_argument('--output-prefix', type=str, required=True,
                        help='Prefix for output files')
    parser.add_argument('--output-dir', type=Path, default=Path('.'),
                        help='Output directory')
    parser.add_argument('--window-5prime', type=int, default=50)
    parser.add_argument('--window-3prime', type=int, default=5)
    parser.add_argument('--min-end-distance', type=int, default=10,
                        help='Skip reads whose 3\' end is within this many bp '
                             'of the isoform boundary (default 10)')
    parser.add_argument('--title-prefix', type=str, default='')
    parser.add_argument('--dpi', type=int, default=300)
    parser.add_argument('--verbose', action='store_true')

    args = parser.parse_args()

    # Parse tool:path pairs
    tool_beds: Dict[str, Path] = {}
    for spec in args.isoform_beds:
        if ':' not in spec:
            print(f"Error: --isoform-beds entries must be tool:path, got '{spec}'",
                  file=sys.stderr)
            sys.exit(1)
        tool, path = spec.split(':', 1)
        tool_beds[tool] = Path(path)

    tool_maps: Dict[str, Path] = {}
    if args.map_files:
        for spec in args.map_files:
            if ':' not in spec:
                continue
            tool, path = spec.split(':', 1)
            tool_maps[tool] = Path(path)

    # Import the analysis module directly, bypassing __init__.py which pulls
    # in heavy deps (scipy via dexseq_ends) that aren't needed here.
    import importlib.util
    _eval_dir = Path(__file__).resolve().parent
    # First load utils and read_analysis so their relative imports resolve
    for _mod_name in ('utils', 'read_analysis'):
        _spec = importlib.util.spec_from_file_location(
            f'evaluation.{_mod_name}', _eval_dir / f'{_mod_name}.py',
            submodule_search_locations=[])
        _mod = importlib.util.module_from_spec(_spec)
        sys.modules[f'evaluation.{_mod_name}'] = _mod
        _spec.loader.exec_module(_mod)
    # Ensure 'evaluation' package exists in sys.modules
    if 'evaluation' not in sys.modules:
        import types
        _pkg = types.ModuleType('evaluation')
        _pkg.__path__ = [str(_eval_dir)]
        sys.modules['evaluation'] = _pkg
    _spec = importlib.util.spec_from_file_location(
        'evaluation.tool_divergence', _eval_dir / 'tool_divergence.py',
        submodule_search_locations=[])
    _td_mod = importlib.util.module_from_spec(_spec)
    sys.modules['evaluation.tool_divergence'] = _td_mod
    _spec.loader.exec_module(_td_mod)
    compute_pairwise_jaccard = _td_mod.compute_pairwise_jaccard
    compute_motif_collapse = _td_mod.compute_motif_collapse

    # 1. Jaccard divergence
    print("Computing pairwise Jaccard indices...")
    jaccard_result = compute_pairwise_jaccard(
        tool_beds,
        window_5prime=args.window_5prime,
        window_3prime=args.window_3prime,
    )

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.output_prefix

    if jaccard_result:
        # Write Jaccard metrics TSV first (downstream processes depend on it)
        import csv
        tsv_path = out_dir / f"{prefix}_divergence_metrics.tsv"
        flat_metrics = {k: v for k, v in jaccard_result.items()
                        if not isinstance(v, (dict, list))}
        with open(tsv_path, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(flat_metrics.keys())
            writer.writerow(flat_metrics.values())
        print(f"Wrote divergence metrics to {tsv_path}")

        plot_jaccard_heatmaps(
            jaccard_result,
            out_dir / f"{prefix}_jaccard_heatmaps.png",
            title_prefix=args.title_prefix, dpi=args.dpi,
        )
        plot_jaccard_comparison_bars(
            jaccard_result,
            out_dir / f"{prefix}_jaccard_comparison.png",
            title_prefix=args.title_prefix, dpi=args.dpi,
        )

    # 2. Motif collapse (for each tool individually, both 5' and 3')
    collapse_results_3p = {}
    collapse_results_5p = {}
    if args.reads_bed and args.genome and tool_maps:
        for tool, bed_path in tool_beds.items():
            map_path = tool_maps.get(tool)
            if map_path and map_path.exists():
                # 3' UTR motif collapse (PAS / ARE / miRNA)
                print(f"Computing 3' motif collapse for {tool}...")
                collapse_3p = compute_motif_collapse(
                    iso_bed=bed_path,
                    reads_bed=args.reads_bed,
                    map_file=map_path,
                    genome_path=args.genome,
                    window_3prime=args.window_3prime,
                    min_end_distance=args.min_end_distance,
                    end_type='tts',
                )
                collapse_results_3p[tool] = collapse_3p

                plot_motif_collapse_summary(
                    collapse_3p,
                    out_dir / f"{prefix}_{tool}_motif_collapse_3p.png",
                    title_prefix=f"{args.title_prefix}{tool}: ",
                    dpi=args.dpi,
                )

                # 5' RNA motif collapse (Inr / TOP)
                print(f"Computing 5' motif collapse for {tool}...")
                collapse_5p = compute_motif_collapse(
                    iso_bed=bed_path,
                    reads_bed=args.reads_bed,
                    map_file=map_path,
                    genome_path=args.genome,
                    window_3prime=args.window_3prime,
                    min_end_distance=args.min_end_distance,
                    end_type='tss',
                )
                collapse_results_5p[tool] = collapse_5p

                plot_motif_collapse_summary(
                    collapse_5p,
                    out_dir / f"{prefix}_{tool}_motif_collapse_5p.png",
                    title_prefix=f"{args.title_prefix}{tool}: ",
                    dpi=args.dpi,
                )

                # Side-by-side 5' vs 3' for this tool
                plot_motif_collapse_paired(
                    collapse_5p,
                    collapse_3p,
                    out_dir / f"{prefix}_{tool}_motif_collapse_paired.png",
                    title_prefix=f"{args.title_prefix}{tool}: ",
                    dpi=args.dpi,
                )

    # 3. Combined dashboard (per-tool, not pooled)
    if jaccard_result:
        plot_tool_divergence_dashboard(
            jaccard_result,
            collapse_results_3p if collapse_results_3p else None,
            collapse_results_5p if collapse_results_5p else None,
            out_dir / f"{prefix}_divergence_dashboard.png",
            title_prefix=args.title_prefix, dpi=args.dpi,
        )

    print("Done.")


if __name__ == '__main__':
    main()
