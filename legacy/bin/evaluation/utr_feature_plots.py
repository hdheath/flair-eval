#!/usr/bin/env python3
"""
5'UTR and 3'UTR feature comparison plots.

Creates multi-panel figures comparing assembler isoform UTR features
(length, GC content, uORF count) against the reference annotation.

Usage:
    python utr_feature_plots.py \\
        --isoform-beds tool1:path1.bed tool2:path2.bed ... \\
        --gtf reference.gtf \\
        --genome genome.fa \\
        --output-prefix test_name \\
        --output-dir . \\
        [--title-prefix "Test: "] [--dpi 300] [--verbose]

Can also be imported and called programmatically.
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from pub_style import style_ax, savefig, PALETTE


def plot_utr_length_comparison(
    comparison_results: Dict[str, dict],
    output_path: Path,
    title_prefix: str = "",
    dpi: int = 300,
) -> bool:
    """Violin/box plots of 5'UTR length: each assembler vs reference.

    One panel per assembler, plus the reference distribution.
    """
    if not HAS_MATPLOTLIB:
        return False

    tools = sorted(comparison_results.keys())
    if not tools:
        return False

    n_panels = len(tools)
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5.5), squeeze=False)
    # (no suptitle — pub-quality)

    for i, tool in enumerate(tools):
        ax = axes[0][i]
        comp = comparison_results[tool]
        asm_lengths = comp.get('_assembler_lengths', [])
        ref_lengths = comp.get('_reference_lengths', [])

        if not asm_lengths and not ref_lengths:
            ax.text(0.5, 0.5, 'No data', transform=ax.transAxes,
                    ha='center', va='center')
            continue

        data = []
        labels = []
        colors = []
        if ref_lengths:
            data.append(ref_lengths)
            labels.append('Reference')
            colors.append('#95a5a6')
        if asm_lengths:
            data.append(asm_lengths)
            labels.append(tool)
            colors.append('#3498db')

        vp = ax.violinplot(data, positions=range(len(data)),
                           showextrema=False, showmedians=False)
        for j, body in enumerate(vp['bodies']):
            body.set_facecolor(colors[j])
            body.set_alpha(0.6)

        bp = ax.boxplot(data, positions=range(len(data)),
                        widths=0.15, patch_artist=False,
                        medianprops=dict(color='red', linewidth=2),
                        whiskerprops=dict(linewidth=0.8),
                        flierprops=dict(markersize=2, alpha=0.3))

        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, fontsize=7)
        style_ax(ax, ylabel='5\u2032UTR Length (bp)', title=tool)

        # Cap y-axis at 99th percentile for readability
        all_vals = asm_lengths + ref_lengths
        if all_vals:
            p99 = sorted(all_vals)[int(len(all_vals) * 0.99)]
            ax.set_ylim(0, min(p99 * 1.2, max(all_vals)))

        # Stats annotation
        asm_med = comp.get('assembler_utr5_length_median', 0)
        ref_med = comp.get('reference_utr5_length_median', 0)
        delta = comp.get('utr5_length_median_delta', 0)
        n_asm = comp.get('assembler_n_isoforms_with_utr', 0)
        n_ref = comp.get('reference_n_transcripts_with_utr', 0)
        txt = (f"Ref median: {ref_med:.0f} bp (n={n_ref})\n"
               f"Asm median: {asm_med:.0f} bp (n={n_asm})\n"
               f"\u0394 median: {delta:+.0f} bp")
        ax.text(0.98, 0.98, txt, transform=ax.transAxes, fontsize=8,
                va='top', ha='right',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.7))

    fig.tight_layout()
    savefig(fig, output_path, dpi=dpi)
    print(f"Saved 5'UTR length comparison to {output_path}")
    return True


def plot_utr_gc_comparison(
    comparison_results: Dict[str, dict],
    output_path: Path,
    title_prefix: str = "",
    dpi: int = 300,
) -> bool:
    """Histogram overlay of 5'UTR GC content: each assembler vs reference."""
    if not HAS_MATPLOTLIB:
        return False

    tools = sorted(comparison_results.keys())
    if not tools:
        return False

    n_panels = len(tools)
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5), squeeze=False)
    # (no suptitle — pub-quality)

    for i, tool in enumerate(tools):
        ax = axes[0][i]
        comp = comparison_results[tool]
        asm_gc = comp.get('_assembler_gc', [])
        ref_gc = comp.get('_reference_gc', [])

        bins = np.linspace(0, 1, 41)

        if ref_gc:
            ax.hist(ref_gc, bins=bins, alpha=0.5, color=PALETTE[7],
                    label='Reference', density=True, edgecolor='none')
        if asm_gc:
            ax.hist(asm_gc, bins=bins, alpha=0.5, color=PALETTE[5],
                    label=tool, density=True, edgecolor='none')

        style_ax(ax, xlabel='GC Content', ylabel='Density', title=tool)
        ax.legend(fontsize=8, frameon=False)

        asm_mean = comp.get('assembler_utr5_gc_mean', 0)
        ref_mean = comp.get('reference_utr5_gc_mean', 0)
        delta = comp.get('utr5_gc_mean_delta', 0)
        txt = (f"Ref GC: {ref_mean:.3f}\n"
               f"Asm GC: {asm_mean:.3f}\n"
               f"\u0394: {delta:+.3f}")
        ax.text(0.98, 0.98, txt, transform=ax.transAxes, fontsize=8,
                va='top', ha='right',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.7))

    fig.tight_layout()
    savefig(fig, output_path, dpi=dpi)
    print(f"Saved 5'UTR GC comparison to {output_path}")
    return True


def plot_utr_uorf_comparison(
    comparison_results: Dict[str, dict],
    output_path: Path,
    title_prefix: str = "",
    dpi: int = 300,
) -> bool:
    """Bar chart of uORF counts: assembler vs reference distribution."""
    if not HAS_MATPLOTLIB:
        return False

    tools = sorted(comparison_results.keys())
    if not tools:
        return False

    n_panels = len(tools)
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5), squeeze=False)
    # (no suptitle — pub-quality)

    for i, tool in enumerate(tools):
        ax = axes[0][i]
        comp = comparison_results[tool]
        asm_uorfs = comp.get('_assembler_uorfs', [])
        ref_uorfs = comp.get('_reference_uorfs', [])

        max_uorf = 15  # cap visualization
        bins = range(0, max_uorf + 2)

        if ref_uorfs:
            ref_capped = [min(u, max_uorf) for u in ref_uorfs]
            ax.hist(ref_capped, bins=bins, alpha=0.5, color=PALETTE[7],
                    label='Reference', density=True, edgecolor='none',
                    rwidth=0.85, align='left')
        if asm_uorfs:
            asm_capped = [min(u, max_uorf) for u in asm_uorfs]
            ax.hist(asm_capped, bins=bins, alpha=0.5, color=PALETTE[2],
                    label=tool, density=True, edgecolor='none',
                    rwidth=0.85, align='left')

        style_ax(ax, xlabel='uORF Count (NUG-initiated)', ylabel='Density', title=tool)
        ax.legend(fontsize=8, frameon=False)

        from matplotlib.ticker import MaxNLocator
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

        asm_mean = comp.get('assembler_utr5_uorf_mean', 0)
        ref_mean = comp.get('reference_utr5_uorf_mean', 0)
        delta = comp.get('utr5_uorf_mean_delta', 0)
        txt = (f"Ref mean: {ref_mean:.1f}\n"
               f"Asm mean: {asm_mean:.1f}\n"
               f"\u0394: {delta:+.1f}")
        ax.text(0.98, 0.98, txt, transform=ax.transAxes, fontsize=8,
                va='top', ha='right',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.7))

    fig.tight_layout()
    savefig(fig, output_path, dpi=dpi)
    print(f"Saved 5'UTR uORF comparison to {output_path}")
    return True


def plot_utr_feature_dashboard(
    comparison_results: Dict[str, dict],
    output_path: Path,
    title_prefix: str = "",
    dpi: int = 300,
) -> bool:
    """Combined 3-panel dashboard: length + GC + uORF across all assemblers.

    Layout: 1 × 3
      Panel 1: Grouped bars of median 5'UTR length (ref vs each assembler)
      Panel 2: Grouped bars of mean GC content
      Panel 3: Grouped bars of mean uORF count
    """
    if not HAS_MATPLOTLIB:
        return False

    tools = sorted(comparison_results.keys())
    if not tools:
        return False

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5.5))

    x = np.arange(len(tools))
    width = 0.35

    # ── Panel 1: median 5'UTR length ──
    ref_lens = [comparison_results[t].get('reference_utr5_length_median', 0) for t in tools]
    asm_lens = [comparison_results[t].get('assembler_utr5_length_median', 0) for t in tools]

    ax1.bar(x - width / 2, ref_lens, width, label='Reference',
            color=PALETTE[7], edgecolor='none', alpha=0.85)
    ax1.bar(x + width / 2, asm_lens, width, label='Assembler',
            color=PALETTE[1], edgecolor='none', alpha=0.85)
    ax1.set_xticks(x)
    ax1.set_xticklabels(tools, fontsize=7, rotation=30, ha='right')
    style_ax(ax1, ylabel='Median 5\u2032UTR Length (bp)', title='5\u2032UTR Length')
    ax1.legend(fontsize=8, frameon=False)

    # Value labels
    for bars in [ref_lens, asm_lens]:
        offset = -width / 2 if bars is ref_lens else width / 2
        for j, v in enumerate(bars):
            ax1.text(x[j] + offset, v + max(max(ref_lens + asm_lens) * 0.02, 1),
                     f'{v:.0f}', ha='center', va='bottom', fontsize=7)

    # ── Panel 2: mean GC content ──
    ref_gcs = [comparison_results[t].get('reference_utr5_gc_mean', 0) for t in tools]
    asm_gcs = [comparison_results[t].get('assembler_utr5_gc_mean', 0) for t in tools]

    ax2.bar(x - width / 2, ref_gcs, width, label='Reference',
            color=PALETTE[7], edgecolor='none', alpha=0.85)
    ax2.bar(x + width / 2, asm_gcs, width, label='Assembler',
            color=PALETTE[5], edgecolor='none', alpha=0.85)
    ax2.set_xticks(x)
    ax2.set_xticklabels(tools, fontsize=7, rotation=30, ha='right')
    style_ax(ax2, ylabel='Mean GC Content', title='5\u2032UTR GC Content')
    ax2.legend(fontsize=8, frameon=False)
    ax2.set_ylim(0, max(max(ref_gcs + asm_gcs, default=0.5) * 1.2, 0.1))

    for bars in [ref_gcs, asm_gcs]:
        offset = -width / 2 if bars is ref_gcs else width / 2
        for j, v in enumerate(bars):
            ax2.text(x[j] + offset, v + max(max(ref_gcs + asm_gcs) * 0.02, 0.005),
                     f'{v:.3f}', ha='center', va='bottom', fontsize=7)

    # ── Panel 3: mean uORF count ──
    ref_uorfs = [comparison_results[t].get('reference_utr5_uorf_mean', 0) for t in tools]
    asm_uorfs = [comparison_results[t].get('assembler_utr5_uorf_mean', 0) for t in tools]

    ax3.bar(x - width / 2, ref_uorfs, width, label='Reference',
            color=PALETTE[7], edgecolor='none', alpha=0.85)
    ax3.bar(x + width / 2, asm_uorfs, width, label='Assembler',
            color=PALETTE[2], edgecolor='none', alpha=0.85)
    ax3.set_xticks(x)
    ax3.set_xticklabels(tools, fontsize=7, rotation=30, ha='right')
    style_ax(ax3, ylabel='Mean uORF Count', title='5\u2032UTR uORFs (NUG-initiated)')
    ax3.legend(fontsize=8, frameon=False)

    for bars in [ref_uorfs, asm_uorfs]:
        offset = -width / 2 if bars is ref_uorfs else width / 2
        for j, v in enumerate(bars):
            ax3.text(x[j] + offset, v + max(max(ref_uorfs + asm_uorfs) * 0.02, 0.05),
                     f'{v:.1f}', ha='center', va='bottom', fontsize=7)

    fig.tight_layout()
    savefig(fig, output_path, dpi=dpi)
    print(f"Saved 5'UTR feature dashboard to {output_path}")
    return True


# ---------------------------------------------------------------------------
# 3'UTR plotting functions
# ---------------------------------------------------------------------------

def plot_utr3_length_comparison(
    comparison_results: Dict[str, dict],
    output_path: Path,
    title_prefix: str = "",
    dpi: int = 300,
) -> bool:
    """Violin/box plots of 3'UTR length: each assembler vs reference."""
    if not HAS_MATPLOTLIB:
        return False

    tools = sorted(comparison_results.keys())
    if not tools:
        return False

    n_panels = len(tools)
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5.5), squeeze=False)

    for i, tool in enumerate(tools):
        ax = axes[0][i]
        comp = comparison_results[tool]
        asm_lengths = comp.get('_assembler_lengths', [])
        ref_lengths = comp.get('_reference_lengths', [])

        if not asm_lengths and not ref_lengths:
            ax.text(0.5, 0.5, 'No data', transform=ax.transAxes,
                    ha='center', va='center')
            continue

        data = []
        labels = []
        colors = []
        if ref_lengths:
            data.append(ref_lengths)
            labels.append('Reference')
            colors.append('#95a5a6')
        if asm_lengths:
            data.append(asm_lengths)
            labels.append(tool)
            colors.append('#e74c3c')

        vp = ax.violinplot(data, positions=range(len(data)),
                           showextrema=False, showmedians=False)
        for j, body in enumerate(vp['bodies']):
            body.set_facecolor(colors[j])
            body.set_alpha(0.6)

        bp = ax.boxplot(data, positions=range(len(data)),
                        widths=0.15, patch_artist=False,
                        medianprops=dict(color='red', linewidth=2),
                        whiskerprops=dict(linewidth=0.8),
                        flierprops=dict(markersize=2, alpha=0.3))

        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, fontsize=7)
        style_ax(ax, ylabel='3\u2032UTR Length (bp)', title=tool)

        all_vals = asm_lengths + ref_lengths
        if all_vals:
            p99 = sorted(all_vals)[int(len(all_vals) * 0.99)]
            ax.set_ylim(0, min(p99 * 1.2, max(all_vals)))

        asm_med = comp.get('assembler_utr3_length_median', 0)
        ref_med = comp.get('reference_utr3_length_median', 0)
        delta = comp.get('utr3_length_median_delta', 0)
        n_asm = comp.get('assembler_n_isoforms_with_utr3', 0)
        n_ref = comp.get('reference_n_transcripts_with_utr3', 0)
        txt = (f"Ref median: {ref_med:.0f} bp (n={n_ref})\n"
               f"Asm median: {asm_med:.0f} bp (n={n_asm})\n"
               f"\u0394 median: {delta:+.0f} bp")
        ax.text(0.98, 0.98, txt, transform=ax.transAxes, fontsize=8,
                va='top', ha='right',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.7))

    fig.tight_layout()
    savefig(fig, output_path, dpi=dpi)
    print(f"Saved 3'UTR length comparison to {output_path}")
    return True


def plot_utr3_gc_comparison(
    comparison_results: Dict[str, dict],
    output_path: Path,
    title_prefix: str = "",
    dpi: int = 300,
) -> bool:
    """Histogram overlay of 3'UTR GC content: each assembler vs reference."""
    if not HAS_MATPLOTLIB:
        return False

    tools = sorted(comparison_results.keys())
    if not tools:
        return False

    n_panels = len(tools)
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5), squeeze=False)

    for i, tool in enumerate(tools):
        ax = axes[0][i]
        comp = comparison_results[tool]
        asm_gc = comp.get('_assembler_gc', [])
        ref_gc = comp.get('_reference_gc', [])

        bins = np.linspace(0, 1, 41)

        if ref_gc:
            ax.hist(ref_gc, bins=bins, alpha=0.5, color=PALETTE[7],
                    label='Reference', density=True, edgecolor='none')
        if asm_gc:
            ax.hist(asm_gc, bins=bins, alpha=0.5, color=PALETTE[5],
                    label=tool, density=True, edgecolor='none')

        style_ax(ax, xlabel='GC Content', ylabel='Density', title=tool)
        ax.legend(fontsize=8, frameon=False)

        asm_mean = comp.get('assembler_utr3_gc_mean', 0)
        ref_mean = comp.get('reference_utr3_gc_mean', 0)
        delta = comp.get('utr3_gc_mean_delta', 0)
        txt = (f"Ref GC: {ref_mean:.3f}\n"
               f"Asm GC: {asm_mean:.3f}\n"
               f"\u0394: {delta:+.3f}")
        ax.text(0.98, 0.98, txt, transform=ax.transAxes, fontsize=8,
                va='top', ha='right',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.7))

    fig.tight_layout()
    savefig(fig, output_path, dpi=dpi)
    print(f"Saved 3'UTR GC comparison to {output_path}")
    return True


def plot_utr3_feature_dashboard(
    comparison_results: Dict[str, dict],
    output_path: Path,
    title_prefix: str = "",
    dpi: int = 300,
) -> bool:
    """Combined 2-panel dashboard: 3'UTR length + GC across all assemblers.

    Layout: 1 × 2
      Panel 1: Grouped bars of median 3'UTR length (ref vs each assembler)
      Panel 2: Grouped bars of mean GC content
    """
    if not HAS_MATPLOTLIB:
        return False

    tools = sorted(comparison_results.keys())
    if not tools:
        return False

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5.5))

    x = np.arange(len(tools))
    width = 0.35

    # ── Panel 1: median 3'UTR length ──
    ref_lens = [comparison_results[t].get('reference_utr3_length_median', 0) for t in tools]
    asm_lens = [comparison_results[t].get('assembler_utr3_length_median', 0) for t in tools]

    ax1.bar(x - width / 2, ref_lens, width, label='Reference',
            color=PALETTE[7], edgecolor='none', alpha=0.85)
    ax1.bar(x + width / 2, asm_lens, width, label='Assembler',
            color=PALETTE[5], edgecolor='none', alpha=0.85)
    ax1.set_xticks(x)
    ax1.set_xticklabels(tools, fontsize=7, rotation=30, ha='right')
    style_ax(ax1, ylabel='Median 3\u2032UTR Length (bp)', title='3\u2032UTR Length')
    ax1.legend(fontsize=8, frameon=False)

    for bars in [ref_lens, asm_lens]:
        offset = -width / 2 if bars is ref_lens else width / 2
        for j, v in enumerate(bars):
            ax1.text(x[j] + offset, v + max(max(ref_lens + asm_lens, default=1) * 0.02, 1),
                     f'{v:.0f}', ha='center', va='bottom', fontsize=7)

    # ── Panel 2: mean GC content ──
    ref_gcs = [comparison_results[t].get('reference_utr3_gc_mean', 0) for t in tools]
    asm_gcs = [comparison_results[t].get('assembler_utr3_gc_mean', 0) for t in tools]

    ax2.bar(x - width / 2, ref_gcs, width, label='Reference',
            color=PALETTE[7], edgecolor='none', alpha=0.85)
    ax2.bar(x + width / 2, asm_gcs, width, label='Assembler',
            color=PALETTE[5], edgecolor='none', alpha=0.85)
    ax2.set_xticks(x)
    ax2.set_xticklabels(tools, fontsize=7, rotation=30, ha='right')
    style_ax(ax2, ylabel='Mean GC Content', title='3\u2032UTR GC Content')
    ax2.legend(fontsize=8, frameon=False)
    ax2.set_ylim(0, max(max(ref_gcs + asm_gcs, default=0.5) * 1.2, 0.1))

    for bars in [ref_gcs, asm_gcs]:
        offset = -width / 2 if bars is ref_gcs else width / 2
        for j, v in enumerate(bars):
            ax2.text(x[j] + offset, v + max(max(ref_gcs + asm_gcs, default=0.1) * 0.02, 0.005),
                     f'{v:.3f}', ha='center', va='bottom', fontsize=7)

    fig.tight_layout()
    savefig(fig, output_path, dpi=dpi)
    print(f"Saved 3'UTR feature dashboard to {output_path}")
    return True


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compute and plot 5'UTR and 3'UTR feature comparisons")

    parser.add_argument('--isoform-beds', nargs='+', required=True,
                        help='tool:path pairs, e.g. flair:isoforms.bed bambu:isoforms.bed')
    parser.add_argument('--gtf', type=Path, required=True,
                        help='Reference annotation GTF (with CDS features)')
    parser.add_argument('--genome', type=Path, required=True,
                        help='Indexed genome FASTA')
    parser.add_argument('--output-prefix', type=str, required=True,
                        help='Prefix for output files')
    parser.add_argument('--output-dir', type=Path, default=Path('.'),
                        help='Output directory')
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

    # Import the UTR features module via importlib to avoid heavy __init__.py deps
    import importlib.util
    _eval_dir = Path(__file__).resolve().parent
    for _mod_name in ('utils',):
        _spec = importlib.util.spec_from_file_location(
            f'evaluation.{_mod_name}', _eval_dir / f'{_mod_name}.py',
            submodule_search_locations=[])
        _mod = importlib.util.module_from_spec(_spec)
        sys.modules[f'evaluation.{_mod_name}'] = _mod
        _spec.loader.exec_module(_mod)
    if 'evaluation' not in sys.modules:
        import types
        _pkg = types.ModuleType('evaluation')
        _pkg.__path__ = [str(_eval_dir)]
        sys.modules['evaluation'] = _pkg
    _spec = importlib.util.spec_from_file_location(
        'evaluation.utr_features', _eval_dir / 'utr_features.py',
        submodule_search_locations=[])
    _utr_mod = importlib.util.module_from_spec(_spec)
    sys.modules['evaluation.utr_features'] = _utr_mod
    _spec.loader.exec_module(_utr_mod)
    extract_5utr_features = _utr_mod.extract_5utr_features
    extract_reference_5utr_features = _utr_mod.extract_reference_5utr_features
    compare_utr_features = _utr_mod.compare_utr_features
    extract_3utr_features = _utr_mod.extract_3utr_features
    extract_reference_3utr_features = _utr_mod.extract_reference_3utr_features
    compare_3utr_features = _utr_mod.compare_3utr_features

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.output_prefix

    # ── 5'UTR analysis ──────────────────────────────────────────────────
    print("Extracting reference 5'UTR features...")
    ref_5utr = extract_reference_5utr_features(args.gtf, args.genome)
    if not ref_5utr:
        print("Warning: No reference 5'UTR features found (no CDS in GTF?)",
              file=sys.stderr)

    comparison_5utr: Dict[str, dict] = {}
    for tool, bed_path in sorted(tool_beds.items()):
        print(f"Extracting 5'UTR features for {tool}...")
        asm_features = extract_5utr_features(
            iso_bed=bed_path, genome_path=args.genome, gtf_path=args.gtf)
        print(f"Comparing {tool} 5'UTR vs reference...")
        comparison_5utr[tool] = compare_utr_features(asm_features, ref_5utr)

    # ── 3'UTR analysis ──────────────────────────────────────────────────
    print("Extracting reference 3'UTR features...")
    ref_3utr = extract_reference_3utr_features(args.gtf, args.genome)
    if not ref_3utr:
        print("Warning: No reference 3'UTR features found (no CDS in GTF?)",
              file=sys.stderr)

    comparison_3utr: Dict[str, dict] = {}
    for tool, bed_path in sorted(tool_beds.items()):
        print(f"Extracting 3'UTR features for {tool}...")
        asm_features = extract_3utr_features(
            iso_bed=bed_path, genome_path=args.genome, gtf_path=args.gtf)
        print(f"Comparing {tool} 3'UTR vs reference...")
        comparison_3utr[tool] = compare_3utr_features(asm_features, ref_3utr)

    # ── Write combined metrics TSV ──────────────────────────────────────
    if comparison_5utr or comparison_3utr:
        tsv_path = out_dir / f"{prefix}_utr_features.tsv"
        with open(tsv_path, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            # Merge 5' and 3' keys
            all_keys = []
            if comparison_5utr:
                first_5 = next(iter(comparison_5utr.values()))
                all_keys += [k for k in first_5 if not k.startswith('_')]
            if comparison_3utr:
                first_3 = next(iter(comparison_3utr.values()))
                all_keys += [k for k in first_3 if not k.startswith('_')]
            writer.writerow(['tool'] + all_keys)
            all_tools = sorted(set(list(comparison_5utr.keys()) + list(comparison_3utr.keys())))
            for tool in all_tools:
                row_vals = [tool]
                comp5 = comparison_5utr.get(tool, {})
                comp3 = comparison_3utr.get(tool, {})
                merged = {**comp5, **comp3}
                row_vals += [merged.get(k, '') for k in all_keys]
                writer.writerow(row_vals)
        print(f"Wrote UTR feature metrics to {tsv_path}")

    # ── Generate 5'UTR plots ───────────────────────────────────────────
    if comparison_5utr:
        plot_utr_length_comparison(
            comparison_5utr, out_dir / f"{prefix}_utr5_length.png",
            title_prefix=args.title_prefix, dpi=args.dpi)
        plot_utr_gc_comparison(
            comparison_5utr, out_dir / f"{prefix}_utr5_gc.png",
            title_prefix=args.title_prefix, dpi=args.dpi)
        plot_utr_uorf_comparison(
            comparison_5utr, out_dir / f"{prefix}_utr5_uorfs.png",
            title_prefix=args.title_prefix, dpi=args.dpi)
        plot_utr_feature_dashboard(
            comparison_5utr, out_dir / f"{prefix}_utr5_dashboard.png",
            title_prefix=args.title_prefix, dpi=args.dpi)

    # ── Generate 3'UTR plots ───────────────────────────────────────────
    if comparison_3utr:
        plot_utr3_length_comparison(
            comparison_3utr, out_dir / f"{prefix}_utr3_length.png",
            title_prefix=args.title_prefix, dpi=args.dpi)
        plot_utr3_gc_comparison(
            comparison_3utr, out_dir / f"{prefix}_utr3_gc.png",
            title_prefix=args.title_prefix, dpi=args.dpi)
        plot_utr3_feature_dashboard(
            comparison_3utr, out_dir / f"{prefix}_utr3_dashboard.png",
            title_prefix=args.title_prefix, dpi=args.dpi)

    print("Done.")


if __name__ == '__main__':
    main()
