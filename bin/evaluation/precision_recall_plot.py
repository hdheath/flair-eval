#!/usr/bin/env python3
"""
Precision-Recall summary: individual plots for cross-mode comparison.

  - 5' P/R scatter (coloured by mode)
  - 3' P/R scatter (coloured by mode)
  - Isoforms per gene (grouped bar)
  - F1 score comparison (5'/3' bars)

Publication-quality formatting via ``pub_style``.
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from pub_style import ModeStyler, style_ax, legend_outside, savefig


def _short_mode(name: str) -> str:
    return (
        name.replace("density-", "d-")
        .replace("_default", "")
        .replace("_pacbio", "")
        .replace("end-scoring-", "es-")
    )



def load_evaluation_files(input_files):
    dfs = []
    for f in input_files:
        try:
            df = pd.read_csv(f, sep='\t')
            df['source_file'] = str(f)
            dfs.append(df)
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}", file=sys.stderr)
    if not dfs:
        return None
    return pd.concat(dfs, ignore_index=True)


def _safe_numeric(df, col):
    if col not in df.columns:
        return pd.Series([np.nan] * len(df), index=df.index)
    return pd.to_numeric(df[col], errors="coerce")


def _resolve_mode_order(df):
    default_order = [
        "default", "density-asymmetric", "density-asymmetric-softclip",
        "density-plain", "density-strict", "k-means", "more-ends",
        "bambu_default", "isoquant_pacbio",
    ]
    present = set(df["transcriptome_mode"].unique())
    order = [m for m in default_order if m in present]
    for m in sorted(present):
        if m not in order:
            order.append(m)
    return order


def _resolve_baseline(mode_order, baseline_arg):
    if not baseline_arg:
        return None
    modes_set = set(mode_order)
    if baseline_arg in modes_set:
        return baseline_arg
    if "baseline" in modes_set:
        return "baseline"
    prefix_matches = [m for m in mode_order if m.startswith(baseline_arg)]
    if len(prefix_matches) == 1:
        return prefix_matches[0]
    contains_matches = [m for m in mode_order if baseline_arg in m]
    if len(contains_matches) == 1:
        return contains_matches[0]
    return None


# ---------------------------------------------------------------------------
# Individual plot functions
# ---------------------------------------------------------------------------

def _plot_pr_scatter(df, output_path, mode_order, styler, baseline_mode,
                     end_type, end_label, title_prefix=""):
    p_col = f'{end_type}_precision'
    r_col = f'{end_type}_recall'

    fig, ax = plt.subplots(figsize=(3.5, 3.5))
    for _, row in df.iterrows():
        m = str(row.get('transcriptome_mode', 'unknown'))
        p_val = pd.to_numeric(row.get(p_col), errors='coerce')
        r_val = pd.to_numeric(row.get(r_col), errors='coerce')
        if pd.isna(p_val) or pd.isna(r_val):
            continue
        ax.scatter(r_val * 100, p_val * 100,
                   s=55, c=styler.color(m), marker=styler.marker(m),
                   edgecolors='white', linewidth=0.5, alpha=0.9, zorder=2)

    style_ax(ax, xlabel=f"{end_label} Recall (%)", ylabel=f"{end_label} Precision (%)",
             title=f"{title_prefix}{end_label} Precision vs Recall")
    ax.set_xlim(0, 105)
    ax.set_ylim(0, 105)
    ax.plot([0, 100], [0, 100], "--", color="#cccccc", linewidth=0.8, zorder=0)

    if baseline_mode:
        bl_sub = df[df['transcriptome_mode'] == baseline_mode]
        bl_p = _safe_numeric(bl_sub, p_col).mean()
        bl_r = _safe_numeric(bl_sub, r_col).mean()
        if not (np.isnan(bl_r) or np.isnan(bl_p)):
            ax.axhline(bl_p * 100, color='#888888', linestyle=':', linewidth=1.0, zorder=1, alpha=0.7)
            ax.axvline(bl_r * 100, color='#888888', linestyle=':', linewidth=1.0, zorder=1, alpha=0.7)

    handles = [styler.legend_handle(m, label=m.replace("_", " "), markersize=7)
               for m in mode_order]
    legend_outside(fig, handles=handles, loc='upper left', bbox_to_anchor=(1.02, 1.0),
                   ncol=1, fontsize=7)
    fig.tight_layout()
    savefig(fig, output_path, dpi=300)


def _plot_f1_bars(df, output_path, mode_order, styler, baseline_mode, title_prefix=""):
    x = np.arange(len(mode_order))
    f1_5p, f1_3p = [], []
    for mode in mode_order:
        subset = df[df['transcriptome_mode'] == mode]
        p5 = _safe_numeric(subset, '5prime_precision').mean()
        r5 = _safe_numeric(subset, '5prime_recall').mean()
        p3 = _safe_numeric(subset, '3prime_precision').mean()
        r3 = _safe_numeric(subset, '3prime_recall').mean()
        f1_5 = 2 * p5 * r5 / (p5 + r5) if (p5 + r5) > 0 and not (np.isnan(p5) or np.isnan(r5)) else 0
        f1_3 = 2 * p3 * r3 / (p3 + r3) if (p3 + r3) > 0 and not (np.isnan(p3) or np.isnan(r3)) else 0
        f1_5p.append(f1_5 * 100)
        f1_3p.append(f1_3 * 100)

    bar_width = 0.35
    fig, ax = plt.subplots(figsize=(max(3.5, len(mode_order) * 0.7), 3.0))
    for i, mode in enumerate(mode_order):
        c = styler.color(mode)
        ax.bar(x[i] - bar_width / 2, f1_5p[i], width=bar_width * 0.9,
               color=c, alpha=0.9, edgecolor='white', linewidth=0.4)
        ax.bar(x[i] + bar_width / 2, f1_3p[i], width=bar_width * 0.9,
               color=c, alpha=0.55, edgecolor='white', linewidth=0.4, hatch='///')
        if f1_5p[i] > 0:
            ax.text(x[i] - bar_width / 2, f1_5p[i] + 0.5, f"{f1_5p[i]:.1f}",
                    ha='center', va='bottom', fontsize=6, color='#333333')
        if f1_3p[i] > 0:
            ax.text(x[i] + bar_width / 2, f1_3p[i] + 0.5, f"{f1_3p[i]:.1f}",
                    ha='center', va='bottom', fontsize=6, color='#333333')

    if baseline_mode:
        bl_sub = df[df['transcriptome_mode'] == baseline_mode]
        bl_p5 = _safe_numeric(bl_sub, '5prime_precision').mean()
        bl_r5 = _safe_numeric(bl_sub, '5prime_recall').mean()
        bl_p3 = _safe_numeric(bl_sub, '3prime_precision').mean()
        bl_r3 = _safe_numeric(bl_sub, '3prime_recall').mean()
        if not (np.isnan(bl_p5) or np.isnan(bl_r5)) and (bl_p5 + bl_r5) > 0:
            bl_f1_5 = 2 * bl_p5 * bl_r5 / (bl_p5 + bl_r5) * 100
            ax.axhline(bl_f1_5, color='#555555', linestyle=':', linewidth=1.0, alpha=0.6, zorder=1,
                       label=f"5' baseline ({bl_f1_5:.1f})")
        if not (np.isnan(bl_p3) or np.isnan(bl_r3)) and (bl_p3 + bl_r3) > 0:
            bl_f1_3 = 2 * bl_p3 * bl_r3 / (bl_p3 + bl_r3) * 100
            ax.axhline(bl_f1_3, color='#555555', linestyle='--', linewidth=1.0, alpha=0.6, zorder=1,
                       label=f"3' baseline ({bl_f1_3:.1f})")

    mode_handles = [styler.legend_handle(m, label=_short_mode(m), markersize=6)
                    for m in mode_order]
    hatch_handles = [
        mpatches.Patch(facecolor='#888888', alpha=0.9, edgecolor='none', label="5\u2032 F1"),
        mpatches.Patch(facecolor='#888888', alpha=0.55, edgecolor='none', hatch='///', label="3\u2032 F1"),
    ]
    ax.set_xticks(x)
    ax.set_xticklabels([_short_mode(m) for m in mode_order], rotation=45, ha='right', fontsize=7)
    ax.set_ylim(0, 105)
    legend_outside(fig, handles=mode_handles + hatch_handles,
                   loc='upper left', bbox_to_anchor=(1.02, 1.0), ncol=1, fontsize=7)
    style_ax(ax, ylabel="F1 Score (%)", faint_y_grid=True)
    fig.tight_layout()
    savefig(fig, output_path, dpi=300)


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------

def create_precision_recall_plots(df, output_dir, title_prefix="", baseline=None):
    """Create individual precision-recall, isoforms/gene, and F1 plots."""
    mode_order = _resolve_mode_order(df)
    if not mode_order:
        print("No modes found - skipping precision/recall plots", file=sys.stderr)
        return False

    styler = ModeStyler(mode_order)
    baseline_mode = _resolve_baseline(mode_order, baseline)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    _plot_pr_scatter(df, output_dir / "pr_5prime_scatter.png", mode_order, styler,
                     baseline_mode, "5prime", "5\u2032 TSS", title_prefix)
    _plot_pr_scatter(df, output_dir / "pr_3prime_scatter.png", mode_order, styler,
                     baseline_mode, "3prime", "3\u2032 TTS", title_prefix)
    _plot_f1_bars(df, output_dir / "f1_score.png", mode_order, styler,
                  baseline_mode, title_prefix)
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Create individual precision-recall, isoforms/gene, and F1 plots"
    )
    parser.add_argument('--input', '-i', nargs='+', required=True,
                        help="Input evaluation TSV file(s)")
    parser.add_argument('--output', '-o', required=True,
                        help="Output directory (individual PNGs will be saved here)")
    parser.add_argument('--title-prefix', default="",
                        help="Optional prefix for the plot title")
    parser.add_argument('--baseline', default=None,
                        help="Baseline mode name for reference lines")
    parser.add_argument('--verbose', '-v', action='store_true')

    args = parser.parse_args()

    df = load_evaluation_files(args.input)
    if df is None or len(df) == 0:
        print("Error: No valid data found in input files", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        modes = df['transcriptome_mode'].unique()
        print(f"Loaded {len(df)} rows, {len(modes)} modes: {list(modes)}")

    success = create_precision_recall_plots(df, args.output, args.title_prefix,
                                            baseline=args.baseline)
    if success:
        print(f"Saved precision-recall plots to {args.output}")
    else:
        print("Failed to create plots", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
