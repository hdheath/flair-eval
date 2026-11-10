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


def _plot_paired_summary(df, output_dir, mode_order, styler, baseline_mode, title_prefix=""):
    """Three paired-end plots: paired precision, paired recall (geom mean), paired F1.

    Paired precision  = paired_dedup_precision (JC-unique (TSS-peak, TTS-peak) pairs where
                        both ends hit a peak / total isoforms emitted by the method)
    Paired recall     = sqrt(5prime_recall * 3prime_recall)  [geometric mean]
    Paired F1         = harmonic mean of paired precision and paired recall
    """
    output_dir = Path(output_dir)

    paired_prec, paired_rec, paired_f1 = [], [], []
    for mode in mode_order:
        subset = df[df['transcriptome_mode'] == mode]
        r5 = _safe_numeric(subset, '5prime_recall').mean()
        r3 = _safe_numeric(subset, '3prime_recall').mean()
        pp = _safe_numeric(subset, 'paired_dedup_precision').mean()

        if not (np.isnan(pp)):
            paired_prec.append(pp * 100)
        else:
            paired_prec.append(np.nan)

        if not (np.isnan(r5) or np.isnan(r3)) and r5 >= 0 and r3 >= 0:
            pr = np.sqrt(r5 * r3)
            paired_rec.append(pr * 100)
        else:
            paired_rec.append(np.nan)
            pr = np.nan

        pp_val = pp if not np.isnan(pp) else np.nan
        pr_val = np.sqrt(r5 * r3) if (not (np.isnan(r5) or np.isnan(r3)) and r5 >= 0 and r3 >= 0) else np.nan
        if not (np.isnan(pp_val) or np.isnan(pr_val)) and (pp_val + pr_val) > 0:
            paired_f1.append(2 * pp_val * pr_val / (pp_val + pr_val) * 100)
        else:
            paired_f1.append(np.nan)

    # Baseline values
    bl_prec = bl_rec = bl_f1 = np.nan
    if baseline_mode:
        bl_sub = df[df['transcriptome_mode'] == baseline_mode]
        bl_r5 = _safe_numeric(bl_sub, '5prime_recall').mean()
        bl_r3 = _safe_numeric(bl_sub, '3prime_recall').mean()
        bl_pp = _safe_numeric(bl_sub, 'paired_dedup_precision').mean()
        bl_prec = bl_pp * 100 if not np.isnan(bl_pp) else np.nan
        if not (np.isnan(bl_r5) or np.isnan(bl_r3)) and bl_r5 >= 0 and bl_r3 >= 0:
            bl_pr = np.sqrt(bl_r5 * bl_r3)
            bl_rec = bl_pr * 100
            if not np.isnan(bl_pp) and (bl_pp + bl_pr) > 0:
                bl_f1 = 2 * bl_pp * bl_pr / (bl_pp + bl_pr) * 100

    # Single grouped bar plot: precision / recall / F1 per mode
    x = np.arange(len(mode_order))
    bar_width = 0.25
    fig, ax = plt.subplots(figsize=(max(4.0, len(mode_order) * 0.9), 3.2))
    for i, mode in enumerate(mode_order):
        c = styler.color(mode)
        for j, (vals, alpha, hatch) in enumerate([
            (paired_prec, 0.90, None),
            (paired_rec,  0.55, '///'),
            (paired_f1,   0.75, 'xx'),
        ]):
            v = vals[i]
            offset = (j - 1) * bar_width
            bar = ax.bar(x[i] + offset, v if not np.isnan(v) else 0,
                         width=bar_width * 0.9, color=c, alpha=alpha,
                         edgecolor='white', linewidth=0.4,
                         hatch=hatch if hatch else '')
            if not np.isnan(v) and v > 0:
                ax.text(bar[0].get_x() + bar[0].get_width() / 2, v + 0.5,
                        f"{v:.1f}", ha='center', va='bottom', fontsize=5, color='#333333')

    for bl_val, ls in [(bl_prec, ':'), (bl_rec, '--'), (bl_f1, '-.')]:
        if not np.isnan(bl_val):
            ax.axhline(bl_val, color='#555555', linestyle=ls, linewidth=1.0, alpha=0.6)

    mode_handles = [
        mpatches.Patch(facecolor=styler.color(m), edgecolor='none',
                       label=_short_mode(m), alpha=0.9)
        for m in mode_order
    ]
    metric_handles = [
        mpatches.Patch(facecolor='#666666', alpha=0.90, edgecolor='none', label="Paired Precision"),
        mpatches.Patch(facecolor='#666666', alpha=0.55, edgecolor='none', hatch='///', label="Paired Recall (\u221ar5\u00b7r3)"),
        mpatches.Patch(facecolor='#666666', alpha=0.75, edgecolor='none', hatch='xx',  label="Paired F1"),
    ]
    ax.set_xticks(x)
    ax.set_xticklabels([_short_mode(m) for m in mode_order], rotation=45, ha='right', fontsize=7)
    ax.set_ylim(0, 105)
    legend_outside(fig, handles=mode_handles + metric_handles,
                   loc='upper left', bbox_to_anchor=(1.02, 1.0), ncol=1, fontsize=7)
    style_ax(ax, ylabel="Score (%)", faint_y_grid=True,
             title=f"{title_prefix}Paired End Precision / Recall / F1")
    fig.tight_layout()
    savefig(fig, output_dir / "paired_summary.png", dpi=300)


def _plot_paired_pr_scatter(df, output_path, mode_order, styler, baseline_mode, title_prefix=""):
    """Scatter: paired recall (geom mean) on X, paired dedup precision on Y."""
    fig, ax = plt.subplots(figsize=(3.5, 3.5))
    for _, row in df.iterrows():
        m = str(row.get('transcriptome_mode', 'unknown'))
        pp  = pd.to_numeric(row.get('paired_dedup_precision'), errors='coerce')
        r5  = pd.to_numeric(row.get('5prime_recall'), errors='coerce')
        r3  = pd.to_numeric(row.get('3prime_recall'), errors='coerce')
        if pd.isna(pp) or pd.isna(r5) or pd.isna(r3) or r5 < 0 or r3 < 0:
            continue
        pr = np.sqrt(r5 * r3)
        ax.scatter(pr * 100, pp * 100,
                   s=55, c=styler.color(m), marker=styler.marker(m),
                   edgecolors='white', linewidth=0.5, alpha=0.9, zorder=2)

    style_ax(ax,
             xlabel="Paired Recall (%) (\u221a(5\u2032\u00d73\u2032 recall))",
             ylabel="Paired Precision (%)",
             title=f"{title_prefix}Paired Precision vs Recall")
    ax.set_xlim(0, 105)
    ax.set_ylim(0, 105)
    ax.plot([0, 100], [0, 100], "--", color="#cccccc", linewidth=0.8, zorder=0)

    if baseline_mode:
        bl_sub = df[df['transcriptome_mode'] == baseline_mode]
        bl_pp = _safe_numeric(bl_sub, 'paired_dedup_precision').mean()
        bl_r5 = _safe_numeric(bl_sub, '5prime_recall').mean()
        bl_r3 = _safe_numeric(bl_sub, '3prime_recall').mean()
        if not (np.isnan(bl_pp) or np.isnan(bl_r5) or np.isnan(bl_r3)):
            bl_pr = np.sqrt(bl_r5 * bl_r3)
            ax.axhline(bl_pp * 100, color='#888888', linestyle=':', linewidth=1.0, zorder=1, alpha=0.7)
            ax.axvline(bl_pr * 100, color='#888888', linestyle=':', linewidth=1.0, zorder=1, alpha=0.7)

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

    # Legend: one patch per mode (matching scatter marker color) + hatch legend
    mode_handles = [
        mpatches.Patch(facecolor=styler.color(m), edgecolor='none',
                       label=_short_mode(m), alpha=0.9)
        for m in mode_order
    ]
    hatch_handles = [
        mpatches.Patch(facecolor='#666666', alpha=0.9, edgecolor='none', label="5\u2032 F1"),
        mpatches.Patch(facecolor='#666666', alpha=0.55, edgecolor='none', hatch='///', label="3\u2032 F1"),
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
# Per-method 5′ vs 3′ scatter grid
# ---------------------------------------------------------------------------

def _plot_5v3_scatter(
    ax,
    df: "pd.DataFrame",
    col_x: str,
    col_y: str,
    xlabel: str,
    ylabel: str,
    title: str,
    mode_order: list,
    styler: "ModeStyler",
    label_points: bool = True,
):
    """Scatter one metric per axis (5′ on X, 3′ on Y), one point per method."""
    for m in mode_order:
        row = df[df["transcriptome_mode"] == m]
        x = _safe_numeric(row, col_x).mean()
        y = _safe_numeric(row, col_y).mean()
        if pd.isna(x) or pd.isna(y):
            continue
        ax.scatter(x * 100, y * 100,
                   s=60, c=styler.color(m), marker=styler.marker(m),
                   edgecolors="white", linewidth=0.5, alpha=0.9, zorder=2)
        if label_points:
            ax.annotate(_short_mode(m), (x * 100, y * 100),
                        textcoords="offset points", xytext=(4, 2),
                        fontsize=5, color="#444444", zorder=3)
    ax.plot([0, 100], [0, 100], "--", color="#cccccc", linewidth=0.8, zorder=0)
    ax.set_xlim(0, 105)
    ax.set_ylim(0, 105)
    style_ax(ax, xlabel=xlabel, ylabel=ylabel, title=title)


def plot_per_method_5v3(
    df_ortho: "pd.DataFrame",
    df_ref: "pd.DataFrame | None",
    output_dir: "Path",
    mode_order: list,
    styler: "ModeStyler",
    title_prefix: str = "",
):
    """Save per-method 5′ vs 3′ scatter grid under output_dir/per_methods/.

    Each plot has one point per method.  Axes are the same metric measured
    at each end (5′ on X, 3′ on Y), so you can see whether methods trade
    off 5′ and 3′ performance symmetrically.

    Orthogonal plots (require df_ortho):
      5prime_recall × 3prime_recall         → ortho_recall_5v3.png
      5prime_precision × 3prime_precision   → ortho_precision_5v3.png
      5prime_f1 × 3prime_f1                 → ortho_f1_5v3.png

    Reference plots (require df_ref):
      same three metrics from GTF TSV       → ref_recall_5v3.png etc.

    Cross plots (require both):
      5prime ortho precision × ref precision → cross_5prime_ortho_vs_ref_prec.png
      3prime ortho precision × ref precision → cross_3prime_ortho_vs_ref_prec.png
      5prime ortho recall    × ref recall    → cross_5prime_ortho_vs_ref_recall.png
      3prime ortho recall    × ref recall    → cross_3prime_ortho_vs_ref_recall.png
    """
    per_methods_dir = Path(output_dir)
    per_methods_dir.mkdir(parents=True, exist_ok=True)
    # Filename convention (all flat in output_dir, captured by *.png glob):
    #   5v3_*        : same metric at 5′ (X) vs 3′ (Y) — one point per method
    #   paired_*     : paired-end precision/recall derived metrics
    #   cross_*      : per-end orthogonal vs reference comparison

    def _f1(df, end):
        """Compute per-mode F1 and return as a synthetic single-column df."""
        rows = []
        for m in mode_order:
            sub = df[df["transcriptome_mode"] == m]
            p = _safe_numeric(sub, f"{end}_precision").mean()
            r = _safe_numeric(sub, f"{end}_recall").mean()
            f1 = 2 * p * r / (p + r) if not (pd.isna(p) or pd.isna(r)) and (p + r) > 0 else np.nan
            rows.append({"transcriptome_mode": m, f"{end}_f1": f1})
        return pd.DataFrame(rows)

    # ── Orthogonal plots ────────────────────────────────────────────────────
    if df_ortho is not None:
        df_f1_5o = _f1(df_ortho, "5prime")
        df_f1_3o = _f1(df_ortho, "3prime")
        df_ortho_f1 = df_f1_5o.merge(df_f1_3o, on="transcriptome_mode")

        specs_ortho = [
            ("5prime_recall",    "3prime_recall",    "5′ Recall (%)",    "3′ Recall (%)",    "Orthogonal recall: 5′ vs 3′",    "5v3_ortho_recall"),
            ("5prime_precision", "3prime_precision", "5′ Precision (%)", "3′ Precision (%)", "Orthogonal precision: 5′ vs 3′", "5v3_ortho_prec"),
            ("5prime_f1",        "3prime_f1",        "5′ F1 (%)",        "3′ F1 (%)",        "Orthogonal F1: 5′ vs 3′",        "5v3_ortho_f1"),
        ]
        for col_x, col_y, xlabel, ylabel, title, fname in specs_ortho:
            fig, ax = plt.subplots(figsize=(3.5, 3.5))
            src = df_ortho_f1 if "f1" in col_x else df_ortho
            _plot_5v3_scatter(ax, src, col_x, col_y, xlabel, ylabel,
                              f"{title_prefix}{title}", mode_order, styler)
            handles = [styler.legend_handle(m, label=_short_mode(m), markersize=7)
                       for m in mode_order]
            legend_outside(fig, handles=handles, loc="upper left",
                           bbox_to_anchor=(1.02, 1.0), ncol=1, fontsize=7)
            fig.tight_layout()
            savefig(fig, per_methods_dir / fname, dpi=300)

    # ── Reference plots ─────────────────────────────────────────────────────
    if df_ref is not None:
        df_f1_5r = _f1(df_ref, "5prime")
        df_f1_3r = _f1(df_ref, "3prime")
        df_ref_f1 = df_f1_5r.merge(df_f1_3r, on="transcriptome_mode")

        specs_ref = [
            ("5prime_recall",    "3prime_recall",    "5′ Recall (%)",    "3′ Recall (%)",    "Reference recall: 5′ vs 3′",    "5v3_ref_recall"),
            ("5prime_precision", "3prime_precision", "5′ Precision (%)", "3′ Precision (%)", "Reference precision: 5′ vs 3′", "5v3_ref_prec"),
            ("5prime_f1",        "3prime_f1",        "5′ F1 (%)",        "3′ F1 (%)",        "Reference F1: 5′ vs 3′",        "5v3_ref_f1"),
        ]
        for col_x, col_y, xlabel, ylabel, title, fname in specs_ref:
            fig, ax = plt.subplots(figsize=(3.5, 3.5))
            src = df_ref_f1 if "f1" in col_x else df_ref
            _plot_5v3_scatter(ax, src, col_x, col_y, xlabel, ylabel,
                              f"{title_prefix}{title}", mode_order, styler)
            handles = [styler.legend_handle(m, label=_short_mode(m), markersize=7)
                       for m in mode_order]
            legend_outside(fig, handles=handles, loc="upper left",
                           bbox_to_anchor=(1.02, 1.0), ncol=1, fontsize=7)
            fig.tight_layout()
            savefig(fig, per_methods_dir / fname, dpi=300)

    # ── Paired plots ────────────────────────────────────────────────────────
    # Paired precision = paired_dedup_precision
    # Paired recall    = sqrt(5prime_recall * 3prime_recall)
    # Paired F1        = harmonic mean of the two
    def _paired_series(df):
        """Return df with columns paired_prec, paired_rec, paired_f1 per mode."""
        rows = []
        for m in mode_order:
            sub = df[df["transcriptome_mode"] == m]
            pp = _safe_numeric(sub, "paired_dedup_precision").mean()
            r5 = _safe_numeric(sub, "5prime_recall").mean()
            r3 = _safe_numeric(sub, "3prime_recall").mean()
            pr = np.sqrt(r5 * r3) if not (pd.isna(r5) or pd.isna(r3)) and r5 >= 0 and r3 >= 0 else np.nan
            f1 = 2 * pp * pr / (pp + pr) if not (pd.isna(pp) or pd.isna(pr)) and (pp + pr) > 0 else np.nan
            rows.append({"transcriptome_mode": m, "paired_prec": pp,
                         "paired_rec": pr, "paired_f1": f1})
        return pd.DataFrame(rows)

    if df_ortho is not None:
        ps_o = _paired_series(df_ortho)
        fig, ax = plt.subplots(figsize=(3.5, 3.5))
        for m in mode_order:
            row = ps_o[ps_o["transcriptome_mode"] == m]
            x = row["paired_rec"].values[0] if len(row) else np.nan
            y = row["paired_prec"].values[0] if len(row) else np.nan
            if pd.isna(x) or pd.isna(y):
                continue
            ax.scatter(x * 100, y * 100,
                       s=60, c=styler.color(m), marker=styler.marker(m),
                       edgecolors="white", linewidth=0.5, alpha=0.9, zorder=2)
            ax.annotate(_short_mode(m), (x * 100, y * 100),
                        textcoords="offset points", xytext=(4, 2),
                        fontsize=5, color="#444444", zorder=3)
        ax.plot([0, 100], [0, 100], "--", color="#cccccc", linewidth=0.8, zorder=0)
        ax.set_xlim(0, 105); ax.set_ylim(0, 105)
        style_ax(ax, xlabel="Paired Recall (%) (√r5·r3)",
                 ylabel="Paired Precision (%)",
                 title=f"{title_prefix}Orthogonal paired precision vs recall")
        handles = [styler.legend_handle(m, label=_short_mode(m), markersize=7) for m in mode_order]
        legend_outside(fig, handles=handles, loc="upper left",
                       bbox_to_anchor=(1.02, 1.0), ncol=1, fontsize=7)
        fig.tight_layout()
        savefig(fig, per_methods_dir / "paired_ortho_pr", dpi=300)

    if df_ref is not None:
        ps_r = _paired_series(df_ref)
        fig, ax = plt.subplots(figsize=(3.5, 3.5))
        for m in mode_order:
            row = ps_r[ps_r["transcriptome_mode"] == m]
            x = row["paired_rec"].values[0] if len(row) else np.nan
            y = row["paired_prec"].values[0] if len(row) else np.nan
            if pd.isna(x) or pd.isna(y):
                continue
            ax.scatter(x * 100, y * 100,
                       s=60, c=styler.color(m), marker=styler.marker(m),
                       edgecolors="white", linewidth=0.5, alpha=0.9, zorder=2)
            ax.annotate(_short_mode(m), (x * 100, y * 100),
                        textcoords="offset points", xytext=(4, 2),
                        fontsize=5, color="#444444", zorder=3)
        ax.plot([0, 100], [0, 100], "--", color="#cccccc", linewidth=0.8, zorder=0)
        ax.set_xlim(0, 105); ax.set_ylim(0, 105)
        style_ax(ax, xlabel="Paired Recall (%) (√r5·r3)",
                 ylabel="Paired Precision (%)",
                 title=f"{title_prefix}Reference paired precision vs recall")
        handles = [styler.legend_handle(m, label=_short_mode(m), markersize=7) for m in mode_order]
        legend_outside(fig, handles=handles, loc="upper left",
                       bbox_to_anchor=(1.02, 1.0), ncol=1, fontsize=7)
        fig.tight_layout()
        savefig(fig, per_methods_dir / "paired_ref_pr", dpi=300)

    if df_ortho is not None and df_ref is not None:
        ps_o = _paired_series(df_ortho)
        ps_r = _paired_series(df_ref)
        for metric, xlabel, ylabel, fname in [
            ("paired_prec", "Reference Paired Precision (%)", "Orthogonal Paired Precision (%)", "paired_cross_prec"),
            ("paired_rec",  "Reference Paired Recall (%)",    "Orthogonal Paired Recall (%)",    "paired_cross_recall"),
            ("paired_f1",   "Reference Paired F1 (%)",        "Orthogonal Paired F1 (%)",        "paired_cross_f1"),
        ]:
            fig, ax = plt.subplots(figsize=(3.5, 3.5))
            handles = []
            for m in mode_order:
                ro = ps_o[ps_o["transcriptome_mode"] == m]
                rr = ps_r[ps_r["transcriptome_mode"] == m]
                if ro.empty or rr.empty:
                    continue
                x = rr[metric].values[0]
                y = ro[metric].values[0]
                if pd.isna(x) or pd.isna(y):
                    continue
                ax.scatter(x * 100, y * 100,
                           s=60, c=styler.color(m), marker=styler.marker(m),
                           edgecolors="white", linewidth=0.5, alpha=0.9, zorder=2)
                ax.annotate(_short_mode(m), (x * 100, y * 100),
                            textcoords="offset points", xytext=(4, 2),
                            fontsize=5, color="#444444", zorder=3)
                handles.append(styler.legend_handle(m, label=_short_mode(m), markersize=7))
            ax.plot([0, 100], [0, 100], "--", color="#cccccc", linewidth=0.8, zorder=0)
            ax.set_xlim(0, 105); ax.set_ylim(0, 105)
            label = metric.replace("paired_", "").capitalize()
            style_ax(ax, xlabel=xlabel, ylabel=ylabel,
                     title=f"{title_prefix}Paired {label}: Ortho vs Ref")
            legend_outside(fig, handles=handles, loc="upper left",
                           bbox_to_anchor=(1.02, 1.0), ncol=1, fontsize=7)
            fig.tight_layout()
            savefig(fig, per_methods_dir / fname, dpi=300)

    # ── Cross plots (orthogonal metric on Y, reference metric on X) ─────────
    if df_ortho is not None and df_ref is not None:
        cross_specs = [
            ("5prime", "precision", "5′ Reference Precision (%)", "5′ Orthogonal Precision (%)"),
            ("3prime", "precision", "3′ Reference Precision (%)", "3′ Orthogonal Precision (%)"),
            ("5prime", "recall",    "5′ Reference Recall (%)",    "5′ Orthogonal Recall (%)"),
            ("3prime", "recall",    "3′ Reference Recall (%)",    "3′ Orthogonal Recall (%)"),
        ]
        for end, metric, xlabel, ylabel in cross_specs:
            col = f"{end}_{metric}"
            fig, ax = plt.subplots(figsize=(3.5, 3.5))
            handles = []
            for m in mode_order:
                row_o = df_ortho[df_ortho["transcriptome_mode"] == m]
                row_r = df_ref[df_ref["transcriptome_mode"] == m]
                if row_o.empty or row_r.empty:
                    continue
                x = _safe_numeric(row_r, col).mean()
                y = _safe_numeric(row_o, col).mean()
                if pd.isna(x) or pd.isna(y):
                    continue
                ax.scatter(x * 100, y * 100,
                           s=60, c=styler.color(m), marker=styler.marker(m),
                           edgecolors="white", linewidth=0.5, alpha=0.9, zorder=2)
                ax.annotate(_short_mode(m), (x * 100, y * 100),
                            textcoords="offset points", xytext=(4, 2),
                            fontsize=5, color="#444444", zorder=3)
                handles.append(styler.legend_handle(m, label=_short_mode(m), markersize=7))
            ax.plot([0, 100], [0, 100], "--", color="#cccccc", linewidth=0.8, zorder=0)
            ax.set_xlim(0, 105)
            ax.set_ylim(0, 105)
            end_short = "5′" if end == "5prime" else "3′"
            style_ax(ax, xlabel=xlabel, ylabel=ylabel,
                     title=f"{title_prefix}{end_short} Ortho vs Ref {metric.capitalize()}")
            legend_outside(fig, handles=handles, loc="upper left",
                           bbox_to_anchor=(1.02, 1.0), ncol=1, fontsize=7)
            fig.tight_layout()
            end_tag = "5p" if end == "5prime" else "3p"
            savefig(fig, per_methods_dir / f"cross_{end_tag}_{metric}", dpi=300)


# ---------------------------------------------------------------------------
# Reference vs Orthogonal precision scatter
# ---------------------------------------------------------------------------

def plot_ref_vs_ortho_precision(
    df_ortho: "pd.DataFrame",
    df_ref: "pd.DataFrame",
    output_path: "Path",
    mode_order: list,
    styler: "ModeStyler",
    title_prefix: str = "",
):
    """Scatter: GTF-reference precision (X) vs orthogonal-signal precision (Y).

    Each point is one method.  Points above the diagonal have better orthogonal
    precision than annotation concordance (they place ends on real signal even
    when that differs from the annotation).  Points below rely more on annotation.
    Separate panels for 5′ and 3′.
    """
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.5))

    for ax, end, end_label in [
        (axes[0], "5prime", "5′ TSS"),
        (axes[1], "3prime", "3′ TTS"),
    ]:
        pcol = f"{end}_precision"
        handles = []
        for m in mode_order:
            row_o = df_ortho[df_ortho["transcriptome_mode"] == m]
            row_r = df_ref[df_ref["transcriptome_mode"] == m]
            if row_o.empty or row_r.empty:
                continue
            p_ortho = _safe_numeric(row_o, pcol).mean()
            p_ref   = _safe_numeric(row_r, pcol).mean()
            if pd.isna(p_ortho) or pd.isna(p_ref):
                continue
            ax.scatter(
                p_ref * 100, p_ortho * 100,
                s=60, c=styler.color(m), marker=styler.marker(m),
                edgecolors="white", linewidth=0.5, alpha=0.9, zorder=2,
            )
            handles.append(styler.legend_handle(m, label=m.replace("_", " "), markersize=7))

        # Diagonal: orthogonal == reference
        ax.plot([0, 100], [0, 100], "--", color="#cccccc", linewidth=0.8, zorder=0)
        ax.set_xlim(0, 105)
        ax.set_ylim(0, 105)
        style_ax(ax,
                 xlabel=f"{end_label} Reference precision (%)",
                 ylabel=f"{end_label} Orthogonal precision (%)",
                 title=f"{title_prefix}{end_label}")

        # Shade quadrants
        ax.axhspan(50, 105, xmin=0, xmax=0.5, alpha=0.03, color="steelblue")
        ax.axhspan(0,  50,  xmin=0.5, xmax=1, alpha=0.03, color="salmon")

        ax.text(2, 98, "better\northogonal", fontsize=5, color="#888888",
                va="top", style="italic")
        ax.text(52, 2, "better\nannotation", fontsize=5, color="#888888",
                va="bottom", style="italic")

    legend_outside(fig, handles=handles, loc="upper left",
                   bbox_to_anchor=(1.02, 1.0), ncol=1, fontsize=7)
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
    _plot_paired_pr_scatter(df, output_dir / "pr_paired_scatter.png", mode_order, styler,
                            baseline_mode, title_prefix)
    _plot_paired_summary(df, output_dir, mode_order, styler, baseline_mode, title_prefix)
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Create individual precision-recall, isoforms/gene, and F1 plots"
    )
    parser.add_argument('--input', '-i', nargs='+', required=True,
                        help="Input evaluation TSV file(s) — orthogonal-signal precision")
    parser.add_argument('--gtf-input', nargs='+', default=None,
                        help="GTF-reference precision TSV file(s) — enables ref-vs-ortho scatter")
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

    # Optional: reference-vs-orthogonal precision scatter + per-method 5′×3′ grid
    if args.gtf_input:
        df_ref = load_evaluation_files(args.gtf_input)
        if df_ref is not None and len(df_ref) > 0:
            mode_order = _resolve_mode_order(df)
            styler = ModeStyler(mode_order)
            output_dir = Path(args.output)
            output_dir.mkdir(parents=True, exist_ok=True)
            plot_ref_vs_ortho_precision(
                df, df_ref,
                output_dir / "pr_ref_vs_ortho_precision.png",
                mode_order, styler, args.title_prefix,
            )
            plot_per_method_5v3(
                df, df_ref, output_dir, mode_order, styler, args.title_prefix,
            )
            if args.verbose:
                print("Saved reference-vs-orthogonal precision scatter and per-method 5v3 grid")
    else:
        # Orthogonal-only per-method plots (no cross plots)
        mode_order = _resolve_mode_order(df)
        styler = ModeStyler(mode_order)
        output_dir = Path(args.output)
        plot_per_method_5v3(
            df, None, output_dir, mode_order, styler, args.title_prefix,
        )
        if args.verbose:
            print("Saved per-method 5v3 orthogonal plots")

    if success:
        print(f"Saved precision-recall plots to {args.output}")
    else:
        print("Failed to create plots", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
