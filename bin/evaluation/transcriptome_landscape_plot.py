#!/usr/bin/env python3
"""
Transcriptome landscape dashboard: cross-mode comparison of key assembly metrics.

Reads merged evaluation TSVs and produces a 6-panel summary figure comparing
transcriptome properties across assembler modes for a single sample:

  Panel 1: % Reads Assigned
  Panel 2: Isoforms per Gene
  Panel 3: SQANTI Category Breakdown (FSM, ISM, NIC, NNC, SEM, SEN)
  Panel 4: 5'/3' Precision (All vs Reference-only)
  Panel 5: 5'/3' Recall (All vs Reference-only)
  Panel 6: Read Support Quality (mean reads/isoform, single-read %, well-supported)

Usage:
    python transcriptome_landscape_plot.py \\
        --input eval1.tsv eval2.tsv ... \\
        --output transcriptome_landscape.png \\
        [--title-prefix "test_name: "] [--verbose]
"""

import argparse
import sys
from pathlib import Path

import numpy as np

try:
    import pandas as pd
except ImportError:
    pd = None

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from pub_style import style_ax, savefig, MODE_COLORS as _PUB_MODE_COLORS, PALETTE


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MODE_COLORS = _PUB_MODE_COLORS

DEFAULT_MODE_ORDER = [
    "default",
    "density-asymmetric",
    "density-asymmetric-softclip",
    "density-plain",
    "density-strict",
    "k-means",
    "more-ends",
    "bambu_default",
    "isoquant_pacbio",
]

# SQANTI category colours (distinct, colorblind-friendly)
SQANTI_COLORS = {
    "FSM": "#009E73",   # bluish green
    "ISM": "#56B4E9",   # sky blue
    "NIC": "#E69F00",   # orange
    "NNC": "#D55E00",   # vermillion
    "SEM": "#CC79A7",   # reddish purple
    "SEN": "#999999",   # grey
}
SQANTI_ORDER = ["FSM", "ISM", "NIC", "NNC", "SEM", "SEN"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _identify_assembler(mode: str) -> str:
    if "bambu" in mode:
        return "bambu"
    if "isoquant" in mode:
        return "isoquant"
    return "flair"


def _short_mode(name: str) -> str:
    """Abbreviate mode names for axis labels."""
    return (
        name.replace("density-", "d-")
        .replace("_default", "")
        .replace("_pacbio", "")
        .replace("end-scoring-", "es-")
    )


def _resolve_mode_order(df):
    """Return ordered list of modes present in the DataFrame."""
    present = set(df["transcriptome_mode"].unique())
    order = [m for m in DEFAULT_MODE_ORDER if m in present]
    for m in sorted(present):
        if m not in order:
            order.append(m)
    return order


def _mode_color(mode: str) -> str:
    """Assign a deterministic colour to each mode, cycling PALETTE for unknowns."""
    if mode in MODE_COLORS:
        return MODE_COLORS[mode]
    # Deterministic fallback from PALETTE based on hash
    idx = hash(mode) % len(PALETTE)
    return PALETTE[idx]


def load_evaluation_files(input_files):
    """Load and concatenate evaluation TSVs into a DataFrame."""
    if pd is None:
        print("pandas required for transcriptome_landscape_plot.py", file=sys.stderr)
        return None
    dfs = []
    for f in input_files:
        try:
            df = pd.read_csv(f, sep="\t")
            if not df.empty:
                dfs.append(df)
        except Exception as e:
            print(f"Warning: could not read {f}: {e}", file=sys.stderr)
    if not dfs:
        return None
    return pd.concat(dfs, ignore_index=True)


def _safe_numeric(df, col):
    """Convert column to numeric, coercing errors to NaN."""
    if col not in df.columns:
        return pd.Series([np.nan] * len(df), index=df.index)
    return pd.to_numeric(df[col], errors="coerce")


# ---------------------------------------------------------------------------
# Individual panel plotters
# ---------------------------------------------------------------------------

def _plot_reads_assigned(ax, df, mode_order):
    """Panel 1: % reads assigned as bar chart."""
    x = np.arange(len(mode_order))
    rates = []
    colors = []
    for mode in mode_order:
        subset = df[df["transcriptome_mode"] == mode]
        val = _safe_numeric(subset, "assignment_rate").mean()
        rates.append(val * 100 if not np.isnan(val) else 0)
        colors.append(_mode_color(mode))

    bars = ax.bar(x, rates, color=colors, alpha=0.85, edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels([_short_mode(m) for m in mode_order], rotation=45, ha="right", fontsize=8)
    ax.set_ylim(0, 105)

    # Annotate bars
    for bar, rate in zip(bars, rates):
        if rate > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                    f"{rate:.1f}%", ha="center", va="bottom", fontsize=7, color="#333333")

    style_ax(ax, ylabel="% Reads Assigned", faint_y_grid=True)


def _plot_isoforms_per_gene(ax, df, mode_order):
    """Panel 2: bar — isoforms per gene (combined)."""
    x = np.arange(len(mode_order))
    width = 0.6

    tpg_all = []
    tpg_na = []   # track which bars have no gene count (e.g. GTF-based tools)
    for mode in mode_order:
        subset = df[df["transcriptome_mode"] == mode]
        total_isos = _safe_numeric(subset, "isoforms_observed").sum()
        total_genes = _safe_numeric(subset, "genes_observed").sum()

        if total_genes > 0 and not np.isnan(total_genes):
            tpg_all.append(total_isos / total_genes)
            tpg_na.append(False)
        else:
            tpg_all.append(0)
            tpg_na.append(True)

    colors = [_mode_color(m) for m in mode_order]
    ax.bar(x, tpg_all, width=width * 0.9, color=colors,
           alpha=0.85, edgecolor="none")

    y_max = max(tpg_all) if tpg_all else 1
    for i, (val, na) in enumerate(zip(tpg_all, tpg_na)):
        if na:
            ax.text(i, y_max * 0.02 + 0.05, "N/A",
                    ha="center", va="bottom", fontsize=6.5, color="#888888")
        elif val > 0:
            ax.text(i, val + 0.03, f"{val:.1f}",
                    ha="center", va="bottom", fontsize=6.5, color="#333333")

    ax.set_xticks(x)
    ax.set_xticklabels([_short_mode(m) for m in mode_order], rotation=45, ha="right", fontsize=8)
    style_ax(ax, ylabel="Isoforms / Gene", faint_y_grid=True)


def _plot_sqanti_breakdown(ax, df, mode_order):
    """Panel 3: stacked bar — SQANTI category percentages."""
    x = np.arange(len(mode_order))
    category_arrs = {cat: [] for cat in SQANTI_ORDER}

    for mode in mode_order:
        subset = df[df["transcriptome_mode"] == mode]
        counts = {}
        for cat in SQANTI_ORDER:
            val = _safe_numeric(subset, cat).sum()
            counts[cat] = val if not np.isnan(val) else 0
        total = sum(counts.values())
        for cat in SQANTI_ORDER:
            pct = (counts[cat] / total * 100) if total > 0 else 0
            category_arrs[cat].append(pct)

    bottom = np.zeros(len(mode_order))
    for cat in SQANTI_ORDER:
        vals = np.array(category_arrs[cat])
        ax.bar(x, vals, bottom=bottom, color=SQANTI_COLORS[cat], alpha=0.85,
               edgecolor="none", label=cat)
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels([_short_mode(m) for m in mode_order], rotation=45, ha="right", fontsize=8)
    ax.set_ylim(0, 105)
    ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1.01, 1.0),
              frameon=False, ncol=1, borderaxespad=0)
    style_ax(ax, ylabel="% of Isoforms", faint_y_grid=True)


def _plot_precision(ax, df, mode_order):
    """Panel 4: grouped bar — 5'/3' precision."""
    x = np.arange(len(mode_order))
    width = 0.3
    metrics = [
        ("5prime_precision", "5\u2032"),
        ("3prime_precision", "3\u2032"),
    ]
    colors = [PALETTE[4], PALETTE[5]]

    for i, (col, label) in enumerate(metrics):
        vals = []
        for mode in mode_order:
            subset = df[df["transcriptome_mode"] == mode]
            v = _safe_numeric(subset, col).mean()
            vals.append(v * 100 if not np.isnan(v) else 0)
        offset = (i - 0.5) * width
        ax.bar(x + offset, vals, width=width * 0.9, color=colors[i], alpha=0.85,
               edgecolor="none", label=label)

    ax.set_xticks(x)
    ax.set_xticklabels([_short_mode(m) for m in mode_order], rotation=45, ha="right", fontsize=8)
    ax.set_ylim(0, 105)
    ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1.01, 1.0),
              frameon=False, ncol=1, borderaxespad=0)
    style_ax(ax, ylabel="Precision (%)", faint_y_grid=True)


def _plot_recall(ax, df, mode_order):
    """Panel 5: grouped bar — 5'/3' recall."""
    x = np.arange(len(mode_order))
    width = 0.3
    metrics = [
        ("5prime_recall", "5\u2032"),
        ("3prime_recall", "3\u2032"),
    ]
    colors = [PALETTE[4], PALETTE[5]]

    for i, (col, label) in enumerate(metrics):
        vals = []
        for mode in mode_order:
            subset = df[df["transcriptome_mode"] == mode]
            v = _safe_numeric(subset, col).mean()
            vals.append(v * 100 if not np.isnan(v) else 0)
        offset = (i - 0.5) * width
        ax.bar(x + offset, vals, width=width * 0.9, color=colors[i], alpha=0.85,
               edgecolor="none", label=label)

    ax.set_xticks(x)
    ax.set_xticklabels([_short_mode(m) for m in mode_order], rotation=45, ha="right", fontsize=8)
    ax.set_ylim(0, 105)
    ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1.01, 1.0),
              frameon=False, ncol=1, borderaxespad=0)
    style_ax(ax, ylabel="Recall (%)", faint_y_grid=True)


def _plot_read_support(ax, df, mode_order):
    """Panel 6: bar — mean reads per isoform with single-read % and well-supported % overlaid.

    Falls back to plotting mean reads/isoform only when the per-isoform breakdown
    columns (single_read_isoform_rate, well_supported_isoforms) are absent from the
    evaluation TSV.
    """
    x = np.arange(len(mode_order))
    width = 0.35

    mean_rpi = []
    single_pct = []
    well_pct = []
    has_breakdown = (
        "single_read_isoform_rate" in df.columns
        and df["single_read_isoform_rate"].notna().any()
    )

    for mode in mode_order:
        subset = df[df["transcriptome_mode"] == mode]
        rpi = _safe_numeric(subset, "reads_per_isoform_mean").mean()
        mean_rpi.append(rpi if not np.isnan(rpi) else 0)
        if has_breakdown:
            sr = _safe_numeric(subset, "single_read_isoform_rate").mean()
            ws = _safe_numeric(subset, "well_supported_isoforms").sum()
            total = _safe_numeric(subset, "isoforms_observed").sum()
            single_pct.append(sr * 100 if not np.isnan(sr) else 0)
            well_pct.append((ws / total * 100) if total > 0 and not np.isnan(ws) else 0)

    colors = [_mode_color(m) for m in mode_order]

    if has_breakdown:
        bars1 = ax.bar(x - width / 2, single_pct, width=width * 0.9, color=PALETTE[5],
                       alpha=0.85, edgecolor="none", label="Single-read %")
        bars2 = ax.bar(x + width / 2, well_pct, width=width * 0.9, color=PALETTE[2],
                       alpha=0.85, edgecolor="none", label="Well-supported %")
        y_max = max(max(single_pct), max(well_pct)) if (single_pct or well_pct) else 1
        ax.set_ylim(0, y_max * 1.15 + 1)
        for i, rpi in enumerate(mean_rpi):
            y_top = max(single_pct[i], well_pct[i])
            ax.text(i, y_top + y_top * 0.03 + 1, f"\u03bc={rpi:.0f}",
                    ha="center", va="bottom", fontsize=6.5, color="#555555", style="italic")
        ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1.01, 1.0),
                  frameon=False, ncol=1, borderaxespad=0)
        style_ax(ax, ylabel="% of Isoforms", faint_y_grid=True)
    else:
        # Fallback: plot mean reads/isoform as a bar chart
        bars = ax.bar(x, mean_rpi, color=colors, alpha=0.85, edgecolor="none")
        for bar, val in zip(bars, mean_rpi):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + bar.get_height() * 0.02 + 0.5,
                        f"\u03bc={val:.0f}", ha="center", va="bottom", fontsize=6.5,
                        color="#333333")
        y_max = max(mean_rpi) if mean_rpi else 1
        ax.set_ylim(0, y_max * 1.2 + 1)
        style_ax(ax, ylabel="Mean Reads / Isoform", faint_y_grid=True)

    ax.set_xticks(x)
    ax.set_xticklabels([_short_mode(m) for m in mode_order], rotation=45, ha="right", fontsize=8)


def _plot_gene_detection(ax, df, mode_order):
    """Panel 7: bar — genes observed per mode.

    gene_detection_rate (genes observed / reference genes) is not written to the
    evaluation TSV, so we plot the raw genes_observed count instead.  Tools that
    output GTF without ENSEMBL-prefixed gene IDs (e.g. IsoQuant) will show 0 — this
    is a data limitation noted in the bar annotation rather than treated as an error.
    """
    x = np.arange(len(mode_order))
    observed = []
    colors = []
    for mode in mode_order:
        subset = df[df["transcriptome_mode"] == mode]
        obs = _safe_numeric(subset, "genes_observed").sum()
        observed.append(int(obs) if not np.isnan(obs) else 0)
        colors.append(_mode_color(mode))

    bars = ax.bar(x, observed, color=colors, alpha=0.85, edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels([_short_mode(m) for m in mode_order], rotation=45, ha="right", fontsize=8)
    y_max = max(observed) if observed else 1
    ax.set_ylim(0, y_max * 1.2 + 1)

    for bar, obs in zip(bars, observed):
        label = str(obs) if obs > 0 else "N/A"
        ax.text(bar.get_x() + bar.get_width() / 2,
                max(bar.get_height(), 0) + y_max * 0.02 + 0.5,
                label, ha="center", va="bottom", fontsize=6.5, color="#333333")

    style_ax(ax, ylabel="Genes Observed", faint_y_grid=True)


def _plot_end_diversity_recovery(ax, df, mode_order):
    """Panel 8: Cleveland-style connected dot plot — end diversity & peak recovery.

    Two connected-dot sub-panels side by side:
      Left:  TSS vs TTS end-diversity calibration
      Right: 5' vs 3' peak recovery rate

    Each method gets a horizontal line connecting its 5'/TSS and 3'/TTS dots.
    """
    pass  # replaced by dedicated functions below


def _plot_cleveland_dot(
    df, mode_order, metric_a, metric_b, label_a, label_b,
    xlabel, output_path, show_ref_line=None,
):
    """Cleveland-style connected dot plot comparing two paired metrics across methods.

    Parameters
    ----------
    metric_a, metric_b : column names in *df*
    label_a, label_b   : legend text for the two dot series
    xlabel              : x-axis label
    show_ref_line       : optional x-value for a vertical reference line (e.g. 100%)
    """
    from pub_style import W1

    methods = list(reversed(mode_order))          # top-first ordering
    y_pos   = np.arange(len(methods))

    vals_a, vals_b = [], []
    for mode in methods:
        subset = df[df["transcriptome_mode"] == mode]
        va = _safe_numeric(subset, metric_a).mean()
        vb = _safe_numeric(subset, metric_b).mean()
        vals_a.append(va * 100 if not np.isnan(va) else np.nan)
        vals_b.append(vb * 100 if not np.isnan(vb) else np.nan)

    fig_h = max(2.0, 0.38 * len(methods) + 0.8)
    fig, ax = plt.subplots(figsize=(W1, fig_h))

    # Reference line
    if show_ref_line is not None:
        ax.axvline(show_ref_line, color="#cccccc", linewidth=0.8, zorder=0)

    # Connecting segments
    for i, (a, b) in enumerate(zip(vals_a, vals_b)):
        if np.isnan(a) or np.isnan(b):
            continue
        ax.plot([a, b], [i, i], color="#bbbbbb", linewidth=0.8, zorder=1)

    # Dots
    ax.scatter(vals_a, y_pos, color=PALETTE[4], s=30, zorder=2, label=label_a,
               edgecolors="none")
    ax.scatter(vals_b, y_pos, color=PALETTE[5], s=30, zorder=2, label=label_b,
               edgecolors="none", marker="D")

    ax.set_yticks(y_pos)
    ax.set_yticklabels([_short_mode(m) for m in methods], fontsize=7)
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0),
              frameon=False, ncol=1, borderaxespad=0)
    style_ax(ax, xlabel=xlabel, faint_y_grid=True)
    fig.tight_layout()
    savefig(fig, Path(output_path), dpi=300)


# ---------------------------------------------------------------------------
# Individual figure wrappers
# ---------------------------------------------------------------------------

def _single_panel(plot_fn, df, mode_order, output_path, figsize=(3.5, 3.0)):
    """Create a standalone figure calling a panel plotter that takes (ax, df, mode_order)."""
    fig, ax = plt.subplots(figsize=figsize)
    plot_fn(ax, df, mode_order)
    fig.tight_layout()
    savefig(fig, Path(output_path), dpi=300)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def create_landscape_plots(df, output_dir, title_prefix=""):
    """Create individual transcriptome landscape figures in output_dir."""
    mode_order = _resolve_mode_order(df)
    n_modes = len(mode_order)

    if n_modes == 0:
        print("No modes found in data — skipping landscape plots.", file=sys.stderr)
        return False

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fw = max(3.5, n_modes * 0.7)

    _single_panel(_plot_reads_assigned, df, mode_order,
                  output_dir / "reads_assigned.png", figsize=(fw, 3.0))
    _single_panel(_plot_isoforms_per_gene, df, mode_order,
                  output_dir / "isoforms_per_gene.png", figsize=(fw, 3.0))
    _single_panel(_plot_sqanti_breakdown, df, mode_order,
                  output_dir / "sqanti_breakdown.png", figsize=(fw, 3.0))
    _single_panel(_plot_precision, df, mode_order,
                  output_dir / "end_precision.png", figsize=(fw, 3.0))
    _single_panel(_plot_recall, df, mode_order,
                  output_dir / "end_recall.png", figsize=(fw, 3.0))
    _single_panel(_plot_read_support, df, mode_order,
                  output_dir / "read_support.png", figsize=(fw, 3.0))
    _single_panel(_plot_gene_detection, df, mode_order,
                  output_dir / "gene_detection.png", figsize=(fw, 3.0))

    # Alternative end recovery — 5'/3' recall as a proxy for end calibration
    # (tss_end_diversity_calibration / tts_end_diversity_calibration are not written
    # to the evaluation TSV; 5prime_recall / 3prime_recall are the best available
    # per-end accuracy metrics)
    _plot_cleveland_dot(
        df, mode_order,
        "5prime_recall", "3prime_recall",
        "5\u2032 Recall", "3\u2032 Recall",
        xlabel="End Recall (%)",
        output_path=output_dir / "alternative_end_recovery.png",
        show_ref_line=100,
    )
    _plot_cleveland_dot(
        df, mode_order,
        "5prime_peak_recovery_rate", "3prime_peak_recovery_rate",
        "5\u2032", "3\u2032",
        xlabel="Peak Recovery Rate (%)",
        output_path=output_dir / "peak_recovery_rate.png",
    )

    print(f"Saved transcriptome landscape plots to {output_dir}")
    return True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Transcriptome landscape: individual cross-mode metric plots",
    )
    parser.add_argument("--input", nargs="+", required=True,
                        help="Evaluation TSV files")
    parser.add_argument("--output", required=True, help="Output directory for individual plots")
    parser.add_argument("--title-prefix", default="", help="Title prefix")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if pd is None:
        print("ERROR: pandas is required", file=sys.stderr)
        sys.exit(1)

    df = load_evaluation_files(args.input)
    if df is None or df.empty:
        print("No data loaded — skipping plot", file=sys.stderr)
        sys.exit(0)

    if "transcriptome_mode" not in df.columns:
        print("No 'transcriptome_mode' column — skipping plot", file=sys.stderr)
        sys.exit(0)

    if args.verbose:
        modes = df["transcriptome_mode"].unique()
        print(f"Loaded {len(df)} rows, {len(modes)} modes: {list(modes)}")

    success = create_landscape_plots(df, args.output, args.title_prefix)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
