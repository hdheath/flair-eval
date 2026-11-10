#!/usr/bin/env python3
"""
Cross-run summary plots for concordance and boundary selection metrics.

Reads merged evaluation TSVs from the test suite and produces a 2x3 panel figure:
  Row 1: 5' (TSS / CAGE)
  Row 2: 3' (TTS / dRNA)
  Col 1: Spearman rho  -- read-support x expression concordance  (Metric 1)
  Col 2: Concordance index -- boundary selection quality          (Metric 2a)
  Col 3: Signal vs end count -- over-segmentation diagnostic      (Metric 2b)

Usage:
    python concordance_plots.py \\
        --input eval1.tsv eval2.tsv ... \\
        --output concordance_summary.png \\
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

ASSEMBLER_COLORS = {
    "flair": PALETTE[4],
    "bambu": PALETTE[0],
    "isoquant": PALETTE[2],
    "unknown": PALETTE[7],
}

# Fine-grained mode colours (FLAIR modes get blue shades, tools get distinct)
MODE_COLORS = _PUB_MODE_COLORS

# Preferred display order — modes not listed here are appended alphabetically
DEFAULT_MODE_ORDER = list(MODE_COLORS.keys())


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
    )


def load_evaluation_files(input_files):
    """Load and concatenate evaluation TSVs into a DataFrame."""
    if pd is None:
        print("pandas required for concordance_plots.py", file=sys.stderr)
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


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _resolve_mode_order(df):
    """Return ordered list of modes present in the DataFrame."""
    present = set(df["transcriptome_mode"].unique())
    order = [m for m in DEFAULT_MODE_ORDER if m in present]
    for m in sorted(present):
        if m not in order:
            order.append(m)
    return order


def _mode_color(mode: str) -> str:
    if mode in MODE_COLORS:
        return MODE_COLORS[mode]
    idx = hash(mode) % len(PALETTE)
    return PALETTE[idx]


def _plot_metric_bars_single(df, output_path, mode_order, col_name, ylabel, title):
    """DEPRECATED — replaced by _plot_concordance_cleveland below."""
    pass


def _plot_concordance_cleveland(
    df, output_path, mode_order, col_5prime, col_3prime, xlabel,
    baseline_mode=None,
):
    """Cleveland-style horizontal dot plot comparing 5' and 3' for one metric.

    Methods are listed on the Y-axis (top-first). 5' and 3' values appear as
    paired dots connected by a thin segment. An optional vertical reference
    line is drawn at the baseline method's value.
    """
    from pub_style import W1

    methods = list(reversed(mode_order))
    y_pos = np.arange(len(methods))

    vals_5, vals_3 = [], []
    for mode in methods:
        subset = df[df["transcriptome_mode"] == mode]
        v5 = pd.to_numeric(subset[col_5prime], errors="coerce").dropna().mean() \
            if col_5prime in subset.columns else np.nan
        v3 = pd.to_numeric(subset[col_3prime], errors="coerce").dropna().mean() \
            if col_3prime in subset.columns else np.nan
        vals_5.append(v5 if not np.isnan(v5) else np.nan)
        vals_3.append(v3 if not np.isnan(v3) else np.nan)

    if all(np.isnan(v) for v in vals_5) and all(np.isnan(v) for v in vals_3):
        return

    fig_h = max(2.0, 0.38 * len(methods) + 0.8)
    fig, ax = plt.subplots(figsize=(W1, fig_h))

    # Baseline reference lines
    if baseline_mode and baseline_mode in methods:
        bl_idx = methods.index(baseline_mode)
        bl_5 = vals_5[bl_idx]
        bl_3 = vals_3[bl_idx]
        if not np.isnan(bl_5):
            ax.axvline(bl_5, color=PALETTE[4], linewidth=0.6, linestyle=":",
                       alpha=0.5, zorder=0)
        if not np.isnan(bl_3):
            ax.axvline(bl_3, color=PALETTE[5], linewidth=0.6, linestyle=":",
                       alpha=0.5, zorder=0)

    # Connecting segments
    for i, (a, b) in enumerate(zip(vals_5, vals_3)):
        if np.isnan(a) or np.isnan(b):
            continue
        ax.plot([a, b], [i, i], color="#bbbbbb", linewidth=0.8, zorder=1)

    # Dots
    ax.scatter(vals_5, y_pos, color=PALETTE[4], s=30, zorder=2,
               label="5\u2032", edgecolors="none")
    ax.scatter(vals_3, y_pos, color=PALETTE[5], s=30, zorder=2,
               label="3\u2032", edgecolors="none", marker="D")

    ax.set_yticks(y_pos)
    ax.set_yticklabels([_short_mode(m) for m in methods], fontsize=7)
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0),
              frameon=False, ncol=1, borderaxespad=0)
    style_ax(ax, xlabel=xlabel, faint_y_grid=True)
    fig.tight_layout()
    savefig(fig, output_path, dpi=300)


def _plot_signal_by_ends_single(df, output_path, mode_order, end_type, title):
    """Standalone line plot: median signal vs number of distinct ends per junction chain."""
    end_labels = ["1", "2", "3+"]
    signal_cols = [
        f"{end_type}_signal_median_1end",
        f"{end_type}_signal_median_2end",
        f"{end_type}_signal_median_3plus",
    ]

    x_pos = np.arange(len(end_labels))
    any_plotted = False

    fig, ax = plt.subplots(figsize=(3.5, 3.0))

    for mode in mode_order:
        subset = df[df["transcriptome_mode"] == mode]
        avgs = []
        for col in signal_cols:
            if col not in subset.columns:
                avgs.append(np.nan)
                continue
            vals = pd.to_numeric(subset[col], errors="coerce").dropna()
            avgs.append(vals.mean() if len(vals) > 0 else np.nan)

        if sum(1 for v in avgs if not np.isnan(v)) < 2:
            continue

        ax.plot(
            x_pos, avgs, marker="o", label=_short_mode(mode),
            color=_mode_color(mode), linewidth=1.5, markersize=5,
        )
        any_plotted = True

    if not any_plotted:
        plt.close(fig)
        return

    ax.set_xticks(x_pos)
    ax.set_xticklabels(end_labels)
    style_ax(ax, xlabel="Distinct ends per junction chain",
             ylabel="Median orthogonal signal")
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0),
             frameon=False, ncol=1, borderaxespad=0)
    fig.tight_layout()
    savefig(fig, output_path, dpi=300)


# ---------------------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------------------

def create_concordance_plots(df, output_dir, title_prefix=""):
    """Create individual concordance summary figures in output_dir."""
    mode_order = _resolve_mode_order(df)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Detect baseline mode for reference lines
    baseline_mode = None
    for m in mode_order:
        if "baseline" in m:
            baseline_mode = m
            break

    # Spearman rho — Cleveland dot plot (5' vs 3')
    col_5 = "5prime_concordance_spearman_rho"
    col_3 = "3prime_concordance_spearman_rho"
    if col_5 in df.columns or col_3 in df.columns:
        _plot_concordance_cleveland(
            df, output_dir / "concordance_spearman.png",
            mode_order, col_5, col_3,
            xlabel="Spearman \u03c1",
            baseline_mode=baseline_mode,
        )

    # Concordance index — Cleveland dot plot (5' vs 3')
    col_5 = "5prime_boundary_concordance_index"
    col_3 = "3prime_boundary_concordance_index"
    if col_5 in df.columns or col_3 in df.columns:
        _plot_concordance_cleveland(
            df, output_dir / "concordance_index.png",
            mode_order, col_5, col_3,
            xlabel="Concordance Index",
            baseline_mode=baseline_mode,
        )

    # Signal vs end count — keep as separate per-end line plots
    for end_type, end_label in [
        ("5prime", "5\u2032 TSS"),
        ("3prime", "3\u2032 TTS"),
    ]:
        _plot_signal_by_ends_single(
            df, output_dir / f"signal_by_ends_{end_type}.png",
            mode_order, end_type,
            title="",  # title suppressed
        )

    print(f"Saved concordance plots to {output_dir}")
    return True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Cross-run concordance & boundary selection summary plots",
    )
    parser.add_argument(
        "--input", nargs="+", required=True, help="Evaluation TSV files",
    )
    parser.add_argument("--output", required=True, help="Output directory for individual plots")
    parser.add_argument("--title-prefix", default="", help="Title prefix")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if pd is None:
        print("ERROR: pandas is required", file=sys.stderr)
        sys.exit(1)

    df = load_evaluation_files(args.input)
    if df is None or df.empty:
        print("No data loaded -- skipping plot", file=sys.stderr)
        sys.exit(0)

    needed = [
        "5prime_concordance_spearman_rho",
        "3prime_concordance_spearman_rho",
    ]
    if not any(c in df.columns for c in needed):
        print(
            "No concordance metrics found in evaluation data -- skipping plot",
            file=sys.stderr,
        )
        sys.exit(0)

    create_concordance_plots(df, args.output, args.title_prefix)


if __name__ == "__main__":
    main()
