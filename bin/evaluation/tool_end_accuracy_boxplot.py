#!/usr/bin/env python3
"""
Tool End-Accuracy Box Plot: cross-tool comparison of 5'/3' precision & recall.

X-axis: assembler tools, each with 4 sub-positions (5'P, 5'R, 3'P, 3'R).
Y-axis: metric value 0–100.
Dots colored by sample/dataset, with box plots overlaid.

Publication-quality formatting via ``pub_style``.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import numpy as np

from pub_style import (
    PALETTE, SAMPLE_PALETTE, ASSEMBLER_COLORS,
    style_ax, legend_outside, savefig, W2, HMAX,
)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ── Tool label mapping ──────────────────────────────────────────────────────

# Use transcriptome_mode directly as the tool label so each distinct mode
# (e.g. "default", "ted-2d", "ted-1d2d") gets its own column in the plot.

_MODE_PALETTE = [
    "#0072B2", "#56B4E9", "#E69F00", "#009E73",
    "#D55E00", "#CC79A7", "#999999", "#F0E442",
    "#332288", "#88CCEE", "#44AA99", "#117733",
]

METRIC_LABELS = ["5'P", "5'R", "3'P", "3'R"]
METRIC_COLS   = ["5prime_precision", "5prime_recall", "3prime_precision", "3prime_recall"]


# ── Data loading ────────────────────────────────────────────────────────────

def load_evaluation_files(input_files):
    dfs = []
    for f in input_files:
        try:
            df = pd.read_csv(f, sep='\t')
            dfs.append(df)
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}", file=sys.stderr)
    if not dfs:
        return None
    return pd.concat(dfs, ignore_index=True)


# ── Plotting ────────────────────────────────────────────────────────────────

def plot_tool_end_accuracy(df, output_dir, verbose=False):
    """Create the tool end-accuracy box + strip plot."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Use transcriptome_mode directly as the tool label (no collapsing)
    df = df.copy()
    df["tool"] = df["transcriptome_mode"]

    # Keep modes in the order they first appear in the data
    tools_present = list(dict.fromkeys(df["tool"]))
    if not tools_present:
        print("No recognised tools found in evaluation data.", file=sys.stderr)
        return

    # Assign colours from palette
    _tool_colors = {t: _MODE_PALETTE[i % len(_MODE_PALETTE)] for i, t in enumerate(tools_present)}

    # No averaging: each (mode, dataset) row is one data point
    agg_df = df[["tool", "dataset"] + METRIC_COLS].copy()

    # Melt to long format: tool, dataset, metric, value
    long = agg_df.melt(
        id_vars=["tool", "dataset"],
        value_vars=METRIC_COLS,
        var_name="metric_col",
        value_name="value",
    )
    long["metric"] = long["metric_col"].map(dict(zip(METRIC_COLS, METRIC_LABELS)))

    # Convert 0–1 fractions to 0–100 percentages
    long["value"] = long["value"] * 100

    # Dataset colour assignment
    datasets = sorted(long["dataset"].unique())
    ds_colors = {ds: SAMPLE_PALETTE[i % len(SAMPLE_PALETTE)] for i, ds in enumerate(datasets)}

    # Layout constants
    n_tools = len(tools_present)
    n_metrics = len(METRIC_LABELS)
    sub_width = 0.18          # width of each sub-position
    sub_gap   = 0.06          # gap between sub-positions within a tool
    tool_gap  = 0.6           # extra gap between tools
    group_span = n_metrics * sub_width + (n_metrics - 1) * sub_gap

    # Compute x-positions
    tool_centers = []
    sub_positions = {}  # (tool, metric_label) → x
    x = 0.0
    for tool in tools_present:
        center = x + group_span / 2
        tool_centers.append(center)
        for j, ml in enumerate(METRIC_LABELS):
            sub_positions[(tool, ml)] = x + j * (sub_width + sub_gap) + sub_width / 2
        x += group_span + tool_gap

    # Figure
    fig_w = max(W2, 0.8 * n_tools)
    fig, ax = plt.subplots(figsize=(fig_w, W2 * 0.55))

    # Box plots per (tool, metric)
    for tool in tools_present:
        for ml in METRIC_LABELS:
            xpos = sub_positions[(tool, ml)]
            vals = long.loc[(long["tool"] == tool) & (long["metric"] == ml), "value"].dropna()
            if vals.empty:
                continue
            bp = ax.boxplot(
                vals,
                positions=[xpos],
                widths=sub_width * 0.8,
                patch_artist=True,
                showfliers=False,
                zorder=2,
            )
            for patch in bp["boxes"]:
                patch.set_facecolor(_tool_colors.get(tool, "#333333"))
                patch.set_alpha(0.15)
                patch.set_edgecolor(_tool_colors.get(tool, "#333333"))
                patch.set_linewidth(0.6)
            for element in ("whiskers", "caps", "medians"):
                for line in bp[element]:
                    line.set_color(_tool_colors.get(tool, "#333333"))
                    line.set_alpha(0.4)
                    line.set_linewidth(0.6)

    # Strip plot (jittered dots)
    rng = np.random.default_rng(42)
    for tool in tools_present:
        for ml in METRIC_LABELS:
            xpos = sub_positions[(tool, ml)]
            subset = long.loc[(long["tool"] == tool) & (long["metric"] == ml)]
            for _, row in subset.iterrows():
                jitter = rng.uniform(-sub_width * 0.25, sub_width * 0.25)
                ax.scatter(
                    xpos + jitter,
                    row["value"],
                    c=ds_colors[row["dataset"]],
                    s=18,
                    edgecolors="white",
                    linewidths=0.3,
                    zorder=3,
                    alpha=0.9,
                )

    # X-axis: tool labels at centre, metric labels below
    ax.set_xticks([sub_positions[(t, ml)] for t in tools_present for ml in METRIC_LABELS])
    ax.set_xticklabels(METRIC_LABELS * n_tools, fontsize=5.5, rotation=0)

    # Secondary tool-name labels via text
    for tool, center in zip(tools_present, tool_centers):
        ax.text(
            center, -0.08, tool,
            transform=ax.get_xaxis_transform(),
            ha="center", va="top",
            fontsize=7, fontweight="bold",
            color=_tool_colors.get(tool, "#333333"),
        )

    # Axis styling
    ax.set_ylim(0, 105)
    style_ax(ax, faint_y_grid=True)

    # Vertical separator lines between tool groups
    for i in range(1, n_tools):
        sep_x = (tool_centers[i - 1] + tool_centers[i]) / 2
        ax.axvline(sep_x, color="#cccccc", linewidth=0.4, linestyle="--", zorder=0)

    # Pad left/right
    all_x = [sub_positions[(t, ml)] for t in tools_present for ml in METRIC_LABELS]
    ax.set_xlim(min(all_x) - sub_width, max(all_x) + sub_width)

    # Legend for datasets
    ds_handles = [
        mpatches.Patch(facecolor=ds_colors[ds], edgecolor="white", linewidth=0.3, label=ds)
        for ds in datasets
    ]
    legend_outside(ax, handles=ds_handles, labels=datasets, ncol=1, fontsize=6)

    # Increase bottom margin for tool labels
    fig.subplots_adjust(bottom=0.18)

    savefig(fig, output_dir / "tool_end_accuracy_boxplot.png", dpi=600)

    if verbose:
        print(f"Saved tool_end_accuracy_boxplot to {output_dir}", file=sys.stderr)
        print(f"  Tools: {tools_present}", file=sys.stderr)
        print(f"  Datasets: {datasets}", file=sys.stderr)
        print(f"  Total data points: {len(long)}", file=sys.stderr)


# ── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Tool end-accuracy box + strip plot across assemblers."
    )
    parser.add_argument(
        "--input", nargs="+", required=True,
        help="Evaluation TSV file(s).",
    )
    parser.add_argument(
        "--output", required=True,
        help="Output directory for plots.",
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Print progress info.",
    )
    args = parser.parse_args()

    df = load_evaluation_files(args.input)
    if df is None or df.empty:
        print("Error: no evaluation data loaded.", file=sys.stderr)
        sys.exit(1)

    plot_tool_end_accuracy(df, args.output, verbose=args.verbose)


if __name__ == "__main__":
    main()
