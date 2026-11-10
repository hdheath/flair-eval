#!/usr/bin/env python3
"""
depth_calibration.py — TED depth-score calibration plots.

Reads TED log files and shows how read support (n_reads per cluster) relates
to acceptance/rejection, and how the depth score component (TED_depth) is
distributed across outcomes.  Answers the question: "Is the depth saturation
point set appropriately — are we rejecting high-read-count clusters or passing
low-read-count ones?"

Produces:
  D1a. depth_read_count_ecdf.png
       Empirical CDF of n_reads for accepted vs rejected clusters, per mode.
       A vertical dashed line marks the depth_saturation value if provided.
       If the two CDFs overlap strongly at low read counts, the saturation is
       pulling depth scores too flat — consider a lower saturation.

  D1b. depth_score_violin.png
       Violin / strip plot of TED_depth score for each outcome
       (pass / reject / noise) per mode.  Reveals whether the depth component
       is informative (good separation) or flat.

  D1c. depth_accept_rate_by_bin.png
       Bar chart: fraction of clusters ACCEPTED in each n_reads bin, per mode.
       A monotone-increasing curve is healthy — if acceptance rate is low even
       at high read counts, the depth threshold may be too strict, or model/
       annot scores are overriding depth.

Usage:
    python depth_calibration.py \\
        --ted-log label1:log1.tsv label2:log2.tsv ... \\
        --depth-saturation 50 \\
        --output output_dir/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside

apply_rc()

# ── Read-count bins ───────────────────────────────────────────────────────────
READ_BINS      = [0, 1, 2, 4, 9, 19, 49, float("inf")]
READ_LABELS    = ["1", "2", "3–4", "5–9", "10–19", "20–49", "50+"]

STATUS_COLORS = {
    "pass":   "#0072B2",   # blue
    "reject": "#D55E00",   # vermillion
    "noise":  "#CC79A7",   # pink
    "merged": "#56B4E9",   # sky blue
}

# ── Loading ───────────────────────────────────────────────────────────────────

def load_ted_log(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    for col in ("n_reads", "TED_depth", "TED_confidence"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_all_logs(entries: List[str]) -> Dict[str, pd.DataFrame]:
    """Parse 'label:path' entries, return {label: dataframe}."""
    result: Dict[str, pd.DataFrame] = {}
    for entry in entries:
        if ":" not in entry:
            print(f"WARNING: skipping malformed entry '{entry}' (expected label:path)",
                  file=sys.stderr)
            continue
        label, path = entry.split(":", 1)
        if not Path(path).exists():
            print(f"WARNING: {path} not found", file=sys.stderr)
            continue
        try:
            result[label] = load_ted_log(path)
        except Exception as e:
            print(f"WARNING: could not read {path}: {e}", file=sys.stderr)
    return result


def _read_bin_idx(n: int) -> int:
    for i in range(len(READ_BINS) - 1):
        if READ_BINS[i] < n <= READ_BINS[i + 1]:
            return i
    return len(READ_LABELS) - 1


# ── Plot D1a: ECDF of n_reads per status ─────────────────────────────────────

def plot_read_count_ecdf(
    logs: Dict[str, pd.DataFrame],
    output_path: Path,
    depth_saturation: Optional[int] = None,
):
    methods = list(logs.keys())
    if not methods:
        return

    n      = len(methods)
    ncols  = min(n, 4)
    nrows  = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(W2, W1 * 0.75 * nrows),
                             squeeze=False)

    for idx, m in enumerate(methods):
        ax  = axes[idx // ncols][idx % ncols]
        df  = logs[m]
        if "n_reads" not in df.columns or "status" not in df.columns:
            ax.set_visible(False)
            continue

        for status in ("pass", "reject"):
            sub = df[df["status"] == status]["n_reads"].dropna()
            if sub.empty:
                continue
            xs   = np.sort(sub.values)
            ecdf = np.arange(1, len(xs) + 1) / len(xs)
            ax.step(xs, ecdf, where="post",
                    color=STATUS_COLORS.get(status, "#888888"),
                    lw=1.2, label=status, alpha=0.85)

        if depth_saturation is not None:
            ax.axvline(depth_saturation, color="#555555",
                       linestyle="--", lw=0.9, alpha=0.7,
                       label=f"sat={depth_saturation}")

        ax.set_xscale("log")
        ax.set_xlim(left=0.8)
        style_ax(ax)
        ax.text(0.97, 0.05, m, transform=ax.transAxes,
                ha="right", va="bottom", fontsize=5, fontweight="bold")
        if idx >= n - ncols:
            ax.set_xlabel("n reads per cluster", fontsize=6)
        if idx % ncols == 0:
            ax.set_ylabel("ECDF", fontsize=6)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    # Shared legend
    handles = [
        mpatches.Patch(facecolor=STATUS_COLORS["pass"],   label="pass"),
        mpatches.Patch(facecolor=STATUS_COLORS["reject"], label="reject"),
    ]
    if depth_saturation is not None:
        from matplotlib.lines import Line2D
        handles.append(Line2D([0], [0], color="#555555", linestyle="--",
                               lw=0.9, label=f"sat={depth_saturation}"))
    fig.legend(handles=handles, loc="lower right",
               bbox_to_anchor=(1.0, 0.0), fontsize=6, frameon=False)
    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── Plot D1b: TED_depth violin by outcome ────────────────────────────────────

def plot_depth_score_violin(
    logs: Dict[str, pd.DataFrame],
    output_path: Path,
):
    methods = list(logs.keys())
    if not methods:
        return

    target_statuses = ["pass", "reject", "noise"]

    n     = len(methods)
    ncols = min(n, 4)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(W2, W1 * 0.75 * nrows),
                             squeeze=False)

    for idx, m in enumerate(methods):
        ax  = axes[idx // ncols][idx % ncols]
        df  = logs[m]
        if "TED_depth" not in df.columns or "status" not in df.columns:
            ax.set_visible(False)
            continue

        parts_data = []
        xtick_labels = []
        colors = []
        for s in target_statuses:
            sub = df[df["status"] == s]["TED_depth"].dropna().values
            if len(sub) < 2:
                continue
            parts_data.append(sub)
            xtick_labels.append(s)
            colors.append(STATUS_COLORS.get(s, "#888888"))

        if not parts_data:
            ax.set_visible(False)
            continue

        vp = ax.violinplot(parts_data, positions=range(len(parts_data)),
                           showmedians=True, showextrema=False, widths=0.7)
        for body, col in zip(vp["bodies"], colors):
            body.set_facecolor(col)
            body.set_alpha(0.6)
            body.set_edgecolor("none")
        vp["cmedians"].set_color("#333333")
        vp["cmedians"].set_linewidth(1.2)

        ax.set_xticks(range(len(xtick_labels)))
        ax.set_xticklabels(xtick_labels, fontsize=6)
        ax.set_ylim(-0.05, 1.05)
        style_ax(ax)
        ax.text(0.97, 0.97, m, transform=ax.transAxes,
                ha="right", va="top", fontsize=5, fontweight="bold")
        if idx % ncols == 0:
            ax.set_ylabel("TED depth score", fontsize=6)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── Plot D1c: accept-rate by n_reads bin ─────────────────────────────────────

def plot_accept_rate_by_bin(
    logs: Dict[str, pd.DataFrame],
    output_path: Path,
):
    methods = list(logs.keys())
    if not methods:
        return
    styler  = ModeStyler(methods)
    n_bins  = len(READ_LABELS)
    x       = np.arange(n_bins)
    width   = 0.8 / max(len(methods), 1)

    fig, ax = plt.subplots(figsize=(W2, W1 * 0.7))
    handles = []
    for i, m in enumerate(methods):
        df = logs[m]
        if "n_reads" not in df.columns or "status" not in df.columns:
            continue
        rates = []
        for b in range(n_bins):
            lo = READ_BINS[b]
            hi = READ_BINS[b + 1]
            sub = df[(df["n_reads"] > lo) & (df["n_reads"] <= hi)]
            if len(sub) == 0:
                rates.append(np.nan)
            else:
                rates.append((sub["status"] == "pass").sum() / len(sub))

        offset = (i - len(methods) / 2 + 0.5) * width
        ax.bar(x + offset, rates, width * 0.9,
               color=styler.color(m), edgecolor="none", alpha=0.85,
               label=m)
        handles.append(styler.legend_handle(m, label=m, markersize=5))

    ax.set_xticks(x)
    ax.set_xticklabels(READ_LABELS, fontsize=7)
    ax.set_ylim(0, 1)
    style_ax(ax,
             xlabel="Reads per cluster (n_reads)",
             ylabel="Fraction of clusters accepted (pass)")
    ax.set_title("TED acceptance rate by read-support bin", fontsize=7)
    legend_outside(fig, handles=handles, loc="upper left",
                   bbox_to_anchor=(1.0, 1.0), ncol=1, fontsize=6)
    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--ted-log", nargs="+", required=True,
                        help="label:path pairs for TED log TSV files")
    parser.add_argument("--depth-saturation", type=int, default=None,
                        help="Depth saturation value used during assembly "
                             "(drawn as reference line in ECDF plot)")
    parser.add_argument("--output", required=True,
                        help="Output directory")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    logs = load_all_logs(args.ted_log)
    if not logs:
        print("No TED log data loaded — skipping", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        for m, df in logs.items():
            statuses = df["status"].value_counts().to_dict() if "status" in df.columns else {}
            print(f"  {m}: {len(df)} rows  {statuses}", file=sys.stderr)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_read_count_ecdf(logs, output_dir / "depth_read_count_ecdf.png",
                         depth_saturation=args.depth_saturation)
    plot_depth_score_violin(logs, output_dir / "depth_score_violin.png")
    plot_accept_rate_by_bin(logs, output_dir / "depth_accept_rate_by_bin.png")

    print(f"Saved depth calibration plots to {args.output}")


if __name__ == "__main__":
    main()
