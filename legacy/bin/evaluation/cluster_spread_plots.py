#!/usr/bin/env python3
"""
cluster_spread_plots.py — TED cluster spread (IQR) diagnostic plots.

Uses the new tss_spread_iqr / tts_spread_iqr columns added to the TED log to
reveal whether accepted clusters are genuinely tight or whether TED is passing
widely scattered read-end groups that orthogonal signal does not support.

Produces:
  D2a. spread_violin_by_status.png
       Violin of TSS and TTS spread IQR split by cluster status
       (pass / reject) per mode.  A well-calibrated mode shows
       clearly lower IQR for accepted clusters.

Usage:
    python cluster_spread_plots.py \\
        --ted-log label1:log1.tsv label2:log2.tsv ... \\
        --output output_dir/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2

apply_rc()

STATUS_COLORS = {
    "pass":   "#0072B2",
    "reject": "#D55E00",
}

# ── Loading ────────────────────────────────────────────────────────────────────

def load_ted_log(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    for col in ("n_reads", "jc_n_reads_total",
                "tss_spread_iqr", "tts_spread_iqr",
                "TED_tss_depth", "TED_tts_depth", "TED_depth",
                "TED_confidence", "threshold_tss", "threshold_tts"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    # Keep only clustered rows (not noise, not merged)
    if "status" in df.columns:
        df = df[df["status"].isin(["pass", "reject"])].copy()
    return df


def load_all_logs(entries: List[str]) -> Dict[str, pd.DataFrame]:
    result: Dict[str, pd.DataFrame] = {}
    for entry in entries:
        if ":" not in entry:
            continue
        label, path = entry.split(":", 1)
        if not Path(path).exists():
            print(f"WARNING: {path} not found", file=sys.stderr)
            continue
        try:
            df = load_ted_log(path)
            if not df.empty:
                result[label] = df
        except Exception as e:
            print(f"WARNING: could not read {path}: {e}", file=sys.stderr)
    return result


# ── Plot D2a: violin of spread IQR by status ──────────────────────────────────

def plot_spread_violin(logs: Dict[str, pd.DataFrame], output_path: Path):
    methods = list(logs.keys())
    if not methods:
        return

    n     = len(methods)
    ncols = min(n, 4)
    nrows = (n + ncols - 1) // ncols

    # Two rows per panel: TSS IQR on top, TTS IQR on bottom
    fig, axes = plt.subplots(
        nrows * 2, ncols,
        figsize=(W2, W1 * 0.7 * nrows),
        squeeze=False,
    )

    for idx, m in enumerate(methods):
        df = logs[m]
        col_r = idx % ncols
        for row_off, (iqr_col, label) in enumerate([
            ("tss_spread_iqr", "TSS IQR (bp)"),
            ("tts_spread_iqr", "TTS IQR (bp)"),
        ]):
            ax = axes[idx // ncols * 2 + row_off][col_r]
            if iqr_col not in df.columns:
                ax.set_visible(False)
                continue

            parts_data, xtick_labels, colors = [], [], []
            for status in ("pass", "reject"):
                sub = df[df["status"] == status][iqr_col].dropna()
                # IQR of -1 means fallback (single read) — exclude
                sub = sub[sub >= 0]
                if len(sub) < 2:
                    continue
                parts_data.append(sub.values)
                xtick_labels.append(f"{status}\n(n={len(sub)})")
                colors.append(STATUS_COLORS[status])

            if not parts_data:
                ax.set_visible(False)
                continue

            vp = ax.violinplot(parts_data, positions=range(len(parts_data)),
                               showmedians=True, showextrema=False, widths=0.65)
            for body, col in zip(vp["bodies"], colors):
                body.set_facecolor(col)
                body.set_alpha(0.6)
                body.set_edgecolor("none")
            vp["cmedians"].set_color("#333333")
            vp["cmedians"].set_linewidth(1.2)

            ax.set_xticks(range(len(xtick_labels)))
            ax.set_xticklabels(xtick_labels, fontsize=5)
            style_ax(ax)
            ax.set_ylabel(label, fontsize=6)
            if row_off == 0:
                ax.text(0.97, 0.97, m, transform=ax.transAxes,
                        ha="right", va="top", fontsize=5, fontweight="bold")

    for idx in range(n, nrows * ncols):
        for row_off in range(2):
            axes[idx // ncols * 2 + row_off][idx % ncols].set_visible(False)

    fig.suptitle("Cluster spread IQR by acceptance status", fontsize=7)
    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── CLI ────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--ted-log", nargs="+", required=True,
                        help="label:path pairs for TED log TSV files")
    parser.add_argument("--output",  required=True)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    logs = load_all_logs(args.ted_log)
    if not logs:
        print("No TED log data loaded — skipping", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        for m, df in logs.items():
            sc = df["status"].value_counts().to_dict() if "status" in df.columns else {}
            print(f"  {m}: {len(df)} rows  {sc}", file=sys.stderr)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_spread_violin(logs, output_dir / "spread_violin_by_status.png")

    print(f"Saved cluster spread plots to {args.output}")


if __name__ == "__main__":
    main()
