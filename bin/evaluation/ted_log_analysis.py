#!/usr/bin/env python3
"""
ted_log_analysis.py — TED log: end spread analysis and CAGE peak width comparison.

Uses the tss_spread_iqr / tts_spread_iqr columns from the updated TED log to
answer two questions:

  E1. ted_end_spread_comparison.png
      Are TSS clusters tighter than TTS clusters for the same mode?
      Are accepted clusters meaningfully tighter than rejected ones?
      Paired violin plot: TSS IQR vs TTS IQR side-by-side per mode,
      split into pass / reject rows.  Also includes a summary scatter
      showing median TSS IQR vs median TTS IQR per mode — modes above
      the diagonal have more variable 5' ends than 3' ends.

  E2a. cage_width_vs_cluster_iqr.png
       For every accepted TSS cluster: CAGE peak width at its centroid (x)
       vs tss_spread_iqr (y), one panel per mode.  Points coloured by
       n_reads.  A tight 1:1 alignment suggests TED clusters faithfully
       track peak width; scatter above the line means TED is accepting
       wider clusters than the CAGE peak justifies.

  E2b. cage_width_vs_cluster_iqr_kde.png
       KDE density version of E2a across all modes combined — shows the
       population-level relationship without per-mode fragmentation.

  E2c. threshold_margin_violin.png
       For each mode, violin of (TED_tss_reality − threshold_tss) and
       (TED_tts_reality − threshold_tts) for pass clusters.  The "margin"
       shows how confidently clusters are accepted — thin margins suggest
       the threshold is near-optimally tight; fat margins suggest there is
       room to raise the threshold without losing many real isoforms.

Usage:
    python ted_log_analysis.py \\
        --ted-log label1:log1.tsv label2:log2.tsv ... \\
        --cage-peaks cage_peaks.bed \\
        --output output_dir/
"""

from __future__ import annotations

import argparse
import sys
from bisect import bisect_left
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.cm as cm

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside

apply_rc()

STATUS_COLORS = {"pass": "#0072B2", "reject": "#D55E00"}
MATCH_WINDOW  = 100   # bp to search for a CAGE peak near TSS centroid


# ── Loading ────────────────────────────────────────────────────────────────────

def load_ted_log(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    num_cols = [
        "n_reads", "jc_n_reads_total",
        "tss_spread_iqr", "tts_spread_iqr",
        "TED_tss_reality", "TED_tts_reality",
        "TED_confidence", "threshold_tss", "threshold_tts",
        "tss_pos", "tts_pos",
    ]
    for col in num_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
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


def _parse_strand(junc_id: str) -> str:
    """Extract strand from junc_id formatted as chrom:start-end:strand:Nj."""
    parts = str(junc_id).split(":")
    return parts[2] if len(parts) > 2 else "+"


def _parse_chrom(junc_id: str) -> str:
    parts = str(junc_id).split(":")
    return parts[0] if parts else "?"


# ── CAGE peak index ────────────────────────────────────────────────────────────

class CagePeakIndex:
    """Sorted BED6 CAGE peaks indexed for O(log n) nearest-peak lookup."""

    def __init__(self, path: str | Path):
        # {(chrom, strand): sorted list of (start, end)}
        self._peaks: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
        with open(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 6:
                    continue
                chrom  = cols[0]
                start  = int(cols[1])
                end    = int(cols[2])
                strand = cols[5]
                key    = (chrom, strand)
                self._peaks.setdefault(key, []).append((start, end))
        for k in self._peaks:
            self._peaks[k].sort()

    def nearest_peak_width(
        self, chrom: str, pos: int, strand: str, window: int = MATCH_WINDOW
    ) -> Optional[int]:
        """Return width (bp) of the nearest CAGE peak within *window* of *pos*."""
        key = (chrom, strand)
        peaks = self._peaks.get(key)
        if not peaks:
            return None
        starts = [p[0] for p in peaks]
        idx = bisect_left(starts, pos)
        best_width = None
        best_dist  = window + 1
        for i in (idx - 1, idx):
            if 0 <= i < len(peaks):
                s, e = peaks[i]
                # distance: 0 if pos inside peak, else distance to nearest edge
                dist = 0 if s <= pos < e else min(abs(pos - s), abs(pos - (e - 1)))
                if dist < best_dist:
                    best_dist  = dist
                    best_width = e - s
        return best_width if best_dist <= window else None


def _add_cage_width(
    df: pd.DataFrame, cage_index: CagePeakIndex
) -> pd.DataFrame:
    """Add cage_peak_width column to accepted-cluster rows."""
    widths = []
    for _, row in df.iterrows():
        junc_id = str(row.get("junc_id", ""))
        chrom   = _parse_chrom(junc_id)
        strand  = _parse_strand(junc_id)
        tss_pos = row.get("tss_pos", -1)
        if pd.isna(tss_pos) or tss_pos < 0:
            widths.append(np.nan)
            continue
        w = cage_index.nearest_peak_width(chrom, int(tss_pos), strand)
        widths.append(float(w) if w is not None else np.nan)
    df = df.copy()
    df["cage_peak_width"] = widths
    return df


# ── Plot E1: TSS vs TTS spread violin + summary scatter ──────────────────────

def plot_end_spread_comparison(
    logs: Dict[str, pd.DataFrame], output_path: Path
):
    methods = list(logs.keys())
    if not methods:
        return

    n     = len(methods)
    ncols = min(n, 4)
    nrows = (n + ncols - 1) // ncols

    # Main panel: paired violins TSS IQR (blue) vs TTS IQR (orange) for pass clusters
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(W2, W1 * 0.75 * nrows),
                             squeeze=False)

    col_tss = "#0072B2"
    col_tts = "#E69F00"

    for idx, m in enumerate(methods):
        ax  = axes[idx // ncols][idx % ncols]
        df  = logs[m]
        sub = df[df["status"] == "pass"].copy()

        tss_iqr = sub["tss_spread_iqr"].dropna() if "tss_spread_iqr" in sub.columns else pd.Series(dtype=float)
        tts_iqr = sub["tts_spread_iqr"].dropna() if "tts_spread_iqr" in sub.columns else pd.Series(dtype=float)
        tss_iqr = tss_iqr[tss_iqr >= 0]
        tts_iqr = tts_iqr[tts_iqr >= 0]

        if len(tss_iqr) < 2 and len(tts_iqr) < 2:
            ax.set_visible(False)
            continue

        positions = []
        parts_data = []
        colors = []
        xtick_pos = []
        xtick_labels = []
        pos = 0
        for vals, col, lbl in [
            (tss_iqr.values, col_tss, "TSS"),
            (tts_iqr.values, col_tts, "TTS"),
        ]:
            if len(vals) >= 2:
                parts_data.append(vals)
                positions.append(pos)
                colors.append(col)
                xtick_pos.append(pos)
                xtick_labels.append(f"{lbl}\nn={len(vals)}")
            pos += 1

        if not parts_data:
            ax.set_visible(False)
            continue

        vp = ax.violinplot(parts_data, positions=positions,
                           showmedians=True, showextrema=False, widths=0.6)
        for body, c in zip(vp["bodies"], colors):
            body.set_facecolor(c)
            body.set_alpha(0.65)
            body.set_edgecolor("none")
        vp["cmedians"].set_color("#333333")
        vp["cmedians"].set_linewidth(1.2)

        ax.set_xticks(xtick_pos)
        ax.set_xticklabels(xtick_labels, fontsize=6)
        style_ax(ax)
        ax.text(0.97, 0.97, m, transform=ax.transAxes,
                ha="right", va="top", fontsize=5, fontweight="bold")
        if idx % ncols == 0:
            ax.set_ylabel("Spread IQR (bp)", fontsize=6)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    handles = [
        mpatches.Patch(facecolor=col_tss, label="TSS IQR"),
        mpatches.Patch(facecolor=col_tts, label="TTS IQR"),
    ]
    fig.legend(handles=handles, loc="lower right",
               bbox_to_anchor=(1.0, 0.0), fontsize=6, frameon=False)
    fig.suptitle("TSS vs TTS cluster spread IQR (accepted clusters)", fontsize=7)
    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


def plot_end_spread_summary_scatter(
    logs: Dict[str, pd.DataFrame], output_path: Path
):
    """Scatter of median TSS IQR vs median TTS IQR per mode."""
    methods = list(logs.keys())
    if not methods:
        return
    styler = ModeStyler(methods)

    fig, ax = plt.subplots(figsize=(W1 * 1.1, W1 * 1.1))

    for m in methods:
        df  = logs[m]
        sub = df[df["status"] == "pass"]
        tss_iqr = sub["tss_spread_iqr"].dropna() if "tss_spread_iqr" in sub.columns else pd.Series(dtype=float)
        tts_iqr = sub["tts_spread_iqr"].dropna() if "tts_spread_iqr" in sub.columns else pd.Series(dtype=float)
        tss_iqr = tss_iqr[tss_iqr >= 0]
        tts_iqr = tts_iqr[tts_iqr >= 0]
        if tss_iqr.empty or tts_iqr.empty:
            continue
        med_tss = float(np.median(tss_iqr))
        med_tts = float(np.median(tts_iqr))
        ax.scatter(med_tss, med_tts,
                   s=40, color=styler.color(m), marker=styler.marker(m),
                   edgecolors="white", linewidth=0.5, zorder=2)
        ax.text(med_tss + 0.3, med_tts,
                m.replace("TED-", "").replace("FLAIR-", ""),
                fontsize=4, va="center", alpha=0.8)

    lim = ax.get_xlim()[1]
    ax.plot([0, lim], [0, lim], "--", color="#cccccc", lw=0.8, zorder=0)
    style_ax(ax, xlabel="Median TSS spread IQR (bp)",
             ylabel="Median TTS spread IQR (bp)")
    ax.set_title("TSS vs TTS spread per mode\n(above diagonal = TSS wider than TTS)",
                 fontsize=6)
    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── Plot E2a/b: CAGE peak width vs TSS cluster IQR ───────────────────────────

def plot_cage_width_vs_cluster_iqr(
    logs: Dict[str, pd.DataFrame], output_path: Path
):
    methods = [m for m in logs if "cage_peak_width" in logs[m].columns]
    if not methods:
        return

    n     = len(methods)
    ncols = min(n, 4)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(W2, W1 * 0.8 * nrows),
                             squeeze=False)

    for idx, m in enumerate(methods):
        ax  = axes[idx // ncols][idx % ncols]
        df  = logs[m]
        sub = df[
            (df["status"] == "pass") &
            (df["tss_spread_iqr"] >= 0) &
            df["cage_peak_width"].notna()
        ].copy()
        if sub.empty:
            ax.set_visible(False)
            continue

        depth = sub["n_reads"].fillna(1).values
        norm  = mcolors.LogNorm(
            vmin=max(depth.min(), 1), vmax=max(depth.max(), 2))
        colors = cm.get_cmap("plasma")(norm(np.maximum(depth, 1)))

        ax.scatter(sub["cage_peak_width"].values,
                   sub["tss_spread_iqr"].values,
                   s=2, c=colors, alpha=0.4, edgecolors="none",
                   rasterized=True)

        # 1:1 reference line
        lim = max(sub["cage_peak_width"].max(), sub["tss_spread_iqr"].max()) * 1.05
        ax.plot([0, lim], [0, lim], "--", color="#cccccc", lw=0.8, zorder=0)

        style_ax(ax)
        ax.text(0.97, 0.97, m, transform=ax.transAxes,
                ha="right", va="top", fontsize=5, fontweight="bold")
        if idx >= n - ncols:
            ax.set_xlabel("CAGE peak width (bp)", fontsize=6)
        if idx % ncols == 0:
            ax.set_ylabel("TSS cluster IQR (bp)", fontsize=6)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    sm = cm.ScalarMappable(cmap="plasma",
                           norm=mcolors.LogNorm(vmin=1, vmax=100))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, shrink=0.5, pad=0.02)
    cbar.set_label("n reads in cluster", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    fig.suptitle("CAGE peak width vs TED TSS cluster IQR\n"
                 "(points above 1:1 = cluster wider than CAGE peak)", fontsize=7)
    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


def plot_cage_width_vs_iqr_kde(
    logs: Dict[str, pd.DataFrame], output_path: Path
):
    """All-mode combined KDE density scatter."""
    # Collect all pass rows that have both fields
    frames = []
    for m, df in logs.items():
        if "cage_peak_width" not in df.columns:
            continue
        sub = df[
            (df["status"] == "pass") &
            (df["tss_spread_iqr"] >= 0) &
            df["cage_peak_width"].notna()
        ][["cage_peak_width", "tss_spread_iqr"]].copy()
        sub["mode"] = m
        frames.append(sub)

    if not frames:
        return

    all_df = pd.concat(frames, ignore_index=True)
    xs = all_df["cage_peak_width"].values
    ys = all_df["tss_spread_iqr"].values

    # Simple 2D histogram as proxy for KDE
    fig, ax = plt.subplots(figsize=(W1 * 1.2, W1 * 1.2))
    h = ax.hist2d(xs, ys, bins=50,
                  norm=mcolors.LogNorm(), cmap="YlOrRd")
    fig.colorbar(h[3], ax=ax, label="Count (log)", shrink=0.7)
    lim = max(xs.max(), ys.max()) * 1.05
    ax.plot([0, lim], [0, lim], "--", color="#333333", lw=0.8, zorder=3)
    style_ax(ax, xlabel="CAGE peak width (bp)",
             ylabel="TSS cluster IQR (bp)")
    ax.set_title("All modes: CAGE width vs cluster IQR", fontsize=7)
    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── Plot E2c: threshold margin violin ────────────────────────────────────────

def plot_threshold_margin_violin(
    logs: Dict[str, pd.DataFrame], output_path: Path
):
    """Violin of (reality − threshold) for accepted clusters, per end, per mode."""
    methods = list(logs.keys())
    if not methods:
        return

    n     = len(methods)
    ncols = min(n, 4)
    nrows = (n + ncols - 1) // ncols

    col_tss = "#0072B2"
    col_tts = "#E69F00"

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(W2, W1 * 0.75 * nrows),
                             squeeze=False)

    for idx, m in enumerate(methods):
        ax  = axes[idx // ncols][idx % ncols]
        df  = logs[m]
        sub = df[df["status"] == "pass"]

        required = ["TED_tss_reality", "TED_tts_reality",
                    "threshold_tss", "threshold_tts"]
        if not all(c in sub.columns for c in required):
            ax.set_visible(False)
            continue

        tss_margin = (sub["TED_tss_reality"] - sub["threshold_tss"]).dropna()
        tts_margin = (sub["TED_tts_reality"] - sub["threshold_tts"]).dropna()

        if len(tss_margin) < 2 and len(tts_margin) < 2:
            ax.set_visible(False)
            continue

        parts_data, positions, colors = [], [], []
        xtick_pos, xtick_labels = [], []
        pos = 0
        for vals, col, lbl in [
            (tss_margin.values, col_tss, "TSS margin"),
            (tts_margin.values, col_tts, "TTS margin"),
        ]:
            if len(vals) >= 2:
                parts_data.append(vals)
                positions.append(pos)
                colors.append(col)
                xtick_pos.append(pos)
                xtick_labels.append(f"{lbl}\nn={len(vals)}")
            pos += 1

        if not parts_data:
            ax.set_visible(False)
            continue

        vp = ax.violinplot(parts_data, positions=positions,
                           showmedians=True, showextrema=False, widths=0.6)
        for body, c in zip(vp["bodies"], colors):
            body.set_facecolor(c)
            body.set_alpha(0.65)
            body.set_edgecolor("none")
        vp["cmedians"].set_color("#333333")
        vp["cmedians"].set_linewidth(1.2)

        ax.axhline(0, color="#888888", lw=0.7, linestyle=":", zorder=0)
        ax.set_xticks(xtick_pos)
        ax.set_xticklabels(xtick_labels, fontsize=6)
        style_ax(ax)
        ax.text(0.97, 0.97, m, transform=ax.transAxes,
                ha="right", va="top", fontsize=5, fontweight="bold")
        if idx % ncols == 0:
            ax.set_ylabel("Reality − threshold (margin)", fontsize=6)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    handles = [
        mpatches.Patch(facecolor=col_tss, label="TSS margin"),
        mpatches.Patch(facecolor=col_tts, label="TTS margin"),
    ]
    fig.legend(handles=handles, loc="lower right",
               bbox_to_anchor=(1.0, 0.0), fontsize=6, frameon=False)
    fig.suptitle("Threshold margin for accepted clusters\n"
                 "(TED_reality − threshold; wider = more headroom to raise threshold)",
                 fontsize=7)
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
    parser.add_argument("--cage-peaks", default=None,
                        help="BED6 CAGE peaks file for peak-width comparison (E2)")
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
            sc = df["status"].value_counts().to_dict() if "status" in df.columns else {}
            print(f"  {m}: {len(df)} rows  {sc}", file=sys.stderr)

    # Enrich with CAGE peak widths if provided
    if args.cage_peaks and Path(args.cage_peaks).exists():
        if args.verbose:
            print(f"  Loading CAGE peaks from {args.cage_peaks}…", file=sys.stderr)
        cage_index = CagePeakIndex(args.cage_peaks)
        for m in list(logs.keys()):
            pass_df = logs[m][logs[m]["status"] == "pass"].copy()
            enriched = _add_cage_width(pass_df, cage_index)
            # Merge cage_peak_width back into full df
            logs[m] = logs[m].copy()
            logs[m]["cage_peak_width"] = np.nan
            logs[m].loc[enriched.index, "cage_peak_width"] = enriched["cage_peak_width"].values
    elif args.cage_peaks:
        print(f"WARNING: {args.cage_peaks} not found — skipping E2 plots",
              file=sys.stderr)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # E1
    plot_end_spread_comparison(logs, output_dir / "ted_end_spread_comparison.png")
    plot_end_spread_summary_scatter(logs,
                                    output_dir / "ted_end_spread_summary_scatter.png")

    # E2
    has_cage = args.cage_peaks and Path(args.cage_peaks).exists()
    if has_cage:
        plot_cage_width_vs_cluster_iqr(logs,
                                       output_dir / "cage_width_vs_cluster_iqr.png")
        plot_cage_width_vs_iqr_kde(logs,
                                   output_dir / "cage_width_vs_cluster_iqr_kde.png")

    # E2c (no CAGE peaks needed — uses reality/threshold from log)
    plot_threshold_margin_violin(logs,
                                  output_dir / "threshold_margin_violin.png")

    print(f"Saved TED log analysis plots to {args.output}")


if __name__ == "__main__":
    main()
