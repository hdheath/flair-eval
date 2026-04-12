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

  D2b. spread_ecdf_pass_vs_reject.png
       Empirical CDF of spread IQR for pass vs reject clusters.
       The further the curves diverge, the more spread is
       discriminating acceptance.

  D2c. spread_vs_signal.png
       Scatter of tss_spread_iqr (x) vs CAGE signal at TSS (y)
       for every accepted cluster, one panel per mode.
       Coloured by jc_n_reads_total (depth-aware).
       Reveals whether wide, low-signal clusters slip through.

  D2d. spread_vs_n_reads.png
       Scatter of n_reads (x, log) vs tss/tts spread IQR (y) for
       pass clusters — shows whether deeper clusters converge.

Usage:
    python cluster_spread_plots.py \\
        --ted-log label1:log1.tsv label2:log2.tsv ... \\
        --cage-plus cage_plus.bg --cage-minus cage_minus.bg \\
        --qs-plus qs_plus.bg   --qs-minus qs_minus.bg \\
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
import matplotlib.colors as mcolors
import matplotlib.cm as cm

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
    from signal_utils import BedGraphTrack, load_signal_tracks
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
    from evaluation.signal_utils import BedGraphTrack, load_signal_tracks

apply_rc()

STATUS_COLORS = {
    "pass":   "#0072B2",
    "reject": "#D55E00",
}
SIG_WINDOW = 50   # bp window for signal query at cluster centroid


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


def _query_signal_at_pos(
    chrom: str, pos: int, strand: str,
    cage_p: BedGraphTrack, cage_m: BedGraphTrack,
    qs_p: BedGraphTrack, qs_m: BedGraphTrack,
    window: int = SIG_WINDOW,
) -> Tuple[float, float]:
    """Return (tss_signal, tts_signal) at genomic position."""
    if strand == "+":
        tss_sig = cage_p.query(chrom, pos - window, pos + window)
        tts_sig = qs_p.query(chrom, pos - window, pos + window)
    else:
        tss_sig = cage_m.query(chrom, pos - window, pos + window)
        tts_sig = qs_m.query(chrom, pos - window, pos + window)
    return tss_sig, tts_sig


def _add_signal(
    df: pd.DataFrame,
    cage_p: BedGraphTrack, cage_m: BedGraphTrack,
    qs_p: BedGraphTrack, qs_m: BedGraphTrack,
) -> pd.DataFrame:
    """Add tss_signal and tts_signal columns by querying signal tracks."""
    if "tss_pos" not in df.columns or "tts_pos" not in df.columns:
        return df
    # Extract chrom and strand from junc_id: chr22:start-end:strand:Nj
    def _parse_junc(junc_id):
        parts = str(junc_id).split(":")
        chrom = parts[0] if parts else "?"
        strand = parts[2] if len(parts) > 2 else "+"
        return chrom, strand

    tss_sigs, tts_sigs = [], []
    for _, row in df.iterrows():
        chrom, strand = _parse_junc(row.get("junc_id", ""))
        tss_pos = int(row.get("tss_pos", -1))
        tts_pos = int(row.get("tts_pos", -1))
        if tss_pos < 0 or tts_pos < 0:
            tss_sigs.append(np.nan)
            tts_sigs.append(np.nan)
            continue
        ts, tt = _query_signal_at_pos(
            chrom, tss_pos, strand, cage_p, cage_m, qs_p, qs_m)
        tss_sigs.append(ts)
        tts_sigs.append(tt)
    df = df.copy()
    df["tss_signal"] = tss_sigs
    df["tts_signal"] = tts_sigs
    return df


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


# ── Plot D2b: ECDF of spread IQR pass vs reject ───────────────────────────────

def plot_spread_ecdf(logs: Dict[str, pd.DataFrame], output_path: Path):
    methods = list(logs.keys())
    if not methods:
        return

    n     = len(methods)
    ncols = min(n, 4)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(W2, W1 * 0.7 * nrows),
                             squeeze=False)

    for idx, m in enumerate(methods):
        ax  = axes[idx // ncols][idx % ncols]
        df  = logs[m]
        col = "tss_spread_iqr"  # use TSS; TTS available but same story
        if col not in df.columns:
            ax.set_visible(False)
            continue

        for status, ls in [("pass", "-"), ("reject", "--")]:
            sub = df[df["status"] == status][col].dropna()
            sub = sub[sub >= 0]
            if sub.empty:
                continue
            xs   = np.sort(sub.values)
            ecdf = np.arange(1, len(xs) + 1) / len(xs)
            ax.step(xs, ecdf, where="post",
                    color=STATUS_COLORS[status], lw=1.1,
                    linestyle=ls, label=status, alpha=0.85)

        style_ax(ax)
        ax.text(0.97, 0.05, m, transform=ax.transAxes,
                ha="right", va="bottom", fontsize=5, fontweight="bold")
        if idx >= n - ncols:
            ax.set_xlabel("TSS spread IQR (bp)", fontsize=6)
        if idx % ncols == 0:
            ax.set_ylabel("ECDF", fontsize=6)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    handles = [
        mpatches.Patch(facecolor=STATUS_COLORS["pass"],   label="pass"),
        mpatches.Patch(facecolor=STATUS_COLORS["reject"], label="reject"),
    ]
    fig.legend(handles=handles, loc="lower right",
               bbox_to_anchor=(1.0, 0.0), fontsize=6, frameon=False)
    fig.suptitle("Cluster TSS spread IQR: pass vs reject (ECDF)", fontsize=7)
    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── Plot D2c: spread vs signal (accepted clusters) ────────────────────────────

def plot_spread_vs_signal(logs: Dict[str, pd.DataFrame], output_path: Path):
    methods = [m for m in logs if "tss_signal" in logs[m].columns]
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
        sub = df[(df["status"] == "pass") &
                 (df["tss_spread_iqr"] >= 0) &
                 df["tss_signal"].notna()].copy()
        if sub.empty:
            ax.set_visible(False)
            continue

        depth = sub["jc_n_reads_total"].fillna(0).values
        norm  = mcolors.LogNorm(
            vmin=max(depth.min(), 1), vmax=max(depth.max(), 2))
        colors = cm.get_cmap("viridis")(norm(np.maximum(depth, 1)))

        ax.scatter(sub["tss_spread_iqr"].values,
                   sub["tss_signal"].values + 1e-4,
                   s=2, c=colors, alpha=0.4, edgecolors="none",
                   rasterized=True)
        ax.set_yscale("log")
        style_ax(ax)
        ax.text(0.97, 0.97, m, transform=ax.transAxes,
                ha="right", va="top", fontsize=5, fontweight="bold")
        if idx >= n - ncols:
            ax.set_xlabel("TSS spread IQR (bp)", fontsize=6)
        if idx % ncols == 0:
            ax.set_ylabel("CAGE signal at TSS", fontsize=6)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    # Colorbar
    sm = cm.ScalarMappable(cmap="viridis",
                           norm=mcolors.LogNorm(vmin=1, vmax=100))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, shrink=0.5, pad=0.02)
    cbar.set_label("JC total reads", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    fig.suptitle("Accepted cluster TSS spread vs CAGE signal", fontsize=7)
    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── Plot D2d: spread vs n_reads (accepted clusters) ──────────────────────────

def plot_spread_vs_nreads(logs: Dict[str, pd.DataFrame], output_path: Path):
    methods = list(logs.keys())
    if not methods:
        return
    styler = ModeStyler(methods)

    fig, axes = plt.subplots(1, 2, figsize=(W2, W1 * 0.7))
    for ax, iqr_col, label in [
        (axes[0], "tss_spread_iqr", "TSS spread IQR (bp)"),
        (axes[1], "tts_spread_iqr", "TTS spread IQR (bp)"),
    ]:
        for m in methods:
            df  = logs[m]
            sub = df[(df["status"] == "pass") &
                     (df[iqr_col] >= 0)].copy() if iqr_col in df.columns else pd.DataFrame()
            if sub.empty:
                continue
            ax.scatter(sub["n_reads"].values + 0.5,
                       sub[iqr_col].values,
                       s=1.5, alpha=0.25, color=styler.color(m),
                       edgecolors="none", rasterized=True)

        ax.set_xscale("log")
        style_ax(ax, xlabel="n reads in cluster (log)", ylabel=label)

    handles = [styler.legend_handle(m, label=m, markersize=4) for m in methods]
    legend_outside(fig, handles=handles, loc="upper right",
                   bbox_to_anchor=(1.0, 1.0), ncol=1, fontsize=5)
    fig.suptitle("Do deeper clusters converge? Spread vs read count (pass only)",
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
    parser.add_argument("--cage-plus",  default=None)
    parser.add_argument("--cage-minus", default=None)
    parser.add_argument("--qs-plus",    default=None)
    parser.add_argument("--qs-minus",   default=None)
    parser.add_argument("--output",     required=True)
    parser.add_argument("--verbose",    action="store_true")
    args = parser.parse_args()

    logs = load_all_logs(args.ted_log)
    if not logs:
        print("No TED log data loaded — skipping", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        for m, df in logs.items():
            sc = df["status"].value_counts().to_dict() if "status" in df.columns else {}
            print(f"  {m}: {len(df)} rows  {sc}", file=sys.stderr)

    # Optionally enrich with signal
    has_signal = all(x is not None for x in
                     [args.cage_plus, args.cage_minus,
                      args.qs_plus, args.qs_minus])
    if has_signal:
        if args.verbose:
            print("  Loading signal tracks…", file=sys.stderr)
        cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
            args.cage_plus, args.cage_minus, args.qs_plus, args.qs_minus)
        for m in list(logs.keys()):
            logs[m] = _add_signal(logs[m], cage_p, cage_m, qs_p, qs_m)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_spread_violin(logs, output_dir / "spread_violin_by_status.png")
    plot_spread_ecdf(logs, output_dir / "spread_ecdf_pass_vs_reject.png")
    plot_spread_vs_nreads(logs, output_dir / "spread_vs_n_reads.png")
    if has_signal:
        plot_spread_vs_signal(logs, output_dir / "spread_vs_signal.png")
    elif args.verbose:
        print("  Skipping spread_vs_signal.png (no signal tracks provided)",
              file=sys.stderr)

    print(f"Saved cluster spread plots to {args.output}")


if __name__ == "__main__":
    main()
