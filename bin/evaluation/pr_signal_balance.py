#!/usr/bin/env python3
"""
pr_signal_balance.py — Precision / recall / boundary-signal balance plots.

Three complementary views of the three-way trade-off between isoform-end
precision, recall, and the fraction of isoforms whose ends land on genuine
signal (CAGE at TSS, dRNA at TTS).

Plot 1 — Three-axis scatter  (pr_signal_scatter.png)
    X = 5′ precision,  Y = 3′ precision.
    Point size ∝ zero-signal fraction at TSS (larger = worse).
    Point colour = F1 (mean of 5′ and 3′ F1), RdYlGn colourmap.
    One point per mode, averaged over datasets.

Plot 2 — Pareto frontier  (pr_signal_pareto.png)
    Two panels side by side:  left = 5′ end,  right = 3′ end.
    X = end precision,  Y = zero-signal fraction at that end.
    Pareto-optimal modes (high precision AND low zero-signal) annotated.
    Colour = mode family via ModeStyler.

Plot 3 — Dead-zone bar + F1 overlay  (pr_signal_deadzone.png)
    Grouped horizontal bars sorted by composite signal quality.
    Left bar  = zero-signal fraction at TSS (solid fill).
    Right bar = zero-signal fraction at TTS (hatched fill).
    Right-axis line = mean F1 (avg of 5′ and 3′), one dot per mode.

Usage:
    python pr_signal_balance.py \\
        --tsv label1:eval1.tsv [label2:eval2.tsv ...] \\
        --bed label1:bed1.bed [label2:bed2.bed ...] \\
        --cage-plus cage.plus.bg --cage-minus cage.minus.bg \\
        --qs-plus qs.plus.bg --qs-minus qs.minus.bg \\
        --output outdir/
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
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

try:
    from pub_style import (apply_rc, style_ax, savefig, legend_outside,
                           W1, W2, ModeStyler)
    from signal_utils import parse_isoforms, load_signal_tracks, isoform_signal
except ImportError:
    from evaluation.pub_style import (apply_rc, style_ax, savefig,
                                      legend_outside, W1, W2, ModeStyler)
    from evaluation.signal_utils import (parse_isoforms, load_signal_tracks,
                                         isoform_signal)

apply_rc()

# ── Constants ─────────────────────────────────────────────────────────────────

# Minimum isoforms for a mode to be included
MIN_ISOFORMS = 20
# Marker size scaling for scatter plot: base + scale * zero_frac
_SZ_BASE  = 12
_SZ_SCALE = 180


# ── Data loading ──────────────────────────────────────────────────────────────

def _parse_label_path(entries: List[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for e in entries:
        if ":" not in e:
            continue
        label, path = e.split(":", 1)
        out[label] = path
    return out


def load_pr_stats(tsv_paths: Dict[str, str]) -> pd.DataFrame:
    """Load and concatenate per-sample combined evaluation TSVs.

    Returns one row per (dataset_label, transcriptome_mode) with averaged
    5′/3′ precision, recall and F1.
    """
    frames = []
    for label, path in tsv_paths.items():
        p = Path(path)
        if not p.exists():
            print(f"WARNING: {path} not found", file=sys.stderr)
            continue
        df = pd.read_csv(p, sep="\t")
        df["_label"] = label
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def compute_zero_signal(
    beds_by_mode: Dict[str, List[dict]],
    cage_p, cage_m, qs_p, qs_m,
) -> Dict[str, Dict[str, float]]:
    """Return {mode: {tss_zero_frac, tts_zero_frac}} from BED isoforms."""
    result: Dict[str, Dict[str, float]] = {}
    for mode, isos in beds_by_mode.items():
        if len(isos) < MIN_ISOFORMS:
            continue
        tss_zero = tts_zero = 0
        for iso in isos:
            tss_sig, tts_sig = isoform_signal(iso, cage_p, cage_m, qs_p, qs_m)
            if tss_sig == 0:
                tss_zero += 1
            if tts_sig == 0:
                tts_zero += 1
        n = len(isos)
        result[mode] = {
            "tss_zero_frac": tss_zero / n,
            "tts_zero_frac": tts_zero / n,
        }
    return result


def merge_stats(
    pr_df: pd.DataFrame,
    zero_by_mode: Dict[str, Dict[str, float]],
) -> pd.DataFrame:
    """Average P/R across datasets, join with zero-signal fractions."""
    if pr_df.empty:
        return pd.DataFrame()

    needed = ["transcriptome_mode", "5prime_precision", "5prime_recall",
              "5prime_f1", "3prime_precision", "3prime_recall", "3prime_f1"]
    for col in needed:
        if col not in pr_df.columns:
            return pd.DataFrame()

    agg = (
        pr_df[needed]
        .groupby("transcriptome_mode")
        .mean(numeric_only=True)
        .reset_index()
    )

    zero_df = pd.DataFrame([
        {"transcriptome_mode": m, **v}
        for m, v in zero_by_mode.items()
    ])

    merged = agg.merge(zero_df, on="transcriptome_mode", how="inner")
    merged["mean_f1"] = (merged["5prime_f1"] + merged["3prime_f1"]) / 2
    # composite for sorting: precision both ends + signal quality
    merged["_composite"] = (
        0.30 * merged["5prime_precision"]
        + 0.30 * merged["3prime_precision"]
        + 0.15 * (1 - merged["tss_zero_frac"])
        + 0.15 * (1 - merged["tts_zero_frac"])
        + 0.10 * merged["mean_f1"]
    )
    return merged.sort_values("_composite", ascending=False).reset_index(drop=True)


# ── Plot 1: Three-axis scatter ────────────────────────────────────────────────

def plot_scatter(df: pd.DataFrame, styler: ModeStyler, output_path: Path) -> None:
    """5′ prec vs 3′ prec; size = TSS zero-signal frac; colour = mean F1."""
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(W1 * 1.55, W1 * 1.35))

    cmap  = plt.get_cmap("RdYlGn")
    f1min = df["mean_f1"].min()
    f1max = df["mean_f1"].max()
    f1rng = f1max - f1min if f1max > f1min else 1e-6

    # Normalise zero-signal to its own range so small differences are visible
    zmin = df["tss_zero_frac"].min()
    zmax = df["tss_zero_frac"].max()
    zrng = zmax - zmin if zmax > zmin else 1e-6

    for _, row in df.iterrows():
        mode  = row["transcriptome_mode"]
        x     = row["5prime_precision"]
        y     = row["3prime_precision"]
        zfrac = row["tss_zero_frac"]
        f1    = row["mean_f1"]

        norm_f1  = (f1 - f1min) / f1rng
        norm_z   = (zfrac - zmin) / zrng   # 0 = best, 1 = worst
        color    = cmap(norm_f1)
        size     = _SZ_BASE + _SZ_SCALE * norm_z

        ax.scatter(x, y, s=size, color=color,
                   edgecolors=styler.color(mode),
                   linewidths=0.8, zorder=3, alpha=0.88,
                   marker=styler.marker(mode))
        # short label offset
        ax.annotate(
            _short(mode), (x, y),
            fontsize=4.5, ha="left", va="bottom",
            xytext=(3, 3), textcoords="offset points",
            color="#333333",
        )

    # Colourbar for F1
    sm = plt.cm.ScalarMappable(
        cmap=cmap,
        norm=mcolors.Normalize(vmin=f1min, vmax=f1max),
    )
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Mean F1 (5′ + 3′) / 2", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    # Size legend (relative zero-signal — scaled to data range)
    for norm, lbl in [(0.0, "lowest"), (0.5, "mid"), (1.0, "highest")]:
        ax.scatter([], [], s=_SZ_BASE + _SZ_SCALE * norm,
                   color="lightgrey", edgecolors="#555555",
                   linewidths=0.6, label=f"{lbl} zero-signal")
    ax.legend(title="TSS zero-signal\n(relative, point size)", title_fontsize=5,
              fontsize=5, loc="lower right", frameon=False,
              handletextpad=0.4, borderpad=0.5)

    style_ax(ax,
             xlabel="5′ precision (TSS)",
             ylabel="3′ precision (TTS)")
    ax.set_xlim(left=max(0, df["5prime_precision"].min() - 0.05))
    ax.set_ylim(bottom=max(0, df["3prime_precision"].min() - 0.05))

    # Diagonal reference (equal 5′/3′ precision)
    lims = [
        max(ax.get_xlim()[0], ax.get_ylim()[0]),
        min(ax.get_xlim()[1], ax.get_ylim()[1]),
    ]
    ax.plot(lims, lims, ls="--", lw=0.5, color="#999999", zorder=0)

    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── Plot 2: Pareto frontier ───────────────────────────────────────────────────

def _pareto_front(xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    """Return boolean mask of Pareto-optimal points (max x, min y)."""
    n = len(xs)
    mask = np.ones(n, dtype=bool)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if xs[j] >= xs[i] and ys[j] <= ys[i] and (xs[j] > xs[i] or ys[j] < ys[i]):
                mask[i] = False
                break
    return mask


def plot_pareto(df: pd.DataFrame, styler: ModeStyler, output_path: Path) -> None:
    """Precision vs zero-signal Pareto for 5′ and 3′ ends."""
    if df.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(W2, W1 * 1.1), sharey=False)

    panels = [
        (axes[0], "5prime_precision", "tss_zero_frac", "5′ precision (TSS)", "Zero-signal fraction at TSS"),
        (axes[1], "3prime_precision", "tts_zero_frac", "3′ precision (TTS)", "Zero-signal fraction at TTS"),
    ]

    for ax, prec_col, zero_col, xlabel, ylabel in panels:
        xs = df[prec_col].values
        ys = df[zero_col].values
        pareto = _pareto_front(xs, ys)

        for i, row in df.iterrows():
            mode  = row["transcriptome_mode"]
            x, y  = row[prec_col], row[zero_col]
            is_p  = pareto[i]
            ax.scatter(x, y,
                       s=22 if is_p else 14,
                       color=styler.color(mode),
                       marker=styler.marker(mode),
                       edgecolors="white" if is_p else "none",
                       linewidths=0.6,
                       zorder=4 if is_p else 3,
                       alpha=0.95 if is_p else 0.70)
            if is_p:
                ax.annotate(
                    _short(mode), (x, y),
                    fontsize=4.5, ha="left", va="top",
                    xytext=(3, -3), textcoords="offset points",
                    color="#222222", fontweight="bold",
                )

        # Shade Pareto-optimal region lightly
        if pareto.any():
            px = xs[pareto]
            py = ys[pareto]
            order = np.argsort(px)
            ax.step(
                np.r_[px[order], px[order][-1]],
                np.r_[py[order][0], py[order]],
                where="post", color="#009E73", lw=0.8,
                ls="--", alpha=0.6, zorder=2,
                label="Pareto frontier",
            )

        ax.set_ylim(bottom=0)
        style_ax(ax, xlabel=xlabel, ylabel=ylabel)
        ax.legend(fontsize=5, frameon=False, loc="upper left")

    # Shared mode legend
    handles = [styler.legend_handle(m, label=_short(m), markersize=5)
               for m in df["transcriptome_mode"]]
    legend_outside(fig, handles=handles, loc="lower center",
                   bbox_to_anchor=(0.5, -0.18), ncol=4, fontsize=5)
    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── Plot 3: Dead-zone bar + F1 overlay ───────────────────────────────────────

def plot_deadzone(df: pd.DataFrame, styler: ModeStyler, output_path: Path) -> None:
    """Horizontal grouped bars: zero-signal fractions; F1 annotated as text.

    Uses a single x-axis for zero-signal fractions only — no twin-axis
    confusion.  Mean F1 is written as a right-aligned label on each bar row.
    """
    if df.empty:
        return

    # Sort by composite (best at top)
    sdf = df.sort_values("_composite", ascending=True).reset_index(drop=True)
    modes  = sdf["transcriptome_mode"].tolist()
    n      = len(modes)
    y      = np.arange(n)
    height = 0.35

    fig, ax = plt.subplots(figsize=(W2 * 0.72, max(W1 * 0.5, n * 0.22 + 0.5)))

    x_max = max(sdf[["tss_zero_frac", "tts_zero_frac"]].values.max() * 1.35 + 0.01,
                0.06)   # ensure room for small values

    for i, row in sdf.iterrows():
        mode  = row["transcriptome_mode"]
        tss_z = row["tss_zero_frac"]
        tts_z = row["tts_zero_frac"]
        f1    = row["mean_f1"]
        col   = styler.color(mode)

        ax.barh(y[i] + height / 2, tss_z, height,
                color=col, alpha=0.90, edgecolor="none")
        ax.barh(y[i] - height / 2, tts_z, height,
                color=col, alpha=0.45, edgecolor=col,
                linewidth=0.4, hatch="///")

        # Annotate F1 at the right edge of the plot area
        ax.text(x_max, y[i], f"F1={f1:.2f}",
                ha="right", va="center", fontsize=5, color="#333333")

    ax.set_yticks(y)
    ax.set_yticklabels([_short(m) for m in modes], fontsize=6)
    ax.set_xlim(0, x_max)
    style_ax(ax, xlabel="Fraction of isoforms with zero boundary signal")

    legend_handles = [
        Patch(facecolor="grey", alpha=0.90, edgecolor="none",
              label="TSS  (CAGE, solid)"),
        Patch(facecolor="grey", alpha=0.45, edgecolor="grey",
              linewidth=0.4, hatch="///",
              label="TTS  (dRNA, hatched)"),
    ]
    ax.legend(handles=legend_handles, fontsize=5, frameon=False,
              loc="lower right")

    fig.tight_layout(pad=0.5)
    savefig(fig, output_path)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _short(mode: str) -> str:
    return (mode
            .replace("TED-", "")
            .replace("FLAIR-", "FL-")
            .replace("isoquant_", "IQ-"))


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--tsv", nargs="+", required=True,
                        help="label:path pairs for combined evaluation TSVs "
                             "(one per dataset/sample)")
    parser.add_argument("--bed", nargs="+", required=True,
                        help="label:path pairs for BED12 isoform files "
                             "(one per transcriptome mode)")
    parser.add_argument("--cage-plus",  required=True)
    parser.add_argument("--cage-minus", required=True)
    parser.add_argument("--qs-plus",    required=True)
    parser.add_argument("--qs-minus",   required=True)
    parser.add_argument("--output",     required=True)
    parser.add_argument("--verbose",    action="store_true")
    args = parser.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    # ── Load P/R stats from TSVs ──────────────────────────────────────────────
    tsv_paths = _parse_label_path(args.tsv)
    pr_df = load_pr_stats(tsv_paths)
    if pr_df.empty:
        print("No P/R data — exiting", file=sys.stderr)
        sys.exit(1)
    if args.verbose:
        print(f"  Loaded P/R for {pr_df['transcriptome_mode'].nunique()} modes "
              f"across {len(tsv_paths)} datasets", file=sys.stderr)

    # ── Load BED isoforms ─────────────────────────────────────────────────────
    bed_paths = _parse_label_path(args.bed)
    beds_by_mode: Dict[str, List[dict]] = {}
    for label, path in bed_paths.items():
        if not Path(path).exists():
            print(f"WARNING: {path} not found", file=sys.stderr)
            continue
        isos = parse_isoforms(path)
        if isos:
            beds_by_mode[label] = isos
            if args.verbose:
                print(f"  {label}: {len(isos)} isoforms", file=sys.stderr)

    if not beds_by_mode:
        print("No BED isoform data — exiting", file=sys.stderr)
        sys.exit(1)

    # ── Load signal tracks ────────────────────────────────────────────────────
    if args.verbose:
        print("  Loading signal tracks...", file=sys.stderr)
    cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
        args.cage_plus, args.cage_minus, args.qs_plus, args.qs_minus,
    )

    # ── Compute zero-signal fractions ─────────────────────────────────────────
    if args.verbose:
        print("  Computing zero-signal fractions...", file=sys.stderr)
    zero_by_mode = compute_zero_signal(beds_by_mode, cage_p, cage_m, qs_p, qs_m)
    if not zero_by_mode:
        print("No signal data computed — exiting", file=sys.stderr)
        sys.exit(1)

    # ── Merge into one table ──────────────────────────────────────────────────
    df = merge_stats(pr_df, zero_by_mode)
    if df.empty:
        print("No modes with both P/R and signal data — exiting", file=sys.stderr)
        sys.exit(1)
    if args.verbose:
        print(f"  Final table: {len(df)} modes", file=sys.stderr)

    styler = ModeStyler(df["transcriptome_mode"].tolist())

    # ── Produce plots ─────────────────────────────────────────────────────────
    plot_scatter(df, styler, out / "pr_signal_scatter.png")
    plot_pareto(df,  styler, out / "pr_signal_pareto.png")
    plot_deadzone(df, styler, out / "pr_signal_deadzone.png")

    print(f"Saved P/R × signal balance plots to {args.output}")


if __name__ == "__main__":
    main()
