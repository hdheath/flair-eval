#!/usr/bin/env python3
"""
sqanti_precision.py — SQANTI-stratified precision/recall dashboard.

Reads the combined evaluation TSV (one row per mode, containing both
precision/recall metrics and SQANTI structural-category counts) and
produces three complementary panels:

  C1a. sqanti_composition.png
       Stacked horizontal bar chart showing the fraction of isoforms in each
       SQANTI category (FSM, ISM, NIC, NNC, SEM, SEN) per mode.
       Modes sorted by NNC fraction (most noise on bottom) so the plot
       immediately shows which method produces the cleanest assemblies.

  C1b. sqanti_precision_scatter.png
       5' and 3' precision/recall scatter, one point per mode.
       Point colour encodes NNC fraction (diverging palette: low=blue,
       high=red) so you can see whether precision gains are real or just NNC
       suppression.

  C1c. sqanti_nnc_vs_precision.png
       NNC fraction (x) vs 5' and 3' precision (y) scatter — direct view of
       the NNC-precision trade-off across parameter variants.

Usage:
    python sqanti_precision.py \\
        --input combined_eval.tsv [combined_eval2.tsv ...] \\
        --output output_dir/ \\
        [--dataset LABEL]  # filter to one dataset if TSV has multiple
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside

apply_rc()

# ── SQANTI category colours (Okabe-Ito / publication-safe) ───────────────────
SQANTI_CATS = ["FSM", "ISM", "NIC", "NNC", "SEM", "SEN"]
SQANTI_COLORS = {
    "FSM": "#0072B2",   # blue  — full splice match (best)
    "ISM": "#56B4E9",   # sky   — incomplete splice match
    "NIC": "#009E73",   # green — novel in catalog (splice variant)
    "NNC": "#D55E00",   # vermillion — novel not in catalog (most suspect)
    "SEM": "#CC79A7",   # pink  — single-exon match
    "SEN": "#E69F00",   # amber — single-exon novel
}

NNC_CMAP = "RdBu_r"   # blue=low NNC, red=high NNC


# ── Data loading ──────────────────────────────────────────────────────────────

def load_data(paths: List[str | Path], dataset: Optional[str] = None) -> pd.DataFrame:
    dfs = []
    for p in paths:
        try:
            df = pd.read_csv(p, sep="\t")
            dfs.append(df)
        except Exception as e:
            print(f"WARNING: could not read {p}: {e}", file=sys.stderr)
    if not dfs:
        return pd.DataFrame()
    out = pd.concat(dfs, ignore_index=True)
    if dataset and "dataset" in out.columns:
        out = out[out["dataset"] == dataset].copy()
    # Ensure numeric columns
    num_cols = ["5prime_precision", "5prime_recall", "5prime_f1",
                 "3prime_precision", "3prime_recall", "3prime_f1",
                 "FSM", "ISM", "NIC", "NNC", "SEM", "SEN"]
    for c in num_cols:
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce")
    return out


def _sqanti_fractions(df: pd.DataFrame) -> pd.DataFrame:
    """Add per-mode SQANTI fraction columns and NNC fraction."""
    cats = [c for c in SQANTI_CATS if c in df.columns]
    total = df[cats].sum(axis=1).replace(0, np.nan)
    for c in cats:
        df[f"{c}_frac"] = df[c] / total
    if "NNC" in cats:
        df["nnc_frac"] = df["NNC"] / total
    else:
        df["nnc_frac"] = 0.0
    return df


# ── Plot C1a: stacked composition bar ────────────────────────────────────────

def plot_sqanti_composition(df: pd.DataFrame, output_path: Path):
    cats = [c for c in SQANTI_CATS if c in df.columns]
    if not cats:
        return

    # Sort modes by NNC fraction ascending (cleanest on top)
    df_sorted = df.sort_values("nnc_frac", ascending=False).reset_index(drop=True)
    modes = df_sorted["transcriptome_mode"].tolist()
    n = len(modes)

    fig, ax = plt.subplots(figsize=(W2 * 0.55, max(W1 * 0.4, n * 0.22 + 0.6)))

    lefts = np.zeros(n)
    handles = []
    for cat in cats:
        col = f"{cat}_frac"
        if col not in df_sorted.columns:
            continue
        vals = df_sorted[col].fillna(0).values
        bars = ax.barh(np.arange(n), vals, left=lefts,
                       color=SQANTI_COLORS[cat], edgecolor="none", alpha=0.9)
        lefts += vals
        handles.append(mpatches.Patch(facecolor=SQANTI_COLORS[cat],
                                       edgecolor="none", label=cat))

    ax.set_yticks(np.arange(n))
    ax.set_yticklabels(modes, fontsize=6)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of isoforms", fontsize=7)
    ax.invert_yaxis()
    style_ax(ax)
    ax.set_title("SQANTI structural-category composition", fontsize=7, pad=4)

    fig.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, 0.0),
               ncol=1, fontsize=6, frameon=False)
    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── Plot C1b: precision/recall scatter coloured by NNC fraction ───────────────

def plot_precision_scatter(df: pd.DataFrame, output_path: Path):
    if df.empty:
        return
    methods = df["transcriptome_mode"].tolist()
    styler = ModeStyler(methods)

    nnc_vals = df["nnc_frac"].fillna(0).values
    norm = mcolors.Normalize(vmin=0, vmax=max(nnc_vals.max(), 0.01))
    cmap = cm.get_cmap(NNC_CMAP)

    fig, axes = plt.subplots(1, 2, figsize=(W2, W1 * 0.75))

    for ax, end, label in [
        (axes[0], "5prime", "5\u2032 TSS"),
        (axes[1], "3prime", "3\u2032 TTS"),
    ]:
        p_col = f"{end}_precision"
        r_col = f"{end}_recall"
        for _, row in df.iterrows():
            p = row.get(p_col)
            r = row.get(r_col)
            nnc = row.get("nnc_frac", 0)
            m = row["transcriptome_mode"]
            if pd.isna(p) or pd.isna(r):
                continue
            color = cmap(norm(nnc))
            ax.scatter(r * 100, p * 100,
                       s=40, color=color, marker=styler.marker(m),
                       edgecolors="white", linewidth=0.5, alpha=0.9, zorder=2)
            ax.text(r * 100 + 0.5, p * 100, m.replace("TED-", "").replace("FLAIR-", ""),
                    fontsize=4, va="center", alpha=0.7)

        ax.plot([0, 100], [0, 100], "--", color="#cccccc", lw=0.8, zorder=0)
        ax.set_xlim(30, 105)
        ax.set_ylim(30, 105)
        style_ax(ax, xlabel=f"{label} Recall (%)", ylabel=f"{label} Precision (%)")
        ax.set_title(label, fontsize=7)

    # Colourbar for NNC fraction
    sm = cm.ScalarMappable(cmap=NNC_CMAP, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, shrink=0.6, pad=0.02)
    cbar.set_label("NNC fraction", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    fig.suptitle("Precision vs Recall coloured by NNC fraction", fontsize=7, y=1.01)
    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── Plot C1c: NNC fraction vs precision direct scatter ───────────────────────

def plot_nnc_vs_precision(df: pd.DataFrame, output_path: Path):
    if df.empty:
        return
    methods = df["transcriptome_mode"].tolist()
    styler = ModeStyler(methods)

    fig, axes = plt.subplots(1, 2, figsize=(W2, W1 * 0.7))

    for ax, end, label in [
        (axes[0], "5prime", "5\u2032 Precision"),
        (axes[1], "3prime", "3\u2032 Precision"),
    ]:
        p_col = f"{end}_precision"
        for _, row in df.iterrows():
            p   = row.get(p_col)
            nnc = row.get("nnc_frac", 0)
            m   = row["transcriptome_mode"]
            if pd.isna(p) or pd.isna(nnc):
                continue
            ax.scatter(nnc * 100, p * 100,
                       s=40, color=styler.color(m), marker=styler.marker(m),
                       edgecolors="white", linewidth=0.5, alpha=0.9, zorder=2)
            ax.text(nnc * 100 + 0.3, p * 100,
                    m.replace("TED-", "").replace("FLAIR-", ""),
                    fontsize=4, va="center", alpha=0.7)

        style_ax(ax, xlabel="NNC fraction (%)", ylabel=label + " (%)")
        ax.set_title(label, fontsize=7)

    # Shared legend
    handles = [styler.legend_handle(m, label=m, markersize=4) for m in methods]
    legend_outside(fig, handles=handles, loc="upper right",
                   bbox_to_anchor=(1.0, 1.0), ncol=1, fontsize=5)
    fig.suptitle("NNC fraction vs Precision (noise-precision trade-off)", fontsize=7, y=1.01)
    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", "-i", nargs="+", required=True,
                        help="Combined evaluation TSV file(s)")
    parser.add_argument("--output", "-o", required=True,
                        help="Output directory")
    parser.add_argument("--dataset", default=None,
                        help="Filter to a single dataset label (dataset column)")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    df = load_data(args.input, dataset=args.dataset)
    if df.empty:
        print("No data loaded — skipping", file=sys.stderr)
        sys.exit(1)

    df = _sqanti_fractions(df)

    if args.verbose:
        modes = df["transcriptome_mode"].unique()
        print(f"  {len(df)} rows, {len(modes)} modes", file=sys.stderr)
        print(f"  NNC range: {df['nnc_frac'].min():.2f}–{df['nnc_frac'].max():.2f}",
              file=sys.stderr)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_sqanti_composition(df, output_dir / "sqanti_composition.png")
    plot_precision_scatter(df, output_dir / "sqanti_precision_scatter.png")
    plot_nnc_vs_precision(df, output_dir / "sqanti_nnc_vs_precision.png")

    print(f"Saved SQANTI precision plots to {args.output}")


if __name__ == "__main__":
    main()
