#!/usr/bin/env python3
"""
signal_read_support.py — Signal at isoform ends across assembly modes.

Asks: do different modes place isoform ends on stronger signal peaks?

Produces:
  signal_distribution_5prime.png
      ECDF of CAGE signal at TSS per mode (non-zero isoforms only),
      plus fraction with zero signal shown in legend.

  signal_distribution_3prime.png
      Same for dRNA signal at TTS.

  signal_zero_fraction.png
      Grouped bar chart: fraction of isoforms with zero signal at TSS and TTS
      per mode — direct cross-mode quality comparison.

  signal_zero_count.png
      Same as signal_zero_fraction.png but Y-axis is raw isoform count,
      making absolute scale differences between methods visible.

  signal_quantiles.png
      Median + IQR (P25–P75) of non-zero signal per mode, dot-and-range plot,
      for TSS and TTS side by side — which mode best concentrates ends on
      real signal peaks?

Usage:
    python signal_read_support.py \\
        --bed label1:bed1.bed label2:bed2.bed ... \\
        --cage-plus cage_plus.bg --cage-minus cage_minus.bg \\
        --qs-plus qs_plus.bg   --qs-minus qs_minus.bg \\
        --output output_dir/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
    from signal_utils import parse_isoforms, load_signal_tracks, isoform_signal
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
    from evaluation.signal_utils import parse_isoforms, load_signal_tracks, isoform_signal

apply_rc()


# ── Data building ─────────────────────────────────────────────────────────────

def build_data(
    beds_by_method: Dict[str, List[dict]],
    cage_p, cage_m, qs_p, qs_m,
) -> Dict[str, dict]:
    """Return {method: {tss: np.ndarray, tts: np.ndarray}} of signal values."""
    data: Dict[str, dict] = {}
    for method, isos in beds_by_method.items():
        tss_vals, tts_vals = [], []
        for iso in isos:
            tss_sig, tts_sig = isoform_signal(iso, cage_p, cage_m, qs_p, qs_m)
            tss_vals.append(tss_sig)
            tts_vals.append(tts_sig)
        data[method] = {
            "tss": np.array(tss_vals, dtype=float),
            "tts": np.array(tts_vals, dtype=float),
        }
    return data


# ── Plot 1: ECDF of signal (non-zero isoforms) ───────────────────────────────

def plot_signal_ecdf(
    data: Dict[str, dict],
    sig_key: str,
    sig_label: str,
    output_path: Path,
    styler: ModeStyler,
):
    methods = [m for m in data if len(data[m][sig_key]) > 0]
    if not methods:
        return

    fig, ax = plt.subplots(figsize=(W2 * 0.55, W1 * 0.75))
    handles = []

    for m in methods:
        vals = data[m][sig_key]
        n_total = len(vals)
        nonzero = vals[vals > 0]
        zero_pct = (vals == 0).sum() / n_total * 100 if n_total > 0 else 0

        if len(nonzero) == 0:
            continue

        sorted_v = np.sort(nonzero)
        ecdf_y   = np.arange(1, len(sorted_v) + 1) / len(sorted_v)

        line, = ax.plot(
            sorted_v, ecdf_y,
            color=styler.color(m),
            linewidth=1.0,
            alpha=0.85,
        )
        handles.append(styler.legend_handle(
            m,
            label=f"{m}  ({zero_pct:.0f}% zero)",
            markersize=5,
        ))

    ax.set_xscale("log")
    style_ax(ax,
             xlabel=f"{sig_label} signal at isoform end (non-zero)",
             ylabel="Cumulative fraction of isoforms")
    ax.set_ylim(0, 1)

    legend_outside(fig, handles=handles, loc="lower right",
                   bbox_to_anchor=(1.0, 0.0), ncol=1, fontsize=5,
                   title="mode  (% zero-signal)")
    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── Plot 2: Zero-signal fraction per mode ────────────────────────────────────

def plot_zero_fraction(
    data: Dict[str, dict],
    output_path: Path,
    styler: ModeStyler,
):
    methods = list(data.keys())
    if not methods:
        return

    tss_fracs = []
    tts_fracs = []
    for m in methods:
        tss = data[m]["tss"]
        tts = data[m]["tts"]
        tss_fracs.append((tss == 0).sum() / len(tss) if len(tss) > 0 else 0)
        tts_fracs.append((tts == 0).sum() / len(tts) if len(tts) > 0 else 0)

    x     = np.arange(len(methods))
    width = 0.35

    fig, ax = plt.subplots(figsize=(W2, W1 * 0.65))

    bars_tss = ax.bar(x - width / 2, tss_fracs, width,
                      color=[styler.color(m) for m in methods],
                      edgecolor="none", alpha=0.9, label="TSS (CAGE)")
    bars_tts = ax.bar(x + width / 2, tts_fracs, width,
                      color=[styler.color(m) for m in methods],
                      edgecolor="none", alpha=0.45, label="TTS (dRNA)")

    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=35, ha="right", fontsize=6)
    ax.set_ylim(0, 1)
    style_ax(ax, ylabel="Fraction of isoforms with zero signal")

    # Hatch legend to distinguish TSS vs TTS bars
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="grey", alpha=0.9,  label="TSS (CAGE)"),
        Patch(facecolor="grey", alpha=0.45, label="TTS (dRNA)"),
    ]
    ax.legend(handles=legend_handles, fontsize=6, frameon=False,
              loc="upper right")

    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── Plot 2b: Raw count of zero-signal isoforms per mode ──────────────────────

def plot_zero_count(
    data: Dict[str, dict],
    output_path: Path,
    styler: ModeStyler,
):
    methods = list(data.keys())
    if not methods:
        return

    tss_counts = []
    tts_counts = []
    for m in methods:
        tss = data[m]["tss"]
        tts = data[m]["tts"]
        tss_counts.append(int((tss == 0).sum()))
        tts_counts.append(int((tts == 0).sum()))

    x     = np.arange(len(methods))
    width = 0.35

    fig, ax = plt.subplots(figsize=(W2, W1 * 0.65))

    ax.bar(x - width / 2, tss_counts, width,
           color=[styler.color(m) for m in methods],
           edgecolor="none", alpha=0.9, label="TSS (CAGE)")
    ax.bar(x + width / 2, tts_counts, width,
           color=[styler.color(m) for m in methods],
           edgecolor="none", alpha=0.45, label="TTS (dRNA)")

    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=35, ha="right", fontsize=6)
    style_ax(ax, ylabel="Number of isoforms with zero signal")

    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="grey", alpha=0.9,  label="TSS (CAGE)"),
        Patch(facecolor="grey", alpha=0.45, label="TTS (dRNA)"),
    ]
    ax.legend(handles=legend_handles, fontsize=6, frameon=False,
              loc="upper right")

    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── Plot 3: Median + IQR of non-zero signal per mode ─────────────────────────

def plot_signal_quantiles(
    data: Dict[str, dict],
    output_path: Path,
    styler: ModeStyler,
):
    methods = list(data.keys())
    if not methods:
        return

    fig, axes = plt.subplots(1, 2, figsize=(W2, W1 * 0.7), sharey=False)

    for ax, sig_key, sig_label in [
        (axes[0], "tss", "CAGE signal at TSS"),
        (axes[1], "tts", "dRNA signal at TTS"),
    ]:
        medians, p25s, p75s = [], [], []
        valid_methods = []
        for m in methods:
            nz = data[m][sig_key]
            nz = nz[nz > 0]
            if len(nz) < 5:
                continue
            medians.append(np.median(nz))
            p25s.append(np.percentile(nz, 25))
            p75s.append(np.percentile(nz, 75))
            valid_methods.append(m)

        if not valid_methods:
            ax.set_visible(False)
            continue

        y = np.arange(len(valid_methods))
        for i, m in enumerate(valid_methods):
            ax.plot([p25s[i], p75s[i]], [i, i],
                    color=styler.color(m), linewidth=2.0, alpha=0.6,
                    solid_capstyle="round")
            ax.scatter([medians[i]], [i],
                       color=styler.color(m), s=20, zorder=3,
                       edgecolors="white", linewidths=0.4)

        ax.set_yticks(y)
        ax.set_yticklabels(valid_methods, fontsize=6)
        ax.set_xscale("log")
        style_ax(ax, xlabel=sig_label + " (non-zero, median ± IQR)")
        ax.invert_yaxis()

    fig.tight_layout(pad=0.4)
    savefig(fig, output_path)


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bed",        nargs="+", required=True,
                        help="label:path pairs for BED12 isoform files")
    parser.add_argument("--cage-plus",  required=True)
    parser.add_argument("--cage-minus", required=True)
    parser.add_argument("--qs-plus",    required=True)
    parser.add_argument("--qs-minus",   required=True)
    parser.add_argument("--output",     required=True)
    parser.add_argument("--verbose",    action="store_true")
    args = parser.parse_args()

    beds_by_method: Dict[str, List[dict]] = {}
    for entry in args.bed:
        if ":" not in entry:
            continue
        label, path = entry.split(":", 1)
        if not Path(path).exists():
            print(f"WARNING: {path} not found", file=sys.stderr)
            continue
        isos = parse_isoforms(path)
        if isos:
            beds_by_method[label] = isos
            if args.verbose:
                print(f"  {label}: {len(isos)} isoforms", file=sys.stderr)

    if not beds_by_method:
        print("No isoform data — skipping", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print("  Loading signal tracks...", file=sys.stderr)
    cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
        args.cage_plus, args.cage_minus, args.qs_plus, args.qs_minus,
    )

    if args.verbose:
        print("  Computing per-isoform signal...", file=sys.stderr)
    data = build_data(beds_by_method, cage_p, cage_m, qs_p, qs_m)

    styler = ModeStyler(list(beds_by_method.keys()))

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_signal_ecdf(data, "tss", "CAGE",
                     output_dir / "signal_distribution_5prime.png", styler)
    plot_signal_ecdf(data, "tts", "dRNA",
                     output_dir / "signal_distribution_3prime.png", styler)
    plot_zero_fraction(data, output_dir / "signal_zero_fraction.png", styler)
    plot_zero_count(data, output_dir / "signal_zero_count.png", styler)
    plot_signal_quantiles(data, output_dir / "signal_quantiles.png", styler)

    print(f"Saved signal distribution plots to {args.output}")


if __name__ == "__main__":
    main()
