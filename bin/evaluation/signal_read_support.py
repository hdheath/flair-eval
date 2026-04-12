#!/usr/bin/env python3
"""
signal_read_support.py — Signal × read support scatter plot.

For each isoform, plots orthogonal signal strength at the TSS (CAGE) and TTS
(dRNA) against read support, revealing whether low-support isoforms are also
low-signal (noise) or have genuine signal (real biology being dropped by
frac_support or threshold tuning).

Produces:
  B1a. signal_x_read_support_5prime.png
       x = read support (log), y = CAGE signal at TSS.
       Points coloured by zero vs non-zero signal.
       One overlaid panel per method (KDE-smoothed density contours optional).

  B1b. signal_x_read_support_3prime.png
       Same, but y = dRNA signal at TTS.

  B1c. signal_x_read_support_zero_pct.png
       Bar chart: fraction of isoforms with zero signal at TSS and TTS
       per read-support bin, per method — shows whether zero-signal isoforms
       are concentrated in low-support bins (noise) or spread across bins (real).

Usage:
    python signal_read_support.py \\
        --bed label1:bed1.bed label2:bed2.bed ... \\
        --read-map label1:map1.txt label2:map2.txt ... \\
        --cage-plus cage_plus.bg --cage-minus cage_minus.bg \\
        --qs-plus qs_plus.bg   --qs-minus qs_minus.bg \\
        --output output_dir/
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
    from signal_utils import parse_isoforms, load_signal_tracks, isoform_signal
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
    from evaluation.signal_utils import parse_isoforms, load_signal_tracks, isoform_signal

apply_rc()

# ── Constants ────────────────────────────────────────────────────────────────

MAX_ISO    = 50_000   # cap for scatter (rasterized)
RNG_SEED   = 42
SUPPORT_BINS   = [1, 2, 5, 10, 20, float("inf")]
SUPPORT_LABELS = ["1", "2–4", "5–9", "10–19", "20+"]

# Colours for zero vs non-zero signal
COL_ZERO    = "#CC79A7"   # reddish purple — zero signal
COL_NONZERO = "#0072B2"   # blue — has signal


# ── Read-map loading ─────────────────────────────────────────────────────────

def load_read_map(path: str | Path) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    n_self = n_total = 0
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            iso_id  = parts[0]
            reads   = parts[1].split(",")
            n_reads = len(reads)
            counts[iso_id] = n_reads
            n_total += 1
            if n_reads == 1 and reads[0] == iso_id:
                n_self += 1
    if n_total > 0 and n_self / n_total > 0.9:
        return {}
    return counts


_ENSG_RE = re.compile(r"_ENSG\d")


def _lookup_count(name: str, rc: Dict[str, int]) -> int:
    if name in rc:
        return rc[name]
    m = _ENSG_RE.search(name)
    if m:
        tid = name[:m.start()]
        if tid in rc:
            return rc[tid]
    idx = name.rfind("_")
    if idx > 0 and name[:idx] in rc:
        return rc[name[:idx]]
    return 0


# ── Data building ─────────────────────────────────────────────────────────────

def build_data(
    beds_by_method: Dict[str, List[dict]],
    read_maps:       Dict[str, Dict[str, int]],
    cage_p, cage_m, qs_p, qs_m,
) -> Dict[str, List[dict]]:
    """Return {method: [{reads, tss_sig, tts_sig}]}."""
    data: Dict[str, List[dict]] = {}
    rng = np.random.default_rng(RNG_SEED)
    for method, isos in beds_by_method.items():
        rc = read_maps.get(method, {})
        rows = []
        for iso in isos:
            reads = _lookup_count(iso["name"], rc) if rc else 0
            if reads == 0:
                continue
            tss_sig, tts_sig = isoform_signal(iso, cage_p, cage_m, qs_p, qs_m)
            rows.append({"reads": reads, "tss_sig": tss_sig, "tts_sig": tts_sig})
        # subsample if very large
        if len(rows) > MAX_ISO:
            idx = rng.choice(len(rows), MAX_ISO, replace=False)
            rows = [rows[i] for i in idx]
        data[method] = rows
    return data


def support_bin(n: int) -> int:
    for i, upper in enumerate(SUPPORT_BINS[1:]):
        if n < upper:
            return i
    return len(SUPPORT_BINS) - 2


# ── Plot B1a/b: scatter panels ────────────────────────────────────────────────

def _scatter_panel(
    data: Dict[str, List[dict]],
    sig_key: str,
    sig_label: str,
    output_path: Path,
):
    methods = [m for m in data if data[m]]
    if not methods:
        return

    n     = len(methods)
    ncols = min(n, 4)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(W2, W1 * 0.82 * nrows / max(nrows, 1)),
        squeeze=False,
    )

    for idx, m in enumerate(methods):
        ax   = axes[idx // ncols][idx % ncols]
        rows = data[m]
        xs   = np.array([r["reads"]   for r in rows], dtype=float)
        ys   = np.array([r[sig_key]   for r in rows], dtype=float)

        zero    = ys == 0
        nonzero = ~zero

        # Non-zero signal first (bottom layer), zero on top
        if nonzero.any():
            ax.scatter(xs[nonzero] + 0.5, ys[nonzero] + 1e-4,
                       s=1.5, alpha=0.3, color=COL_NONZERO,
                       edgecolors="none", rasterized=True)
        if zero.any():
            # jitter y slightly for visibility
            jitter = np.random.default_rng(42).uniform(-0.3, 0.3, zero.sum())
            ax.scatter(xs[zero] + 0.5, np.full(zero.sum(), 1e-4) + 10 ** jitter * 1e-5,
                       s=1.0, alpha=0.2, color=COL_ZERO,
                       edgecolors="none", rasterized=True)

        ax.set_xscale("log")
        ax.set_yscale("log")
        style_ax(ax)
        pct_zero = zero.sum() / len(rows) * 100 if rows else 0
        ax.text(0.04, 0.97, m, transform=ax.transAxes,
                ha="left", va="top", fontsize=5, fontweight="bold")
        ax.text(0.04, 0.88, f"{pct_zero:.0f}% zero-signal",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=5, color=COL_ZERO)

        if idx >= n - ncols:
            ax.set_xlabel("Read support", fontsize=6)
        if idx % ncols == 0:
            ax.set_ylabel(sig_label, fontsize=6)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    # Shared legend
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_NONZERO,
               markersize=5, label="Signal > 0"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_ZERO,
               markersize=5, label="Zero signal"),
    ]
    fig.legend(handles=handles, loc="lower right", fontsize=6,
               frameon=False, bbox_to_anchor=(1.0, 0.0))
    fig.tight_layout(pad=0.3)
    savefig(fig, output_path)


# ── Plot B1c: zero-signal fraction by read-support bin ───────────────────────

def plot_zero_by_support_bin(
    data: Dict[str, List[dict]],
    output_dir: Path,
):
    methods = list(data.keys())
    if not methods:
        return
    styler = ModeStyler(methods)
    n_bins = len(SUPPORT_LABELS)

    for sig_key, sig_label, fname in [
        ("tss_sig", "TSS (CAGE)",  "signal_x_read_support_zero_pct_5prime.png"),
        ("tts_sig", "TTS (dRNA)",  "signal_x_read_support_zero_pct_3prime.png"),
    ]:
        x     = np.arange(n_bins)
        width = 0.8 / max(len(methods), 1)

        fig, ax = plt.subplots(figsize=(W2, W1 * 0.65))
        handles = []
        for i, m in enumerate(methods):
            rows  = data[m]
            fracs = []
            for b in range(n_bins):
                bin_rows = [r for r in rows if support_bin(r["reads"]) == b]
                if bin_rows:
                    fracs.append(sum(r[sig_key] == 0 for r in bin_rows) / len(bin_rows))
                else:
                    fracs.append(0.0)
            offset = (i - len(methods) / 2 + 0.5) * width
            ax.bar(x + offset, fracs, width * 0.9,
                   color=styler.color(m), edgecolor="none", alpha=0.85)
            handles.append(styler.legend_handle(m, label=m, markersize=5))

        ax.set_xticks(x)
        ax.set_xticklabels(SUPPORT_LABELS, fontsize=7)
        ax.set_ylim(0, 1)
        style_ax(ax,
                 xlabel="Read support (reads per isoform)",
                 ylabel=f"Fraction with zero {sig_label} signal")
        legend_outside(fig, handles=handles, loc="upper right",
                       bbox_to_anchor=(1.0, 1.0), ncol=1, fontsize=6)
        fig.tight_layout(pad=0.3)
        savefig(fig, output_dir / fname)


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bed",        nargs="+", required=True,
                        help="label:path pairs for BED12 isoform files")
    parser.add_argument("--read-map",   nargs="+", default=[],
                        help="label:path pairs for isoform read-map files")
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

    read_maps: Dict[str, Dict[str, int]] = {}
    for entry in args.read_map:
        if ":" not in entry:
            continue
        label, path = entry.split(":", 1)
        if not Path(path).exists():
            continue
        rc = load_read_map(path)
        if rc:
            read_maps[label] = rc

    if args.verbose:
        print("  Loading signal tracks...", file=sys.stderr)
    cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
        args.cage_plus, args.cage_minus, args.qs_plus, args.qs_minus,
    )

    if args.verbose:
        print("  Computing per-isoform signal...", file=sys.stderr)
    data = build_data(beds_by_method, read_maps, cage_p, cage_m, qs_p, qs_m)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    _scatter_panel(data, "tss_sig", "CAGE signal at TSS",
                   output_dir / "signal_x_read_support_5prime.png")
    _scatter_panel(data, "tts_sig", "dRNA signal at TTS",
                   output_dir / "signal_x_read_support_3prime.png")
    plot_zero_by_support_bin(data, output_dir)

    print(f"Saved signal × read support plots to {args.output}")


if __name__ == "__main__":
    main()
