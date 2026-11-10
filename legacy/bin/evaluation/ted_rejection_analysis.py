#!/usr/bin/env python3
"""
ted_rejection_analysis.py — Two views on why TED clusters are rejected:

  1. Rejection reason breakdown (per config)
     Stacked bar showing what fraction of rejected candidates failed on:
       - TSS only     (tss_below_threshold)
       - TTS only     (tts_below_threshold)
       - Both ends    (tss_below_threshold+tts_below_threshold)
       - Snap filter  (dist=Xbp proximity collapse)
       - HDBSCAN noise

  2. Spliced-length stratification of end hit rates
     Bins candidates by spliced RNA length (sum of exon block sizes from
     the firstpass BED12) and plots how joint TP rate (both ends hit
     orthogonal peaks), 5' TSS hit rate, and 3' TTS hit rate change
     across length bins — one line per config.

Usage
-----
  python ted_rejection_analysis.py \\
      --ted-log label:ted_log.tsv:firstpass.bed [label:...] \\
      --cage-peaks cage.bed \\
      --drna-peaks drna.bed \\
      --output outdir/ \\
      [--chrom chr22] \\
      [--window 50] \\
      [--n-bins 10]
"""

from __future__ import annotations

import argparse
import logging
import sys
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s")
log = logging.getLogger(__name__)

WINDOW = 50

try:
    from matplotlib import rcParams
    rcParams.update({
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
        "xtick.labelsize": 6, "ytick.labelsize": 6,
        "legend.fontsize": 6, "figure.dpi": 150,
        "axes.spines.top": False, "axes.spines.right": False,
    })
except Exception:
    pass

PALETTE = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
           "#0072B2", "#D55E00", "#CC79A7", "#999999"]

REJECT_COLORS = {
    "pass":     "#009E73",
    "tss_only": "#0072B2",
    "tts_only": "#E69F00",
    "both":     "#D55E00",
    "snap":     "#CC79A7",
    "noise":    "#bbbbbb",
    "other":    "#555555",
}


# ── Peak helpers ─────────────────────────────────────────────────────────────

def _parse_peaks(path: str, chrom: str = None
                 ) -> Dict[Tuple[str, str], List[Tuple[int, int]]]:
    peaks: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            c, start, end, strand = cols[0], int(cols[1]), int(cols[2]), cols[5]
            if chrom and c != chrom:
                continue
            peaks[(c, strand)].append((start, end))
    for k in peaks:
        peaks[k].sort()
    return peaks


def _hit_peak(pos: int, intervals: List[Tuple[int, int]]) -> bool:
    if not intervals:
        return False
    starts = [iv[0] for iv in intervals]
    idx = bisect_left(starts, pos)
    for i in (idx - 1, idx):
        if 0 <= i < len(intervals):
            s, e = intervals[i]
            d = 0 if s <= pos < e else min(abs(pos - s), abs(pos - (e - 1)))
            if d <= WINDOW:
                return True
    return False


# ── Firstpass BED12 ──────────────────────────────────────────────────────────

def load_firstpass_bed(path: str, chrom: str = None) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, dtype=str,
                     names=list(range(22)))
    # BED12 cols: 0=chrom, 1=start (genomic), 2=end (genomic), 5=strand, 10=block_sizes
    df = df.rename(columns={0: "chrom", 1: "bed_start", 2: "bed_end",
                             5: "strand", 10: "block_sizes"})
    df["bed_start"] = pd.to_numeric(df["bed_start"], errors="coerce")
    df["bed_end"]   = pd.to_numeric(df["bed_end"],   errors="coerce")
    if chrom:
        df = df[df["chrom"] == chrom]

    def _spliced(s):
        try:
            return sum(int(x) for x in str(s).rstrip(",").split(",") if x)
        except Exception:
            return np.nan

    df["spliced_length"] = df["block_sizes"].apply(_spliced)
    df = df[["chrom", "bed_start", "bed_end", "strand", "spliced_length"]].dropna()
    df["bed_start"] = df["bed_start"].astype(np.int64)
    df["bed_end"]   = df["bed_end"].astype(np.int64)

    # tss_pos = 5' end, tts_pos = 3' end — strand-aware
    df["tss_pos"] = np.where(df["strand"] == "+", df["bed_start"], df["bed_end"])
    df["tts_pos"] = np.where(df["strand"] == "+", df["bed_end"],   df["bed_start"])

    return df[["chrom", "tss_pos", "tts_pos", "strand", "spliced_length"]].reset_index(drop=True)


# ── TED log loading ──────────────────────────────────────────────────────────

def _chrom(jid): return str(jid).split(":")[0]
def _strand(jid):
    p = str(jid).split(":")
    return p[2] if len(p) > 2 else "+"


def load_ted_log(path: str, chrom: str = None) -> pd.DataFrame:
    """Load a flair_ted.ted_log.tsv, dropping --ted_global partition-summary
    rows and adding `chrom`/`strand` columns derived from junc_id.

    Delegates the schema-aware filtering to ted_log_loader.load_ted_log so
    every diagnostic that consumes ted_log.tsv applies the same global-mode
    filters. Without that filter, partition-level `global_peak` rows
    (junc_id like `chr*:None-None:strand:0j`, with one of tss_pos/tts_pos
    set to the -1 sentinel) leak through as phantom per-isoform rejections
    and corrupt the rejection-reason and length-stratified plots downstream.
    """
    # Sibling import — the shim has no relative imports so it works whether
    # this script is invoked directly (`python ted_rejection_analysis.py`)
    # or imported as part of the `evaluation` package.
    from ted_log_loader import load_ted_log as _shared_load_ted_log

    df = _shared_load_ted_log(path)
    df["chrom"]  = df["junc_id"].apply(_chrom)
    df["strand"] = df["junc_id"].apply(_strand)
    if chrom:
        df = df[df["chrom"] == chrom]
    return df.reset_index(drop=True)


# ── TP labelling ─────────────────────────────────────────────────────────────

def assign_hits(df: pd.DataFrame, cage_peaks, drna_peaks) -> pd.DataFrame:
    """Add tss_hit, tts_hit, joint_tp columns to pass/reject rows."""
    mask = df["status"].isin(["pass", "reject"]) & df["tss_pos"].notna()
    sub  = df[mask].copy()
    tss_hit = np.zeros(len(sub), dtype=np.int8)
    tts_hit = np.zeros(len(sub), dtype=np.int8)
    chroms  = sub["chrom"].values
    strands = sub["strand"].values
    tss_pos = sub["tss_pos"].values.astype(np.int64)
    tts_pos = sub["tts_pos"].values.astype(np.int64)
    for i in range(len(sub)):
        key = (chroms[i], strands[i])
        tss_hit[i] = int(_hit_peak(int(tss_pos[i]), cage_peaks.get(key, [])))
        tts_hit[i] = int(_hit_peak(int(tts_pos[i]), drna_peaks.get(key, [])))
    sub["tss_hit"]  = tss_hit
    sub["tts_hit"]  = tts_hit
    sub["joint_tp"] = (tss_hit & tts_hit).astype(np.int8)
    return sub


# ── Rejection reason classification ─────────────────────────────────────────

def classify_drop(reason) -> str:
    if pd.isna(reason) or str(reason) in ("nan", ""):
        return "pass"
    r = str(reason).lower()
    if r.startswith("hdbscan"):
        return "noise"
    if r.startswith("dist="):
        return "snap"
    if "tss" in r and "tts" in r:
        return "both"
    if "tss" in r:
        return "tss_only"
    if "tts" in r:
        return "tts_only"
    return "other"


# ── Plot 1: rejection reason breakdown ───────────────────────────────────────

def plot_rejection_reasons(datasets: Dict[str, pd.DataFrame], outdir: Path):
    categories = ["pass", "tss_only", "tts_only", "both", "snap", "noise", "other"]
    labels = list(datasets.keys())
    counts = {c: [] for c in categories}

    for label, df in datasets.items():
        df["_class"] = df["drop_reason"].apply(classify_drop)
        vc = df["_class"].value_counts()
        for c in categories:
            counts[c].append(vc.get(c, 0))

    x = np.arange(len(labels))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(max(7, len(labels) * 1.4), 3.5))

    # Absolute counts
    bottoms = np.zeros(len(labels))
    handles = []
    for cat in categories:
        vals = np.array(counts[cat], dtype=float)
        ax1.bar(x, vals, bottom=bottoms, color=REJECT_COLORS[cat],
                width=0.6, edgecolor="white", linewidth=0.3,
                label=cat.replace("_", " "))
        if vals.sum() > 0:
            handles.append(mpatches.Patch(facecolor=REJECT_COLORS[cat],
                                          label=cat.replace("_", " ")))
        bottoms += vals

    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=40, ha="right", fontsize=6)
    ax1.set_ylabel("Candidate count")
    ax1.set_title("Rejection reason (counts)")
    ax1.legend(handles=handles, frameon=False, fontsize=6)

    # % of rejected only
    reject_cats = ["tss_only", "tts_only", "both", "snap", "noise", "other"]
    totals = np.array([sum(counts[c][i] for c in reject_cats)
                       for i in range(len(labels))], dtype=float)
    totals[totals == 0] = 1
    bottoms2 = np.zeros(len(labels))
    for cat in reject_cats:
        vals = np.array(counts[cat], dtype=float) / totals * 100
        ax2.bar(x, vals, bottom=bottoms2, color=REJECT_COLORS[cat],
                width=0.6, edgecolor="white", linewidth=0.3)
        bottoms2 += vals

    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, rotation=40, ha="right", fontsize=6)
    ax2.set_ylabel("% of rejected candidates")
    ax2.set_title("Rejection reason (% of rejects)")
    ax2.set_ylim(0, 105)

    fig.tight_layout()
    fig.savefig(outdir / "rejection_reasons.png", bbox_inches="tight", dpi=150)
    plt.close(fig)

    # TSV
    rows = []
    for i, label in enumerate(labels):
        total = sum(counts[c][i] for c in categories)
        row = {"config": label}
        for c in categories:
            row[c] = counts[c][i]
            row[f"{c}_pct"] = f"{counts[c][i]/total*100:.1f}" if total else "0.0"
        rows.append(row)
    pd.DataFrame(rows).to_csv(outdir / "rejection_reasons.tsv", sep="\t", index=False)
    log.info("Saved rejection_reasons.png + .tsv")


# ── Plot 2: end hit rate vs spliced length ────────────────────────────────────

def plot_length_stratification(
    datasets: Dict[str, Tuple[pd.DataFrame, pd.DataFrame]],
    outdir: Path,
    n_bins: int = 10,
):
    """
    datasets: {label: (ted_log_df_with_hits, firstpass_bed_df)}
    Uses genomic span (|tts_pos - tss_pos|) from the TED log as a length proxy
    for binning — the exact-coordinate merge with the BED rarely succeeds because
    TED log tss_pos/tts_pos are cluster centroids, not isoform endpoints.
    """
    merged_dfs = {}
    all_lengths = []
    for label, (ted_df, _bed_df) in datasets.items():
        m = ted_df.copy()
        m["spliced_length"] = (m["tts_pos"] - m["tss_pos"]).abs()
        m = m[m["spliced_length"].notna() & (m["spliced_length"] > 0)].copy()
        merged_dfs[label] = m
        all_lengths.extend(m["spliced_length"].values.tolist())

    if not all_lengths:
        log.warning("No spliced length data available — skipping length plot")
        return

    bin_edges = np.unique(np.percentile(all_lengths, np.linspace(0, 100, n_bins + 1)))
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    metrics = [
        ("tss_hit",  "5′ TSS hit rate"),
        ("tts_hit",  "3′ TTS hit rate"),
        ("joint_tp", "Paired TP rate (both ends)"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(14, 3.5))

    for ax, (metric, ylabel) in zip(axes, metrics):
        for ci, (label, merged) in enumerate(merged_dfs.items()):
            color = PALETTE[ci % len(PALETTE)]
            rates, counts_, mids = [], [], []
            for j in range(len(bin_edges) - 1):
                mask = ((merged["spliced_length"] >= bin_edges[j]) &
                        (merged["spliced_length"] <  bin_edges[j + 1]))
                sub = merged[mask]
                if len(sub) < 5:
                    continue
                rates.append(sub[metric].mean())
                counts_.append(len(sub))
                mids.append(bin_centers[j])

            ax.plot(mids, rates, marker="o", markersize=4,
                    linewidth=1.2, color=color, label=label)
            sizes = [max(10, 5 * np.log1p(c)) for c in counts_]
            ax.scatter(mids, rates, s=sizes, color=color, zorder=3, alpha=0.6)

        ax.set_xlabel("Genomic span (nt)")
        ax.set_ylabel(ylabel)
        ax.set_ylim(0, 1.05)
        ax.set_title(ylabel)
        ax.legend(frameon=False, fontsize=6)

    fig.tight_layout()
    fig.savefig(outdir / "end_hit_rate_vs_spliced_length.png",
                bbox_inches="tight", dpi=150)
    plt.close(fig)
    log.info("Saved end_hit_rate_vs_spliced_length.png")


# ── Plot 3: rejection reason stratified by spliced length ────────────────────

def plot_rejection_by_length(
    datasets: Dict[str, Tuple[pd.DataFrame, pd.DataFrame]],
    outdir: Path,
    n_bins: int = 8,
):
    """Stacked bar showing rejection reason breakdown across genomic-span bins."""
    reject_cats = ["tss_only", "tts_only", "both", "snap", "noise"]
    labels = list(datasets.keys())

    # Use genomic span as length proxy (same as plot_length_stratification)
    all_lengths = []
    merged_dfs = {}
    for label, (ted_df, _bed_df) in datasets.items():
        m = ted_df.copy()
        m["spliced_length"] = (m["tts_pos"] - m["tss_pos"]).abs()
        m = m[m["spliced_length"].notna() & (m["spliced_length"] > 0)].copy()
        m["_class"] = m["drop_reason"].apply(classify_drop)
        merged_dfs[label] = m
        all_lengths.extend(m["spliced_length"].values.tolist())

    if not all_lengths:
        return

    bin_edges = np.unique(np.percentile(all_lengths, np.linspace(0, 100, n_bins + 1)))
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_xlabels = [f"{int(e/1000)}k" if e >= 1000 else str(int(e))
                   for e in bin_centers]

    n_labels = len(labels)
    fig, axes = plt.subplots(1, n_labels,
                             figsize=(max(5, n_labels * 3.5), 3.5),
                             squeeze=False)

    for col, label in enumerate(labels):
        ax = axes[0][col]
        merged = merged_dfs[label]
        x = np.arange(len(bin_centers))

        cat_pcts = {c: [] for c in reject_cats}
        for j in range(len(bin_edges) - 1):
            mask = ((merged["spliced_length"] >= bin_edges[j]) &
                    (merged["spliced_length"] <  bin_edges[j + 1]))
            sub = merged[mask]
            reject_sub = sub[sub["_class"] != "pass"]
            total = max(len(reject_sub), 1)
            for c in reject_cats:
                cat_pcts[c].append((reject_sub["_class"] == c).sum() / total * 100)

        bottoms = np.zeros(len(bin_centers))
        for cat in reject_cats:
            vals = np.array(cat_pcts[cat])
            ax.bar(x, vals, bottom=bottoms, color=REJECT_COLORS[cat],
                   width=0.7, edgecolor="white", linewidth=0.3,
                   label=cat.replace("_", " "))
            bottoms += vals

        ax.set_xticks(x)
        ax.set_xticklabels(bin_xlabels, rotation=40, ha="right", fontsize=5)
        ax.set_xlabel("Genomic span (nt)")
        ax.set_ylim(0, 105)
        ax.set_title(label, fontsize=6)
        if col == 0:
            ax.set_ylabel("% of rejected candidates")
            handles = [mpatches.Patch(facecolor=REJECT_COLORS[c],
                                      label=c.replace("_", " "))
                       for c in reject_cats]
            ax.legend(handles=handles, frameon=False, fontsize=5,
                      loc="upper right")

    fig.suptitle("Rejection reason by spliced RNA length", fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "rejection_reason_by_length.png",
                bbox_inches="tight", dpi=150)
    plt.close(fig)
    log.info("Saved rejection_reason_by_length.png")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--ted-log", nargs="+", required=True,
                        help="label:ted_log.tsv:firstpass.bed triplets")
    parser.add_argument("--cage-peaks", required=True)
    parser.add_argument("--drna-peaks", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--chrom", default=None)
    parser.add_argument("--window", type=int, default=50)
    parser.add_argument("--n-bins", type=int, default=10)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    global WINDOW
    WINDOW = args.window

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    cage_peaks = _parse_peaks(args.cage_peaks, chrom=args.chrom)
    drna_peaks = _parse_peaks(args.drna_peaks, chrom=args.chrom)
    log.info("%d CAGE peaks, %d dRNA peaks",
             sum(len(v) for v in cage_peaks.values()),
             sum(len(v) for v in drna_peaks.values()))

    raw_logs: Dict[str, pd.DataFrame] = {}      # full log (all statuses) for rejection plot
    length_datasets: Dict[str, Tuple] = {}       # (ted_with_hits, bed) for length plots

    for entry in args.ted_log:
        parts = entry.split(":")
        if len(parts) != 3:
            log.error("Expected label:ted_log:firstpass_bed, got: %s", entry)
            sys.exit(1)
        label, log_path, bed_path = parts

        if not Path(log_path).exists():
            log.warning("TED log not found: %s", log_path); continue
        if not Path(bed_path).exists():
            log.warning("Firstpass BED not found: %s", bed_path); continue

        full_df = load_ted_log(log_path, chrom=args.chrom)
        bed_df  = load_firstpass_bed(bed_path, chrom=args.chrom)

        # Label hits on pass/reject subset
        hits_df = assign_hits(full_df, cage_peaks, drna_peaks)

        raw_logs[label] = full_df
        length_datasets[label] = (hits_df, bed_df)

        log.info("%s: %d total rows, %d bed entries (pass=%d reject=%d)",
                 label, len(full_df), len(bed_df),
                 (full_df["status"] == "pass").sum(),
                 (full_df["status"] == "reject").sum())

    if raw_logs:
        plot_rejection_reasons(raw_logs, outdir)

    if length_datasets:
        plot_length_stratification(length_datasets, outdir, n_bins=args.n_bins)
        plot_rejection_by_length(length_datasets, outdir, n_bins=args.n_bins)

    log.info("Done — output in %s", outdir)


if __name__ == "__main__":
    main()
