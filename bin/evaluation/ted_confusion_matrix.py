#!/usr/bin/env python3
"""
ted_confusion_matrix.py — Per-TED-config confusion matrices from the ted_log.

Reads one or more TED logs (with their --threshold_tss / --threshold_tts already
recorded per row) and the orthogonal CAGE / dRNA peak BEDs.  For every TED log:

  * Labels each cluster as a joint TP iff its TSS is within ±WINDOW of a CAGE
    peak AND its TTS is within ±WINDOW of a dRNA peak.
  * Crosses with the cluster's pass/reject status to build a 2x2 matrix of
    {pass, reject} x {joint TP, joint FP}.
  * Writes one PNG per TED config plus a multi-config comparison figure and a
    summary TSV.

This intentionally does NOT do JC-level deduplication — it's a per-cluster
diagnostic.  For dedup'd precision/recall use the existing precision_recall
TSVs from TedEndPrecision.

Usage
-----
    python ted_confusion_matrix.py \
        --ted-log label:path [label:path ...] \
        --cage-peaks cage.bed \
        --drna-peaks drna.bed \
        --output outdir/ \
        [--window 50] [--verbose]
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


# ── peak helpers ───────────────────────────────────────────────────────────

def _parse_peaks(path: str) -> Dict[Tuple[str, str], List[Tuple[int, int]]]:
    peaks: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            c, start, end, strand = cols[0], int(cols[1]), int(cols[2]), cols[5]
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


# ── TED log loading ────────────────────────────────────────────────────────

def _parse_chrom(junc_id: str) -> str:
    return str(junc_id).split(":")[0]


def _parse_strand(junc_id: str) -> str:
    parts = str(junc_id).split(":")
    return parts[2] if len(parts) > 2 else "+"


def load_ted_log(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    for col in ("tss_pos", "tts_pos", "n_reads"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[df["status"].isin(["pass", "reject"])].copy()
    df = df[df["tss_pos"].notna() & (df["tss_pos"] >= 0)].copy()
    df["chrom"] = df["junc_id"].apply(_parse_chrom)
    df["strand"] = df["junc_id"].apply(_parse_strand)
    return df.reset_index(drop=True)


def label_joint_tp(df: pd.DataFrame, cage: dict, drna: dict) -> pd.DataFrame:
    df = df.copy()
    tss = df["tss_pos"].values.astype(np.int64)
    tts = df["tts_pos"].values.astype(np.int64)
    chroms = df["chrom"].values
    strands = df["strand"].values
    tss_hit = np.zeros(len(df), dtype=np.int8)
    tts_hit = np.zeros(len(df), dtype=np.int8)
    for i in range(len(df)):
        key = (chroms[i], strands[i])
        tss_hit[i] = int(_hit_peak(int(tss[i]), cage.get(key, [])))
        tts_hit[i] = int(_hit_peak(int(tts[i]), drna.get(key, [])))
    df["tss_hit"] = tss_hit
    df["tts_hit"] = tts_hit
    df["joint_tp"] = (tss_hit & tts_hit).astype(np.int8)
    return df


def confusion_counts(df: pd.DataFrame) -> Dict[str, int]:
    """Return {tp, fp, fn, tn} counts.

    pass × joint_tp=1 → TP        (kept a real both-end isoform)
    pass × joint_tp=0 → FP        (kept a not-both-end isoform)
    reject × joint_tp=1 → FN      (rejected a real both-end isoform)
    reject × joint_tp=0 → TN      (correctly rejected a not-both-end isoform)
    """
    pass_mask = df["status"] == "pass"
    tp = int(((df["joint_tp"] == 1) & pass_mask).sum())
    fp = int(((df["joint_tp"] == 0) & pass_mask).sum())
    fn = int(((df["joint_tp"] == 1) & ~pass_mask).sum())
    tn = int(((df["joint_tp"] == 0) & ~pass_mask).sum())
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn}


def _plot_one(ax, counts: dict, label: str):
    tp, fp, fn, tn = counts["tp"], counts["fp"], counts["fn"], counts["tn"]
    mat = np.array([[tp, fp], [fn, tn]], dtype=float)
    im = ax.imshow(mat, cmap="Blues", aspect="auto")
    for r in range(2):
        for c in range(2):
            v = int(mat[r, c])
            ax.text(c, r, f"{v:,}", ha="center", va="center",
                    fontsize=8,
                    color="white" if mat[r, c] > mat.max() * 0.5 else "black")
    ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
    ax.set_xticklabels(["Joint TP", "Joint FP"])
    ax.set_yticklabels(["Pass", "Reject"])
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    rec  = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1   = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
    fpr  = fp / (fp + tn) if (fp + tn) > 0 else 0.0
    ax.set_title(
        f"{label}\nP={prec:.3f}  R={rec:.3f}  F1={f1:.3f}  FPR={fpr:.3f}",
        fontsize=7,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ted-log", nargs="+", required=True,
                   help="label:path TED log entries")
    p.add_argument("--cage-peaks", required=True)
    p.add_argument("--drna-peaks", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--window", type=int, default=50)
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args()

    global WINDOW
    WINDOW = args.window
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("Loading peaks")
    cage = _parse_peaks(args.cage_peaks)
    drna = _parse_peaks(args.drna_peaks)

    rows = []
    config_data: List[Tuple[str, Dict[str, int]]] = []
    for entry in args.ted_log:
        if ":" not in entry:
            log.warning("Bad entry %s — expected label:path", entry); continue
        label, path = entry.split(":", 1)
        if not Path(path).exists() or Path(path).stat().st_size == 0:
            log.warning("Skipping missing/empty %s", path); continue

        log.info("== %s ==", label)
        df = load_ted_log(path)
        if df.empty:
            log.warning("  empty after filtering, skipping"); continue
        df = label_joint_tp(df, cage, drna)
        c = confusion_counts(df)

        log.info("  pass=%d  reject=%d  TP=%d  FP=%d  FN=%d  TN=%d",
                 (df["status"] == "pass").sum(),
                 (df["status"] == "reject").sum(),
                 c["tp"], c["fp"], c["fn"], c["tn"])

        # Per-config single-panel
        fig, ax = plt.subplots(figsize=(3.6, 3.4))
        _plot_one(ax, c, label)
        fig.tight_layout()
        fig.savefig(outdir / f"confusion_{label}.png", bbox_inches="tight", dpi=150)
        plt.close(fig)
        log.info("  → confusion_%s.png", label)

        config_data.append((label, c))
        prec = c["tp"] / (c["tp"] + c["fp"]) if (c["tp"] + c["fp"]) else 0.0
        rec  = c["tp"] / (c["tp"] + c["fn"]) if (c["tp"] + c["fn"]) else 0.0
        f1   = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
        rows.append(dict(
            config=label, n_clusters=len(df),
            n_pass=int((df["status"] == "pass").sum()),
            n_reject=int((df["status"] == "reject").sum()),
            **c,
            precision=round(prec, 4), recall=round(rec, 4), f1=round(f1, 4),
        ))

    if not config_data:
        log.error("No usable TED logs"); sys.exit(1)

    # Multi-config comparison panel
    n = len(config_data)
    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(3.4 * ncols, 3.2 * nrows),
                             squeeze=False)
    for i, (label, c) in enumerate(config_data):
        ax = axes[i // ncols][i % ncols]
        _plot_one(ax, c, label)
    for i in range(n, nrows * ncols):
        axes[i // ncols][i % ncols].set_visible(False)
    fig.suptitle("TED confusion matrices: rows = pass/reject, cols = joint TP / joint FP",
                 fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "confusion_all_configs.png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    log.info("→ confusion_all_configs.png")

    # Summary TSV
    pd.DataFrame(rows).to_csv(outdir / "confusion_summary.tsv",
                              sep="\t", index=False)
    log.info("→ confusion_summary.tsv")
    log.info("Done — outputs in %s", outdir)


if __name__ == "__main__":
    main()
