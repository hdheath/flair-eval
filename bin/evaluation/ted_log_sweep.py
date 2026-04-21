#!/usr/bin/env python3
"""
ted_log_sweep.py — Threshold / weight sweep over the full TED candidate pool.

The standard TED component diagnostic optimises weights on POST-filter isoforms
(only what TED accepted into the BED file).  That is biased: the rejected
candidates never appear, so the sweep F1 looks much higher than the P/R the
evaluation pipeline actually reports.

This script uses the TED log (--ted_log), which records EVERY candidate cluster
(status = pass OR reject) with raw component scores.  For each point in the
weight / threshold grid it:

  1. Recomputes TSS and TTS reality scores from raw components.
  2. Applies the threshold pair (t_tss, t_tts) — both ends must pass.
  3. Groups survivors by junction chain (junc_id field).
  4. Deduplicates ends within each junction-chain group (same logic as
     compute_jc_deduplicated_precision_recall in ted_end_precision.py).
  5. Evaluates orthogonal-peak precision + peak-coverage recall → F1.

This produces a sweep F1 that matches what the pipeline reports, so the
optimal (weights, thresholds) found here can be used directly as TED flags.

Outputs
-------
  weight_sweep_5prime.png      — ternary weight surface for TSS
  weight_sweep_3prime.png      — ternary weight surface for TTS
  joint_threshold_sweep.png    — 2-D (t_tss × t_tts) F1 heatmap
  summary_stats.tsv            — optimal per-end and joint settings

Usage
-----
  python ted_log_sweep.py \\
      --ted-log label1:log1.tsv [label2:log2.tsv ...] \\
      --cage-peaks cage.bed \\
      --drna-peaks drna.bed \\
      --output outdir/
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, PALETTE
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, PALETTE

apply_rc()
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s")
log = logging.getLogger(__name__)

WINDOW = 50   # bp match window (same as evaluation pipeline)


# ── Peak index ──────────────────────────────────────────────────────────────

def _parse_peaks(path: str, chrom: str = None) -> Dict[Tuple[str, str], List[Tuple[int, int]]]:
    """BED6 → {(chrom, strand): sorted [(start, end)]}.  Optionally filter to one chrom."""
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


def _hit_peak(pos: int, intervals: List[Tuple[int, int]]) -> Optional[Tuple[int, int]]:
    """Return nearest peak interval within WINDOW bp of pos, or None."""
    if not intervals:
        return None
    starts = [iv[0] for iv in intervals]
    idx = bisect_left(starts, pos)
    best_iv, best_d = None, WINDOW + 1
    for i in (idx - 1, idx):
        if 0 <= i < len(intervals):
            s, e = intervals[i]
            d = 0 if s <= pos < e else min(abs(pos - s), abs(pos - (e - 1)))
            if d < best_d:
                best_d, best_iv = d, intervals[i]
    return best_iv if best_d <= WINDOW else None


# ── TED log loading ─────────────────────────────────────────────────────────

def _parse_chrom(junc_id: str) -> str:
    return str(junc_id).split(":")[0]


def _parse_strand(junc_id: str) -> str:
    parts = str(junc_id).split(":")
    return parts[2] if len(parts) > 2 else "+"


def _parse_njunc(junc_id: str) -> int:
    """Extract junction count from junc_id suffix like '6j'."""
    parts = str(junc_id).split(":")
    if parts:
        suf = parts[-1]
        if suf.endswith("j"):
            try:
                return int(suf[:-1])
            except ValueError:
                pass
    return 0


def load_ted_log(path: str, chrom: str = None) -> pd.DataFrame:
    """Load TED log, keep only scoreable pass/reject rows.

    If chrom is given, restricts to rows where junc_id starts with that chromosome.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    num_cols = [
        "tss_pos", "tts_pos", "n_reads",
        "TED_depth", "TED_tss_model", "TED_tts_model",
        "TED_tss_annot", "TED_tts_annot",
        "TED_tss_reality", "TED_tts_reality",
        "threshold_tss", "threshold_tts",
    ]
    for col in num_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df[df["status"].isin(["pass", "reject"])].copy()
    df = df[df["tss_pos"] >= 0].copy()
    df = df[df["TED_depth"].notna()].copy()

    if chrom:
        df = df[df["junc_id"].str.startswith(chrom + ":")].copy()

    return df


def load_all_logs(entries: List[str], chrom: str = None) -> Dict[str, pd.DataFrame]:
    result: Dict[str, pd.DataFrame] = {}
    for entry in entries:
        if ":" not in entry:
            continue
        label, path = entry.split(":", 1)
        if not Path(path).exists():
            log.warning("Not found: %s", path)
            continue
        df = load_ted_log(path, chrom=chrom)
        if not df.empty:
            result[label] = df
            log.info("Loaded %d candidates for %s (pass=%d, reject=%d)%s",
                     len(df), label,
                     (df["status"] == "pass").sum(),
                     (df["status"] == "reject").sum(),
                     f" [chrom={chrom}]" if chrom else "")
    return result


# ── Core metric replay ──────────────────────────────────────────────────────
# Strategy: precompute per-row peak hit IDs once (expensive), then all sweep
# iterations work on integer arrays with no Python loops over rows.

def _assign_peak_ids(
    df: pd.DataFrame,
    cage_peaks: Dict[Tuple[str, str], List[Tuple[int, int]]],
    drna_peaks: Dict[Tuple[str, str], List[Tuple[int, int]]],
) -> pd.DataFrame:
    """Add integer columns to df (done once per config, not per sweep point).

    tss_peak_id  — integer ID of matched CAGE peak interval, or -1
    tts_peak_id  — integer ID of matched dRNA peak interval, or -1
    cage_peak_uid — unique (chrom, strand, peak_idx) encoded as int for recall counting
    drna_peak_uid — same for dRNA

    Also adds: chrom, strand, junc_id (string) columns for grouping.
    """
    # Build flat peak lists with integer IDs per (chrom, strand)
    cage_flat: Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray, int]] = {}
    for key, ivs in cage_peaks.items():
        starts = np.array([iv[0] for iv in ivs], dtype=np.int64)
        ends   = np.array([iv[1] for iv in ivs], dtype=np.int64)
        cage_flat[key] = (starts, ends)

    drna_flat: Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray]] = {}
    for key, ivs in drna_peaks.items():
        starts = np.array([iv[0] for iv in ivs], dtype=np.int64)
        ends   = np.array([iv[1] for iv in ivs], dtype=np.int64)
        drna_flat[key] = (starts, ends)

    # Build a global peak → uid mapping so recall just counts unique uids
    cage_uid_base: Dict[Tuple[str, str], int] = {}
    base = 0
    for key in cage_peaks:
        cage_uid_base[key] = base
        base += len(cage_peaks[key])
    total_cage = base

    drna_uid_base: Dict[Tuple[str, str], int] = {}
    base = 0
    for key in drna_peaks:
        drna_uid_base[key] = base
        base += len(drna_peaks[key])
    total_drna = base

    def _lookup(pos_arr, chrom_arr, strand_arr, flat, uid_base):
        """Vectorised peak lookup. Returns (peak_id array, uid array), both int64, -1=no hit."""
        n = len(pos_arr)
        peak_ids = np.full(n, -1, dtype=np.int64)
        uids     = np.full(n, -1, dtype=np.int64)
        for key, (starts, ends) in flat.items():
            c, s = key
            mask = (chrom_arr == c) & (strand_arr == s)
            if not mask.any():
                continue
            idxs = np.where(mask)[0]
            pos  = pos_arr[idxs]
            ins  = np.searchsorted(starts, pos)
            for local_i, (gi, p, ins_i) in enumerate(zip(idxs, pos, ins)):
                best_pid, best_d = -1, WINDOW + 1
                for ci in (ins_i - 1, ins_i):
                    if 0 <= ci < len(starts):
                        d = 0 if starts[ci] <= p < ends[ci] else min(
                            abs(p - starts[ci]), abs(p - (ends[ci] - 1)))
                        if d < best_d:
                            best_d, best_pid = d, ci
                if best_d <= WINDOW:
                    peak_ids[gi] = best_pid
                    uids[gi] = uid_base[key] + best_pid
        return peak_ids, uids

    junc_ids = df["junc_id"].astype(str).values
    chrom_arr  = np.array([_parse_chrom(j) for j in junc_ids])
    strand_arr = np.array([_parse_strand(j) for j in junc_ids])
    tss_pos = df["tss_pos"].values.astype(np.int64)
    tts_pos = df["tts_pos"].values.astype(np.int64)

    tss_pid, tss_uid = _lookup(tss_pos, chrom_arr, strand_arr, cage_flat, cage_uid_base)
    tts_pid, tts_uid = _lookup(tts_pos, chrom_arr, strand_arr, drna_flat, drna_uid_base)

    df = df.copy()
    df["_tss_peak_id"] = tss_pid   # local peak index within (chrom,strand), -1=miss
    df["_tts_peak_id"] = tts_pid
    df["_tss_uid"]     = tss_uid   # global uid for recall counting
    df["_tts_uid"]     = tts_uid
    df["_chrom"]       = chrom_arr
    df["_strand"]      = strand_arr
    df["_total_cage"]  = total_cage
    df["_total_drna"]  = total_drna
    return df


def _dedup_pr_fast(
    mask: np.ndarray,
    junc_ids: np.ndarray,
    tss_pid: np.ndarray,
    tts_pid: np.ndarray,
    tss_uid: np.ndarray,
    tts_uid: np.ndarray,
    total_cage: int,
    total_drna: int,
) -> Tuple[float, float, float, float, float, float]:
    """Compute dedup P/R from boolean survivor mask + precomputed arrays.

    No Python loops over rows — uses numpy groupby via structured arrays.
    """
    if not mask.any():
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    jids  = junc_ids[mask]
    t_pid = tss_pid[mask]
    s_pid = tts_pid[mask]
    t_uid = tss_uid[mask]
    s_uid = tts_uid[mask]

    # Precision: per junc_id group, count unique (tss_peak_id, tts_peak_id) pairs
    # Encode pair as single int: use (tss_pid+1)*100003 + (tts_pid+1) (both >= -1)
    pair_enc = (t_pid.astype(np.int64) + 1) * 100003 + (s_pid.astype(np.int64) + 1)

    # Sort by junc_id then pair_enc to find unique pairs per group efficiently
    order = np.lexsort((pair_enc, jids))
    jids_s    = jids[order]
    t_pid_s   = t_pid[order]
    s_pid_s   = s_pid[order]
    pair_s    = pair_enc[order]

    # Find group boundaries
    boundaries = np.concatenate(([True], jids_s[1:] != jids_s[:-1], [True]))
    group_start = np.where(boundaries[:-1])[0]
    group_end   = np.where(boundaries[1:])[0]

    total_pairs = 0
    total_5p_tp = 0
    total_3p_tp = 0
    for gs, ge in zip(group_start, group_end):
        sl = slice(gs, ge + 1)
        # unique pairs within group
        up = np.unique(pair_s[sl])
        total_pairs += len(up)
        # reconstruct tss_pid from encoding: tss_pid = (enc // 100003) - 1
        up_tss = (up // 100003) - 1
        up_tts = (up %  100003) - 1
        total_5p_tp += int((up_tss >= 0).sum())
        total_3p_tp += int((up_tts >= 0).sum())

    prec_5 = total_5p_tp / total_pairs if total_pairs > 0 else 0.0
    prec_3 = total_3p_tp / total_pairs if total_pairs > 0 else 0.0

    # Recall: unique peak uids touched by any survivor
    rec_5 = len(np.unique(t_uid[t_uid >= 0])) / total_cage if total_cage > 0 else 0.0
    rec_3 = len(np.unique(s_uid[s_uid >= 0])) / total_drna if total_drna > 0 else 0.0

    def _f1(p, r):
        return 2 * p * r / (p + r) if (p + r) > 0 else 0.0

    return prec_5, rec_5, _f1(prec_5, rec_5), prec_3, rec_3, _f1(prec_3, rec_3)


def _make_arrays(df: pd.DataFrame):
    """Extract numpy arrays from an _assign_peak_ids-enriched DataFrame."""
    return (
        df["junc_id"].astype(str).values,
        df["_tss_peak_id"].values.astype(np.int64),
        df["_tts_peak_id"].values.astype(np.int64),
        df["_tss_uid"].values.astype(np.int64),
        df["_tts_uid"].values.astype(np.int64),
        int(df["_total_cage"].iloc[0]),
        int(df["_total_drna"].iloc[0]),
    )


def _score_arrays(df: pd.DataFrame):
    """Return raw component score arrays (float64) for vectorised threshold application."""
    return (
        df["TED_depth"].values.astype(np.float64),
        df["TED_tss_model"].values.astype(np.float64),
        df["TED_tss_annot"].values.astype(np.float64),
        df["TED_tts_model"].values.astype(np.float64),
        df["TED_tts_annot"].values.astype(np.float64),
    )


# ── Weight sweep ─────────────────────────────────────────────────────────────

def _sweep_weights(
    df: pd.DataFrame,
    end: str,
    n_steps: int = 11,
    fixed_other_weights: dict = None,
    fixed_other_thresh: float = 0.5,
) -> dict:
    """Sweep (w_depth, w_model, w_annot) for one end, holding the other fixed.

    Uses precomputed peak arrays — no per-row Python loops.
    Returns arrays: w_depth, w_model, w_annot, best_f1, best_thresh.
    """
    if fixed_other_weights is None:
        fixed_other_weights = {"wd": 0.33, "wm": 0.34, "wa": 0.33}

    jids, tss_pid, tts_pid, tss_uid, tts_uid, total_cage, total_drna = _make_arrays(df)
    depth, tss_model, tss_annot, tts_model, tts_annot = _score_arrays(df)

    fo = fixed_other_weights
    fixed_tss = fo["wd"] * depth + fo["wm"] * tss_model + fo["wa"] * tss_annot
    fixed_tts = fo["wd"] * depth + fo["wm"] * tts_model + fo["wa"] * tts_annot

    grid = np.linspace(0, 1, n_steps)
    thresh_grid = np.linspace(0.05, 0.95, 19)
    w_ds, w_ms, w_as, f1s, thresholds = [], [], [], [], []

    for wd in grid:
        for wm in grid:
            wa = 1.0 - wd - wm
            if wa < -0.01:
                continue
            wa = max(0.0, wa)

            sweep_score = wd * depth + wm * (tss_model if end == "tss" else tts_model) \
                                     + wa * (tss_annot if end == "tss" else tts_annot)

            best_f1, best_t = 0.0, 0.5
            for t in thresh_grid:
                if end == "tss":
                    mask = (sweep_score >= t) & (fixed_tts >= fixed_other_thresh)
                    _, _, f1, _, _, _ = _dedup_pr_fast(
                        mask, jids, tss_pid, tts_pid, tss_uid, tts_uid, total_cage, total_drna)
                else:
                    mask = (fixed_tss >= fixed_other_thresh) & (sweep_score >= t)
                    _, _, _, _, _, f1 = _dedup_pr_fast(
                        mask, jids, tss_pid, tts_pid, tss_uid, tts_uid, total_cage, total_drna)
                if f1 > best_f1:
                    best_f1, best_t = f1, float(t)

            w_ds.append(wd)
            w_ms.append(wm)
            w_as.append(wa)
            f1s.append(best_f1)
            thresholds.append(best_t)

    return {
        "w_depth":    np.asarray(w_ds),
        "w_model":    np.asarray(w_ms),
        "w_annot":    np.asarray(w_as),
        "best_f1":    np.asarray(f1s),
        "best_thresh": np.asarray(thresholds),
    }


def _joint_threshold_sweep(
    df: pd.DataFrame,
    opt_5: dict, opt_3: dict,
    n_steps: int = 21,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """2D (t_tss × t_tts) sweep using optimal weights. All numpy, no iterrows."""
    jids, tss_pid, tts_pid, tss_uid, tts_uid, total_cage, total_drna = _make_arrays(df)
    depth, tss_model, tss_annot, tts_model, tts_annot = _score_arrays(df)

    tss_score = (opt_5["w_depth"] * depth
                 + opt_5["w_model"] * tss_model
                 + opt_5["w_annot"] * tss_annot)
    tts_score = (opt_3["w_depth"] * depth
                 + opt_3["w_model"] * tts_model
                 + opt_3["w_annot"] * tts_annot)

    grid = np.linspace(0.05, 0.95, n_steps)
    f1_5_grid  = np.zeros((n_steps, n_steps))
    f1_3_grid  = np.zeros((n_steps, n_steps))
    prec5_grid = np.zeros((n_steps, n_steps))
    prec3_grid = np.zeros((n_steps, n_steps))
    rec5_grid  = np.zeros((n_steps, n_steps))
    rec3_grid  = np.zeros((n_steps, n_steps))

    for ti, t_tss in enumerate(grid):
        tss_mask = tss_score >= t_tss
        for tj, t_tts in enumerate(grid):
            mask = tss_mask & (tts_score >= t_tts)
            p5, r5, f5, p3, r3, f3 = _dedup_pr_fast(
                mask, jids, tss_pid, tts_pid, tss_uid, tts_uid, total_cage, total_drna)
            f1_5_grid[ti, tj]  = f5
            f1_3_grid[ti, tj]  = f3
            prec5_grid[ti, tj] = p5
            prec3_grid[ti, tj] = p3
            rec5_grid[ti, tj]  = r5
            rec3_grid[ti, tj]  = r3

    return grid, f1_5_grid, f1_3_grid, prec5_grid, prec3_grid, rec5_grid, rec3_grid


# ── Plots ────────────────────────────────────────────────────────────────────

def plot_weight_sweep(sw: dict, label: str, end_label: str, outdir: Path) -> dict:
    """Ternary-style 2D heatmap: x=w_model, y=w_depth, colour=best_F1."""
    best_idx = int(np.argmax(sw["best_f1"]))
    opt = {
        "w_depth":   float(sw["w_depth"][best_idx]),
        "w_model":   float(sw["w_model"][best_idx]),
        "w_annot":   float(sw["w_annot"][best_idx]),
        "f1":        float(sw["best_f1"][best_idx]),
        "threshold": float(sw["best_thresh"][best_idx]),
    }

    fig, ax = plt.subplots(figsize=(W1 * 1.1, W1 * 1.1))
    sc = ax.scatter(sw["w_model"], sw["w_depth"], c=sw["best_f1"],
                    cmap="viridis", s=18, edgecolors="none",
                    vmin=0, vmax=max(float(sw["best_f1"].max()), 0.5))
    ax.scatter([opt["w_model"]], [opt["w_depth"]], marker="*",
               s=120, color="red", edgecolors="black", linewidths=0.5, zorder=5)
    ax.text(0.03, 0.03,
            f"Best F1={opt['f1']:.3f}\n"
            f"d={opt['w_depth']:.2f} m={opt['w_model']:.2f} a={opt['w_annot']:.2f}\n"
            f"t={opt['threshold']:.2f}",
            transform=ax.transAxes, fontsize=5, va="bottom",
            bbox=dict(facecolor="white", alpha=0.85, pad=1, edgecolor="none"))
    ax.fill_between([0, 1], [1, 0], [1.05, 1.05], color="grey", alpha=0.15)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    fig.colorbar(sc, ax=ax, shrink=0.8, pad=0.02, label="Best F1 (dedup)")
    style_ax(ax, xlabel="w_model", ylabel="w_depth",
             title=f"{label} — {end_label}\n(w_annot = 1 − w_depth − w_model)")
    fig.tight_layout()
    savefig(fig, outdir / f"weight_sweep_{end_label}_{label}")
    return opt


def plot_joint_sweep(
    grid, f1_5, f1_3, prec5, prec3, rec5, rec3,
    label: str, outdir: Path,
) -> Tuple[dict, dict]:
    """Two heatmaps side-by-side: joint F1 for 5' and 3' ends."""
    fig, (ax5, ax3) = plt.subplots(1, 2, figsize=(W2 * 0.75, W1 * 1.1))

    def _draw(ax, grid, f1_grid, title):
        best_flat = int(np.argmax(f1_grid))
        bi, bj = np.unravel_index(best_flat, f1_grid.shape)
        im = ax.imshow(f1_grid.T, origin="lower", aspect="auto",
                       extent=[grid[0], grid[-1], grid[0], grid[-1]],
                       cmap="viridis", vmin=0, vmax=max(float(f1_grid.max()), 0.5))
        cs = ax.contour(grid, grid, f1_grid.T, levels=6,
                        colors="white", linewidths=0.4, alpha=0.6)
        ax.clabel(cs, fmt="%.2f", fontsize=4, inline=True)
        ax.scatter([grid[bi]], [grid[bj]], marker="*", s=120,
                   color="red", edgecolors="black", linewidths=0.5, zorder=5)
        ax.text(0.03, 0.97,
                f"F1={f1_grid[bi,bj]:.3f}\n"
                f"t_tss={grid[bi]:.2f}  t_tts={grid[bj]:.2f}",
                transform=ax.transAxes, fontsize=5, va="top",
                bbox=dict(facecolor="white", alpha=0.85, pad=1, edgecolor="none"))
        style_ax(ax, xlabel="t_tss", ylabel="t_tts", title=title)
        return im, int(bi), int(bj)

    im5, bi5, bj5 = _draw(ax5, grid, f1_5, f"5′ F1 — {label}")
    im3, bi3, bj3 = _draw(ax3, grid, f1_3, f"3′ F1 — {label}")
    fig.colorbar(im3, ax=[ax5, ax3], shrink=0.7, pad=0.02, label="F1 (dedup)")
    fig.suptitle("Joint threshold sweep (optimal weights, dedup P/R metric)",
                 fontsize=7, y=1.02)
    fig.tight_layout()
    savefig(fig, outdir / f"joint_threshold_sweep_{label}")

    opt_5 = {"t_tss": float(grid[bi5]), "t_tts": float(grid[bj5]),
              "f1": float(f1_5[bi5, bj5]),
              "precision": float(prec5[bi5, bj5]), "recall": float(rec5[bi5, bj5])}
    opt_3 = {"t_tss": float(grid[bi3]), "t_tts": float(grid[bj3]),
              "f1": float(f1_3[bi3, bj3]),
              "precision": float(prec3[bi3, bj3]), "recall": float(rec3[bi3, bj3])}
    return opt_5, opt_3


# ── Summary TSV ──────────────────────────────────────────────────────────────

def write_summary(
    label: str,
    opt_w5: dict, opt_w3: dict,
    opt_j5: dict, opt_j3: dict,
    outdir: Path,
):
    path = outdir / "summary_stats.tsv"
    first_write = not path.exists()
    with open(path, "a", newline="") as f:
        fields = [
            "config",
            "opt_5prime_w_depth",  "opt_5prime_w_model",  "opt_5prime_w_annot",
            "opt_5prime_f1_wsweep", "opt_5prime_thresh_wsweep",
            "opt_3prime_w_depth",  "opt_3prime_w_model",  "opt_3prime_w_annot",
            "opt_3prime_f1_wsweep", "opt_3prime_thresh_wsweep",
            "joint_5prime_t_tss",  "joint_5prime_t_tts",
            "joint_5prime_f1",     "joint_5prime_precision", "joint_5prime_recall",
            "joint_3prime_t_tss",  "joint_3prime_t_tts",
            "joint_3prime_f1",     "joint_3prime_precision", "joint_3prime_recall",
        ]
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        if first_write:
            writer.writeheader()
        row = {
            "config":                   label,
            "opt_5prime_w_depth":       f"{opt_w5['w_depth']:.2f}",
            "opt_5prime_w_model":       f"{opt_w5['w_model']:.2f}",
            "opt_5prime_w_annot":       f"{opt_w5['w_annot']:.2f}",
            "opt_5prime_f1_wsweep":     f"{opt_w5['f1']:.4f}",
            "opt_5prime_thresh_wsweep": f"{opt_w5['threshold']:.2f}",
            "opt_3prime_w_depth":       f"{opt_w3['w_depth']:.2f}",
            "opt_3prime_w_model":       f"{opt_w3['w_model']:.2f}",
            "opt_3prime_w_annot":       f"{opt_w3['w_annot']:.2f}",
            "opt_3prime_f1_wsweep":     f"{opt_w3['f1']:.4f}",
            "opt_3prime_thresh_wsweep": f"{opt_w3['threshold']:.2f}",
            "joint_5prime_t_tss":       f"{opt_j5['t_tss']:.2f}",
            "joint_5prime_t_tts":       f"{opt_j5['t_tts']:.2f}",
            "joint_5prime_f1":          f"{opt_j5['f1']:.4f}",
            "joint_5prime_precision":   f"{opt_j5['precision']:.4f}",
            "joint_5prime_recall":      f"{opt_j5['recall']:.4f}",
            "joint_3prime_t_tss":       f"{opt_j3['t_tss']:.2f}",
            "joint_3prime_t_tts":       f"{opt_j3['t_tts']:.2f}",
            "joint_3prime_f1":          f"{opt_j3['f1']:.4f}",
            "joint_3prime_precision":   f"{opt_j3['precision']:.4f}",
            "joint_3prime_recall":      f"{opt_j3['recall']:.4f}",
        }
        writer.writerow(row)
    log.info("Summary → %s", path)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--ted-log", nargs="+", required=True,
                        help="label:path pairs for TED log TSVs")
    parser.add_argument("--cage-peaks", required=True,
                        help="BED6 CAGE peaks (5′ orthogonal signal)")
    parser.add_argument("--drna-peaks", required=True,
                        help="BED6 dRNA peaks (3′ orthogonal signal)")
    parser.add_argument("--output", required=True,
                        help="Output directory")
    parser.add_argument("--weight-steps", type=int, default=11,
                        help="Grid steps for weight sweep (default: 11, ~66 points)")
    parser.add_argument("--threshold-steps", type=int, default=21,
                        help="Grid steps for threshold sweep (default: 21)")
    parser.add_argument("--chrom", default=None,
                        help="Restrict analysis to this chromosome (e.g. chr22). "
                             "Strongly recommended to keep runtime manageable.")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    chrom = args.chrom
    if chrom:
        log.info("Restricting analysis to chromosome: %s", chrom)

    logs = load_all_logs(args.ted_log, chrom=chrom)
    if not logs:
        log.error("No TED log data loaded")
        sys.exit(1)

    log.info("Loading CAGE peaks from %s", args.cage_peaks)
    cage_peaks = _parse_peaks(args.cage_peaks, chrom=chrom)
    log.info("Loading dRNA peaks from %s", args.drna_peaks)
    drna_peaks = _parse_peaks(args.drna_peaks, chrom=chrom)

    total_cage = sum(len(v) for v in cage_peaks.values())
    total_drna = sum(len(v) for v in drna_peaks.values())
    log.info("  %d CAGE peaks, %d dRNA peaks", total_cage, total_drna)

    for label, df in logs.items():
        log.info("── %s (%d candidates) ──", label, len(df))

        # Precompute peak hit IDs once — all sweep iterations reuse these arrays
        log.info("  Precomputing peak hits …")
        df = _assign_peak_ids(df, cage_peaks, drna_peaks)

        # ── Step 1: weight sweep for 5' (TSS) ──────────────────────────
        log.info("  Weight sweep 5′ …")
        sw5 = _sweep_weights(
            df, "tss",
            n_steps=args.weight_steps,
        )
        opt_w5 = plot_weight_sweep(sw5, label, "5prime", outdir)
        log.info("  Opt 5′: d=%.2f m=%.2f a=%.2f  F1=%.4f  t=%.2f",
                 opt_w5["w_depth"], opt_w5["w_model"], opt_w5["w_annot"],
                 opt_w5["f1"], opt_w5["threshold"])

        # ── Step 2: weight sweep for 3' (TTS), using opt 5' weights ────
        log.info("  Weight sweep 3′ …")
        sw3 = _sweep_weights(
            df, "tts",
            n_steps=args.weight_steps,
            fixed_other_weights={
                "wd": opt_w5["w_depth"],
                "wm": opt_w5["w_model"],
                "wa": opt_w5["w_annot"],
            },
            fixed_other_thresh=opt_w5["threshold"],
        )
        opt_w3 = plot_weight_sweep(sw3, label, "3prime", outdir)
        log.info("  Opt 3′: d=%.2f m=%.2f a=%.2f  F1=%.4f  t=%.2f",
                 opt_w3["w_depth"], opt_w3["w_model"], opt_w3["w_annot"],
                 opt_w3["f1"], opt_w3["threshold"])

        # ── Step 3: joint threshold sweep ──────────────────────────────
        log.info("  Joint threshold sweep …")
        grid, f1_5, f1_3, p5, p3, r5, r3 = _joint_threshold_sweep(
            df,
            {"w_depth": opt_w5["w_depth"], "w_model": opt_w5["w_model"], "w_annot": opt_w5["w_annot"]},
            {"w_depth": opt_w3["w_depth"], "w_model": opt_w3["w_model"], "w_annot": opt_w3["w_annot"]},
            n_steps=args.threshold_steps,
        )
        opt_j5, opt_j3 = plot_joint_sweep(
            grid, f1_5, f1_3, p5, p3, r5, r3, label, outdir)
        log.info("  Joint opt 5′: t_tss=%.2f t_tts=%.2f  F1=%.4f",
                 opt_j5["t_tss"], opt_j5["t_tts"], opt_j5["f1"])
        log.info("  Joint opt 3′: t_tss=%.2f t_tts=%.2f  F1=%.4f",
                 opt_j3["t_tss"], opt_j3["t_tts"], opt_j3["f1"])

        write_summary(label, opt_w5, opt_w3, opt_j5, opt_j3, outdir)

    log.info("Done — output in %s", outdir)


if __name__ == "__main__":
    main()
