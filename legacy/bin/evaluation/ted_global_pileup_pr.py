#!/usr/bin/env python3
"""
ted_global_pileup_pr.py - Precision/recall audit for --ted_global pileup peaks.

Global TED emits two useful stages into <prefix>.ted_log.tsv:

  global_peak        Partition-level TSS/TTS pileup peaks after find_peaks,
                     min-read filtering, merge/rescue/drop logic, and TED
                     validation. Rows have one real end coordinate and one
                     sentinel -1 coordinate.
  global_assign      Per-junction-chain groups snapped to approved global
                     TSS/TTS peaks. These are downstream of the peak gate.

This script compares both stages to orthogonal CAGE/dRNA BED6 peaks so the
pipeline can separate:

  raw pileup peak precision/recall,
  validated peak precision/recall,
  assigned snapped-end precision/recall,
  and read-weighted precision for each stage.
"""

from __future__ import annotations

import argparse
import logging
import sys
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ted_log_loader import load_ted_log
# Shared, contract-compliant matcher: collapse peaks to one summit each and
# match called ends one-to-one (greedy nearest).  Same module the logreg
# fit uses, so "TP" means the same thing.  (R3 fix, 2026-05-28 memo §5.2.)
from peak_matching import (
    match_positions_one_to_one, count_summits,
    load_summit_index_from_bed,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s")
log = logging.getLogger(__name__)
logging.getLogger("matplotlib").setLevel(logging.WARNING)

# A (chrom, strand) -> sorted list of summit positions (one per peak).
PeakIndex = Dict[Tuple[str, str], List[int]]
Candidate = Dict[str, object]

STAGE_ORDER = ["global_peak_all", "global_peak_pass", "global_assign_pass"]
STAGE_LABEL = {
    "global_peak_all": "raw peaks",
    "global_peak_pass": "validated peaks",
    "global_assign_pass": "assigned ends",
}
END_LABEL = {"tss": "5prime", "tts": "3prime"}


def _parse_peaks(path: str, chroms: Optional[Iterable[str]] = None) -> PeakIndex:
    """Parse BED6 peaks into a (chrom, strand) -> sorted summit index.

    Each peak interval is collapsed to ONE summit (peak centre; no signal
    bedGraph wired in here) so downstream one-to-one matching can credit at
    most one called end per peak.  (R3 fix.)

    `chroms`, when given, restricts the truth peak set to those chromosomes.
    On a scoped run (candidates only exist on one chromosome) this keeps the
    recall denominator (total summits) on the same chromosome as the
    candidates, otherwise recall/F1 are deflated by the genome/chrom ratio.
    """
    return load_summit_index_from_bed(path, chroms=chroms)


def _parse_chrom(junc_id: str) -> str:
    return str(junc_id).split(":")[0]


def _parse_strand(junc_id: str) -> str:
    parts = str(junc_id).split(":")
    return parts[2] if len(parts) > 2 else "+"


def _as_int(value, default: int = -1) -> int:
    try:
        if pd.isna(value):
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _as_float(value, default: float = np.nan) -> float:
    try:
        if pd.isna(value):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _nearest_peak(pos: int, summits: List[int]
                  ) -> Tuple[Optional[int], Optional[int]]:
    """Return nearest summit position and its absolute distance.

    Diagnostic only — the TP/precision/recall accounting uses one-to-one
    matching in _score_candidates.  Operates on the summit index (one
    scalar position per peak).
    """
    if not summits:
        return None, None
    idx = bisect_left(summits, pos)
    best_pos = None
    best_dist = None
    for i in (idx - 1, idx):
        if 0 <= i < len(summits):
            dist = abs(pos - summits[i])
            if best_dist is None or dist < best_dist:
                best_dist = dist
                best_pos = summits[i]
    return best_pos, best_dist


def _candidate(label: str, stage: str, end_type: str, row: pd.Series,
               pos: int) -> Candidate:
    return {
        "label": label,
        "stage": stage,
        "end": end_type,
        "chrom": _parse_chrom(row.get("junc_id", "")),
        "strand": _parse_strand(row.get("junc_id", "")),
        "pos": pos,
        "status": str(row.get("status", "")),
        "n_reads": _as_int(row.get("n_reads", 0), 0),
        "TED_confidence": _as_float(row.get("TED_confidence")),
        "TED_depth": _as_float(row.get("TED_depth")),
        "TED_model": _as_float(row.get("TED_tss_model" if end_type == "tss"
                                        else "TED_tts_model")),
        "TED_annot": _as_float(row.get("TED_tss_annot" if end_type == "tss"
                                        else "TED_tts_annot")),
        "TED_reality": _as_float(row.get("TED_tss_reality" if end_type == "tss"
                                         else "TED_tts_reality")),
        "drop_reason": str(row.get("drop_reason", "")),
    }


def _load_candidates(label: str, path: str) -> Dict[Tuple[str, str], List[Candidate]]:
    df = load_ted_log(path, drop_partition_summary=False)
    buckets: Dict[Tuple[str, str], List[Candidate]] = defaultdict(list)
    if "stage" not in df.columns:
        return buckets

    global_peak = df[df["stage"] == "global_peak"].copy()
    if global_peak.empty:
        return buckets

    for _, row in global_peak.iterrows():
        tss_pos = _as_int(row.get("tss_pos"))
        tts_pos = _as_int(row.get("tts_pos"))
        if tss_pos >= 0:
            cand = _candidate(label, "global_peak_all", "tss", row, tss_pos)
            buckets[("global_peak_all", "tss")].append(cand)
            if cand["status"] == "pass":
                buckets[("global_peak_pass", "tss")].append(
                    {**cand, "stage": "global_peak_pass"})
        if tts_pos >= 0:
            cand = _candidate(label, "global_peak_all", "tts", row, tts_pos)
            buckets[("global_peak_all", "tts")].append(cand)
            if cand["status"] == "pass":
                buckets[("global_peak_pass", "tts")].append(
                    {**cand, "stage": "global_peak_pass"})

    assigned = df[(df["stage"] == "global_assign") & (df["status"] == "pass")]
    for _, row in assigned.iterrows():
        tss_pos = _as_int(row.get("tss_pos"))
        tts_pos = _as_int(row.get("tts_pos"))
        if tss_pos >= 0:
            buckets[("global_assign_pass", "tss")].append(
                _candidate(label, "global_assign_pass", "tss", row, tss_pos))
        if tts_pos >= 0:
            buckets[("global_assign_pass", "tts")].append(
                _candidate(label, "global_assign_pass", "tts", row, tts_pos))
    return buckets


def _score_candidates(candidates: List[Candidate], peaks: PeakIndex,
                      window: int) -> Tuple[dict, List[dict]]:
    """Score candidate ends against truth summits with one-to-one matching.

    R3 fix: candidates are bucketed by (chrom, strand) and matched to the
    summit index with the shared greedy-nearest one-to-one matcher, so each
    summit backs at most one candidate TP.  Precision = #candidates matched
    / #candidates; recall = #summits matched / #summits.  The
    nearest-peak-distance diagnostic columns (median/p90/within-Nbp) are
    kept as informational distance histograms, not TP accounting.
    """
    total_peaks = count_summits(peaks)
    detail_rows: List[dict] = []
    reads_total = 0
    reads_hit = 0
    distances = []
    unique_positions = {}

    # ── One-to-one matching, bucketed by (chrom, strand) ──────────────
    buckets: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    for i, cand in enumerate(candidates):
        key = (str(cand["chrom"]), str(cand["strand"]))
        buckets[key].append(i)
    matched_summit_pos: List[Optional[int]] = [None] * len(candidates)
    matched_peaks = set()
    for key, idxs in buckets.items():
        summits = peaks.get(key, [])
        if not summits:
            continue
        positions = [int(candidates[i]["pos"]) for i in idxs]
        assign = match_positions_one_to_one(positions, summits, window)
        for local_i, ci in enumerate(idxs):
            j = assign[local_i]
            if j >= 0:
                matched_summit_pos[ci] = summits[j]
                matched_peaks.add((key[0], key[1], summits[j]))

    candidate_hits = 0
    for i, cand in enumerate(candidates):
        key = (str(cand["chrom"]), str(cand["strand"]))
        # Diagnostic nearest-summit distance (independent of one-to-one hit).
        nearest_iv, dist = _nearest_peak(int(cand["pos"]), peaks.get(key, []))
        hit = matched_summit_pos[i] is not None
        if hit:
            candidate_hits += 1
        n_reads = int(cand.get("n_reads", 0) or 0)
        reads_total += n_reads
        if hit:
            reads_hit += n_reads
        if dist is not None:
            distances.append(dist)
        pos_key = (key[0], key[1], int(cand["pos"]))
        unique_positions[pos_key] = unique_positions.get(pos_key, False) or hit
        detail_rows.append({
            **cand,
            "orthogonal_hit": int(hit),
            "matched_summit": matched_summit_pos[i] if hit else -1,
            "nearest_peak_start": nearest_iv if nearest_iv is not None else -1,
            "nearest_peak_end": nearest_iv if nearest_iv is not None else -1,
            "nearest_peak_distance": dist if dist is not None else -1,
        })

    n_candidates = len(candidates)
    precision = candidate_hits / n_candidates if n_candidates else np.nan
    recall = len(matched_peaks) / total_peaks if total_peaks else np.nan
    weighted_precision = reads_hit / reads_total if reads_total else np.nan
    unique_hits = sum(1 for hit in unique_positions.values() if hit)
    unique_precision = (unique_hits / len(unique_positions)
                        if unique_positions else np.nan)
    f1 = (2 * precision * recall / (precision + recall)
          if precision == precision and recall == recall and (precision + recall) else np.nan)

    row = {
        "n_candidates": n_candidates,
        "n_candidate_hit": candidate_hits,
        "precision": precision,
        "n_unique_positions": len(unique_positions),
        "n_unique_position_hit": unique_hits,
        "unique_position_precision": unique_precision,
        "reads_total": reads_total,
        "reads_hit": reads_hit,
        "weighted_precision": weighted_precision,
        "n_orthogonal_peaks_hit": len(matched_peaks),
        "n_orthogonal_peaks_total": total_peaks,
        "recall": recall,
        "f1": f1,
        "median_nearest_peak_distance": float(np.median(distances)) if distances else np.nan,
        "p90_nearest_peak_distance": float(np.percentile(distances, 90)) if distances else np.nan,
        "within_0bp": sum(1 for d in distances if d == 0),
        "within_10bp": sum(1 for d in distances if d <= 10),
        "within_25bp": sum(1 for d in distances if d <= 25),
        "within_50bp": sum(1 for d in distances if d <= 50),
    }
    return row, detail_rows


def _plot_precision_recall(metrics: pd.DataFrame, outdir: Path) -> None:
    if metrics.empty:
        return
    labels = list(metrics["label"].drop_duplicates())
    colors = plt.cm.tab10(np.linspace(0, 1, max(1, len(labels))))
    color_map = dict(zip(labels, colors))
    markers = {"global_peak_all": "o", "global_peak_pass": "s",
               "global_assign_pass": "^"}

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.4), sharex=True, sharey=True)
    for ax, end_type in zip(axes, ["tss", "tts"]):
        sub = metrics[metrics["end"] == end_type]
        for _, row in sub.iterrows():
            ax.scatter(row["recall"], row["precision"],
                       s=max(25, min(220, np.log10(row["n_candidates"] + 1) * 55)),
                       marker=markers.get(row["stage"], "o"),
                       color=color_map[row["label"]],
                       edgecolor="black", linewidth=0.4, alpha=0.85)
            ax.text(row["recall"], row["precision"],
                    STAGE_LABEL.get(row["stage"], row["stage"]),
                    fontsize=5, ha="left", va="bottom")
        ax.set_title(END_LABEL[end_type])
        ax.set_xlabel("Orthogonal peak recall")
        ax.grid(True, linewidth=0.3, alpha=0.35)
    axes[0].set_ylabel("Candidate precision")
    handles = [plt.Line2D([0], [0], marker="o", color="w", label=label,
                          markerfacecolor=color_map[label], markeredgecolor="black",
                          markersize=6)
               for label in labels]
    fig.legend(handles=handles, loc="lower center", ncol=min(3, len(handles)),
               fontsize=6, frameon=False)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(outdir / "global_pileup_precision_recall.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)


def _plot_weighted_precision(metrics: pd.DataFrame, outdir: Path) -> None:
    if metrics.empty:
        return
    metrics = metrics.copy()
    metrics["x_label"] = metrics["label"] + "\n" + metrics["stage"].map(STAGE_LABEL)
    stages = [s for s in STAGE_ORDER if s in set(metrics["stage"])]
    labels = list(metrics["label"].drop_duplicates())
    x_keys = [(label, stage) for label in labels for stage in stages]
    x = np.arange(len(x_keys))
    width = 0.38

    fig, axes = plt.subplots(2, 1, figsize=(max(7, len(x_keys) * 0.8), 5.0),
                             sharex=True)
    for ax, metric, ylabel in [
        (axes[0], "weighted_precision", "Read-weighted precision"),
        (axes[1], "n_candidates", "Candidate count"),
    ]:
        for offset, end_type in [(-width / 2, "tss"), (width / 2, "tts")]:
            vals = []
            for label, stage in x_keys:
                sub = metrics[(metrics["label"] == label)
                              & (metrics["stage"] == stage)
                              & (metrics["end"] == end_type)]
                vals.append(float(sub[metric].iloc[0]) if not sub.empty else np.nan)
            ax.bar(x + offset, vals, width=width, label=END_LABEL[end_type])
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", linewidth=0.3, alpha=0.35)
    axes[0].legend(frameon=False, fontsize=6)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([f"{label}\n{STAGE_LABEL[stage]}" for label, stage in x_keys],
                            rotation=45, ha="right", fontsize=6)
    fig.tight_layout()
    fig.savefig(outdir / "global_pileup_weighted_precision_counts.png",
                dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ted-log", nargs="+", required=True,
                        help="label:path TED log entries")
    parser.add_argument("--cage-peaks", required=True)
    parser.add_argument("--drna-peaks", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--window", type=int, default=50)
    parser.add_argument(
        "--chrom", "--region", dest="chrom", default=None,
        help="Restrict truth CAGE/dRNA peaks to this chromosome (e.g. chr22). "
             "On scoped runs the global_peak candidates only exist on the "
             "scoped chromosome, so without this the recall denominator spans "
             "the whole genome and deflates recall/F1. Default: no restriction "
             "(genome-wide, preserves prior behaviour).")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)
    chroms = {args.chrom} if args.chrom else None
    if chroms:
        log.info("Restricting truth peaks to chromosome: %s", args.chrom)
    peaks_by_end = {
        "tss": _parse_peaks(args.cage_peaks, chroms=chroms),
        "tts": _parse_peaks(args.drna_peaks, chroms=chroms),
    }

    metric_rows = []
    detail_rows = []
    for entry in args.ted_log:
        if ":" not in entry:
            log.warning("Bad --ted-log entry %s; expected label:path", entry)
            continue
        label, path = entry.split(":", 1)
        if not Path(path).exists() or Path(path).stat().st_size == 0:
            log.warning("Skipping missing/empty TED log: %s", path)
            continue
        buckets = _load_candidates(label, path)
        if not buckets:
            log.info("No global_peak rows for %s; skipping", label)
            continue
        for stage in STAGE_ORDER:
            for end_type in ("tss", "tts"):
                candidates = buckets.get((stage, end_type), [])
                if not candidates:
                    continue
                row, details = _score_candidates(candidates, peaks_by_end[end_type],
                                                 args.window)
                row.update({"label": label, "stage": stage, "end": end_type,
                            "window": args.window})
                metric_rows.append(row)
                detail_rows.extend(details)

    if not metric_rows:
        (outdir / "global_pileup_metrics.tsv").write_text(
            "label\tstage\tend\tn_candidates\n")
        log.warning("No global TED pileup candidates found")
        return

    metrics = pd.DataFrame(metric_rows)
    metrics = metrics[["label", "stage", "end", "window",
                       "n_candidates", "n_candidate_hit", "precision",
                       "n_unique_positions", "n_unique_position_hit",
                       "unique_position_precision", "reads_total", "reads_hit",
                       "weighted_precision", "n_orthogonal_peaks_hit",
                       "n_orthogonal_peaks_total", "recall", "f1",
                       "median_nearest_peak_distance", "p90_nearest_peak_distance",
                       "within_0bp", "within_10bp", "within_25bp", "within_50bp"]]
    metrics.to_csv(outdir / "global_pileup_metrics.tsv", sep="\t",
                   index=False, float_format="%.6g")

    details = pd.DataFrame(detail_rows)
    if not details.empty:
        details.to_csv(outdir / "global_pileup_candidates.tsv", sep="\t",
                       index=False, float_format="%.6g")

    _plot_precision_recall(metrics, outdir)
    _plot_weighted_precision(metrics, outdir)
    log.info("Done - outputs in %s", outdir)


if __name__ == "__main__":
    main()
