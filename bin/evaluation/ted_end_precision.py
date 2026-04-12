#!/usr/bin/env python3
"""
eval_end_precision.py — Junction-chain-deduplicated end precision/recall.

For each junction chain (group of isoforms sharing all splice junctions):
  - Multiple isoform ends mapping to the SAME annotation window count as
    ONE true positive (deduplicated precision).
  - Different junction chains CAN independently match the same annotation
    window (each gets credit).

This prevents inflated precision from over-segmented ends within a locus
while still rewarding legitimate alternative splicing that shares an end.

Outputs:
  <outdir>/precision_recall_summary.tsv   — per-mode summary metrics
  <outdir>/per_junction_chain.tsv         — per-JC-group detail
"""

import argparse
import csv
import logging
import math
import sys
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

try:
    from signal_utils import parse_bed12, tss_tts
except ImportError:
    from evaluation.signal_utils import parse_bed12, tss_tts

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s")
log = logging.getLogger(__name__)


# ── Parse inputs ────────────────────────────────────────────────────────────

def parse_gtf_ends(gtf_path: str, region_chrom: str = None,
                   region_start: int = None, region_end: int = None
                   ) -> Dict[str, Dict[str, List[int]]]:
    """Extract annotated TSS/TTS positions from GTF, optionally filtered to a region."""
    ends: Dict[str, Dict[str, set]] = defaultdict(lambda: {"tss": set(), "tts": set()})
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            chrom = cols[0]
            if region_chrom and chrom != region_chrom:
                continue
            start = int(cols[3]) - 1  # GTF is 1-based
            end = int(cols[4])
            strand = cols[6]
            if region_start is not None and region_end is not None:
                if end < region_start or start > region_end:
                    continue
            if strand == "+":
                ends[chrom]["tss"].add(start)
                ends[chrom]["tts"].add(end)
            else:
                ends[chrom]["tss"].add(end)
                ends[chrom]["tts"].add(start)
    return {c: {"tss": sorted(ends[c]["tss"]), "tts": sorted(ends[c]["tts"])}
            for c in ends}


def parse_gtf_transcripts(gtf_path: str, region_chrom: str = None,
                          region_start: int = None, region_end: int = None
                          ) -> List[dict]:
    """Parse GTF transcript+exon features to get per-transcript junction chains and ends."""
    tx_exons: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    tx_info: Dict[str, Tuple[str, str, int, int]] = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom = cols[0]
            if region_chrom and chrom != region_chrom:
                continue
            feature = cols[2]
            if feature not in ("transcript", "exon"):
                continue
            start = int(cols[3]) - 1  # GTF is 1-based
            end = int(cols[4])
            strand = cols[6]
            if region_start is not None and region_end is not None:
                if end < region_start or start > region_end:
                    continue
            # Extract transcript_id
            tx_id = None
            for attr in cols[8].split(";"):
                attr = attr.strip()
                if attr.startswith("transcript_id"):
                    parts = attr.split('"')
                    if len(parts) >= 2:
                        tx_id = parts[1]
                    break
            if tx_id is None:
                continue
            if feature == "transcript":
                tx_info[tx_id] = (chrom, strand, start, end)
            elif feature == "exon":
                tx_exons[tx_id].append((start, end))
    transcripts = []
    for tx_id, exons in tx_exons.items():
        if tx_id not in tx_info:
            continue
        chrom, strand, tx_start, tx_end = tx_info[tx_id]
        exons_sorted = sorted(exons)
        introns = tuple((exons_sorted[i][1], exons_sorted[i + 1][0])
                        for i in range(len(exons_sorted) - 1))
        n_exons = len(exons_sorted)
        tss, tts = tss_tts(tx_start, tx_end, strand)
        transcripts.append({
            "tx_id": tx_id, "chrom": chrom, "strand": strand,
            "junctions": introns, "tss": tss, "tts": tts,
            # Fields expected by compute_jc_deduplicated_precision_recall:
            "n_exons": n_exons,
            "start": tx_start,
            "end": tx_end,
        })
    return transcripts


def parse_isoforms_bed(bed_path: str) -> List[dict]:
    """Parse BED12 isoforms into list of dicts with junction chains."""
    return parse_bed12(bed_path)


def parse_isoforms_gtf(gtf_path: str, region_chrom: str = None,
                       region_start: int = None, region_end: int = None
                       ) -> List[dict]:
    """Parse assembler GTF (Bambu/IsoQuant/StringTie2/etc.) into same format as parse_bed12.

    Reuses parse_gtf_transcripts — already reads exon features and builds
    junction chains.  Returns a list of transcript dicts compatible with
    compute_jc_deduplicated_precision_recall.
    """
    return parse_gtf_transcripts(gtf_path, region_chrom, region_start, region_end)


def parse_peaks_bed(path: str, region_chrom: str = None,
                    region_start: int = None, region_end: int = None
                    ) -> Dict[Tuple[str, str], List[Tuple[int, int]]]:
    """Parse BED6 peaks into {(chrom, strand): sorted (start, end) intervals}.

    Intervals use BED half-open convention [start, end).  Sorted by start so
    binary search can locate nearby peaks efficiently.
    """
    peaks: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            chrom, start, end, strand = cols[0], int(cols[1]), int(cols[2]), cols[5]
            if region_chrom and chrom != region_chrom:
                continue
            if region_start is not None and region_end is not None:
                if end < region_start or start > region_end:
                    continue
            peaks[(chrom, strand)].append((start, end))
    for k in peaks:
        peaks[k].sort()
    return peaks


# ── Core metric ─────────────────────────────────────────────────────────────

def _nearest_pos(pos: int, sorted_positions: List[int], window: int
                 ) -> Optional[int]:
    """Return the nearest scalar position within *window*, or None.

    Used for GTF-based fallback matching where only integer end positions
    are available (no interval width).
    """
    if not sorted_positions:
        return None
    idx = bisect_left(sorted_positions, pos)
    best_pos = None
    best_dist = window + 1
    for i in (idx - 1, idx):
        if 0 <= i < len(sorted_positions):
            d = abs(pos - sorted_positions[i])
            if d < best_dist:
                best_dist = d
                best_pos = sorted_positions[i]
    if best_dist <= window:
        return best_pos
    return None


def _nearest_annot(pos: int, sorted_intervals: List[Tuple[int, int]], window: int
                   ) -> Optional[Tuple[int, int]]:
    """Return the nearest peak interval within *window* bp of *pos*, or None.

    Distance is the minimum distance from *pos* to any point inside the
    half-open interval [start, end).  Matches interval-edge semantics of
    bedtools closest used in ted_core.py — more permissive than midpoint
    matching for wide peaks (e.g. CAGE).  Returns the (start, end) tuple
    as a hashable peak identity key.
    """
    if not sorted_intervals:
        return None
    starts = [iv[0] for iv in sorted_intervals]
    idx = bisect_left(starts, pos)
    best_iv = None
    best_dist = window + 1
    for i in (idx - 1, idx):
        if 0 <= i < len(sorted_intervals):
            s, e = sorted_intervals[i]
            d = 0 if s <= pos < e else min(abs(pos - s), abs(pos - (e - 1)))
            if d < best_dist:
                best_dist = d
                best_iv = sorted_intervals[i]
    if best_dist <= window:
        return best_iv
    return None


def compute_jc_deduplicated_precision_recall(
    isoforms: List[dict],
    annotated_ends: Dict[str, Dict[str, List[int]]],
    annot_transcripts: List[dict],
    window: int = 50,
    peaks_5prime: Dict[Tuple[str, str], List[int]] = None,
    peaks_3prime: Dict[Tuple[str, str], List[int]] = None,
) -> dict:
    """
    Compute junction-chain-deduplicated precision and recall.

    Precision (pair-aware deduplicated):
      For each junction chain group, compute the (tss_peak, tts_peak) pair
      for every isoform.  Deduplicate by pair — an isoform is only redundant
      if another isoform in the same JC group matches the SAME pair of peaks.
      Then:
        5' dedup TP = unique pairs with tss_peak matched
        3' dedup TP = unique pairs with tts_peak matched
      This means same-TSS-different-TTS isoforms both count as 5' TPs,
      and same-TTS-different-TSS both count as 3' TPs.

    Recall:
      When peaks are provided: fraction of peaks matched by any isoform.
      When using GTF: JC-matched recall (only annotation ends whose
      junction chain matches a predicted isoform).

    Returns dict with per-end-type and aggregate metrics, plus per-JC detail.
    """
    # Group multi-exon isoforms by junction chain.
    # Group single-exon isoforms by genomic overlap on the same strand:
    # overlapping single-exon isoforms are alternative-end variants of the
    # same locus and must be deduplicated the same way as JC groups.
    jc_groups: Dict[Tuple, List[dict]] = defaultdict(list)
    for iso in isoforms:
        if iso["n_exons"] >= 2 and iso["junctions"]:
            key = (iso["chrom"], iso["strand"], iso["junctions"])
            jc_groups[key].append(iso)
        else:
            iso["_se"] = True  # mark for sweep-line pass

    # Sweep-line overlap grouping for single-exon isoforms.
    se_isos = [iso for iso in isoforms if iso.get("_se")]
    if se_isos:
        se_isos_sorted = sorted(se_isos, key=lambda x: (x["chrom"], x["strand"], x["start"]))
        group_id = 0
        current_end = -1
        current_chrom = None
        current_strand = None
        for iso in se_isos_sorted:
            chrom, strand = iso["chrom"], iso["strand"]
            start, end = iso["start"], iso["end"]
            if chrom != current_chrom or strand != current_strand or start >= current_end:
                group_id += 1
                current_chrom = chrom
                current_strand = strand
                current_end = end
            else:
                current_end = max(current_end, end)
            jc_groups[(chrom, strand, "SE", group_id)].append(iso)

    # Build set of predicted JC keys for recall filtering (multi-exon only)
    predicted_jc_keys = set(k for k in jc_groups if len(k) == 3)

    # Collect annotation ends only from transcripts whose JC matches a predicted one
    relevant_annot_ends: Dict[str, Dict[str, Set[int]]] = defaultdict(
        lambda: {"tss": set(), "tts": set()}
    )
    for tx in annot_transcripts:
        jc_key = (tx["chrom"], tx["strand"], tx["junctions"])
        if jc_key in predicted_jc_keys:
            relevant_annot_ends[tx["chrom"]]["tss"].add(tx["tss"])
            relevant_annot_ends[tx["chrom"]]["tts"].add(tx["tts"])

    # Clean up temporary marker
    for iso in isoforms:
        iso.pop("_se", None)

    results = {"per_jc": []}

    # Build sorted lists of JC-relevant annotation ends for recall matching
    relevant_annot_sorted: Dict[str, Dict[str, List[int]]] = {}
    for chrom in relevant_annot_ends:
        relevant_annot_sorted[chrom] = {
            "tss": sorted(relevant_annot_ends[chrom]["tss"]),
            "tts": sorted(relevant_annot_ends[chrom]["tts"]),
        }

    # ── Pair-aware per-JC precision ─────────────────────────────────────
    # Process both ends together per JC group so dedup uses the full pair.
    total_5p_ends = 0
    total_3p_ends = 0
    total_5p_dedup_tp = 0
    total_3p_dedup_tp = 0
    total_5p_redundant = 0
    total_3p_redundant = 0
    total_paired_isoforms = 0
    total_paired_unique = 0
    total_paired_both_hit = 0
    matched_annots_global_tss: Set[Tuple[str, int]] = set()
    matched_annots_global_tts: Set[Tuple[str, int]] = set()

    for jc_key, members in jc_groups.items():
        chrom, strand = jc_key[0], jc_key[1]
        is_se_group = (len(jc_key) == 4 and jc_key[2] == "SE")
        junctions = jc_key[2] if (len(jc_key) == 3) else ()

        # When orthogonal peaks provided, use interval-edge matching (_nearest_annot).
        # GTF fallback uses scalar positions (_nearest_pos).
        peaks_5ref  = peaks_5prime.get((chrom, strand), []) if peaks_5prime is not None else None
        peaks_3ref  = peaks_3prime.get((chrom, strand), []) if peaks_3prime is not None else None
        gtf_tss_ref = annotated_ends.get(chrom, {}).get("tss", [])
        gtf_tts_ref = annotated_ends.get(chrom, {}).get("tts", [])

        tss_relevant = relevant_annot_sorted.get(chrom, {}).get("tss", [])
        tts_relevant = relevant_annot_sorted.get(chrom, {}).get("tts", [])

        end_pairs: List[Tuple] = []
        tss_positions = []
        tts_positions = []
        for iso in members:
            tss_match = (_nearest_annot(iso["tss"], peaks_5ref, window)
                         if peaks_5ref is not None
                         else _nearest_pos(iso["tss"], gtf_tss_ref, window))
            tts_match = (_nearest_annot(iso["tts"], peaks_3ref, window)
                         if peaks_3ref is not None
                         else _nearest_pos(iso["tts"], gtf_tts_ref, window))
            end_pairs.append((tss_match, tts_match))
            tss_positions.append(iso["tss"])
            tts_positions.append(iso["tts"])

            # GTF-based recall (only for multi-exon, scalar positions)
            if peaks_5prime is None and not is_se_group:
                rel_m = _nearest_pos(iso["tss"], tss_relevant, window)
                if rel_m is not None:
                    matched_annots_global_tss.add((chrom, rel_m))
            if peaks_3prime is None and not is_se_group:
                rel_m = _nearest_pos(iso["tts"], tts_relevant, window)
                if rel_m is not None:
                    matched_annots_global_tts.add((chrom, rel_m))

        unique_pairs = set(end_pairs)
        n_members = len(members)

        n_5p_tp = sum(1 for p in unique_pairs if p[0] is not None)
        n_3p_tp = sum(1 for p in unique_pairs if p[1] is not None)
        n_both_hit = sum(1 for p in unique_pairs if p[0] is not None and p[1] is not None)

        n_unique = len(unique_pairs)
        n_5p_unique_pos = len(set(tss_positions))
        n_3p_unique_pos = len(set(tts_positions))
        n_5p_redundant = max(0, n_members - max(1, n_unique))
        n_3p_redundant = max(0, n_members - max(1, n_unique))

        total_5p_ends += n_members
        total_3p_ends += n_members
        total_5p_dedup_tp += n_5p_tp
        total_3p_dedup_tp += n_3p_tp
        total_5p_redundant += n_5p_redundant
        total_3p_redundant += n_3p_redundant
        total_paired_isoforms += n_members
        total_paired_unique += n_unique
        total_paired_both_hit += n_both_hit

        results["per_jc"].append({
            "chrom": chrom, "strand": strand,
            "jc_hash": f"SE{jc_key[3]:04d}" if is_se_group else f"{hash(junctions) & 0xFFFF:04x}",
            "n_isoforms": n_members,
            "n_junctions": 0 if is_se_group else len(junctions),
            "5prime_unique_pos": n_5p_unique_pos,
            "5prime_dedup_tp": n_5p_tp,
            "5prime_redundant": n_5p_redundant,
            "3prime_unique_pos": n_3p_unique_pos,
            "3prime_dedup_tp": n_3p_tp,
            "3prime_redundant": n_3p_redundant,
            "paired_unique": n_unique,
            "paired_both_hit": n_both_hit,
            "paired_redundant": n_members - n_unique,
        })

    # ── Per-end precision ───────────────────────────────────────────────
    for end_label, total_ends, total_tp, total_red in [
        ("5prime", total_5p_ends, total_5p_dedup_tp, total_5p_redundant),
        ("3prime", total_3p_ends, total_3p_dedup_tp, total_3p_redundant),
    ]:
        dedup_precision = total_tp / total_ends if total_ends > 0 else None
        naive_precision = (total_tp + total_red) / total_ends if total_ends > 0 else None

        # ── Recall ──────────────────────────────────────────────────────
        end_type = "tss" if end_label == "5prime" else "tts"
        current_peaks = peaks_5prime if end_type == "tss" else peaks_3prime
        use_peaks = current_peaks is not None
        matched_annots_global = matched_annots_global_tss if end_type == "tss" else matched_annots_global_tts

        if use_peaks:
            # Peak-based recall: all isoforms (incl. single-exon) vs all peaks.
            # Uses interval-edge matching consistent with the precision dedup loop.
            all_matched_peaks: Set[Tuple[str, str, int, int]] = set()
            for iso in isoforms:
                pos = iso[end_type]
                peak_ivs = current_peaks.get((iso["chrom"], iso["strand"]), [])
                m = _nearest_annot(pos, peak_ivs, window)
                if m is not None:
                    all_matched_peaks.add((iso["chrom"], iso["strand"], m[0], m[1]))
            total_peaks = sum(len(v) for v in current_peaks.values())
            recall = len(all_matched_peaks) / total_peaks if total_peaks > 0 else None
            n_annot_matched = len(all_matched_peaks)
            n_annot_total = total_peaks
        else:
            # GTF-based: JC-matched recall
            relevant = set()
            for chrom in relevant_annot_ends:
                for pos in relevant_annot_ends[chrom].get(end_type, set()):
                    relevant.add((chrom, pos))
            recall = len(matched_annots_global) / len(relevant) if relevant else None
            n_annot_matched = len(matched_annots_global)
            n_annot_total = len(relevant)

        # F1
        f1 = None
        if dedup_precision is not None and recall is not None and (dedup_precision + recall) > 0:
            f1 = 2 * dedup_precision * recall / (dedup_precision + recall)

        results[f"{end_label}_n_isoform_ends"] = total_ends
        results[f"{end_label}_dedup_tp"] = total_tp
        results[f"{end_label}_redundant_calls"] = total_red
        results[f"{end_label}_dedup_precision"] = dedup_precision
        results[f"{end_label}_naive_precision"] = naive_precision
        results[f"{end_label}_recall"] = recall
        results[f"{end_label}_dedup_f1"] = f1
        results[f"{end_label}_n_annot_matched"] = n_annot_matched
        results[f"{end_label}_n_annot_total"] = n_annot_total

    n_se_groups = sum(1 for k in jc_groups if len(k) == 4 and k[2] == "SE")
    n_se_isos   = sum(len(v) for k, v in jc_groups.items() if len(k) == 4 and k[2] == "SE")
    results["n_jc_groups"]   = len(jc_groups)
    results["n_se_groups"]   = n_se_groups
    results["n_single_exon"] = n_se_isos
    results["n_isoforms_total"] = len(isoforms)

    # ── Paired-end dedup precision ──────────────────────────────────────
    results["paired_n_isoforms"] = total_paired_isoforms
    results["paired_unique_pairs"] = total_paired_unique
    results["paired_both_hit"] = total_paired_both_hit
    results["paired_dedup_precision"] = (
        total_paired_both_hit / total_paired_isoforms
        if total_paired_isoforms > 0 else None
    )
    results["paired_naive_precision"] = (
        total_paired_unique / total_paired_isoforms
        if total_paired_isoforms > 0 else None
    )

    return results


# ── Output ──────────────────────────────────────────────────────────────────

def write_summary_tsv(results: dict, outpath: Path, mode_label: str = ""):
    # Column names use the shared vocabulary of precision_recall_plot.py:
    #   transcriptome_mode, 5prime_precision, 3prime_precision, 5prime_recall, 3prime_recall
    # The "dedup" prefix is dropped here since this IS the dedup metric; raw (non-dedup)
    # values are kept under the *_naive_precision names for reference.
    fields = [
        "transcriptome_mode", "n_isoforms_total", "n_jc_groups", "n_single_exon",
        "5prime_n_isoform_ends", "5prime_dedup_tp", "5prime_redundant_calls",
        "5prime_precision", "5prime_naive_precision", "5prime_recall", "5prime_f1",
        "5prime_n_annot_matched", "5prime_n_annot_total",
        "3prime_n_isoform_ends", "3prime_dedup_tp", "3prime_redundant_calls",
        "3prime_precision", "3prime_naive_precision", "3prime_recall", "3prime_f1",
        "3prime_n_annot_matched", "3prime_n_annot_total",
        "paired_n_isoforms", "paired_unique_pairs", "paired_both_hit",
        "paired_dedup_precision", "paired_naive_precision",
    ]
    # Internal results dict still uses _dedup_ names; remap for output
    _remap = {
        "5prime_dedup_precision": "5prime_precision",
        "5prime_dedup_f1":        "5prime_f1",
        "3prime_dedup_precision": "3prime_precision",
        "3prime_dedup_f1":        "3prime_f1",
    }
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        row = {"transcriptome_mode": mode_label}
        for col in fields[1:]:
            # Look up via remapped key if needed
            src_key = {v: k for k, v in _remap.items()}.get(col, col)
            v = results.get(src_key)
            if isinstance(v, float):
                row[col] = f"{v:.6f}"
            elif v is not None:
                row[col] = v
            else:
                row[col] = ""
        w.writerow(row)
    log.info(f"  → {outpath.name}")


def write_per_jc_tsv(per_jc: list, outpath: Path):
    if not per_jc:
        return
    fields = ["chrom", "strand", "jc_hash", "n_isoforms", "n_junctions",
              "5prime_unique_pos", "5prime_dedup_tp", "5prime_redundant",
              "3prime_unique_pos", "3prime_dedup_tp", "3prime_redundant",
              "paired_unique", "paired_both_hit", "paired_redundant"]
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for row in per_jc:
            w.writerow(row)
    log.info(f"  → {outpath.name}")


# ── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Junction-chain-deduplicated end precision/recall."
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--isoforms-bed",
                        help="Isoforms BED12 (FLAIR output)")
    input_group.add_argument("--isoforms-gtf",
                        help="Isoforms GTF (Bambu/IsoQuant/StringTie2/etc.)")
    parser.add_argument("--gtf", required=True,
                        help="Reference annotation GTF (full or partitioned)")
    parser.add_argument("--peaks-5prime", default=None,
                        help="BED6 file of CAGE peaks for TSS evaluation")
    parser.add_argument("--peaks-3prime", default=None,
                        help="BED6 file of dRNA peaks for TTS evaluation")
    parser.add_argument("--window", type=int, default=50,
                        help="Max distance (bp) for end matching (default: 50)")
    parser.add_argument("--region", default=None,
                        help="Restrict annotation to region (e.g. chr7:116000000-117000000)")
    parser.add_argument("--mode", default="",
                        help="Label for the transcriptome mode (for summary TSV)")
    parser.add_argument("--outdir", required=True,
                        help="Output directory")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Parse region filter
    region_chrom, region_start, region_end = None, None, None
    if args.region:
        if ":" in args.region:
            region_chrom, coords = args.region.split(":", 1)
            region_start, region_end = [int(x) for x in coords.split("-")]
        else:
            region_chrom = args.region

    log.info(f"Parsing GTF ends (region={args.region or 'all'})...")
    annotated_ends = parse_gtf_ends(args.gtf, region_chrom, region_start, region_end)
    n_tss = sum(len(v["tss"]) for v in annotated_ends.values())
    n_tts = sum(len(v["tts"]) for v in annotated_ends.values())
    log.info(f"  {n_tss} annotated TSS, {n_tts} annotated TTS (region total)")

    log.info("Parsing GTF transcripts for JC-matched recall...")
    annot_transcripts = parse_gtf_transcripts(args.gtf, region_chrom, region_start, region_end)
    log.info(f"  {len(annot_transcripts)} annotation transcripts")

    # Parse orthogonal peaks if provided
    peaks_5prime = None
    peaks_3prime = None
    if args.peaks_5prime:
        log.info(f"Parsing 5' peaks (CAGE): {args.peaks_5prime}")
        peaks_5prime = parse_peaks_bed(args.peaks_5prime, region_chrom, region_start, region_end)
        n5 = sum(len(v) for v in peaks_5prime.values())
        log.info(f"  {n5} peaks in region")
    if args.peaks_3prime:
        log.info(f"Parsing 3' peaks (dRNA): {args.peaks_3prime}")
        peaks_3prime = parse_peaks_bed(args.peaks_3prime, region_chrom, region_start, region_end)
        n3 = sum(len(v) for v in peaks_3prime.values())
        log.info(f"  {n3} peaks in region")

    if args.isoforms_bed:
        log.info("Parsing isoforms BED...")
        isoforms = parse_isoforms_bed(args.isoforms_bed)
    else:
        log.info("Parsing isoforms GTF...")
        isoforms = parse_isoforms_gtf(args.isoforms_gtf, region_chrom, region_start, region_end)
    log.info(f"  {len(isoforms)} isoforms")

    log.info("Computing JC-deduplicated precision/recall...")
    results = compute_jc_deduplicated_precision_recall(
        isoforms, annotated_ends, annot_transcripts, window=args.window,
        peaks_5prime=peaks_5prime, peaks_3prime=peaks_3prime,
    )

    log.info(f"Results: {results['n_jc_groups']} JC groups, "
             f"{results['n_isoforms_total']} isoforms")
    for end in ("5prime", "3prime"):
        dp = results.get(f"{end}_dedup_precision")
        np_ = results.get(f"{end}_naive_precision")
        r = results.get(f"{end}_recall")
        f1 = results.get(f"{end}_dedup_f1")
        log.info(f"  {end}: dedup_P={dp:.4f}, naive_P={np_:.4f}, R={r:.4f}, F1={f1:.4f}"
                 if all(v is not None for v in (dp, np_, r, f1)) else f"  {end}: insufficient data")

    write_summary_tsv(results, outdir / "precision_recall_summary.tsv",
                      mode_label=args.mode)
    write_per_jc_tsv(results["per_jc"], outdir / "per_junction_chain.tsv")
    log.info("Done.")


if __name__ == "__main__":
    main()
