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

# Shared, contract-compliant matcher (summit dedup + one-to-one greedy
# nearest).  Used for peak-based RECALL so a fat cluster of called ends
# near one truth peak yields exactly one matched summit, not one per call
# (the R3 recall inflator).  Same module the logreg fit + pileup PR use.
try:
    from peak_matching import (
        build_summit_index, match_positions_one_to_one,
        is_canonical_chrom, _match_distance, _peak_lo, _peak_hi)
except ImportError:
    from evaluation.peak_matching import (
        build_summit_index, match_positions_one_to_one,
        is_canonical_chrom, _match_distance, _peak_lo, _peak_hi)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s")
log = logging.getLogger(__name__)

Region = Tuple[str, Optional[int], Optional[int]]


def parse_region_spec(region: str) -> Region:
    """Parse chr or chr:start-end into a normalized region tuple."""
    if ":" not in region:
        return region, None, None
    chrom, coords = region.split(":", 1)
    start_s, end_s = coords.split("-", 1)
    start, end = int(start_s), int(end_s)
    if start > end:
        start, end = end, start
    return chrom, start, end


def parse_region_values(region_values: Optional[List[str]]) -> List[Region]:
    """Parse optional argparse --region values."""
    if not region_values:
        return []
    return [parse_region_spec(r) for r in region_values if r]


def _overlaps_regions(chrom: str, start: int, end: int,
                      regions: Optional[List[Region]] = None) -> bool:
    """Return True when an interval overlaps any requested region."""
    if not regions:
        return True
    for r_chrom, r_start, r_end in regions:
        if chrom != r_chrom:
            continue
        if r_start is None or r_end is None:
            return True
        if not (end < r_start or start > r_end):
            return True
    return False


def parse_feature_attributes(attr_string: str) -> dict:
    """Parse GTF/GFF3 attributes into a dict."""
    attrs = {}
    for part in attr_string.strip().rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        if ' "' in part:
            key, value = part.split(' "', 1)
            attrs[key.strip()] = value.rstrip('"')
        elif "=" in part:
            key, value = part.split("=", 1)
            attrs[key.strip()] = value.strip().strip('"')
        elif " " in part:
            key, value = part.split(" ", 1)
            attrs[key.strip()] = value.strip().strip('"')
    return attrs


def _first_attr_value(attrs: dict, key: str) -> Optional[str]:
    value = attrs.get(key)
    if not value:
        return None
    return str(value).split(",", 1)[0]


def _transcript_id_from_attrs(attrs: dict, feature: str) -> Optional[str]:
    """Get transcript ID from GTF or GFF3 attributes."""
    if attrs.get("transcript_id"):
        return attrs["transcript_id"]
    if feature in ("transcript", "mRNA") and attrs.get("ID"):
        return attrs["ID"]
    parent = _first_attr_value(attrs, "Parent")
    if parent:
        return parent
    return None


# ── Parse inputs ────────────────────────────────────────────────────────────

def parse_gtf_ends(gtf_path: str, regions: Optional[List[Region]] = None
                   ) -> Dict[str, Dict[str, List[int]]]:
    """Extract annotated TSS/TTS positions from GTF/GFF, optionally filtered to regions."""
    ends: Dict[str, Dict[str, set]] = defaultdict(lambda: {"tss": set(), "tts": set()})
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] not in ("transcript", "mRNA"):
                continue
            chrom = cols[0]
            start = int(cols[3]) - 1  # GTF is 1-based
            end = int(cols[4])
            strand = cols[6]
            if not _overlaps_regions(chrom, start, end, regions):
                continue
            if strand == "+":
                ends[chrom]["tss"].add(start)
                ends[chrom]["tts"].add(end)
            else:
                ends[chrom]["tss"].add(end)
                ends[chrom]["tts"].add(start)
    return {c: {"tss": sorted(ends[c]["tss"]), "tts": sorted(ends[c]["tts"])}
            for c in ends}


def parse_gtf_transcripts(gtf_path: str, regions: Optional[List[Region]] = None
                          ) -> List[dict]:
    """Parse GTF/GFF transcript+exon features to get junction chains and ends."""
    tx_exons: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    tx_info: Dict[str, Tuple[str, str]] = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom = cols[0]
            feature = cols[2]
            if feature not in ("transcript", "mRNA", "exon"):
                continue
            start = int(cols[3]) - 1  # GTF is 1-based
            end = int(cols[4])
            strand = cols[6]
            if not _overlaps_regions(chrom, start, end, regions):
                continue
            attrs = parse_feature_attributes(cols[8])
            tx_id = _transcript_id_from_attrs(attrs, feature)
            if tx_id is None:
                continue
            if feature in ("transcript", "mRNA"):
                tx_info[tx_id] = (chrom, strand)
            elif feature == "exon":
                tx_info.setdefault(tx_id, (chrom, strand))
                tx_exons[tx_id].append((start, end))
    transcripts = []
    for tx_id, exons in tx_exons.items():
        if tx_id not in tx_info:
            continue
        chrom, strand = tx_info[tx_id]
        exons_sorted = sorted(exons)
        tx_start = min(start for start, _ in exons_sorted)
        tx_end = max(end for _, end in exons_sorted)
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


def parse_isoforms_bed(bed_path: str, regions: Optional[List[Region]] = None) -> List[dict]:
    """Parse BED12 isoforms into list of dicts with junction chains."""
    isoforms = parse_bed12(bed_path)
    if not regions:
        return isoforms
    return [iso for iso in isoforms
            if _overlaps_regions(iso["chrom"], iso["start"], iso["end"], regions)]


def parse_isoforms_gtf(gtf_path: str, region_chrom: str = None,
                       region_start: int = None, region_end: int = None,
                       regions: Optional[List[Region]] = None
                       ) -> List[dict]:
    """Parse assembler GTF (Bambu/IsoQuant/StringTie2/etc.) into same format as parse_bed12.

    Reuses parse_gtf_transcripts — already reads exon features and builds
    junction chains.  Returns a list of transcript dicts compatible with
    compute_jc_deduplicated_precision_recall.
    """
    if regions is None and isinstance(region_chrom, list):
        regions = region_chrom
    elif regions is None and region_chrom:
        regions = [(region_chrom, region_start, region_end)]
    return parse_gtf_transcripts(gtf_path, regions)


def load_supported_ids(counts_path: Optional[str], min_support: int) -> Optional[Set[str]]:
    """Return transcript IDs with count >= min_support, or None if unavailable.

    Supports IsoQuant two-column counts, Bambu TXNAME/GENEID/sample counts,
    and FLAIR isoform count tables by treating column 1 as the transcript ID
    and summing numeric fields to its right.
    """
    if not counts_path:
        return None
    path = Path(counts_path)
    if not path.exists():
        return None
    supported: Set[str] = set()
    n_parsed = 0
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            numeric_values = []
            for value in parts[1:]:
                try:
                    numeric_values.append(float(value))
                except ValueError:
                    continue
            if not numeric_values:
                continue
            n_parsed += 1
            if sum(numeric_values) >= min_support:
                supported.add(parts[0])
    return supported if n_parsed else None


def _matches_supported_id(tx_id: str, supported_ids: Set[str]) -> bool:
    """Match direct IDs plus the FLAIR `transcript_gene` naming convention."""
    if tx_id in supported_ids:
        return True
    if "_" in tx_id and tx_id.split("_", 1)[0] in supported_ids:
        return True
    return False


def filter_isoforms_by_counts(
    isoforms: List[dict],
    supported_ids: Optional[Set[str]],
) -> List[dict]:
    """Filter parsed BED/GTF isoforms to supported transcript IDs."""
    if supported_ids is None:
        return isoforms
    filtered = []
    for iso in isoforms:
        tx_id = str(iso.get("tx_id") or iso.get("name") or "")
        if tx_id and _matches_supported_id(tx_id, supported_ids):
            filtered.append(iso)
    return filtered


def parse_peaks_bed(path: str, regions: Optional[List[Region]] = None,
                    canonical_only: bool = True
                    ) -> Dict[Tuple[str, str], List[Tuple[int, int]]]:
    """Parse BED6 peaks into {(chrom, strand): sorted (start, end) intervals}.

    Intervals use BED half-open convention [start, end).  Sorted by start so
    binary search can locate nearby peaks efficiently.

    When `canonical_only` is True (default) peaks on non-canonical contigs
    (chrM and unplaced/alt scaffolds — anything that is not chr1-22/X/Y) are
    dropped, so the recall denominator and precision are computed only over
    the canonical human chromosomes.  A `--region` allow-list (regions) is
    applied on top.
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
            if canonical_only and not is_canonical_chrom(chrom):
                continue
            if not _overlaps_regions(chrom, start, end, regions):
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


def _interval_centre(start: int, end: int) -> int:
    """Summit of a peak with no signal track — the interval centre.

    Mirrors peak_matching._interval_summit's no-signal branch so the
    paired-dedup path uses the SAME summit definition as the one-to-one /
    single-end paths.
    """
    return (int(start) + int(end)) // 2


def _nearest_annot(pos: int, sorted_intervals: List[Tuple[int, int]], window: int
                   ) -> Optional[Tuple[int, int, int]]:
    """Return the matching peak (summit, start, end) under the OR predicate, or None.

    THE contract predicate (shared with peak_matching.match_positions_one_to_one):

        call c matches peak p  iff  (|c - p.summit| <= W)  OR  (p.start <= c < p.end)

    A call is credited if it is within ±W of the peak SUMMIT (interval centre,
    no signal track) *or* falls INSIDE the peak interval [start, end).  The
    interior case gives broad-promoter credit; for peaks with width <= 2W the
    interval is a subset of summit±W so the OR adds nothing (backward-compat).

    Among competing peaks the one with the smallest match-distance wins
    (0 when interior, else |c - summit|), tie-broken by |c - summit|, so an
    interior call still claims the peak it sits inside.  Returns a hashable
    (summit, start, end) identity so paired-end dedup collapses isoforms that
    hit the SAME peak.

    A peak's interval may extend beyond ±W from its summit, so the summit
    bisection bound is widened by the max interval half-extent to never miss
    an interior match.
    """
    if not sorted_intervals:
        return None
    # sorted_intervals is sorted by start (parse_peaks_bed); compute summits
    # and the max half-extent so we can bound the scan and never miss an
    # interior-but-summit-far peak.
    starts = [iv[0] for iv in sorted_intervals]
    max_left = 0   # max (summit - start)
    max_right = 0  # max (end - summit)
    for s, e in sorted_intervals:
        c = _interval_centre(s, e)
        if c - s > max_left:
            max_left = c - s
        if e - c > max_right:
            max_right = e - c
    lo_bound = pos - max(window, max_right)
    hi_bound = pos + max(window, max_left)
    best_key = None
    best_match_dist = None
    best_summit_dist = None
    # Walk intervals whose start could still bring the peak into range.  Start
    # a little left of lo_bound's start-position via bisect on starts.
    i = bisect_left(starts, lo_bound - max_left - 1)
    n = len(sorted_intervals)
    while i < n:
        s, e = sorted_intervals[i]
        if s > hi_bound:
            break
        summit = _interval_centre(s, e)
        summit_dist = abs(pos - summit)
        inside = s <= pos < e
        if inside or summit_dist <= window:
            match_dist = 0 if inside else summit_dist
            if (best_match_dist is None
                    or match_dist < best_match_dist
                    or (match_dist == best_match_dist
                        and summit_dist < best_summit_dist)):
                best_match_dist = match_dist
                best_summit_dist = summit_dist
                best_key = (summit, s, e)
        i += 1
    return best_key


def _summit_index_from_intervals(
    peaks: Dict[Tuple[str, str], List[Tuple[int, int]]]
) -> Dict[Tuple[str, str], List[int]]:
    """Collapse an interval peak dict to a (chrom,strand)->summit index.

    One summit per peak interval (peak centre; no signal bedGraph wired
    in).  Used to give peak-based recall the contract's one-to-one,
    summit-deduped accounting.
    """
    rows: List[Tuple[str, int, int, str]] = []
    for (chrom, strand), ivs in peaks.items():
        for s, e in ivs:
            rows.append((chrom, s, e, strand))
    return build_summit_index(rows)


def _one_to_one_match_counts(
    isoforms: List[dict],
    summit_index: Dict[Tuple[str, str], List[int]],
    end_type: str,
    window: int,
) -> Tuple[int, int, int]:
    """One-to-one (greedy nearest) peak matching — shared TP for P AND R.

    Buckets all isoform ends by (chrom, strand), matches each bucket to its
    summits with the shared greedy-nearest matcher (each summit consumed at
    most once), and returns ``(tp, n_calls, n_summits)``.

    Because matching is one-to-one, the number of matched calls equals the
    number of matched summits, so a single TP count feeds BOTH sides of the
    F1 (the W1 unification, mirroring ted_core.one_to_one_match_counts):

        precision = tp / n_calls      recall = tp / n_summits

    A fat cluster of called ends near one summit therefore contributes
    exactly ONE true positive on both sides — removing the many-to-one
    precision inflator that ``_nearest_annot`` produced (each ref peak could
    be credited once per nearby call).
    """
    buckets: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    n_calls = 0
    for iso in isoforms:
        buckets[(iso["chrom"], iso["strand"])].append(int(iso[end_type]))
        n_calls += 1
    matched = 0
    for key, positions in buckets.items():
        summits = summit_index.get(key, [])
        if not summits:
            continue
        assign = match_positions_one_to_one(positions, summits, window)
        matched += sum(1 for j in assign if j >= 0)
    total = sum(len(v) for v in summit_index.values())
    return matched, n_calls, total


def _many_to_one_match_counts(
    isoforms: List[dict],
    summit_index: Dict[Tuple[str, str], List[int]],
    end_type: str,
    window: int,
) -> Tuple[int, int, int, int]:
    """MANY-TO-ONE 'near a peak' counts for SINGLE-END precision/recall.

    This is the DEFINITION FOR SINGLE 5'/3' ENDS (hdheath 2026-06-15): a
    promoter (or cleavage region) that produces several isoforms is EXPECTED,
    so every isoform end within `window` of a peak is credited — duplication
    is allowed. Uniqueness / dedup is enforced ONLY in the PAIRED metric
    (paired_dedup_precision), never here.

    Mirrors bed_utils.near_match_counts (the function ted_core.py /
    combined_evaluation.tsv use) but on this module's summit_index/isoform
    structures, so the plotted single-end numbers and combined_evaluation.tsv
    use the SAME rule and cannot diverge.

      ends_near       = # isoform ends within `window` of ANY summit (precision num)
      peaks_recovered = # summits with ANY isoform end within `window` (recall num)
    Returns (ends_near, n_calls, peaks_recovered, n_summits).
    """
    ends_by_key: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    n_calls = 0
    for iso in isoforms:
        ends_by_key[(iso["chrom"], iso["strand"])].append(int(iso[end_type]))
        n_calls += 1
    n_summits = sum(len(v) for v in summit_index.values())

    def _peak_half_extents(peaks) -> Tuple[int, int]:
        """(max summit-to-start, max end-to-summit) over a peak bucket."""
        ml = mr = 0
        for p in peaks:
            s = int(p)
            ml = max(ml, s - _peak_lo(p))
            mr = max(mr, _peak_hi(p) - s)
        return ml, mr

    def _end_near_any_peak(e: int, sorted_peaks, ml: int, mr: int) -> bool:
        """CORRECTED predicate (2026-06-21): e within W of ANY part of the peak
        INTERVAL, i.e. peak_matching._match_distance(e, peak) <= window.

        This is the documented eval-contract predicate (supported zone
        [lo-W, (hi-1)+W]); for a wide peak it ALSO credits a call within W of an
        EDGE, not only within W of the summit.  For width<=2W peaks (e.g. the
        width-1 dRNA cleavage peaks) the interval is a subset of summit±W so the
        result is identical to the previous summit-OR-interior predicate — only
        wide CAGE promoters (median 45 bp, 14% >100 bp on chr22) shift, and they
        shift toward the contract.  This is the SAME primitive stratified_pr.py /
        end_dist / peak_matching.match_positions_one_to_one use, so the single-
        end headline 5prime/3prime precision matches across all eval scripts.
        """
        # Scan bound: a peak's summit can sit up to (W + mr) left of e (e past
        # the right edge) or (W + ml) right of e (e before the left edge) and
        # still match within W of an edge.  Mirrors
        # peak_matching.match_positions_one_to_one.
        lo_bound = e - (window + mr)
        hi_bound = e + (window + ml)
        j = bisect_left(sorted_peaks, lo_bound)
        while j < len(sorted_peaks):
            p = sorted_peaks[j]
            s = int(p)
            if s > hi_bound:
                break
            if _match_distance(e, p) <= window:
                return True
            j += 1
        return False

    # ── precision numerator: # ends near ANY peak (interval-OR-summit) ──
    ends_near = 0
    for key, ends in ends_by_key.items():
        peaks = summit_index.get(key, [])
        if not peaks:
            continue
        sp = sorted(peaks)
        ml, mr = _peak_half_extents(sp)
        for e in ends:
            if _end_near_any_peak(e, sp, ml, mr):
                ends_near += 1

    # ── recall numerator: # peaks recovered by ANY end (SAME corrected predicate) ──
    peaks_recovered = 0
    for key, peaks in summit_index.items():
        ends = ends_by_key.get(key, [])
        if not ends:
            continue
        es = sorted(ends)
        for p in peaks:
            s = int(p)
            lo = _peak_lo(p)
            hi = _peak_hi(p)
            # A peak is recovered if any end is within W of ANY part of the
            # interval, i.e. _match_distance(e, p) <= window (supported zone
            # [lo-W, (hi-1)+W]).  Scan ends across that zone.
            scan_lo = min(lo, s) - window
            scan_hi = max(hi, s + 1) + window
            i = bisect_left(es, scan_lo)
            recovered = False
            while i < len(es):
                e = es[i]
                if e >= scan_hi:
                    break
                if _match_distance(e, p) <= window:
                    recovered = True
                    break
                i += 1
            if recovered:
                peaks_recovered += 1
    return ends_near, n_calls, peaks_recovered, n_summits


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
        paired dedup TP = unique pairs with BOTH peaks matched
      This means same-TSS-different-TTS isoforms both count as 5' TPs,
      and same-TTS-different-TSS both count as 3' TPs.

      Denominator for ALL dedup precisions = total isoforms emitted by the
      method (n_isoforms_total).  An isoform that fails (off-peak end, or
      same-peak-pair as a JC sibling) shows up in the denominator but not
      the numerator.  This rewards methods for emitting many distinct
      peak-pair hits while penalizing both off-peak emissions and same-pair
      duplicates within a JC.

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

    # Clean up temporary marker
    for iso in isoforms:
        iso.pop("_se", None)

    results = {"per_jc": []}

    # GTF-based recall now uses every distinct annotated TSS/TTS as the
    # denominator (see the recall block below).  The previous JC-filtered
    # denominator made annotation-passthrough tools score ~100% trivially.

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
    total_paired_both_hit_raw = 0
    total_5p_raw_tp = 0
    total_3p_raw_tp = 0
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

        unique_pairs = set(end_pairs)
        n_members = len(members)

        n_5p_tp = sum(1 for p in unique_pairs if p[0] is not None)
        n_3p_tp = sum(1 for p in unique_pairs if p[1] is not None)
        n_both_hit = sum(1 for p in unique_pairs if p[0] is not None and p[1] is not None)
        n_5p_raw_tp = sum(1 for p in end_pairs if p[0] is not None)
        n_3p_raw_tp = sum(1 for p in end_pairs if p[1] is not None)
        n_both_hit_raw = sum(1 for p in end_pairs if p[0] is not None and p[1] is not None)

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
        total_paired_both_hit_raw += n_both_hit_raw
        total_5p_raw_tp += n_5p_raw_tp
        total_3p_raw_tp += n_3p_raw_tp

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
            "paired_both_hit_raw": n_both_hit_raw,
            "paired_redundant": n_members - n_unique,
        })

    # ── Per-end precision ───────────────────────────────────────────────
    for end_label, total_ends, total_tp, total_raw_tp, total_red in [
        ("5prime", total_5p_ends, total_5p_dedup_tp, total_5p_raw_tp, total_5p_redundant),
        ("3prime", total_3p_ends, total_3p_dedup_tp, total_3p_raw_tp, total_3p_redundant),
    ]:
        # JC-deduplicated (many-to-one via _nearest_annot) and naive (raw)
        # precision are kept as DESCRIPTIVE columns.  In the peak-based path
        # the HEADLINE precision below is overridden with the one-to-one value
        # so it shares the same TP as recall (W1; see ted_core.py ~211-233).
        dedup_precision = total_tp / total_ends if total_ends > 0 else None
        naive_precision = total_raw_tp / total_ends if total_ends > 0 else None

        # ── Recall ──────────────────────────────────────────────────────
        end_type = "tss" if end_label == "5prime" else "tts"
        current_peaks = peaks_5prime if end_type == "tss" else peaks_3prime
        use_peaks = current_peaks is not None

        if use_peaks:
            # SINGLE-END headline = MANY-TO-ONE (hdheath 2026-06-15).
            # An isoform end is credited if it's within `window` of ANY peak,
            # and a peak is recovered if ANY end is near it — duplication is
            # ALLOWED, because one promoter/cleavage region legitimately yields
            # several isoforms. This is the SAME rule as ted_core.near_match_counts
            # (combined_evaluation.tsv), so the plotted single-end numbers match
            # those tables. The one-to-one (unique-consumption) logic belongs
            # ONLY to the PAIRED metric (paired_dedup_precision) below — it is
            # kept here purely as a labelled diagnostic (*_onetoone_*).
            summit_index = _summit_index_from_intervals(current_peaks)
            ends_near, n_calls, peaks_rec, total_peaks = _many_to_one_match_counts(
                isoforms, summit_index, end_type, window)
            precision = ends_near / n_calls if n_calls > 0 else None
            recall = peaks_rec / total_peaks if total_peaks > 0 else None
            dedup_precision = precision  # -> headline {end}_precision (many-to-one)
            n_annot_matched = peaks_rec
            n_annot_total = total_peaks
            # one-to-one kept as a diagnostic only (NOT plotted as the headline)
            oto_tp, _oto_calls, _oto_summits = _one_to_one_match_counts(
                isoforms, summit_index, end_type, window)
            results[f"{end_label}_onetoone_precision"] = (
                oto_tp / n_calls if n_calls > 0 else None)
            results[f"{end_label}_onetoone_recall"] = (
                oto_tp / total_peaks if total_peaks > 0 else None)
        else:
            # GTF-based recall: every distinct annotated TSS/TTS in the
            # reference GTF is a target — same denominator for every tool,
            # matching the orthogonal-peak pattern above.  Previous behavior
            # filtered the denominator to "annotated transcripts whose JC was
            # also predicted by this tool", which made the denominator
            # tool-dependent and trivially gave annotation-passthrough tools
            # ~100% recall (their predicted JCs are reference JCs and their
            # predicted ends are reference ends).
            all_matched_annots: Set[Tuple[str, int]] = set()
            for iso in isoforms:
                pos = iso[end_type]
                annot_ref = annotated_ends.get(iso["chrom"], {}).get(end_type, [])
                m = _nearest_pos(pos, annot_ref, window)
                if m is not None:
                    all_matched_annots.add((iso["chrom"], m))
            total_annot_ends = sum(
                len(annotated_ends[c].get(end_type, []))
                for c in annotated_ends
            )
            recall = (len(all_matched_annots) / total_annot_ends
                      if total_annot_ends > 0 else None)
            n_annot_matched = len(all_matched_annots)
            n_annot_total = total_annot_ends

        # F1
        f1 = None
        if dedup_precision is not None and recall is not None and (dedup_precision + recall) > 0:
            f1 = 2 * dedup_precision * recall / (dedup_precision + recall)

        results[f"{end_label}_n_isoform_ends"] = total_ends
        results[f"{end_label}_dedup_tp"] = total_tp
        results[f"{end_label}_raw_tp"] = total_raw_tp
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
    results["paired_both_hit_raw"] = total_paired_both_hit_raw
    # paired_dedup_precision:
    #   numerator   = JC-unique (TSS-peak, TTS-peak) pairs where BOTH ends hit
    #                 a peak within tolerance.  An isoform is "redundant" with
    #                 a JC sibling only when it lands on the SAME (TSS-peak,
    #                 TTS-peak) combination — single-end peak collisions or
    #                 off-peak emissions do not collapse.
    #   denominator = total isoforms emitted by the method.
    # This penalizes both genuinely-wrong calls (off-peak) AND same-peak-pair
    # JC sibling duplicates, while crediting every distinct peak-pair hit.
    # paired_naive_precision: both-hit isoforms (raw, no dedup) / total isoforms
    results["paired_dedup_precision"] = (
        total_paired_both_hit / total_paired_isoforms
        if total_paired_isoforms > 0 else None
    )
    results["paired_naive_precision"] = (
        total_paired_both_hit_raw / total_paired_isoforms
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
        "paired_both_hit_raw", "paired_dedup_precision", "paired_naive_precision",
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
    # Always write the file (header-only when empty) so downstream consumers —
    # e.g. the Nextflow module's unconditional `mv per_junction_chain.tsv` — never
    # fail on tools that yield 0 isoforms in the region (e.g. isoseq on chr22).
    fields = ["chrom", "strand", "jc_hash", "n_isoforms", "n_junctions",
              "5prime_unique_pos", "5prime_dedup_tp", "5prime_redundant",
              "3prime_unique_pos", "3prime_dedup_tp", "3prime_redundant",
              "paired_unique", "paired_both_hit", "paired_both_hit_raw",
              "paired_redundant"]
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
    parser.add_argument("--region", nargs="+", default=None,
                        help=("Restrict evaluation to one or more regions "
                              "(e.g. chr7:116000000-117000000 chr11:64000000-69000000)"))
    parser.add_argument("--mode", default="",
                        help="Label for the transcriptome mode (for summary TSV)")
    parser.add_argument("--counts", default=None,
                        help="Optional per-transcript counts table for support filtering")
    parser.add_argument("--min-support", type=int, default=1,
                        help="Minimum transcript count retained when --counts is provided")
    parser.add_argument("--outdir", required=True,
                        help="Output directory")
    parser.add_argument("--all-contigs", action="store_true",
                        help="Keep peaks/isoforms on ALL contigs (chrM + "
                             "GL*/KI* scaffolds).  Default restricts matching "
                             "to the canonical human chromosomes chr1-22,X,Y.")
    args = parser.parse_args()
    canonical_only = not args.all_contigs

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    regions = parse_region_values(args.region)
    region_label = " ".join(args.region) if args.region else "all"

    log.info(f"Parsing GTF ends (region={region_label})...")
    annotated_ends = parse_gtf_ends(args.gtf, regions)
    n_tss = sum(len(v["tss"]) for v in annotated_ends.values())
    n_tts = sum(len(v["tts"]) for v in annotated_ends.values())
    log.info(f"  {n_tss} annotated TSS, {n_tts} annotated TTS (region total)")

    log.info("Parsing GTF transcripts for JC-matched recall...")
    annot_transcripts = parse_gtf_transcripts(args.gtf, regions)
    log.info(f"  {len(annot_transcripts)} annotation transcripts")

    # Parse orthogonal peaks if provided
    peaks_5prime = None
    peaks_3prime = None
    if args.peaks_5prime:
        log.info(f"Parsing 5' peaks (CAGE): {args.peaks_5prime}")
        peaks_5prime = parse_peaks_bed(args.peaks_5prime, regions,
                                       canonical_only=canonical_only)
        n5 = sum(len(v) for v in peaks_5prime.values())
        log.info(f"  {n5} peaks in region")
    if args.peaks_3prime:
        log.info(f"Parsing 3' peaks (dRNA): {args.peaks_3prime}")
        peaks_3prime = parse_peaks_bed(args.peaks_3prime, regions,
                                       canonical_only=canonical_only)
        n3 = sum(len(v) for v in peaks_3prime.values())
        log.info(f"  {n3} peaks in region")

    if args.isoforms_bed:
        log.info("Parsing isoforms BED...")
        isoforms = parse_isoforms_bed(args.isoforms_bed, regions)
    else:
        log.info("Parsing isoforms GTF...")
        isoforms = parse_gtf_transcripts(args.isoforms_gtf, regions)
    if canonical_only:
        n_before = len(isoforms)
        isoforms = [i for i in isoforms if is_canonical_chrom(i["chrom"])]
        if len(isoforms) != n_before:
            log.info("Canonical-chrom filter: kept %d/%d isoforms "
                     "(dropped chrM/scaffolds)", len(isoforms), n_before)
    supported_ids = load_supported_ids(args.counts, args.min_support)
    if supported_ids is not None:
        n_before = len(isoforms)
        isoforms = filter_isoforms_by_counts(isoforms, supported_ids)
        log.info(
            "Read-support filter: kept %d/%d isoforms "
            "(count >= %s in %s)",
            len(isoforms), n_before, args.min_support, Path(args.counts).name,
        )
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
