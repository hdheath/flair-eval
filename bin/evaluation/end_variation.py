"""
Alternative end-site analysis for shared splice junction chains.

Detects groups of multi-exon isoforms that share the same junction chain
but differ at transcript start/end positions (APU/APA-style variation).
"""

from bisect import bisect_right
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from signal_utils import group_by_junction_chain as _group_by_junction_chain, gene_from_name as _gene_from_name
except ImportError:
    from evaluation.signal_utils import group_by_junction_chain as _group_by_junction_chain, gene_from_name as _gene_from_name


def _load_peaks_by_chrom_strand(peaks_path: Path) -> Dict[Tuple[str, str], List[Tuple[int, int, str]]]:
    """Load BED6 peaks grouped by (chrom, strand).

    Returns dict: (chrom, strand) -> sorted list of (start, end, peak_id).
    Peak IDs are 'chrom_start_end'.
    """
    result: Dict[Tuple[str, str], List[Tuple[int, int, str]]] = defaultdict(list)
    with open(peaks_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 6:
                continue
            chrom, start, end, strand = cols[0], int(cols[1]), int(cols[2]), cols[5]
            pid = f"{chrom}_{start}_{end}"
            result[(chrom, strand)].append((start, end, pid))
    # Sort by start position for binary search
    for key in result:
        result[key].sort()
    return result


def _find_peak_for_position(
    pos: int,
    peaks: List[Tuple[int, int, str]],
    window: int,
) -> Optional[str]:
    """Find the closest peak (by gap distance) within *window* bp of *pos*.

    Uses bedtools-consistent gap distance: if pos falls inside [start, end),
    distance is 0; otherwise it's the gap to the nearest edge.

    Returns peak_id or None if no peak within window.
    """
    if not peaks:
        return None
    # Binary search for candidate peaks
    # peaks are sorted by start; we search for peaks whose start <= pos + window
    right_idx = bisect_right(peaks, (pos + window + 1,)) 
    best_pid: Optional[str] = None
    best_dist = window + 1
    # Check candidates from right_idx backward
    for i in range(right_idx - 1, -1, -1):
        p_start, p_end, pid = peaks[i]
        if p_start > pos + window:
            continue
        if p_end < pos - window:
            break  # No more candidates (sorted by start, so earlier peaks are even further left)
        # Compute gap distance
        if pos < p_start:
            dist = p_start - pos
        elif pos >= p_end:
            dist = pos - p_end + 1  # +1 because BED is half-open [start, end)
        else:
            dist = 0  # Inside peak
        if dist <= window and dist < best_dist:
            best_dist = dist
            best_pid = pid
    return best_pid


def summarize_junction_chain_end_variation(
    isoforms: Dict[str, dict],
    min_group_size: int = 1,
    peaks_5prime_path: Optional[Path] = None,
    peaks_3prime_path: Optional[Path] = None,
    window_5prime: int = 50,
    window_3prime: int = 50,
) -> Tuple[Dict[str, int], Dict[str, List[int]]]:
    """
    Summarize alternative-end events for shared junction chains.

    When peak files are provided, also computes how many "unique" end
    positions within a group actually map to distinct peaks.  Ends that
    map to the same peak are considered *redundant*.

    Args:
        isoforms: Output of parse_isoform_ends()
        min_group_size: Minimum isoforms required in a junction-chain group
        peaks_5prime_path: CAGE peaks BED6 (optional)
        peaks_3prime_path: dRNA peaks BED6 (optional)
        window_5prime: Window for TSS → peak matching
        window_3prime: Window for TTS → peak matching

    Returns:
        metrics:
            - junction_chain_groups_total
            - junction_chain_groups_alt_end
            - junction_chain_groups_apu_only
            - junction_chain_groups_apa_only
            - junction_chain_groups_both_ends
            - max_unique_tss_per_junction_chain
            - max_unique_tts_per_junction_chain
            - max_isoforms_per_junction_chain
            - (if peaks provided) unique_tss_after_dedup_total, redundant_tss_total, etc.
        distributions:
            - unique_tss_count_distribution
            - unique_tts_count_distribution
            - tss_span_distribution_bp
            - tts_span_distribution_bp
            - dedup_tss_count_distribution (unique peaks per group)
            - dedup_tts_count_distribution
            - redundant_tss_count_distribution (unique_positions - unique_peaks per group)
            - redundant_tts_count_distribution
    """
    groups = _group_by_junction_chain(isoforms)
    valid_groups = [txs for txs in groups.values() if len(txs) >= min_group_size]

    # Load peaks if provided
    peaks_5 = (
        _load_peaks_by_chrom_strand(peaks_5prime_path)
        if peaks_5prime_path and peaks_5prime_path.exists()
        else {}
    )
    peaks_3 = (
        _load_peaks_by_chrom_strand(peaks_3prime_path)
        if peaks_3prime_path and peaks_3prime_path.exists()
        else {}
    )

    tss_count_dist: List[int] = []
    tts_count_dist: List[int] = []
    tss_span_dist: List[int] = []
    tts_span_dist: List[int] = []
    # Peak-deduplicated distributions
    dedup_tss_count_dist: List[int] = []
    dedup_tts_count_dist: List[int] = []
    redundant_tss_count_dist: List[int] = []
    redundant_tts_count_dist: List[int] = []

    alt_end_groups = 0
    apu_only_groups = 0
    apa_only_groups = 0
    both_groups = 0

    max_unique_tss = 0
    max_unique_tts = 0
    max_group_size = 0

    total_unique_tss = 0
    total_unique_tts = 0
    total_dedup_tss = 0
    total_dedup_tts = 0

    # Retrieve group keys to access chrom/strand for peak lookup
    group_keys_and_txs = [
        (key, txs)
        for key, txs in groups.items()
        if len(txs) >= min_group_size
    ]

    for (chrom, strand, _), txs in group_keys_and_txs:
        tss_vals = sorted({int(t["tss"]) for t in txs})
        tts_vals = sorted({int(t["tts"]) for t in txs})

        n_tss = len(tss_vals)
        n_tts = len(tts_vals)
        tss_count_dist.append(n_tss)
        tts_count_dist.append(n_tts)

        max_unique_tss = max(max_unique_tss, n_tss)
        max_unique_tts = max(max_unique_tts, n_tts)
        max_group_size = max(max_group_size, len(txs))

        if n_tss > 1:
            tss_span_dist.append(tss_vals[-1] - tss_vals[0])
        if n_tts > 1:
            tts_span_dist.append(tts_vals[-1] - tts_vals[0])

        has_apu = n_tss > 1
        has_apa = n_tts > 1
        if has_apu or has_apa:
            alt_end_groups += 1
            if has_apu and has_apa:
                both_groups += 1
            elif has_apu:
                apu_only_groups += 1
            else:
                apa_only_groups += 1

        # ── Peak-dedup analysis ──
        # For each unique end position, find which peak it maps to.
        # Count distinct peaks hit → that's the "deduplicated" count.
        cs_key = (chrom, strand)

        if peaks_5:
            cs_peaks_5 = peaks_5.get(cs_key, [])
            tss_peak_ids = set()
            for pos in tss_vals:
                pid = _find_peak_for_position(pos, cs_peaks_5, window_5prime)
                if pid is not None:
                    tss_peak_ids.add(pid)
            n_dedup_tss = max(len(tss_peak_ids), 1) if n_tss > 0 else 0
            # Positions that didn't hit any peak still count as unique
            n_no_peak_tss = sum(
                1 for pos in tss_vals
                if _find_peak_for_position(pos, cs_peaks_5, window_5prime) is None
            )
            n_dedup_tss = len(tss_peak_ids) + n_no_peak_tss
            dedup_tss_count_dist.append(n_dedup_tss)
            redundant_tss_count_dist.append(max(0, n_tss - n_dedup_tss))
            total_unique_tss += n_tss
            total_dedup_tss += n_dedup_tss
        else:
            dedup_tss_count_dist.append(n_tss)
            redundant_tss_count_dist.append(0)
            total_unique_tss += n_tss
            total_dedup_tss += n_tss

        if peaks_3:
            cs_peaks_3 = peaks_3.get(cs_key, [])
            tts_peak_ids = set()
            for pos in tts_vals:
                pid = _find_peak_for_position(pos, cs_peaks_3, window_3prime)
                if pid is not None:
                    tts_peak_ids.add(pid)
            n_no_peak_tts = sum(
                1 for pos in tts_vals
                if _find_peak_for_position(pos, cs_peaks_3, window_3prime) is None
            )
            n_dedup_tts = len(tts_peak_ids) + n_no_peak_tts
            dedup_tts_count_dist.append(n_dedup_tts)
            redundant_tts_count_dist.append(max(0, n_tts - n_dedup_tts))
            total_unique_tts += n_tts
            total_dedup_tts += n_dedup_tts
        else:
            dedup_tts_count_dist.append(n_tts)
            redundant_tts_count_dist.append(0)
            total_unique_tts += n_tts
            total_dedup_tts += n_tts

    metrics = {
        "junction_chain_groups_total": len(valid_groups),
        "junction_chain_groups_alt_end": alt_end_groups,
        "junction_chain_groups_apu_only": apu_only_groups,
        "junction_chain_groups_apa_only": apa_only_groups,
        "junction_chain_groups_both_ends": both_groups,
        "max_unique_tss_per_junction_chain": max_unique_tss,
        "max_unique_tts_per_junction_chain": max_unique_tts,
        "max_isoforms_per_junction_chain": max_group_size,
        "unique_tss_total": total_unique_tss,
        "unique_tts_total": total_unique_tts,
        "dedup_tss_total": total_dedup_tss,
        "dedup_tts_total": total_dedup_tts,
        "redundant_tss_total": total_unique_tss - total_dedup_tss,
        "redundant_tts_total": total_unique_tts - total_dedup_tts,
    }

    distributions = {
        "unique_tss_count_distribution": tss_count_dist,
        "unique_tts_count_distribution": tts_count_dist,
        "tss_span_distribution_bp": tss_span_dist,
        "tts_span_distribution_bp": tts_span_dist,
        "dedup_tss_count_distribution": dedup_tss_count_dist,
        "dedup_tts_count_distribution": dedup_tts_count_dist,
        "redundant_tss_count_distribution": redundant_tss_count_dist,
        "redundant_tts_count_distribution": redundant_tts_count_dist,
    }
    return metrics, distributions


def classify_genes_by_variation(
    isoforms: List[dict],
    min_exons: int = 2,
) -> Dict[str, int]:
    """Classify genes by their type of isoform variation.

    For each gene (extracted via ``gene_from_name``), multi-exon isoforms are
    grouped by junction chain.  The gene is then placed in exactly one
    category:

      - **single_isoform**: only one multi-exon isoform
      - **alt_ends_only**: all isoforms share one SJC but differ in TSS/TTS
      - **alt_splicing_only**: multiple SJCs, each with a single TSS+TTS combo
      - **alt_splicing_and_ends**: multiple SJCs *and* at least one SJC has
        >1 unique TSS or TTS

    Returns a dict of category → count.
    """
    genes: Dict[str, List[dict]] = defaultdict(list)
    for iso in isoforms:
        if int(iso.get("n_exons", 1)) < min_exons:
            continue
        g = _gene_from_name(iso["name"])
        genes[g].append(iso)

    counts = {
        "single_isoform": 0,
        "alt_ends_only": 0,
        "alt_splicing_only": 0,
        "alt_splicing_and_ends": 0,
    }

    for gene_id, gene_isos in genes.items():
        if len(gene_isos) == 1:
            counts["single_isoform"] += 1
            continue

        # Group isoforms by junction chain
        sjc_groups: Dict[tuple, List[dict]] = defaultdict(list)
        for iso in gene_isos:
            key = (iso["chrom"], iso["strand"], iso["junctions"])
            sjc_groups[key].append(iso)

        n_sjcs = len(sjc_groups)

        # Check if any SJC has alternative ends
        has_alt_ends = False
        for sjc_isos in sjc_groups.values():
            n_tss = len({iso["tss"] for iso in sjc_isos})
            n_tts = len({iso["tts"] for iso in sjc_isos})
            if n_tss > 1 or n_tts > 1:
                has_alt_ends = True
                break

        if n_sjcs == 1:
            # Single SJC across all isoforms — variation is ends only
            counts["alt_ends_only"] += 1
        elif has_alt_ends:
            counts["alt_splicing_and_ends"] += 1
        else:
            counts["alt_splicing_only"] += 1

    return counts

