"""
Read-end entropy calculations.

Provides functions for computing Shannon entropy of read end positions
relative to isoform model endpoints.
"""

import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

from .utils import get_logger
from .read_analysis import parse_isoform_ends, parse_read_map, parse_reads_bed_ends

logger = get_logger()


def shannon_entropy(values: List[int], bin_size: int = 10) -> float:
    """Compute Shannon entropy (bits) of a list of integer values, binned."""
    if not values:
        return 0.0
    # Bin the values
    binned = [v // bin_size for v in values]
    counts = defaultdict(int)
    for b in binned:
        counts[b] += 1
    n = len(binned)
    entropy = 0.0
    for c in counts.values():
        p = c / n
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy


def compute_read_end_entropy(
    iso_bed: Path,
    map_path: Path,
    reads_bed: Path,
    bin_size: int = 10,
) -> dict:
    """
    Compute per-isoform entropy of read-end positions relative to the isoform model.

    For each isoform, collects the 5' and 3' end positions of all assigned reads,
    computes the offset from the isoform model's TSS/TTS, and measures Shannon entropy
    of these offsets (binned). High entropy = reads are spread out; low = tightly clustered.

    Returns dict with:
        - tss_entropy_per_isoform: list of (isoform_id, entropy, n_reads)
        - tts_entropy_per_isoform: list of (isoform_id, entropy, n_reads)
        - all_tss_offsets: list of all read TSS offsets from model TSS
        - all_tts_offsets: list of all read TTS offsets from model TTS
        - summary metrics (mean/median entropy)
    """
    isoforms = parse_isoform_ends(iso_bed)
    iso_to_reads = parse_read_map(map_path)
    read_ends = parse_reads_bed_ends(reads_bed)

    tss_entropy_list = []  # (iso_id, entropy, n_reads)
    tts_entropy_list = []
    all_tss_offsets = []
    all_tts_offsets = []

    for iso_id, read_ids in iso_to_reads.items():
        if iso_id not in isoforms:
            continue
        iso = isoforms[iso_id]
        model_tss = iso['tss']
        model_tts = iso['tts']
        strand = iso['strand']

        tss_offsets = []
        tts_offsets = []

        for rid in read_ids:
            if rid not in read_ends:
                continue
            r = read_ends[rid]
            # Signed offset: positive = downstream of model end
            if strand == '+':
                tss_offsets.append(r['tss'] - model_tss)
                tts_offsets.append(r['tts'] - model_tts)
            else:
                # For minus strand, flip sign so positive = downstream (towards 5')
                tss_offsets.append(model_tss - r['tss'])
                tts_offsets.append(model_tts - r['tts'])

        if tss_offsets:
            tss_ent = shannon_entropy(tss_offsets, bin_size=bin_size)
            tss_entropy_list.append((iso_id, tss_ent, len(tss_offsets)))
            all_tss_offsets.extend(tss_offsets)

        if tts_offsets:
            tts_ent = shannon_entropy(tts_offsets, bin_size=bin_size)
            tts_entropy_list.append((iso_id, tts_ent, len(tts_offsets)))
            all_tts_offsets.extend(tts_offsets)

    # Summary metrics
    tss_entropies = [e for _, e, _ in tss_entropy_list]
    tts_entropies = [e for _, e, _ in tts_entropy_list]

    result = {
        "tss_entropy_per_isoform": tss_entropy_list,
        "tts_entropy_per_isoform": tts_entropy_list,
        "all_tss_offsets": all_tss_offsets,
        "all_tts_offsets": all_tts_offsets,
        "tss_entropy_mean": statistics.mean(tss_entropies) if tss_entropies else None,
        "tss_entropy_median": statistics.median(tss_entropies) if tss_entropies else None,
        "tts_entropy_mean": statistics.mean(tts_entropies) if tts_entropies else None,
        "tts_entropy_median": statistics.median(tts_entropies) if tts_entropies else None,
    }

    logger.info(f"Read-end entropy: {len(tss_entropy_list)} isoforms with TSS data, "
                f"{len(tts_entropy_list)} with TTS data")

    return result
