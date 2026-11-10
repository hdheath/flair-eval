"""
Truncation pattern analysis functions.

Provides utilities for characterizing read start distribution patterns
relative to experimental peaks, which can indicate RT truncation.
"""

import statistics
from typing import Dict, List, Tuple


def characterize_truncation_pattern(
    read_starts: List[int],
    peak_pos: int,
    strand: str,
    sharp_threshold: int = 20,
    trailing_ratio_threshold: float = 0.3,
) -> Tuple[str, dict]:
    """
    Characterize the shape of read start distribution relative to peak.

    For 5' ends, examines whether reads cluster at the peak or trail off,
    which can indicate RT truncation patterns.

    Patterns:
    - "sharp": reads cluster tightly at peak (good signal, low truncation)
    - "trailing": reads trail off in 5' direction (truncation gradient)
    - "bimodal": two distinct clusters (possible alternative TSS)
    - "sparse": too few reads to characterize

    Args:
        read_starts: List of read start positions (5' ends)
        peak_pos: Position of the experimental peak
        strand: Strand of the peak ('+' or '-')
        sharp_threshold: Max distance from peak for "sharp" classification
        trailing_ratio_threshold: Ratio of upstream to downstream reads for "trailing"

    Returns:
        Tuple of (pattern_name, details_dict)
    """
    if len(read_starts) < 5:
        return "sparse", {"n_reads": len(read_starts)}

    # Calculate signed distances (positive = downstream of peak in transcription direction)
    if strand == '+':
        distances = [pos - peak_pos for pos in read_starts]
    else:
        distances = [peak_pos - pos for pos in read_starts]

    # Upstream (5' of peak, negative values) vs downstream (3' of peak, positive values)
    upstream = [d for d in distances if d < 0]
    downstream = [d for d in distances if d >= 0]
    at_peak = [d for d in distances if abs(d) <= sharp_threshold]

    details = {
        "n_reads": len(read_starts),
        "n_upstream": len(upstream),
        "n_downstream": len(downstream),
        "n_at_peak": len(at_peak),
        "mean_distance": statistics.mean(distances) if distances else 0,
        "median_distance": statistics.median(distances) if distances else 0,
        "std_distance": statistics.stdev(distances) if len(distances) > 1 else 0,
    }

    # Check for sharp clustering at peak
    at_peak_ratio = len(at_peak) / len(read_starts)
    if at_peak_ratio >= 0.7:
        return "sharp", details

    # Check for trailing pattern (many reads upstream/5' of peak = truncation)
    if len(upstream) > 0 and len(downstream) > 0:
        upstream_ratio = len(upstream) / len(read_starts)
        if upstream_ratio >= trailing_ratio_threshold:
            # Further check: is there a gradient (more reads closer to peak)?
            if upstream:
                upstream_dists = sorted([abs(d) for d in upstream])
                # Check if closer reads are more numerous
                close_upstream = sum(1 for d in upstream if abs(d) <= 100)
                far_upstream = sum(1 for d in upstream if abs(d) > 100)
                if close_upstream > far_upstream:
                    details["trailing_gradient"] = True
            return "trailing", details

    # Check for bimodal distribution
    if len(distances) >= 10:
        sorted_dists = sorted(distances)
        # Simple bimodal check: look for gap in middle
        mid_start = len(sorted_dists) // 3
        mid_end = 2 * len(sorted_dists) // 3
        mid_range = sorted_dists[mid_end] - sorted_dists[mid_start]
        full_range = sorted_dists[-1] - sorted_dists[0]
        if full_range > 0:
            # If middle third spans less than 20% of range, might be bimodal
            if mid_range / full_range < 0.2 and full_range > 100:
                return "bimodal", details

    # Default: dispersed pattern
    return "dispersed", details
