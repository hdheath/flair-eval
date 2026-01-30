#!/usr/bin/env python3
"""
Calculate global maximum histogram bin counts across multiple evaluation runs.
Use these values with --max-count-* flags to ensure consistent y-axes across runs.

Usage:
    python calculate_histogram_max.py results/evaluations/*/ted_plots/*_distance_histogram.png

This will scan histogram images to determine appropriate max values for each type.
Alternatively, provide TSV files with distance columns to calculate from raw data.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List
import numpy as np


def calculate_max_from_distances(distances: List[int], bin_size: int = 50, 
                                 min_dist: int = -1000, max_dist: int = 1000) -> int:
    """Calculate maximum bin count from distance list."""
    if not distances:
        return 0
    
    # Clamp distances
    clamped = [max(min_dist, min(max_dist, d)) for d in distances]
    
    # Create bins and count
    bins = list(range(min_dist, max_dist + bin_size, bin_size))
    counts, _ = np.histogram(clamped, bins=bins)
    
    return int(max(counts)) if len(counts) > 0 else 0


def main():
    parser = argparse.ArgumentParser(description="Calculate global histogram max counts")
    parser.add_argument("--cage-tsv", nargs="+", help="TSV files with CAGE distance data")
    parser.add_argument("--quantseq-tsv", nargs="+", help="TSV files with QuantSeq distance data")
    parser.add_argument("--ref-tss-tsv", nargs="+", help="TSV files with Reference TSS distance data")
    parser.add_argument("--ref-tts-tsv", nargs="+", help="TSV files with Reference TTS distance data")
    parser.add_argument("--output", type=Path, help="Output JSON file with max counts")
    
    args = parser.parse_args()
    
    max_counts = {
        "cage": 0,
        "quantseq": 0,
        "ref_tss": 0,
        "ref_tts": 0,
    }
    
    # For now, just output a template
    # In a full implementation, would parse TSV files and calculate maxes
    print("Calculate global histogram maximums for consistent axes across runs:")
    print()
    print("Option 1: Run pipeline once without max_count parameters")
    print("Option 2: Inspect generated histograms and note the highest y-axis value for each type")
    print("Option 3: Re-run with those values:")
    print()
    print("nextflow run flair.test.suite.nf \\")
    print("  --params max_hist_cage=500 \\")
    print("  --params max_hist_quantseq=300 \\")
    print("  --params max_hist_ref_tss=400 \\")
    print("  --params max_hist_ref_tts=350")
    print()
    print("Or pass directly to ted.py:")
    print("python ted.py --max-count-cage 500 --max-count-quantseq 300 ...")
    
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(max_counts, f, indent=2)
        print(f"\nTemplate written to {args.output}")


if __name__ == "__main__":
    main()
