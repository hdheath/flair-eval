#!/usr/bin/env python3
"""
Synthesize TED and FLAIR evaluation results into unified summary tables.

This script combines TED metrics (.tsv files) and FLAIR evaluation results (.tsv files)
into comprehensive summary tables with consistent formatting.

This is the CLI entry point. All logic is implemented in the evaluation package.
"""

import argparse
from pathlib import Path

from evaluation import (
    process_evaluation_files,
    write_tsv_merged,
)


def main():
    parser = argparse.ArgumentParser(description="Synthesize TED and FLAIR evaluation results")
    parser.add_argument('--ted-files', nargs='+', required=True,
                       help='TED evaluation TSV files')
    parser.add_argument('--flair-files', nargs='+', required=True,
                       help='FLAIR evaluation TSV files')
    parser.add_argument('--output', required=True,
                       help='Output TSV file for unified results')
    parser.add_argument('--test-name', required=True,
                       help='Test set name for labeling')

    args = parser.parse_args()

    # Dictionary to hold results keyed by sample identifier
    results_by_sample = {}

    # Process both TED and FLAIR files using unified function
    process_evaluation_files(args.ted_files, 'TED', args.test_name, results_by_sample)
    process_evaluation_files(args.flair_files, 'FLAIR', args.test_name, results_by_sample)

    if not results_by_sample:
        print("No valid results found!")
        return 1

    # Convert to list of rows
    all_results = list(results_by_sample.values())

    # Write unified TSV
    success = write_tsv_merged(all_results, args.output, args.test_name)

    if success:
        print(f"Saved unified evaluation summary to {args.output}")
        print(f"Total merged samples: {len(all_results)}")
        return 0
    else:
        print("Failed to write output file")
        return 1


if __name__ == '__main__':
    exit(main())
