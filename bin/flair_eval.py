#!/usr/bin/env python3
"""
FLAIR Evaluation Script

Performs comprehensive evaluation of FLAIR isoform predictions including:
- Genic region analysis
- Splice junction evaluation
- Transcript classification (FSM, ISM, NIC, NNC, SEM, SEN)

Usage:
    python flair_eval.py --reads-bed <reads.bed> --isoforms-bed <isoforms.bed> --gtf <annotation.gtf> --output <output.txt>

This is the CLI entry point. All logic is implemented in the evaluation package.
"""

import argparse
import sys

from evaluation import (
    get_chromtoint,
    get_regions,
    get_intersect_count,
    extract_sj_info,
    parse_gtf_transcripts,
    build_reference_structures,
    classify_transcripts,
)


def main():
    parser = argparse.ArgumentParser(description='Evaluate FLAIR isoform predictions')
    parser.add_argument('--reads-bed', required=True, help='Input reads BED file')
    parser.add_argument('--isoforms-bed', required=True, help='FLAIR isoforms BED file')
    parser.add_argument('--gtf', required=True, help='Reference annotation GTF file')
    parser.add_argument('--output', required=True, help='Output evaluation summary file (TSV format)')
    parser.add_argument('--verbose', action='store_true', help='Print verbose output')

    # Metadata arguments
    parser.add_argument('--test-name', help='Test set name')
    parser.add_argument('--dataset-name', help='Dataset name')
    parser.add_argument('--align-mode', help='Alignment mode')
    parser.add_argument('--partition-mode', help='Partition mode')
    parser.add_argument('--pipeline-mode', help='Pipeline mode (e.g., collapse_with-gtf_default)')
    parser.add_argument('--stage', help='Stage name (e.g., collapse, transcriptome)')

    args = parser.parse_args()

    if args.verbose:
        print(f"Evaluating: {args.isoforms_bed}")
        print(f"Against reads: {args.reads_bed}")
        print(f"Using annotation: {args.gtf}")

    # EVALUATE GENIC REGIONS
    chromtoint, totreads = get_chromtoint(args.reads_bed)
    reads_intervals_file = args.reads_bed.replace('.bed', '.intervals.bed')
    totregions = get_regions(chromtoint, reads_intervals_file)

    chromtoint_iso, _ = get_chromtoint(args.isoforms_bed)
    iso_intervals_file = args.isoforms_bed.replace('.bed', '.intervals.bed')
    get_regions(chromtoint_iso, iso_intervals_file)

    foundregions = get_intersect_count(reads_intervals_file, iso_intervals_file)
    genicreads = get_intersect_count(args.reads_bed, iso_intervals_file)

    # EVALUATE SPLICE JUNCTIONS
    read_sjc, read_se_ends = extract_sj_info(args.reads_bed)
    found_sjc, found_se_ends = extract_sj_info(args.isoforms_bed)

    found_subsets = {}
    for cs in found_sjc:
        found_subsets[cs] = set()
        for sjc in found_sjc[cs]:
            for slen in range(len(sjc)-1, 0, -1):
                for i in range(0, len(sjc)-slen+1):
                    found_subsets[cs].add(sjc[i:i+slen])

    tot_sjc, sup_sjc, tot_se, sup_se = 0, 0, 0, 0
    subset_sjc = 0
    for cs in read_sjc:
        if cs in found_sjc:
            for sjc in read_sjc[cs]:
                if sjc in found_sjc[cs]:
                    sup_sjc += read_sjc[cs][sjc]
                elif sjc in found_subsets[cs]:
                    subset_sjc += read_sjc[cs][sjc]
                tot_sjc += read_sjc[cs][sjc]
        else:
            for sjc in read_sjc[cs]:
                tot_sjc += read_sjc[cs][sjc]

    for cs in read_se_ends:
        if cs in found_se_ends:
            for se in read_se_ends[cs]:
                if se in found_se_ends[cs]:
                    sup_se += read_se_ends[cs][se]
                tot_se += read_se_ends[cs][se]
        else:
            for se in read_se_ends[cs]:
                tot_se += read_se_ends[cs][se]

    # EVALUATE TRANSCRIPT CLASSIFICATION
    transcripttoexons = parse_gtf_transcripts(args.gtf)
    refjuncs, refjuncchains, refseends = build_reference_structures(transcripttoexons)
    fsm, ism, nic, nnc, sem, sen, tot = classify_transcripts(args.isoforms_bed, refjuncs, refjuncchains, refseends)

    # Write output as simple TSV (one header row, one data row)
    with open(args.output, 'w') as outfile:
        # Build header with metadata columns first
        header = []
        if args.test_name:
            header.append('test_name')
        if args.dataset_name:
            header.append('dataset')
        if args.align_mode:
            header.append('align_mode')
        if args.partition_mode:
            header.append('partition_mode')

        # Extract transcriptome_mode from pipeline_mode if provided
        # pipeline_mode format is typically "transcriptome_<mode>" or "collapse_<mode>"
        transcriptome_mode = None
        if args.pipeline_mode:
            if args.pipeline_mode.startswith('transcriptome_'):
                transcriptome_mode = args.pipeline_mode.replace('transcriptome_', '', 1)
            elif args.pipeline_mode.startswith('collapse_'):
                transcriptome_mode = args.pipeline_mode.replace('collapse_', '', 1)
            else:
                transcriptome_mode = args.pipeline_mode
            header.append('transcriptome_mode')

        # Add metrics columns
        header.extend(['total_read_regions', 'found_regions', 'genic_reads',
                      'total_sjc', 'supported_sjc', 'subset_sjc', 'total_se', 'supported_se',
                      'FSM', 'ISM', 'NIC', 'NNC', 'SEM', 'SEN'])
        outfile.write('\t'.join(header) + '\n')

        # Build data row with metadata values first
        values = []
        if args.test_name:
            values.append(args.test_name)
        if args.dataset_name:
            values.append(args.dataset_name)
        if args.align_mode:
            values.append(args.align_mode)
        if args.partition_mode:
            values.append(args.partition_mode)
        if transcriptome_mode is not None:
            values.append(transcriptome_mode)

        # Add metrics values
        values.extend([totregions, foundregions, genicreads,
                      tot_sjc, sup_sjc, subset_sjc, tot_se, sup_se,
                      fsm, ism, nic, nnc, sem, sen])
        outfile.write('\t'.join(str(v) for v in values) + '\n')

    if args.verbose:
        print(f"Evaluation complete. Results written to {args.output}")
        print(f"Total transcripts classified: {tot}")

if __name__ == "__main__":
    main()
