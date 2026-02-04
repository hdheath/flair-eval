#!/usr/bin/env python3
"""
TED (Transcriptome End Distribution) - Standard Library Version
Evaluates TSS/TTS accuracy and read assignment metrics for FLAIR outputs.
Uses only Python standard library + bedtools/samtools.

This is the CLI entry point. All logic is implemented in the evaluation package.
"""

import argparse
import csv
import logging
from pathlib import Path

from evaluation import (
    calculate_ted_metrics,
    get_timing_data,
    write_timing_report,
)
from evaluation.utils import get_logger

logger = get_logger()


def main():
    parser = argparse.ArgumentParser(
        description="TED: Transcriptome End Distribution metrics for FLAIR outputs"
    )
    # Input options - either BED+map or GTF-only
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--isoforms-bed", type=Path, help="Isoforms BED file (*.isoforms.bed)")
    input_group.add_argument("--gtf-input", type=Path, help="Isoforms GTF file (for Bambu/IsoQuant outputs)")
    
    parser.add_argument("--map-file", type=Path, help="Read map file (*.isoform.read.map.txt) - required unless --skip-read-metrics")
    parser.add_argument("--bam", type=Path, help="BAM file for alignment metrics")
    parser.add_argument("--corrected-bed", type=Path, help="Corrected BED (for collapse stage)")
    parser.add_argument("--reads-bed", type=Path, help="Reads BED12 file (for read-end entropy analysis)")
    parser.add_argument("--prime5-peaks", type=Path, help="Experimental 5' peaks (CAGE)")
    parser.add_argument("--prime3-peaks", type=Path, help="Experimental 3' peaks (QuantSeq/PolyA)")
    parser.add_argument("--ref-prime5-peaks", type=Path, help="Reference 5' peaks")
    parser.add_argument("--ref-prime3-peaks", type=Path, help="Reference 3' peaks")
    parser.add_argument("--window", type=int, default=50, help="Distance window for TSS/TTS matching (default: 50)")
    parser.add_argument("--stage", choices=["collapse", "transcriptome"], default="collapse", help="Pipeline stage")
    parser.add_argument("--output", type=Path, required=True, help="Output TSV file")
    parser.add_argument("--plot-output-dir", type=Path, help="Directory to save distance histogram plots (CAGE and QuantSeq)")
    parser.add_argument("--test-regions-dir", type=Path, help="Directory to save recoverable peak BED files")
    parser.add_argument("--genome", type=Path, help="Genome FASTA file (for motif analysis)")
    parser.add_argument("--timing-output", type=Path, help="Output file for performance timing report")
    parser.add_argument("--max-count-cage", type=int, help="Fixed y-axis limit for CAGE histograms across runs (optional)")
    parser.add_argument("--max-count-quantseq", type=int, help="Fixed y-axis limit for QuantSeq histograms across runs (optional)")
    parser.add_argument("--max-count-ref-tss", type=int, help="Fixed y-axis limit for Reference TSS histograms across runs (optional)")
    parser.add_argument("--max-count-ref-tts", type=int, help="Fixed y-axis limit for Reference TTS histograms across runs (optional)")
    # Simplified evaluation mode (for Bambu/IsoQuant - no read-level metrics)
    parser.add_argument("--skip-read-metrics", action="store_true", 
                        help="Skip read-level metrics (entropy, motifs, truncation). Use for Bambu/IsoQuant evaluation.")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    # Metadata arguments for result tracking
    parser.add_argument("--test-name", type=str, help="Test set name")
    parser.add_argument("--dataset-name", type=str, help="Dataset name")
    parser.add_argument("--align-mode", type=str, help="Alignment mode")
    parser.add_argument("--partition-mode", type=str, help="Partition mode")
    parser.add_argument("--pipeline-mode", type=str, help="Pipeline mode (process label)")

    args = parser.parse_args()

    # Validate arguments
    if not args.skip_read_metrics and not args.map_file:
        parser.error("--map-file is required unless --skip-read-metrics is specified")

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Handle GTF input - convert to BED12 internally
    iso_bed = args.isoforms_bed
    if args.gtf_input:
        import tempfile
        import subprocess
        # Convert GTF to BED12 using our converter script
        temp_bed = Path(tempfile.mktemp(suffix='.bed'))
        logger.info(f"Converting GTF to BED12: {args.gtf_input} -> {temp_bed}")
        from pathlib import Path as P
        gtf_to_bed12_script = P(__file__).parent / "gtf_to_bed12.py"
        subprocess.run([
            "python", str(gtf_to_bed12_script),
            "--gtf", str(args.gtf_input),
            "--output", str(temp_bed),
            "--verbose" if args.verbose else ""
        ], check=True)
        iso_bed = temp_bed

    # Build plot prefix from metadata
    plot_prefix_parts = []
    if args.dataset_name:
        plot_prefix_parts.append(args.dataset_name)
    if args.align_mode:
        plot_prefix_parts.append(args.align_mode)
    if args.partition_mode:
        plot_prefix_parts.append(args.partition_mode)
    if args.pipeline_mode:
        plot_prefix_parts.append(args.pipeline_mode)
    if args.stage:
        plot_prefix_parts.append(args.stage)
    plot_prefix = "_".join(plot_prefix_parts) if plot_prefix_parts else "ted"

    # Calculate metrics
    input_file = args.gtf_input if args.gtf_input else args.isoforms_bed
    logger.info(f"Calculating TED metrics for {input_file}")
    metrics = calculate_ted_metrics(
        iso_bed=iso_bed,
        map_file=args.map_file,
        bam_file=args.bam,
        corrected_bed=args.corrected_bed,
        reads_bed=args.reads_bed,
        prime5_peaks=args.prime5_peaks,
        prime3_peaks=args.prime3_peaks,
        ref_prime5_peaks=args.ref_prime5_peaks,
        ref_prime3_peaks=args.ref_prime3_peaks,
        window=args.window,
        stage=args.stage,
        plot_output_dir=args.plot_output_dir,
        plot_prefix=plot_prefix,
        test_regions_dir=args.test_regions_dir,
        genome_path=args.genome,
        max_count_cage=args.max_count_cage,
        max_count_quantseq=args.max_count_quantseq,
        max_count_ref_tss=args.max_count_ref_tss,
        max_count_ref_tts=args.max_count_ref_tts,
        skip_read_metrics=args.skip_read_metrics,
    )

    # Clean up temp file if we created one
    if args.gtf_input and iso_bed.exists():
        iso_bed.unlink()

    # Add metadata to metrics if provided
    if args.test_name:
        metrics = {'test_name': args.test_name, **metrics}
    if args.dataset_name:
        metrics = {'dataset': args.dataset_name, **metrics}
    if args.align_mode:
        metrics = {'align_mode': args.align_mode, **metrics}
    if args.partition_mode:
        metrics = {'partition_mode': args.partition_mode, **metrics}

    # Extract transcriptome_mode from pipeline_mode if provided
    # pipeline_mode format is typically "transcriptome_<mode>" or "collapse_<mode>"
    if args.pipeline_mode:
        if args.pipeline_mode.startswith('transcriptome_'):
            transcriptome_mode = args.pipeline_mode.replace('transcriptome_', '', 1)
        elif args.pipeline_mode.startswith('collapse_'):
            transcriptome_mode = args.pipeline_mode.replace('collapse_', '', 1)
        else:
            transcriptome_mode = args.pipeline_mode
        metrics = {'transcriptome_mode': transcriptome_mode, **metrics}

    # Write output as TSV
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=metrics.keys(), delimiter='\t')
        writer.writeheader()
        writer.writerow(metrics)

    logger.info(f"Wrote metrics to {args.output}")

    # Write timing report if requested
    if args.timing_output:
        _TIMING_DATA = get_timing_data()
        with open(args.timing_output, 'w') as tf:
            tf.write("TED Performance Timing Report\n")
            tf.write("=" * 80 + "\n\n")

            # Sort by total time (sum of all calls)
            sorted_sections = sorted(
                _TIMING_DATA.items(),
                key=lambda x: sum(x[1]),
                reverse=True
            )

            total_time = sum(sum(times) for times in _TIMING_DATA.values())
            tf.write(f"Total instrumented time: {total_time:.2f}s\n\n")

            tf.write(f"{'Section':<50} {'Calls':>8} {'Total (s)':>12} {'Avg (s)':>12} {'% of Total':>12}\n")
            tf.write("-" * 100 + "\n")

            for section, times in sorted_sections:
                n_calls = len(times)
                total = sum(times)
                avg = total / n_calls
                pct = (total / total_time * 100) if total_time > 0 else 0
                tf.write(f"{section:<50} {n_calls:>8} {total:>12.2f} {avg:>12.2f} {pct:>11.1f}%\n")

        logger.info(f"Wrote timing report to {args.timing_output}")

    # Print summary
    logger.info(f"Isoforms: {metrics['isoforms_observed']}")
    logger.info(f"Genes: {metrics['genes_observed']}")
    if metrics['assignment_rate']:
        logger.info(f"Assignment rate: {metrics['assignment_rate']:.2%}")
    else:
        logger.info("Assignment rate: N/A")


if __name__ == "__main__":
    main()
