#!/usr/bin/env python3
"""
FLAIR Evaluation Script

Computes SQANTI-style transcript classification (FSM, ISM, NIC, NNC, SEM, SEN)
for an isoform set against a reference GTF.

Previously also computed read-based metrics (genic-region overlap, splice-chain
support stats), but those columns were never visualized in any downstream plot
and cost ~1.5h on full-genome ONT data due to repeated multi-pass scans of the
3GB reads BED.  They have been removed; the columns remain in the output TSV
header for back-compat but are written as empty strings.

Usage:
    python flair_eval.py --isoforms-bed <isoforms.bed> --gtf <annotation.gtf> --output <output.txt>

This is the CLI entry point. All logic is implemented in the evaluation package.
"""

import argparse
import sys
import tempfile
import subprocess
from pathlib import Path

from evaluation import (
    parse_gtf_transcripts,
    build_reference_structures,
    classify_transcripts_per_isoform,
)


def load_supported_ids(counts_path: str, min_support: int) -> set[str] | None:
    """Parse supported transcript IDs from IsoQuant/Bambu/FLAIR count tables."""
    if not counts_path:
        return None
    path = Path(counts_path)
    if not path.exists():
        return None
    supported: set[str] = set()
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


def is_supported_bed_name(name: str, supported_ids: set[str]) -> bool:
    """Match direct IDs and FLAIR-style transcript_gene BED names."""
    if name in supported_ids:
        return True
    if "_" in name and name.split("_", 1)[0] in supported_ids:
        return True
    return False


def filter_bed_by_supported_ids(
    bed_path: str,
    supported_ids: set[str],
    *,
    min_support: int,
    counts_path: str,
    verbose: bool,
) -> str:
    """Write a temporary BED containing only count-supported transcripts."""
    filtered_bed = tempfile.mktemp(suffix=f".support_ge{min_support}.bed")
    n_pre = 0
    n_post = 0
    with open(bed_path) as fin, open(filtered_bed, "w") as fout:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            n_pre += 1
            if is_supported_bed_name(parts[3], supported_ids):
                fout.write(line)
                n_post += 1
    if verbose:
        print(
            f"Read-support filter: kept {n_post}/{n_pre} isoforms "
            f"(count >= {min_support} in {Path(counts_path).name})",
            flush=True,
        )
    return filtered_bed


def main():
    parser = argparse.ArgumentParser(description='Evaluate FLAIR isoform predictions')
    # --reads-bed kept for CLI back-compat; no longer used.
    parser.add_argument('--reads-bed', help='[DEPRECATED] Input reads BED file '
                                            '(no longer used; kept for CLI compatibility)')
    # Input options - either BED or GTF
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--isoforms-bed', help='FLAIR isoforms BED file')
    input_group.add_argument('--gtf-input', help='Isoforms GTF/GFF file (for GTF/GFF-based assemblers)')

    parser.add_argument('--gtf', required=True, help='Reference annotation GTF file')
    parser.add_argument('--output', required=True, help='Output evaluation summary file (TSV format)')
    parser.add_argument('--verbose', action='store_true', help='Print verbose output')

    # Metadata arguments
    parser.add_argument('--test-name', help='Test set name')
    parser.add_argument('--dataset-name', help='Dataset name')
    parser.add_argument('--library-type', help='Library type (e.g., pacbio_cDNA, ont_cDNA, ont_dRNA)')
    parser.add_argument('--align-mode', help='Alignment mode')
    parser.add_argument('--partition-mode', help='Partition mode')
    parser.add_argument('--pipeline-mode', help='Pipeline mode (e.g., collapse_with-gtf_default)')
    parser.add_argument('--stage', help='Stage name (e.g., collapse, transcriptome)')
    parser.add_argument('--plot-output-dir', help='Directory to save structural evaluation plots')
    parser.add_argument('--plot-prefix', default='', help='Prefix for plot filenames')
    parser.add_argument('--categories-output', help='Path to write per-isoform category TSV (name, category)')
    parser.add_argument('--counts', help='Optional per-transcript counts table for support filtering')
    parser.add_argument('--min-support', type=int, default=1,
                        help='Minimum transcript count retained when --counts is provided')

    args = parser.parse_args()

    # Handle GTF input - convert to BED12
    isoforms_bed = args.isoforms_bed
    temp_bed = None
    if args.gtf_input:
        # Convert GTF to BED12 using our converter script
        temp_bed = tempfile.mktemp(suffix='.bed')
        if args.verbose:
            print(f"Converting GTF to BED12: {args.gtf_input} -> {temp_bed}", flush=True)
        gtf_to_bed12_script = Path(__file__).parent / "gtf_to_bed12.py"
        cmd = ["python", str(gtf_to_bed12_script), "--gtf", args.gtf_input, "--output", temp_bed]
        if args.verbose:
            cmd.append("--verbose")
        subprocess.run(cmd, check=True)
        isoforms_bed = temp_bed

    supported_ids = load_supported_ids(args.counts, args.min_support)
    filtered_bed = None
    if supported_ids is not None:
        filtered_bed = filter_bed_by_supported_ids(
            isoforms_bed,
            supported_ids,
            min_support=args.min_support,
            counts_path=args.counts,
            verbose=args.verbose,
        )
        isoforms_bed = filtered_bed

    if args.verbose:
        print(f"Evaluating: {isoforms_bed}", flush=True)
        print(f"Using annotation: {args.gtf}", flush=True)

    # Read-based metrics (genic regions, splice-chain support) are removed.
    # The columns stay in the output for back-compat but are empty strings.
    totregions = ''
    foundregions = ''
    genicreads = ''
    tot_sjc = ''
    sup_sjc = ''
    subset_sjc = ''
    tot_se = ''
    sup_se = ''

    # EVALUATE TRANSCRIPT CLASSIFICATION
    if args.verbose:
        print(f"Parsing reference GTF: {args.gtf}", flush=True)
    transcripttoexons = parse_gtf_transcripts(args.gtf)
    if args.verbose:
        print(f"Building reference structures", flush=True)
    refjuncs, refjuncchains, refseends = build_reference_structures(transcripttoexons)

    # Use per-isoform classification so we can both count totals AND emit labels
    iso_classifications = classify_transcripts_per_isoform(
        isoforms_bed, refjuncs, refjuncchains, refseends
    )
    fsm = sum(1 for r in iso_classifications if r["category"] == "FSM")
    ism = sum(1 for r in iso_classifications if r["category"] == "ISM")
    nic = sum(1 for r in iso_classifications if r["category"] == "NIC")
    nnc = sum(1 for r in iso_classifications if r["category"] == "NNC")
    sem = sum(1 for r in iso_classifications if r["category"] == "SEM")
    sen = sum(1 for r in iso_classifications if r["category"] == "SEN")
    tot = len(iso_classifications)

    # Write per-isoform category TSV if requested (consumed downstream by sqanti_precision.py)
    if args.categories_output:
        with open(args.categories_output, 'w') as cat_out:
            cat_out.write("isoform_name\tcategory\n")
            for r in iso_classifications:
                cat_out.write(f"{r['name']}\t{r['category']}\n")

    # Write output as simple TSV (one header row, one data row)
    with open(args.output, 'w') as outfile:
        # Build header with metadata columns first
        header = []
        if args.test_name:
            header.append('test_name')
        if args.dataset_name:
            header.append('dataset')
        if args.library_type:
            header.append('library_type')
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
        if args.library_type:
            values.append(args.library_type)
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

    # Note: transcript_classification and splice_junction_support plots removed.
    # The classification counts are still emitted to the evaluation TSV above.

    # Clean up temp file if we created one
    if filtered_bed and Path(filtered_bed).exists():
        Path(filtered_bed).unlink()
    if temp_bed and Path(temp_bed).exists():
        Path(temp_bed).unlink()

    if args.verbose:
        print(f"Evaluation complete. Results written to {args.output}", flush=True)
        print(f"Total transcripts classified: {tot}", flush=True)

if __name__ == "__main__":
    main()
