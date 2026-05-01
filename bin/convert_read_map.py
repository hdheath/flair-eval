#!/usr/bin/env python3
"""
Convert read-to-isoform assignment files from various assemblers into FLAIR's
isoform.read.map.txt format:

    isoform_id\tread1,read2,read3,...

Supported input formats:
  --isoquant-model-reads: IsoQuant transcript_model_reads.tsv(.gz)
      Two columns: #read_id  transcript_id
      Unassigned reads have '*' for transcript_id.

  --bambu-read-map: Bambu read-to-transcript TSV (from trackReads export)
      Two columns: read_id  transcript_id
      One row per assignment (a read may appear multiple times for ambiguous).

  --isoseq-read-stat: IsoSeq read_stat.txt from `isoseq collapse`
      Columns include: id, length, is_fl, stat, pbid
      Reads with a valid pbid are grouped by cluster.

  --flames-realign-bam: FLAMES realign2transcript BAM from BulkPipeline
      BAM where reads are aligned to transcript sequences; reference name is
      the transcript_id and query name is the read_id.

  --stringtie2-gtf: StringTie2 output GTF
      Creates an empty map because StringTie2 does not provide per-read
      transcript assignments.

Usage:
  python convert_read_map.py --isoquant-model-reads model_reads.tsv.gz --output read.map.txt
  python convert_read_map.py --bambu-read-map bambu_reads.tsv --output read.map.txt
  python convert_read_map.py --isoseq-read-stat read_stat.txt --output read.map.txt
  python convert_read_map.py --flames-realign-bam sample_realign2transcript.bam --output read.map.txt
  python convert_read_map.py --stringtie2-gtf assembled.gtf --output read.map.txt
"""

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path


def convert_isoquant_model_reads(input_path: Path, output_path: Path) -> int:
    """Convert IsoQuant transcript_model_reads.tsv(.gz) to FLAIR read map format.

    IsoQuant format: read_id<TAB>transcript_id (one row per read assignment).
    Unassigned reads have '*' as transcript_id — these are skipped.

    Returns number of isoforms written.
    """
    iso_to_reads = defaultdict(list)
    opener = gzip.open if str(input_path).endswith('.gz') else open

    with opener(input_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            read_id, transcript_id = parts[0], parts[1]
            if transcript_id == '*' or not transcript_id:
                continue
            iso_to_reads[transcript_id].append(read_id)

    with open(output_path, 'w') as out:
        for iso_id in sorted(iso_to_reads.keys()):
            reads = iso_to_reads[iso_id]
            out.write(f"{iso_id}\t{','.join(reads)}\n")

    return len(iso_to_reads)


def convert_bambu_read_map(input_path: Path, output_path: Path) -> int:
    """Convert Bambu read-to-transcript assignment TSV to FLAIR read map format.

    Bambu format: read_id<TAB>transcript_id (one row per assignment).
    A read can appear multiple times if ambiguously assigned — we keep ALL
    assignments (same behavior as FLAIR's map which lists all reads per isoform).

    Returns number of isoforms written.
    """
    iso_to_reads = defaultdict(list)
    opener = gzip.open if str(input_path).endswith('.gz') else open

    with opener(input_path, 'rt') as f:
        for line in f:
            if line.startswith('#') or line.startswith('read_id') or line.startswith('readId'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            read_id, transcript_id = parts[0], parts[1]
            if not transcript_id or transcript_id in ('NA', '.', '*'):
                continue
            iso_to_reads[transcript_id].append(read_id)

    with open(output_path, 'w') as out:
        for iso_id in sorted(iso_to_reads.keys()):
            reads = iso_to_reads[iso_id]
            out.write(f"{iso_id}\t{','.join(reads)}\n")

    return len(iso_to_reads)


def convert_isoseq_read_stat(input_path: Path, output_path: Path) -> int:
    """Convert IsoSeq read_stat.txt from `isoseq collapse` to FLAIR read map format.

    The read_stat.txt file has columns: id, length, is_fl, stat, pbid
    where 'pbid' is the collapsed transcript cluster ID (e.g., PB.1.1).
    Reads with stat != 'unique' or missing pbid are skipped.

    Returns number of isoforms written.
    """
    iso_to_reads = defaultdict(list)
    opener = gzip.open if str(input_path).endswith('.gz') else open

    with opener(input_path, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            read_id = row.get('id', '').strip()
            pbid = row.get('pbid', '').strip()
            if not read_id or not pbid or pbid in ('NA', '*', ''):
                continue
            iso_to_reads[pbid].append(read_id)

    with open(output_path, 'w') as out:
        for iso_id in sorted(iso_to_reads.keys()):
            reads = iso_to_reads[iso_id]
            out.write(f"{iso_id}\t{','.join(reads)}\n")

    return len(iso_to_reads)


def convert_flames_realign_bam(input_path: Path, output_path: Path) -> int:
    """Convert FLAMES realign2transcript BAM to FLAIR read map format.

    FLAMES' BulkPipeline produces a realign2transcript BAM where reads are
    aligned to transcript sequences. The BAM reference name is the transcript_id
    and the query name is the read_id. Unmapped reads and secondary/supplementary
    alignments are skipped.

    Returns number of isoforms written.
    """
    import pysam

    iso_to_reads = defaultdict(list)

    with pysam.AlignmentFile(str(input_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            transcript_id = read.reference_name
            read_id = read.query_name
            if not transcript_id or not read_id:
                continue
            iso_to_reads[transcript_id].append(read_id)

    with open(output_path, 'w') as out:
        for iso_id in sorted(iso_to_reads.keys()):
            reads = iso_to_reads[iso_id]
            out.write(f"{iso_id}\t{','.join(reads)}\n")

    return len(iso_to_reads)


def convert_stringtie2_gtf(input_path: Path, output_path: Path) -> int:
    """Create an empty read map for StringTie2 output.

    StringTie2 does not natively produce per-read transcript assignments.
    A transcript_id -> transcript_id placeholder looks like a valid map but
    corrupts assignment-rate and read-end entropy metrics, so downstream
    evaluation should see an existing zero-byte file and skip read metrics.
    """
    output_path.write_text("")
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Convert assembler read-to-isoform assignments to FLAIR read map format"
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--isoquant-model-reads", type=Path,
        help="IsoQuant transcript_model_reads.tsv(.gz)"
    )
    input_group.add_argument(
        "--bambu-read-map", type=Path,
        help="Bambu read-to-transcript TSV (from trackReads export)"
    )
    input_group.add_argument(
        "--isoseq-read-stat", type=Path,
        help="IsoSeq read_stat.txt from isoseq collapse"
    )
    input_group.add_argument(
        "--flames-realign-bam", type=Path,
        help="FLAMES realign2transcript BAM from BulkPipeline"
    )
    input_group.add_argument(
        "--stringtie2-gtf", type=Path,
        help="StringTie2 output GTF (transcript-level placeholder map)"
    )
    parser.add_argument("--output", "-o", type=Path, required=True,
                        help="Output file in FLAIR isoform.read.map.txt format")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    if args.isoquant_model_reads:
        if not args.isoquant_model_reads.exists():
            print(f"ERROR: Input file not found: {args.isoquant_model_reads}", file=sys.stderr)
            sys.exit(1)
        n = convert_isoquant_model_reads(args.isoquant_model_reads, args.output)
        if args.verbose:
            print(f"Converted IsoQuant model reads: {n} isoforms -> {args.output}")
    elif args.bambu_read_map:
        if not args.bambu_read_map.exists():
            print(f"ERROR: Input file not found: {args.bambu_read_map}", file=sys.stderr)
            sys.exit(1)
        n = convert_bambu_read_map(args.bambu_read_map, args.output)
        if args.verbose:
            print(f"Converted Bambu read map: {n} isoforms -> {args.output}")
    elif args.isoseq_read_stat:
        if not args.isoseq_read_stat.exists():
            print(f"ERROR: Input file not found: {args.isoseq_read_stat}", file=sys.stderr)
            sys.exit(1)
        n = convert_isoseq_read_stat(args.isoseq_read_stat, args.output)
        if args.verbose:
            print(f"Converted IsoSeq read stat: {n} isoforms -> {args.output}")
    elif args.flames_realign_bam:
        if not args.flames_realign_bam.exists():
            print(f"ERROR: Input file not found: {args.flames_realign_bam}", file=sys.stderr)
            sys.exit(1)
        n = convert_flames_realign_bam(args.flames_realign_bam, args.output)
        if args.verbose:
            print(f"Converted FLAMES realign BAM: {n} isoforms -> {args.output}")
    elif args.stringtie2_gtf:
        if not args.stringtie2_gtf.exists():
            print(f"ERROR: Input file not found: {args.stringtie2_gtf}", file=sys.stderr)
            sys.exit(1)
        n = convert_stringtie2_gtf(args.stringtie2_gtf, args.output)
        if args.verbose:
            print(f"Created empty StringTie2 read map: {n} assignments -> {args.output}")


if __name__ == "__main__":
    main()
