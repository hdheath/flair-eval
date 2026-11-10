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

  --lraa-tracking: LRAA <prefix>.quant.tracking[.gz] from LRAA bulk runs
      TSV with header columns including read_name, gene_id, transcript_id,
      transcript_splice_hash_code, frac_assigned. Rows with empty
      transcript_id or non-positive frac_assigned are skipped.

  --lraa-quant-expr / --counts-output: LRAA <prefix>.quant.expr
      TSV with header columns including transcript_id and all_reads. Emits
      a `#feature_id<TAB>count` TSV (IsoQuant convention) consumed by
      ted.py --counts. Use --counts-output (not --output) for this mode.

Usage:
  python convert_read_map.py --isoquant-model-reads model_reads.tsv.gz --output read.map.txt
  python convert_read_map.py --bambu-read-map bambu_reads.tsv --output read.map.txt
  python convert_read_map.py --isoseq-read-stat read_stat.txt --output read.map.txt
  python convert_read_map.py --flames-realign-bam sample_realign2transcript.bam --output read.map.txt
  python convert_read_map.py --stringtie2-gtf assembled.gtf --output read.map.txt
  python convert_read_map.py --lraa-tracking sample.quant.tracking.gz --output read.map.txt
  python convert_read_map.py --lraa-quant-expr sample.quant.expr --counts-output transcript_counts.tsv
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


def convert_lraa_tracking(input_path: Path, output_path: Path) -> int:
    """Convert LRAA <prefix>.quant.tracking[.gz] to FLAIR read map format.

    LRAA's tracking TSV (see util/sc/singlecell_tracking_to_sparse_matrix.py
    upstream) has header columns including:
      read_name, gene_id, transcript_id, transcript_splice_hash_code, frac_assigned

    Reads with empty transcript_id, transcript_id == '*' or '.' / 'NA', or
    frac_assigned <= 0 are skipped. A read may appear multiple times if LRAA
    assigned it fractionally to several transcripts; all positive
    assignments are kept (same convention as Bambu).

    Returns number of isoforms written.
    """
    iso_to_reads = defaultdict(list)
    opener = gzip.open if str(input_path).endswith('.gz') else open

    with opener(input_path, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            transcript_id = (row.get('transcript_id') or '').strip()
            read_id = (row.get('read_name') or row.get('read_id') or '').strip()
            if not transcript_id or transcript_id in ('NA', '.', '*'):
                continue
            if not read_id:
                continue
            try:
                frac = float(row.get('frac_assigned', '1') or '1')
            except ValueError:
                frac = 1.0
            if frac <= 0:
                continue
            iso_to_reads[transcript_id].append(read_id)

    with open(output_path, 'w') as out:
        for iso_id in sorted(iso_to_reads.keys()):
            reads = iso_to_reads[iso_id]
            out.write(f"{iso_id}\t{','.join(reads)}\n")

    return len(iso_to_reads)


def convert_lraa_quant_expr(input_path: Path, output_path: Path) -> int:
    """Convert LRAA <prefix>.quant.expr to a ted.py --counts compatible TSV.

    LRAA's quant.expr is a tab-separated table with header columns including
    `transcript_id` and `all_reads` (the integer read count per transcript;
    also has `TPM` and `RPM_total_reads`). We emit a 2-column TSV with
    `#feature_id<TAB>count` (matching IsoQuant's convention parsed by
    _load_supported_ids in bin/evaluation/ted_core.py).

    Returns number of transcripts written.
    """
    opener = gzip.open if str(input_path).endswith('.gz') else open

    rows = []
    with opener(input_path, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if reader.fieldnames is None:
            raise RuntimeError(f"No header found in LRAA quant.expr: {input_path}")
        if 'all_reads' not in reader.fieldnames:
            raise RuntimeError(
                f"Expected 'all_reads' column in LRAA quant.expr, got: {reader.fieldnames}"
            )
        # transcript_id is the canonical key; some shards use 'transcript' as a fallback.
        tid_key = 'transcript_id' if 'transcript_id' in reader.fieldnames else (
            'transcript' if 'transcript' in reader.fieldnames else None
        )
        if tid_key is None:
            raise RuntimeError(
                f"Expected 'transcript_id' column in LRAA quant.expr, got: {reader.fieldnames}"
            )
        for row in reader:
            tid = (row.get(tid_key) or '').strip()
            if not tid:
                continue
            try:
                count = float(row.get('all_reads', '0') or '0')
            except ValueError:
                continue
            rows.append((tid, count))

    with open(output_path, 'w') as out:
        out.write("#feature_id\tcount\n")
        for tid, count in rows:
            # ted.py's _load_supported_ids parses any numeric value, so we
            # preserve the float-string form (LRAA may emit fractional EM
            # estimates). Use a compact format with up to 3 decimals.
            if count == int(count):
                out.write(f"{tid}\t{int(count)}\n")
            else:
                out.write(f"{tid}\t{count:.3f}\n")

    return len(rows)


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
    input_group.add_argument(
        "--lraa-tracking", type=Path,
        help="LRAA <prefix>.quant.tracking[.gz] (per-read transcript assignments)"
    )
    input_group.add_argument(
        "--lraa-quant-expr", type=Path,
        help="LRAA <prefix>.quant.expr (per-transcript counts; pair with --counts-output)"
    )
    parser.add_argument("--output", "-o", type=Path, required=False,
                        help="Output file in FLAIR isoform.read.map.txt format (read-map modes)")
    parser.add_argument("--counts-output", type=Path, required=False,
                        help="Output TSV in #feature_id<TAB>count form (used with --lraa-quant-expr)")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    # --output is required for read-map modes; --counts-output is required
    # for --lraa-quant-expr.  Enforce here so the mutually-exclusive group
    # can still keep a single CLI shape.
    if args.lraa_quant_expr:
        if not args.counts_output:
            parser.error("--counts-output is required with --lraa-quant-expr")
    else:
        if not args.output:
            parser.error("--output is required for read-map conversion modes")

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
    elif args.lraa_tracking:
        if not args.lraa_tracking.exists():
            print(f"ERROR: Input file not found: {args.lraa_tracking}", file=sys.stderr)
            sys.exit(1)
        n = convert_lraa_tracking(args.lraa_tracking, args.output)
        if args.verbose:
            print(f"Converted LRAA tracking: {n} isoforms -> {args.output}")
    elif args.lraa_quant_expr:
        if not args.lraa_quant_expr.exists():
            print(f"ERROR: Input file not found: {args.lraa_quant_expr}", file=sys.stderr)
            sys.exit(1)
        n = convert_lraa_quant_expr(args.lraa_quant_expr, args.counts_output)
        if args.verbose:
            print(f"Converted LRAA quant.expr: {n} transcripts -> {args.counts_output}")


if __name__ == "__main__":
    main()
