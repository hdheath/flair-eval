#!/usr/bin/env python3
"""
generate_transcriptome_fa.py
============================
Extract spliced cDNA transcript sequences from a genome FASTA + GTF annotation.
Produces a FASTA suitable as Badread --reference input, with per-transcript
depth headers so Badread controls relative abundance.

Usage:
    python generate_transcriptome_fa.py \
        --genome genome.fa \
        --gtf annotation.gtf \
        --output transcriptome.fa \
        --depth 10              # uniform read depth per transcript
        --depth-map tx_depths.tsv  # per-transcript depth (optional)

Output FASTA headers:
    >TX1a depth=10
    ATCGATCG...
"""

import argparse
import sys
from pathlib import Path

try:
    import pysam
except ImportError:
    sys.exit("ERROR: pysam is required. Install with: pip install pysam")


def parse_gtf_transcripts(gtf_path: Path) -> dict:
    """Parse GTF into {transcript_id: {chrom, strand, exons: [(start, end), ...]}}."""
    transcripts = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom = fields[0]
            feature = fields[2]
            start = int(fields[3]) - 1  # convert to 0-based
            end = int(fields[4])
            strand = fields[6]
            attrs = fields[8]

            # Extract transcript_id
            tid = None
            for token in attrs.split(';'):
                token = token.strip()
                if token.startswith('transcript_id'):
                    # Handle both transcript_id "X" and transcript_id X formats
                    if '"' in token:
                        tid = token.split('"')[1]
                    else:
                        tid = token.split()[-1]
                    break
            if not tid:
                continue

            if feature == 'transcript':
                transcripts[tid] = {
                    'chrom': chrom,
                    'strand': strand,
                    'exons': [],
                }
            elif feature == 'exon' and tid in transcripts:
                transcripts[tid]['exons'].append((start, end))

    # Sort exons by position
    for tid in transcripts:
        transcripts[tid]['exons'].sort()

    return transcripts


def reverse_complement(seq: str) -> str:
    """Reverse complement a DNA sequence."""
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(comp)[::-1]


def extract_sequences(genome: pysam.FastaFile, transcripts: dict,
                      default_depth: int, depth_map: dict) -> list:
    """Extract spliced cDNA sequences from genome for each transcript."""
    results = []
    genome_chroms = set(genome.references)

    for tid, info in transcripts.items():
        if not info['exons']:
            continue
        if info['chrom'] not in genome_chroms:
            print(f"WARNING: {info['chrom']} not in genome FASTA, skipping {tid}",
                  file=sys.stderr)
            continue

        # Extract and concatenate exon sequences
        seq_parts = []
        for start, end in info['exons']:
            chrom_len = genome.get_reference_length(info['chrom'])
            # Clamp to chromosome bounds
            s = max(0, start)
            e = min(chrom_len, end)
            if e > s:
                seq_parts.append(genome.fetch(info['chrom'], s, e))

        if not seq_parts:
            continue

        cdna_seq = ''.join(seq_parts).upper()

        # For minus-strand transcripts, reverse complement
        if info['strand'] == '-':
            cdna_seq = reverse_complement(cdna_seq)

        depth = depth_map.get(tid, default_depth)
        results.append((tid, depth, cdna_seq, info['strand'], len(info['exons'])))

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Extract spliced cDNA sequences for Badread simulation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--genome', required=True, type=Path,
                        help='Reference genome FASTA (must be indexed)')
    parser.add_argument('--gtf', required=True, type=Path,
                        help='GTF annotation file')
    parser.add_argument('--output', required=True, type=Path,
                        help='Output transcriptome FASTA')
    parser.add_argument('--depth', type=int, default=10,
                        help='Uniform read depth per transcript (default: 10)')
    parser.add_argument('--depth-map', type=Path, default=None,
                        help='TSV file with transcript_id<TAB>depth overrides')
    parser.add_argument('--min-length', type=int, default=50,
                        help='Minimum transcript length to include (default: 50)')
    parser.add_argument('--region', type=str, nargs='+', default=None,
                        help='Restrict to genomic region(s) (e.g. chr1 chr3:48000000-53000000)')
    args = parser.parse_args()

    # Load optional per-transcript depth map
    depth_map = {}
    if args.depth_map and args.depth_map.exists():
        with open(args.depth_map) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    depth_map[parts[0]] = int(parts[1])

    # Parse GTF
    transcripts = parse_gtf_transcripts(args.gtf)
    print(f"Parsed {len(transcripts)} transcripts from {args.gtf.name}",
          file=sys.stderr)

    # Filter by region(s) if specified
    if args.region:
        filtered = {}
        for region in args.region:
            if ':' in region:
                chrom, coords = region.split(':', 1)
                start, end = [int(x.replace(',', '')) for x in coords.split('-')]
                for tid, info in transcripts.items():
                    if info['chrom'] == chrom and any(e[0] >= start and e[1] <= end for e in info['exons']):
                        filtered[tid] = info
            else:
                # Just a chromosome name
                for tid, info in transcripts.items():
                    if info['chrom'] == region:
                        filtered[tid] = info
        transcripts = filtered
        print(f"After region filter '{args.region}': {len(transcripts)} transcripts",
              file=sys.stderr)

    # Extract sequences
    genome = pysam.FastaFile(str(args.genome))
    results = extract_sequences(genome, transcripts, args.depth, depth_map)
    genome.close()

    # Write output FASTA
    n_written = 0
    with open(args.output, 'w') as fh:
        for tid, depth, seq, strand, n_exons in results:
            if len(seq) < args.min_length:
                continue
            fh.write(f">{tid} depth={depth} strand={strand} exons={n_exons} length={len(seq)}\n")
            # Write sequence in 80-char lines
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + '\n')
            n_written += 1

    print(f"Wrote {n_written} transcript sequences to {args.output.name}",
          file=sys.stderr)


if __name__ == '__main__':
    main()
