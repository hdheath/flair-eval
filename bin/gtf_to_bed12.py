#!/usr/bin/env python3
"""
Convert GTF file to BED12 format.

This script converts transcript annotations from GTF format to BED12 format,
which is required for evaluation scripts that expect isoforms in BED12 format.

BED12 format columns:
1. chrom - Chromosome name
2. chromStart - Start position (0-based)
3. chromEnd - End position
4. name - Transcript ID
5. score - Score (set to 0)
6. strand - Strand (+/-)
7. thickStart - Same as chromStart for transcripts
8. thickEnd - Same as chromEnd for transcripts
9. itemRgb - RGB color (set to 0)
10. blockCount - Number of exons
11. blockSizes - Comma-separated exon sizes
12. blockStarts - Comma-separated exon starts (relative to chromStart)

Usage:
    python gtf_to_bed12.py --gtf input.gtf --output output.bed
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path


def parse_gtf_attributes(attr_string: str) -> dict:
    """Parse GTF attribute string into a dictionary."""
    attrs = {}
    # Handle both formats: key "value"; and key=value;
    parts = attr_string.strip().rstrip(';').split(';')
    for part in parts:
        part = part.strip()
        if not part:
            continue
        # Try key "value" format first (standard GTF)
        if ' "' in part:
            key, value = part.split(' "', 1)
            attrs[key.strip()] = value.rstrip('"')
        # Try key=value format (GFF3-like)
        elif '=' in part:
            key, value = part.split('=', 1)
            attrs[key.strip()] = value.strip('"')
        # Try key value format (space-separated)
        elif ' ' in part:
            key, value = part.split(' ', 1)
            attrs[key.strip()] = value.strip().strip('"')
    return attrs


def gtf_to_bed12(gtf_path: Path, output_path: Path, source_filter: str = None, verbose: bool = False):
    """
    Convert GTF to BED12 format.
    
    Args:
        gtf_path: Path to input GTF file
        output_path: Path to output BED12 file
        source_filter: Optional source field filter (e.g., 'Bambu', 'isoquant')
        verbose: Print progress information
    """
    # Collect exons by transcript
    # Structure: {transcript_id: {'chrom': str, 'strand': str, 'gene_id': str, 'exons': [(start, end), ...]}}
    transcripts = defaultdict(lambda: {'chrom': None, 'strand': None, 'gene_id': None, 'exons': []})
    
    if verbose:
        print(f"Reading GTF: {gtf_path}", file=sys.stderr)
    
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            
            # Apply source filter if specified
            if source_filter and source != source_filter:
                continue
            
            # Only process exon features
            if feature != 'exon':
                continue
            
            attrs = parse_gtf_attributes(attributes)
            transcript_id = attrs.get('transcript_id')
            gene_id = attrs.get('gene_id', transcript_id)
            
            if not transcript_id:
                continue
            
            # Convert to 0-based coordinates
            start = int(start) - 1
            end = int(end)
            
            tx = transcripts[transcript_id]
            tx['chrom'] = chrom
            tx['strand'] = strand
            tx['gene_id'] = gene_id
            tx['exons'].append((start, end))
    
    if verbose:
        print(f"Found {len(transcripts)} transcripts", file=sys.stderr)
    
    # Write BED12 output
    n_written = 0
    with open(output_path, 'w') as out:
        for tx_id, tx_data in transcripts.items():
            if not tx_data['exons']:
                continue
            
            chrom = tx_data['chrom']
            strand = tx_data['strand']
            gene_id = tx_data['gene_id']
            
            # Sort exons by start position
            exons = sorted(tx_data['exons'], key=lambda x: x[0])
            
            # Calculate transcript boundaries
            tx_start = exons[0][0]
            tx_end = exons[-1][1]
            
            # Calculate block sizes and starts
            block_count = len(exons)
            block_sizes = ','.join(str(end - start) for start, end in exons)
            block_starts = ','.join(str(start - tx_start) for start, end in exons)
            
            # Build transcript name (include gene_id if different from tx_id)
            name = tx_id
            if gene_id and gene_id != tx_id:
                name = f"{gene_id}_{tx_id}"
            
            # Write BED12 line
            bed_fields = [
                chrom,           # chrom
                str(tx_start),   # chromStart
                str(tx_end),     # chromEnd
                name,            # name
                '0',             # score
                strand,          # strand
                str(tx_start),   # thickStart
                str(tx_end),     # thickEnd
                '0',             # itemRgb
                str(block_count),# blockCount
                block_sizes,     # blockSizes
                block_starts,    # blockStarts
            ]
            out.write('\t'.join(bed_fields) + '\n')
            n_written += 1
    
    if verbose:
        print(f"Wrote {n_written} transcripts to {output_path}", file=sys.stderr)
    
    return n_written


def main():
    parser = argparse.ArgumentParser(
        description='Convert GTF file to BED12 format for isoform evaluation'
    )
    parser.add_argument('--gtf', '-g', required=True, type=Path,
                        help='Input GTF file')
    parser.add_argument('--output', '-o', required=True, type=Path,
                        help='Output BED12 file')
    parser.add_argument('--source', '-s', type=str, default=None,
                        help='Filter by source field (e.g., "Bambu", "isoquant")')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print progress information')
    
    args = parser.parse_args()
    
    if not args.gtf.exists():
        print(f"Error: GTF file not found: {args.gtf}", file=sys.stderr)
        sys.exit(1)
    
    n_transcripts = gtf_to_bed12(args.gtf, args.output, args.source, args.verbose)
    
    if n_transcripts == 0:
        print("Warning: No transcripts written to output", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
