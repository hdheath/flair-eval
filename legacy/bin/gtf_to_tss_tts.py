#!/usr/bin/env python3
"""
Extract single-position TSS and TTS coordinates from GTF file.
Creates BED files suitable for TED reference peak comparison.
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict
import gzip

def open_file(path):
    """Open regular or gzipped file."""
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def parse_gtf_attributes(attr_string):
    """Parse GTF attribute string into dict."""
    attrs = {}
    for item in attr_string.strip().split(';'):
        item = item.strip()
        if not item:
            continue
        if ' ' not in item:
            continue
        key, value = item.split(' ', 1)
        # Remove quotes
        value = value.strip('"').strip("'")
        attrs[key] = value
    return attrs


def extract_transcripts_from_gtf(gtf_path, feature_type='transcript'):
    """
    Extract transcript coordinates from GTF.
    
    Returns dict: transcript_id -> (chrom, start, end, strand, gene_id)
    """
    transcripts = {}
    
    # with open_file(gtf_path) as fh:
    #     for line in fh:
    #         if line.startswith('#'):
    #             continue
            
    #         parts = line.strip().split('\t')
    #         if len(parts) < 9:
    #             continue
            
    #         chrom, source, feature, start, end, score, strand, frame, attributes = parts
            
    #         # Only process transcript features
    #         if feature != feature_type:
    #             continue
            
    #         attrs = parse_gtf_attributes(attributes)
    #         transcript_id = attrs.get('transcript_id')
    #         gene_id = attrs.get('gene_id', 'unknown')
            
    #         if not transcript_id:
    #             continue
            
    #         # GTF is 1-based, convert to 0-based for BED
    #         start_pos = int(start) - 1
    #         end_pos = int(end)
            
    #         transcripts[transcript_id] = (chrom, start_pos, end_pos, strand, gene_id)
    
    # return transcripts

    with open_file(gtf_path) as fh:
        for line in fh:
            if not line.startswith('#'):
            
                parts = line.strip().split('\t')
                if len(parts) >= 9:
                
                    chrom, source, feature, start, end, score, strand, frame, attributes = parts
                    
                    # Only process transcript features
                    if feature == feature_type:
                    
                        attrs = parse_gtf_attributes(attributes)
                        transcript_id = attrs.get('transcript_id')
                        gene_id = attrs.get('gene_id', 'unknown')
                        
                        if transcript_id:
                        
                            # GTF is 1-based, convert to 0-based for BED
                            start_pos = int(start) - 1
                            end_pos = int(end)
                            
                            transcripts[transcript_id] = (chrom, start_pos, end_pos, strand, gene_id)
    
    return transcripts


def create_tss_tts_beds(transcripts, output_prefix, deduplicate=True):
    """
    Create TSS and TTS BED files from transcript coordinates.
    
    For + strand: TSS = start, TTS = end
    For - strand: TSS = end, TTS = start
    
    Args:
        transcripts: dict of transcript_id -> (chrom, start, end, strand, gene_id)
        output_prefix: prefix for output files
        deduplicate: if True, keep only unique positions per strand
    """
    tss_coords = []  # (chrom, pos, strand, transcript_id)
    tts_coords = []
    
    for tx_id, (chrom, start, end, strand, gene_id) in transcripts.items():
        if strand == '+':
            tss_pos = start
            tts_pos = end
        elif strand == '-':
            tss_pos = end
            tts_pos = start
        else:
            print(f"Warning: Unknown strand '{strand}' for {tx_id}, skipping", file=sys.stderr)
            continue
        
        tss_coords.append((chrom, tss_pos, strand, tx_id, gene_id))
        tts_coords.append((chrom, tts_pos, strand, tx_id, gene_id))
    
    # Optionally deduplicate
    if deduplicate:
        # Keep unique (chrom, pos, strand) combinations
        tss_unique = {}
        for chrom, pos, strand, tx_id, gene_id in tss_coords:
            key = (chrom, pos, strand)
            if key not in tss_unique:
                tss_unique[key] = (chrom, pos, strand, tx_id, gene_id)
        tss_coords = list(tss_unique.values())
        
        tts_unique = {}
        for chrom, pos, strand, tx_id, gene_id in tts_coords:
            key = (chrom, pos, strand)
            if key not in tts_unique:
                tts_unique[key] = (chrom, pos, strand, tx_id, gene_id)
        tts_coords = list(tts_unique.values())
    
    # Sort by chrom, position
    tss_coords.sort(key=lambda x: (x[0], x[1]))
    tts_coords.sort(key=lambda x: (x[0], x[1]))
    
    # Write TSS BED (single position, 1bp window)
    tss_path = Path(f"{output_prefix}_tss.bed")
    with open(tss_path, 'w') as fh:
        for chrom, pos, strand, tx_id, gene_id in tss_coords:
            # BED format: chrom, start, end, name, score, strand
            fh.write(f"{chrom}\t{pos}\t{pos+1}\t{tx_id}\t0\t{strand}\n")
    
    # Write TTS BED (single position, 1bp window)
    tts_path = Path(f"{output_prefix}_tts.bed")
    with open(tts_path, 'w') as fh:
        for chrom, pos, strand, tx_id, gene_id in tts_coords:
            fh.write(f"{chrom}\t{pos}\t{pos+1}\t{tx_id}\t0\t{strand}\n")
    
    print(f"Processed {len(transcripts)} transcripts")
    print(f"TSS sites: {len(tss_coords)} ({'unique' if deduplicate else 'total'})")
    print(f"TTS sites: {len(tts_coords)} ({'unique' if deduplicate else 'total'})")
    print(f"Created: {tss_path}")
    print(f"Created: {tts_path}")
    
    return tss_path, tts_path


def create_gene_level_beds(transcripts, output_prefix):
    """
    Create gene-level TSS and TTS BED files.
    Aggregates all transcript starts/ends per gene.
    """
    gene_tss = defaultdict(set)  # gene_id -> set of (chrom, pos, strand)
    gene_tts = defaultdict(set)
    
    for tx_id, (chrom, start, end, strand, gene_id) in transcripts.items():
        if strand == '+':
            tss_pos = start
            tts_pos = end
        elif strand == '-':
            tss_pos = end
            tts_pos = start
        else:
            continue
        
        gene_tss[gene_id].add((chrom, tss_pos, strand))
        gene_tts[gene_id].add((chrom, tts_pos, strand))
    
    # Write gene-level TSS BED
    tss_coords = []
    for gene_id, positions in gene_tss.items():
        for chrom, pos, strand in positions:
            tss_coords.append((chrom, pos, strand, gene_id))
    tss_coords.sort(key=lambda x: (x[0], x[1]))
    
    tss_path = Path(f"{output_prefix}_gene_tss.bed")
    with open(tss_path, 'w') as fh:
        for chrom, pos, strand, gene_id in tss_coords:
            fh.write(f"{chrom}\t{pos}\t{pos+1}\t{gene_id}\t0\t{strand}\n")
    
    # Write gene-level TTS BED
    tts_coords = []
    for gene_id, positions in gene_tts.items():
        for chrom, pos, strand in positions:
            tts_coords.append((chrom, pos, strand, gene_id))
    tts_coords.sort(key=lambda x: (x[0], x[1]))
    
    tts_path = Path(f"{output_prefix}_gene_tts.bed")
    with open(tts_path, 'w') as fh:
        for chrom, pos, strand, gene_id in tts_coords:
            fh.write(f"{chrom}\t{pos}\t{pos+1}\t{gene_id}\t0\t{strand}\n")
    
    print(f"\nGene-level aggregation:")
    print(f"Genes: {len(gene_tss)}")
    print(f"Gene TSS sites: {len(tss_coords)}")
    print(f"Gene TTS sites: {len(tts_coords)}")
    print(f"Created: {tss_path}")
    print(f"Created: {tts_path}")
    
    return tss_path, tts_path


def main():
    parser = argparse.ArgumentParser(
        description="Extract single-position TSS/TTS coordinates from GTF file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract transcript-level TSS/TTS
  %(prog)s --gtf gencode.gtf --output-prefix ref_transcript
  
  # Extract with deduplication (unique positions only)
  %(prog)s --gtf gencode.gtf --output-prefix ref_transcript --deduplicate
  
  # Extract gene-level aggregated TSS/TTS
  %(prog)s --gtf gencode.gtf --output-prefix ref --gene-level
  
  # Extract both transcript and gene level
  %(prog)s --gtf gencode.gtf --output-prefix ref --gene-level

Output files:
  {prefix}_tss.bed        - Transcript TSS positions (1bp windows)
  {prefix}_tts.bed        - Transcript TTS positions (1bp windows)
  {prefix}_gene_tss.bed   - Gene-level TSS positions (with --gene-level)
  {prefix}_gene_tts.bed   - Gene-level TTS positions (with --gene-level)
        """
    )
    
    parser.add_argument(
        '--gtf',
        required=True,
        type=Path,
        help='Input GTF file (can be gzipped)'
    )
    
    parser.add_argument(
        '--output-prefix',
        required=True,
        type=str,
        help='Prefix for output BED files'
    )
    
    parser.add_argument(
        '--feature-type',
        default='transcript',
        choices=['transcript', 'exon', 'CDS'],
        help='GTF feature type to extract (default: transcript)'
    )
    
    parser.add_argument(
        '--deduplicate',
        action='store_true',
        help='Keep only unique positions per strand (recommended)'
    )
    
    parser.add_argument(
        '--gene-level',
        action='store_true',
        help='Also create gene-level aggregated TSS/TTS files'
    )
    
    parser.add_argument(
        '--chrom-filter',
        type=str,
        help='Only process this chromosome (e.g., chr22)'
    )
    
    args = parser.parse_args()
    
    if not args.gtf.exists():
        print(f"Error: GTF file not found: {args.gtf}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Reading GTF file: {args.gtf}")
    print(f"Feature type: {args.feature_type}")
    
    # Extract transcripts
    transcripts = extract_transcripts_from_gtf(args.gtf, args.feature_type)
    
    if not transcripts:
        print("Warning: No transcripts found in GTF file - creating empty output files", file=sys.stderr)
        # Create empty output files
        tss_file = Path(f"{args.output_prefix}_tss.bed")
        tts_file = Path(f"{args.output_prefix}_tts.bed")
        tss_file.touch()
        tts_file.touch()
        print(f"Created empty TSS file: {tss_file}")
        print(f"Created empty TTS file: {tts_file}")
        return
    
    # Filter by chromosome if requested
    if args.chrom_filter:
        print(f"Filtering to chromosome: {args.chrom_filter}")
        transcripts = {
            tx_id: coords 
            for tx_id, coords in transcripts.items() 
            if coords[0] == args.chrom_filter
        }
        print(f"Transcripts after filter: {len(transcripts)}")
    
    # Create transcript-level TSS/TTS BEDs
    create_tss_tts_beds(transcripts, args.output_prefix, args.deduplicate)
    
    # Optionally create gene-level TSS/TTS BEDs
    if args.gene_level:
        create_gene_level_beds(transcripts, args.output_prefix)
    
    print("\nDone!")


if __name__ == '__main__':
    main()
