#!/usr/bin/env python3
"""
Simple partition script for FLAIR outputs.
Just extracts reads from a specific genomic region - no complex metadata tracking.
"""

import argparse
import subprocess
import sys
from pathlib import Path


def run_command(cmd):
    """Run a shell command and check for errors."""
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    return result.stdout


def partition_bam(input_bam, output_bam, region):
    """Extract reads from BAM file for specified region."""
    # Index input BAM if needed
    bai_file = Path(str(input_bam) + ".bai")
    if not bai_file.exists():
        print(f"Indexing input BAM file: {input_bam}")
        run_command(["samtools", "index", str(input_bam)])
    
    # Extract region - samtools accepts both 'chr1' and 'chr1:1000-2000' formats
    print(f"Extracting region {region} from {input_bam}")
    run_command([
        "samtools", "view", "-b", 
        "-o", str(output_bam),
        str(input_bam), region
    ])
    
    # Index the output BAM file (required for FLAIR transcriptome)
    print(f"Indexing output BAM file: {output_bam}")
    run_command(["samtools", "index", str(output_bam)])





def parse_region(region_str):
    """Parse region string like 'chr1:1000-2000' or just 'chr1' into components."""
    # Handle case where only chromosome is specified (e.g., 'chr1')
    if ':' not in region_str:
        # Just chromosome name - return None for start/end to indicate whole chromosome
        return region_str, None, None
    
    # Handle case with coordinates (e.g., 'chr1:1000-2000')
    if '-' not in region_str:
        raise ValueError(f"Region with coordinates must be in format 'chr:start-end', got: {region_str}")
    
    chrom, pos_range = region_str.split(':', 1)
    start_str, end_str = pos_range.split('-', 1)
    
    try:
        start = int(start_str)
        end = int(end_str)
    except ValueError:
        raise ValueError(f"Start and end positions must be integers: {region_str}")
    
    if start > end:
        start, end = end, start  # Swap if backwards
    
    return chrom, start, end


def partition_bed_file(input_file, output_file, chrom, start, end, create_empty_if_missing=False):
    """Filter any BED-like file to only include entries in the specified region.

    Args:
        input_file: Path to input BED file
        output_file: Path to output BED file
        chrom: Chromosome to filter to
        start: Start position (0-based, inclusive)
        end: End position (0-based, exclusive)
        create_empty_if_missing: If True, create empty output file even if input doesn't exist

    Returns:
        True if output file was created, False otherwise
    """
    if not input_file.exists():
        if create_empty_if_missing:
            print(f"Input file not found, creating empty output: {input_file} -> {output_file}")
            # Create empty file
            output_file.touch()
            return True
        else:
            print(f"Warning: Optional file not found, skipping: {input_file}")
            return False
    
    # Format region string for logging
    if start is None or end is None:
        region_str = chrom  # Whole chromosome
    else:
        region_str = f"{chrom}:{start}-{end}"
        
    print(f"Filtering {input_file} for {region_str}")
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 3:
                outfile.write(line)  # Keep non-standard lines
                continue
                
            bed_chrom = fields[0]
            try:
                bed_start = int(fields[1])
                bed_end = int(fields[2])
            except ValueError:
                outfile.write(line)  # Keep lines with non-numeric coordinates
                continue
            
            # Check if this entry is in the target chromosome
            if bed_chrom == chrom:
                # If no coordinates specified, include all entries from this chromosome
                if start is None or end is None:
                    outfile.write(line)
                # Otherwise check if fully contained within specified region
                elif bed_start >= start and bed_end <= end:
                    outfile.write(line)
    
    print(f"Created: {output_file}")
    return True


def partition_gtf_file(input_file, output_file, chrom, start, end):
    """Filter GTF file to only include transcripts fully contained in the region."""
    if not input_file.exists():
        print(f"Warning: GTF file not found, skipping: {input_file}")
        return False
    
    # Format region string for logging
    if start is None or end is None:
        region_str = chrom  # Whole chromosome
    else:
        region_str = f"{chrom}:{start}-{end}"
        
    print(f"Filtering GTF {input_file} for {region_str}")
    
    # First pass: identify which transcripts are fully contained in the region
    valid_transcripts = set()
    transcript_bounds = {}  # transcript_id -> (start, end)
    
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            gtf_chrom = fields[0]
            feature_type = fields[2]
            gtf_start = int(fields[3])  # GTF is 1-based
            gtf_end = int(fields[4])    # GTF end is inclusive
            attributes = fields[8]
            
            # Only look at chromosome of interest
            if gtf_chrom != chrom:
                continue
            
            # Extract transcript_id from attributes
            transcript_id = None
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    transcript_id = attr.split('"')[1]
                    break
            
            if not transcript_id:
                continue
            
            # Track transcript boundaries (using actual coordinates, not 0-based)
            if transcript_id not in transcript_bounds:
                transcript_bounds[transcript_id] = [gtf_start, gtf_end]
            else:
                transcript_bounds[transcript_id][0] = min(transcript_bounds[transcript_id][0], gtf_start)
                transcript_bounds[transcript_id][1] = max(transcript_bounds[transcript_id][1], gtf_end)
    
    # Determine which transcripts are fully contained
    if start is None or end is None:
        # No region filter - keep all transcripts on this chromosome
        valid_transcripts = set(transcript_bounds.keys())
    else:
        # GTF uses 1-based inclusive coordinates
        # Our region uses 0-based half-open (start is 0-based, end is exclusive in BED convention)
        # So we need: transcript_start >= start+1 AND transcript_end <= end
        for transcript_id, (trans_start, trans_end) in transcript_bounds.items():
            if trans_start >= start + 1 and trans_end <= end:
                valid_transcripts.add(transcript_id)
    
    print(f"Found {len(valid_transcripts)} transcripts fully contained in {region_str}")
    
    # Second pass: write only lines belonging to valid transcripts
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                outfile.write(line)  # Keep non-standard lines
                continue
            
            gtf_chrom = fields[0]
            attributes = fields[8]
            
            # Only process lines from target chromosome
            if gtf_chrom != chrom:
                continue
            
            # Extract transcript_id
            transcript_id = None
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    transcript_id = attr.split('"')[1]
                    break
            
            # Write line if transcript is valid
            if transcript_id and transcript_id in valid_transcripts:
                outfile.write(line)
    
    print(f"Created: {output_file}")
    return True


def main():
    parser = argparse.ArgumentParser(description="Simple partition script for FLAIR BAM/BED files")
    parser.add_argument("--bam", required=True, help="Input BAM file from FLAIR align")
    parser.add_argument("--bed", required=True, help="Input BED file from FLAIR align") 
    parser.add_argument("--region", help="Region to extract (e.g., chr1:1000-2000)")
    parser.add_argument("--all", action="store_true", help="Pass through all data (no filtering)")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output files")
    
    # Optional files
    parser.add_argument("--gtf", help="GTF annotation file (optional)")
    parser.add_argument("--genome", help="Genome FASTA file (optional)")
    parser.add_argument("--cage-peaks", help="CAGE peaks BED file (optional)")
    parser.add_argument("--quantseq-peaks", help="QuantSeq peaks BED file (optional)")
    parser.add_argument("--junctions", help="Junction file (STAR SJ.out.tab or BED format)")
    parser.add_argument("--target-regions", help="Target regions BED file (optional)")
    
    args = parser.parse_args()
    
    # Check arguments
    if not args.all and not args.region:
        print("Error: Must specify either --all or --region", file=sys.stderr)
        sys.exit(1)
    
    if args.all and args.region:
        print("Error: Cannot specify both --all and --region", file=sys.stderr)
        sys.exit(1)
    
    # Parse inputs
    input_bam = Path(args.bam)
    input_bed = Path(args.bed)
    
    if args.region:
        chrom, start, end = parse_region(args.region)
    else:
        chrom, start, end = None, None, None  # For --all mode
    
    # Check inputs exist
    if not input_bam.exists():
        print(f"Error: BAM file not found: {input_bam}", file=sys.stderr)
        sys.exit(1)
    if not input_bed.exists():
        print(f"Error: BED file not found: {input_bed}", file=sys.stderr)
        sys.exit(1)
    
    # Define outputs
    output_bam = Path(f"{args.output_prefix}.bam")
    output_bed = Path(f"{args.output_prefix}.bed")
    
    if args.all:
        print("Pass-through mode: creating symlinks to original files")
        print(f"Input BAM: {input_bam}")
        print(f"Input BED: {input_bed}")
        print(f"Output BAM: {output_bam}")
        print(f"Output BED: {output_bed}")
        
        # Create symlinks instead of copying
        output_bam.symlink_to(input_bam.resolve())
        output_bed.symlink_to(input_bed.resolve())
        
        # For BAM index - symlink if exists, create if doesn't
        input_bai = Path(str(input_bam) + ".bai")
        output_bai = Path(str(output_bam) + ".bai")
        if input_bai.exists():
            output_bai.symlink_to(input_bai.resolve())
        else:
            print(f"Creating BAM index for: {output_bam}")
            run_command(["samtools", "index", str(output_bam)])
    else:
        # Format region string for logging
        if start is None or end is None:
            region_display = chrom
        else:
            region_display = f"{chrom}:{start}-{end}"
        
        print(f"Partitioning to region: {region_display}")
        print(f"Input BAM: {input_bam}")
        print(f"Input BED: {input_bed}")
        print(f"Output BAM: {output_bam}")
        print(f"Output BED: {output_bed}")
        
        # Do the work - core BAM/BED files
        partition_bam(input_bam, output_bam, args.region)
        partition_bed_file(input_bed, output_bed, chrom, start, end)
    
    created_files = [str(output_bam), str(output_bed)]
    
    # Handle optional files
    optional_files = {
        'gtf': args.gtf,
        'genome': args.genome,
        'cage': args.cage_peaks,
        'quantseq': args.quantseq_peaks,
        'junctions': args.junctions,
        'targets': args.target_regions
    }

    for file_type, file_path in optional_files.items():
        # Special handling for CAGE and QuantSeq: always create output files (even if empty)
        # This satisfies Nextflow's output requirements while allowing evaluation to handle missing data
        if file_type in ['cage', 'quantseq'] and not file_path:
            output_file = Path(f"{args.output_prefix}_{file_type}.bed")
            print(f"Input {file_type} file not provided, creating empty output: {output_file}")
            output_file.touch()
            created_files.append(str(output_file))
            continue

        if file_path:
            input_path = Path(file_path)
            
            if args.all and input_path.exists():
                # Pass-through mode: create symlinks
                if file_type == 'genome':
                    output_file = Path(f"{args.output_prefix}_genome.fa")
                elif file_type == 'gtf':
                    output_file = Path(f"{args.output_prefix}_annotation.gtf")
                else:
                    output_file = Path(f"{args.output_prefix}_{file_type}.bed")
                
                output_file.symlink_to(input_path.resolve())
                created_files.append(str(output_file))
                
            elif file_type == 'genome' and input_path.exists():
                # For genome FASTA, create symlink to preserve coordinates
                # DO NOT extract sequence - FLAIR needs full genome for coordinate matching
                output_genome = Path(f"{args.output_prefix}_genome.fa")
                print(f"Creating symlink to full genome (preserving coordinates): {input_path}")
                output_genome.symlink_to(input_path.resolve())
                created_files.append(str(output_genome))
                    
            elif file_type == 'gtf' and input_path.exists():
                # For GTF, filter annotations to only include transcripts fully in this region
                output_gtf = Path(f"{args.output_prefix}_annotation.gtf") 
                if partition_gtf_file(input_path, output_gtf, chrom, start, end):
                    created_files.append(str(output_gtf))
                    
            elif input_path.exists():
                # For BED-like files (CAGE, QuantSeq, junctions, targets)
                output_file = Path(f"{args.output_prefix}_{file_type}.bed")
                # For CAGE and QuantSeq, always create output file (even if empty) to satisfy Nextflow
                # The evaluation script handles empty/missing files gracefully
                create_empty = (file_type in ['cage', 'quantseq'])
                if partition_bed_file(input_path, output_file, chrom, start, end, create_empty_if_missing=create_empty):
                    created_files.append(str(output_file))
            elif file_type in ['cage', 'quantseq']:
                # Input file doesn't exist but we need to create empty output for Nextflow
                output_file = Path(f"{args.output_prefix}_{file_type}.bed")
                print(f"Input {file_type} file not provided, creating empty output: {output_file}")
                output_file.touch()
                created_files.append(str(output_file))
    
    print("\nPartitioning complete!")
    print("Created files:")
    for file_path in created_files:
        print(f"  {file_path}")


if __name__ == "__main__":
    main()