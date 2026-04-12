#!/usr/bin/env python3
"""
Simple partition script for FLAIR outputs.
Just extracts reads from a specific genomic region - no complex metadata tracking.
"""

import argparse
import subprocess
import sys
from pathlib import Path

import pysam


def run_command(cmd):
    """Run a shell command and check for errors."""
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    return result.stdout


def partition_bam(input_bam, output_bam, region):
    """Extract reads from BAM file for specified region(s).

    Args:
        input_bam: path to input BAM
        output_bam: path to output BAM
        region: a single region string OR a list of region strings
    """
    # Index input BAM if needed
    bai_file = Path(str(input_bam) + ".bai")
    if not bai_file.exists():
        print(f"Indexing input BAM file: {input_bam}")
        run_command(["samtools", "index", str(input_bam)])

    regions = region if isinstance(region, list) else [region]

    if len(regions) == 1:
        # Single region: direct extraction
        print(f"Extracting region {regions[0]} from {input_bam}")
        run_command([
            "samtools", "view", "-b",
            "-o", str(output_bam),
            str(input_bam), regions[0]
        ])
    else:
        # Multiple regions: samtools view accepts multiple region args
        print(f"Extracting {len(regions)} regions from {input_bam}: {regions}")
        cmd = ["samtools", "view", "-b", "-o", str(output_bam), str(input_bam)] + regions
        run_command(cmd)

    # Index the output BAM file (required for FLAIR transcriptome)
    print(f"Indexing output BAM file: {output_bam}")
    run_command(["samtools", "index", str(output_bam)])


def bam_to_bed12(input_bam, output_bed):
    """Convert BAM file to BED12 format using pysam."""
    print(f"Converting BAM to BED12: {input_bam} -> {output_bed}")

    with pysam.AlignmentFile(str(input_bam), "rb") as bam, open(str(output_bed), "w") as out:
        for read in bam.fetch():
            if read.is_unmapped:
                continue
            chrom = read.reference_name
            start = read.reference_start
            end = read.reference_end
            name = read.query_name
            score = read.mapping_quality
            strand = "-" if read.is_reverse else "+"

            blocks = read.get_blocks()
            if not blocks:
                continue
            block_count = len(blocks)
            block_sizes = ",".join(str(b[1] - b[0]) for b in blocks)
            block_starts = ",".join(str(b[0] - start) for b in blocks)

            out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t"
                       f"{start}\t{end}\t0\t{block_count}\t{block_sizes}\t{block_starts}\n")

    print(f"Created BED12: {output_bed}")





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
    return _partition_gtf_multi(input_file, output_file, [(chrom, start, end)])


def _partition_bed_multi(input_file, output_file, parsed_regions):
    """Filter BED file to include entries in ANY of the specified regions.

    Args:
        input_file: Path to input BED file
        output_file: Path to output BED file
        parsed_regions: list of (chrom, start, end) tuples
    """
    input_file = Path(input_file)
    output_file = Path(output_file)
    if not input_file.exists():
        print(f"Warning: File not found, skipping: {input_file}")
        return False

    region_strs = []
    for chrom, start, end in parsed_regions:
        if start is None:
            region_strs.append(chrom)
        else:
            region_strs.append(f"{chrom}:{start}-{end}")
    print(f"Filtering {input_file} for regions: {region_strs}")

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            fields = line.strip().split('\t')
            if len(fields) < 3:
                outfile.write(line)
                continue
            bed_chrom = fields[0]
            try:
                bed_start = int(fields[1])
                bed_end = int(fields[2])
            except ValueError:
                outfile.write(line)
                continue

            for r_chrom, r_start, r_end in parsed_regions:
                if bed_chrom != r_chrom:
                    continue
                if r_start is None or r_end is None:
                    outfile.write(line)
                    break
                elif bed_start >= r_start and bed_end <= r_end:
                    outfile.write(line)
                    break

    print(f"Created: {output_file}")
    return True


def _partition_gtf_multi(input_file, output_file, parsed_regions):
    """Filter GTF to include transcripts fully contained in ANY of the specified regions.

    Args:
        input_file: Path to input GTF file
        output_file: Path to output GTF file
        parsed_regions: list of (chrom, start, end) tuples
    """
    input_file = Path(input_file)
    output_file = Path(output_file)
    if not input_file.exists():
        print(f"Warning: GTF file not found, skipping: {input_file}")
        return False

    region_strs = []
    for chrom, start, end in parsed_regions:
        if start is None:
            region_strs.append(chrom)
        else:
            region_strs.append(f"{chrom}:{start}-{end}")
    print(f"Filtering GTF {input_file} for regions: {region_strs}")

    # Build set of target chromosomes for quick filtering
    target_chroms = set(chrom for chrom, _, _ in parsed_regions)

    # First pass: identify which transcripts are fully contained in any region
    transcript_bounds = {}  # transcript_id -> (chrom, start, end)

    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            gtf_chrom = fields[0]
            if gtf_chrom not in target_chroms:
                continue
            gtf_start = int(fields[3])  # 1-based
            gtf_end = int(fields[4])    # inclusive
            attributes = fields[8]
            transcript_id = None
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    transcript_id = attr.split('"')[1]
                    break
            if not transcript_id:
                continue
            if transcript_id not in transcript_bounds:
                transcript_bounds[transcript_id] = [gtf_chrom, gtf_start, gtf_end]
            else:
                transcript_bounds[transcript_id][1] = min(transcript_bounds[transcript_id][1], gtf_start)
                transcript_bounds[transcript_id][2] = max(transcript_bounds[transcript_id][2], gtf_end)

    # Determine which transcripts are fully contained in at least one region
    valid_transcripts = set()
    for tid, (t_chrom, t_start, t_end) in transcript_bounds.items():
        for r_chrom, r_start, r_end in parsed_regions:
            if t_chrom != r_chrom:
                continue
            if r_start is None or r_end is None:
                valid_transcripts.add(tid)
                break
            elif t_start >= r_start + 1 and t_end <= r_end:
                valid_transcripts.add(tid)
                break

    print(f"Found {len(valid_transcripts)} transcripts fully contained in target regions")

    # Second pass: write valid transcript lines
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                outfile.write(line)
                continue
            gtf_chrom = fields[0]
            if gtf_chrom not in target_chroms:
                continue
            attributes = fields[8]
            transcript_id = None
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    transcript_id = attr.split('"')[1]
                    break
            if transcript_id and transcript_id in valid_transcripts:
                outfile.write(line)

    print(f"Created: {output_file}")
    return True


def main():
    parser = argparse.ArgumentParser(description="Simple partition script for FLAIR BAM/BED files")
    parser.add_argument("--bam", required=True, help="Input BAM file from FLAIR align")
    parser.add_argument("--bed", help="Input BED file from FLAIR align (optional if --generate-bed is used)")
    parser.add_argument("--generate-bed", action="store_true",
                        help="Generate BED12 from BAM after partitioning (saves time/storage for region partitions)")
    parser.add_argument("--region", nargs='+',
                        help="Region(s) to extract (e.g., chr1:1000-2000 chr3:48000000-53000000). "
                             "Multiple regions are merged into one output.")
    parser.add_argument("--all", action="store_true", help="Pass through all data (no filtering)")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output files")
    
    # Optional files
    parser.add_argument("--gtf", help="GTF annotation file (optional)")
    parser.add_argument("--genome", help="Genome FASTA file (optional)")
    parser.add_argument("--cage-peaks", help="CAGE peaks BED file (optional)")
    parser.add_argument("--drna-peaks", help="dRNA peaks BED file (optional)")
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
    
    # Validate BED input requirements
    if not args.generate_bed and not args.bed:
        print("Error: Must specify either --bed or --generate-bed", file=sys.stderr)
        sys.exit(1)

    # Parse inputs
    input_bam = Path(args.bam)
    input_bed = Path(args.bed) if args.bed else None

    # Parse regions — may be a list of 1+ regions
    parsed_regions = []  # list of (chrom, start, end) tuples
    if args.region:
        for r in args.region:
            parsed_regions.append(parse_region(r))
    # chrom/start/end shortcuts for single-region backward compat
    if len(parsed_regions) == 1:
        chrom, start, end = parsed_regions[0]
    elif len(parsed_regions) > 1:
        chrom, start, end = None, None, None  # multi-region mode
    else:
        chrom, start, end = None, None, None  # --all mode

    # Check inputs exist
    if not input_bam.exists():
        print(f"Error: BAM file not found: {input_bam}", file=sys.stderr)
        sys.exit(1)
    if input_bed and not input_bed.exists():
        print(f"Error: BED file not found: {input_bed}", file=sys.stderr)
        sys.exit(1)
    
    # Define outputs
    output_bam = Path(f"{args.output_prefix}.bam")
    output_bed = Path(f"{args.output_prefix}.bed")
    
    if args.all:
        print("Pass-through mode: creating symlinks to original files")
        print(f"Input BAM: {input_bam}")
        print(f"Output BAM: {output_bam}")
        print(f"Output BED: {output_bed}")

        # Create symlink for BAM
        output_bam.symlink_to(input_bam.resolve())

        # For BAM index - symlink if exists, create if doesn't
        input_bai = Path(str(input_bam) + ".bai")
        output_bai = Path(str(output_bam) + ".bai")
        if input_bai.exists():
            output_bai.symlink_to(input_bai.resolve())
        else:
            print(f"Creating BAM index for: {output_bam}")
            run_command(["samtools", "index", str(output_bam)])

        # Handle BED: either symlink existing or generate from BAM
        if args.generate_bed:
            print("Generating BED12 from BAM (--generate-bed mode)")
            bam_to_bed12(input_bam, output_bed)
        elif input_bed:
            print(f"Input BED: {input_bed}")
            output_bed.symlink_to(input_bed.resolve())
        else:
            print("Error: No BED source available in --all mode", file=sys.stderr)
            sys.exit(1)
    else:
        print(f"Partitioning to region(s): {[r for r in args.region]}")
        print(f"Input BAM: {input_bam}")
        print(f"Output BAM: {output_bam}")
        print(f"Output BED: {output_bed}")

        # Partition BAM first (fast with samtools — supports multiple regions natively)
        partition_bam(input_bam, output_bam, args.region)

        # Handle BED: either partition existing or generate from partitioned BAM
        if args.generate_bed:
            # Convert the partitioned BAM to BED12 (more efficient than partitioning full BED)
            print("Generating BED12 from partitioned BAM (--generate-bed mode)")
            bam_to_bed12(output_bam, output_bed)
        elif input_bed:
            print(f"Input BED: {input_bed}")
            _partition_bed_multi(input_bed, output_bed, parsed_regions)
        else:
            print("Error: No BED source available", file=sys.stderr)
            sys.exit(1)
    
    created_files = [str(output_bam), str(output_bed)]
    
    # Handle optional files
    optional_files = {
        'gtf': args.gtf,
        'genome': args.genome,
        'cage': args.cage_peaks,
        'drna': args.drna_peaks,
        'junctions': args.junctions,
        'targets': args.target_regions
    }

    for file_type, file_path in optional_files.items():
        # Special handling for CAGE and dRNA: always create output files (even if empty)
        # This satisfies Nextflow's output requirements while allowing evaluation to handle missing data
        if file_type in ['cage', 'drna'] and not file_path:
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
                # For GTF, filter annotations to only include transcripts fully in target region(s)
                output_gtf = Path(f"{args.output_prefix}_annotation.gtf") 
                if _partition_gtf_multi(input_path, output_gtf, parsed_regions):
                    created_files.append(str(output_gtf))
                    
            elif input_path.exists():
                # For BED-like files (CAGE, dRNA, junctions, targets)
                output_file = Path(f"{args.output_prefix}_{file_type}.bed")
                # For CAGE and dRNA, always create output file (even if empty) to satisfy Nextflow
                # The evaluation script handles empty/missing files gracefully
                create_empty = (file_type in ['cage', 'drna'])
                if _partition_bed_multi(input_path, output_file, parsed_regions):
                    created_files.append(str(output_file))
                elif create_empty:
                    output_file.touch()
                    created_files.append(str(output_file))
            elif file_type in ['cage', 'drna']:
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