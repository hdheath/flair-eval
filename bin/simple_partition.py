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
    # Index BAM if needed
    bai_file = Path(str(input_bam) + ".bai")
    if not bai_file.exists():
        print(f"Indexing BAM file: {input_bam}")
        run_command(["samtools", "index", str(input_bam)])
    
    # Extract region
    print(f"Extracting region {region} from {input_bam}")
    run_command([
        "samtools", "view", "-b", 
        "-o", str(output_bam),
        str(input_bam), region
    ])





def parse_region(region_str):
    """Parse region string like 'chr1:1000-2000' into components."""
    if ':' not in region_str or '-' not in region_str:
        raise ValueError(f"Region must be in format 'chr:start-end', got: {region_str}")
    
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


def partition_bed_file(input_file, output_file, chrom, start, end):
    """Filter any BED-like file to only include entries in the specified region."""
    if not input_file.exists():
        print(f"Warning: Optional file not found, skipping: {input_file}")
        return False
        
    print(f"Filtering {input_file} for {chrom}:{start}-{end}")
    
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
            
            # Check if this entry overlaps with our region
            if (bed_chrom == chrom and 
                not (bed_end < start or bed_start > end)):
                outfile.write(line)
    
    print(f"Created: {output_file}")
    return True


def main():
    parser = argparse.ArgumentParser(description="Simple partition script for FLAIR BAM/BED files")
    parser.add_argument("--bam", required=True, help="Input BAM file from FLAIR align")
    parser.add_argument("--bed", required=True, help="Input BED file from FLAIR align") 
    parser.add_argument("--region", required=True, help="Region to extract (e.g., chr1:1000-2000)")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output files")
    
    # Optional files
    parser.add_argument("--gtf", help="GTF annotation file (optional)")
    parser.add_argument("--genome", help="Genome FASTA file (optional)")
    parser.add_argument("--cage-peaks", help="CAGE peaks BED file (optional)")
    parser.add_argument("--quantseq-peaks", help="QuantSeq peaks BED file (optional)")
    parser.add_argument("--junctions", help="Junction file (STAR SJ.out.tab or BED format)")
    parser.add_argument("--target-regions", help="Target regions BED file (optional)")
    
    args = parser.parse_args()
    
    # Parse inputs
    input_bam = Path(args.bam)
    input_bed = Path(args.bed)
    chrom, start, end = parse_region(args.region)
    
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
    
    print(f"Partitioning to region: {chrom}:{start}-{end}")
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
        if file_path:
            input_path = Path(file_path)
            
            if file_type == 'genome' and input_path.exists():
                # For genome FASTA, extract the sequence for this region
                output_genome = Path(f"{args.output_prefix}_genome.fa")
                print(f"Extracting genome sequence for {args.region}")
                try:
                    run_command([
                        "samtools", "faidx", str(input_path), args.region,
                        "-o", str(output_genome)
                    ])
                    created_files.append(str(output_genome))
                except:
                    print(f"Warning: Could not extract genome sequence (samtools faidx failed)")
                    
            elif file_type == 'gtf' and input_path.exists():
                # For GTF, filter annotations in this region
                output_gtf = Path(f"{args.output_prefix}_annotation.gtf") 
                if partition_bed_file(input_path, output_gtf, chrom, start, end):
                    created_files.append(str(output_gtf))
                    
            elif input_path.exists():
                # For BED-like files (CAGE, QuantSeq, junctions, targets)
                output_file = Path(f"{args.output_prefix}_{file_type}.bed")
                if partition_bed_file(input_path, output_file, chrom, start, end):
                    created_files.append(str(output_file))
    
    print("\nPartitioning complete!")
    print("Created files:")
    for file_path in created_files:
        print(f"  {file_path}")


if __name__ == "__main__":
    main()