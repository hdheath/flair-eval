#!/usr/bin/env python3
"""
Synthesize TED and FLAIR evaluation results into unified summary tables.

This script combines TED metrics (.tsv files) and FLAIR evaluation results (.tsv files)
into comprehensive summary tables with consistent formatting.
"""

import argparse
import csv
from pathlib import Path
from collections import OrderedDict

# Metadata fields to exclude from metric columns
METADATA_FIELDS = {'test_name', 'dataset', 'align_mode', 'partition_mode', 'stage', 'pipeline_mode'}

# Minimum lines required for a valid TSV file (header + at least one data row)
MIN_TSV_LINES = 2


def parse_filename_metadata(filename):
    """
    Extract metadata from standardized filename patterns.
    
    Expected patterns:
    - TED: dataset_align_partition_process_stage_ted.tsv
    - FLAIR: dataset_align_partition_process_stage_flair_eval.tsv
    
    Example: a549_chr5_default_all_transcriptome_with_gtf_transcriptome_ted.tsv
    -> dataset=a549_chr5, align=default, partition=all, process=transcriptome_with_gtf
    """
    base_name = Path(filename).stem
    
    # Remove common suffixes (order matters - do flair_eval before ted to avoid partial match)
    base_name = base_name.replace('_flair_eval', '').replace('_ted', '').replace('.eval_summary', '')
    
    # Split by underscores
    parts = base_name.split('_')
    
    # Dataset is always first part + second part if second part starts with 'chr'
    # e.g., a549_chr5 or a549_chr1 or a549_chr22
    if len(parts) > 1 and parts[1].startswith('chr'):
        dataset = f"{parts[0]}_{parts[1]}"
        remaining_parts = parts[2:]
    else:
        dataset = parts[0]
        remaining_parts = parts[1:]
    
    # Next part is align_mode
    align_mode = remaining_parts[0] if len(remaining_parts) > 0 else 'unknown'
    
    # Next part is partition_mode  
    partition_mode = remaining_parts[1] if len(remaining_parts) > 1 else 'unknown'
    
    # Remaining parts define the pipeline stage and mode
    process_parts = remaining_parts[2:] if len(remaining_parts) > 2 else []
    
    # Determine pipeline_stage and pipeline_mode from the remaining parts
    if 'transcriptome' in process_parts:
        pipeline_stage = 'transcriptome'
        # Remove 'transcriptome' and join the rest as mode
        mode_parts = [p for p in process_parts if p != 'transcriptome']
        pipeline_mode = '_'.join(mode_parts) if mode_parts else 'default'
    elif 'collapse' in process_parts:
        pipeline_stage = 'collapse'
        # Format: collapse_correct_mode_collapse_mode_collapse
        # e.g., ['collapse', 'with', 'gtf', 'default', 'collapse']
        # Remove first and last 'collapse' to get the modes
        collapse_indices = [i for i, p in enumerate(process_parts) if p == 'collapse']
        
        if len(collapse_indices) >= 2:
            # Standard format with two 'collapse' markers
            first_idx = collapse_indices[0]
            last_idx = collapse_indices[-1]
            
            # Everything between first and last collapse
            middle_parts = process_parts[first_idx+1:last_idx]
            
            # The last element before final 'collapse' is the collapse mode
            # Everything else is the correct mode
            if middle_parts:
                if len(middle_parts) == 1:
                    # Just collapse mode, use default for correct
                    correct_mode = 'default'
                    collapse_mode = middle_parts[0]
                else:
                    # Last part is collapse mode, rest is correct mode
                    collapse_mode = middle_parts[-1]
                    correct_mode = '_'.join(middle_parts[:-1])
                pipeline_mode = f"{correct_mode}+{collapse_mode}"
            else:
                pipeline_mode = 'default'
        elif len(collapse_indices) == 1:
            # Single collapse marker - everything before it is the mode
            collapse_idx = collapse_indices[0]
            if collapse_idx > 0:
                pipeline_mode = '_'.join(process_parts[:collapse_idx])
            else:
                # Just the stage, get everything after collapse
                remaining = process_parts[collapse_idx+1:]
                pipeline_mode = '_'.join(remaining) if remaining else 'default'
        else:
            # No collapse marker found (shouldn't happen)
            pipeline_mode = '_'.join(process_parts)
    else:
        pipeline_stage = 'unknown'
        pipeline_mode = '_'.join(process_parts) if process_parts else 'unknown'
    
    metadata = {
        'dataset': dataset,
        'align_mode': align_mode,
        'partition_mode': partition_mode,
        'stage': pipeline_stage,
        'pipeline_mode': pipeline_mode,
        'source_file': filename
    }
    
    return metadata


def parse_tsv_file(filepath):
    """
    Parse evaluation TSV file (TED or FLAIR) and return list of row dictionaries.
    
    Args:
        filepath: Path object pointing to TSV file
        
    Returns:
        List of OrderedDict objects, one per data row with metadata prepended
    """
    try:
        results = []
        metadata = parse_filename_metadata(filepath.name)
        
        with open(filepath, 'r') as f:
            # Use csv.DictReader for more robust TSV parsing
            reader = csv.DictReader(f, delimiter='\t')
            
            for row_dict in reader:
                row = OrderedDict()
                
                # Add metadata first (exclude source_file)
                for key, value in metadata.items():
                    if key != 'source_file':
                        row[key] = value
                
                # Add metric columns
                for key, value in row_dict.items():
                    row[key] = value
                
                results.append(row)
        
        if not results:
            print(f"Warning: No data rows found in {filepath}")
        
        return results
        
    except FileNotFoundError:
        print(f"Error: File not found: {filepath}")
        return []
    except csv.Error as e:
        print(f"Error parsing CSV/TSV file {filepath}: {e}")
        return []
    except Exception as e:
        print(f"Unexpected error parsing {filepath}: {e}")
        return []


def process_evaluation_files(file_list, file_type, test_name, results_by_sample):
    """
    Process a list of evaluation files and merge into results dictionary.
    
    Args:
        file_list: List of file paths to process
        file_type: String description for logging (e.g., 'TED', 'FLAIR')
        test_name: Test set name for labeling
        results_by_sample: Dictionary to populate with results, keyed by sample tuple
    """
    print(f"Processing {len(file_list)} {file_type} files...")
    
    for filepath_str in file_list:
        filepath = Path(filepath_str)
        if not filepath.exists():
            print(f"Warning: {file_type} file not found: {filepath_str}")
            continue
        
        rows = parse_tsv_file(filepath)
        for row in rows:
            # Create sample key from metadata
            sample_key = (row['dataset'], row['align_mode'], row['partition_mode'], 
                         row['stage'], row['pipeline_mode'])
            
            # Initialize sample if not exists
            if sample_key not in results_by_sample:
                results_by_sample[sample_key] = OrderedDict()
                results_by_sample[sample_key]['test_name'] = test_name
                results_by_sample[sample_key]['dataset'] = row['dataset']
                results_by_sample[sample_key]['align_mode'] = row['align_mode']
                results_by_sample[sample_key]['partition_mode'] = row['partition_mode']
                results_by_sample[sample_key]['stage'] = row['stage']
                results_by_sample[sample_key]['pipeline_mode'] = row['pipeline_mode']
            
            # Add metrics (exclude metadata fields)
            for key, value in row.items():
                if key not in METADATA_FIELDS:
                    results_by_sample[sample_key][key] = value


def main():
    parser = argparse.ArgumentParser(description="Synthesize TED and FLAIR evaluation results")
    parser.add_argument('--ted-files', nargs='+', required=True, 
                       help='TED evaluation TSV files')
    parser.add_argument('--flair-files', nargs='+', required=True,
                       help='FLAIR evaluation TSV files')
    parser.add_argument('--output', required=True,
                       help='Output TSV file for unified results')
    parser.add_argument('--test-name', required=True,
                       help='Test set name for labeling')
    
    args = parser.parse_args()
    
    # Dictionary to hold results keyed by sample identifier
    results_by_sample = {}
    
    # Process both TED and FLAIR files using unified function
    process_evaluation_files(args.ted_files, 'TED', args.test_name, results_by_sample)
    process_evaluation_files(args.flair_files, 'FLAIR', args.test_name, results_by_sample)
    
    if not results_by_sample:
        print("No valid results found!")
        return 1
    
    # Convert to list of rows
    all_results = list(results_by_sample.values())
    
    # Write unified TSV
    success = write_tsv_merged(all_results, args.output, args.test_name)
    
    if success:
        print(f"Saved unified evaluation summary to {args.output}")
        print(f"Total merged samples: {len(all_results)}")
        return 0
    else:
        print("Failed to write output file")
        return 1


def write_tsv_merged(rows, output_path, test_name):
    """
    Write merged rows to TSV file with TED and FLAIR metrics combined.
    
    Args:
        rows: List of OrderedDict objects containing merged data
        output_path: Path to output TSV file
        test_name: Test set name (unused but kept for API compatibility)
        
    Returns:
        True if successful, False otherwise
    """
    if not rows:
        print("No rows to write!")
        return False
    
    try:
        # Collect all columns present in data
        all_columns = OrderedDict()
        for row in rows:
            for col in row.keys():
                all_columns[col] = None
        
        # Define metadata column order
        metadata_cols = ['test_name', 'dataset', 'align_mode', 'partition_mode', 
                         'stage', 'pipeline_mode']
        
        # Define known TED metrics (from TED TSV files) in logical grouping order
        known_ted_metrics = [
            # Isoform counts
            'isoforms_observed', 'genes_observed',
            # Read assignment counts
            'assigned_unique_read_ids', 'assigned_primary_alignments', 
            'assigned_supplementary_alignments', 'assigned_total_alignments',
            # Per-isoform statistics
            'reads_per_isoform_mean', 'reads_per_isoform_median', 
            'reads_per_isoform_min', 'reads_per_isoform_max',
            # Input counts
            'input_molecules', 'input_primary_alignments', 
            'input_supplementary_alignments', 'input_total_alignments',
            # Utilization rates
            'assignment_rate', 'primary_alignment_utilization', 'total_alignment_utilization',
            # 5' end metrics
            '5prime_precision', '5prime_recall', '5prime_f1',
            # 3' end metrics
            '3prime_precision', '3prime_recall', '3prime_f1',
            # Reference 5' end metrics
            'ref5prime_precision', 'ref5prime_recall', 'ref5prime_f1',
            # Reference 3' end metrics
            'ref3prime_precision', 'ref3prime_recall', 'ref3prime_f1'
        ]
        
        # Define known FLAIR metrics (from FLAIR eval files)
        known_flair_metrics = [
            # Genic region metrics
            'total_read_regions', 'found_regions', 'genic_reads',
            # Splice junction metrics
            'total_sjc', 'supported_sjc', 'subset_sjc', 'total_se', 'supported_se',
            # SQANTI-like classification
            'FSM', 'ISM', 'NIC', 'NNC', 'SEM', 'SEN'
        ]
        
        # Filter to only include columns that exist in the data
        ted_cols = [c for c in known_ted_metrics if c in all_columns]
        flair_cols = [c for c in known_flair_metrics if c in all_columns]
        
        # Identify any unexpected columns not in known lists
        known_all = set(metadata_cols + known_ted_metrics + known_flair_metrics)
        other_cols = sorted([c for c in all_columns.keys() if c not in known_all])
        
        # Create final column order: metadata, TED metrics, FLAIR metrics, other
        final_columns = metadata_cols + ted_cols + flair_cols + other_cols
        
        # Write TSV using csv module for proper escaping
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            
            # Write header
            writer.writerow(final_columns)
            
            # Write data rows
            for row in rows:
                values = [row.get(col, '') for col in final_columns]
                writer.writerow(values)
        
        return True
        
    except IOError as e:
        print(f"Error writing output file {output_path}: {e}")
        return False
    except Exception as e:
        print(f"Unexpected error writing output: {e}")
        return False


if __name__ == '__main__':
    exit(main())