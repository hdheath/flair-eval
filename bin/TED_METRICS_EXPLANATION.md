# TED Option 4: Comprehensive Read Assignment Metrics

## Overview
The TED script now implements **Option 4** for tracking read assignments with multiple denominators and detailed breakdowns of primary vs supplementary alignment usage.

## New Metrics Reported

### 1. Assignment Counts
- **`assigned_unique_read_ids`**: Number of unique read IDs assigned to isoforms (from `*.isoform.read.map.txt`)
  - Represents the number of **molecules** (unique reads) that were assigned
  - Each read ID counted once, even if it has multiple alignments

- **`assigned_primary_alignments`**: Number of primary alignments for assigned reads
  - Extracted from BAM file using read IDs from map file
  - Excludes secondary and unmapped reads

- **`assigned_supplementary_alignments`**: Number of supplementary alignments for assigned reads
  - Includes alignments marked with supplementary flag (0x800)
  - Useful for tracking fusion transcripts and complex splicing

- **`assigned_total_alignments`**: Total alignments used (primary + supplementary)
  - Shows all alignment records that contributed to isoform discovery

### 2. Input Totals
- **`input_molecules`**: Number of input molecules to the process
  - For **collapse**: line count from corrected BED (`*_all_corrected.bed`)
  - For **transcriptome**: primary alignments in input BAM
  
- **`input_primary_alignments`**: Total primary alignments in input BAM
  - Counted with `samtools view -c -F 2308` (excludes unmapped, secondary, supplementary)

- **`input_supplementary_alignments`**: Total supplementary alignments in input BAM
  - Counted with `samtools view -c -f 2048 -F 260` (includes only supplementary)

- **`input_total_alignments`**: Total alignments (primary + supplementary)

### 3. Assignment Percentages

#### **`assignment_rate`** (Primary Metric)
```
assignment_rate = assigned_unique_read_ids / input_molecules
```
- **Interpretation**: What percentage of input molecules were successfully assigned to isoforms?
- **Denominator**:
  - Collapse: corrected reads (filtered input)
  - Transcriptome: primary alignments (all aligned reads)
- **Use case**: Primary QC metric for isoform reconstruction success

#### **`primary_alignment_utilization`**
```
primary_alignment_utilization = assigned_primary_alignments / input_primary_alignments
```
- **Interpretation**: What percentage of primary alignments contributed to isoforms?
- **Always uses primary alignments as denominator**
- **Use case**: Understanding alignment usage efficiency

#### **`total_alignment_utilization`**
```
total_alignment_utilization = assigned_total_alignments / input_total_alignments
```
- **Interpretation**: What percentage of ALL alignments (including supplementary) were used?
- **Use case**: Detecting fusion transcript contribution, complex splicing patterns

## Why Multiple Metrics?

### Different Questions, Different Metrics

1. **"How well did isoform reconstruction work?"**
   → Use `assignment_rate` (molecules in → molecules assigned)

2. **"Are supplementary alignments helping?"**
   → Compare `primary_alignment_utilization` vs `total_alignment_utilization`
   - If similar: supplementary alignments not contributing much
   - If different: supplementary alignments important for complex isoforms

3. **"How many fusion/chimeric transcripts?"**
   → Look at `assigned_supplementary_alignments`
   - High counts suggest complex transcript structures

## Implementation Details

### New Functions Added

#### `_get_assigned_read_ids(map_path: Path) -> set`
- Extracts unique read IDs from `*.isoform.read.map.txt`
- Returns set of read IDs for downstream counting

#### `_count_assigned_reads_by_type(bam_path: Path, assigned_read_ids: set) -> dict`
- Iterates through BAM file
- Counts primary vs supplementary alignments for assigned reads only
- Requires `pysam`

#### `_count_total_alignments_bam(bam: Path) -> dict`
- Uses `samtools view` to count alignments
- Separates primary and supplementary counts
- Fast (no Python BAM iteration)

### BAM Flag Reference
- `0x4` (4): Unmapped
- `0x100` (256): Secondary alignment
- `0x800` (2048): Supplementary alignment
- `2308` = 4 + 256 + 2048 (exclude all non-primary)

## Output Format

### TED.tsv Columns (Relevant to Assignment)
```
assigned_unique_read_ids                 | 850
assigned_primary_alignments              | 900
assigned_supplementary_alignments        | 50
assigned_total_alignments                | 950
input_molecules                          | 1000
input_primary_alignments                 | 1000
input_supplementary_alignments           | 80
input_total_alignments                   | 1080
assignment_rate                          | 0.85 (85%)
primary_alignment_utilization            | 0.90 (90%)
total_alignment_utilization              | 0.88 (88%)
```

### Interpretation Example
- 85% of input molecules assigned to isoforms ✓
- 90% of primary alignments used (good efficiency)
- 88% total alignment utilization (supplementary reads contributing)
- 50 supplementary alignments used (5.3% of assigned alignments)
  → Some complex isoforms with split alignments

## Performance Considerations

- **Fast**: `samtools view -c` for totals (no Python iteration)
- **Moderate**: Reading map file (text parsing)
- **Slower**: BAM iteration for assigned read counting (requires `pysam`)
  - Only iterates once per BAM
  - Skips reads not in assigned set (early filtering)

## Backward Compatibility

**Breaking Change**: Column names changed
- Old: `assigned_primary_reads`, `align_primary_total`, `assigned_pct`
- New: `assigned_unique_read_ids`, `input_molecules`, `assignment_rate`, etc.

Scripts parsing TED.tsv will need updates to use new column names.

## Future Enhancements

Possible additions:
1. Per-gene assignment rates
2. Multi-mapping read analysis
3. Strandedness metrics
4. Insert size distributions for paired-end data
5. Read quality stratification

## References

- FLAIR documentation: https://github.com/BrooksLabUCSC/flair
- SAM flag specification: https://samtools.github.io/hts-specs/SAMv1.pdf
- Supplementary alignments in long-read data: https://lh3.github.io/minimap2/minimap2.html
