# TED Standalone - Integration Guide

## What Was Removed

### 1. Pipeline-Specific Plumbing
- ❌ `from . import register, write_metrics` - Pipeline registration system
- ❌ `@register("TED")` decorators - Auto-discovery for other pipeline
- ❌ `PathBuilder` dependencies - Custom path management
- ❌ `upstreams` parameter - Inter-stage dependencies
- ❌ Config object parsing (`_cfg_get`, config traversal)
- ❌ Auto-discovery of files (`_dirs_with_file`, `_find_*` functions)
- ❌ `_build_region_metrics_index` - Region metadata indexing
- ❌ Sidecar metrics JSON generation
- ❌ Audit trail TSV generation
- ❌ Visualization generation (scatter plots, transcriptome browser)

### 2. Regionalization Logic
- ❌ Region pattern detection (`_is_region_tag`)
- ❌ Multi-region iteration
- ❌ Sliced peak file discovery
- ❌ Expected gene/isoform counts from region index

### 3. Complex File Discovery
- ❌ Searching through output directories for corrected BED
- ❌ Searching for align BAM files
- ❌ Dynamic resolution of peak files

## What Was Kept

### Core Functionality ✅
- TSS/TTS precision/recall/F1 calculation
- Read assignment metrics (Option 4 implementation)
- Primary vs supplementary alignment tracking
- BED file parsing and sorting
- bedtools integration
- samtools integration
- All metric calculations

### Simplified To
- **Single function**: `calculate_ted_metrics()` - does all calculations
- **Command-line interface**: Takes explicit file paths as arguments
- **One sample at a time**: Nextflow handles parallelization
- **No config files**: All parameters as CLI args

## How to Use in Nextflow

### Option A: As Python Script (Recommended)

```groovy
process FlairTED {
    publishDir "results/ted/${test_name}", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    tuple val(test_name), val(dataset_name), val(align_label), val(partition_label), val(process_label),
          path(isoforms_bed), path(isoform_map), path(bam), path(corrected_bed),
          path(cage_peaks), path(quantseq_peaks), path(ref_tss_peaks), path(ref_tts_peaks)
    
    output:
    tuple val(test_name), val(dataset_name), val(align_label), val(partition_label), val(process_label),
          path("${dataset_name}_${align_label}_${partition_label}_${process_label}.ted.tsv"), emit: ted_results
    
    script:
    def stage = process_label.contains("collapse") ? "collapse" : "transcriptome"
    def bam_arg = bam.name != "NO_BAM" ? "--bam ${bam}" : ""
    def corrected_arg = corrected_bed.name != "NO_CORRECTED" ? "--corrected-bed ${corrected_bed}" : ""
    def cage_arg = cage_peaks.name != "NO_CAGE" ? "--prime5-peaks ${cage_peaks}" : ""
    def qs_arg = quantseq_peaks.name != "NO_QUANTSEQ" ? "--prime3-peaks ${quantseq_peaks}" : ""
    def ref_tss_arg = ref_tss_peaks.name != "NO_REF_TSS" ? "--ref-prime5-peaks ${ref_tss_peaks}" : ""
    def ref_tts_arg = ref_tts_peaks.name != "NO_REF_TTS" ? "--ref-prime3-peaks ${ref_tts_peaks}" : ""
    
    """
    python3 ${projectDir}/bin/ted_standalone.py \\
        --isoforms-bed ${isoforms_bed} \\
        --map-file ${isoform_map} \\
        ${bam_arg} \\
        ${corrected_arg} \\
        ${cage_arg} \\
        ${qs_arg} \\
        ${ref_tss_arg} \\
        ${ref_tts_arg} \\
        --window 50 \\
        --stage ${stage} \\
        --output ${dataset_name}_${align_label}_${partition_label}_${process_label}.ted.tsv
    """
}
```

### Option B: Embedded in Process (Like FlairEval)

You could also embed the `calculate_ted_metrics()` function directly into a Nextflow process script block (similar to how FlairEval works), but Option A is cleaner.

## File Requirements

### Always Required
- `--isoforms-bed`: `*.isoforms.bed` file
- `--map-file`: `*.isoform.read.map.txt` file
- `--output`: Output TSV path

### Optional (for complete metrics)
- `--bam`: BAM file (for alignment utilization metrics)
- `--corrected-bed`: Corrected BED (for collapse stage denominator)
- `--prime5-peaks`: CAGE peaks (for 5' precision/recall)
- `--prime3-peaks`: QuantSeq/PolyA peaks (for 3' precision/recall)
- `--ref-prime5-peaks`: Reference TSS (for ref 5' metrics)
- `--ref-prime3-peaks`: Reference TTS (for ref 3' metrics)

### Parameters
- `--window`: Distance threshold (default 50bp)
- `--stage`: "collapse" or "transcriptome"
- `--verbose`: Enable debug logging

## Output Format

Single-row TSV with columns:
```
stage
isoforms_observed
genes_observed
assigned_unique_read_ids
assigned_primary_alignments
assigned_supplementary_alignments
assigned_total_alignments
input_molecules
input_primary_alignments
input_supplementary_alignments
input_total_alignments
assignment_rate
primary_alignment_utilization
total_alignment_utilization
5prime_precision
5prime_recall
5prime_f1
3prime_precision
3prime_recall
3prime_f1
ref5prime_precision
ref5prime_recall
ref5prime_f1
ref3prime_precision
ref3prime_recall
ref3prime_f1
```

## Dependencies

- Python 3.7+
- pandas
- numpy
- pysam (optional, for supplementary alignment counting)
- bedtools (command-line tool)
- samtools (command-line tool)

## Example Usage

```bash
python3 ted_standalone.py \\
    --isoforms-bed sample.isoforms.bed \\
    --map-file sample.isoform.read.map.txt \\
    --bam sample.bam \\
    --corrected-bed sample_all_corrected.bed \\
    --prime5-peaks cage_peaks.bed \\
    --prime3-peaks quantseq_peaks.bed \\
    --ref-prime5-peaks reference_tss.bed \\
    --ref-prime3-peaks reference_tts.bed \\
    --window 50 \\
    --stage collapse \\
    --output sample.ted.tsv
```

## Advantages of Standalone Version

1. **Explicit dependencies**: No hidden file lookups
2. **Nextflow-friendly**: One input → one output
3. **Testable**: Can run independently
4. **Cacheable**: Nextflow can cache results
5. **Parallelizable**: Nextflow handles parallelization
6. **Debuggable**: Clear input/output contract
7. **Portable**: No pipeline-specific imports

## Migration from Original TED

If you have existing TED configs, map them as follows:

```python
# Old (config-based)
[run]
experiment_5_prime_regions_bed_file = "path/to/cage.bed"

# New (CLI args)
--prime5-peaks path/to/cage.bed
```

## Next Steps

1. Make script executable: `chmod +x bin/ted_standalone.py`
2. Add to Nextflow process (see Option A above)
3. Wire up inputs from FlairTranscriptome/FlairCollapse outputs
4. Add Dataset fields for peak files if not already present
