#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.json.JsonSlurper

// =============================================================================
// UTILITY FUNCTIONS
// =============================================================================

// Sanitize mode names to use hyphens instead of underscores
// This ensures underscores can be used as delimiters in filenames
def sanitizeModeName(String name) {
    return name.replaceAll('_', '-')
}

// Sanitize all keys in a mode map
def sanitizeModeMap(Map modes) {
    if (modes == null || modes.isEmpty()) {
        return modes
    }
    def result = [:]
    modes.each { k, v ->
        result[sanitizeModeName(k)] = v
    }
    return result
}

// =============================================================================
// CLASS DEFINITIONS - Define dataset and test configuration objects
// =============================================================================

class Dataset {
    String name
    def reads  // Can be String (single file) or List (multiple files)
    String genome
    String gtf
    String cage
    String quantseq
    String bam  // Optional: pre-aligned BAM file
    String bai  // Optional: BAM index file (auto-constructed from bam path)
    String junction_tab  // Optional: short-read junction file for flair transcriptome

    // Constructor
    Dataset(String name, Map config) {
        this.name = name
        this.reads = config.reads
        this.genome = config.genome
        this.gtf = config.gtf
        this.cage = config.cage
        this.quantseq = config.quantseq
        this.bam = config.bam
        // Auto-construct BAI path from BAM path if BAM is provided
        this.bai = config.bam ? "${config.bam}.bai" : null
        this.junction_tab = config.junction_tab
    }

    // Helper methods
    boolean hasCage() { return cage != null }
    boolean hasQuantseq() { return quantseq != null }
    boolean hasBam() { return bam != null }
    boolean hasJunctionTab() { return junction_tab != null }

    // Get reads as a list (handles both single file and multiple files)
    List<String> getReadsList() {
        return reads instanceof List ? reads : [reads]
    }

    String description() {
        return "Dataset[${name}]"
    }
}

class TestSet {
    String name
    Dataset dataset
    Map alignModes
    Map partitionModes
    Map transcriptomeModes

    // Constructor with simplified parameter structure
    TestSet(String name, Dataset dataset, Map modes) {
        this.name = name
        this.dataset = dataset
        this.alignModes = modes.align ?: [:]
        this.partitionModes = modes.partition ?: [:]
        this.transcriptomeModes = modes.transcriptome ?: [:]
    }

    // Calculate total number of jobs this test set will produce
    int totalJobs() {
        def align_count = alignModes.size() ?: 1
        def partition_count = partitionModes.size() ?: 1

        return align_count * partition_count
    }

    String description() {
        return "TestSet[${name}, ${totalJobs()} jobs]"
    }
}

// =============================================================================
// PROCESS DEFINITIONS
// =============================================================================

process FlairAlign {
    // Aligns long-read RNA-seq data to the genome using minimap2 (via FLAIR)
    // Produces BAM, BAI, and BED12 files for downstream analysis
    publishDir "results/align", mode: 'symlink'
    publishDir "results/logs/align", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${it}" }
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}"

    input:
    // test_name: Passed through for downstream tracking (not used in script)
    // dataset_name: Used for output file naming
    // reads: FASTQ file(s) with long-read RNA-seq data
    // align_mode: Mode identifier for this alignment configuration
    // align_args: Additional command-line arguments for flair align
    // genome: Reference genome FASTA file
    tuple val(test_name), val(dataset_name), path(reads), val(align_mode), val(align_args), path(genome)

    output:
    // Produces aligned BAM with index and BED12 representation of reads
    tuple val(test_name), val(dataset_name), val(align_mode),
          path("${dataset_name}_${align_mode}.bam"),
          path("${dataset_name}_${align_mode}.bam.bai"),
          path("${dataset_name}_${align_mode}.bed"), emit: alignments

    script:
    """
    # FLAIR align wraps minimap2 for spliced alignment of long reads
    # Automatically generates BAM, BAI, and BED12 files
    flair align ${align_args} \\
        -r ${reads} \\
        -g ${genome} \\
        -o ${dataset_name}_${align_mode}
    """
}

process BamToBed {
    // Converts pre-aligned BAM files to BED12 format for downstream FLAIR processes
    // Used when user provides pre-aligned BAM instead of raw reads
    publishDir "results/align", mode: 'symlink'
    tag "${dataset_name}_${align_mode}"

    input:
    // test_name: Passed through for downstream tracking (not used in script)
    // dataset_name: Used for output file naming
    // align_mode: Mode identifier for this alignment
    // bam: Pre-aligned BAM file
    // bai: BAM index file (required implicitly for efficient BAM access)
    tuple val(test_name), val(dataset_name), val(align_mode), path(bam), path(bai)

    output:
    // Returns original BAM/BAI plus newly generated BED12 file
    tuple val(test_name), val(dataset_name), val(align_mode),
          path(bam), path(bai), path("${dataset_name}_${align_mode}.bed"), emit: alignments

    script:
    """
    # Convert BAM alignments to BED12 format (preserves splice junctions)
    bedtools bamtobed -bed12 -i ${bam} > ${dataset_name}_${align_mode}.bed
    """
}

process FlairPartition {
    // Partitions data to a specific genomic region or subset for testing
    // Extracts alignments, genome sequence, annotations, and experimental peak files for the specified region
    publishDir "results/partition/${test_name}", mode: 'symlink'
    tag "${dataset_name}_${align_mode}_${partition_mode}"

    input:
    // test_name: Used for publishDir organization (not used in script)
    // dataset_name, align_mode, partition_mode: Used for output file naming
    // bam, bai: Input alignment files (bai required implicitly by samtools)
    // bed: BED12 representation of alignments
    // partition_mode: Mode identifier for this partition configuration
    // partition_args: Arguments specifying region (e.g., --region chr19:1000-2000 or --all)
    // genome: Full reference genome FASTA
    // gtf: Full reference annotation GTF
    // cage_peaks: CAGE-seq TSS peaks (optional, may be NO_CAGE placeholder)
    // quantseq_peaks: QuantSeq TTS peaks (optional, may be NO_QUANTSEQ placeholder)
    tuple val(test_name), val(dataset_name), val(align_mode), path(bam), path(bai), path(bed),
          val(partition_mode), val(partition_args), path(genome), path(gtf),
          path(cage_peaks), path(quantseq_peaks)

    output:
    // Produces subsetted BAM, BED, genome, GTF, and peak files for the specified partition
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}.bam"),
          path("${dataset_name}_${align_mode}_${partition_mode}.bam.bai"),
          path("${dataset_name}_${align_mode}_${partition_mode}.bed"),
          path("${dataset_name}_${align_mode}_${partition_mode}_genome.fa"),
          path("${dataset_name}_${align_mode}_${partition_mode}_annotation.gtf"),
          path("${dataset_name}_${align_mode}_${partition_mode}_cage.bed", optional: true),
          path("${dataset_name}_${align_mode}_${partition_mode}_quantseq.bed", optional: true), emit: partitioned

    script:
    def output_prefix = "${dataset_name}_${align_mode}_${partition_mode}"
    def cage_arg = cage_peaks.name != 'NO_CAGE' ? "--cage-peaks ${cage_peaks}" : ""
    def quantseq_arg = quantseq_peaks.name != 'NO_QUANTSEQ' ? "--quantseq-peaks ${quantseq_peaks}" : ""

    """
    # Custom script to partition data to a specific genomic region
    # Subsets BAM, BED, genome sequence, GTF annotation, and experimental peak files to the target region
    python ${projectDir}/bin/simple_partition.py \\
        --bam ${bam} \\
        --bed ${bed} \\
        --genome ${genome} \\
        --gtf ${gtf} \\
        ${cage_arg} \\
        ${quantseq_arg} \\
        --output-prefix ${output_prefix} \\
        ${partition_args}
    """
}

process FlairTranscriptome {
    publishDir "results/transcriptome/${test_name}", mode: 'symlink'
    publishDir "results/logs/${test_name}", mode: 'copy', pattern: '*.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_transcriptome.${it.tokenize('.')[-1]}" }
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome"

    input:
    // test_name: Used for publishDir organization
    // dataset_name, align_mode, partition_mode: Used for output file naming
    // partition_args, transcriptome_mode: Passed through for downstream tracking
    // bam, bai: Input alignment files (bai required implicitly by samtools)
    // genome: Reference genome FASTA for isoform sequence extraction
    // gtf: Reference annotation for guided assembly
    // transcriptome_args: Additional command-line arguments for flair transcriptome
    //                     If args contain '--junction_tab' flag, will use junction_tab file from sample
    // junction_tab: Optional short-read junction file (may be NO_JUNCTION_TAB placeholder)
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), val(partition_args),
          path(bam), path(bai), path(genome), path(gtf),
          val(transcriptome_mode), val(transcriptome_args), path(junction_tab)

    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), val(partition_args), val(transcriptome_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome.isoforms.bed", optional: true),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome.isoforms.gtf", optional: true),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome.isoforms.fa", optional: true),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome.isoform.counts.txt", optional: true),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome.isoform.read.map.txt", optional: true), emit: transcriptome

    script:
    // Check if transcriptome_args requests junction_tab AND sample has a valid junction_tab file
    // The '--junction_tab' in args acts as a flag to enable junction_tab usage
    def use_junction_tab = transcriptome_args.contains('--junction_tab') && junction_tab.name != 'NO_JUNCTION_TAB'
    // Build the junction_tab argument (only if both conditions met)
    def junction_tab_arg = use_junction_tab ? "--junction_tab ${junction_tab}" : ""
    // Remove the '--junction_tab' flag from args (it's just a marker, actual path is added separately)
    def cleaned_args = transcriptome_args.replaceAll('--junction_tab\\s*', '')

    """
    # Run FLAIR transcriptome assembly to generate isoform models
    # Creates BED, GTF, FASTA, counts, and read-to-isoform mapping files
    flair transcriptome \\
        -b ${bam} \\
        --genome ${genome} \\
        -f ${gtf} \\
        ${junction_tab_arg} \\
        ${cleaned_args} \\
        -o ${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome
    """
}

process FlairEvaluation {
    // Evaluates isoform predictions using two complementary metrics:
    // 1. TED (Transcript End Distance) - measures TSS/TTS accuracy
    // 2. FLAIR eval - measures structural accuracy against reference annotation
    publishDir "results/evaluations/individual/${test_name}", mode: 'symlink', pattern: '*.tsv'
    publishDir "results/evaluations/ted_plots", mode: 'copy', pattern: 'ted_plots/*.png'
    publishDir "results/evaluations/ted_test_regions", mode: 'copy', pattern: 'test_regions/*.bed'
    publishDir "results/logs/${test_name}", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_${process_label}_${stage}_${it}" }
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_${process_label}_${stage}"

    input:
    // Metadata for tracking and output naming
    // test_name, dataset_name, align_mode, partition_mode: Identify this specific analysis
    // partition_args: Partition configuration (e.g., genomic region)
    // process_label: Pipeline stage being evaluated (e.g., "transcriptome_with-gtf")
    // stage: Evaluation stage identifier
    //
    // Data files for evaluation
    // isoforms_bed: Predicted isoforms from FLAIR
    // isoform_read_map: Mapping of reads to isoforms
    // bam, bai: Alignment files (used for read coverage analysis)
    // reads_bed: Original aligned reads in BED12 format
    // corrected_bed: Optional corrected reads (may be NO_CORRECTED placeholder)
    // genome, gtf: Reference files for comparison
    //
    // Peak files for TED evaluation (optional, may be NO_* placeholders)
    // cage_peaks: CAGE-seq TSS peaks (5' ends)
    // quantseq_peaks: QuantSeq TTS peaks (3' ends)
    // ref_tss, ref_tts: Reference TSS/TTS from annotation
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(partition_args), val(process_label), val(stage),
          path(isoforms_bed), path(isoform_read_map), path(bam), path(bai),
          path(reads_bed), path(corrected_bed), path(genome), path(gtf),
          path(cage_peaks), path(quantseq_peaks),
          path(ref_tss), path(ref_tts)

    output:
    // Two TSV files with evaluation metrics
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(partition_args), val(process_label), val(stage),
          path("${dataset_name}_${align_mode}_${partition_mode}_${process_label}_${stage}_ted.tsv"),
          path("${dataset_name}_${align_mode}_${partition_mode}_${process_label}_${stage}_flair_eval.tsv"), emit: evaluation_results
    // Distance histogram plots (CAGE, QuantSeq, and Reference TSS/TTS)
    path "ted_plots/*_cage_distance_histogram.png", optional: true, emit: cage_plots
    path "ted_plots/*_quantseq_distance_histogram.png", optional: true, emit: quantseq_plots
    path "ted_plots/*_ref_tss_distance_histogram.png", optional: true, emit: ref_tss_plots
    path "ted_plots/*_ref_tts_distance_histogram.png", optional: true, emit: ref_tts_plots
    path "ted_plots/*_read_tss_offset_histogram.png", optional: true, emit: read_tss_offset_plots
    path "ted_plots/*_read_tts_offset_histogram.png", optional: true, emit: read_tts_offset_plots
    path "ted_plots/*_tss_entropy_distribution.png", optional: true, emit: tss_entropy_plots
    path "ted_plots/*_tts_entropy_distribution.png", optional: true, emit: tts_entropy_plots
    path "ted_plots/*_cage_peak_read_support.png", optional: true, emit: cage_read_support_plots
    path "ted_plots/*_quantseq_peak_read_support.png", optional: true, emit: quantseq_read_support_plots
    // Recoverable peak BED files (peaks with at least one long read end within window)
    path "test_regions/*_recoverable_cage_peaks.bed", optional: true, emit: recoverable_cage_peaks
    path "test_regions/*_recoverable_quantseq_peaks.bed", optional: true, emit: recoverable_quantseq_peaks

    script:
    // Build optional arguments - only include if files are not placeholders and not empty
    // Note: Empty CAGE/QuantSeq files may be created by partition script to satisfy Nextflow output requirements
    // We check file size to avoid passing empty files to evaluation (ted.py handles None gracefully)
    def cage_arg = (cage_peaks.name != 'NO_CAGE' && cage_peaks.size() > 0) ? "--prime5-peaks ${cage_peaks}" : ""
    def quantseq_arg = (quantseq_peaks.name != 'NO_QUANTSEQ' && quantseq_peaks.size() > 0) ? "--prime3-peaks ${quantseq_peaks}" : ""
    def ref_tss_arg = ref_tss.name != 'NO_REF_TSS' ? "--ref-prime5-peaks ${ref_tss}" : ""
    def ref_tts_arg = ref_tts.name != 'NO_REF_TTS' ? "--ref-prime3-peaks ${ref_tts}" : ""
    def corrected_arg = corrected_bed.name != 'NO_CORRECTED' ? "--corrected-bed ${corrected_bed}" : ""

    """
    # Create output directories for TED plots and test regions
    mkdir -p ted_plots
    mkdir -p test_regions

    # TED (Transcript End Distance) evaluation
    # Measures distance between predicted and actual TSS/TTS positions
    # Uses CAGE-seq (5' ends) and QuantSeq (3' ends) data if available
    # Also generates distance histogram plots showing transcript end accuracy
    python ${projectDir}/bin/ted.py \\
        --isoforms-bed ${isoforms_bed} \\
        --map-file ${isoform_read_map} \\
        --bam ${bam} \\
        --reads-bed ${reads_bed} \\
        ${corrected_arg} \\
        ${cage_arg} \\
        ${quantseq_arg} \\
        ${ref_tss_arg} \\
        ${ref_tts_arg} \\
        --window 50 \\
        --stage ${stage} \\
        --test-name ${test_name} \\
        --dataset-name ${dataset_name} \\
        --align-mode ${align_mode} \\
        --partition-mode ${partition_mode} \\
        --pipeline-mode ${process_label} \\
        --plot-output-dir ted_plots \\
        --test-regions-dir test_regions \\
        --output ${dataset_name}_${align_mode}_${partition_mode}_${process_label}_${stage}_ted.tsv \\
        --verbose

    # FLAIR evaluation metrics
    # Compares predicted isoforms against reference annotation
    # Measures sensitivity, precision, and structural accuracy
    python ${projectDir}/bin/flair_eval.py \\
        --reads-bed ${reads_bed} \\
        --isoforms-bed ${isoforms_bed} \\
        --gtf ${gtf} \\
        --test-name ${test_name} \\
        --dataset-name ${dataset_name} \\
        --align-mode ${align_mode} \\
        --partition-mode ${partition_mode} \\
        --pipeline-mode ${process_label} \\
        --stage ${stage} \\
        --output ${dataset_name}_${align_mode}_${partition_mode}_${process_label}_${stage}_flair_eval.tsv \\
        --verbose
    """
}

process PrepareReferencePeaks {
    // Extracts TSS (transcription start sites) and TTS (transcription termination sites)
    // from reference annotation GTF for use as ground truth in TED evaluation
    publishDir "results/reference_peaks/${test_name}", mode: 'symlink'
    tag "${dataset_name}_${align_mode}_${partition_mode}"

    input:
    // test_name, dataset_name, align_mode, partition_mode: Used for tracking and naming
    // gtf: Reference annotation to extract TSS/TTS positions from
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), path(gtf)

    output:
    // Two BED files: one with TSS positions, one with TTS positions
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_ref_tss.bed"),
          path("${dataset_name}_${align_mode}_${partition_mode}_ref_tts.bed"), emit: ref_peaks

    script:
    """
    # Extract TSS and TTS positions from GTF annotation
    # --deduplicate removes redundant positions from overlapping transcripts
    python ${projectDir}/bin/gtf_to_tss_tts.py \\
        --gtf ${gtf} \\
        --output-prefix ${dataset_name}_${align_mode}_${partition_mode}_ref \\
        --deduplicate
    """
}

process SynthesizeEvaluations {
    // Aggregates all individual evaluation results into a single comprehensive summary
    // Combines TED and FLAIR evaluation metrics across all test configurations
    publishDir "results/evaluations/summary", mode: 'symlink'
    tag "${test_name}_evaluation_summary"

    input:
    // test_name: Overall test suite identifier
    // ted_files: All TED evaluation TSV files from FlairEvaluation
    // eval_files: All FLAIR evaluation TSV files from FlairEvaluation
    tuple val(test_name), path(ted_files), path(eval_files)

    output:
    // Single comprehensive TSV with all evaluation metrics
    path "${test_name}_evaluation_summary.tsv", emit: evaluation_summary

    script:
    """
    # Synthesize all evaluation results into one summary table
    # Merges metrics from multiple runs for easy comparison
    python ${projectDir}/bin/synthesize_evaluations.py \\
        --ted-files ${ted_files.join(' ')} \\
        --flair-files ${eval_files.join(' ')} \\
        --output ${test_name}_evaluation_summary.tsv \\
        --test-name ${test_name}
    """
}

process PlotIsoforms {
    // Generates visualization plots of isoform structures and read assignments
    // Only runs for partitions ≤ 400kb to avoid generating overly complex plots
    publishDir "results/isoform_plots/${test_name}", mode: 'copy'
    publishDir "results/logs/${test_name}", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_plot_${it}" }
    conda '/private/home/hdheath/miniforge3/envs/nextflow_env'  // TODO: Make this configurable
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}"

    input:
    // test_name: Used for publishDir organization (not used in script)
    // dataset_name, align_mode, partition_mode: Used for output file naming
    // bam, bai: Alignment files for read coverage visualization
    // isoforms_bed: Predicted isoforms to plot
    // isoform_read_map: Read-to-isoform assignments for coloring/grouping
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path(bam), path(bai), path(isoforms_bed), path(isoform_read_map)

    output:
    // PNG image visualizing isoform structures and read support
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_isoform_plot.png"), emit: plots

    script:
    """
    # Generate isoform structure plot with read coverage
    # Shows splice junctions, exons, and read assignments
    python ${projectDir}/bin/isoform_plot_test.py \\
        --bam ${bam} \\
        --readmap ${isoform_read_map} \\
        --isoforms ${isoforms_bed} \\
        --output ${dataset_name}_${align_mode}_${partition_mode}_isoform_plot
    """
}

workflow {
    // Validate input parameters
    if (!params.input) {
        error "ERROR: Please provide a samplesheet via"
    }
    if (!params.params_file) {
        error "ERROR: Please provide a parameters configuration file via --params_file"
    }

    // Display test name warning if using default
    if (params.test_name == 'flair_test_suite') {
        log.warn "WARNING: Using default test name 'flair_test_suite'. Specify a custom name with --test_name"
    }

    // Read and parse pipeline parameters configuration from JSON file
    def jsonSlurper = new JsonSlurper()
    def pipeline_config = jsonSlurper.parse(file(params.params_file))

    // Extract mode maps from JSON config (with defaults if not specified)
    def align_modes = pipeline_config.align ?: [default: '']
    def partition_modes = pipeline_config.partition ?: [all: '--all']
    def transcriptome_modes = pipeline_config.transcriptome ?: [default: '']

    // Read and parse samplesheet CSV into a list
    def test_sets_list = []
    file(params.input).withReader { reader ->
        def header = reader.readLine().split(',')
        reader.eachLine { line ->
            // Skip empty lines
            if (!line.trim()) return

            // Split and trim values, handling trailing commas
            def values = line.split(',', -1)  // -1 keeps trailing empty strings

            // Ensure values array matches header length
            if (values.size() > header.size()) {
                values = values[0..<header.size()]  // Truncate extra values
            } else if (values.size() < header.size()) {
                // Pad with empty strings if needed
                values = values + [''] * (header.size() - values.size())
            }

            def row = [header, values].transpose().collectEntries()

            // Create Dataset object from CSV row
            def dataset = new Dataset(row.sample_id, ([
                genome: file(row.genome),
                gtf: file(row.gtf),
                bam: row.bam && row.bam != '' ? file(row.bam) : null,
                reads: row.reads && row.reads != '' ? file(row.reads) : null,
                cage: row.cage && row.cage != '' ? file(row.cage) : null,
                quantseq: row.quantseq && row.quantseq != '' ? file(row.quantseq) : null,
                junction_tab: row.junction_tab && row.junction_tab != '' ? file(row.junction_tab) : null
            ]))

            // Create TestSet with the specified test_name and modes from JSON config
            test_sets_list.add(new TestSet("${params.test_name}_${dataset.name}", dataset, ([
                align: sanitizeModeMap(align_modes),
                partition: sanitizeModeMap(partition_modes),
                transcriptome: sanitizeModeMap(transcriptome_modes)
            ])))
        }
    }

    // Use test_name for the comprehensive summary
    def test_set_name = params.test_name

    // Print summary
    println("\n=== FLAIR Test Suite ===")
    println("Input samplesheet: ${params.input}")
    println("Test name: ${test_set_name}")
    println("===================================")
    test_sets_list.each { test ->
        println("${test.name}: ${test.totalJobs()} jobs using dataset '${test.dataset.name}'")
        if (test.dataset.hasBam()) {
            println("  ⚠ WARNING: Pre-aligned BAM detected for '${test.dataset.name}' - FlairAlign will be SKIPPED")
            if (test.dataset.reads != null) {
                println("  ⚠ Note: 'reads' field will be ignored in favor of provided BAM files")
            }
        }
    }
    println("Total jobs: ${test_sets_list.sum { it.totalJobs() }}")
    println("===================================\n")

    // =============================================================================
    // CHANNEL CONSTRUCTION - Build input channels for processes
    // =============================================================================

    // Create main datasets channel from test sets list
    // Each element: [test_name, dataset, align_modes, partition_modes, transcriptome_modes]
    datasets_ch = Channel.from(test_sets_list)
        .map { test_set ->
            [test_set.name, test_set.dataset, test_set.alignModes, test_set.partitionModes,
             test_set.transcriptomeModes]
        }

    // Split datasets into two branches based on whether BAM files are provided
    // Branch 1: Pre-aligned BAM files provided (skip FlairAlign)
    datasets_with_bam = datasets_ch.filter { test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes ->
        dataset.hasBam()
    }

    // Branch 2: Raw reads provided (run FlairAlign)
    datasets_without_bam = datasets_ch.filter { test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes ->
        !dataset.hasBam()
    }

    // Prepare inputs for FlairAlign process
    // Expands datasets into individual alignment jobs (one per read file per align mode)
    align_inputs = datasets_without_bam.flatMap { test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes ->
        // For each alignment mode configuration
        ds_align_modes.collectMany { align_mode, align_args ->
            // For each reads file (handles both single file and multi-file datasets)
            dataset.getReadsList().collect { reads_file ->
                [test_name, dataset.name, file(reads_file), align_mode, align_args, file(dataset.genome)]
            }
        }
    }

    // Prepare inputs for BamToBed process (pre-aligned data)
    // Creates one job per align mode for each pre-aligned BAM file
    prealigned_bam_inputs = datasets_with_bam.flatMap { test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes ->
        ds_align_modes.collect { align_mode, align_args ->
            [test_name, dataset.name, align_mode, file(dataset.bam), file(dataset.bai)]
        }
    }

    // =============================================================================
    // ALIGNMENT STAGE - Run alignment or convert pre-aligned BAMs
    // =============================================================================

    FlairAlign(align_inputs)
    BamToBed(prealigned_bam_inputs)

    // Merge alignment outputs from both branches (FlairAlign and BamToBed)
    // Both produce: [test_name, dataset_name, align_mode, bam, bai, bed]
    // Use concat instead of mix for better empty channel handling
    all_alignments = FlairAlign.out.alignments.concat(BamToBed.out.alignments)

    // =============================================================================
    // PARTITION STAGE - Subset data to specific genomic regions
    // =============================================================================

    // Prepare inputs for FlairPartition by combining alignments with partition modes
    // For each alignment, create one partition job per partition mode
    partition_inputs = all_alignments.combine(datasets_ch, by: 0).flatMap {
        test_name, dataset_name, align_mode, bam, bai, bed, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes ->
        // Prepare CAGE and QuantSeq files (use placeholders if not provided)
        def cage_file = dataset.cage ? file(dataset.cage) : file("${workflow.workDir}/NO_CAGE", checkIfExists: false)
        def quantseq_file = dataset.quantseq ? file(dataset.quantseq) : file("${workflow.workDir}/NO_QUANTSEQ", checkIfExists: false)
        // Generate tuple for each partition mode configuration
        ds_partition_modes.collect { partition_mode, partition_args ->
            [test_name, dataset_name, align_mode, bam, bai, bed, partition_mode, partition_args,
             file(dataset.genome), file(dataset.gtf), cage_file, quantseq_file]
        }
    }

    FlairPartition(partition_inputs)

    // =============================================================================
    // TRANSCRIPTOME STAGE - Generate isoform predictions
    // =============================================================================

    // Prepare inputs for FlairTranscriptome by combining partitioned data with transcriptome modes
    // For each partitioned output, generate one transcriptome job per transcriptome mode
    transcriptome_inputs = FlairPartition.out.partitioned.combine(datasets_ch, by: 0).flatMap {
        test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, quantseq_peaks,
        dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes ->
        // Retrieve partition arguments for this specific partition mode
        def partition_args = ds_partition_modes[partition_mode] ?: ''
        // Prepare junction_tab file (use placeholder if not provided)
        def junction_tab_file = dataset.junction_tab ? file(dataset.junction_tab) : file("${workflow.workDir}/NO_JUNCTION_TAB", checkIfExists: false)
        // Generate a tuple for each transcriptome mode configuration
        ds_transcriptome_modes.collect { transcriptome_mode, transcriptome_args ->
            // Note: bed, cage_peaks, quantseq_peaks are not used by flair transcriptome
            [test_name, dataset_name, align_mode, partition_mode, partition_args, bam, bai, genome, gtf,
             transcriptome_mode, transcriptome_args, junction_tab_file]
        }
    }

    FlairTranscriptome(transcriptome_inputs)

    // =============================================================================
    // EVALUATION SETUP - Prepare reference data and organize inputs
    // =============================================================================

    // Extract GTF files from partitioned data to generate reference TSS/TTS peaks
    ref_peak_inputs = FlairPartition.out.partitioned.map {
        test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, quantseq_peaks ->
        [test_name, dataset_name, align_mode, partition_mode, gtf]
    }

    PrepareReferencePeaks(ref_peak_inputs)

    // Prepare alignment data for evaluation (keep only necessary fields)
    align_for_eval = all_alignments.map { test_name, dataset_name, align_mode, bam, bai, bed ->
        [test_name, dataset_name, align_mode, bam, bai, bed]
    }

    // =============================================================================
    // EVALUATION STAGE - Combine all data sources for comprehensive evaluation
    // =============================================================================
    //
    // This is the most complex channel operation in the pipeline. It combines:
    // 1. Transcriptome outputs (isoform predictions)
    // 2. Partitioned data (BAM, BED, genome, GTF, CAGE, QuantSeq)
    // 3. Reference TSS/TTS peaks
    //
    // The combine operations use 'by' parameter to match on common keys:
    // - by: [0, 1, 2, 3] matches on [test_name, dataset_name, align_mode, partition_mode]

    eval_transcriptome_inputs = FlairTranscriptome.out.transcriptome
        // Step 1: Extract isoform data and add process labels
        .map { test_name, dataset_name, align_mode, partition_mode, partition_args, transcriptome_mode,
               isoforms_bed, isoforms_gtf, isoforms_fa, isoform_counts, isoform_read_map ->
            // Add descriptive process_label and stage identifier
            // Use placeholder for corrected_bed (not used in transcriptome-only evaluation)
            [test_name, dataset_name, align_mode, partition_mode, partition_args, "transcriptome_${transcriptome_mode}",
             "transcriptome", isoforms_bed, isoform_read_map, file("${workflow.workDir}/NO_CORRECTED", checkIfExists: false)]
        }
        // Step 2: Add partitioned BAM/BED/genome/GTF/CAGE/QuantSeq (match on test_name, dataset_name, align_mode, partition_mode)
        .combine(FlairPartition.out.partitioned.map { test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, quantseq_peaks ->
            [test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, quantseq_peaks]
        }, by: [0, 1, 2, 3])
        // Step 3: Add reference TSS/TTS peaks (match on test_name, dataset_name, align_mode, partition_mode)
        .combine(PrepareReferencePeaks.out.ref_peaks, by: [0, 1, 2, 3])
        // Step 4: Reorganize into final evaluation input tuple
        .map { test_name, dataset_name, align_mode, partition_mode, partition_args, process_label, stage,
               isoforms_bed, isoform_read_map, corrected_file, bam, bai, reads_bed, genome, gtf, cage_peaks, quantseq_peaks,
               ref_tss, ref_tts ->
            // Final tuple with all inputs needed for FlairEvaluation
            // Use partitioned CAGE/QuantSeq files or placeholders
            def cage_file = cage_peaks ?: file("${workflow.workDir}/NO_CAGE", checkIfExists: false)
            def quantseq_file = quantseq_peaks ?: file("${workflow.workDir}/NO_QUANTSEQ", checkIfExists: false)
            [test_name, dataset_name, align_mode, partition_mode, partition_args, process_label, stage,
             isoforms_bed, isoform_read_map, bam, bai, reads_bed, corrected_file, genome, gtf,
             cage_file, quantseq_file, ref_tss, ref_tts]
        }

    FlairEvaluation(eval_transcriptome_inputs)

    // =============================================================================
    // VISUALIZATION - Generate isoform plots for small regions
    // =============================================================================

    // Prepare inputs for plotting, but only for small regions (≤ 400kb)
    // Larger regions would produce overly complex plots
    plot_inputs = FlairTranscriptome.out.transcriptome
        .map { test_name, dataset_name, align_mode, partition_mode, partition_args, transcriptome_mode,
               isoforms_bed, isoforms_gtf, isoforms_fa, isoform_counts, isoform_read_map ->
            // Extract region size from partition_args if present
            // Looks for patterns like: --region chr19:10000-20000
            def region_size = 0
            if (partition_args && partition_args.contains('--region')) {
                def matcher = partition_args =~ /chr\w+:(\d+)-(\d+)/
                if (matcher.find()) {
                    def start = matcher.group(1).toLong()
                    def end = matcher.group(2).toLong()
                    region_size = end - start
                }
            }
            [test_name, dataset_name, align_mode, partition_mode, partition_args, region_size,
             isoforms_bed, isoform_read_map]
        }
        // Filter: only include regions ≤ 400kb for plotting
        .filter { test_name, dataset_name, align_mode, partition_mode, partition_args, region_size,
                  isoforms_bed, isoform_read_map ->
            region_size > 0 && region_size <= 400000
        }
        // Add BAM files from partition output (match on test_name, dataset_name, align_mode, partition_mode)
        .combine(FlairPartition.out.partitioned.map { test_name, dataset_name, align_mode, partition_mode,
                                                        bam, bai, bed, genome, gtf, cage_peaks, quantseq_peaks ->
            [test_name, dataset_name, align_mode, partition_mode, bam, bai]
        }, by: [0, 1, 2, 3])
        // Remove region_size (no longer needed) and reorganize for PlotIsoforms input
        .map { test_name, dataset_name, align_mode, partition_mode, partition_args, region_size,
               isoforms_bed, isoform_read_map, bam, bai ->
            [test_name, dataset_name, align_mode, partition_mode, bam, bai, isoforms_bed, isoform_read_map]
        }

    PlotIsoforms(plot_inputs)

    // =============================================================================
    // SUMMARY GENERATION - Aggregate all evaluation results
    // =============================================================================

    // Collect all TED evaluation files from all runs
    // .collect() waits for all FlairEvaluation processes to complete and gathers files into a list
    all_ted_files = FlairEvaluation.out.evaluation_results
        .map { test_name, dataset_name, align_mode, partition_mode, partition_args, process_label, stage, ted_file, flair_file ->
            ted_file
        }
        .collect()

    // Collect all FLAIR evaluation files from all runs
    all_flair_files = FlairEvaluation.out.evaluation_results
        .map { test_name, dataset_name, align_mode, partition_mode, partition_args, process_label, stage, ted_file, flair_file ->
            flair_file
        }
        .collect()

    // Combine TED and FLAIR file lists into a single channel for synthesis
    // Use the overall test_set_name (not individual test names) for the comprehensive summary
    all_results_input = all_ted_files
        .map { ted_files -> [test_set_name, ted_files] }
        .join(all_flair_files.map { flair_files -> [test_set_name, flair_files] })

    // Generate one comprehensive summary TSV with all evaluation metrics
    // This produces a single file combining results from all configurations
    SynthesizeEvaluations(all_results_input)
}
