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
    Map bambuModes      // Optional: Bambu transcriptome assembler modes
    Map isoquantModes   // Optional: IsoQuant transcriptome assembler modes

    // Constructor with simplified parameter structure
    TestSet(String name, Dataset dataset, Map modes) {
        this.name = name
        this.dataset = dataset
        this.alignModes = modes.align ?: [:]
        this.partitionModes = modes.partition ?: [:]
        this.transcriptomeModes = modes.transcriptome ?: [:]
        this.bambuModes = modes.bambu ?: [:]
        this.isoquantModes = modes.isoquant ?: [:]
    }

    // Helper methods for checking if alternative assemblers are enabled
    boolean hasBambu() { return bambuModes != null && !bambuModes.isEmpty() }
    boolean hasIsoquant() { return isoquantModes != null && !isoquantModes.isEmpty() }

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

process FlairPartition {
    // Partitions data to a specific genomic region or subset for testing
    // Extracts alignments, genome sequence, annotations, and experimental peak files for the specified region
    // Supports two modes:
    //   1. With pre-computed BED file (from FlairAlign)
    //   2. With NO_BED placeholder (from pre-aligned BAM) - generates BED from partitioned BAM
    publishDir "results/partition/${test_name}", mode: 'symlink'
    tag "${dataset_name}_${align_mode}_${partition_mode}"

    input:
    // test_name: Used for publishDir organization (not used in script)
    // dataset_name, align_mode, partition_mode: Used for output file naming
    // bam, bai: Input alignment files (bai required implicitly by samtools)
    // bed: BED12 representation of alignments (or NO_BED placeholder for pre-aligned BAM)
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
    // If BED is a placeholder (NO_BED), use --generate-bed to create BED from partitioned BAM
    // This avoids expensive full BAM->BED conversion for pre-aligned data
    def bed_arg = bed.name != 'NO_BED' ? "--bed ${bed}" : "--generate-bed"

    """
    # Custom script to partition data to a specific genomic region
    # Subsets BAM, BED, genome sequence, GTF annotation, and experimental peak files to the target region
    # When --generate-bed is used, BED is created from the partitioned BAM (more efficient for large files)
    python ${projectDir}/bin/simple_partition.py \\
        --bam ${bam} \\
        ${bed_arg} \\
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

// =============================================================================
// ALTERNATIVE TRANSCRIPTOME ASSEMBLERS
// =============================================================================

process BambuAssembly {
    // Bambu: R-based transcript discovery and quantification
    // Only runs if bambu modes are defined in params.json
    //
    // IMPORTANT: XGBoost version compatibility issue
    // Bambu's pre-trained XGBoost model (for novel transcript discovery) is incompatible
    // with XGBoost >= 2.1.0. If you see "xgb.Booster object is corrupted" error:
    //   Option 1: Downgrade XGBoost to 1.7.x (r-xgboost=1.7.6 in conda)
    //   Option 2: Use discovery=FALSE in bambu_args to disable novel discovery
    //   Option 3: Update to latest Bambu (3.5.1+) with R >= 4.4
    //
    // Key parameters (set via bambu_args in JSON config):
    //   NDR=<0-1>: Novel Discovery Rate threshold (higher = more novel transcripts)
    //   discovery=TRUE/FALSE: Enable/disable novel transcript discovery
    //   quant=TRUE/FALSE: Enable/disable quantification
    publishDir "results/bambu/${test_name}", mode: 'symlink'
    publishDir "results/logs/${test_name}", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_bambu_${bambu_mode}_${it}" }
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_bambu_${bambu_mode}"

    input:
    // test_name, dataset_name, align_mode, partition_mode: Tracking and naming
    // bam, bai: Input alignment files
    // genome: Reference genome FASTA
    // gtf: Reference annotation GTF
    // bambu_mode: Mode identifier for this bambu configuration
    // bambu_args: Additional arguments for bambu (R format, e.g., "NDR=0.1, discovery=FALSE")
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path(bam), path(bai), path(genome), path(gtf),
          val(bambu_mode), val(bambu_args)

    output:
    // GTF with transcript models - primary output for evaluation
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(bambu_mode), path("${dataset_name}_${align_mode}_${partition_mode}_bambu_${bambu_mode}.gtf"),
          emit: bambu_gtf
    // Optional: counts and other outputs
    path "${dataset_name}_${align_mode}_${partition_mode}_bambu_${bambu_mode}_counts_transcript.txt", optional: true, emit: transcript_counts
    path "${dataset_name}_${align_mode}_${partition_mode}_bambu_${bambu_mode}_counts_gene.txt", optional: true, emit: gene_counts

    script:
    def output_prefix = "${dataset_name}_${align_mode}_${partition_mode}_bambu_${bambu_mode}"
    // Parse bambu_args - if empty, use defaults
    def bambu_params = bambu_args ?: ""
    """
    #!/usr/bin/env Rscript

    library(bambu)

    # Prepare annotations (recommended for multiple runs)
    annotations <- prepareAnnotations("${gtf}")

    # Run bambu with specified parameters
    # Note: If XGBoost version incompatibility error occurs, add discovery=FALSE
    # to bambu_args in your params JSON file, or downgrade XGBoost to 1.7.x
    se <- bambu(
        reads = "${bam}",
        annotations = annotations,
        genome = "${genome}",
        ncore = ${task.cpus}${bambu_params ? ", ${bambu_params}" : ""}
    )

    # Write GTF output
    # writeBambuOutput creates: {prefix}extended_annotations.gtf (note: no underscore before "extended")
    writeBambuOutput(se, path = ".", prefix = "${output_prefix}")

    # Rename the extended_annotations output file to our expected naming convention
    file.rename("${output_prefix}extended_annotations.gtf", "${output_prefix}.gtf")
    """
}

process IsoQuantAssembly {
    // IsoQuant: Python-based transcript reconstruction from long reads
    // Only runs if isoquant modes are defined in params.json
    // Key parameters (set via isoquant_args in JSON config):
    //   --data_type: nanopore, pacbio_ccs, or assembly (REQUIRED in args)
    //   --complete_genedb: Use for complete gene databases like GENCODE
    //   --model_construction_strategy: sensitive, default, or fl_pacbio
    publishDir "results/isoquant/${test_name}", mode: 'symlink'
    publishDir "results/logs/${test_name}", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}_${it}" }
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}"

    input:
    // test_name, dataset_name, align_mode, partition_mode: Tracking and naming
    // bam, bai: Input alignment files
    // genome: Reference genome FASTA
    // gtf: Reference annotation GTF
    // isoquant_mode: Mode identifier for this isoquant configuration
    // isoquant_args: Additional arguments for isoquant (e.g., "--data_type nanopore --complete_genedb")
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path(bam), path(bai), path(genome), path(gtf),
          val(isoquant_mode), val(isoquant_args)

    output:
    // GTF with transcript models - primary output for evaluation
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(isoquant_mode), path("${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}.gtf"),
          emit: isoquant_gtf
    // Optional: read assignments and counts
    path "${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}_read_assignments.tsv", optional: true, emit: read_assignments
    path "${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}_transcript_counts.tsv", optional: true, emit: transcript_counts

    script:
    def output_prefix = "${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}"
    """
    # Run IsoQuant
    isoquant.py \\
        --reference ${genome} \\
        --genedb ${gtf} \\
        --bam ${bam} \\
        --output isoquant_out \\
        --threads ${task.cpus} \\
        --prefix ${output_prefix} \\
        ${isoquant_args}

    # Copy outputs to working directory with standardized names
    cp isoquant_out/${output_prefix}/${output_prefix}.transcript_models.gtf ${output_prefix}.gtf

    # Copy optional files if they exist
    if [ -f "isoquant_out/${output_prefix}/${output_prefix}.read_assignments.tsv" ]; then
        cp isoquant_out/${output_prefix}/${output_prefix}.read_assignments.tsv ${output_prefix}_read_assignments.tsv
    fi
    if [ -f "isoquant_out/${output_prefix}/${output_prefix}.transcript_counts.tsv" ]; then
        cp isoquant_out/${output_prefix}/${output_prefix}.transcript_counts.tsv ${output_prefix}_transcript_counts.tsv
    fi
    """
}

process Evaluation {
    // Unified evaluation process for all assemblers (FLAIR, Bambu, IsoQuant)
    // Automatically detects assembler type and runs appropriate evaluation:
    //   - FLAIR: Full evaluation with read-level metrics (entropy, motifs, truncation patterns)
    //   - Bambu/IsoQuant: Simplified evaluation (TSS/TTS precision/recall, structural metrics only)
    //
    // Outputs:
    //   1. TED (Transcript End Distance) - measures TSS/TTS accuracy
    //   2. FLAIR eval - measures structural accuracy against reference annotation
    //   3. Synthesized TSV combining all metrics
    publishDir "results/evaluations/${test_name}", mode: 'symlink', pattern: '*_evaluation.tsv'
    publishDir "results/evaluations/${test_name}/ted_plots", mode: 'copy', pattern: 'ted_plots/*.png'
    publishDir "results/evaluations/${test_name}/test_regions", mode: 'copy', pattern: 'test_regions/*.bed'
    publishDir "results/evaluations/${test_name}/test_regions", mode: 'copy', pattern: 'test_regions/*.csv'
    publishDir "results/evaluations/${test_name}/test_regions", mode: 'copy', pattern: 'test_regions/*.tsv'
    publishDir "results/evaluations/${test_name}/timing", mode: 'copy', pattern: '*_ted_timing.txt'
    publishDir "results/logs/${test_name}", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_${it}" }
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}"

    input:
    // Metadata for tracking and output naming
    // test_name, dataset_name, align_mode, partition_mode: Identify this specific analysis
    // transcriptome_mode: Assembler and mode (e.g., "default", "bambu_default", "isoquant_pacbio")
    //
    // Data files for evaluation
    // isoforms_bed: Predicted isoforms in BED12 format (from FLAIR, or NO_ISOFORMS_BED placeholder for Bambu/IsoQuant)
    // isoforms_gtf: Predicted isoforms in GTF format (from Bambu/IsoQuant, or NO_ISOFORMS_GTF placeholder for FLAIR)
    // isoform_read_map: Mapping of reads to isoforms (from FLAIR, or NO_READ_MAP placeholder for Bambu/IsoQuant)
    // bam, bai: Alignment files (used for read coverage analysis in FLAIR mode)
    // reads_bed: Original aligned reads in BED12 format
    // genome, gtf: Reference files for comparison
    //
    // Peak files for TED evaluation (optional, may be NO_* placeholders)
    // cage_peaks: CAGE-seq TSS peaks (5' ends)
    // quantseq_peaks: QuantSeq TTS peaks (3' ends)
    // ref_tss, ref_tts: Reference TSS/TTS from annotation
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(transcriptome_mode),
          path(isoforms_bed), path(isoforms_gtf), path(isoform_read_map),
          path(bam), path(bai), path(reads_bed), path(genome), path(gtf),
          path(cage_peaks), path(quantseq_peaks),
          path(ref_tss), path(ref_tts)

    output:
    // Single synthesized TSV file with all evaluation metrics
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(transcriptome_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_evaluation.tsv"), emit: evaluation_results
    // Distance histogram plots (CAGE, QuantSeq, and Reference TSS/TTS)
    path "ted_plots/*_cage_distance_histogram.png", optional: true, emit: cage_plots
    path "ted_plots/*_quantseq_distance_histogram.png", optional: true, emit: quantseq_plots
    path "ted_plots/*_ref_tss_distance_histogram.png", optional: true, emit: ref_tss_plots
    path "ted_plots/*_ref_tts_distance_histogram.png", optional: true, emit: ref_tts_plots
    // FLAIR-only outputs (read-level analysis)
    path "ted_plots/*_read_tss_offset_histogram.png", optional: true, emit: read_tss_offset_plots
    path "ted_plots/*_read_tts_offset_histogram.png", optional: true, emit: read_tts_offset_plots
    path "ted_plots/*_tss_entropy_distribution.png", optional: true, emit: tss_entropy_plots
    path "ted_plots/*_tts_entropy_distribution.png", optional: true, emit: tts_entropy_plots
    path "ted_plots/*_cage_peak_read_support.png", optional: true, emit: cage_read_support_plots
    path "ted_plots/*_quantseq_peak_read_support.png", optional: true, emit: quantseq_read_support_plots
    path "ted_plots/*_cage_read_classification.png", optional: true, emit: cage_read_classification_plots
    path "ted_plots/*_quantseq_read_classification.png", optional: true, emit: quantseq_read_classification_plots
    path "ted_plots/*_cage_truncation_patterns.png", optional: true, emit: cage_truncation_plots
    path "ted_plots/*_cage_missed_sj_support.png", optional: true, emit: cage_missed_sj_support_plots
    path "ted_plots/*_quantseq_missed_sj_support.png", optional: true, emit: quantseq_missed_sj_support_plots
    // Peak recovery by expression/width
    path "ted_plots/*_cage_recovery_by_expression.png", optional: true, emit: cage_recovery_expression_plots
    path "ted_plots/*_cage_recovery_by_width.png", optional: true, emit: cage_recovery_width_plots
    path "ted_plots/*_quantseq_recovery_by_expression.png", optional: true, emit: quantseq_recovery_expression_plots
    path "ted_plots/*_tss_motif_logo.png", optional: true, emit: tss_motif_plots
    path "ted_plots/*_tts_motif_logo.png", optional: true, emit: tts_motif_plots
    // Structural evaluation plots (all assemblers)
    path "ted_plots/*_transcript_classification.png", optional: true, emit: transcript_classification_plots
    path "ted_plots/*_splice_junction_support.png", optional: true, emit: splice_junction_support_plots
    // Test region outputs (FLAIR-only)
    path "test_regions/*_recoverable_cage_peaks.bed", optional: true, emit: recoverable_cage_peaks
    path "test_regions/*_recoverable_quantseq_peaks.bed", optional: true, emit: recoverable_quantseq_peaks
    path "test_regions/*_missed_cage_peaks_annotated.csv", optional: true, emit: missed_cage_peaks_annotated
    path "test_regions/*_missed_quantseq_peaks_annotated.csv", optional: true, emit: missed_quantseq_peaks_annotated
    path "test_regions/*_missed_cage_peaks.tsv", optional: true, emit: missed_cage_peaks_tsv
    path "test_regions/*_missed_quantseq_peaks.tsv", optional: true, emit: missed_quantseq_peaks_tsv
    // Performance timing report
    path "*_ted_timing.txt", optional: true, emit: timing_reports

    script:
    // Determine assembler type from transcriptome_mode
    def is_bambu = transcriptome_mode.contains('bambu')
    def is_isoquant = transcriptome_mode.contains('isoquant')
    def is_simplified = is_bambu || is_isoquant

    // Build optional arguments
    def cage_arg = (cage_peaks.name != 'NO_CAGE' && cage_peaks.size() > 0) ? "--prime5-peaks ${cage_peaks}" : ""
    def quantseq_arg = (quantseq_peaks.name != 'NO_QUANTSEQ' && quantseq_peaks.size() > 0) ? "--prime3-peaks ${quantseq_peaks}" : ""
    def ref_tss_arg = ref_tss.name != 'NO_REF_TSS' ? "--ref-prime5-peaks ${ref_tss}" : ""
    def ref_tts_arg = ref_tts.name != 'NO_REF_TTS' ? "--ref-prime3-peaks ${ref_tts}" : ""
    def output_prefix = "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}"

    if (is_simplified)
    """
    # Simplified evaluation for Bambu/IsoQuant (no read-level metrics)
    mkdir -p ted_plots

    # TED evaluation with --gtf-input and --skip-read-metrics
    python ${projectDir}/bin/ted.py \\
        --gtf-input ${isoforms_gtf} \\
        --skip-read-metrics \\
        ${cage_arg} \\
        ${quantseq_arg} \\
        ${ref_tss_arg} \\
        ${ref_tts_arg} \\
        --window 50 \\
        --stage transcriptome \\
        --test-name ${test_name} \\
        --dataset-name ${dataset_name} \\
        --align-mode ${align_mode} \\
        --partition-mode ${partition_mode} \\
        --pipeline-mode ${transcriptome_mode} \\
        --plot-output-dir ted_plots \\
        --output ${output_prefix}_ted.tsv \\
        --verbose

    # FLAIR structural evaluation (using GTF input)
    python ${projectDir}/bin/flair_eval.py \\
        --reads-bed ${reads_bed} \\
        --gtf-input ${isoforms_gtf} \\
        --gtf ${gtf} \\
        --test-name ${test_name} \\
        --dataset-name ${dataset_name} \\
        --align-mode ${align_mode} \\
        --partition-mode ${partition_mode} \\
        --pipeline-mode ${transcriptome_mode} \\
        --stage transcriptome \\
        --plot-output-dir ted_plots \\
        --plot-prefix ${output_prefix} \\
        --output ${output_prefix}_flair_eval.tsv \\
        --verbose

    # Synthesize results
    python ${projectDir}/bin/synthesize_evaluations.py \\
        --ted-files ${output_prefix}_ted.tsv \\
        --flair-files ${output_prefix}_flair_eval.tsv \\
        --output ${output_prefix}_evaluation.tsv \\
        --test-name ${test_name}

    # Clean up intermediate files
    rm -f ${output_prefix}_ted.tsv ${output_prefix}_flair_eval.tsv
    """

    else
    """
    # Full FLAIR evaluation with read-level metrics
    mkdir -p ted_plots
    mkdir -p test_regions

    # TED evaluation with full read-level analysis
    python ${projectDir}/bin/ted.py \\
        --isoforms-bed ${isoforms_bed} \\
        --map-file ${isoform_read_map} \\
        --bam ${bam} \\
        --reads-bed ${reads_bed} \\
        --genome ${genome} \\
        ${cage_arg} \\
        ${quantseq_arg} \\
        ${ref_tss_arg} \\
        ${ref_tts_arg} \\
        --window 50 \\
        --stage transcriptome \\
        --test-name ${test_name} \\
        --dataset-name ${dataset_name} \\
        --align-mode ${align_mode} \\
        --partition-mode ${partition_mode} \\
        --pipeline-mode ${transcriptome_mode} \\
        --plot-output-dir ted_plots \\
        --test-regions-dir test_regions \\
        --timing-output ${output_prefix}_ted_timing.txt \\
        --output ${output_prefix}_ted.tsv \\
        --verbose

    # FLAIR structural evaluation
    python ${projectDir}/bin/flair_eval.py \\
        --reads-bed ${reads_bed} \\
        --isoforms-bed ${isoforms_bed} \\
        --gtf ${gtf} \\
        --test-name ${test_name} \\
        --dataset-name ${dataset_name} \\
        --align-mode ${align_mode} \\
        --partition-mode ${partition_mode} \\
        --pipeline-mode ${transcriptome_mode} \\
        --stage transcriptome \\
        --plot-output-dir ted_plots \\
        --plot-prefix ${output_prefix} \\
        --output ${output_prefix}_flair_eval.tsv \\
        --verbose

    # Synthesize results
    python ${projectDir}/bin/synthesize_evaluations.py \\
        --ted-files ${output_prefix}_ted.tsv \\
        --flair-files ${output_prefix}_flair_eval.tsv \\
        --output ${output_prefix}_evaluation.tsv \\
        --test-name ${test_name}

    # Clean up intermediate files
    rm -f ${output_prefix}_ted.tsv ${output_prefix}_flair_eval.tsv
    """
}

process SummaryPlots {
    // Creates summary plots comparing all assemblers after evaluations complete
    // Currently generates:
    //   - 3-panel precision-recall comparison plot
    //     Panel 1: CAGE (5') Precision vs Recall
    //     Panel 2: QuantSeq (3') Precision vs Recall
    //     Panel 3: Reference 5' Precision vs Reference 3' Precision
    //   Points are colored by assembler (FLAIR, Bambu, IsoQuant)
    publishDir "results/evaluations/${test_name}/ted_plots", mode: 'copy'
    tag "${test_name}"

    input:
    // test_name: Used for output naming and publishDir
    // evaluation_files: All evaluation TSV files from Evaluation process
    tuple val(test_name), path(evaluation_files)

    output:
    // PNG image with 3-panel precision-recall comparison
    path "${test_name}_precision_recall_comparison.png", emit: precision_recall_plot

    script:
    """
    # Generate precision-recall comparison plot across all assemblers
    python ${projectDir}/bin/evaluation/precision_recall_plot.py \\
        --input ${evaluation_files} \\
        --output ${test_name}_precision_recall_comparison.png \\
        --title-prefix "${test_name}: " \\
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
    // Alternative transcriptome assemblers (optional - only run if defined in JSON)
    def bambu_modes = pipeline_config.bambu ?: [:]
    def isoquant_modes = pipeline_config.isoquant ?: [:]

    // Log alternative assembler status
    if (!bambu_modes.isEmpty()) {
        log.info "Bambu modes enabled: ${bambu_modes.keySet().join(', ')}"
    }
    if (!isoquant_modes.isEmpty()) {
        log.info "IsoQuant modes enabled: ${isoquant_modes.keySet().join(', ')}"
    }

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
                transcriptome: sanitizeModeMap(transcriptome_modes),
                bambu: sanitizeModeMap(bambu_modes),
                isoquant: sanitizeModeMap(isoquant_modes)
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
    // Each element: [test_name, dataset, align_modes, partition_modes, transcriptome_modes, bambu_modes, isoquant_modes]
    datasets_ch = Channel.from(test_sets_list)
        .map { test_set ->
            [test_set.name, test_set.dataset, test_set.alignModes, test_set.partitionModes,
             test_set.transcriptomeModes, test_set.bambuModes, test_set.isoquantModes]
        }

    // Split datasets into two branches based on whether BAM files are provided
    // Branch 1: Pre-aligned BAM files provided (skip FlairAlign)
    datasets_with_bam = datasets_ch.filter { test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes ->
        dataset.hasBam()
    }

    // Branch 2: Raw reads provided (run FlairAlign)
    datasets_without_bam = datasets_ch.filter { test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes ->
        !dataset.hasBam()
    }

    // Prepare inputs for FlairAlign process
    // Expands datasets into individual alignment jobs (one per read file per align mode)
    align_inputs = datasets_without_bam.flatMap { test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes ->
        // For each alignment mode configuration
        ds_align_modes.collectMany { align_mode, align_args ->
            // For each reads file (handles both single file and multi-file datasets)
            dataset.getReadsList().collect { reads_file ->
                [test_name, dataset.name, file(reads_file), align_mode, align_args, file(dataset.genome)]
            }
        }
    }

    // Prepare inputs for partition process from pre-aligned BAM (no BED conversion)
    // Uses NO_BED placeholder - FlairPartition will generate BED from partitioned BAM
    // This avoids expensive full BAM->BED conversion for large pre-aligned files
    prealigned_partition_inputs = datasets_with_bam.flatMap { test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes ->
        // Prepare CAGE and QuantSeq files (use placeholders if not provided)
        def cage_file = dataset.cage ? file(dataset.cage) : file("${workflow.workDir}/NO_CAGE", checkIfExists: false)
        def quantseq_file = dataset.quantseq ? file(dataset.quantseq) : file("${workflow.workDir}/NO_QUANTSEQ", checkIfExists: false)
        // Use NO_BED placeholder - FlairPartition will use --generate-bed flag
        def bed_placeholder = file("${workflow.workDir}/NO_BED", checkIfExists: false)
        // Generate tuples for each align_mode x partition_mode combination
        ds_align_modes.collectMany { align_mode, align_args ->
            ds_partition_modes.collect { partition_mode, partition_args ->
                [test_name, dataset.name, align_mode, file(dataset.bam), file(dataset.bai), bed_placeholder,
                 partition_mode, partition_args, file(dataset.genome), file(dataset.gtf),
                 cage_file, quantseq_file]
            }
        }
    }

    // =============================================================================
    // ALIGNMENT STAGE - Run alignment (only for datasets without pre-aligned BAM)
    // =============================================================================

    FlairAlign(align_inputs)

    // =============================================================================
    // PARTITION STAGE - Subset data to specific genomic regions
    // =============================================================================

    // For FlairAlign outputs: standard partition (BED already exists from alignment)
    partition_inputs_from_align = FlairAlign.out.alignments.combine(datasets_ch, by: 0).flatMap {
        test_name, dataset_name, align_mode, bam, bai, bed, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes ->
        // Prepare CAGE and QuantSeq files (use placeholders if not provided)
        def cage_file = dataset.cage ? file(dataset.cage) : file("${workflow.workDir}/NO_CAGE", checkIfExists: false)
        def quantseq_file = dataset.quantseq ? file(dataset.quantseq) : file("${workflow.workDir}/NO_QUANTSEQ", checkIfExists: false)
        // Generate tuple for each partition mode configuration
        ds_partition_modes.collect { partition_mode, partition_args ->
            [test_name, dataset_name, align_mode, bam, bai, bed, partition_mode, partition_args,
             file(dataset.genome), file(dataset.gtf), cage_file, quantseq_file]
        }
    }

    // Merge partition inputs from both branches (FlairAlign with BED, pre-aligned with NO_BED)
    // FlairPartition handles both cases - uses --bed or --generate-bed based on placeholder
    all_partition_inputs = partition_inputs_from_align.concat(prealigned_partition_inputs)

    FlairPartition(all_partition_inputs)

    // Use FlairPartition output directly (no merging needed anymore)
    all_partitioned = FlairPartition.out.partitioned

    // =============================================================================
    // TRANSCRIPTOME STAGE - Generate isoform predictions
    // =============================================================================

    // Prepare inputs for FlairTranscriptome by combining partitioned data with transcriptome modes
    // For each partitioned output, generate one transcriptome job per transcriptome mode
    transcriptome_inputs = all_partitioned.combine(datasets_ch, by: 0).flatMap {
        test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, quantseq_peaks,
        dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes ->
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
    // ALTERNATIVE ASSEMBLERS - Bambu and IsoQuant (conditional on mode definition)
    // =============================================================================

    // Prepare inputs for BambuAssembly
    // Only generates inputs if bambu modes are defined in params.json
    bambu_inputs = all_partitioned.combine(datasets_ch, by: 0).flatMap {
        test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, quantseq_peaks,
        dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes ->
        // Only generate inputs if bambu modes are defined
        if (ds_bambu_modes.isEmpty()) {
            return []
        }
        // Generate a tuple for each bambu mode configuration
        ds_bambu_modes.collect { bambu_mode, bambu_args ->
            [test_name, dataset_name, align_mode, partition_mode, bam, bai, genome, gtf,
             bambu_mode, bambu_args]
        }
    }

    // Run Bambu if there are any inputs
    BambuAssembly(bambu_inputs)

    // Prepare inputs for IsoQuantAssembly
    // Only generates inputs if isoquant modes are defined in params.json
    isoquant_inputs = all_partitioned.combine(datasets_ch, by: 0).flatMap {
        test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, quantseq_peaks,
        dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes ->
        // Only generate inputs if isoquant modes are defined
        if (ds_isoquant_modes.isEmpty()) {
            return []
        }
        // Generate a tuple for each isoquant mode configuration
        ds_isoquant_modes.collect { isoquant_mode, isoquant_args ->
            [test_name, dataset_name, align_mode, partition_mode, bam, bai, genome, gtf,
             isoquant_mode, isoquant_args]
        }
    }

    // Run IsoQuant if there are any inputs
    IsoQuantAssembly(isoquant_inputs)

    // =============================================================================
    // EVALUATION SETUP - Prepare reference data and organize inputs
    // =============================================================================

    // Extract GTF files from partitioned data to generate reference TSS/TTS peaks
    ref_peak_inputs = all_partitioned.map {
        test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, quantseq_peaks ->
        [test_name, dataset_name, align_mode, partition_mode, gtf]
    }

    PrepareReferencePeaks(ref_peak_inputs)

    // =============================================================================
    // EVALUATION STAGE - Unified evaluation for all assemblers
    // =============================================================================
    //
    // All assembler outputs (FLAIR, Bambu, IsoQuant) are processed through the
    // unified Evaluation process. The process automatically detects the assembler
    // type and runs appropriate evaluation (full for FLAIR, simplified for others).
    //
    // Input tuple format for Evaluation:
    // [test_name, dataset_name, align_mode, partition_mode, transcriptome_mode,
    //  isoforms_bed, isoforms_gtf, isoform_read_map,
    //  bam, bai, reads_bed, genome, gtf,
    //  cage_peaks, quantseq_peaks, ref_tss, ref_tts]

    // Prepare FLAIR transcriptome outputs for evaluation
    flair_eval_inputs = FlairTranscriptome.out.transcriptome
        .map { test_name, dataset_name, align_mode, partition_mode, partition_args, transcriptome_mode,
               isoforms_bed, isoforms_gtf, isoforms_fa, isoform_counts, isoform_read_map ->
            // FLAIR outputs: use isoforms_bed and isoform_read_map, placeholder for isoforms_gtf
            [test_name, dataset_name, align_mode, partition_mode, transcriptome_mode,
             isoforms_bed, file("${workflow.workDir}/NO_ISOFORMS_GTF", checkIfExists: false), isoform_read_map]
        }
        .combine(all_partitioned.map { test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, quantseq_peaks ->
            [test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, quantseq_peaks]
        }, by: [0, 1, 2, 3])
        .combine(PrepareReferencePeaks.out.ref_peaks, by: [0, 1, 2, 3])
        .map { test_name, dataset_name, align_mode, partition_mode, transcriptome_mode,
               isoforms_bed, isoforms_gtf_placeholder, isoform_read_map,
               bam, bai, reads_bed, genome, gtf, cage_peaks, quantseq_peaks,
               ref_tss, ref_tts ->
            def cage_file = cage_peaks ?: file("${workflow.workDir}/NO_CAGE", checkIfExists: false)
            def quantseq_file = quantseq_peaks ?: file("${workflow.workDir}/NO_QUANTSEQ", checkIfExists: false)
            [test_name, dataset_name, align_mode, partition_mode, transcriptome_mode,
             isoforms_bed, isoforms_gtf_placeholder, isoform_read_map,
             bam, bai, reads_bed, genome, gtf,
             cage_file, quantseq_file, ref_tss, ref_tts]
        }

    // Prepare Bambu outputs for evaluation
    bambu_eval_inputs = BambuAssembly.out.bambu_gtf
        .combine(all_partitioned.map { test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, quantseq_peaks ->
            [test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, quantseq_peaks]
        }, by: [0, 1, 2, 3])
        .combine(PrepareReferencePeaks.out.ref_peaks, by: [0, 1, 2, 3])
        .map { test_name, dataset_name, align_mode, partition_mode, bambu_mode, isoforms_gtf,
               bam, bai, reads_bed, genome, ref_gtf, cage_peaks, quantseq_peaks, ref_tss, ref_tts ->
            def cage_file = cage_peaks ?: file("${workflow.workDir}/NO_CAGE", checkIfExists: false)
            def quantseq_file = quantseq_peaks ?: file("${workflow.workDir}/NO_QUANTSEQ", checkIfExists: false)
            // Bambu outputs: use isoforms_gtf, placeholders for isoforms_bed and isoform_read_map
            [test_name, dataset_name, align_mode, partition_mode, "bambu_${bambu_mode}",
             file("${workflow.workDir}/NO_ISOFORMS_BED", checkIfExists: false), isoforms_gtf,
             file("${workflow.workDir}/NO_READ_MAP", checkIfExists: false),
             bam, bai, reads_bed, genome, ref_gtf,
             cage_file, quantseq_file, ref_tss, ref_tts]
        }

    // Prepare IsoQuant outputs for evaluation
    isoquant_eval_inputs = IsoQuantAssembly.out.isoquant_gtf
        .combine(all_partitioned.map { test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, quantseq_peaks ->
            [test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, quantseq_peaks]
        }, by: [0, 1, 2, 3])
        .combine(PrepareReferencePeaks.out.ref_peaks, by: [0, 1, 2, 3])
        .map { test_name, dataset_name, align_mode, partition_mode, isoquant_mode, isoforms_gtf,
               bam, bai, reads_bed, genome, ref_gtf, cage_peaks, quantseq_peaks, ref_tss, ref_tts ->
            def cage_file = cage_peaks ?: file("${workflow.workDir}/NO_CAGE", checkIfExists: false)
            def quantseq_file = quantseq_peaks ?: file("${workflow.workDir}/NO_QUANTSEQ", checkIfExists: false)
            // IsoQuant outputs: use isoforms_gtf, placeholders for isoforms_bed and isoform_read_map
            [test_name, dataset_name, align_mode, partition_mode, "isoquant_${isoquant_mode}",
             file("${workflow.workDir}/NO_ISOFORMS_BED", checkIfExists: false), isoforms_gtf,
             file("${workflow.workDir}/NO_READ_MAP", checkIfExists: false),
             bam, bai, reads_bed, genome, ref_gtf,
             cage_file, quantseq_file, ref_tss, ref_tts]
        }

    // Merge all evaluation inputs and run unified Evaluation process
    all_eval_inputs = flair_eval_inputs.concat(bambu_eval_inputs).concat(isoquant_eval_inputs)
    Evaluation(all_eval_inputs)

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
        .combine(all_partitioned.map { test_name, dataset_name, align_mode, partition_mode,
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
    // SUMMARY PLOTS - Compare all assemblers
    // =============================================================================

    // Collect all evaluation TSVs and group by test_name for summary plots
    all_evaluations = Evaluation.out.evaluation_results
        .map { test_name, dataset_name, align_mode, partition_mode, transcriptome_mode, eval_tsv ->
            [test_name, eval_tsv]
        }
        .groupTuple()

    SummaryPlots(all_evaluations)

}
