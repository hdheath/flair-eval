#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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
    String junctions
    String cage
    String quantseq
    String bam  // Optional: pre-aligned BAM file
    String bai  // Optional: BAM index file (auto-constructed from bam path)

    // Constructor
    Dataset(String name, Map config) {
        this.name = name
        this.reads = config.reads
        this.genome = config.genome
        this.gtf = config.gtf
        this.junctions = config.junctions
        this.cage = config.cage
        this.quantseq = config.quantseq
        this.bam = config.bam
        // Auto-construct BAI path from BAM path if BAM is provided
        this.bai = config.bam ? "${config.bam}.bai" : null
    }

    // Helper methods
    boolean hasJunctions() { return junctions != null }
    boolean hasCage() { return cage != null }
    boolean hasQuantseq() { return quantseq != null }
    boolean hasBam() { return bam != null }

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
// PARAMETER DEFINITIONS - global
// =============================================================================

def standard_align_modes = [
    default: ''
]

def standard_partition_modes = [
    //'SMARCA4': '--region chr19:10900001-11100000'
    'all': '--all'
]

def standard_transcriptome_modes = [
    'with-gtf': ''
]

// =============================================================================
// PROCESS DEFINITIONS
// =============================================================================

process FlairAlign {
    publishDir "results/align", mode: 'symlink'
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}"

    input:
    tuple val(test_name), val(dataset_name), path(reads), val(align_mode), val(align_args), path(genome)

    output:
    tuple val(test_name), val(dataset_name), val(align_mode),
          path("${dataset_name}_${align_mode}.bam"),
          path("${dataset_name}_${align_mode}.bam.bai"),
          path("${dataset_name}_${align_mode}.bed"), emit: alignments

    script:
    """
    flair align ${align_args} \\
        -r ${reads} \\
        -g ${genome} \\
        -o ${dataset_name}_${align_mode}
    """
}

process BamToBed {
    publishDir "results/align", mode: 'symlink'
    tag "${dataset_name}_${align_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), path(bam), path(bai)

    output:
    tuple val(test_name), val(dataset_name), val(align_mode),
          path(bam), path(bai), path("${dataset_name}_${align_mode}.bed"), emit: alignments

    script:
    """
    # Generate BED file from BAM
    bedtools bamtobed -bed12 -i ${bam} > ${dataset_name}_${align_mode}.bed
    """
}

process FlairPartition {
    publishDir "results/partition/${test_name}", mode: 'symlink'
    tag "${dataset_name}_${align_mode}_${partition_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), path(bam), path(bai), path(bed),
          val(partition_mode), val(partition_args), path(genome), path(gtf)

    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}.bam"),
          path("${dataset_name}_${align_mode}_${partition_mode}.bam.bai"),
          path("${dataset_name}_${align_mode}_${partition_mode}.bed"),
          path("${dataset_name}_${align_mode}_${partition_mode}_genome.fa"),
          path("${dataset_name}_${align_mode}_${partition_mode}_annotation.gtf"), emit: partitioned

    script:
    def output_prefix = "${dataset_name}_${align_mode}_${partition_mode}"

    """
    python ${projectDir}/bin/simple_partition.py \\
        --bam ${bam} \\
        --bed ${bed} \\
        --genome ${genome} \\
        --gtf ${gtf} \\
        --output-prefix ${output_prefix} \\
        ${partition_args}
    """
}

process FlairTranscriptome {
    publishDir "results/transcriptome/${test_name}", mode: 'symlink'
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_transcriptome"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), val(partition_args),
          path(bam), path(bai), path(bed), path(genome), path(gtf),
          val(transcriptome_mode), val(transcriptome_args)

    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), val(partition_args), val(transcriptome_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_transcriptome.isoforms.bed"),
          path("${dataset_name}_${align_mode}_${partition_mode}_transcriptome.isoforms.gtf"),
          path("${dataset_name}_${align_mode}_${partition_mode}_transcriptome.isoforms.fa"),
          path("${dataset_name}_${align_mode}_${partition_mode}_transcriptome.isoform.counts.txt"),
          path("${dataset_name}_${align_mode}_${partition_mode}_transcriptome.isoform.read.map.txt"), emit: transcriptome

    script:
    def extra_args = transcriptome_args ? "${transcriptome_args} \\" : ""
    """
    flair transcriptome \\
        -b ${bam} \\
        --genome ${genome} \\
        -f ${gtf} \\
        ${extra_args}        -o ${dataset_name}_${align_mode}_${partition_mode}_transcriptome
    """
}

process FlairEvaluation {
    publishDir "results/evaluations/individual/${test_name}", mode: 'symlink', pattern: '*.tsv'
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_${process_label}_${stage}"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(partition_args), val(process_label), val(stage),
          path(isoforms_bed), path(isoform_read_map), path(bam), path(bai),
          path(reads_bed), path(corrected_bed), path(genome), path(gtf),
          path(cage_peaks), path(quantseq_peaks),
          path(ref_tss), path(ref_tts)

    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(partition_args), val(process_label), val(stage),
          path("${dataset_name}_${align_mode}_${partition_mode}_${process_label}_${stage}_ted.tsv"),
          path("${dataset_name}_${align_mode}_${partition_mode}_${process_label}_${stage}_flair_eval.tsv"), emit: evaluation_results

    script:
    // Build optional arguments for TED
    def cage_arg = cage_peaks.name != 'NO_CAGE' ? "--prime5-peaks ${cage_peaks}" : ""
    def quantseq_arg = quantseq_peaks.name != 'NO_QUANTSEQ' ? "--prime3-peaks ${quantseq_peaks}" : ""
    def ref_tss_arg = ref_tss.name != 'NO_REF_TSS' ? "--ref-prime5-peaks ${ref_tss}" : ""
    def ref_tts_arg = ref_tts.name != 'NO_REF_TTS' ? "--ref-prime3-peaks ${ref_tts}" : ""
    def corrected_arg = corrected_bed.name != 'NO_CORRECTED' ? "--corrected-bed ${corrected_bed}" : ""

    """
    # Run TED evaluation
    python ${projectDir}/bin/ted.py \\
        --isoforms-bed ${isoforms_bed} \\
        --map-file ${isoform_read_map} \\
        --bam ${bam} \\
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
        --output ${dataset_name}_${align_mode}_${partition_mode}_${process_label}_${stage}_ted.tsv \\
        --verbose

    # Run FLAIR evaluation
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
    publishDir "results/reference_peaks/${test_name}", mode: 'symlink'
    tag "${dataset_name}_${align_mode}_${partition_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), path(gtf)

    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_ref_tss.bed"),
          path("${dataset_name}_${align_mode}_${partition_mode}_ref_tts.bed"), emit: ref_peaks

    script:
    """
    python ${projectDir}/bin/gtf_to_tss_tts.py \\
        --gtf ${gtf} \\
        --output-prefix ${dataset_name}_${align_mode}_${partition_mode}_ref \\
        --deduplicate
    """
}

process SynthesizeEvaluations {
    publishDir "results/evaluations/summary", mode: 'symlink'
    tag "${test_name}_evaluation_summary"

    input:
    tuple val(test_name), path(ted_files), path(eval_files)

    output:
    path "${test_name}_evaluation_summary.tsv", emit: evaluation_summary

    script:
    """
    python ${projectDir}/bin/synthesize_evaluations.py \\
        --ted-files ${ted_files.join(' ')} \\
        --flair-files ${eval_files.join(' ')} \\
        --output ${test_name}_evaluation_summary.tsv \\
        --test-name ${test_name}
    """
}

process PlotIsoforms {
    publishDir "results/isoform_plots/${test_name}", mode: 'copy'
    conda '/private/home/hdheath/miniforge3/envs/nextflow_env'
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path(bam), path(bai), path(isoforms_bed), path(isoform_read_map)

    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_isoform_plot.png"), emit: plots

    script:
    """
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
        error "ERROR: Please provide a samplesheet via --input"
    }

    // Display test name warning if using default
    if (params.test_name == 'flair_test_suite') {
        log.warn "WARNING: Using default test name 'flair_test_suite'. Specify a custom name with --test_name"
    }

    // Read and parse samplesheet CSV into a list
    def test_sets_list = []
    file(params.input).withReader { reader ->
        def header = reader.readLine().split(',')
        reader.eachLine { line ->
            def values = line.split(',')
            def row = [header, values].transpose().collectEntries()

            // Create Dataset object from CSV row
            def dataset = new Dataset(row.sample_id, ([
                genome: file(row.genome),
                gtf: file(row.gtf),
                bam: row.bam && row.bam != '' ? file(row.bam) : null,
                reads: row.reads && row.reads != '' ? file(row.reads) : null,
                cage: row.cage && row.cage != '' ? file(row.cage) : null,
                quantseq: row.quantseq && row.quantseq != '' ? file(row.quantseq) : null
            ]))

            // Create TestSet with the specified test_name
            test_sets_list.add(new TestSet("${params.test_name}_${dataset.name}", dataset, ([
                align: sanitizeModeMap(standard_align_modes),
                partition: sanitizeModeMap(standard_partition_modes),
                transcriptome: sanitizeModeMap(standard_transcriptome_modes)
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

    // Build channels
    datasets_ch = Channel.from(test_sets_list)
        .map { test_set -> 
            [test_set.name, test_set.dataset, test_set.alignModes, test_set.partitionModes,
             test_set.transcriptomeModes]
        }

    datasets_with_bam = datasets_ch.filter { test_name, dataset, align_modes, partition_modes, transcriptome_modes ->
        dataset.hasBam()
    }

    datasets_without_bam = datasets_ch.filter { test_name, dataset, align_modes, partition_modes, transcriptome_modes ->
        !dataset.hasBam()
    }

    align_inputs = datasets_without_bam.flatMap { test_name, dataset, align_modes, partition_modes, transcriptome_modes ->
        align_modes.collectMany { align_mode, align_args ->
            dataset.getReadsList().collect { reads_file ->
                [test_name, dataset.name, file(reads_file), align_mode, align_args, file(dataset.genome)]
            }
        }
    }

    prealigned_bam_inputs = datasets_with_bam.flatMap { test_name, dataset, align_modes, partition_modes, transcriptome_modes ->
        align_modes.collect { align_mode, align_args ->
            [test_name, dataset.name, align_mode, file(dataset.bam), file(dataset.bai)]
        }
    }

    // Run processes
    FlairAlign(align_inputs)
    BamToBed(prealigned_bam_inputs)

    // Mix outputs - use concat instead of mix for better empty channel handling
    all_alignments = FlairAlign.out.alignments.concat(BamToBed.out.alignments)

    partition_inputs = all_alignments.combine(datasets_ch, by: 0).flatMap { 
        test_name, dataset_name, align_mode, bam, bai, bed, dataset, align_modes, partition_modes, transcriptome_modes ->
        partition_modes.collect { partition_mode, partition_args ->
            [test_name, dataset_name, align_mode, bam, bai, bed, partition_mode, partition_args, file(dataset.genome), file(dataset.gtf)]
        }
    }

    FlairPartition(partition_inputs)

    transcriptome_inputs = FlairPartition.out.partitioned.combine(datasets_ch, by: 0).flatMap { 
        test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, 
        dataset, align_modes, partition_modes, transcriptome_modes ->
        def partition_args = partition_modes[partition_mode] ?: ''
        transcriptome_modes.collect { transcriptome_mode, transcriptome_args ->
            [test_name, dataset_name, align_mode, partition_mode, partition_args, bam, bai, bed, genome, gtf, 
             transcriptome_mode, transcriptome_args]
        }
    }

    FlairTranscriptome(transcriptome_inputs)

    // Evaluation logic
    ref_peak_inputs = FlairPartition.out.partitioned.map { 
        test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf ->
        [test_name, dataset_name, align_mode, partition_mode, gtf]
    }

    PrepareReferencePeaks(ref_peak_inputs)

    align_for_eval = all_alignments.map { test_name, dataset_name, align_mode, bam, bai, bed ->
        [test_name, dataset_name, align_mode, bam, bai, bed]
    }

    cage_quantseq_by_dataset = datasets_ch.map { 
        test_name, dataset, align_modes, partition_modes, transcriptome_modes ->
        def cage_file = dataset.cage ? file(dataset.cage) : file("${workflow.workDir}/NO_CAGE", checkIfExists: false)
        def quantseq_file = dataset.quantseq ? file(dataset.quantseq) : file("${workflow.workDir}/NO_QUANTSEQ", checkIfExists: false)
        [test_name, dataset.name, cage_file, quantseq_file]
    }

    eval_transcriptome_inputs = FlairTranscriptome.out.transcriptome
        .map { test_name, dataset_name, align_mode, partition_mode, partition_args, transcriptome_mode,
               isoforms_bed, isoforms_gtf, isoforms_fa, isoform_counts, isoform_read_map ->
            [test_name, dataset_name, align_mode, partition_mode, partition_args, "transcriptome_${transcriptome_mode}", 
             "transcriptome", isoforms_bed, isoform_read_map, file("${workflow.workDir}/NO_CORRECTED", checkIfExists: false)]
        }
        .combine(FlairPartition.out.partitioned.map { test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf ->
            [test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf]
        }, by: [0, 1, 2, 3])
        .combine(datasets_ch, by: 0)
        .combine(cage_quantseq_by_dataset, by: [0, 1])
        .combine(PrepareReferencePeaks.out.ref_peaks, by: [0, 1, 2, 3])
        .map { test_name, dataset_name, align_mode, partition_mode, partition_args, process_label, stage,
               isoforms_bed, isoform_read_map, corrected_file, bam, bai, reads_bed, genome, gtf,
               dataset, align_modes, partition_modes, transcriptome_modes,
               cage_file, quantseq_file, ref_tss, ref_tts ->
            [test_name, dataset_name, align_mode, partition_mode, partition_args, process_label, stage,
             isoforms_bed, isoform_read_map, bam, bai, reads_bed, corrected_file, genome, gtf,
             cage_file, quantseq_file, ref_tss, ref_tts]
        }

    FlairEvaluation(eval_transcriptome_inputs)

    // Prepare inputs for plotting (filter for regions ≤ 400kb)
    plot_inputs = FlairTranscriptome.out.transcriptome
        .map { test_name, dataset_name, align_mode, partition_mode, partition_args, transcriptome_mode,
               isoforms_bed, isoforms_gtf, isoforms_fa, isoform_counts, isoform_read_map ->
            // Extract region size from partition_args if present
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
        .filter { test_name, dataset_name, align_mode, partition_mode, partition_args, region_size,
                  isoforms_bed, isoform_read_map ->
            // Only plot if region is 400kb or less
            region_size > 0 && region_size <= 400000
        }
        .combine(FlairPartition.out.partitioned.map { test_name, dataset_name, align_mode, partition_mode, 
                                                        bam, bai, bed, genome, gtf ->
            [test_name, dataset_name, align_mode, partition_mode, bam, bai]
        }, by: [0, 1, 2, 3])
        .map { test_name, dataset_name, align_mode, partition_mode, partition_args, region_size,
               isoforms_bed, isoform_read_map, bam, bai ->
            [test_name, dataset_name, align_mode, partition_mode, bam, bai, isoforms_bed, isoform_read_map]
        }

    PlotIsoforms(plot_inputs)

    // Collect all TED and FLAIR evaluation files
    all_ted_files = FlairEvaluation.out.evaluation_results
        .map { test_name, dataset_name, align_mode, partition_mode, partition_args, process_label, stage, ted_file, flair_file ->
            ted_file
        }
        .collect()

    all_flair_files = FlairEvaluation.out.evaluation_results
        .map { test_name, dataset_name, align_mode, partition_mode, partition_args, process_label, stage, ted_file, flair_file ->
            flair_file
        }
        .collect()

    // Combine all files for a single comprehensive summary
    all_results_input = all_ted_files
        .map { ted_files -> [test_set_name, ted_files] }
        .join(all_flair_files.map { flair_files -> [test_set_name, flair_files] })

    // Generate one comprehensive summary with all TED and FLAIR evaluations
    SynthesizeEvaluations(all_results_input)
}
