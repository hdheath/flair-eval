#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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
    
    // Constructor
    Dataset(String name, Map config) {
        this.name = name
        this.reads = config.reads
        this.genome = config.genome
        this.gtf = config.gtf
        this.junctions = config.junctions
        this.cage = config.cage
        this.quantseq = config.quantseq
    }
    
    // Helper methods
    boolean hasJunctions() { return junctions != null }
    boolean hasCage() { return cage != null }
    boolean hasQuantseq() { return quantseq != null }
    
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
    Map correctModes
    Map collapseModes
    Map transcriptomeModes
    
    // Constructor with simplified parameter structure
    TestSet(String name, Dataset dataset, Map modes) {
        this.name = name
        this.dataset = dataset
        this.alignModes = modes.align ?: [:]
        this.partitionModes = modes.partition ?: [:]
        this.correctModes = modes.correct ?: [:]
        this.collapseModes = modes.collapse ?: [:]
        this.transcriptomeModes = modes.transcriptome ?: [:]
    }
    
    // Calculate total number of jobs this test set will produce
    int totalJobs() {
        def align_count = alignModes.size() ?: 1
        def partition_count = partitionModes.size() ?: 1
        def correct_count = correctModes.size() ?: 1
        def collapse_count = collapseModes.size() ?: 1
        
        return align_count * partition_count * correct_count * collapse_count
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
    all: '--all',
    chr1: '--region chr1',
    chr5: '--region chr5',
    chr22: '--region chr22'
]

def standard_correct_modes = [
    with_gtf: '',
]

def standard_collapse_modes = [
    default: '',
]

def standard_transcriptome_modes = [
    with_gtf: ''
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

process FlairCorrect {
    publishDir "results/correct/${test_name}", mode: 'symlink'
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_${correct_mode}"
    
    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), 
          path(bam), path(bai), path(bed), path(genome), path(gtf),
          val(correct_mode), val(correct_args)
    
    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), val(correct_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_${correct_mode}_all_corrected.bed", optional: true),
          path("${dataset_name}_${align_mode}_${partition_mode}_${correct_mode}_all_inconsistent.bed", optional: true), emit: corrected
    
    script:
    """
    flair correct ${correct_args} \\
        -q ${bed} \\
        -f ${gtf} \\
        -o ${dataset_name}_${align_mode}_${partition_mode}_${correct_mode}
    """
}

process FlairCollapse {
    publishDir "results/collapse/${test_name}", mode: 'symlink'
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_${correct_mode}_${collapse_mode}"
    
    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), val(correct_mode),
          path(corrected_bed), path(inconsistent_bed), path(reads),
          path(genome), path(gtf), val(collapse_mode), val(collapse_args)
    
    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), val(correct_mode), val(collapse_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_${correct_mode}_${collapse_mode}.isoforms.bed"),
          path("${dataset_name}_${align_mode}_${partition_mode}_${correct_mode}_${collapse_mode}.isoforms.fa"),
          path("${dataset_name}_${align_mode}_${partition_mode}_${correct_mode}_${collapse_mode}.isoforms.gtf"),
          path("${dataset_name}_${align_mode}_${partition_mode}_${correct_mode}_${collapse_mode}.isoform.read.map.txt"),
          path("${dataset_name}_${align_mode}_${partition_mode}_${correct_mode}_${collapse_mode}.isoform.counts.txt"), emit: collapsed
    
    script:
    """
    flair collapse ${collapse_args} \\
        -r ${reads} \\
        -q ${corrected_bed} \\
        -g ${genome} \\
        -f ${gtf} \\
        --generate_map \\
        -o ${dataset_name}_${align_mode}_${partition_mode}_${correct_mode}_${collapse_mode}
    """
}

process FlairTranscriptome {
    publishDir "results/transcriptome/${test_name}", mode: 'symlink'
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_transcriptome"
    
    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path(bam), path(bai), path(bed), path(genome), path(gtf),
          val(transcriptome_mode), val(transcriptome_args)
    
    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), val(transcriptome_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_transcriptome.isoforms.bed"),
          path("${dataset_name}_${align_mode}_${partition_mode}_transcriptome.isoforms.gtf"),
          path("${dataset_name}_${align_mode}_${partition_mode}_transcriptome.isoforms.fa"),
          path("${dataset_name}_${align_mode}_${partition_mode}_transcriptome.isoform.counts.txt"),
          path("${dataset_name}_${align_mode}_${partition_mode}_transcriptome.isoform.read.map.txt"), emit: transcriptome
    
    script:
    """
    flair transcriptome \\
        -b ${bam} \\
        --genome ${genome} \\
        -f ${gtf} \\
        ${transcriptome_args} \\
        -o ${dataset_name}_${align_mode}_${partition_mode}_transcriptome
    """
}

process FlairEval {
    publishDir "results/eval/${test_name}", mode: 'symlink'
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_${process_label}"
    
    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), val(process_label),
          path(isoforms_bed), path(reads_bed), path(gtf)
    
    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), val(process_label),
          path("${dataset_name}_${align_mode}_${partition_mode}_${process_label}.eval_summary.txt"), emit: eval_results
    
    script:
    """
    python ${projectDir}/bin/flair_eval.py \\
        --reads-bed ${reads_bed} \\
        --isoforms-bed ${isoforms_bed} \\
        --gtf ${gtf} \\
        --output ${dataset_name}_${align_mode}_${partition_mode}_${process_label}.eval_summary.txt \\
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

process FlairTED {
    publishDir "results/ted/${test_name}/${dataset_name}", mode: 'symlink', 
               pattern: "*.tsv"
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_${process_label}_${stage}"
    
    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), 
          val(process_label), val(stage),
          path(isoforms_bed), path(isoform_read_map), path(bam), path(bai),
          path(corrected_bed), path(cage_peaks), path(quantseq_peaks), 
          path(ref_tss), path(ref_tts)
    
    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), 
          val(process_label), val(stage),
          path("${dataset_name}_${align_mode}_${partition_mode}_${process_label}_${stage}_ted.tsv"), emit: ted_metrics
    
    script:
    // Build optional arguments
    def cage_arg = cage_peaks.name != 'NO_CAGE' ? "--prime5-peaks ${cage_peaks}" : ""
    def quantseq_arg = quantseq_peaks.name != 'NO_QUANTSEQ' ? "--prime3-peaks ${quantseq_peaks}" : ""
    def ref_tss_arg = ref_tss.name != 'NO_REF_TSS' ? "--ref-prime5-peaks ${ref_tss}" : ""
    def ref_tts_arg = ref_tts.name != 'NO_REF_TTS' ? "--ref-prime3-peaks ${ref_tts}" : ""
    def corrected_arg = corrected_bed.name != 'NO_CORRECTED' ? "--corrected-bed ${corrected_bed}" : ""
    
    """
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
        --output ${dataset_name}_${align_mode}_${partition_mode}_${process_label}_${stage}_ted.tsv \\
        --verbose
    """
}

process SynthesizeEvaluations {
    publishDir "results/summary", mode: 'symlink'
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

workflow {
    // =============================================================================
    // DATASET DEFINITIONS
    // =============================================================================
    
    def a549_chr1_dataset = new Dataset('a549_chr1', [
        reads: '/private/groups/brookslab/hdheath/projects/flair-eval/data/A549_chr1_reads_300kb.fasta',
        genome: '/private/groups/brookslab/hdheath/projects/flair-eval/data/A549_chr1_reference.fasta',
        gtf: '/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/gencode.v48.annotation.gtf'
    ])

    def a549_chr5_dataset = new Dataset('a549_chr5', [
        reads: '/private/groups/brookslab/hdheath/projects/flair-eval/data/A549_chr5_reads_300kb.fasta',
        genome: '/private/groups/brookslab/hdheath/projects/flair-eval/data/A549_chr5_reference.fasta',
        gtf: '/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/gencode.v48.annotation.gtf'
    ])

    def a549_chr22_dataset = new Dataset('a549_chr22', [
        reads: '/private/groups/brookslab/cafelton/testflairanyvcf/simNMD/smallchr22test/WTC11.pacBio.chr22smalllocus.fasta',
        genome: '/private/groups/brookslab/cafelton/testflairanyvcf/simNMD/smallchr22test/GRCh38.chr22.genome.fa',
        gtf: '/private/groups/brookslab/cafelton/testflairanyvcf/simNMD/smallchr22test/gencode.v38.annotation.chr22.gtf'
    ])
    
    // =============================================================================
    // TEST SET DEFINITIONS 1 - Single Dataset, individual testsets, custom options
    // =============================================================================
    
    // // Quick test - minimal options
    // def quick_test = new TestSet('test_w_colette', test_w_colette, [
    //     align: [default: ''],
    //     partition: [chr22_100k: '--region chr22:31000000-32000000'],
    //     correct: [with_gtf: ''],
    //     collapse: [default: ''],
    //     transcriptome: [
    //         with_gtf: '',
    //         max_ends: '--max_ends 3'
    //     ]
    // ])

    // def quick_test2 = new TestSet('test_w_colette', test_w_colette, [
    //     align: [default: '']
    // ])

    // // Partition comparison test - focused on partition modes
    // def partition_test = new TestSet('partition_test', a549_dataset, [
    //     align: [default: ''],
    //     partition: [
    //         all: '--all',
    //         chr1_100k: '--region chr1:1-100000'
    //     ],
    //     correct: [with_gtf: ''],
    //     collapse: [default: ''],
    //     transcriptome: [with_gtf: '']
    // ])

    // =============================================================================
    // TEST SET DEFINITIONS 2 - Single Dataset, individual testsets, global options
    // =============================================================================
    
    // // Comprehensive tests - uses standard modes for full matrix across multiple datasets
    // def comprehensive_chr1_test = new TestSet('comprehensive_chr1', a549_chr1_dataset, [
    //     align: standard_align_modes,
    //     partition: standard_partition_modes,
    //     correct: standard_correct_modes,
    //     collapse: standard_collapse_modes,
    //     transcriptome: standard_transcriptome_modes
    // ])
    
    // def comprehensive_chr5_test = new TestSet('comprehensive_chr5', a549_chr5_dataset, [
    //     align: standard_align_modes,
    //     partition: standard_partition_modes,
    //     correct: standard_correct_modes,
    //     collapse: standard_collapse_modes,
    //     transcriptome: standard_transcriptome_modes
    // ])
    
    // def comprehensive_chr22_test = new TestSet('comprehensive_chr22', a549_chr22_dataset, [
    //     align: standard_align_modes,
    //     partition: standard_partition_modes,
    //     correct: standard_correct_modes,
    //     collapse: standard_collapse_modes,
    //     transcriptome: standard_transcriptome_modes
    // ])

    // =============================================================================
    // TEST SET DEFINITIONS 3 - Multi Dataset, using global options
    // =============================================================================
    
    // Helper function to create a testset of multiple datasets with global options
    def createMultiDatasetRun = { datasets ->
        return datasets.collect { dataset ->
            new TestSet("comprehensive_${dataset.name}", dataset, [
                align: standard_align_modes,
                partition: standard_partition_modes,
                correct: standard_correct_modes,
                collapse: standard_collapse_modes,
                transcriptome: standard_transcriptome_modes
            ])
        }
    }
    
    // Create comprehensive tests for all chromosome datasets
    def all_comprehensive_tests = createMultiDatasetRun([
        a549_chr1_dataset, 
        a549_chr5_dataset, 
        a549_chr22_dataset
    ])
    
    // =============================================================================
    // COLLECT ALL TEST SETS TO RUN
    // =============================================================================
    
    def test_sets = all_comprehensive_tests
    
    // Alternative ways to configure test_sets:
    // def test_sets = [quick_test]
    // def test_sets = [quick_test, partition_test]
    // def test_sets = [comprehensive_chr1_test, comprehensive_chr5_test, comprehensive_chr22_test]
    // def test_sets = [quick_test] + all_comprehensive_tests
    
    // Print summary
    println("\n=== FLAIR Test Suite ===")
    test_sets.each { test ->
        println("${test.name}: ${test.totalJobs()} jobs using dataset '${test.dataset.name}'")
    }
    println("Total jobs: ${test_sets.sum { it.totalJobs() }}")
    println("===================================\n")
    
    // =============================================================================
    // BUILD CHANNELS
    // =============================================================================
    
    // Create datasets channel
    datasets_ch = Channel.from(test_sets)
        .map { test_set -> 
            [test_set.name, test_set.dataset, test_set.alignModes, test_set.partitionModes, 
             test_set.correctModes, test_set.collapseModes, test_set.transcriptomeModes]
        }
    
    // Generate align combinations (include dataset-specific genome)
    align_inputs = datasets_ch
        .flatMap { test_name, dataset, align_modes, partition_modes, correct_modes, collapse_modes, transcriptome_modes ->
            align_modes.collectMany { align_mode, align_args ->
                dataset.getReadsList().collect { reads_file ->
                    [test_name, dataset.name, file(reads_file), align_mode, align_args, file(dataset.genome)]
                }
            }
        }
    
    // Run alignment (genome files are now included in align_inputs)
    FlairAlign(align_inputs)
    
    // Generate partition combinations using combine (include dataset-specific genome)
    partition_inputs = FlairAlign.out.alignments
        .combine(datasets_ch, by: 0)  // Join by test_name
        .flatMap { test_name, dataset_name, align_mode, bam, bai, bed, dataset, align_modes, partition_modes, correct_modes, collapse_modes, transcriptome_modes ->
            partition_modes.collect { partition_mode, partition_args ->
                [test_name, dataset_name, align_mode, bam, bai, bed, partition_mode, partition_args, file(dataset.genome), file(dataset.gtf)]
            }
        }
    
    // Run partition (genome and GTF files are now included in partition_inputs)
    FlairPartition(partition_inputs)
    
    // Generate transcriptome combinations
    transcriptome_inputs = FlairPartition.out.partitioned
        .combine(datasets_ch, by: 0)  // Join by test_name
        .flatMap { test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, 
                   dataset, align_modes, partition_modes, correct_modes, collapse_modes, transcriptome_modes ->
            transcriptome_modes.collect { transcriptome_mode, transcriptome_args ->
                [test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, 
                 transcriptome_mode, transcriptome_args]
            }
        }
    
    // Run transcriptome
    FlairTranscriptome(transcriptome_inputs)
    
    // Generate correct combinations
    correct_inputs = FlairPartition.out.partitioned
        .combine(datasets_ch, by: 0)  // Join by test_name  
        .flatMap { test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf,
                   dataset, align_modes, partition_modes, correct_modes, collapse_modes, transcriptome_modes ->
            correct_modes.collect { correct_mode, correct_args ->
                [test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf,
                 correct_mode, correct_args]
            }
        }
    
    // Run correct
    FlairCorrect(correct_inputs)
    
    // Generate collapse combinations
    collapse_inputs = FlairCorrect.out.corrected
        .filter { test_name, dataset_name, align_mode, partition_mode, correct_mode, corrected_bed, inconsistent_bed ->
            corrected_bed != null && corrected_bed.size() > 0
        }
        .combine(datasets_ch, by: 0)  // Join by test_name
        .flatMap { test_name, dataset_name, align_mode, partition_mode, correct_mode, corrected_bed, inconsistent_bed,
                   dataset, align_modes, partition_modes, correct_modes, collapse_modes, transcriptome_modes ->
            collapse_modes.collect { collapse_mode, collapse_args ->
                [test_name, dataset_name, align_mode, partition_mode, correct_mode, corrected_bed, inconsistent_bed,
                 file(dataset.reads), file(dataset.genome), file(dataset.gtf), collapse_mode, collapse_args]
            }
        }
    
    // Run collapse
    FlairCollapse(collapse_inputs)
    
    // =============================================================================
    // EVALUATION LOGIC
    // =============================================================================
    
    // Prepare align outputs for evaluation (get reads BED files)
    align_for_eval = FlairAlign.out.alignments
        .map { test_name, dataset_name, align_mode, bam, bai, bed ->
            [test_name, dataset_name, align_mode, bed]
        }
    
    // Prepare transcriptome outputs for evaluation
    eval_transcriptome_inputs = FlairTranscriptome.out.transcriptome
        .map { test_name, dataset_name, align_mode, partition_mode, transcriptome_mode, 
               isoforms_bed, isoforms_gtf, isoforms_fa, isoform_counts, isoform_read_map ->
            def process_label = "transcriptome_${transcriptome_mode}"
            [test_name, dataset_name, align_mode, partition_mode, process_label, isoforms_bed]
        }
        .combine(align_for_eval, by: [0, 1, 2])  // Join by test_name, dataset_name, align_mode
        .combine(datasets_ch, by: 0)  // Join by test_name to get dataset info
        .map { test_name, dataset_name, align_mode, partition_mode, process_label, isoforms_bed, reads_bed,
               dataset, align_modes, partition_modes, correct_modes, collapse_modes, transcriptome_modes ->
            [test_name, dataset_name, align_mode, partition_mode, process_label, isoforms_bed, reads_bed, file(dataset.gtf)]
        }
    
    // Prepare collapse outputs for evaluation
    eval_collapse_inputs = FlairCollapse.out.collapsed
        .map { test_name, dataset_name, align_mode, partition_mode, correct_mode, collapse_mode,
               isoforms_bed, isoforms_fa, isoforms_gtf, isoform_read_map, isoform_counts ->
            def process_label = "collapse_${correct_mode}_${collapse_mode}"
            [test_name, dataset_name, align_mode, partition_mode, process_label, isoforms_bed]
        }
        .combine(align_for_eval, by: [0, 1, 2])  // Join by test_name, dataset_name, align_mode
        .combine(datasets_ch, by: 0)  // Join by test_name to get dataset info
        .map { test_name, dataset_name, align_mode, partition_mode, process_label, isoforms_bed, reads_bed,
               dataset, align_modes, partition_modes, correct_modes, collapse_modes, transcriptome_modes ->
            [test_name, dataset_name, align_mode, partition_mode, process_label, isoforms_bed, reads_bed, file(dataset.gtf)]
        }
    
    // Combine both evaluation inputs and run FlairEval
    eval_all_inputs = eval_transcriptome_inputs.mix(eval_collapse_inputs)
    FlairEval(eval_all_inputs)
    
    // =============================================================================
    // TED METRICS LOGIC
    // =============================================================================
    
    // Generate reference peaks from partitioned GTFs
    ref_peak_inputs = FlairPartition.out.partitioned
        .map { test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf ->
            [test_name, dataset_name, align_mode, partition_mode, gtf]
        }
    
    PrepareReferencePeaks(ref_peak_inputs)
    
    // Prepare align outputs for TED (get BAM files)
    align_for_ted = FlairAlign.out.alignments
        .map { test_name, dataset_name, align_mode, bam, bai, bed ->
            [test_name, dataset_name, align_mode, bam, bai]
        }
    
    // Get cage and quantseq files by dataset
    cage_quantseq_by_dataset = datasets_ch
        .map { test_name, dataset, align_modes, partition_modes, correct_modes, collapse_modes, transcriptome_modes ->
            def cage_file = dataset.cage ? file(dataset.cage) : file('NO_CAGE')
            def quantseq_file = dataset.quantseq ? file(dataset.quantseq) : file('NO_QUANTSEQ')
            [test_name, dataset.name, cage_file, quantseq_file]
        }
    
    // Prepare transcriptome outputs for TED
    ted_transcriptome_inputs = FlairTranscriptome.out.transcriptome
        .map { test_name, dataset_name, align_mode, partition_mode, transcriptome_mode,
               isoforms_bed, isoforms_gtf, isoforms_fa, isoform_counts, isoform_read_map ->
            def process_label = "transcriptome_${transcriptome_mode}"
            def stage = "transcriptome"
            def corrected_file = file('NO_CORRECTED')  // No corrected BED for transcriptome
            [test_name, dataset_name, align_mode, partition_mode, process_label, stage, 
             isoforms_bed, isoform_read_map, corrected_file]
        }
        .combine(align_for_ted, by: [0, 1, 2])  // Join by test_name, dataset_name, align_mode
        .combine(cage_quantseq_by_dataset, by: [0, 1])  // Join by test_name, dataset_name
        .combine(PrepareReferencePeaks.out.ref_peaks, by: [0, 1, 2, 3])  // Join by test_name, dataset_name, align_mode, partition_mode
        .map { test_name, dataset_name, align_mode, partition_mode, process_label, stage,
               isoforms_bed, isoform_read_map, corrected_file, bam, bai, cage_file, quantseq_file,
               ref_tss, ref_tts ->
            [test_name, dataset_name, align_mode, partition_mode, process_label, stage,
             isoforms_bed, isoform_read_map, bam, bai, corrected_file, cage_file, quantseq_file, ref_tss, ref_tts]
        }
    
    // Prepare collapse outputs for TED
    ted_collapse_inputs = FlairCollapse.out.collapsed
        .map { test_name, dataset_name, align_mode, partition_mode, correct_mode, collapse_mode,
               isoforms_bed, isoforms_fa, isoforms_gtf, isoform_read_map, isoform_counts ->
            def process_label = "collapse_${correct_mode}_${collapse_mode}"
            def stage = "collapse"
            [test_name, dataset_name, align_mode, partition_mode, correct_mode, process_label, stage,
             isoforms_bed, isoform_read_map]
        }
        .combine(align_for_ted, by: [0, 1, 2])  // Join by test_name, dataset_name, align_mode
        .combine(
            FlairCorrect.out.corrected.map { test_name, dataset_name, align_mode, partition_mode, correct_mode, 
                                              corrected_bed, inconsistent_bed ->
                [test_name, dataset_name, align_mode, partition_mode, correct_mode, corrected_bed]
            },
            by: [0, 1, 2, 3, 4]  // Join by test_name, dataset_name, align_mode, partition_mode, correct_mode
        )
        .combine(cage_quantseq_by_dataset, by: [0, 1])  // Join by test_name, dataset_name
        .combine(PrepareReferencePeaks.out.ref_peaks, by: [0, 1, 2, 3])  // Join by test_name, dataset_name, align_mode, partition_mode
        .map { test_name, dataset_name, align_mode, partition_mode, correct_mode, process_label, stage,
               isoforms_bed, isoform_read_map, bam, bai, corrected_bed, cage_file, quantseq_file,
               ref_tss, ref_tts ->
            [test_name, dataset_name, align_mode, partition_mode, process_label, stage,
             isoforms_bed, isoform_read_map, bam, bai, corrected_bed, cage_file, quantseq_file, ref_tss, ref_tts]
        }
    
    // Combine transcriptome and collapse TED inputs
    ted_all_inputs = ted_transcriptome_inputs.mix(ted_collapse_inputs)
    FlairTED(ted_all_inputs)
    
    // =============================================================================
    // SYNTHESIS AND SUMMARY LOGIC
    // =============================================================================
    
    // Group TED results by test name for synthesis
    ted_by_test = FlairTED.out.ted_metrics
        .map { test_name, dataset_name, align_mode, partition_mode, process_label, stage, ted_file ->
            [test_name, ted_file]
        }
        .groupTuple(by: 0)
    
    // Group FLAIR eval results by test name for synthesis
    flair_by_test = FlairEval.out.eval_results
        .map { test_name, dataset_name, align_mode, partition_mode, process_label, eval_file ->
            [test_name, eval_file]
        }
        .groupTuple(by: 0)
    
    // Combine both TED and FLAIR results for comprehensive synthesis
    all_results_input = ted_by_test.join(flair_by_test, by: 0)
    
    // Generate comprehensive summary with both TED and FLAIR evaluations
    SynthesizeEvaluations(all_results_input)
}