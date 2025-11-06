#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =============================================================================
// PARAMETERS - Define configurable options
// =============================================================================

params {
    // Input data parameters - for single dataset mode
    input_reads = null               // Path to reads file (FASTA/FASTQ)
    genome = null                    // Path to reference genome (FASTA)  
    gtf = null                       // Path to GTF annotation file
    
    // FLAIR ALIGN parameters (based on flair align -h)
    // Only set defaults for required parameters, let users decide everything else
    threads = null                   // minimap2 number of threads, null = use FLAIR default (4)
    junction_bed = null              // annotated isoforms/junctions bed file
    nvrna = null                     // use native-RNA alignment parameters, null = use FLAIR default (false)
    quality = null                   // minimum MAPQ of read alignment (0-60), null = use FLAIR default (0)
    minfragmentsize = null           // minimum size of alignment kept
    maxintronlen = null              // maximum intron length in genomic alignment
    filtertype = null                // 'removesup', 'separate', 'keysup', null = use FLAIR default ('removesup')
    remove_internal_priming = null   // remove reads with internal priming, null = use FLAIR default (false)
    intprimingthreshold = null       // bases that are at least 75% As for internal priming
    intprimingfracAs = null          // fraction As required for internal priming
    remove_singleexon = null         // remove unspliced reads, null = use FLAIR default (false)
    quiet = null                     // suppress minimap progress statements, null = use FLAIR default (false)
    
    // ALIGN MODE PRESETS - minimal presets, only override FLAIR defaults when necessary
    align_presets = [
        'default': [
            // Use all FLAIR defaults - no overrides
        ],
        'nvrna': [
            nvrna: true  // Only override what's needed for native RNA
        ],
        'stringent': [
            quality: 20,           // Higher quality threshold
            minfragmentsize: 200   // Longer reads only
        ],
        'sensitive': [
            quality: 5,              // Lower quality threshold  
            maxintronlen: 2000000    // Allow very long introns
        ]
    ]
    
    // FLAIR TRANSCRIPTOME parameters (based on flair transcriptome -h)
    // Only set defaults for required parameters, let users decide everything else
    junction_tab = null              // short-read junctions in SJ.out.tab format
    junction_support = null          // minimum junction support required
    ss_window = null                 // window size for correcting splice sites, null = use FLAIR default (15)
    end_window = null                // window size for comparing TSS/TES, null = use FLAIR default (100)
    sjc_support = null               // minimum reads for spliced isoform, null = use FLAIR default (1)
    se_support = null                // minimum reads for single exon isoform, null = use FLAIR default (3)
    frac_support = null              // minimum fraction of gene locus support, null = use FLAIR default (0.05)
    no_stringent = null              // don't require full-length reads, null = use FLAIR default (false)
    no_check_splice = null           // don't enforce accurate splice site alignment, null = use FLAIR default (false)
    no_align_to_annot = null         // don't align to annotation first, null = use FLAIR default (false)
    no_redundant = null              // 'none', 'longest', 'best_only', null = use FLAIR default ('none')
    max_ends = null                  // maximum TSS/TES per isoform, null = use FLAIR default (1)
    filter = null                    // 'nosubset', 'bysupport', 'comprehensive', 'ginormous'
    predict_cds = null               // predict CDS of final isoforms, null = use FLAIR default (false)
    keep_intermediate = null         // keep intermediate files, null = use FLAIR default (false)
    keep_sup = null                  // keep supplementary alignments, null = use FLAIR default (false)
    
    // TRANSCRIPTOME MODE PRESETS - minimal presets, only override FLAIR defaults when necessary
    transcriptome_presets = [
        'with_gtf': [
            // Use all FLAIR defaults - no overrides
        ],
        'denovo': [
            // Use all FLAIR defaults - no overrides
        ],
        'stringent': [
            sjc_support: 3,      // Higher support requirements
            se_support: 5,
            frac_support: 0.1,   // Require 10% of gene locus support
            end_window: 50       // Tighter end clustering
        ],
        'permissive': [
            sjc_support: 1,        // Lower support requirements
            se_support: 1,
            frac_support: 0.01,    // Only require 1% of gene locus support
            end_window: 200,       // Looser end clustering
            max_ends: 3,           // Allow more TSS/TES variants
            no_stringent: true     // Don't require full-length reads
        ]
    ]
    
    // Preset modes for easy configuration
    align_mode = 'default'           // 'default', 'nvrna', 'stringent', 'sensitive'
    transcriptome_mode = 'with_gtf'  // 'with_gtf', 'denovo', 'stringent', 'permissive'
    
    // Output parameters
    outdir = 'results'
    
    // Test parameters - for when running predefined test sets
    run_test_sets = true
    test_name = 'quick_test'  // 'quick_test', 'align_comparison', 'transcriptome_comparison'
}

// =============================================================================
// CLASS DEFINITIONS - Define dataset and test configuration objects
// =============================================================================

class Dataset {
    String name
    def reads  // Can be String (single file) or List (multiple files)
    String genome
    String gtf
    
    // Optional files for FLAIR processes
    String junction_bed
    String junction_tab
    
    // Constructor
    Dataset(String name, Map config) {
        this.name = name
        this.reads = config.reads
        this.genome = config.genome
        this.gtf = config.gtf
        this.junction_bed = config.junction_bed
        this.junction_tab = config.junction_tab
    }
    
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
    List<List> alignOptions
    List<List> transcriptomeOptions
    
    // Constructor
    TestSet(String name, Dataset dataset, Map options) {
        this.name = name
        this.dataset = dataset
        this.alignOptions = options.align
        this.transcriptomeOptions = options.transcriptome ?: []
    }
    
    // Calculate total number of jobs this test set will produce
    int totalJobs() {
        return alignOptions.size() * transcriptomeOptions.size()
    }
    
    String description() {
        return "TestSet[${name}, ${totalJobs()} jobs]"
    }
}

// =============================================================================
// PROCESS DEFINITIONS
// =============================================================================

process FlairAlign {
    publishDir "${params.outdir}/align", mode: 'symlink'
    errorStrategy 'ignore'  // Allow other test sets to continue if this fails
    conda 'flair-dev'
    
    input:
    tuple val(test_name), val(dataset_name), path(reads), val(align_mode), path(junction_bed, stageAs: 'junction.bed')
    path genome
    
    output:
    tuple val(test_name), val(dataset_name), val(align_mode), 
          path("${dataset_name}_${align_mode}.bam"), 
          path("${dataset_name}_${align_mode}.bam.bai"),
          path("${dataset_name}_${align_mode}.bed"), emit: alignments
    
    script:
    // Get mode-specific defaults, then override with user parameters
    def mode_defaults = params.align_presets[align_mode]
    if (mode_defaults == null) {
        error "Unknown align mode '${align_mode}'. Available modes: ${params.align_presets.keySet()}"
    }
    
    // Resolve final parameter values: user params override mode defaults
    def final_quality = params.quality != null ? params.quality : mode_defaults.quality
    def final_minfragmentsize = params.minfragmentsize != null ? params.minfragmentsize : mode_defaults.minfragmentsize  
    def final_maxintronlen = params.maxintronlen != null ? params.maxintronlen : mode_defaults.maxintronlen
    def final_nvrna = params.nvrna || (mode_defaults.nvrna ?: false)
    
    // Build alignment arguments based on resolved parameters
    def align_args = []
    
    // Only add parameters that are explicitly set (not null)
    if (params.threads != null) align_args << "--threads ${params.threads}"
    if (final_quality != null) align_args << "--quality ${final_quality}"
    if (final_minfragmentsize != null) align_args << "--minfragmentsize ${final_minfragmentsize}"
    if (final_maxintronlen != null) align_args << "--maxintronlen ${final_maxintronlen}"
    if (final_nvrna) align_args << "--nvrna"
    
    // Add other optional parameters only if set
    if (junction_bed.name != 'NO_JUNCTION_BED') align_args << "--junction_bed ${junction_bed}"
    if (params.filtertype != null) align_args << "--filtertype ${params.filtertype}"
    if (params.quiet != null && params.quiet) align_args << "--quiet"
    if (params.remove_singleexon != null && params.remove_singleexon) align_args << "--remove_singleexon"
    
    // Handle internal priming removal
    if (params.remove_internal_priming) {
        align_args << "--remove_internal_priming"
        if (params.gtf) align_args << "-f ${params.gtf}"
        if (params.intprimingthreshold) align_args << "--intprimingthreshold ${params.intprimingthreshold}"
        if (params.intprimingfracAs) align_args << "--intprimingfracAs ${params.intprimingfracAs}"
    }
    
    def args_string = align_args.join(' ')
    
    """
    flair align ${args_string} \\
        -r ${reads} \\
        -g ${genome} \\
        -o ${dataset_name}_${align_mode}
    """
}

process FlairTranscriptome {
    publishDir "${params.outdir}/transcriptome/${test_name}", mode: 'symlink'
    errorStrategy 'ignore'  // Allow other alignments to continue if transcriptome fails
    conda 'flair-dev'
    
    input:
    tuple val(test_name), val(dataset_name), val(align_mode), path(bam), path(bai), path(bed), val(transcriptome_mode), path(junction_tab, stageAs: 'junction.tab'), path(junction_bed, stageAs: 'junction.bed')
    path genome
    path gtf
    
    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(transcriptome_mode),
          path("${dataset_name}_${align_mode}_${transcriptome_mode}.transcriptome.isoforms.bed"),
          path("${dataset_name}_${align_mode}_${transcriptome_mode}.transcriptome.isoforms.gtf"),
          path("${dataset_name}_${align_mode}_${transcriptome_mode}.transcriptome.isoforms.fa"),
          path("${dataset_name}_${align_mode}_${transcriptome_mode}.transcriptome.isoform.counts.txt"),
          path("${dataset_name}_${align_mode}_${transcriptome_mode}.transcriptome.isoform.read.map.txt"), emit: transcriptome
    
    script:
    // Get mode-specific defaults, then override with user parameters
    def mode_defaults = params.transcriptome_presets[transcriptome_mode]
    if (mode_defaults == null) {
        error "Unknown transcriptome mode '${transcriptome_mode}'. Available modes: ${params.transcriptome_presets.keySet()}"
    }
    
    // Resolve final parameter values: user params override mode defaults
    def final_ss_window = params.ss_window != null ? params.ss_window : mode_defaults.ss_window
    def final_end_window = params.end_window != null ? params.end_window : mode_defaults.end_window
    def final_sjc_support = params.sjc_support != null ? params.sjc_support : mode_defaults.sjc_support
    def final_se_support = params.se_support != null ? params.se_support : mode_defaults.se_support
    def final_frac_support = params.frac_support != null ? params.frac_support : mode_defaults.frac_support
    def final_max_ends = params.max_ends != null ? params.max_ends : mode_defaults.max_ends
    def final_no_stringent = params.no_stringent != null ? params.no_stringent : mode_defaults.no_stringent
    
    // Build transcriptome arguments based on resolved parameters
    def transcriptome_args = []
    
    // Only add threads if explicitly set
    if (params.threads != null) transcriptome_args << "--threads ${params.threads}"
    
    // Handle GTF annotation
    if (transcriptome_mode == 'with_gtf' && gtf) {
        transcriptome_args << "-f ${gtf}"
    }
    
    // Add junction support only if provided  
    if (junction_tab.name != 'NO_JUNCTION_TAB') transcriptome_args << "--junction_tab ${junction_tab}"
    if (junction_bed.name != 'NO_JUNCTION_BED') transcriptome_args << "--junction_bed ${junction_bed}"
    if (params.junction_support) transcriptome_args << "--junction_support ${params.junction_support}"
    
    // Add resolved window and support parameters only if set
    if (final_ss_window != null) transcriptome_args << "--ss_window ${final_ss_window}"
    if (final_end_window != null) transcriptome_args << "--end_window ${final_end_window}"
    if (final_sjc_support != null) transcriptome_args << "--sjc_support ${final_sjc_support}"
    if (final_se_support != null) transcriptome_args << "--se_support ${final_se_support}"
    if (final_frac_support != null) transcriptome_args << "--frac_support ${final_frac_support}"
    if (final_max_ends != null) transcriptome_args << "--max_ends ${final_max_ends}"
    
    // Add boolean flags only if explicitly set to true
    if (final_no_stringent != null && final_no_stringent) transcriptome_args << "--no_stringent"
    if (params.no_check_splice != null && params.no_check_splice) transcriptome_args << "--no_check_splice"
    if (params.no_align_to_annot != null && params.no_align_to_annot) transcriptome_args << "--no_align_to_annot"
    if (params.predict_cds != null && params.predict_cds) transcriptome_args << "--predict_cds"
    if (params.keep_intermediate != null && params.keep_intermediate) transcriptome_args << "--keep_intermediate"
    if (params.keep_sup != null && params.keep_sup) transcriptome_args << "--keep_sup"
    
    // Add filter and redundancy options only if set
    if (params.no_redundant != null && params.no_redundant != 'none') transcriptome_args << "--no_redundant ${params.no_redundant}"
    if (params.filter) transcriptome_args << "--filter ${params.filter}"
    
    def args_string = transcriptome_args.join(' ')
    
    """
    flair transcriptome \\
        -b ${bam} \\
        -g ${genome} \\
        ${args_string} \\
        -o ${dataset_name}_${align_mode}_${transcriptome_mode}.transcriptome
    """
}

workflow {
    // =============================================================================
    // DEFINE DATASETS AS OBJECTS
    // =============================================================================
    
    def human_dataset = new Dataset('human', [
        reads: '/private/groups/brookslab/hdheath/projects/flair-eval/cache/datasets/human/WTC11.100reads.fasta',
        genome: '/private/groups/brookslab/reference_genomes/gencode_human/GRCh38.primary_assembly.genome.fa',
        gtf: '/private/groups/brookslab/reference_genomes/gencode_human/gencode.v47.annotation.gtf',
        junction_bed: '/private/groups/brookslab/hdheath/projects/flair-eval/cache/datasets/human/junctions.bed',
        junction_tab: null  // No STAR junctions for this dataset
    ])
    
    def a549_dataset = new Dataset('A549_cDNA', [
        reads: '/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/WTC11.100reads.fasta',
        genome: '/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/GRCh38.primary_assembly.genome.fa',
        gtf: '/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/gencode.v48.annotation.gtf',
        junction_bed: null,  // No junction bed for this dataset
        junction_tab: '/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/WTC11_all.SJ.out.tab'
    ])
    
    // =============================================================================
    // DEFINE STANDARD/GLOBAL OPTIONS
    // =============================================================================
    
    def standard_align_modes = ['default', 'nvrna', 'stringent', 'sensitive']
    def standard_transcriptome_modes = ['with_gtf', 'denovo', 'stringent', 'permissive']
    
    // =============================================================================
    // DEFINE TEST SETS
    // =============================================================================
    
    // Test 1: Quick validation on A549 dataset
    def quick_test = new TestSet('quick_test', a549_dataset, [
        align: [['default']],
        transcriptome: [['with_gtf']]
    ])
    // 1 * 1 = 1 job
    
    // Test 2: Align comparison on human dataset  
    def align_comparison_test = new TestSet('align_comparison', human_dataset, [
        align: standard_align_modes.collect { [it] },
        transcriptome: [['with_gtf']]
    ])
    // 4 * 1 = 4 jobs (default, nvrna, stringent, sensitive)
    
    // Test 3: Transcriptome comparison
    def transcriptome_comparison_test = new TestSet('transcriptome_comparison', a549_dataset, [
        align: [['default']],
        transcriptome: standard_transcriptome_modes.collect { [it] }
    ])
    // 1 * 4 = 4 jobs (with_gtf, denovo, stringent, permissive)
    
    // Test 4: Full matrix test (commented out by default to avoid too many jobs)
    // def full_matrix_test = new TestSet('full_matrix', human_dataset, [
    //     align: standard_align_modes.collect { [it] },
    //     transcriptome: standard_transcriptome_modes.collect { [it] }
    // ])
    // 4 * 4 = 16 jobs (all combinations!)
    
    // =============================================================================
    // COLLECT ALL TEST SETS TO RUN
    // =============================================================================
    
    def test_sets = [
        quick_test,
        // align_comparison_test,
        // transcriptome_comparison_test,
        // full_matrix_test
    ]
    
    // Print summary
    println("\n=== FLAIR Align + Transcriptome Test Suite ===")
    test_sets.each { test ->
        println("${test.name}: ${test.totalJobs()} jobs using dataset '${test.dataset.name}'")
    }
    println("Total jobs: ${test_sets.sum { it.totalJobs() }}")
    println("==============================================\n")
    
    // =============================================================================
    // BUILD CHANNELS FROM TEST SETS
    // =============================================================================
    
    // Create align inputs from all test sets
    align_inputs = Channel.from(test_sets)
        .flatMap { test_set ->
            // Get reads as list (handles both single and multiple files)
            def reads_list = test_set.dataset.getReadsList()
            
            // Create one alignment job per read file × align mode
            reads_list.collectMany { reads_file ->
                test_set.alignOptions.collect { align_option ->
                    // Handle optional junction bed file
                    def junction_bed_file = test_set.dataset.junction_bed ? 
                        file(test_set.dataset.junction_bed) : file('NO_JUNCTION_BED')
                    
                    tuple(
                        test_set.name,              // test_name
                        test_set.dataset.name,      // dataset_name
                        file(reads_file),           // reads (single file)
                        align_option[0],            // align_mode
                        junction_bed_file           // junction_bed (or placeholder)
                    )
                }
            }
        }
    
    // Get unique genomes (in case multiple datasets use same genome)
    unique_genome = Channel.from(test_sets)
        .map { it.dataset.genome }
        .unique()
        .map { file(it) }
        .first()
    
    // Run alignment
    FlairAlign(align_inputs, unique_genome)
    
    // Create transcriptome inputs from alignment outputs
    transcriptome_inputs = FlairAlign.out.alignments
        .flatMap { test_name, dataset_name, align_mode, bam, bai, bed ->
            // Find the test set that produced this alignment
            def test_set = test_sets.find { it.name == test_name }
            
            // If no transcriptome options, skip this alignment
            if (!test_set.transcriptomeOptions || test_set.transcriptomeOptions.isEmpty()) {
                return []
            }
            
            test_set.transcriptomeOptions.collect { transcriptome_option ->
                // Handle optional junction files
                def junction_tab_file = test_set.dataset.junction_tab ? 
                    file(test_set.dataset.junction_tab) : file('NO_JUNCTION_TAB')
                def junction_bed_file = test_set.dataset.junction_bed ? 
                    file(test_set.dataset.junction_bed) : file('NO_JUNCTION_BED')
                
                tuple(
                    test_name,
                    dataset_name,
                    align_mode,
                    bam,
                    bai,
                    bed,
                    transcriptome_option[0],    // transcriptome_mode
                    junction_tab_file,          // junction_tab (or placeholder)
                    junction_bed_file           // junction_bed (or placeholder)
                )
            }
        }
    
    // Get unique GTF files
    unique_gtf = Channel.from(test_sets)
        .map { it.dataset.gtf }
        .unique()
        .map { file(it) }
        .first()
    
    // Run transcriptome generation
    FlairTranscriptome(transcriptome_inputs, unique_genome, unique_gtf)
    
    // =============================================================================
    // OUTPUT SUMMARY
    // =============================================================================
    
    FlairAlign.out.alignments.view { test_name, dataset_name, align_mode, bam, bai, bed ->
        "ALIGN COMPLETE: ${test_name}/${dataset_name}/${align_mode} -> ${bed.name}"
    }
    
    FlairTranscriptome.out.transcriptome.view { test_name, dataset_name, align_mode, transcriptome_mode, 
                                                 isoforms_bed, isoforms_gtf, isoforms_fa, isoform_counts, isoform_read_map ->
        "TRANSCRIPTOME COMPLETE: ${test_name}/${dataset_name}/${align_mode}/${transcriptome_mode} -> ${isoforms_bed.name}"
    }
}