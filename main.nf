#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =============================================================================
// FLAIR Evaluation Pipeline — Modular Entrypoint
// =============================================================================
// Compares long-read RNA-seq isoform assemblers (FLAIR, Bambu, IsoQuant)
// using TED (Transcript End Distance) evaluation and concordance metrics.
//
// Usage:
//   nextflow run main.nf --input data/samplesheet.csv \
//       --params_file data/params.json --test_name my_test
//
// See REFACTORING_PLAN.md for architecture details.
// =============================================================================

// ---------------------------------------------------------------------------
// Module includes
// ---------------------------------------------------------------------------
include { FlairAlign           } from './modules/align/main'
include { FlairPartition       } from './modules/partition/main'
// include { PlotIsoforms         } from './modules/visualization/plot_isoforms/main'  // detached; re-enable when ready

// ---------------------------------------------------------------------------
// Subworkflow includes
// ---------------------------------------------------------------------------
include { ASSEMBLE_AND_EVAL    } from './subworkflows/assemble_and_eval'
include { SUMMARY_AND_VIZ      } from './subworkflows/summary_and_viz'

// ---------------------------------------------------------------------------
// Preflight validation
// ---------------------------------------------------------------------------
include { PreflightValidateFlair     } from './modules/preflight/main'
include { PreflightValidateIsoquant  } from './modules/preflight/main'

// =============================================================================
// WORKFLOW
// =============================================================================

workflow {
    // -------------------------------------------------------------------------
    // 1. Validate input parameters
    // -------------------------------------------------------------------------
    if (!params.input) {
        error "ERROR: Please provide a samplesheet via --input"
    }
    if (!params.params_file) {
        error "ERROR: Please provide a parameters configuration file via --params_file"
    }
    if (params.test_name == 'flair_test_suite') {
        log.warn "WARNING: Using default test name 'flair_test_suite'. Specify a custom name with --test_name"
    }

    // -------------------------------------------------------------------------
    // 2. Parse samplesheet + JSON config  (replaces ~80 lines of inline CSV parsing)
    // -------------------------------------------------------------------------
    def test_sets_list = Utils.parseSamplesheet(
        file(params.input), file(params.params_file), params.test_name
    )

    // Log alternative assembler status
    test_sets_list.each { test ->
        if (!test.bambuModes.isEmpty()) {
            log.info "Bambu modes enabled for ${test.name}: ${test.bambuModes.keySet().join(', ')}"
        }
        if (!test.isoquantModes.isEmpty()) {
            log.info "IsoQuant modes enabled for ${test.name}: ${test.isoquantModes.keySet().join(', ')}"
        }
        if (!test.isoseqModes.isEmpty()) {
            log.info "IsoSeq modes enabled for ${test.name}: ${test.isoseqModes.keySet().join(', ')}"
        }
        if (!test.flamesModes.isEmpty()) {
            log.info "FLAMES modes enabled for ${test.name}: ${test.flamesModes.keySet().join(', ')}"
        }
        if (!test.stringtie2Modes.isEmpty()) {
            log.info "StringTie2 modes enabled for ${test.name}: ${test.stringtie2Modes.keySet().join(', ')}"
        }
    }

    // Print summary
    println("\n=== FLAIR Test Suite ===")
    println("Input samplesheet: ${params.input}")
    println("Test name: ${params.test_name}")
    println("===================================")
    test_sets_list.each { test ->
        println("${test.name}: ${test.totalJobs()} jobs using dataset '${test.dataset.name}'")
        if (test.dataset.hasBam()) {
            println("  ⚠ WARNING: Pre-aligned BAM detected for '${test.dataset.name}' - FlairAlign will be SKIPPED")
            if (test.dataset.reads != null) {
                println("  ⚠ Note: 'reads' field will be ignored in favor of provided BAM files")
            }
        }
        if (!test.dataset.hasBam() && test.dataset.reads == null) {
            println("  🧬 SIMULATION: No BAM or reads — Badread will generate synthetic reads")
        }
    }
    println("Total jobs: ${test_sets_list.sum { it.totalJobs() }}")
    println("===================================\n")

    // =========================================================================
    // PREFLIGHT VALIDATION
    // Validate all assembler arg strings before submitting any real work.
    // Each condition's args are tested against the tool's own argument parser;
    // unknown flags cause an immediate hard failure with a clear error message.
    // =========================================================================

    // Collect every unique (condition_name, transcriptome_args) pair across all test sets
    def flair_preflight_entries = test_sets_list
        .collectMany { ts -> ts.transcriptomeModes.collect { mode, args -> [mode, args] } }
        .unique { it[0] }

    def isoquant_preflight_entries = test_sets_list
        .collectMany { ts -> ts.isoquantModes.collect { mode, args -> [mode, args] } }
        .unique { it[0] }

    // Pass all entries as a single list to one preflight job per assembler
    PreflightValidateFlair(Channel.of(flair_preflight_entries))
    if (!isoquant_preflight_entries.isEmpty()) {
        PreflightValidateIsoquant(Channel.of(isoquant_preflight_entries))
    }

    // Gate: each preflight process emits a single val(true) on success.
    // ASSEMBLE_AND_EVAL will not start until every preflight check has passed.
    def preflight_gate
    if (!isoquant_preflight_entries.isEmpty()) {
        preflight_gate = PreflightValidateFlair.out
            .concat(PreflightValidateIsoquant.out)
            .collect()
            .map { true }
    } else {
        preflight_gate = PreflightValidateFlair.out
    }

    // =========================================================================
    // CHANNEL CONSTRUCTION
    // =========================================================================

    // Create actual empty placeholder files on disk so Nextflow's HashBuilder
    // can resolve file attributes and staging doesn't fail with NoSuchFileException.
    ["NO_BED", "NO_CAGE", "NO_DRNA", "NO_ISOFORMS_BED", "NO_ISOFORMS_GTF", "NO_JUNCTION_TAB"].each {
        def f = new File("${workflow.workDir}/${it}")
        if (!f.exists()) { f.parentFile?.mkdirs(); f.createNewFile() }
    }
    def NO_BED          = file("${workflow.workDir}/NO_BED")
    def NO_CAGE         = file("${workflow.workDir}/NO_CAGE")
    def NO_DRNA     = file("${workflow.workDir}/NO_DRNA")
    def NO_ISOFORMS_BED = file("${workflow.workDir}/NO_ISOFORMS_BED")
    def NO_ISOFORMS_GTF = file("${workflow.workDir}/NO_ISOFORMS_GTF")
    def NO_JUNCTION_TAB = file("${workflow.workDir}/NO_JUNCTION_TAB")

    // Master datasets channel: [test_name, dataset, align_modes, partition_modes,
    //                           transcriptome_modes, bambu_modes, isoquant_modes,
    //                           isoseq_modes, flames_modes, stringtie2_modes]
    datasets_ch = Channel.from(test_sets_list)
        .map { test_set ->
            [test_set.name, test_set.dataset, test_set.alignModes, test_set.partitionModes,
             test_set.transcriptomeModes, test_set.bambuModes, test_set.isoquantModes,
             test_set.isoseqModes, test_set.flamesModes, test_set.stringtie2Modes]
        }

    // Branch: datasets with pre-aligned BAM (skip FlairAlign)
    datasets_with_bam = datasets_ch.filter {
        test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes ->
        dataset.hasBam()
    }

    // Branch: datasets needing alignment
    datasets_without_bam = datasets_ch.filter {
        test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes ->
        !dataset.hasBam() && dataset.reads != null
    }

    // Expand into per-read-file × per-align-mode jobs
    align_inputs = datasets_without_bam.flatMap {
        test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes ->
        ds_align_modes.collectMany { align_mode, align_args ->
            dataset.getReadsList().collect { reads_file ->
                [test_name, dataset.name, file(reads_file), align_mode, align_args, file(dataset.genome)]
            }
        }
    }

    // Pre-aligned BAM → partition (NO_BED placeholder; FlairPartition will generate BED)
    prealigned_partition_inputs = datasets_with_bam.flatMap {
        test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes ->
        def cage_file     = dataset.cage     ? file(dataset.cage)     : NO_CAGE
        def drna_file = dataset.drna ? file(dataset.drna) : NO_DRNA
        ds_align_modes.collectMany { align_mode, align_args ->
            ds_partition_modes.collect { partition_mode, partition_args ->
                [test_name, dataset.name, align_mode, file(dataset.bam), file(dataset.bai), NO_BED,
                 partition_mode, partition_args, file(dataset.genome), file(dataset.gtf),
                 cage_file, drna_file]
            }
        }
    }

    // =========================================================================
    // ALIGNMENT
    // =========================================================================
    FlairAlign(align_inputs)

    // =========================================================================
    // PARTITION
    // =========================================================================

    // Partition inputs from FlairAlign outputs (BED exists from alignment)
    partition_inputs_from_align = FlairAlign.out.alignments
        .combine(datasets_ch, by: 0)
        .flatMap {
            test_name, dataset_name, align_mode, bam, bai, bed,
            dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes ->
            def cage_file     = dataset.cage     ? file(dataset.cage)     : NO_CAGE
            def drna_file = dataset.drna ? file(dataset.drna) : NO_DRNA
            ds_partition_modes.collect { partition_mode, partition_args ->
                [test_name, dataset_name, align_mode, bam, bai, bed, partition_mode, partition_args,
                 file(dataset.genome), file(dataset.gtf), cage_file, drna_file]
            }
        }

    all_partition_inputs = partition_inputs_from_align
        .concat(prealigned_partition_inputs)
    FlairPartition(all_partition_inputs)

    // Canonical reference to all partitioned data
    all_partitioned = FlairPartition.out.partitioned

    // =========================================================================
    // ASSEMBLE, EVALUATE, AND SUMMARIZE (subworkflows)
    // =========================================================================

    // Per-dataset signal bedGraph paths + library type metadata
    dataset_signal_ch = datasets_ch.map {
        test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes ->
        [test_name, dataset.library_type ?: 'unknown',
         dataset.cage_signal_plus ?: '', dataset.cage_signal_minus ?: '',
         dataset.drna_signal_plus ?: '', dataset.drna_signal_minus ?: '']
    }

    // Bundle placeholder files for subworkflows
    def placeholders = [
        NO_ISOFORMS_BED: NO_ISOFORMS_BED,
        NO_ISOFORMS_GTF: NO_ISOFORMS_GTF,
        NO_JUNCTION_TAB: NO_JUNCTION_TAB,
        NO_CAGE:         NO_CAGE,
        NO_DRNA:     NO_DRNA,
    ]

    // --- Transcriptome assembly + evaluation ---
    // Gated on preflight: the combine() holds all_partitioned until the
    // preflight_gate value is emitted (i.e., all arg checks have passed).
    gated_partitioned = all_partitioned.combine(preflight_gate).map { items ->
        // Drop the trailing gate value; restore the original tuple
        items[0..-2]
    }
    ASSEMBLE_AND_EVAL(gated_partitioned, datasets_ch, dataset_signal_ch, placeholders)

    // --- Summary plots + divergence + UTR features ---
    SUMMARY_AND_VIZ(
        ASSEMBLE_AND_EVAL.out.evaluation_results,
        ASSEMBLE_AND_EVAL.out.cage_peak_reason_tsvs,
        ASSEMBLE_AND_EVAL.out.drna_peak_reason_tsvs,
        ASSEMBLE_AND_EVAL.out.all_eval_inputs,
        dataset_signal_ch,
        ASSEMBLE_AND_EVAL.out.ted_precision_metrics
    )

    // =========================================================================
    // VISUALIZATION — isoform structure plots (temporarily disabled)
    // =========================================================================
    // PlotIsoforms is detached for now; re-enable when ready.
    // plot_inputs = ASSEMBLE_AND_EVAL.out.flair_transcriptome ...
    // PlotIsoforms(plot_inputs)
}
