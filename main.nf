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
include { ASSEMBLE             } from './subworkflows/assemble'
include { EVALUATE             } from './subworkflows/evaluate'
include { SUMMARY_AND_VIZ      } from './subworkflows/summary_and_viz'

// ---------------------------------------------------------------------------
// Preflight validation
// ---------------------------------------------------------------------------
include { PreflightValidateFlair     } from './modules/preflight/main'
include { PreflightValidateIsoquant  } from './modules/preflight/main'

include { WriteAssemblyManifest      } from './modules/manifest/main'

// =============================================================================
// WORKFLOW
// =============================================================================

// -----------------------------------------------------------------------------
// STAGE 1 — ASSEMBLY
// -----------------------------------------------------------------------------
// preflight -> align -> partition -> ALL assemblers. STOPS after assemblers.
// Publishes per-tool outputs (via the assembler modules' own publishDir) and
// writes results/assemblies/<test_name>/assembly_manifest.csv.
//
// Emits normalized channels so the DEFAULT (monolith) workflow can drive
// EVALUATE in-process without a manifest round-trip.
// -----------------------------------------------------------------------------
workflow ASSEMBLY {

    main:
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
        if (!test.lraaModes.isEmpty()) {
            log.info "LRAA modes enabled for ${test.name}: ${test.lraaModes.keySet().join(', ')}"
        }
    }

    // -------------------------------------------------------------------------
    // Counts-filter advisory
    // -------------------------------------------------------------------------
    // ted.py --counts filters isoforms_observed to transcripts with read
    // support >= 1.  Assemblers without a native per-transcript count file
    // (currently: StringTie2) skip the filter and report every emitted
    // transcript.  Log this once at startup so it's transparent in the run
    // log alongside the per-test mode summaries above.
    def assemblers_without_counts = ['stringtie2']
    log.info "Assemblers WITHOUT native --counts support (filter skipped): ${assemblers_without_counts.join(', ')}"

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
    //                           isoseq_modes, flames_modes, stringtie2_modes,
    //                           lraa_modes]
    datasets_ch = Channel.from(test_sets_list)
        .map { test_set ->
            [test_set.name, test_set.dataset, test_set.alignModes, test_set.partitionModes,
             test_set.transcriptomeModes, test_set.bambuModes, test_set.isoquantModes,
             test_set.isoseqModes, test_set.flamesModes, test_set.stringtie2Modes,
             test_set.lraaModes]
        }

    // Branch: datasets with pre-aligned BAM (skip FlairAlign)
    datasets_with_bam = datasets_ch.filter {
        test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes, ds_lraa_modes ->
        dataset.hasBam()
    }

    // Branch: datasets needing alignment
    datasets_without_bam = datasets_ch.filter {
        test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes, ds_lraa_modes ->
        !dataset.hasBam() && dataset.reads != null
    }

    // Expand into per-read-file × per-align-mode jobs
    align_inputs = datasets_without_bam.flatMap {
        test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes, ds_lraa_modes ->
        ds_align_modes.collectMany { align_mode, align_args ->
            dataset.getReadsList().collect { reads_file ->
                [test_name, dataset.name, file(reads_file), align_mode, align_args, file(dataset.genome)]
            }
        }
    }

    // Pre-aligned BAM → partition (NO_BED placeholder; FlairPartition will generate BED)
    prealigned_partition_inputs = datasets_with_bam.flatMap {
        test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes, ds_lraa_modes ->
        def cage_file     = dataset.cage     ? file(dataset.cage)     : NO_CAGE
        def drna_file = dataset.drna ? file(dataset.drna) : NO_DRNA
        // genome_orig: original canonical genome path threaded through partition
        // output as a val (never staged/copied). See partition/main.nf output note.
        def genome_orig = file(dataset.genome).toString()
        ds_align_modes.collectMany { align_mode, align_args ->
            ds_partition_modes.collect { partition_mode, partition_args ->
                [test_name, dataset.name, align_mode, file(dataset.bam), file(dataset.bai), NO_BED,
                 partition_mode, partition_args, file(dataset.genome), genome_orig, file(dataset.gtf),
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
            dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes, ds_lraa_modes ->
            def cage_file     = dataset.cage     ? file(dataset.cage)     : NO_CAGE
            def drna_file = dataset.drna ? file(dataset.drna) : NO_DRNA
            // genome passed BOTH as a staged path (partition reads it, symlink,
            // no copy) AND as its original path string (genome_orig val) that we
            // thread through the partition output so downstream assemblers read
            // the canonical genome directly without ever staging/copying 2.8 GB.
            def genome_orig = file(dataset.genome).toString()
            ds_partition_modes.collect { partition_mode, partition_args ->
                [test_name, dataset_name, align_mode, bam, bai, bed, partition_mode, partition_args,
                 file(dataset.genome), genome_orig, file(dataset.gtf), cage_file, drna_file]
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
        test_name, dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes, ds_lraa_modes ->
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

    // --- Transcriptome assembly (STAGE 1 stops here) ---
    // Gated on preflight: the combine() holds all_partitioned until the
    // preflight_gate value is emitted (i.e., all arg checks have passed).
    gated_partitioned = all_partitioned.combine(preflight_gate).map { items ->
        // Drop the trailing gate value; restore the original tuple
        items[0..-2]
    }
    ASSEMBLE(gated_partitioned, datasets_ch, dataset_signal_ch, placeholders)

    // --- Write the assembly manifest (one row per produced assembly) ---
    manifest_by_test = ASSEMBLE.out.manifest_rows
        .map { row -> [row.split(',', 2)[0], row] }   // key by test_name (col 0)
        .groupTuple()
    WriteAssemblyManifest(manifest_by_test)

    emit:
        isoform_tuples    = ASSEMBLE.out.isoform_tuples
        eval_context_seed = ASSEMBLE.out.eval_context_seed
        flair_firstpass   = ASSEMBLE.out.flair_firstpass
        flair_cds_bed     = ASSEMBLE.out.flair_cds_bed
        dataset_signal    = dataset_signal_ch
        manifest          = WriteAssemblyManifest.out.manifest
}

// -----------------------------------------------------------------------------
// STAGE 2 — EVALUATION
// -----------------------------------------------------------------------------
// Takes the STAGE-1 manifest CSV as --input. For each row, LOADS the assembly +
// orthogonal data and runs PrepareReferencePeaks + Evaluation + TedEndPrecision
// + FirstpassComparison + ALL of SUMMARY_AND_VIZ. No assembler/align/partition.
// Cheap, frequently-re-run stage for threshold selection + logreg + plots.
// -----------------------------------------------------------------------------
workflow EVALUATION {

    main:
    if (!params.input) {
        error "ERROR: EVALUATION requires the STAGE-1 manifest CSV via --input"
    }
    if (!params.params_file) {
        error "ERROR: Please provide a parameters configuration file via --params_file (needed for partition-region diversity plots)"
    }

    def manifest_rows = Utils.parseManifest(file(params.input))
    if (manifest_rows.isEmpty()) {
        error "ERROR: manifest '${params.input}' contained no assembly rows"
    }
    println("\n=== FLAIR EVALUATION (STAGE 2) ===")
    println("Manifest      : ${params.input}")
    println("Assembly rows : ${manifest_rows.size()}")
    println("Test name     : ${params.test_name}")
    println("===================================\n")

    // Placeholder files (mirror STAGE-1 channel construction)
    ["NO_CAGE", "NO_DRNA", "NO_ISOFORMS_BED", "NO_ISOFORMS_GTF",
     "NO_TED_LOG", "NO_COUNTS", "NO_FIRSTPASS", "NO_CDS"].each {
        def f = new File("${workflow.workDir}/${it}")
        if (!f.exists()) { f.parentFile?.mkdirs(); f.createNewFile() }
    }
    def NO_CAGE         = file("${workflow.workDir}/NO_CAGE")
    def NO_DRNA         = file("${workflow.workDir}/NO_DRNA")
    def NO_ISOFORMS_BED = file("${workflow.workDir}/NO_ISOFORMS_BED")
    def NO_ISOFORMS_GTF = file("${workflow.workDir}/NO_ISOFORMS_GTF")
    def NO_TED_LOG      = file("${workflow.workDir}/NO_TED_LOG")
    def NO_COUNTS       = file("${workflow.workDir}/NO_COUNTS")
    def NO_FIRSTPASS    = file("${workflow.workDir}/NO_FIRSTPASS")
    def NO_CDS          = file("${workflow.workDir}/NO_CDS")
    def placeholders = [
        NO_ISOFORMS_BED: NO_ISOFORMS_BED, NO_ISOFORMS_GTF: NO_ISOFORMS_GTF,
        NO_CAGE: NO_CAGE, NO_DRNA: NO_DRNA,
    ]

    // Helper: resolve a manifest path field to a file() or a placeholder.
    def resolveFile = { String val, placeholder ->
        (val == null || val.trim().isEmpty()) ? placeholder : file(val)
    }

    // ---- Reconstruct the normalized channels from manifest rows ----
    def iso_list = []
    def ctx_list = []
    def fp_list  = []
    def cds_list = []
    def sig_seen = [:]
    def sig_list = []

    manifest_rows.each { r ->
        def key4 = [r.test_name, r.dataset_name, r.align_mode, r.partition_mode]
        // isoform tuple (10-field)
        iso_list << (key4 + [
            r.transcriptome_mode,
            resolveFile(r.isoforms_bed, NO_ISOFORMS_BED),
            resolveFile(r.isoforms_gtf, NO_ISOFORMS_GTF),
            resolveFile(r.isoform_read_map, NO_COUNTS),
            resolveFile(r.ted_log, NO_TED_LOG),
            resolveFile(r.isoform_counts, NO_COUNTS),
        ])
        // eval_context_seed (12-field): keys + partition outputs + partition_args
        // genome carried as a val string (canonical absolute path).
        ctx_list << (key4 + [
            resolveFile(r.bam, NO_CAGE),
            resolveFile(r.bai, NO_CAGE),
            resolveFile(r.reads_bed, NO_CAGE),
            (r.genome ?: '').toString(),
            resolveFile(r.gtf, NO_CAGE),
            resolveFile(r.cage_peaks, NO_CAGE),
            resolveFile(r.drna_peaks, NO_DRNA),
            (r.partition_args ?: '').toString(),
        ])
        // flair firstpass (7-field) — only for FLAIR rows with a real firstpass
        if (r.tool == 'flair' && r.firstpass && r.firstpass.trim() &&
            !r.firstpass.endsWith('NO_FIRSTPASS')) {
            fp_list << (key4 + [r.partition_args ?: '', r.transcriptome_mode,
                                resolveFile(r.firstpass, NO_FIRSTPASS)])
        }
        // flair cds (8-field) — only for FLAIR rows with a real CDS bed
        if (r.tool == 'flair' && r.cds_bed && r.cds_bed.trim() &&
            !r.cds_bed.endsWith('NO_CDS')) {
            cds_list << (key4 + [r.partition_args ?: '', r.transcriptome_mode,
                                 resolveFile(r.cds_bed, NO_CDS),
                                 resolveFile(r.cds_info, NO_CDS)])
        }
        // dataset signal (one per test_name) — first row wins
        if (!sig_seen.containsKey(r.test_name)) {
            sig_seen[r.test_name] = true
            sig_list << [r.test_name, r.library_type ?: 'unknown',
                         r.cage_signal_plus ?: '', r.cage_signal_minus ?: '',
                         r.drna_signal_plus ?: '', r.drna_signal_minus ?: '']
        }
    }

    isoform_tuples    = Channel.fromList(iso_list)
    eval_context_seed = Channel.fromList(ctx_list)
    flair_firstpass   = Channel.fromList(fp_list)
    flair_cds_bed     = Channel.fromList(cds_list)
    dataset_signal_ch = Channel.fromList(sig_list)

    // ---- Core evaluation ----
    EVALUATE(isoform_tuples, eval_context_seed, flair_firstpass,
             dataset_signal_ch, placeholders)

    // ---- Summary + viz (unchanged surface) ----
    SUMMARY_AND_VIZ(
        EVALUATE.out.evaluation_results,
        EVALUATE.out.isoform_categories,
        EVALUATE.out.cage_peak_reason_tsvs,
        EVALUATE.out.drna_peak_reason_tsvs,
        EVALUATE.out.all_eval_inputs,
        dataset_signal_ch,
        EVALUATE.out.ted_precision_metrics,
        EVALUATE.out.ted_gtf_precision,
        flair_firstpass,
        flair_cds_bed
    )
}

// -----------------------------------------------------------------------------
// DEFAULT workflow — back-compat monolith.
// Runs STAGE 1 (ASSEMBLY) then drives EVALUATE + SUMMARY_AND_VIZ IN-PROCESS
// from ASSEMBLE.out (no manifest round-trip needed). The manifest is still
// written as a side effect so a later EVALUATION run can reuse the assemblies.
// -----------------------------------------------------------------------------
workflow {
    ASSEMBLY()

    EVALUATE(
        ASSEMBLY.out.isoform_tuples,
        ASSEMBLY.out.eval_context_seed,
        ASSEMBLY.out.flair_firstpass,
        ASSEMBLY.out.dataset_signal,
        placeholders_value()
    )

    SUMMARY_AND_VIZ(
        EVALUATE.out.evaluation_results,
        EVALUATE.out.isoform_categories,
        EVALUATE.out.cage_peak_reason_tsvs,
        EVALUATE.out.drna_peak_reason_tsvs,
        EVALUATE.out.all_eval_inputs,
        ASSEMBLY.out.dataset_signal,
        EVALUATE.out.ted_precision_metrics,
        EVALUATE.out.ted_gtf_precision,
        ASSEMBLY.out.flair_firstpass,
        ASSEMBLY.out.flair_cds_bed
    )
}

// Shared placeholder map for the default in-process EVALUATE driver. The
// files are created by ASSEMBLY's channel-construction block; here we only
// need the map of file() handles for EVALUATE's NO_CAGE/NO_DRNA fallbacks.
def placeholders_value() {
    ["NO_ISOFORMS_BED", "NO_ISOFORMS_GTF", "NO_CAGE", "NO_DRNA"].each {
        def f = new File("${workflow.workDir}/${it}")
        if (!f.exists()) { f.parentFile?.mkdirs(); f.createNewFile() }
    }
    return [
        NO_ISOFORMS_BED: file("${workflow.workDir}/NO_ISOFORMS_BED"),
        NO_ISOFORMS_GTF: file("${workflow.workDir}/NO_ISOFORMS_GTF"),
        NO_CAGE:         file("${workflow.workDir}/NO_CAGE"),
        NO_DRNA:         file("${workflow.workDir}/NO_DRNA"),
    ]
}

// =============================================================================
// EMAIL NOTIFICATIONS
// Fires once the workflow finishes (success OR failure). Uses the system
// `sendmail` binary when no `mail {}` SMTP block is configured.
// See: https://docs.seqera.io/nextflow/notifications
// =============================================================================
workflow.onComplete {
    def status     = workflow.success ? 'SUCCESS' : 'FAILED'
    def emoji      = workflow.success ? '✅' : '❌'
    def recipient  = 'hdheath@ucsc.edu'

    def body = """\
        ${emoji} FLAIR Evaluation Pipeline — ${status}

        Test name      : ${params.test_name}
        Started        : ${workflow.start}
        Completed      : ${workflow.complete}
        Duration       : ${workflow.duration}
        Exit status    : ${workflow.exitStatus}
        Command line   : ${workflow.commandLine}
        Work dir       : ${workflow.workDir}
        Output dir     : ${params.outdir}
        Project dir    : ${workflow.projectDir}
        Launch dir     : ${workflow.launchDir}
        Nextflow ver   : ${workflow.nextflow.version}
        Run name       : ${workflow.runName}
        Session ID     : ${workflow.sessionId}

        ${workflow.success ? '' : 'Error message : ' + (workflow.errorMessage ?: 'n/a')}
        ${workflow.success ? '' : 'Error report  : ' + (workflow.errorReport  ?: 'n/a')}
        """.stripIndent()

    // Attach the HTML execution report if it was generated (it usually is,
    // but very-early failures may abort before `report{}` writes the file).
    def report_html = file("${params.outdir}/reports/execution_report.html")
    def attachments = report_html.exists() ? [report_html] : []

    try {
        sendMail(
            to:      recipient,
            from:    'hdheath@ucsc.edu',
            subject: "[Nextflow] flair4_ted ${status}: ${params.test_name} (${workflow.runName})",
            body:    body,
            attach:  attachments
        )
    } catch (Exception e) {
        log.warn "Could not send completion email to ${recipient}: ${e.message}"
    }
}

workflow.onError {
    log.error "Pipeline execution stopped with message: ${workflow.errorMessage}"
}
