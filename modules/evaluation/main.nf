// Module: Evaluation
// Unified evaluation process for all assemblers (FLAIR, Bambu, IsoQuant).
// Automatically detects assembler type and runs appropriate evaluation:
//   - FLAIR: Full evaluation with --isoforms-bed and read-level metrics
//   - GTF/GFF-based assemblers: evaluation with --gtf-input; read metrics only when a real read map exists
//
// Outputs:
//   1. TED (Transcript End Distance) — measures TSS/TTS accuracy
//   2. FLAIR eval — measures structural accuracy against reference annotation
//   3. Synthesized TSV combining all metrics

process Evaluation {
    // --- Per-method plots: category folder, method name in filename ---
    // motif_logos: tss/tts_motif_logo.png
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/per_method/motif_logos", mode: 'copy', pattern: 'ted_plots/*motif_logo*.png', saveAs: { it.toString().tokenize('/').last() }
    // peak_recovery: recovery_by_expression/width, peak_read_support
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/per_method/peak_recovery", mode: 'copy', pattern: 'ted_plots/*recovery_by_*.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/per_method/peak_recovery", mode: 'copy', pattern: 'ted_plots/*peak_read_support.png', saveAs: { it.toString().tokenize('/').last() }
    // distance_histograms: cage/drna distance + read offset histograms
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/per_method/distance_histograms", mode: 'copy', pattern: 'ted_plots/*distance_histogram.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/per_method/distance_histograms", mode: 'copy', pattern: 'ted_plots/*offset_histogram.png', saveAs: { it.toString().tokenize('/').last() }
    // entropy: TSS/TTS entropy distributions
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/per_method/entropy", mode: 'copy', pattern: 'ted_plots/*entropy*.png', saveAs: { it.toString().tokenize('/').last() }
    // --- Logs ---
    publishDir "${params.outdir}/logs/${test_name}", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_${it}" }
    tag "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(transcriptome_mode),
          path(isoforms_bed), path(isoforms_gtf), path(isoform_read_map),
          path(ted_log),
          path(isoform_counts),
          path(bam), path(bai), path(reads_bed), path(genome), path(gtf),
          path(cage_peaks), path(drna_peaks),
          path(ref_tss), path(ref_tts),
          val(library_type),
          val(cage_signal_plus), val(cage_signal_minus),
          val(drna_signal_plus), val(drna_signal_minus)

    output:
    // Core evaluation result
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(transcriptome_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_evaluation.tsv"), emit: evaluation_results
    // Per-isoform SQANTI category labels (consumed by SqantiPrecision to avoid re-classifying)
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(transcriptome_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_isoform_categories.tsv"), emit: isoform_categories
    // All per-run plots (published but not consumed downstream)
    path "ted_plots/*.png", optional: true, emit: all_plots
    // Per-peak reason TSVs (consumed downstream by PeakReasonHeatmap)
    tuple val(test_name), path("ted_plots/*_cage_peak_reasons.tsv"), optional: true, emit: cage_peak_reason_tsvs
    tuple val(test_name), path("ted_plots/*_drna_peak_reasons.tsv"), optional: true, emit: drna_peak_reason_tsvs
    // Test region outputs (published only)
    path "test_regions/*", optional: true, emit: test_regions
    // Performance timing report
    path "*_ted_timing.txt", optional: true, emit: timing_reports

    script:
    // v2: pair-aware dedup precision (2026-04-02)
    // Determine evaluation path: GTF-based (simplified) vs BED-based (FLAIR).
    // Any assembler that provides a GTF/GFF instead of isoforms BED goes through
    // the simplified path (Bambu, IsoQuant, IsoSeq, FLAMES, StringTie2, etc.)
    def is_simplified = (isoforms_bed.name == 'NO_ISOFORMS_BED')

    // Build optional arguments
    def cage_arg = (cage_peaks.name != 'NO_CAGE' && cage_peaks.size() > 0) ? "--prime5-peaks ${cage_peaks}" : ""
    def drna_arg = (drna_peaks.name != 'NO_DRNA' && drna_peaks.size() > 0) ? "--prime3-peaks ${drna_peaks}" : ""
    def ref_tss_arg = ref_tss.name != 'NO_REF_TSS' ? "--ref-prime5-peaks ${ref_tss}" : ""
    def ref_tts_arg = ref_tts.name != 'NO_REF_TTS' ? "--ref-prime3-peaks ${ref_tts}" : ""
    // Signal arguments for ted.py
    def cage_signal_plus_arg_ted = cage_signal_plus ? "--cage-signal-plus ${cage_signal_plus}" : ""
    def cage_signal_minus_arg_ted = cage_signal_minus ? "--cage-signal-minus ${cage_signal_minus}" : ""
    def drna_signal_plus_arg_ted = drna_signal_plus ? "--drna-signal-plus ${drna_signal_plus}" : ""
    def drna_signal_minus_arg_ted = drna_signal_minus ? "--drna-signal-minus ${drna_signal_minus}" : ""
    def library_type_arg = library_type && library_type != 'unknown' ? "--library-type ${library_type}" : ""
    def supports_read_metrics = !transcriptome_mode.toLowerCase().startsWith('stringtie2')
    // TED internal decision log (only exists for FLAIR --ted runs)
    def ted_log_arg = (ted_log.name != 'NO_TED_LOG' && ted_log.size() > 0) ? "--ted-log ${ted_log}" : ""
    // Per-transcript counts file (FLAIR / IsoQuant / Bambu emit one; others pass NO_COUNTS).
    // When present, ted.py filters isoforms_observed to transcripts with count >= 1,
    // preventing zero-support reference carry-throughs from inflating the count
    // (notably IsoQuant: ~30% of GTF transcripts have zero supporting reads).
    def counts_arg = (isoform_counts.name != 'NO_COUNTS' && isoform_counts.size() > 0) ? "--counts ${isoform_counts}" : ""
    def output_prefix = "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}"

    if (is_simplified)
    """
    # v2: pair-aware dedup precision (2026-04-02)
    mkdir -p ted_plots
    mkdir -p test_regions

    READ_MAP_ARG=""
    BAM_ARG=""
    READS_BED_ARG=""
    SKIP_ARG=""
    if [ "${supports_read_metrics}" = "true" ] && [ -s "${isoform_read_map}" ]; then
        READ_MAP_ARG="--map-file ${isoform_read_map}"
        BAM_ARG="--bam ${bam}"
        READS_BED_ARG="--reads-bed ${reads_bed}"
    else
        SKIP_ARG="--skip-read-metrics"
    fi

    python ${projectDir}/bin/ted.py \\
        --gtf-input ${isoforms_gtf} \\
        \$READ_MAP_ARG \\
        \$BAM_ARG \\
        \$READS_BED_ARG \\
        \$SKIP_ARG \\
        --genome ${genome} \\
        --gtf ${gtf} \\
        ${cage_arg} \\
        ${drna_arg} \\
        ${ref_tss_arg} \\
        ${ref_tts_arg} \\
        ${cage_signal_plus_arg_ted} \\
        ${cage_signal_minus_arg_ted} \\
        ${drna_signal_plus_arg_ted} \\
        ${drna_signal_minus_arg_ted} \\
        --window 50 \\
        --window-5prime 50 \\
        --window-3prime 50 \\
        --stage transcriptome \\
        --test-name ${test_name} \\
        --dataset-name ${dataset_name} \\
        --align-mode ${align_mode} \\
        --partition-mode ${partition_mode} \\
        --pipeline-mode ${transcriptome_mode} \\
        ${library_type_arg} \\
        --plot-output-dir ted_plots \\
        --plot-mode both \\
        --plot-dpi 300 \\
        --test-regions-dir test_regions \\
        --timing-output ${output_prefix}_ted_timing.txt \\
        --output ${output_prefix}_ted.tsv \\
        ${ted_log_arg} \\
        ${counts_arg} \\
        --verbose

    python ${projectDir}/bin/flair_eval.py \\
        --reads-bed ${reads_bed} \\
        --gtf-input ${isoforms_gtf} \\
        --gtf ${gtf} \\
        --test-name ${test_name} \\
        --dataset-name ${dataset_name} \\
        --align-mode ${align_mode} \\
        --partition-mode ${partition_mode} \\
        --pipeline-mode ${transcriptome_mode} \\
        ${library_type_arg} \\
        --stage transcriptome \\
        --plot-output-dir ted_plots \\
        --plot-prefix ${output_prefix} \\
        --categories-output ${output_prefix}_isoform_categories.tsv \\
        --output ${output_prefix}_flair_eval.tsv \\
        ${counts_arg} \\
        --verbose

    python ${projectDir}/bin/synthesize_evaluations.py \\
        --ted-files ${output_prefix}_ted.tsv \\
        --flair-files ${output_prefix}_flair_eval.tsv \\
        --output ${output_prefix}_evaluation.tsv \\
        --test-name ${test_name}

    rm -f ${output_prefix}_ted.tsv ${output_prefix}_flair_eval.tsv
    """

    else
    """
    # v2: pair-aware dedup precision (2026-04-02)
    mkdir -p ted_plots
    mkdir -p test_regions

    python ${projectDir}/bin/ted.py \\
        --isoforms-bed ${isoforms_bed} \\
        --map-file ${isoform_read_map} \\
        --bam ${bam} \\
        --reads-bed ${reads_bed} \\
        --genome ${genome} \\
        --gtf ${gtf} \\
        ${cage_arg} \\
        ${drna_arg} \\
        ${ref_tss_arg} \\
        ${ref_tts_arg} \\
        ${cage_signal_plus_arg_ted} \\
        ${cage_signal_minus_arg_ted} \\
        ${drna_signal_plus_arg_ted} \\
        ${drna_signal_minus_arg_ted} \\
        --window 50 \\
        --window-5prime 50 \\
        --window-3prime 50 \\
        --stage transcriptome \\
        --test-name ${test_name} \\
        --dataset-name ${dataset_name} \\
        --align-mode ${align_mode} \\
        --partition-mode ${partition_mode} \\
        --pipeline-mode ${transcriptome_mode} \\
        ${library_type_arg} \\
        --plot-output-dir ted_plots \\
        --plot-mode both \\
        --plot-dpi 300 \\
        --test-regions-dir test_regions \\
        --timing-output ${output_prefix}_ted_timing.txt \\
        --output ${output_prefix}_ted.tsv \\
        ${ted_log_arg} \\
        ${counts_arg} \\
        --verbose

    python ${projectDir}/bin/flair_eval.py \\
        --reads-bed ${reads_bed} \\
        --isoforms-bed ${isoforms_bed} \\
        --gtf ${gtf} \\
        --test-name ${test_name} \\
        --dataset-name ${dataset_name} \\
        --align-mode ${align_mode} \\
        --partition-mode ${partition_mode} \\
        --pipeline-mode ${transcriptome_mode} \\
        ${library_type_arg} \\
        --stage transcriptome \\
        --plot-output-dir ted_plots \\
        --plot-prefix ${output_prefix} \\
        --categories-output ${output_prefix}_isoform_categories.tsv \\
        --output ${output_prefix}_flair_eval.tsv \\
        ${counts_arg} \\
        --verbose

    python ${projectDir}/bin/synthesize_evaluations.py \\
        --ted-files ${output_prefix}_ted.tsv \\
        --flair-files ${output_prefix}_flair_eval.tsv \\
        --output ${output_prefix}_evaluation.tsv \\
        --test-name ${test_name}

    rm -f ${output_prefix}_ted.tsv ${output_prefix}_flair_eval.tsv
    """
}


/*
 * SqantiPrecisionRecall
 * ---------------------
 * Per-SQANTI-category precision/recall and end-redundancy analysis.
 *
/*
 * TedEndPrecision
 * ---------------
 * Junction-chain-deduplicated end precision for TED transcriptome modes.
 *
 * For each junction chain group, counts unique orthogonal-data peaks matched
 * by isoform ends (same peak = 1 TP, not N).  This detects over-segmentation
 * where multiple endvars redundantly map to the same CAGE or dRNA peak.
 *
 * Only runs when at least one orthogonal peak file is provided.
 *
 * Outputs:
 *   precision_recall_summary.tsv — per-mode dedup precision, recall, F1
 *   per_junction_chain.tsv       — per-JC detail
 */
process TedEndPrecision {
    tag "${dataset_name}_${transcriptome_mode}"

    input:
    // isoforms_bed is NO_ISOFORMS_BED for GTF-only assemblers (Bambu/IsoQuant/etc.)
    // isoforms_gtf is NO_ISOFORMS_GTF for FLAIR (BED-only assembler)
    tuple val(test_name), val(dataset_name), val(transcriptome_mode),
          path(isoforms_bed), path(isoforms_gtf), path(annotation_gtf),
          path(isoform_counts),
          path(cage_peaks), path(drna_peaks),
          val(partition_args)

    output:
    tuple val(test_name), val(dataset_name), val(transcriptome_mode),
          path("${transcriptome_mode}_precision_recall_summary.tsv"),
          path("${transcriptome_mode}_per_junction_chain.tsv"), emit: metrics
    // GTF-only precision (no orthogonal peaks) — used for reference-vs-orthogonal scatter
    tuple val(test_name), val(dataset_name), val(transcriptome_mode),
          path("${transcriptome_mode}_gtf_precision_recall_summary.tsv"), emit: gtf_metrics

    script:
    def cage_arg = (cage_peaks.name != 'NO_CAGE' && cage_peaks.size() > 0) ? "--peaks-5prime ${cage_peaks}" : ""
    def drna_arg = (drna_peaks.name != 'NO_DRNA' && drna_peaks.size() > 0) ? "--peaks-3prime ${drna_peaks}" : ""
    // partition_args is e.g. "--region chr22:16000000-26000000 chr3:48000000-53000000" or empty.
    def region_match = (partition_args =~ /--region\s+(.+?)(?:\s+--|$)/)
    def region_arg = region_match ? "--region ${region_match[0][1].trim()}" : ""
    def isoforms_arg = (isoforms_bed.name != 'NO_ISOFORMS_BED') ? "--isoforms-bed ${isoforms_bed}" : "--isoforms-gtf ${isoforms_gtf}"
    def counts_arg = (isoform_counts.name != 'NO_COUNTS' && isoform_counts.size() > 0) ? "--counts ${isoform_counts} --min-support 1" : ""
    """
    # v3: GTF-based recall now uses every distinct annotated TSS/TTS as the
    # denominator (was JC-filtered, which made annotation-passthrough tools
    # score ~100% trivially).
    # --- Orthogonal-signal precision (peaks, primary metric) ---
    python ${projectDir}/bin/evaluation/ted_end_precision.py \\
        ${isoforms_arg} \\
        --gtf ${annotation_gtf} \\
        ${cage_arg} \\
        ${drna_arg} \\
        ${region_arg} \\
        ${counts_arg} \\
        --mode ${transcriptome_mode} \\
        --window 50 \\
        --outdir .

    mv precision_recall_summary.tsv ${transcriptome_mode}_precision_recall_summary.tsv
    mv per_junction_chain.tsv ${transcriptome_mode}_per_junction_chain.tsv

    # --- Reference-GTF precision (no peaks) — annotation concordance ---
    python ${projectDir}/bin/evaluation/ted_end_precision.py \\
        ${isoforms_arg} \\
        --gtf ${annotation_gtf} \\
        ${region_arg} \\
        ${counts_arg} \\
        --mode ${transcriptome_mode} \\
        --window 50 \\
        --outdir gtf_only

    mv gtf_only/precision_recall_summary.tsv ${transcriptome_mode}_gtf_precision_recall_summary.tsv
    """
}


/*
 * FirstpassComparison
 * -------------------
 * Compares firstpass (pre-TED) and final (post-TED) isoform BEDs using
 * junction-chain-deduplicated end precision/recall.  Shows the effect of
 * TED filtering on isoform end accuracy.
 *
 * Outputs:
 *   firstpass_vs_final.tsv — side-by-side P/R/F1 for both stages
 *   firstpass_vs_final.png — grouped bar chart comparison
 */
process FirstpassComparison {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/per_method/firstpass_comparison", mode: 'copy', saveAs: { "${transcriptome_mode}_${it}" }
    tag "${dataset_name}_${transcriptome_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(transcriptome_mode),
          path(firstpass_bed), path(final_bed), path(annotation_gtf),
          path(cage_peaks), path(drna_peaks),
          val(partition_args)

    output:
    tuple val(test_name), val(dataset_name), val(transcriptome_mode),
          path("firstpass_vs_final.tsv"), path("firstpass_vs_final.png"), emit: comparison

    script:
    def cage_arg = (cage_peaks.name != 'NO_CAGE' && cage_peaks.size() > 0) ? "--peaks-5prime ${cage_peaks}" : ""
    def drna_arg = (drna_peaks.name != 'NO_DRNA' && drna_peaks.size() > 0) ? "--peaks-3prime ${drna_peaks}" : ""
    def region_match = (partition_args =~ /--region\s+(.+?)(?:\s+--|$)/)
    def region_arg = region_match ? "--region ${region_match[0][1].trim()}" : ""
    """
    python ${projectDir}/bin/evaluation/firstpass_vs_final.py \\
        --firstpass-bed ${firstpass_bed} \\
        --final-bed ${final_bed} \\
        --gtf ${annotation_gtf} \\
        ${cage_arg} \\
        ${drna_arg} \\
        ${region_arg} \\
        --mode ${transcriptome_mode} \\
        --window 50 \\
        --outdir .
    """
}
