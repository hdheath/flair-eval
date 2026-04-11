// Module: Evaluation
// Unified evaluation process for all assemblers (FLAIR, Bambu, IsoQuant).
// Automatically detects assembler type and runs appropriate evaluation:
//   - FLAIR: Full evaluation with --isoforms-bed and read-level metrics
//   - Bambu/IsoQuant: Full evaluation with --gtf-input and read-level metrics from converted read maps
//
// Outputs:
//   1. TED (Transcript End Distance) — measures TSS/TTS accuracy
//   2. FLAIR eval — measures structural accuracy against reference annotation
//   3. Synthesized TSV combining all metrics

process Evaluation {
    // --- Data files ---
    publishDir "${params.outdir}/evaluations/${test_name}/data", mode: 'symlink', pattern: '*_evaluation.tsv'
    publishDir "${params.outdir}/evaluations/${test_name}/data", mode: 'copy', pattern: 'ted_plots/*_peak_reasons.tsv', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/data/test_regions", mode: 'copy', pattern: 'test_regions/*.bed'
    publishDir "${params.outdir}/evaluations/${test_name}/data/test_regions", mode: 'copy', pattern: 'test_regions/*.csv'
    publishDir "${params.outdir}/evaluations/${test_name}/data/test_regions", mode: 'copy', pattern: 'test_regions/*.tsv'
    publishDir "${params.outdir}/evaluations/${test_name}/data", mode: 'copy', pattern: '*_ted_timing.txt'
    // --- Per-method plots routed into category subdirectories ---
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/signal_profiles", mode: 'copy', pattern: 'ted_plots/*signal_meta_profile.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/boundary_signal", mode: 'copy', pattern: 'ted_plots/*boundary_signal_*.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/motif_logos", mode: 'copy', pattern: 'ted_plots/*motif_logo*.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/peak_recovery", mode: 'copy', pattern: 'ted_plots/*recovery_by_*.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/peak_recovery", mode: 'copy', pattern: 'ted_plots/*peak_width_histogram.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/peak_recovery", mode: 'copy', pattern: 'ted_plots/*peak_recovery_dashboard.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/peak_recovery", mode: 'copy', pattern: 'ted_plots/*peak_read_support.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/peak_diagnostics", mode: 'copy', pattern: 'ted_plots/*missed_peak_*.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/peak_diagnostics", mode: 'copy', pattern: 'ted_plots/*peak_signal_context_dashboard.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/peak_diagnostics", mode: 'copy', pattern: 'ted_plots/*missed_sj_support.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/peak_diagnostics", mode: 'copy', pattern: 'ted_plots/*missed_truncation_patterns.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/distance_histograms", mode: 'copy', pattern: 'ted_plots/*distance_histogram.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/distance_histograms", mode: 'copy', pattern: 'ted_plots/*offset_histogram.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/read_analysis", mode: 'copy', pattern: 'ted_plots/*read_classification.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/read_analysis", mode: 'copy', pattern: 'ted_plots/*read_frequency*.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/read_analysis", mode: 'copy', pattern: 'ted_plots/*_all_truncation_patterns.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/read_analysis", mode: 'copy', pattern: 'ted_plots/*signal_vs_support.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/entropy", mode: 'copy', pattern: 'ted_plots/*entropy*.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/end_structure", mode: 'copy', pattern: 'ted_plots/*end_distance_dashboard.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/end_structure", mode: 'copy', pattern: 'ted_plots/*end_structure_dashboard.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/end_structure", mode: 'copy', pattern: 'ted_plots/*junction_chain_end_variation*.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/internal_priming", mode: 'copy', pattern: 'ted_plots/*internal_priming*.png', saveAs: { it.toString().tokenize('/').last() }
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/internal_priming", mode: 'copy', pattern: 'ted_plots/*proximal_apa*.png', saveAs: { it.toString().tokenize('/').last() }
    // --- Logs ---
    publishDir "${params.outdir}/logs/${test_name}", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_${it}" }
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(transcriptome_mode),
          path(isoforms_bed), path(isoforms_gtf), path(isoform_read_map),
          path(ted_log),
          path(bam), path(bai), path(reads_bed), path(genome), path(gtf),
          path(cage_peaks), path(quantseq_peaks),
          path(ref_tss), path(ref_tts),
          val(library_type),
          val(cage_signal_plus), val(cage_signal_minus),
          val(quantseq_signal_plus), val(quantseq_signal_minus)

    output:
    // Core evaluation result
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(transcriptome_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_evaluation.tsv"), emit: evaluation_results
    // All per-run plots (published but not consumed downstream)
    path "ted_plots/*.png", optional: true, emit: all_plots
    // Per-peak reason TSVs (consumed downstream by PeakReasonHeatmap)
    tuple val(test_name), path("ted_plots/*_cage_peak_reasons.tsv"), optional: true, emit: cage_peak_reason_tsvs
    tuple val(test_name), path("ted_plots/*_quantseq_peak_reasons.tsv"), optional: true, emit: quantseq_peak_reason_tsvs
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
    def quantseq_arg = (quantseq_peaks.name != 'NO_QUANTSEQ' && quantseq_peaks.size() > 0) ? "--prime3-peaks ${quantseq_peaks}" : ""
    def ref_tss_arg = ref_tss.name != 'NO_REF_TSS' ? "--ref-prime5-peaks ${ref_tss}" : ""
    def ref_tts_arg = ref_tts.name != 'NO_REF_TTS' ? "--ref-prime3-peaks ${ref_tts}" : ""
    // Signal arguments for ted.py
    def cage_signal_plus_arg_ted = cage_signal_plus ? "--cage-signal-plus ${cage_signal_plus}" : ""
    def cage_signal_minus_arg_ted = cage_signal_minus ? "--cage-signal-minus ${cage_signal_minus}" : ""
    def quantseq_signal_plus_arg_ted = quantseq_signal_plus ? "--quantseq-signal-plus ${quantseq_signal_plus}" : ""
    def quantseq_signal_minus_arg_ted = quantseq_signal_minus ? "--quantseq-signal-minus ${quantseq_signal_minus}" : ""
    def library_type_arg = library_type && library_type != 'unknown' ? "--library-type ${library_type}" : ""
    // TED internal decision log (only exists for FLAIR --ted runs)
    def ted_log_arg = (ted_log.name != 'NO_TED_LOG' && ted_log.size() > 0) ? "--ted-log ${ted_log}" : ""
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
    if [ -s "${isoform_read_map}" ]; then
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
        ${quantseq_arg} \\
        ${ref_tss_arg} \\
        ${ref_tts_arg} \\
        ${cage_signal_plus_arg_ted} \\
        ${cage_signal_minus_arg_ted} \\
        ${quantseq_signal_plus_arg_ted} \\
        ${quantseq_signal_minus_arg_ted} \\
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
        --output ${output_prefix}_flair_eval.tsv \\
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
        ${quantseq_arg} \\
        ${ref_tss_arg} \\
        ${ref_tts_arg} \\
        ${cage_signal_plus_arg_ted} \\
        ${cage_signal_minus_arg_ted} \\
        ${quantseq_signal_plus_arg_ted} \\
        ${quantseq_signal_minus_arg_ted} \\
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
        --output ${output_prefix}_flair_eval.tsv \\
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
 * where multiple endvars redundantly map to the same CAGE or QuantSeq peak.
 *
 * Only runs when at least one orthogonal peak file is provided.
 *
 * Outputs:
 *   precision_recall_summary.tsv — per-mode dedup precision, recall, F1
 *   per_junction_chain.tsv       — per-JC detail
 */
process TedEndPrecision {
    publishDir "${params.outdir}/evaluations/${test_name}/data/ted_precision/${transcriptome_mode}", mode: 'copy', pattern: '*.tsv'
    errorStrategy 'ignore'
    tag "${dataset_name}_${transcriptome_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(transcriptome_mode),
          path(isoforms_bed), path(annotation_gtf),
          path(cage_peaks), path(quantseq_peaks),
          val(partition_args)

    output:
    tuple val(test_name), val(dataset_name), val(transcriptome_mode),
          path("precision_recall_summary.tsv"), path("per_junction_chain.tsv"), emit: metrics

    script:
    def cage_arg = (cage_peaks.name != 'NO_CAGE' && cage_peaks.size() > 0) ? "--peaks-5prime ${cage_peaks}" : ""
    def quantseq_arg = (quantseq_peaks.name != 'NO_QUANTSEQ' && quantseq_peaks.size() > 0) ? "--peaks-3prime ${quantseq_peaks}" : ""
    // partition_args is e.g. "--region chr22:16000000-26000000" or empty; extract region value
    def region_match = (partition_args =~ /--region\s+(\S+)/)
    def region_arg = region_match ? "--region ${region_match[0][1]}" : ""
    """
    python ${projectDir}/bin/evaluation/ted_end_precision.py \\
        --isoforms-bed ${isoforms_bed} \\
        --gtf ${annotation_gtf} \\
        ${cage_arg} \\
        ${quantseq_arg} \\
        ${region_arg} \\
        --mode ${transcriptome_mode} \\
        --window 50 \\
        --outdir .
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
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/firstpass_comparison/${transcriptome_mode}", mode: 'copy'
    errorStrategy 'ignore'
    tag "${dataset_name}_${transcriptome_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(transcriptome_mode),
          path(firstpass_bed), path(final_bed), path(annotation_gtf),
          path(cage_peaks), path(quantseq_peaks),
          val(partition_args)

    output:
    tuple val(test_name), val(dataset_name), val(transcriptome_mode),
          path("firstpass_vs_final.tsv"), path("firstpass_vs_final.png"), emit: comparison

    script:
    def cage_arg = (cage_peaks.name != 'NO_CAGE' && cage_peaks.size() > 0) ? "--peaks-5prime ${cage_peaks}" : ""
    def quantseq_arg = (quantseq_peaks.name != 'NO_QUANTSEQ' && quantseq_peaks.size() > 0) ? "--peaks-3prime ${quantseq_peaks}" : ""
    def region_match = (partition_args =~ /--region\s+(\S+)/)
    def region_arg = region_match ? "--region ${region_match[0][1]}" : ""
    """
    python ${projectDir}/bin/evaluation/firstpass_vs_final.py \\
        --firstpass-bed ${firstpass_bed} \\
        --final-bed ${final_bed} \\
        --gtf ${annotation_gtf} \\
        ${cage_arg} \\
        ${quantseq_arg} \\
        ${region_arg} \\
        --mode ${transcriptome_mode} \\
        --window 50 \\
        --outdir .
    """
}
