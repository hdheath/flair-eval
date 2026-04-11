// Module: CrossSampleSummary
// Cross-sample summary processes that combine results from ALL samples
// into unified publication-quality figures.

process CrossSamplePrecisionRecall {
    publishDir "${params.outdir}/summary/${params.test_name}/precision_recall", mode: 'copy'
    tag "cross_sample_precision_recall"

    input:
    tuple val(test_name), path(evaluation_files)

    output:
    path "precision_recall/*.{png,svg}", emit: cross_sample_pr_plot, optional: true

    script:
    """
    # v2: pair-aware dedup precision (2026-04-02)
    python ${projectDir}/bin/evaluation/precision_recall_plot.py \\
        --input ${evaluation_files} \\
        --output precision_recall \\
        --verbose
    """
}

process CrossSampleConcordance {
    publishDir "${params.outdir}/summary/${params.test_name}/concordance", mode: 'copy'
    tag "cross_sample_concordance"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(evaluation_files)

    output:
    path "concordance/*.{png,svg}", emit: cross_sample_concordance_plot, optional: true

    script:
    """
    python ${projectDir}/bin/evaluation/concordance_plots.py \\
        --input ${evaluation_files} \\
        --output concordance \\
        --verbose || true
    """
}

process CrossSampleSignalSupport {
    publishDir "${params.outdir}/summary/${params.test_name}/signal_support", mode: 'copy'
    tag "cross_sample_signal_support"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(evaluation_files)

    output:
    path "signal_support/*.{png,svg}", emit: cross_sample_signal_support_plot, optional: true

    script:
    """
    python ${projectDir}/bin/evaluation/signal_support_summary_plot.py \\
        --input ${evaluation_files} \\
        --output signal_support \\
        --verbose || true
    """
}

process CrossSampleLandscape {
    publishDir "${params.outdir}/summary/${params.test_name}/landscape", mode: 'copy'
    tag "cross_sample_landscape"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(evaluation_files)

    output:
    path "landscape/*.{png,svg}", emit: cross_sample_landscape_plot, optional: true

    script:
    """
    python ${projectDir}/bin/evaluation/transcriptome_landscape_plot.py \\
        --input ${evaluation_files} \\
        --output landscape \\
        --verbose || true
    """
}

process CrossSampleToolEndAccuracy {
    publishDir "${params.outdir}/summary/${params.test_name}/tool_end_accuracy", mode: 'copy'
    tag "cross_sample_tool_end_accuracy"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(evaluation_files)

    output:
    path "tool_end_accuracy/*.{png,svg}", emit: tool_end_accuracy_plot, optional: true

    script:
    """
    # v2: use transcriptome_mode labels directly (2026-04-02)
    python ${projectDir}/bin/evaluation/tool_end_accuracy_boxplot.py \\
        --input ${evaluation_files} \\
        --output tool_end_accuracy \\
        --verbose
    """
}

// -----------------------------------------------------------------
// Cross-sample signal-stratified peak recovery curves.
// Aggregates peak_reason TSVs across all datasets; produces mean ± SD
// curves and per-dataset curves for CAGE and dRNA signal types.
// -----------------------------------------------------------------
process CrossSamplePeakRoc {
    publishDir "${params.outdir}/summary/${params.test_name}/peak_roc", mode: 'copy'
    tag "cross_sample_peak_roc"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(reason_tsvs)

    output:
    path "signal_recall_*.{png,svg}", emit: cross_sample_peak_roc_plot, optional: true

    script:
    """
    python ${projectDir}/bin/evaluation/peak_roc_curves.py \\
        --input ${reason_tsvs} \\
        --output . \\
        --dataset "${test_name}" \\
        --verbose || true
    """
}

