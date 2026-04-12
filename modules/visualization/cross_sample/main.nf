// Module: CrossSampleSummary
// Cross-sample summary processes that combine results from ALL samples
// into unified publication-quality figures.

process CrossSamplePrecisionRecall {
    publishDir "${params.outdir}/summary/${params.test_name}/end_accuracy/precision_recall", mode: 'copy'
    tag "cross_sample_precision_recall"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(precision_recall_tsvs)

    output:
    path "*.png", emit: cross_sample_pr_plot, optional: true

    script:
    """
    python ${projectDir}/bin/evaluation/precision_recall_plot.py \\
        --input ${precision_recall_tsvs} \\
        --output . \\
        --verbose
    """
}

process CrossSampleConcordance {
    publishDir "${params.outdir}/summary/${params.test_name}/comparison/concordance", mode: 'copy'
    tag "cross_sample_concordance"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(evaluation_files)

    output:
    path "*.png", emit: cross_sample_concordance_plot, optional: true

    script:
    """
    python ${projectDir}/bin/evaluation/concordance_plots.py \\
        --input ${evaluation_files} \\
        --output . \\
        --verbose || true
    """
}

process CrossSampleSignalSupport {
    publishDir "${params.outdir}/summary/${params.test_name}/comparison/signal_support", mode: 'copy'
    tag "cross_sample_signal_support"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(evaluation_files)

    output:
    path "*.png", emit: cross_sample_signal_support_plot, optional: true

    script:
    """
    python ${projectDir}/bin/evaluation/signal_support_summary_plot.py \\
        --input ${evaluation_files} \\
        --output . \\
        --verbose || true
    """
}

process CrossSampleLandscape {
    publishDir "${params.outdir}/summary/${params.test_name}/comparison/landscape", mode: 'copy'
    tag "cross_sample_landscape"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(evaluation_files)

    output:
    path "*.png", emit: cross_sample_landscape_plot, optional: true

    script:
    """
    python ${projectDir}/bin/evaluation/transcriptome_landscape_plot.py \\
        --input ${evaluation_files} \\
        --output . \\
        --verbose || true
    """
}

process CrossSampleToolEndAccuracy {
    publishDir "${params.outdir}/summary/${params.test_name}/end_accuracy/tool_end_accuracy", mode: 'copy'
    tag "cross_sample_tool_end_accuracy"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(evaluation_files)

    output:
    path "*.png", emit: tool_end_accuracy_plot, optional: true

    script:
    """
    # v2: use transcriptome_mode labels directly (2026-04-02)
    python ${projectDir}/bin/evaluation/tool_end_accuracy_boxplot.py \\
        --input ${evaluation_files} \\
        --output . \\
        --verbose
    """
}

// -----------------------------------------------------------------
// Cross-sample signal-stratified peak recovery + ROC curves.
// Aggregates peak_reason TSVs across all datasets.
// Produces mean ± SD curves for raw, normalised recall, and TPR-FPR ROC.
// -----------------------------------------------------------------
process CrossSamplePeakRoc {
    publishDir "${params.outdir}/summary/${params.test_name}/end_accuracy/peak_roc", mode: 'copy'
    tag "cross_sample_peak_roc"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(reason_tsvs)

    output:
    path "signal_recall_*.png", emit: cross_sample_peak_roc_plot, optional: true
    path "signal_roc_*.png",    optional: true

    script:
    """
    python ${projectDir}/bin/evaluation/peak_roc_curves.py \\
        --input ${reason_tsvs} \\
        --output . \\
        --dataset "${test_name}" \\
        --verbose || true
    """
}

