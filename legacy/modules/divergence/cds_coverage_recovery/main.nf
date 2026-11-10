// Module: CDSCoverageRecovery
//
// Per-gene productivity-class diff vs FLAIR baseline. Surfaces the
// PRO -> PTC and PRO -> NGO transitions caused by TED's TSS/TTS choices,
// plus the genes where productivity is retained.
//
// Requires CDS-aware BEDs (*.isoforms.CDS.bed) from FlairTranscriptome's
// --predict_cds path. Without those, the analysis logs a warning and
// produces no output.

process CDSCoverageRecovery {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/cds_analysis", mode: 'copy'
    publishDir "${params.outdir}/summary/${params.test_name}/cds_analysis", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files)

    output:
    path "${test_name}_cds_coverage_recovery.tsv", optional: true
    path "${test_name}_cds_coverage_recovery_summary.tsv", optional: true
    path "${test_name}_cds_coverage_recovery_panel.png", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    def baseline = params.gene_recovery_baseline ?: 'FLAIR'
    """
    python ${projectDir}/bin/evaluation/cds_coverage_recovery.py \\
        --bed ${bed_args.join(' ')} \\
        --baseline-label ${baseline} \\
        --output-prefix ${test_name} \\
        --output-dir . \\
        --verbose || true
    """
}
