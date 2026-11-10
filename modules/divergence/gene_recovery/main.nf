// Module: GeneRecoveryAnalysis
//
// For each non-baseline transcriptome mode, diff its per-gene isoform
// count against a baseline (default: FLAIR). Surfaces:
//   * genes that lost ALL isoforms ("dropped")
//   * genes that went from N>=2 -> 1 isoform ("collapsed")
//   * partial / unchanged / expanded counts for context
//
// Outputs land in `evaluations/per_sample/<test>/cds_analysis/` so
// CrossSampleGeneRecovery (Track 4) can pick them up alongside the
// CDS / UTR / uORF analyses.
//
// Input shape mirrors IsoformsPerGeneHist: per-method isoform files
// staged together with parallel `bed_labels` for the analysis to know
// which file is which method. No genome / GTF needed — gene_id comes
// from each BED's name field, which FLAIR encodes as <iso_id>_<gene_id>.

process GeneRecoveryAnalysis {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/cds_analysis", mode: 'copy'
    publishDir "${params.outdir}/summary/${params.test_name}/cds_analysis", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files)

    output:
    path "${test_name}_gene_recovery_dropped.tsv", optional: true
    path "${test_name}_gene_recovery_collapsed.tsv", optional: true
    path "${test_name}_gene_recovery_summary.tsv", optional: true
    path "${test_name}_gene_recovery_panel.png", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    // Baseline label is the FLAIR mode (no --ted) by convention. If a
    // test set uses a different name, override via params.gene_recovery_baseline.
    def baseline = params.gene_recovery_baseline ?: 'FLAIR'
    """
    python ${projectDir}/bin/evaluation/gene_recovery_analysis.py \\
        --bed ${bed_args.join(' ')} \\
        --baseline-label ${baseline} \\
        --output-prefix ${test_name} \\
        --output-dir . \\
        --verbose || true
    """
}
