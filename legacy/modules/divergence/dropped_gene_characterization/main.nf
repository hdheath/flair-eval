// Module: DroppedGeneCharacterization
//
// Companion to GeneRecoveryAnalysis. For each per-sample dropped gene
// list (FLAIR baseline has isos, TED mode dropped them all), characterize
// the dropped genes by:
//   - gene biotype (protein_coding vs lncRNA vs pseudogene vs unknown)
//   - annotated transcript length distribution (vs kept genes)
//   - n annotated transcripts per gene (vs kept genes)
//
// Answers: are TED-dropped genes systematically (a) too long for reads
// to span [length], (b) sparsely-annotated non-coding [biotype], or
// (c) something else?
//
// On the 05_08_26_chr22 run, the answer across all 3 samples was clear:
// it's biotype, not length. lncRNAs and unknown-biotype genes drop at
// 50–95%; protein_coding genes drop at 14–19%; transcript length is
// indistinguishable between dropped and kept sets.

process DroppedGeneCharacterization {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/cds_analysis", mode: 'copy'
    publishDir "${params.outdir}/summary/${params.test_name}/cds_analysis", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files), path(gtf)

    output:
    path "${test_name}_dropped_gene_chars.tsv", optional: true
    path "${test_name}_dropped_gene_summary.tsv", optional: true
    path "${test_name}_dropped_gene_panel.png", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    def baseline = params.gene_recovery_baseline ?: 'FLAIR'
    """
    python ${projectDir}/bin/evaluation/dropped_gene_characterization.py \\
        --bed ${bed_args.join(' ')} \\
        --gtf ${gtf} \\
        --baseline-label ${baseline} \\
        --output-prefix ${test_name} \\
        --output-dir . \\
        --verbose || true
    """
}
