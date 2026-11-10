// Module: AltEndCDSAnalysis
//
// For each junction-chain (SJC) group with >=2 isoforms, compute UTR
// length deltas and productivity-class transitions between sibling
// isoforms. Surfaces:
//   * does TED's alt-TSS choice change CDS / 5'UTR length?
//   * does alt-TTS introduce premature termination (PRO->PTC)?
//
// Requires CDS-aware BEDs from FlairTranscriptome's predict_cds path
// (the *.isoforms.CDS.bed files), where predictProductivity has
// populated thickStart/thickEnd and embedded productivity (PRO/PTC/
// NST/NGO) in the name field.
//
// On BEDs without CDS (no --predict_cds), the script logs a warning
// and skips the file gracefully — no failure.

process AltEndCDSAnalysis {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/cds_analysis", mode: 'copy'
    publishDir "${params.outdir}/summary/${params.test_name}/cds_analysis", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files)

    output:
    path "${test_name}_alt_end_cds.tsv", optional: true
    path "${test_name}_alt_end_cds_pairs.tsv", optional: true
    path "${test_name}_alt_end_cds_panel.png", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/alt_end_cds_analysis.py \\
        --bed ${bed_args.join(' ')} \\
        --output-prefix ${test_name} \\
        --output-dir . \\
        --verbose || true
    """
}
