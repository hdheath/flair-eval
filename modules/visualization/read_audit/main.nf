// Module: ReadAudit + ReadAuditSummary
// Per-mode read classification audit (BED9 for IGV + TSV detail),
// followed by cross-mode summary comparison plot.

process ReadAudit {
    publishDir "${params.outdir}/evaluations/${test_name}/per_method/read_audit", mode: 'copy'
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_read_audit"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(transcriptome_mode),
          path(isoforms_bed), path(isoform_read_map),
          path(bam), path(bai)

    output:
    tuple val(test_name), val(transcriptome_mode),
          path("${output_prefix}.read_audit.bed"),
          path("${output_prefix}.read_audit.tsv"), emit: audit_results

    script:
    output_prefix = "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}"
    """
    python ${projectDir}/bin/evaluation/read_audit.py \\
        --bam ${bam} \\
        --read-map ${isoform_read_map} \\
        --isoforms-bed ${isoforms_bed} \\
        --output ${output_prefix} \\
        --tolerance 50 \\
        --verbose
    """
}

process ReadAuditSummary {
    publishDir "${params.outdir}/evaluations/${test_name}/summary/read_audit", mode: 'copy'
    errorStrategy 'ignore'
    tag "${test_name}_read_audit_summary"

    input:
    tuple val(test_name), val(modes), path(audit_tsvs)

    output:
    path "read_audit/*.{png,svg}", emit: read_audit_plot, optional: true

    script:
    // Build mode:path arguments
    def input_args = []
    for (int i = 0; i < modes.size(); i++) {
        input_args << "${modes[i]}:${audit_tsvs[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/read_audit_plot.py \\
        --input ${input_args.join(' ')} \\
        --output read_audit \\
        --verbose || true
    """
}
