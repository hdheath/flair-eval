// Module: PrepareReferencePeaks
// Extracts TSS and TTS from reference annotation GTF for use as ground truth in TED evaluation.

process PrepareReferencePeaks {
    publishDir "${params.outdir}/reference_peaks/${test_name}", mode: 'symlink'
    tag "${dataset_name}_${align_mode}_${partition_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), path(gtf)

    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_ref_tss.bed"),
          path("${dataset_name}_${align_mode}_${partition_mode}_ref_tts.bed"), emit: ref_peaks

    script:
    """
    python ${projectDir}/bin/gtf_to_tss_tts.py \\
        --gtf ${gtf} \\
        --output-prefix ${dataset_name}_${align_mode}_${partition_mode}_ref \\
        --deduplicate
    """
}
