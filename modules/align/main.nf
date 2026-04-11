// Module: FlairAlign
// Aligns long-read RNA-seq data to the genome using minimap2 (via FLAIR)

process FlairAlign {
    publishDir "${params.outdir}/align", mode: 'symlink'
    publishDir "${params.outdir}/logs/align", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${it}" }
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}"

    input:
    tuple val(test_name), val(dataset_name), path(reads), val(align_mode), val(align_args), path(genome)

    output:
    tuple val(test_name), val(dataset_name), val(align_mode),
          path("${dataset_name}_${align_mode}.bam"),
          path("${dataset_name}_${align_mode}.bam.bai"),
          path("${dataset_name}_${align_mode}.bed"), emit: alignments

    script:
    """
    flair align ${align_args} \\
        -r ${reads} \\
        -g ${genome} \\
        -o ${dataset_name}_${align_mode}
    """
}
