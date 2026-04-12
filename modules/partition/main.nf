// Module: FlairPartition
// Partitions data to a specific genomic region or subset for testing.
// Supports two modes:
//   1. With pre-computed BED file (from FlairAlign)
//   2. With NO_BED placeholder (from pre-aligned BAM) — generates BED from partitioned BAM

process FlairPartition {
    publishDir "${params.outdir}/partition/${test_name}", mode: 'symlink'
    tag "${dataset_name}_${align_mode}_${partition_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), path(bam), path(bai), path(bed),
          val(partition_mode), val(partition_args), path(genome), path(gtf),
          path(cage_peaks), path(drna_peaks)

    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}.bam"),
          path("${dataset_name}_${align_mode}_${partition_mode}.bam.bai"),
          path("${dataset_name}_${align_mode}_${partition_mode}.bed"),
          path("${dataset_name}_${align_mode}_${partition_mode}_genome.fa"),
          path("${dataset_name}_${align_mode}_${partition_mode}_annotation.gtf"),
          path("${dataset_name}_${align_mode}_${partition_mode}_cage.bed", optional: true),
          path("${dataset_name}_${align_mode}_${partition_mode}_drna.bed", optional: true), emit: partitioned

    script:
    def output_prefix = "${dataset_name}_${align_mode}_${partition_mode}"
    def cage_arg = cage_peaks.name != 'NO_CAGE' ? "--cage-peaks ${cage_peaks}" : ""
    def drna_arg = drna_peaks.name != 'NO_DRNA' ? "--drna-peaks ${drna_peaks}" : ""
    def bed_arg = bed.name != 'NO_BED' ? "--bed ${bed}" : "--generate-bed"

    """
    python ${projectDir}/bin/simple_partition.py \\
        --bam ${bam} \\
        ${bed_arg} \\
        --genome ${genome} \\
        --gtf ${gtf} \\
        ${cage_arg} \\
        ${drna_arg} \\
        --output-prefix ${output_prefix} \\
        ${partition_args}
    """
}
