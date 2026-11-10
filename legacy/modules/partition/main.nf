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
          val(partition_mode), val(partition_args), path(genome), val(genome_orig), path(gtf),
          path(cage_peaks), path(drna_peaks)

    output:
    // NOTE: the reference genome is NOT a materialized output. It is the same
    // 2.8 GB canonical reference for every sample/partition and is FLAIR-read
    // directly. Emitting it as a path() output caused Nextflow's scratch
    // staging (`nxf_fs_copy`, fire_test scratch=/data/tmp) to DEREFERENCE the
    // symlink and copy the full 2.8 GB genome back into the workDir for EVERY
    // partition task → blew the /private/groups disk quota (the WTC11
    // FlairPartition "Disk quota exceeded" failures). Instead we pass the
    // ORIGINAL genome path through as a `val` so it is never staged or copied;
    // downstream consumers take it as `val(genome)` and read the full genome
    // (correct coordinates) directly from its canonical location.
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}.bam"),
          path("${dataset_name}_${align_mode}_${partition_mode}.bam.bai"),
          path("${dataset_name}_${align_mode}_${partition_mode}.bed"),
          val(genome_orig),
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
        --threads ${task.cpus} \\
        --output-prefix ${output_prefix} \\
        ${partition_args}
    """
}
