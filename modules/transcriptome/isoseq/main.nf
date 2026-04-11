// Module: IsoSeqAssembly
// PacBio IsoSeq: collapse aligned HiFi/CCS reads into transcript isoforms.
// Uses `isoseq collapse` on an already-aligned BAM to produce a GFF and
// a read_stat.txt mapping reads to collapsed isoform clusters.
//
// Key parameters (set via isoseq_args in JSON config):
//   --do-not-collapse-extra-5exons: Keep 5' degraded transcripts separate
//   --min-aln-coverage <float>:     Minimum alignment coverage (default 0.99)
//   --min-aln-identity <float>:     Minimum alignment identity (default 0.95)

process IsoSeqAssembly {
    publishDir "${params.outdir}/isoseq/${test_name}", mode: 'symlink'
    publishDir "${params.outdir}/logs/${test_name}", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_isoseq_${isoseq_mode}_${it}" }
    errorStrategy 'terminate'
    tag "${dataset_name}_${align_mode}_${partition_mode}_isoseq_${isoseq_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path(bam), path(bai), path(genome), path(gtf),
          val(isoseq_mode), val(isoseq_args)

    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(isoseq_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_isoseq_${isoseq_mode}.gff"),
          path("${dataset_name}_${align_mode}_${partition_mode}_isoseq_${isoseq_mode}_read_map.txt"),
          emit: isoseq_gff

    script:
    def output_prefix = "${dataset_name}_${align_mode}_${partition_mode}_isoseq_${isoseq_mode}"
    """
    # IsoSeq collapse: aligned BAM → collapsed isoforms GFF + read stat
    isoseq collapse \\
        ${bam} \\
        ${output_prefix}.gff \\
        --do-not-collapse-extra-5exons \\
        -j ${task.cpus} \\
        ${isoseq_args}

    # Convert the read_stat.txt to FLAIR-format read map
    python ${projectDir}/bin/convert_read_map.py \\
        --isoseq-read-stat ${output_prefix}.read_stat.txt \\
        --output ${output_prefix}_read_map.txt \\
        --verbose
    """
}
