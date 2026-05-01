// Module: StringTie2Assembly
// StringTie2: efficient transcript assembly from long and short reads.
// Uses StringTie2's long-read mode (-L) with reference-guided annotation.
//
// Key parameters (set via stringtie2_args in JSON config):
//   -c <float>:  Minimum reads per bp coverage (default 1)
//   -f <float>:  Minimum isoform abundance fraction (default 0.01)
//   -m <int>:    Minimum transcript length (default 200)
//   -a <int>:    Minimum anchor length for junctions (default 10)
//   -j <float>:  Minimum junction coverage (default 1)
//   -M <float>:  Fraction for multi-hit reads (default 1.0)

process StringTie2Assembly {
    publishDir "${params.outdir}/assemblers/stringtie2/${test_name}", mode: 'symlink'
    publishDir "${params.outdir}/logs/${test_name}", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_stringtie2_${stringtie2_mode}_${it}" }
    errorStrategy 'terminate'
    tag "${dataset_name}_${align_mode}_${partition_mode}_stringtie2_${stringtie2_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path(bam), path(bai), path(genome), path(gtf),
          val(stringtie2_mode), val(stringtie2_args)

    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(stringtie2_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_stringtie2_${stringtie2_mode}.gtf"),
          path("${dataset_name}_${align_mode}_${partition_mode}_stringtie2_${stringtie2_mode}_read_map.txt"),
          emit: stringtie2_gtf

    script:
    def output_prefix = "${dataset_name}_${align_mode}_${partition_mode}_stringtie2_${stringtie2_mode}"
    """
    # StringTie2 long-read mode: -L for long reads, -G for guided annotation
    stringtie ${bam} \\
        -L \\
        -G ${gtf} \\
        -o ${output_prefix}.gtf \\
        -p ${task.cpus} \\
        ${stringtie2_args}

    # StringTie2 does not emit per-read transcript assignments.
    # Keep the expected file path, but leave it empty so Evaluation skips
    # read-assignment and read-end entropy metrics instead of using a
    # transcript_id -> transcript_id placeholder map.
    : > ${output_prefix}_read_map.txt
    """
}
