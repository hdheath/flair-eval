// Module: IsoQuantAssembly
// IsoQuant: Python-based transcript reconstruction from long reads.
// Key parameters (set via isoquant_args in JSON config):
//   --data_type: nanopore, pacbio_ccs, or assembly (REQUIRED)
//   --complete_genedb: Use for complete gene databases like GENCODE
//   --model_construction_strategy: sensitive, default, or fl_pacbio

process IsoQuantAssembly {
    publishDir "${params.outdir}/assemblers/isoquant/${test_name}", mode: 'symlink'
    publishDir "${params.outdir}/logs/${test_name}", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}_${it}" }
    errorStrategy 'terminate'
    tag "${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path(bam), path(bai), path(genome), path(gtf),
          val(isoquant_mode), val(isoquant_args)

    output:
    // transcript_counts is included in the main isoquant_gtf tuple (alongside
    // gtf + read_map) so the evaluation subworkflow can forward it to ted.py's
    // --counts arg. Without this filter, every GTF transcript gets counted —
    // including ~30% IsoQuant carries forward from reference with zero reads.
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(isoquant_mode), path("${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}.gtf"),
          path("${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}_read_map.txt"),
          path("${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}_transcript_counts.tsv", optional: true),
          emit: isoquant_gtf
    path "${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}_read_assignments.tsv", optional: true, emit: read_assignments
    path "${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}_transcript_counts.tsv", optional: true, emit: transcript_counts

    script:
    def output_prefix = "${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}"
    // Absolute path to the isoquant env's python (2026-05-14): the conda
    // activate emitted by Nextflow's process directive works on the login
    // node but on some SLURM compute nodes does NOT prepend the env bin to
    // PATH, so `isoquant.py` runs under the system python (which lacks
    // gffutils). Pinning to the env's python bypasses the activation.
    def isoquant_py = '/private/home/hdheath/miniforge3/envs/isoquant/bin/python3.8 /private/home/hdheath/miniforge3/envs/isoquant/bin/isoquant.py'
    """
    ${isoquant_py} \\
        --reference ${genome} \\
        --genedb ${gtf} \\
        --bam ${bam} \\
        --output isoquant_out \\
        --threads ${task.cpus} \\
        --prefix ${output_prefix} \\
        ${isoquant_args}

    cp isoquant_out/${output_prefix}/${output_prefix}.transcript_models.gtf ${output_prefix}.gtf

    READ_MAP_CREATED=false
    for ext in .transcript_model_reads.tsv.gz .transcript_model_reads.tsv; do
        if [ -f "isoquant_out/${output_prefix}/${output_prefix}\${ext}" ]; then
            python ${projectDir}/bin/convert_read_map.py \\
                --isoquant-model-reads "isoquant_out/${output_prefix}/${output_prefix}\${ext}" \\
                --output ${output_prefix}_read_map.txt \\
                --verbose
            READ_MAP_CREATED=true
            break
        fi
    done
    if [ "\$READ_MAP_CREATED" = false ]; then
        echo "WARNING: No transcript_model_reads file found, creating empty read map" >&2
        touch ${output_prefix}_read_map.txt
    fi

    if [ -f "isoquant_out/${output_prefix}/${output_prefix}.read_assignments.tsv" ]; then
        cp isoquant_out/${output_prefix}/${output_prefix}.read_assignments.tsv ${output_prefix}_read_assignments.tsv
    fi
    if [ -f "isoquant_out/${output_prefix}/${output_prefix}.transcript_counts.tsv" ]; then
        cp isoquant_out/${output_prefix}/${output_prefix}.transcript_counts.tsv ${output_prefix}_transcript_counts.tsv
    fi
    """
}
