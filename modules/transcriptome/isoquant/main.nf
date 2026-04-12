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
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(isoquant_mode), path("${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}.gtf"),
          path("${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}_read_map.txt"),
          emit: isoquant_gtf
    path "${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}_read_assignments.tsv", optional: true, emit: read_assignments
    path "${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}_transcript_counts.tsv", optional: true, emit: transcript_counts

    script:
    def output_prefix = "${dataset_name}_${align_mode}_${partition_mode}_isoquant_${isoquant_mode}"
    """
    isoquant.py \\
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
