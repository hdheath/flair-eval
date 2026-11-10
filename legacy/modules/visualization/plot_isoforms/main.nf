// Module: PlotIsoforms
// Generates visualization plots of isoform structures and read assignments.
// Only runs for partitions <= 400kb to avoid overly complex plots.

process PlotIsoforms {
    publishDir "${params.outdir}/isoform_plots/${test_name}/${transcriptome_mode}", mode: 'copy'
    publishDir "${params.outdir}/logs/${test_name}", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_plot_${it}" }
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), val(transcriptome_mode),
          path(bam), path(bai), val(region_spec), path(cage_peaks), path(drna_peaks),
          val(cage_signal_plus), val(cage_signal_minus),
          val(drna_signal_plus), val(drna_signal_minus),
          path(isoforms_bed), path(isoform_read_map)

    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), val(transcriptome_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_isoform_plot.png"), emit: plots
    path "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_sj_chain_end_scatter_*.png", optional: true, emit: chain_scatter_plots

    script:
    def cage_peaks_arg = cage_peaks.name != 'NO_CAGE' ? "--cage-peaks ${cage_peaks}" : ""
    def drna_peaks_arg = drna_peaks.name != 'NO_DRNA' ? "--drna-peaks ${drna_peaks}" : ""
    def cage_signal_plus_arg = cage_signal_plus ? "--cage-signal-plus ${cage_signal_plus}" : ""
    def cage_signal_minus_arg = cage_signal_minus ? "--cage-signal-minus ${cage_signal_minus}" : ""
    def drna_signal_plus_arg = drna_signal_plus ? "--drna-signal-plus ${drna_signal_plus}" : ""
    def drna_signal_minus_arg = drna_signal_minus ? "--drna-signal-minus ${drna_signal_minus}" : ""
    def region_arg = region_spec ? "--region ${region_spec}" : ""
    """
    python ${projectDir}/bin/isoform_plot_test.py \\
        --bam ${bam} \\
        --readmap ${isoform_read_map} \\
        --isoforms ${isoforms_bed} \\
        ${region_arg} \\
        ${cage_peaks_arg} \\
        ${drna_peaks_arg} \\
        ${cage_signal_plus_arg} \\
        ${cage_signal_minus_arg} \\
        ${drna_signal_plus_arg} \\
        ${drna_signal_minus_arg} \\
        --max-redundant-pattern-reads 50 \\
        --pattern-end-tolerance 1 \\
        --output ${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_isoform_plot

    if [ -n "${region_spec}" ]; then
        python ${projectDir}/bin/read_end_chain_scatter.py \\
            --bam ${bam} \\
            --readmap ${isoform_read_map} \\
            --isoforms ${isoforms_bed} \\
            --region ${region_spec} \\
            ${cage_peaks_arg} \\
            ${drna_peaks_arg} \\
            ${cage_signal_plus_arg} \\
            ${cage_signal_minus_arg} \\
            ${drna_signal_plus_arg} \\
            ${drna_signal_minus_arg} \\
            --max-chains 10 \\
            --focus-padding 500 \\
            --max-focus-span 100000 \\
            --max-region-size 500000 \\
            --signal-bins 300 \\
            --output-prefix ${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}
    fi
    """
}
