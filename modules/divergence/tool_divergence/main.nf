// Module: ToolDivergence
// Section 1.1: Quantify divergence between assemblers.
// Computes pairwise Jaccard for splice junctions vs transcript ends,
// and quantifies regulatory motif collapse per tool.

process ToolDivergence {
    publishDir "${params.outdir}/evaluations/${test_name}/divergence", mode: 'copy'
    publishDir "${params.outdir}/summary/${params.test_name}/divergence", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(isoform_files), val(tool_names),
          path(map_files), val(map_names),
          path(reads_bed), path(genome)

    output:
    path "${test_name}_divergence_dashboard.png", optional: true
    path "${test_name}_jaccard_heatmaps.png", optional: true
    path "${test_name}_jaccard_comparison.png", optional: true
    path "${test_name}*_motif_collapse*.png", optional: true
    path "${test_name}_divergence_metrics.tsv", optional: true, emit: divergence_metrics

    script:
    """
    # Convert any GTF files to BED12 for uniform processing
    for f in *.gtf; do
        if [ -f "\$f" ]; then
            bed_name="\${f%.gtf}.bed"
            python ${projectDir}/bin/gtf_to_bed12.py --gtf "\$f" --output "\$bed_name" || true
        fi
    done

    # Build --isoform-beds arguments from tool_names
    ISOFORM_ARGS=""
    for pair in ${tool_names}; do
        tool=\$(echo "\$pair" | cut -d: -f1)
        fname=\$(echo "\$pair" | cut -d: -f2)
        if [[ "\$fname" == *.gtf ]]; then
            bed_fname="\${fname%.gtf}.bed"
            if [ -f "\$bed_fname" ]; then
                ISOFORM_ARGS="\$ISOFORM_ARGS \$tool:\$bed_fname"
            fi
        else
            ISOFORM_ARGS="\$ISOFORM_ARGS \$tool:\$fname"
        fi
    done

    # Build --map-files arguments from map_names
    MAP_ARGS=""
    for pair in ${map_names}; do
        tool=\$(echo "\$pair" | cut -d: -f1)
        mfname=\$(echo "\$pair" | cut -d: -f2)
        if [ -f "\$mfname" ] && [ -s "\$mfname" ]; then
            MAP_ARGS="\$MAP_ARGS \$tool:\$mfname"
        fi
    done

    python ${projectDir}/bin/evaluation/tool_divergence_plots.py \\
        --isoform-beds \$ISOFORM_ARGS \\
        --map-files \$MAP_ARGS \\
        --reads-bed ${reads_bed} \\
        --genome ${genome} \\
        --output-prefix ${test_name} \\
        --output-dir . \\
        --window-5prime 50 \\
        --window-3prime 5 \\
        --title-prefix "${test_name}: " \\
        --dpi 300 \\
        --verbose || true
    """
}
