// Module: UTRFeatures
// Section 1.2: Compare 5'UTR and 3'UTR features (length, GC, uORFs) across
// assemblers vs reference annotation. Inspired by Weber et al. 2023 (Oncogene).

process UTRFeatures {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/utr", mode: 'copy'
    publishDir "${params.outdir}/summary/${params.test_name}/utr", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(isoform_files), val(tool_names),
          path(genome), path(gtf)

    output:
    path "${test_name}_utr5_dashboard.png", optional: true
    path "${test_name}_utr5_length.png", optional: true
    path "${test_name}_utr5_gc.png", optional: true
    path "${test_name}_utr5_uorfs.png", optional: true
    path "${test_name}_utr3_dashboard.png", optional: true
    path "${test_name}_utr3_length.png", optional: true
    path "${test_name}_utr3_gc.png", optional: true
    path "${test_name}_utr_features.tsv", optional: true

    script:
    """
    # Convert any GTF isoform files to BED12 for uniform processing
    for f in *.gtf; do
        if [ -f "\$f" ] && [ "\$f" != "${gtf.name}" ]; then
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

    python ${projectDir}/bin/evaluation/utr_feature_plots.py \\
        --isoform-beds \$ISOFORM_ARGS \\
        --gtf ${gtf} \\
        --genome ${genome} \\
        --output-prefix ${test_name} \\
        --output-dir . \\
        --title-prefix "${test_name}: " \\
        --dpi 300 \\
        --verbose || true
    """
}
