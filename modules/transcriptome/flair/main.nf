// Module: FlairTranscriptome
// Runs FLAIR transcriptome assembly to generate isoform models.

process FlairTranscriptome {
    publishDir "${params.outdir}/assemblers/flair/${test_name}", mode: 'symlink'
    publishDir "${params.outdir}/assemblers/flair/${test_name}/firstpass", mode: 'symlink', pattern: '*.firstpass*.bed'
    publishDir "${params.outdir}/logs/${test_name}", mode: 'copy', pattern: '*.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_transcriptome.${it.tokenize('.')[-1]}" }
    errorStrategy 'terminate'
    tag "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), val(partition_args),
          path(bam), path(bai), val(genome), path(gtf),
          val(transcriptome_mode), val(transcriptome_args), path(junction_tab)

    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), val(partition_args), val(transcriptome_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome.isoforms.bed", optional: true),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome.isoforms.gtf", optional: true),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome.isoforms.fa", optional: true),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome.isoform.counts.txt", optional: true),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome.isoform.read.map.txt", optional: true),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome.ted_log.tsv", optional: true), emit: transcriptome
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), val(partition_args), val(transcriptome_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome.firstpass.bed", optional: true), emit: firstpass
    // CDS-aware BED produced by --predict_cds (predictProductivity).
    // File path: <out_prefix>.isoforms.CDS.bed; the original .isoforms.bed
    // is left untouched. Optional because some pipeline modes might not
    // run predict_cds (e.g. FLAMES, IsoSeq don't generate this file).
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode), val(partition_args), val(transcriptome_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome.isoforms.CDS.bed", optional: true),
          path("${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome.isoforms.CDS.info.tsv", optional: true), emit: cds_bed

    script:
    // Auto-detect junction file format: .bed → --junction_bed, .tab → --junction_tab
    def has_junction_flag = transcriptome_args.contains('--junction_tab') && junction_tab.name != 'NO_JUNCTION_TAB'
    def junction_flag = junction_tab.name.endsWith('.bed') ? '--junction_bed' : '--junction_tab'
    def junction_tab_arg = has_junction_flag ? "${junction_flag} ${junction_tab.name}" : ""
    def cleaned_args = transcriptome_args.replaceAll('--junction_tab\\s*', '').replaceAll('--ted_log\\s*', '').replaceAll('--predict_cds\\s*', '')
    def is_ted = (transcriptome_args =~ /(?:^|\s)--ted(?:\s|$)/).find()
    def ted_log_flag = is_ted ? '--ted_log' : ''
    def ted_log_file = "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome.ted_log.tsv"
    // Always run predictProductivity so downstream CDS/UTR/uORF analyses
    // (GeneRecoveryAnalysis, AltEndCDSAnalysis, CDSCoverageRecovery)
    // have CDS coordinates + productivity annotation in the final BED.
    // NOTE: the new dev FLAIR (flair-ted-52826) makes CDS/productivity
    // prediction DEFAULT/inline and REMOVED the --predict_cds flag, so we no
    // longer pass it (doing so errors: "unrecognized arguments: --predict_cds").
    // Strip it if a mode's transcriptome_args still carries it, and never add it.
    def predict_cds_flag = ''

    // The new dev FLAIR (flair-ted-52826) REQUIRES --sample_name (old flair4 did
    // not have it). Use the same base name the -o output uses.
    def out_base = "${dataset_name}_${align_mode}_${partition_mode}_${transcriptome_mode}_transcriptome"
    """
    flair transcriptome \\
        -b ${bam} \\
        --genome ${genome} \\
        -f ${gtf} \\
        -t ${task.cpus} \\
        --sample_name ${out_base} \\
        --keep_intermediate \\
        ${junction_tab_arg} \\
        ${ted_log_flag} \\
        ${predict_cds_flag} \\
        ${cleaned_args} \\
        -o ${out_base}
    touch ${ted_log_file}
    """
}
