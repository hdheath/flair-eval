// Module: LraaAssembly
// LRAA (Long Read Alignment Assembler) — isoform discovery + quantification
// from long-read RNA-seq. Distributed only as a Docker image:
//   us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest
// We pull it into a singularity .sif (params.lraa_sif) and run it via
// Nextflow's container directive, mirroring the FlamesAssembly pattern.
//
// Default invocation: reference-guided isoform discovery + quantification.
//   LRAA --genome <fa> --bam <bam> --gtf <reference.gtf> [--HiFi] ...
//
// PacBio HiFi/CCS samples need --HiFi (error rate < 2%); ONT does not.
// The subworkflow (assemble_and_eval.nf) auto-injects --HiFi when the
// dataset's library_type starts with 'pacbio'.
//
// Outputs we consume (per LRAA WDL output set):
//   <prefix>.gtf                   assembled isoform GTF
//   <prefix>.quant.tracking[.gz]   per-read assignments
//   <prefix>.quant.expr            per-transcript counts (col `all_reads`)

process LraaAssembly {
    publishDir "${params.outdir}/assemblers/lraa/${test_name}", mode: 'symlink'
    publishDir "${params.outdir}/logs/${test_name}", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_lraa_${lraa_mode}_${it}" }
    // Benchmark assembler: a broken tool (e.g. FLAMES container, missing dep) must
    // NOT abort the whole comparison run. Ignore failures so the other tools + FLAIR/TED proceed.
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_lraa_${lraa_mode}"
    // The genome is passed as a val (original canonical path, never staged —
    // see partition/main.nf output note) so its directory is NOT auto-mounted
    // into the singularity container. Bind-mount the genome's parent dir so
    // `--genome <abs path>` resolves inside the container.
    containerOptions { "--bind ${file(genome).parent}" }

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path(bam), path(bai), val(genome), path(gtf),
          val(lraa_mode), val(lraa_args)

    output:
    // transcript_counts is included in the main lraa_gtf tuple (alongside
    // gtf + read_map) so the evaluation subworkflow can forward it to
    // ted.py's --counts arg. LRAA's reference-guided mode can carry forward
    // reference transcripts the same way Bambu/IsoQuant do, so the filter
    // (count >= 1) is meaningful here.
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(lraa_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_lraa_${lraa_mode}.gtf"),
          path("${dataset_name}_${align_mode}_${partition_mode}_lraa_${lraa_mode}_read_map.txt"),
          path("${dataset_name}_${align_mode}_${partition_mode}_lraa_${lraa_mode}_transcript_counts.tsv", optional: true),
          emit: lraa_gtf
    path "${dataset_name}_${align_mode}_${partition_mode}_lraa_${lraa_mode}.quant.expr", optional: true, emit: quant_expr
    path "${dataset_name}_${align_mode}_${partition_mode}_lraa_${lraa_mode}.quant.tracking*", optional: true, emit: quant_tracking

    script:
    def output_prefix = "${dataset_name}_${align_mode}_${partition_mode}_lraa_${lraa_mode}"
    def lraa_params = lraa_args ?: ""
    """
    # LRAA writes outputs in CWD with its --output_prefix; we rename to the
    # pipeline's canonical underscore-delimited prefix below.
    LRAA \\
        --genome ${genome} \\
        --bam ${bam} \\
        --gtf ${gtf} \\
        --output_prefix ${output_prefix} \\
        --num_parallel_contigs ${task.cpus} \\
        ${lraa_params}

    # Defensive: LRAA can emit either ${output_prefix}.quant.tracking or
    # ${output_prefix}.quant.tracking.gz depending on version; convert_read_map.py
    # auto-detects .gz so we just point at whichever exists.
    TRACKING=""
    for cand in ${output_prefix}.quant.tracking ${output_prefix}.quant.tracking.gz; do
        if [ -s "\$cand" ]; then TRACKING="\$cand"; break; fi
    done
    if [ -n "\$TRACKING" ]; then
        python ${projectDir}/bin/convert_read_map.py \\
            --lraa-tracking "\$TRACKING" \\
            --output ${output_prefix}_read_map.txt \\
            --verbose
    else
        echo "WARNING: No LRAA quant.tracking file produced; writing empty read map" >&2
        : > ${output_prefix}_read_map.txt
    fi

    # Per-transcript counts file in the IsoQuant convention
    # (#feature_id<TAB>count). Consumed by ted.py --counts to filter
    # isoforms_observed to transcripts with read support >= 1.
    if [ -s "${output_prefix}.quant.expr" ]; then
        python ${projectDir}/bin/convert_read_map.py \\
            --lraa-quant-expr "${output_prefix}.quant.expr" \\
            --counts-output ${output_prefix}_transcript_counts.tsv \\
            --verbose
    else
        echo "WARNING: No LRAA quant.expr file produced; skipping counts TSV" >&2
    fi
    """
}
