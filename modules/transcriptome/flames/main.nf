// Module: FlamesAssembly
// FLAMES (Full-Length Analysis of Mutations and Splicing): isoform detection
// and quantification from long-read RNA-seq data.
//
// FLAMES is an R/Bioconductor package that runs its own alignment pipeline.
// We use Nextflow's native Singularity support via the 'container' directive
// since the authors warn against conda installs. The workflow:
//   1. Convert input BAM to FASTQ (samtools fastq, available in the container)
//   2. Run FLAMES bulk_long_pipeline via Rscript
//   3. Collect isoform_annotated.gff3 and build read map from realign BAM

process FlamesAssembly {
    publishDir "${params.outdir}/assemblers/flames/${test_name}", mode: 'symlink'
    publishDir "${params.outdir}/logs/${test_name}", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_flames_${flames_mode}_${it}" }
    errorStrategy 'terminate'
    tag "${dataset_name}_${align_mode}_${partition_mode}_flames_${flames_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path(bam), path(bai), path(genome), path(gtf),
          val(flames_mode), val(flames_args)

    output:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(flames_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_flames_${flames_mode}.gff3"),
          path("${dataset_name}_${align_mode}_${partition_mode}_flames_${flames_mode}_read_map.txt"),
          emit: flames_gtf

    script:
    def output_prefix = "${dataset_name}_${align_mode}_${partition_mode}_flames_${flames_mode}"
    """
    # Step 1: Convert BAM to FASTQ for FLAMES input
    samtools fastq -@ ${task.cpus} ${bam} > reads.fastq

    # Step 2: Set environment for container execution
    # - HOME/TMPDIR: writable locations for cache/temp files
    # - BASILISK_USE_SYSTEM_DIR: use the container's pre-installed Python
    #   packages instead of building a fresh venv (avoids oarfish build issues)
    export HOME=\$PWD
    export XDG_CACHE_HOME=\$PWD/.cache
    export R_USER_CACHE_DIR=\$PWD/.cache
    export TMPDIR=\$PWD/tmp
    export BASILISK_USE_SYSTEM_DIR=TRUE
    mkdir -p \$PWD/tmp \$PWD/.cache

    # Step 3: Run FLAMES bulk_long_pipeline
    Rscript -e '
    library(FLAMES)
    outdir <- "flames_out"
    dir.create(outdir, showWarnings = FALSE)
    bulk_long_pipeline(
        annotation = "${gtf}",
        fastq = "reads.fastq",
        outdir = outdir,
        genome_fa = "${genome}"
    )
    '

    # Step 3: Collect isoform annotation
    cp flames_out/isoform_annotated.gff3 ${output_prefix}.gff3

    # Step 4: Convert realign2transcript BAM to FLAIR read map format
    # (samtools+awk: extract primary alignments, group reads by transcript)
    samtools view -F 0x904 flames_out/realign2transcript.bam | \
        awk -F'\t' '{print \$3 "\t" \$1}' | \
        sort -t'\t' -k1,1 | \
        awk -F'\t' '{
            if (\$1 != prev) {
                if (prev != "") print prev "\t" reads;
                prev = \$1; reads = \$2
            } else {
                reads = reads "," \$2
            }
        } END { if (prev != "") print prev "\t" reads }' \
        > ${output_prefix}_read_map.txt
    """
}
