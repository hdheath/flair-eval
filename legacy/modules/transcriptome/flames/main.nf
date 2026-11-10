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
    // Benchmark assembler: a broken tool (e.g. FLAMES container, missing dep) must
    // NOT abort the whole comparison run. Ignore failures so the other tools + FLAIR/TED proceed.
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_flames_${flames_mode}"
    // The genome is passed as a val (original canonical path, never staged —
    // see partition/main.nf output note) so its directory is NOT auto-mounted
    // into the singularity container. Bind-mount the genome's parent dir so
    // `genome_fa = <abs path>` resolves inside the container.
    containerOptions { "--bind ${file(genome).parent}" }

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path(bam), path(bai), val(genome), path(gtf),
          val(flames_mode), val(flames_args)

    output:
    // transcript_counts is included in the main flames_gtf tuple (alongside
    // gtf + read_map) so the evaluation subworkflow can forward it to
    // ted.py's --counts arg.  FLAMES' bulk_long_pipeline emits
    // `transcript_count.csv.gz` keyed on transcript_id; the GFF3 may include
    // transcripts with zero re-quantified reads, so the count filter is
    // genuinely useful here (mirrors Bambu/IsoQuant rationale).
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(flames_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_flames_${flames_mode}.gff3"),
          path("${dataset_name}_${align_mode}_${partition_mode}_flames_${flames_mode}_read_map.txt"),
          path("${dataset_name}_${align_mode}_${partition_mode}_flames_${flames_mode}_transcript_counts.tsv", optional: true),
          emit: flames_gtf
    path "${dataset_name}_${align_mode}_${partition_mode}_flames_${flames_mode}_transcript_counts.tsv", optional: true, emit: transcript_counts

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

    # Step 5: Build a per-transcript counts TSV in IsoQuant's
    # `#feature_id<TAB>count` convention (consumed by ted.py --counts).
    # FLAMES emits transcript_count.csv(.gz) with header
    # `transcript_id,gene_id,<sample1>[,<sample2>...]`.  We sum the numeric
    # sample columns and skip the gene_id column, matching the Bambu path.
    COUNT_CSV=""
    for cand in flames_out/transcript_count.csv.gz flames_out/transcript_count.csv; do
        if [ -s "\$cand" ]; then COUNT_CSV="\$cand"; break; fi
    done
    if [ -n "\$COUNT_CSV" ]; then
        python - <<PYEOF
import csv, gzip, sys
src = "\$COUNT_CSV"
opener = gzip.open if src.endswith('.gz') else open
out = "${output_prefix}_transcript_counts.tsv"
n = 0
with opener(src, 'rt') as fh, open(out, 'w') as ofh:
    reader = csv.reader(fh)
    header = next(reader, None)
    if header is None:
        sys.exit(0)
    # Header: transcript_id, gene_id, <one or more sample cols>
    sample_cols = list(range(2, len(header)))
    ofh.write("#feature_id\\tcount\\n")
    for row in reader:
        if not row:
            continue
        tid = row[0].strip()
        if not tid:
            continue
        total = 0.0
        for i in sample_cols:
            if i >= len(row):
                continue
            try:
                total += float(row[i])
            except ValueError:
                continue
        ofh.write(f"{tid}\\t{int(total) if total == int(total) else total}\\n")
        n += 1
print(f"FLAMES counts: wrote {n} transcripts -> {out}", file=sys.stderr)
PYEOF
    else
        echo "WARNING: No FLAMES transcript_count.csv[.gz] found; skipping counts TSV" >&2
    fi
    """
}
