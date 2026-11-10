// Module: IsoSeqAssembly
// PacBio IsoSeq: collapse aligned HiFi/CCS reads into transcript isoforms.
// Uses `isoseq collapse` on an already-aligned BAM to produce a GFF and
// a read_stat.txt mapping reads to collapsed isoform clusters.
//
// Key parameters (set via isoseq_args in JSON config):
//   --do-not-collapse-extra-5exons: Keep 5' degraded transcripts separate
//   --min-aln-coverage <float>:     Minimum alignment coverage (default 0.99)
//   --min-aln-identity <float>:     Minimum alignment identity (default 0.95)

process IsoSeqAssembly {
    publishDir "${params.outdir}/assemblers/isoseq/${test_name}", mode: 'symlink'
    publishDir "${params.outdir}/logs/${test_name}", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_isoseq_${isoseq_mode}_${it}" }
    // Benchmark assembler: a broken tool (e.g. FLAMES container, missing dep) must
    // NOT abort the whole comparison run. Ignore failures so the other tools + FLAIR/TED proceed.
    errorStrategy 'ignore'
    tag "${dataset_name}_${align_mode}_${partition_mode}_isoseq_${isoseq_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path(bam), path(bai), val(genome), path(gtf),
          val(isoseq_mode), val(isoseq_args)

    output:
    // transcript_counts is included in the main isoseq_gff tuple (alongside
    // gff + read_map) so the evaluation subworkflow can forward it to
    // ted.py's --counts arg.  `isoseq collapse` does not carry forward
    // zero-support reference transcripts, so the filter is mostly a no-op,
    // but wiring it keeps behavior symmetric with the other assemblers and
    // makes --min-support > 1 meaningful if ever used.
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(isoseq_mode),
          path("${dataset_name}_${align_mode}_${partition_mode}_isoseq_${isoseq_mode}.gff"),
          path("${dataset_name}_${align_mode}_${partition_mode}_isoseq_${isoseq_mode}_read_map.txt"),
          path("${dataset_name}_${align_mode}_${partition_mode}_isoseq_${isoseq_mode}_transcript_counts.tsv", optional: true),
          emit: isoseq_gff
    path "${dataset_name}_${align_mode}_${partition_mode}_isoseq_${isoseq_mode}_transcript_counts.tsv", optional: true, emit: transcript_counts

    script:
    def output_prefix = "${dataset_name}_${align_mode}_${partition_mode}_isoseq_${isoseq_mode}"
    """
    # IsoSeq collapse: aligned BAM → collapsed isoforms GFF + read stat
    isoseq collapse \\
        ${bam} \\
        ${output_prefix}.gff \\
        --do-not-collapse-extra-5exons \\
        -j ${task.cpus} \\
        ${isoseq_args}

    # Convert the read_stat.txt to FLAIR-format read map
    python ${projectDir}/bin/convert_read_map.py \\
        --isoseq-read-stat ${output_prefix}.read_stat.txt \\
        --output ${output_prefix}_read_map.txt \\
        --verbose

    # Per-transcript counts in IsoQuant's `#feature_id<TAB>count` convention
    # (consumed by ted.py --counts).  read_stat.txt rows already represent
    # one read assignment each, so a pbid count is just `wc -l` per pbid.
    if [ -s "${output_prefix}.read_stat.txt" ]; then
        python - <<PYEOF
import csv, sys
counts = {}
src = "${output_prefix}.read_stat.txt"
with open(src) as fh:
    reader = csv.DictReader(fh, delimiter='\\t')
    for row in reader:
        pbid = (row.get('pbid') or '').strip()
        if not pbid or pbid in ('NA', '*'):
            continue
        counts[pbid] = counts.get(pbid, 0) + 1
out = "${output_prefix}_transcript_counts.tsv"
with open(out, 'w') as ofh:
    ofh.write("#feature_id\\tcount\\n")
    for pbid in sorted(counts.keys()):
        ofh.write(f"{pbid}\\t{counts[pbid]}\\n")
print(f"IsoSeq counts: wrote {len(counts)} transcripts -> {out}", file=sys.stderr)
PYEOF
    else
        echo "WARNING: No IsoSeq read_stat.txt found; skipping counts TSV" >&2
    fi
    """
}
