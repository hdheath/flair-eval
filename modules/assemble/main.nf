/*
 * ASSEMBLE: flair transcriptome (>=3) or correct+collapse (2.x) through the adapter, inside the
 * run's FLAIR env. Writes the normalised files of CONTRACTS section 2 into assemble/.
 */
process ASSEMBLE {
    tag "${meta.id}"
    label 'process_high'

    conda     { flair_env.conda_prefix ?: null }
    container { flair_env.container ?: null }

    input:
    tuple val(meta), path(bam), path(bai), path(reads), path(genome), path(gtf), val(flair_env)
    path(regions)

    output:
    tuple val(meta), path('assemble', type: 'dir'), emit: dir

    script:
    def prefix = flair_env.conda_prefix
    def args   = meta.run.transcriptome_args ?: ''
    """
    ${prefix ? "export PATH=\"${prefix}/bin:\$PATH\"" : ''}
    flair_adapter.py assemble \\
        --bam ${bam} \\
        ${reads ? "--reads ${reads}" : ''} \\
        --genome ${genome} \\
        --gtf ${gtf} \\
        --args '${args.replace("'", "'\\''")}' \\
        --threads ${task.cpus} \\
        ${regions ? "--regions ${regions}" : ''} \\
        --out assemble

    for f in isoforms.bed isoforms.gtf isoforms.fa; do
        if [ ! -s "assemble/\$f" ]; then
            echo "ASSEMBLE: the adapter did not write assemble/\$f (see CONTRACTS section 2)" >&2
            exit 1
        fi
    done
    """

    stub:
    """
    flair_adapter.py --help > /dev/null
    mkdir -p assemble
    printf 'chr1\\t100\\t900\\tiso1_stub\\t1000\\t+\\t100\\t900\\t0\\t2\\t200,300,\\t0,500,\\n' > assemble/isoforms.bed
    printf 'chr1\\tFLAIR\\ttranscript\\t101\\t900\\t.\\t+\\t.\\tgene_id "novel_gene_1"; transcript_id "iso1_stub";\\n' > assemble/isoforms.gtf
    printf '>iso1_stub\\nACGTACGTAC\\n' > assemble/isoforms.fa
    printf 'read1\\tiso1_stub\\n' > assemble/read_map.tsv
    cat > assemble/adapter_log.json <<'JSON'
${groovy.json.JsonOutput.toJson([step: 'assemble', stub: true, commands: [], args: meta.run.transcriptome_args])}
JSON
    """
}
