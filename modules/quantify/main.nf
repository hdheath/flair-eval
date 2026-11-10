/*
 * QUANTIFY: flair quantify through the adapter (CONTRACTS section 2 and section 7 A6-A9).
 * Core passes everything it has: reads (when the sample has them), the job's BAM, the genome,
 * isoforms.fa/.bed and FLAIR's native BED (A8). The adapter picks what the version needs and
 * writes skipped.json instead of counts.tsv when something required is missing.
 */
include { adapterOptFlags } from '../adapter_opts'

process QUANTIFY {
    tag "${meta.id}"
    label 'process_medium'

    conda     { flair_env.conda_prefix ?: null }
    container { flair_env.container ?: null }

    input:
    tuple val(meta), path(assemble), path(reads), path(bam), path(bai), path(genome), val(flair_env)

    output:
    tuple val(meta), path('quantify', type: 'dir'), emit: dir

    script:
    def prefix = flair_env.conda_prefix
    def opts = adapterOptFlags('quantify', [
        ['--isoforms-bed', "${assemble}/isoforms.bed", 'file'],
        ['--flair-bed', "${assemble}/isoforms.flair.bed", 'file'],
        ['--bam', bam ? bam.toString() : null],
        ['--genome', genome.toString()],
        ['--platform', meta.platform],
        ['--library', meta.library],
    ])
    """
    ${prefix ? "export PATH=\"${prefix}/bin:\$PATH\"" : ''}
    ${opts}
    flair_adapter.py quantify \\
        --isoforms-fa ${assemble}/isoforms.fa \\
        ${reads ? "--reads ${reads}" : ''} \\
        --sample ${meta.sample} \\
        --threads ${task.cpus} \\
        --out quantify \\
        \${OPT[@]+"\${OPT[@]}"}
    """

    stub:
    """
    flair_adapter.py --help > /dev/null
    mkdir -p quantify
    printf 'isoform_id\\tcount\\tcpm\\niso1_stub\\t3\\t750000\\niso2_stub\\t1\\t250000\\n' > quantify/counts.tsv
    printf '{"step":"quantify","stub":true,"commands":[]}\\n' > quantify/adapter_log.json
    """
}
