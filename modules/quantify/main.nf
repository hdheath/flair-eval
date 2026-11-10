/*
 * QUANTIFY: flair quantify through the adapter (CONTRACTS section 2). The adapter may write
 * skipped.json instead of counts.tsv when the FLAIR version cannot quantify.
 */
process QUANTIFY {
    tag "${meta.id}"
    label 'process_medium'

    conda     { flair_env.conda_prefix ?: null }
    container { flair_env.container ?: null }

    input:
    tuple val(meta), path(assemble), path(reads), val(flair_env)

    output:
    tuple val(meta), path('quantify', type: 'dir'), emit: dir

    script:
    def prefix = flair_env.conda_prefix
    """
    ${prefix ? "export PATH=\"${prefix}/bin:\$PATH\"" : ''}
    flair_adapter.py quantify \\
        --isoforms-fa ${assemble}/isoforms.fa \\
        --reads ${reads} \\
        --sample ${meta.sample} \\
        --threads ${task.cpus} \\
        --out quantify
    """

    stub:
    """
    flair_adapter.py --help > /dev/null
    mkdir -p quantify
    printf 'isoform_id\\tcount\\tcpm\\niso1_stub\\t1\\t1000000\\n' > quantify/counts.tsv
    printf '{"step":"quantify","stub":true,"commands":[]}\\n' > quantify/adapter_log.json
    """
}
