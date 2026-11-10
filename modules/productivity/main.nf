/*
 * PRODUCTIVITY: flair predictProductivity through the adapter (CONTRACTS section 2). The adapter
 * may write skipped.json instead of productivity.tsv.
 */
process PRODUCTIVITY {
    tag "${meta.id}"
    label 'process_single'

    conda     { flair_env.conda_prefix ?: null }
    container { flair_env.container ?: null }

    input:
    tuple val(meta), path(assemble), path(genome), path(gtf), val(flair_env)

    output:
    tuple val(meta), path('productivity', type: 'dir'), emit: dir

    script:
    def prefix = flair_env.conda_prefix
    """
    ${prefix ? "export PATH=\"${prefix}/bin:\$PATH\"" : ''}
    flair_adapter.py productivity \\
        --isoforms-bed ${assemble}/isoforms.bed \\
        --genome ${genome} \\
        --gtf ${gtf} \\
        --out productivity
    """

    stub:
    """
    flair_adapter.py --help > /dev/null
    mkdir -p productivity
    printf 'isoform_id\\tcds_start\\tcds_end\\tclass\\niso1_stub\\t150\\t800\\tPRO\\n' > productivity/productivity.tsv
    printf '{"step":"productivity","stub":true,"commands":[]}\\n' > productivity/adapter_log.json
    """
}
