/*
 * ALIGN: flair align through the adapter, inside a FLAIR env (CONTRACTS section 2).
 * Used twice: ALIGN_SHARED (once per sample, FLAIR = params.align_flair) and
 * ALIGN_OWN (per job, for runs with own_align: true, in that run's env).
 */
include { adapterOptFlags } from '../adapter_opts'

process ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda     { flair_env.conda_prefix ?: null }
    container { flair_env.container ?: null }

    input:
    tuple val(meta), path(reads), path(genome), path(gtf), val(flair_env)

    output:
    tuple val(meta), path('align/aligned.bam'), path('align/aligned.bam.bai'), emit: bam
    tuple val(meta), path('align', type: 'dir'),                               emit: dir

    script:
    def prefix = flair_env.conda_prefix
    def opts = adapterOptFlags('align', [['--platform', meta.platform], ['--library', meta.library]])
    """
    ${prefix ? "export PATH=\"${prefix}/bin:\$PATH\"" : ''}
    ${opts}
    flair_adapter.py align \\
        --reads ${reads} \\
        --genome ${genome} \\
        --gtf ${gtf} \\
        --threads ${task.cpus} \\
        --out align \\
        \${OPT[@]+"\${OPT[@]}"}
    """

    stub:
    """
    flair_adapter.py --help > /dev/null
    mkdir -p align
    touch align/aligned.bam align/aligned.bam.bai
    printf '{"step":"align","stub":true,"commands":[],"flair_version":"%s"}\\n' '${flair_env.flair_version}' > align/adapter_log.json
    """
}
