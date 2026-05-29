// Module: WriteAssemblyManifest
// Writes the STAGE-1 assembly manifest CSV (header + one row per produced
// assembly) to a stable, copied location so STAGE 2 (EVALUATION) can consume
// it as --input. Rows are pre-formatted CSV strings emitted by the ASSEMBLE
// subworkflow; the header is injected here from Utils.manifestHeader().

process WriteAssemblyManifest {
    publishDir "${params.outdir}/assemblies/${test_name}", mode: 'copy'
    tag "${test_name}_manifest"

    input:
    tuple val(test_name), val(rows)

    output:
    tuple val(test_name), path("assembly_manifest.csv"), emit: manifest

    exec:
    def header = Utils.manifestHeader()
    def outFile = task.workDir.resolve("assembly_manifest.csv")
    outFile.text = ([header] + rows).join('\n') + '\n'
}
