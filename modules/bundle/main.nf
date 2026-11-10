/*
 * BUNDLE: assemble one job's results folder, <outdir>/bundles/<sample>/<run_id>/ (CONTRACTS
 * section 3), and fill bundle.json 'compute', 'pipeline', 'inputs' and 'run_index' (A17).
 * sections.quant / sections.productivity come from `python3 -m flair_eval_lib.sections` in the
 * pipeline's bin/ (A19), plus sections.quant.input from the adapter log (A21).
 * Pure file shuffling, so the stub runs the same code as the real script.
 */
process BUNDLE {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(env_json, stageAs: 'in/flair_env.json'), path(step_dirs, stageAs: 'in/*'), val(info)

    output:
    tuple val(meta), path("${meta.run_id}", type: 'dir'), emit: bundle

    script:
    """
    cat > info.json <<'JSON'
${groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info))}
JSON
    flair_eval_bundle.py --info info.json --env-json ${env_json} --steps in --out ${meta.run_id} \\
        --lib-dir "\$(dirname "\$(command -v structure_metrics.py)")"
    """

    stub:
    """
    cat > info.json <<'JSON'
${groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info))}
JSON
    flair_eval_bundle.py --info info.json --env-json ${env_json} --steps in --out ${meta.run_id} \\
        --lib-dir "\$(dirname "\$(command -v structure_metrics.py)")"
    """
}
