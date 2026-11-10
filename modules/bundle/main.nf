/*
 * BUNDLE: assemble one job's results folder, <outdir>/bundles/<sample>/<run_id>/ (CONTRACTS
 * section 3), and fill bundle.json 'compute', 'pipeline' and 'inputs'.
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
    flair_eval_bundle.py --info info.json --env-json ${env_json} --steps in --out ${meta.run_id}
    """

    stub:
    """
    cat > info.json <<'JSON'
${groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info))}
JSON
    flair_eval_bundle.py --info info.json --env-json ${env_json} --steps in --out ${meta.run_id}
    """
}
