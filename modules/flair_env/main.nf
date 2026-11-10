// PLACEHOLDER owned by the 'adapter' stream (docs/CONTRACTS.md section 2). Keep this signature.
process FLAIR_ENV {
    tag "${run.id}"

    input:
    val(run)

    output:
    tuple val(run.id), path("flair_env.json")

    script:
    def run_json = groovy.json.JsonOutput.toJson(run)
    """
    cat > run.json <<'JSON'
${run_json}
JSON
    flair_env.py resolve --run-json run.json --cache '${params.env_cache}' --out flair_env.json
    """

    stub:
    """
    printf '{"run_id":"%s","source":"stub","requested":"stub","resolved_commit":null,"flair_version":"stub","conda_prefix":null,"container":null,"python":null,"pip_freeze_sha256":null,"created":"1970-01-01T00:00:00Z"}\\n' '${run.id}' > flair_env.json
    """
}
