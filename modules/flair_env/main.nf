// FLAIR_ENV (adapter stream, docs/CONTRACTS.md section 2): resolve a run's `flair:` source to a
// conda prefix or container and write flair_env.json. Signature is part of the contract.
//
// - release/git sources are built once into params.env_cache (keyed release-<v> / git-<sha12>);
//   concurrent tasks for the same key wait on a file lock there, so this process needs no maxForks.
// - A first build downloads packages (PyPI/conda/GitHub) and takes ~5-15 min; a cache hit or an
//   env/container source takes seconds. The task needs network, `git`, and mamba/conda on PATH
//   (or FLAIR_EVAL_CONDA). It must NOT run inside a FLAIR env.
// - Core sizes it via the 'flair_env' label (1 CPU, ~4 GB is plenty; allow ~1 h for a cold build).
process FLAIR_ENV {
    tag "${run.id}"
    label 'flair_env'

    input:
    val(run)

    output:
    tuple val(run.id), path("flair_env.json")

    script:
    def run_json = groovy.json.JsonOutput.toJson(run)
    """
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 XGBOOST_NTHREAD=1
    cat > run.json <<'JSON'
${run_json}
JSON
    flair_env.py resolve --run-json run.json --cache '${params.env_cache}' --out flair_env.json
    """

    stub:
    // CONTRACTS A11: `source` is the run's flair key (release|git|env|container), never "stub".
    def kinds     = ['release', 'git', 'env', 'container']
    def src       = kinds.find { run.flair instanceof Map && run.flair.containsKey(it) } ?: 'env'
    def val       = run.flair instanceof Map ? (run.flair[src] ?: '') : ''
    def requested = (src == 'git' ? "${val}@${run.flair.ref ?: ''}" : "${val}").toString()
    def stub_json = groovy.json.JsonOutput.toJson([
        run_id           : run.id.toString(),
        source           : src,
        requested        : requested,
        resolved_commit  : null,
        flair_version    : 'FLAIR stub',
        conda_prefix     : null,     // stub tasks must not activate a real env
        container        : null,
        python           : null,
        pip_freeze_sha256: null,
        created          : '1970-01-01T00:00:00+00:00'
    ])
    """
    cat > flair_env.json <<'JSON'
${stub_json}
JSON
    """
}
