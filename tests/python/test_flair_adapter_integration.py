"""Integration test: run every adapter step with real FLAIR envs on a tiny data set.

Opt-in (skipped unless FLAIR_EVAL_IT_ENVS is set). Run it on a compute node, e.g.
  FLAIR_EVAL_IT_ENVS=/envs/flair-2.2.0,/envs/flair-3.0.1,/envs/flair-dev \
  FLAIR_EVAL_IT_READS=reads.fa FLAIR_EVAL_IT_GENOME=genome.fa FLAIR_EVAL_IT_GTF=annotation.gtf \
  FLAIR_EVAL_IT_OUT=/scratch/dir [FLAIR_EVAL_IT_REGIONS=regions.bed] [FLAIR_EVAL_IT_ALIGN_ENV=/envs/flair-3.0.1] \
  python -m pytest tests/python/test_flair_adapter_integration.py -v
The shared BAM is aligned once with FLAIR_EVAL_IT_ALIGN_ENV (default: the first env that has
`flair align`), like the pipeline's shared align step. docs/flair_versions.md lists the data used.
"""
import json
import os
import subprocess

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ADAPTER = os.path.join(ROOT, "bin", "flair_adapter.py")
ENVS = [e for e in os.environ.get("FLAIR_EVAL_IT_ENVS", "").split(",") if e]
pytestmark = pytest.mark.skipif(not ENVS, reason="set FLAIR_EVAL_IT_ENVS to run the integration test")
THREADS = os.environ.get("FLAIR_EVAL_IT_THREADS", "4")


def cfg(name, required=True):
    v = os.environ.get("FLAIR_EVAL_IT_" + name)
    if required and not v:
        pytest.skip("set FLAIR_EVAL_IT_%s" % name)
    return v


def adapter(env_prefix, *args):
    env = dict(os.environ)
    env["PATH"] = os.path.join(env_prefix, "bin") + os.pathsep + env["PATH"]
    for k in ("PYTHONPATH", "PYTHONHOME"):
        env.pop(k, None)
    py = os.path.join(env_prefix, "bin", "python")
    p = subprocess.run([py, ADAPTER] + [str(a) for a in args], env=env, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, universal_newlines=True)
    return p


def env_id(prefix):
    return os.path.basename(os.path.normpath(prefix))


def tsv(path):
    with open(path) as fh:
        return [l.rstrip("\n").split("\t") for l in fh if l.strip()]


def log_of(d):
    return json.load(open(os.path.join(d, "adapter_log.json")))


@pytest.fixture(scope="module")
def shared_bam():
    out = os.path.join(cfg("OUT"), "shared_align")
    align_env = os.environ.get("FLAIR_EVAL_IT_ALIGN_ENV") or ENVS[0]
    p = adapter(align_env, "align", "--reads", cfg("READS"), "--genome", cfg("GENOME"),
                "--threads", THREADS, "--out", out)
    assert p.returncode == 0, p.stderr
    if os.path.exists(os.path.join(out, "skipped.json")):
        pytest.skip("FLAIR_EVAL_IT_ALIGN_ENV has no align; set it to a version that does")
    return os.path.join(out, "aligned.bam")


def check_assemble(d):
    bed = tsv(os.path.join(d, "isoforms.bed"))
    assert bed, "no isoforms"
    assert all(len(f) == 12 for f in bed)
    names = [f[3] for f in bed]
    assert len(set(names)) == len(names)
    for f in bed:
        n, sizes, starts = int(f[9]), f[10].split(","), f[11].split(",")
        assert len(sizes) == n == len(starts)
        assert int(f[1]) + int(starts[-1]) + int(sizes[-1]) == int(f[2])
    tx = [f for f in tsv(os.path.join(d, "isoforms.gtf")) if f[2] == "transcript"]
    gtf_ids = {f[8].split('transcript_id "')[1].split('"')[0] for f in tx}
    assert gtf_ids == set(names)
    assert all('gene_id "' in f[8] for f in tx)
    fa = [l[1:].strip() for l in open(os.path.join(d, "isoforms.fa")) if l.startswith(">")]
    assert fa == names
    if os.path.exists(os.path.join(d, "read_map.tsv")):
        rm = tsv(os.path.join(d, "read_map.tsv"))
        assert rm[0] == ["read_id", "isoform_id"]
        assert {r[1] for r in rm[1:]} <= set(names)
    assert log_of(d)["status"] == "ok"
    return names


@pytest.mark.parametrize("prefix", ENVS, ids=[env_id(e) for e in ENVS])
def test_version_end_to_end(prefix, shared_bam):
    base = os.path.join(cfg("OUT"), env_id(prefix))
    reads, genome, gtf = cfg("READS"), cfg("GENOME"), cfg("GTF")
    results = {}

    d = os.path.join(base, "align")
    p = adapter(prefix, "align", "--reads", reads, "--genome", genome, "--threads", THREADS, "--out", d)
    assert p.returncode == 0, p.stderr
    if os.path.exists(os.path.join(d, "skipped.json")):
        results["align"] = "skipped"
    else:
        assert os.path.getsize(os.path.join(d, "aligned.bam")) > 0
        assert os.path.exists(os.path.join(d, "aligned.bam.bai"))
        results["align"] = "ok"

    d = os.path.join(base, "assemble")
    p = adapter(prefix, "assemble", "--bam", shared_bam, "--reads", reads, "--genome", genome, "--gtf", gtf,
                "--args", "", "--threads", THREADS, "--out", d)
    assert p.returncode == 0, p.stderr
    names = check_assemble(d)
    results["assemble"] = len(names)

    # A6: BAM-only sample (no --reads); 2.x extracts the reads from the BAM for collapse
    db = os.path.join(base, "assemble_bam_only")
    p = adapter(prefix, "assemble", "--bam", shared_bam, "--genome", genome, "--gtf", gtf,
                "--args", "", "--threads", THREADS, "--out", db)
    assert p.returncode == 0, p.stderr
    results["assemble_bam_only"] = len(check_assemble(db))

    regions = cfg("REGIONS", required=False)
    if regions:
        dr = os.path.join(base, "assemble_regions")
        p = adapter(prefix, "assemble", "--bam", shared_bam, "--reads", reads, "--genome", genome,
                    "--gtf", gtf, "--args", "", "--threads", THREADS, "--regions", regions, "--out", dr)
        assert p.returncode == 0, p.stderr
        results["assemble_regions"] = len(check_assemble(dr))

    # contract call (reads only): dev quantify needs a BAM, so it must skip cleanly
    dq = os.path.join(base, "quantify_contract")
    p = adapter(prefix, "quantify", "--isoforms-fa", os.path.join(d, "isoforms.fa"), "--reads", reads,
                "--sample", "it_sample", "--threads", THREADS, "--out", dq)
    assert p.returncode == 0, p.stderr
    results["quantify_contract"] = "skipped" if os.path.exists(os.path.join(dq, "skipped.json")) else "ok"

    # A7/A8: everything core has, FLAIR's own BED via --flair-bed (dev quantifies from the BAM)
    full = ["--sample", "it_sample", "--bam", shared_bam, "--genome", genome,
            "--isoforms-bed", os.path.join(d, "isoforms.bed"), "--flair-bed", os.path.join(d, "isoforms.flair.bed"),
            "--threads", THREADS]
    dq = os.path.join(base, "quantify")
    p = adapter(prefix, "quantify", "--isoforms-fa", os.path.join(d, "isoforms.fa"), "--reads", reads,
                *full, "--out", dq)
    assert p.returncode == 0, p.stderr
    assert not os.path.exists(os.path.join(dq, "skipped.json")), open(os.path.join(dq, "skipped.json")).read()
    rows = tsv(os.path.join(dq, "counts.tsv"))
    assert rows[0] == ["isoform_id", "count", "cpm"]
    assert [r[0] for r in rows[1:]] == names
    total = sum(float(r[1]) for r in rows[1:])
    assert total > 0
    assert abs(sum(float(r[2]) for r in rows[1:]) - 1e6) < 1
    results["quantify"] = total
    results["quantify_mode"] = log_of(dq)["profile"]["quantify_mode"]

    # BAM-only sample: only versions that count from the BAM (dev) produce counts
    dqb = os.path.join(base, "quantify_bam_only")
    p = adapter(prefix, "quantify", "--isoforms-fa", os.path.join(d, "isoforms.fa"), *full, "--out", dqb)
    assert p.returncode == 0, p.stderr
    skipped = os.path.exists(os.path.join(dqb, "skipped.json"))
    assert skipped == (results["quantify_mode"] == "reads_manifest")
    results["quantify_bam_only"] = "skipped" if skipped else "ok"

    dp = os.path.join(base, "productivity")
    p = adapter(prefix, "productivity", "--isoforms-bed", os.path.join(d, "isoforms.bed"), "--genome", genome,
                "--gtf", gtf, "--flair-bed", os.path.join(d, "isoforms.flair.bed"), "--out", dp)
    assert p.returncode == 0, p.stderr
    if os.path.exists(os.path.join(dp, "skipped.json")):
        results["productivity"] = "skipped"
    else:
        rows = tsv(os.path.join(dp, "productivity.tsv"))
        assert rows[0] == ["isoform_id", "cds_start", "cds_end", "class"]
        assert [r[0] for r in rows[1:]] == names
        assert {r[3] for r in rows[1:]} <= {"PRO", "PTC", "NST", "NGO"}
        results["productivity"] = log_of(dp)["outputs"]["method"]

    results["flair_version"] = log_of(d)["flair_version"]
    with open(os.path.join(base, "integration_summary.json"), "w") as fh:
        json.dump(results, fh, indent=2)
    print(env_id(prefix), json.dumps(results))
