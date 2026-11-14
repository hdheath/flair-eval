"""frontend/builder/logic.js must agree with this package on every shared function.

Runs logic.js under Node on the same cases and compares. Skipped when Node or the builder
page is not available (e.g. the CLI installed on its own).
"""
import json
import shutil
import subprocess

import pytest

from flair_eval import config, presets
from flair_eval.command import Options, as_cli_run, build_command, format_command
from flair_eval.tokens import TokenError, parse_token, parse_tokens, token_for_flair
from flair_eval.validate import validate_runs, validate_samples

from .test_tokens import CONTRACT_TABLE

GIT = "https://github.com/BrooksLabUCSC/flair"
S1 = {"sample": "HG002", "reads": "/data/HG002.fq.gz", "genome": "/ref/hg38.fa",
      "gtf": "/ref/gencode.gtf"}
R1 = {"id": "flair-3.0.1", "flair": {"release": "3.0.1"}, "transcriptome_args": ""}
R2 = {"id": "flair-dev", "flair": {"git": GIT, "ref": "dev"}, "transcriptome_args": ""}

TOKENS = [t for t, _, _ in CONTRACT_TABLE] + [
    "git:git@github.com:me/flair.git@abc123", "git:/home/me/flair@dev", "git:feature/ends",
    "env:/x/myenv/", "container:docker://org/flair@sha256:abcd", "container:docker://flair:3",
    "  2.2.0 ", "release:v3", "", " ", "master", "release:", "git:", "env:", "container:",
    "git:https://github.com/me/flair", "git:git@github.com:me/flair.git", "env:/", "3.0.1rc1",
    "git:a@b", "git:https://x@host/r@v1", "env:relative/env", "container:image.sif",
    "release:3 0", "git:we ird", "Dev", "1",
]
TOKEN_LISTS = ["3.0.1,dev", "3.0.1,release:3.0.1", "3.0.1,,dev", "", "dev, main ,git:x"]
FLAIRS = [f for _, f, _ in CONTRACT_TABLE] + [
    {"release": 3.0}, {"release": "v3"}, {"release": "3.0.1", "env": "/x"}, {"env": "/a,b"},
    {}, {"git": GIT, "ref": "x"}, {"git": GIT, "ref": "feature/x"}, {"container": " x"},
    {"git": "git@github.com:me/f.git", "ref": "dev"}, {"ref": "dev"},
]
PATHS = ["/x/HG002.fq.gz", "HG002.rep1.fastq", "/a/b/s1.bam", "reads.FA.GZ", "s.fasta",
         "weird", ".fq.gz", "/dir/", "a.b.c.fq.bz2", "x.cram"]
SAMPLE_SETS = [
    {"rows": [S1]},
    {"rows": []},
    {"rows": [S1, dict(S1)]},
    {"rows": [{"sample": "s", "genome": "g", "gtf": "a"}]},
    {"rows": [{"sample": "a b", "reads": "r", "gtf": "a", "platform": "nanopore",
               "library": "cdna", "colour": "x"}]},
    {"rows": [dict(S1, extra="1")], "header": ["sample", "reads", "genome", "gtf", "extra"]},
    {"rows": [{"sample": "s", "bam": "x.bam", "genome": "g", "gtf": "a", "platform": "ont",
               "library": "dRNA", "cage": "c", "pas": "p", "junctions": "j",
               "spikein_gtf": "sirv.gtf"}]},
    {"rows": [{"sample": "s", "reads": "r", "genome": "hg38"}]},             # preset: gtf optional
    {"rows": [{"sample": "s", "reads": "r", "genome": "hg38", "gtf": "own.gtf"}]},
    {"rows": [{"sample": "s", "reads": "r", "genome": "HG38"}]},             # not a preset
]
RUN_DOCS = [
    {"runs": [R1, R2]}, {}, {"runs": []}, [R1], "x", {"runs": [], "extra": 1},
    {"runs": [{"id": "bad id", "flair": {"release": "3"}, "quantify": "yes", "tags": {"n": 3},
               "colour": 1}, {"flair": {"release": "3"}}, R1, dict(R1)]},
    {"runs": [{"id": "a", "flair": f} for f in [
        {"git": "x"}, {"release": "3.0.1", "env": "/x"}, {"ref": "dev"}, {"release": 3.0},
        {"release": ""}, {"git": "u", "ref": ""}, {"version": "3"}, {}, "3.0.1",
        {"git": "u", "ref": 5}, {"env": True}]]},
    {"runs": [dict(R1, label="L", tags={"a": "b"}, own_align=True, quantify=False,
                   productivity=False, transcriptome_args="--ted")]},
]


def cmd(samples, runs, **opts):
    return {"samples": samples, "runs": runs, "opts": opts}


COMMANDS = [
    cmd([S1], [R1, R2]),
    cmd([S1], [dict(R1, transcriptome_args="--stringent --check_splice"),
               dict(R2, transcriptome_args="--stringent --check_splice")],
        ref="2.0.0", profile="slurm", outdir="/big/out dir", resume=True, work_dir="/w",
        params={"baseline": "flair-dev", "report_title": "It's mine", "align_flair": "2.2.0",
                "regions": "/r.bed", "bundles": "/old/bundles,/older/bundles",
                "env_cache": "/cache"}),
    cmd([S1], [R1, R2], params={"align_flair": "3.0.1", "bundles": ""}),  # defaults left off
    cmd([{"sample": "s1", "bam": "/b/s1.bam", "genome": "g", "gtf": "a", "platform": "pacbio",
          "library": "cDNA"}], [R1]),
    cmd([S1, dict(S1, sample="s2", reads='/d/a,"b".fq')], [R1, R2]),
    cmd([dict(S1, cage="c.bed", junctions="SJ.out.tab")], [R1]),
    cmd([dict(S1, bam="x.bam")], [R1]),
    cmd([S1], [dict(R1, id="mine"), R2]),
    cmd([S1], [dict(R1, label="FLAIR 3"), dict(R2, label="flair-dev")]),
    cmd([S1], [dict(R1, tags={"version": "3.0.1", "k": 'q"uote'}), R2]),
    cmd([S1], [dict(R1, quantify=False, productivity=False, own_align=True), R2]),
    cmd([S1], [dict(R1, transcriptome_args="--stringent"), R2]),
    cmd([S1], []),
    cmd([S1], [R1, R2], prefer_files=True, samples_path="cfg/s.csv", runs_path="cfg/r.yaml"),
    cmd([S1], [{"id": "env-ted", "flair": {"env": "/envs/ted"}, "transcriptome_args": "--ted"},
               {"id": "ctr-flair", "flair": {"container": "docker://org/flair:3.0.1"},
                "transcriptome_args": "--ted"},
               {"id": "flair-x", "flair": {"git": "https://github.com/me/flair", "ref": "x"},
                "transcriptome_args": "--ted"}]),
    cmd([S1], [R1], params={"report_title": "-dash first"}),
    cmd([{"sample": "s1", "reads": "/d/s1.fq.gz", "genome": "hg38"}], [R1, R2],
        outdir="/private/groups/brookslab/me/flair-eval/10_11_26_x"),
    cmd([{"sample": "s1", "reads": "/d/s1.fq.gz", "genome": "hg38"}], [R1, R2],
        prefer_files=True),
]


def python_side():
    def tok(t):
        try:
            r = parse_token(t)
            return {"id": r.id, "flair": r.flair}
        except TokenError:
            return {"error": True}

    def toks(t):
        try:
            return [r.id for r in parse_tokens(t)]
        except TokenError:
            return {"error": True}

    def iss(issues):
        return [[list(i.path), i.code, i.message] for i in issues]

    commands = []
    for k in COMMANDS:
        c = build_command(k["samples"], k["runs"], Options(**k["opts"]))
        commands.append({"argv": c.argv, "mode": c.mode, "reasons": c.reasons,
                         "pretty": format_command(c.argv), "line": str(c),
                         "cli": format_command(as_cli_run(c.argv)),
                         "csv": config.samples_csv(k["samples"]),
                         "yaml": config.runs_yaml(k["runs"])})
    return {
        "tokens": [tok(t) for t in TOKENS],
        "token_lists": [toks(t) for t in TOKEN_LISTS],
        "token_for_flair": [token_for_flair(f) for f in FLAIRS],
        "sample_names": [config.default_sample_name(p) for p in PATHS],
        "validate_samples": [iss(validate_samples(s["rows"], s.get("header")))
                             for s in SAMPLE_SETS],
        "validate_runs": [iss(validate_runs(d)) for d in RUN_DOCS],
        "commands": commands,
        "presets": presets.PRESETS,
    }


def run_js(repo):
    """Run logic.js on all cases under Node; None when Node or the page is not available."""
    node = shutil.which("node") or shutil.which("nodejs")
    logic = repo / "frontend" / "builder" / "logic.js"
    if node is None or not logic.exists():
        return None
    cases = {"tokens": TOKENS, "token_lists": TOKEN_LISTS, "flairs": FLAIRS, "paths": PATHS,
             "sample_sets": SAMPLE_SETS, "run_docs": RUN_DOCS, "commands": COMMANDS}
    p = subprocess.run([node, str(repo / "cli" / "tests" / "js_runner.js"), str(logic),
                        str(repo / "frontend" / "builder" / "schemas.js")],
                       input=json.dumps(cases), capture_output=True, text=True, timeout=60)
    assert p.returncode == 0, p.stderr
    return json.loads(p.stdout)


@pytest.fixture(scope="module")
def js(repo):
    out = run_js(repo)
    if out is None:
        pytest.skip("node or frontend/builder not available")
    return out


@pytest.fixture(scope="module")
def py():
    return json.loads(json.dumps(python_side()))  # tuples -> lists, like the JS side


def _sort(issues):
    return sorted(issues, key=lambda i: ([str(x) for x in i[0]], i[1]))


@pytest.mark.parametrize("key", ["tokens", "token_lists", "token_for_flair", "sample_names",
                                 "presets"])
def test_same_results(js, py, key):
    for i, (a, b) in enumerate(zip(py[key], js[key])):
        assert a == b, f"{key}[{i}]"
    assert len(py[key]) == len(js[key])


@pytest.mark.parametrize("key", ["validate_samples", "validate_runs"])
def test_same_issues(js, py, key):
    for i, (a, b) in enumerate(zip(py[key], js[key])):
        pa = [x[:2] for x in _sort(a)]
        ja = [x[:2] for x in _sort(b)]
        assert pa == ja, f"{key}[{i}]"
        # Messages match too, except where Python repr() quoting differs from JS, and the
        # quoting hint for YAML numbers (JS has no float 3.0; the page never makes numbers).
        for (p, j) in zip(_sort(a), _sort(b)):
            if p[1] not in ("pattern", "enum", "duplicate", "const") \
                    and "put it in quotes" not in p[2]:
                assert p[2] == j[2], f"{key}[{i}] {p[0]}"


@pytest.mark.parametrize("i", range(len(COMMANDS)))
def test_same_command_and_files(js, py, i):
    assert js["commands"][i] == py["commands"][i]


def test_same_genome_presets_and_outdir_rules(js):
    from flair_eval import phoenix
    from flair_eval.command import OUTDIR_PLACEHOLDER
    meta = js["meta"]
    assert meta["genome_presets"] == phoenix.GENOME_PRESETS
    assert meta["outdir_placeholder"] == OUTDIR_PLACEHOLDER
    assert meta["outdir_example"] == phoenix.OUTDIR_EXAMPLE
