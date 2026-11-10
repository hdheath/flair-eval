"""The front ends must only write params the pipeline has, and must cover all of them.

Reads the pipeline's own nextflow_schema.json (skipped when the CLI is tested outside the
repository). For the CLI and for the builder page (logic.js under Node):
- every `--param` they can emit exists in nextflow_schema.json with a compatible type,
  the same enum values, and the same default where the front end leaves defaults off;
- every schema param is emitted, or listed in that front end's NOT_EXPOSED with a reason,
  or is hidden in the schema.
"""
import json
import re

import pytest

from flair_eval import cli, command
from flair_eval.command import PARAM_DEFAULTS, Options, build_command, quick_command, \
    report_command

from .test_builder_parity import COMMANDS, run_js


@pytest.fixture(scope="module")
def nf_params(repo):
    p = repo / "nextflow_schema.json"
    if not p.exists():
        pytest.skip("not in the pipeline repository")
    s = json.loads(p.read_text())
    out = {}
    for group in (s.get("$defs") or s.get("definitions") or {}).values():
        out.update(group.get("properties", {}))
    out.update(s.get("properties", {}))
    return out


def params_in(argv):
    """Pipeline params (--name or --name=value, not single-dash Nextflow options) with values."""
    got, i = {}, 0
    while i < len(argv):
        a = argv[i]
        m = re.fullmatch(r"--([A-Za-z0-9_]+)(?:=(.*))?", a, re.S)
        if m:
            if m.group(2) is not None:
                got.setdefault(m.group(1), []).append(m.group(2))
            else:
                got.setdefault(m.group(1), []).append(argv[i + 1])
                i += 1
        i += 1
    return got


ALL_PARAMS = {"bundles": "/old/bundles", "baseline": "flair-dev", "report_title": "T",
              "align_flair": "2.2.0", "regions": "/r.bed", "env_cache": "/cache"}
SAMPLE = {"sample": "s1", "reads": "/r.fq.gz", "genome": "/g.fa", "gtf": "/a.gtf",
          "platform": "ont", "library": "dRNA"}
RUNS = [{"id": "flair-3.0.1", "flair": {"release": "3.0.1"}, "transcriptome_args": "--stringent"},
        {"id": "flair-dev", "flair": {"git": "https://github.com/BrooksLabUCSC/flair",
                                       "ref": "dev"}, "transcriptome_args": "--stringent"}]


def cli_emitted():
    opts = Options(params=dict(ALL_PARAMS), outdir="/out")
    argvs = [
        build_command([SAMPLE], RUNS, opts).argv,                         # one-command form
        build_command([SAMPLE], RUNS, Options(params=dict(ALL_PARAMS), prefer_files=True)).argv,
        build_command([dict(SAMPLE, reads=None, bam="/b.bam")], RUNS, opts).argv,
        quick_command(["--bam", "/b.bam", "--genome", "g", "--gtf", "a", "--sample", "s"],
                      "3.0.1", "--x", opts).argv,
        report_command(["a", "b"], opts).argv,
    ]
    got = {}
    for argv in argvs:
        for k, v in params_in(argv).items():
            got.setdefault(k, []).extend(v)
    return got


def check_emitted(emitted, nf, who):
    problems = []
    for name, values in sorted(emitted.items()):
        spec = nf.get(name)
        if spec is None:
            problems.append(f"{who} writes --{name}, which nextflow_schema.json does not have")
            continue
        t = spec.get("type", "string")
        for v in values:
            if t == "boolean" and v not in ("true", "false"):
                problems.append(f"{who} writes --{name} {v!r}; the pipeline wants true/false")
            if t in ("integer", "number") and not re.fullmatch(r"-?[0-9.]+", v):
                problems.append(f"{who} writes --{name} {v!r}; the pipeline wants a number")
            if "enum" in spec and v not in spec["enum"]:
                problems.append(f"{who} writes --{name} {v!r}; allowed: {spec['enum']}")
    return problems


def test_cli_params_exist_with_compatible_types(nf_params):
    emitted = cli_emitted()
    assert check_emitted(emitted, nf_params, "the CLI") == []
    # the CLI's own flags must cover what it can emit
    parser_dests = {a.dest for p in cli.make_parser()._subparsers._group_actions[0].choices.values()
                    for a in p._actions}
    for name in emitted:
        assert name in parser_dests or name in ("input", "report_only", "flair_args"), name


def test_every_pipeline_param_is_covered_by_the_cli(nf_params):
    emitted = set(cli_emitted())
    missing = [n for n, spec in nf_params.items()
               if n not in emitted and n not in command.NOT_EXPOSED and not spec.get("hidden")]
    assert missing == [], f"add these to the CLI or to command.NOT_EXPOSED: {missing}"
    stale = [n for n in command.NOT_EXPOSED if n not in nf_params]
    assert stale == [], f"command.NOT_EXPOSED lists params the pipeline no longer has: {stale}"


def test_defaults_match(nf_params):
    """The front ends leave a param off when it equals the default, so defaults must agree."""
    wrong = {name: (ours, nf_params.get(name, {}).get("default", "<no such param>"))
             for name, ours in PARAM_DEFAULTS.items()
             if ours is not None and nf_params.get(name, {}).get("default") != ours}
    assert wrong == {}, f"(front-end default, pipeline default): {wrong}"


def test_choices_match(nf_params, repo):
    from flair_eval.validate import schema
    sub = cli.make_parser()._subparsers._group_actions[0].choices["run"]
    choices = {a.dest: list(a.choices) for a in sub._actions if a.choices}
    for name in ("platform", "library"):
        assert choices[name] == nf_params[name]["enum"], name
        assert schema("samples")["properties"][name]["enum"] == nf_params[name]["enum"], name


@pytest.fixture(scope="module")
def js_out(repo):
    out = run_js(repo)
    if out is None:
        pytest.skip("node or frontend/builder not available")
    return out


def test_builder_params_exist_with_compatible_types(nf_params, js_out):
    emitted = {}
    for c in js_out["commands"]:
        for k, v in params_in(c["argv"]).items():
            emitted.setdefault(k, []).extend(v)
    assert check_emitted(emitted, nf_params, "the builder") == []
    meta = js_out["meta"]
    assert meta["param_defaults"] == PARAM_DEFAULTS
    assert meta["param_order"] == command.PARAM_ORDER


def test_every_pipeline_param_is_covered_by_the_builder(nf_params, js_out, repo):
    emitted = set()
    for c in js_out["commands"]:
        emitted |= set(params_in(c["argv"]))
    meta = js_out["meta"]
    missing = [n for n, spec in nf_params.items()
               if n not in emitted and n not in meta["not_exposed"] and not spec.get("hidden")]
    assert missing == [], f"add these to the builder or to NOT_EXPOSED in logic.js: {missing}"
    # every option field on the page is a pipeline param or a declared Nextflow launch option
    html = (repo / "frontend" / "builder" / "index.html").read_text()
    for opt in re.findall(r'data-opt="([a-z_]+)"', html):
        assert opt in nf_params or opt in meta["launch_options"], opt


def test_parity_cases_set_every_param():
    """COMMANDS (shared with the parity test) must exercise every param the page can set."""
    set_somewhere = set()
    for c in COMMANDS:
        set_somewhere |= {k for k, v in c["opts"].get("params", {}).items() if v}
    assert set(command.PARAM_ORDER) <= set_somewhere
