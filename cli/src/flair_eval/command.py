"""Build the `nextflow run hdheath/flair-eval ...` command line.

`build_command()` picks the one-command form (`--reads ... --flair 3.0.1,dev`) when the
samples and runs can be expressed that way exactly, and otherwise the file form
(`--input samples.csv --runs runs.yaml`). frontend/builder/logic.js mirrors this module;
the tests feed both the same cases and compare the output.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from . import DEFAULT_REF, PIPELINE
from .config import TRACK_COLUMNS, run_is_default
from .tokens import token_for_run

#: Pipeline params (CONTRACTS section 1) and their pipeline defaults. A param is only put on
#: the command line when it differs from this default. `outdir` is always written.
PARAM_DEFAULTS = {
    "outdir": "./flair-eval-results",
    "report_title": "FLAIR comparison",
    "baseline": None,
    "align_flair": "3.0.1",
    "regions": None,
    "bundles": None,
    "env_cache": None,  # pipeline default is "${launchDir}/flair-env-cache"
}
PARAM_ORDER = ["bundles", "baseline", "report_title", "align_flair", "regions", "env_cache"]

_SAFE = re.compile(r"^[A-Za-z0-9_@%+=:,./-]+$")


def _sq(s: str) -> str:
    if s and _SAFE.fullmatch(s):
        return s
    return "'" + s.replace("'", "'\"'\"'") + "'"


def quote(arg: str) -> str:
    """POSIX shell quoting, same rule as logic.js: bare when safe, else single quotes.
    `--name=value` keeps the name bare and quotes only the value."""
    s = str(arg)
    m = re.fullmatch(r"(--[A-Za-z0-9_-]+=)(.*)", s, re.S)
    if m:
        return m.group(1) + _sq(m.group(2))
    return _sq(s)


def opt(name: str, value) -> list:
    """argv for one `--name value` pipeline param. Nextflow reads a value that starts with
    `-` as another option (`--flair_args --stringent` sets flair_args=true), so such values
    are passed as `--name=value`."""
    v = str(value)
    return [f"{name}={v}"] if v.startswith("-") else [name, v]


@dataclass
class Options:
    ref: str = DEFAULT_REF
    profile: str | None = None
    outdir: str = PARAM_DEFAULTS["outdir"]
    params: dict = field(default_factory=dict)  # other pipeline params from PARAM_ORDER
    resume: bool = False
    work_dir: str | None = None
    samples_path: str = "samples.csv"
    runs_path: str = "runs.yaml"
    prefer_files: bool = False
    extra: list = field(default_factory=list)  # passed through verbatim


@dataclass
class Command:
    argv: list
    mode: str  # "quick" or "files"
    reasons: list  # why the one-command form was not used (empty in quick mode)

    def pretty(self) -> str:
        return format_command(self.argv)

    def __str__(self) -> str:
        return " ".join(quote(a) for a in self.argv)


def quick_form(samples: list[dict], runs: list[dict]) -> tuple[list | None, list[str]]:
    """Return (args, []) for the one-command form, or (None, reasons) if it cannot express
    these samples and runs exactly."""
    reasons: list[str] = []
    if len(samples) != 1:
        reasons.append(f"{len(samples)} samples (one-command mode takes exactly one)")
    else:
        s = samples[0]
        used = [c for c in TRACK_COLUMNS if s.get(c)]
        if used:
            reasons.append(f"extra tracks ({', '.join(used)}) need samples.csv")
        if s.get("reads") and s.get("bam"):
            reasons.append("the sample has both reads and bam")
    tokens = []
    for r in runs:
        rid = r.get("id", "?")
        tok = token_for_run(r)
        if tok is None:
            reasons.append(f"run {rid}: its id or FLAIR source has no --flair token")
        tokens.append(tok)
        if r.get("label") and r.get("label") != r.get("id"):
            reasons.append(f"run {rid}: labels need runs.yaml")
        if r.get("tags"):
            reasons.append(f"run {rid}: tags need runs.yaml")
        for k in ("quantify", "productivity", "own_align"):
            if not run_is_default(r, k):
                reasons.append(f"run {rid}: {k} is not the default")
    arg_sets = {r.get("transcriptome_args", "") or "" for r in runs}
    if len(arg_sets) > 1:
        reasons.append("runs use different transcriptome_args")
    if not runs:
        reasons.append("no runs")
    if reasons:
        return None, reasons

    s = samples[0]
    args = []
    if s.get("reads"):
        args += ["--reads", s["reads"]]
    else:
        args += ["--bam", s["bam"]]
    args += ["--genome", s["genome"], "--gtf", s["gtf"], "--sample", s["sample"]]
    for k in ("platform", "library"):
        if s.get(k):
            args += [f"--{k}", s[k]]
    args += ["--flair", ",".join(tokens)]
    fa = arg_sets.pop()
    if fa:
        args += opt("--flair_args", fa)
    return args, []


def _tail(opts: Options) -> list:
    args = []
    for k in PARAM_ORDER:
        v = opts.params.get(k)
        if v is None or v == "" or v == PARAM_DEFAULTS.get(k):
            continue
        args += opt(f"--{k}", "true" if v is True else v)
    args += ["--outdir", opts.outdir or PARAM_DEFAULTS["outdir"]]
    if opts.work_dir:
        args += ["-work-dir", opts.work_dir]
    if opts.resume:
        args.append("-resume")
    return args + list(opts.extra)


def _head(opts: Options) -> list:
    argv = ["nextflow", "run", PIPELINE, "-r", opts.ref or DEFAULT_REF]
    if opts.profile:
        argv += ["-profile", opts.profile]
    return argv


def build_command(samples: list[dict], runs: list[dict], opts: Options | None = None) -> Command:
    opts = opts or Options()
    quick, reasons = quick_form(samples, runs)
    if quick is not None and not opts.prefer_files:
        return Command(_head(opts) + quick + _tail(opts), "quick", [])
    if opts.prefer_files and quick is not None:
        reasons = ["files requested"]
    argv = _head(opts) + ["--input", opts.samples_path, "--runs", opts.runs_path] + _tail(opts)
    return Command(argv, "files", reasons)


def files_command(opts: Options | None = None) -> Command:
    opts = opts or Options()
    argv = _head(opts) + ["--input", opts.samples_path, "--runs", opts.runs_path] + _tail(opts)
    return Command(argv, "files", [])


def quick_command(sample_args: list, flair: str, flair_args: str = "",
                  opts: Options | None = None) -> Command:
    """One-command form straight from flags (used by `flair-eval run --reads ...`)."""
    opts = opts or Options()
    args = list(sample_args) + ["--flair", flair]
    if flair_args:
        args += opt("--flair_args", flair_args)
    return Command(_head(opts) + args + _tail(opts), "quick", [])


def report_command(bundles: list[str], opts: Options | None = None) -> Command:
    opts = opts or Options()
    params = dict(opts.params)
    params["bundles"] = ",".join(bundles)
    o = Options(**{**opts.__dict__, "params": params})
    argv = _head(o) + ["--report_only", "true"] + _tail(o)
    return Command(argv, "report", [])


NO_VALUE = {"-resume"}
VALUE_FLAGS = {"-r", "-profile", "-work-dir"}


def format_command(argv: list) -> str:
    """Multi-line display: `nextflow run X -r REF` on the first line, then one option
    (with its value) per line, joined with ` \\`. Same layout as logic.js."""
    q = [quote(a) for a in argv]
    lines = [" ".join(q[:5])]
    i = 5
    while i < len(argv):
        a = argv[i]
        takes_value = a not in NO_VALUE and "=" not in a and (a.startswith("--")
                                                             or a in VALUE_FLAGS)
        if takes_value and i + 1 < len(argv):
            lines.append("  " + q[i] + " " + q[i + 1])
            i += 2
        else:
            lines.append("  " + q[i])
            i += 1
    return " \\\n".join(lines)


def split_extra(argv: list) -> tuple[list, list]:
    """Split ['a', '--', 'b'] into (['a'], ['b'])."""
    if "--" in argv:
        i = argv.index("--")
        return argv[:i], argv[i + 1:]
    return argv, []
