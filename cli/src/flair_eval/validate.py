"""Validate samples.csv / runs.yaml / bundle.json against the bundled JSON schemas,
with plain-English messages.

Every problem is an `Issue(path, code, message, level)`. `path` is a tuple such as
("samples", 0, "genome") or ("runs", 1, "flair"). `code` is a short machine name; the
config-builder page (frontend/builder/logic.js) produces the same (path, code) pairs for
the same input, and the tests check that.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from functools import lru_cache
from importlib import resources
from pathlib import Path

from jsonschema import Draft202012Validator

from . import config

NAME_HINT = "start with a letter or digit, then use only letters, digits, '.', '_' or '-'"
FLAIR_KEYS = ("release", "git", "env", "container")


@dataclass(frozen=True, order=True)
class Issue:
    path: tuple
    code: str
    message: str = field(compare=False)
    level: str = field(default="error", compare=False)

    def where(self) -> str:
        return format_path(self.path)

    def __str__(self) -> str:
        w = self.where()
        return f"{self.level}: {w}: {self.message}" if w else f"{self.level}: {self.message}"


def format_path(path: tuple) -> str:
    if not path:
        return ""
    head, rest = path[0], list(path[1:])
    if head == "samples":
        out = "samples.csv"
        if rest and isinstance(rest[0], int):
            out += f" row {rest.pop(0) + 1}"
    elif head == "runs":
        out = "runs.yaml"
        if rest and isinstance(rest[0], int):
            out += f" run {rest.pop(0) + 1}"
    else:
        out = str(head)
    if rest:
        out += ": " + ".".join(str(p) for p in rest)
    return out


@lru_cache(maxsize=None)
def schema(name: str) -> dict:
    """Bundled schema: name is 'samples', 'runs' or 'bundle'."""
    text = resources.files("flair_eval").joinpath("schemas", f"{name}.schema.json").read_text()
    return json.loads(text)


def _type_name(v) -> str:
    if isinstance(v, bool):
        return "true/false"
    if isinstance(v, (int, float)):
        return "a number"
    if v is None:
        return "empty"
    if isinstance(v, list):
        return "a list"
    if isinstance(v, dict):
        return "a mapping"
    return "text"


def _issues_from_schema(instance, sch: dict, prefix: tuple, kind: str) -> list[Issue]:
    out: list[Issue] = []
    v = Draft202012Validator(sch)
    for e in v.iter_errors(instance):
        path = prefix + tuple(e.absolute_path)
        inst = e.instance
        kw = e.validator
        if kw == "required":
            m = re.match(r"^'(.+)' is a required property", e.message)
            fld = m.group(1) if m else "?"
            out.append(Issue(path + (fld,), "required", f"missing required {kind} '{fld}'"))
        elif kw == "anyOf" and isinstance(inst, dict) and kind == "column":
            out.append(Issue(path + ("reads",), "either_reads_bam",
                             "give reads (FASTQ/FASTA) or bam (sorted, indexed BAM)"))
        elif kw == "oneOf" and path[-1:] == ("flair",):
            out.append(_flair_issue(path, inst))
        elif kw == "additionalProperties":
            allowed = set(e.schema.get("properties", {}))
            for k in sorted(set(inst) - allowed, key=str):
                msg = f"unknown {kind} '{k}'"
                if kind == "key" and len(path) == 2 and path[0] == "runs":
                    msg += ("; allowed keys are id, flair, transcriptome_args, quantify, "
                            "productivity, own_align, label, tags")
                elif kind == "column":
                    msg += "; allowed columns are " + ", ".join(config.ALL_SAMPLE_COLUMNS)
                out.append(Issue(path + (k,), "unknown_key", msg))
        elif kw == "pattern":
            out.append(Issue(path, "pattern", f"{inst!r} is not a valid name: {NAME_HINT}"))
        elif kw == "enum":
            opts = ", ".join(str(x) for x in e.validator_value)
            out.append(Issue(path, "enum", f"{inst!r} is not allowed; use one of {opts}"))
        elif kw == "type":
            want = e.validator_value
            want = want if isinstance(want, str) else "/".join(want)
            nice = {"string": "text", "boolean": "true or false", "object": "a mapping",
                    "array": "a list", "integer": "a whole number", "number": "a number"}
            msg = f"should be {nice.get(want, want)}, not {_type_name(inst)}"
            if want == "string" and isinstance(inst, (int, float)) and not isinstance(inst, bool):
                msg += f' (put it in quotes, e.g. "{inst}")'
            out.append(Issue(path, "type", msg))
        elif kw == "minItems":
            out.append(Issue(path, "min_items", "needs at least one run"))
        elif kw == "minLength":
            out.append(Issue(path, "min_length", "must not be empty"))
        elif kw == "const":
            out.append(Issue(path, "const", f"should be {e.validator_value!r}"))
        else:  # pragma: no cover - schema keywords not used today
            out.append(Issue(path, kw, e.message))
    return out


def _flair_issue(path: tuple, inst) -> Issue:
    msg = "give exactly one FLAIR source: release, git (with ref), env or container"
    if isinstance(inst, dict):
        given = [k for k in FLAIR_KEYS if k in inst]
        extra = sorted(set(inst) - set(FLAIR_KEYS) - {"ref"})
        if len(given) > 1:
            msg = f"give only one FLAIR source, not {' and '.join(given)}"
        elif given == ["git"] and "ref" not in inst:
            msg = "git needs a ref: a branch, tag or commit SHA (e.g. ref: dev)"
        elif "ref" in inst and given != ["git"]:
            msg = "ref only goes with git"
        elif extra:
            msg = f"unknown key '{extra[0]}' in flair; {msg}"
        elif given and not (isinstance(inst[given[0]], str) and inst[given[0]].strip()):
            v = inst[given[0]]
            msg = f"{given[0]} should be non-empty text"
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                msg += f' (put it in quotes, e.g. {given[0]}: "{v}")'
        elif given == ["git"] and not (isinstance(inst["ref"], str) and inst["ref"].strip()):
            msg = "ref should be non-empty text"
            if isinstance(inst["ref"], (int, float)):
                msg += f' (put it in quotes, e.g. ref: "{inst["ref"]}")'
    return Issue(path, "flair_source", msg)


def _duplicates(items: list[dict], key: str, prefix: str, what: str) -> list[Issue]:
    out, first = [], {}
    for i, it in enumerate(items):
        v = it.get(key) if isinstance(it, dict) else None
        if not isinstance(v, str):
            continue
        if v in first:
            out.append(Issue((prefix, i, key), "duplicate",
                             f"{what} {v!r} is already used by {prefix[:-1]} {first[v] + 1}"))
        else:
            first[v] = i
    return out


# ---------------------------------------------------------------- public API
def validate_samples(rows: list[dict], header: list[str] | None = None) -> list[Issue]:
    """Validate sample rows (as dicts with empty cells removed)."""
    out: list[Issue] = []
    if header is not None:
        for c in header:
            if c and c not in config.ALL_SAMPLE_COLUMNS:
                out.append(Issue(("samples", "header", c), "unknown_key",
                                 f"unknown column '{c}'; allowed columns are "
                                 + ", ".join(config.ALL_SAMPLE_COLUMNS)))
    if not rows:
        out.append(Issue(("samples",), "min_items", "needs at least one sample row"))
    sch = schema("samples")
    for i, r in enumerate(rows):
        issues = _issues_from_schema(r, sch, ("samples", i), "column")
        if header is not None:  # reported once on the header instead
            issues = [x for x in issues if x.code != "unknown_key"]
        out += issues
    out += _duplicates(rows, "sample", "samples", "sample name")
    return out


def validate_runs(doc) -> list[Issue]:
    """Validate a parsed runs.yaml document ({"runs": [...]})."""
    if isinstance(doc, list):
        return [Issue(("runs",), "type",
                      "runs.yaml should start with a 'runs:' key that holds the list of runs")]
    if not isinstance(doc, dict):
        return [Issue(("runs",), "type", "runs.yaml should be a mapping with a 'runs:' list")]
    out = _issues_from_schema(doc, schema("runs"), (), "key")
    runs = doc.get("runs")
    if isinstance(runs, list):
        out += _duplicates(runs, "id", "runs", "run id")
    return out


def validate_bundle(doc) -> list[Issue]:
    return _issues_from_schema(doc, schema("bundle"), ("bundle",), "key")


def check_paths(rows: list[dict], base: Path | None = None) -> list[Issue]:
    """Warnings for local input files that do not exist (remote URIs are skipped)."""
    out: list[Issue] = []
    base = base or Path.cwd()
    for i, r in enumerate(rows):
        for col in ("reads", "bam", "genome", "gtf") + tuple(config.TRACK_COLUMNS):
            p = r.get(col)
            if not p or re.match(r"^[a-z][a-z0-9+.-]*://", p):
                continue
            fp = Path(p).expanduser()
            fp = fp if fp.is_absolute() else base / fp
            if any(ch in p for ch in "*?["):
                continue
            if not fp.exists():
                out.append(Issue(("samples", i, col), "missing_file",
                                 f"file not found: {p}", "warning"))
            elif col == "bam" and not (Path(str(fp) + ".bai").exists()
                                       or fp.with_suffix(".bai").exists()):
                out.append(Issue(("samples", i, col), "missing_index",
                                 f"no index next to {p} (run: samtools index {p})", "warning"))
    return out


def check_files(samples: str | Path | None = None, runs: str | Path | None = None,
                files: bool = True) -> list[Issue]:
    """Validate config files on disk. Missing local inputs are warnings."""
    out: list[Issue] = []
    if samples is not None:
        sp = Path(samples)
        if not sp.exists():
            out.append(Issue(("samples",), "missing_file", f"{sp} not found"))
        else:
            header, rows = config.read_samples(sp)
            out += validate_samples(rows, header)
            if files:
                out += check_paths(rows)  # Nextflow resolves relative paths from launchDir
    if runs is not None:
        rp = Path(runs)
        if not rp.exists():
            out.append(Issue(("runs",), "missing_file", f"{rp} not found"))
        else:
            try:
                out += validate_runs(config.read_runs(rp))
            except config.YamlError as e:
                out.append(Issue(("runs",), "yaml", str(e)))
    return out


def errors(issues: list[Issue]) -> list[Issue]:
    return [i for i in issues if i.level == "error"]
