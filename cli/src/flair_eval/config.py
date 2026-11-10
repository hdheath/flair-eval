"""Read and write samples.csv and runs.yaml (docs/CONTRACTS.md section 1).

The writers emit a fixed, readable layout. frontend/builder/logic.js emits the same bytes,
and the tests check that the two agree.
"""

from __future__ import annotations

import csv
import io
import json
import re
from pathlib import Path

import yaml

#: Columns written to every samples.csv, in this order.
SAMPLE_COLUMNS = ["sample", "reads", "bam", "genome", "gtf", "platform", "library"]
#: Optional track columns, written only when some sample uses them.
TRACK_COLUMNS = ["cage", "pas", "junctions", "spikein_gtf"]
ALL_SAMPLE_COLUMNS = SAMPLE_COLUMNS + TRACK_COLUMNS

RUN_DEFAULTS = {"transcriptome_args": "", "quantify": True, "productivity": True,
                "own_align": False}

_READ_EXT = re.compile(r"(\.(fastq|fq|fasta|fa|fna|bam|sam|cram))?(\.(gz|bgz|bz2|xz|zst))?$",
                       re.IGNORECASE)


def default_sample_name(path: str) -> str:
    """Sample name used when --sample is not given: the reads/BAM basename without its
    read-file extensions (`/x/HG002.rep1.fq.gz` -> `HG002.rep1`)."""
    base = str(path).rstrip("/").rsplit("/", 1)[-1]
    name = _READ_EXT.sub("", base, count=1)
    return name or base


# ---------------------------------------------------------------- samples.csv
def clean_sample(row: dict) -> dict:
    """Drop empty cells and strip whitespace, so a row validates like a JSON object."""
    out = {}
    for k, v in row.items():
        if k is None:
            continue
        k = str(k).strip()
        if v is None:
            continue
        v = str(v).strip()
        if v != "":
            out[k] = v
    return out


def read_samples(path: str | Path) -> tuple[list[str], list[dict]]:
    """Return (header, rows). Rows have empty cells removed."""
    text = Path(path).read_text(encoding="utf-8-sig")
    reader = csv.DictReader(io.StringIO(text))
    header = [h.strip() for h in (reader.fieldnames or [])]
    rows = [clean_sample(r) for r in reader]
    rows = [r for r in rows if r]  # skip blank lines
    return header, rows


def samples_csv(samples: list[dict]) -> str:
    cols = SAMPLE_COLUMNS + [c for c in TRACK_COLUMNS if any(s.get(c) for s in samples)]
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(cols)
    for s in samples:
        w.writerow([s.get(c, "") or "" for c in cols])
    return buf.getvalue()


# ---------------------------------------------------------------- runs.yaml
class YamlError(ValueError):
    pass


def read_runs(path: str | Path) -> dict:
    try:
        data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    except yaml.YAMLError as e:
        raise YamlError(f"{path}: not valid YAML: {e}") from e
    return data if data is not None else {}


def _q(s) -> str:
    """A YAML double-quoted scalar (JSON string syntax is valid YAML)."""
    return json.dumps(str(s), ensure_ascii=False)


def runs_yaml(runs: list[dict]) -> str:
    lines = ["# flair-eval runs (schema flair-eval/runs/1). One entry per FLAIR configuration.",
             "runs:"]
    for r in runs:
        lines.append(f"  - id: {_q(r['id'])}")
        if r.get("label"):
            lines.append(f"    label: {_q(r['label'])}")
        lines.append("    flair:")
        for k in ("release", "git", "ref", "env", "container"):
            if k in r.get("flair", {}):
                lines.append(f"      {k}: {_q(r['flair'][k])}")
        lines.append(f"    transcriptome_args: {_q(r.get('transcriptome_args', ''))}")
        for k in ("quantify", "productivity", "own_align"):
            v = r.get(k, RUN_DEFAULTS[k])
            lines.append(f"    {k}: {'true' if v else 'false'}")
        tags = r.get("tags") or {}
        if tags:
            lines.append("    tags:")
            for k, v in tags.items():
                lines.append(f"      {_q(k)}: {_q(v)}")
    return "\n".join(lines) + "\n"


def run_is_default(run: dict, key: str) -> bool:
    return run.get(key, RUN_DEFAULTS[key]) == RUN_DEFAULTS[key]
