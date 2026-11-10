#!/usr/bin/env python3
"""Build one bundle folder (docs/CONTRACTS.md section 3) from a job's step outputs.

Usage:
  flair_eval_bundle.py --info info.json --env-json flair_env.json --steps DIR --out BUNDLE_DIR

DIR holds the step folders written by the pipeline: assemble/ and eval/ (required), and
quantify/, productivity/, align/ (own_align runs) when they ran. info.json is written by the
BUNDLE process: sample, run_id, label, tags, inputs, pipeline, compute, trace rows.
Stdlib only.
"""
import argparse
import datetime as dt
import json
import shutil
import sys
from pathlib import Path

REQUIRED_KEYS = ["schema", "sample", "run_id", "flair", "pipeline", "summary", "categories",
                 "fsm_subcategories", "junctions", "histograms", "sections"]
# CONTRACTS section 7, A15. Raw units: realtime/duration ms, %cpu percent, peak_rss/rss bytes.
TRACE_COLUMNS = ["task_id", "name", "status", "exit", "realtime", "%cpu", "peak_rss", "rss", "duration"]
# sections.quant cpm_hist: fixed edges so runs compare bin by bin; last bin is open-ended.
CPM_EDGES = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
STEP_ORDER = ["align", "assemble", "quantify", "productivity", "eval"]


def die(msg):
    sys.exit(f"flair_eval_bundle.py: {msg}")


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    return str(v)


def read_tsv(path: Path):
    """Rows of a one-header TSV as dicts (CONTRACTS section 7, A10)."""
    with open(path) as fh:
        lines = [l.rstrip("\n").rstrip("\r") for l in fh if l.strip()]
    if not lines:
        return []
    head = lines[0].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[1:]]


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def histogram(values, edges):
    counts = [0] * len(edges)
    for v in values:
        i = 0
        while i + 1 < len(edges) and v >= edges[i + 1]:
            i += 1
        counts[i] += 1
    return counts


def quant_section(counts_tsv: Path):
    """sections.quant (CONTRACTS section 8) from counts.tsv: isoform_id, count, cpm."""
    rows = read_tsv(counts_tsv)
    counts = [num(r.get("count")) for r in rows]
    cpms = [num(r.get("cpm")) for r in rows]
    assigned = sum(c for c in counts if c is not None)
    return {
        "isoforms_quantified": sum(1 for c in counts if c is not None and c > 0),
        "reads_assigned": int(assigned) if float(assigned).is_integer() else round(assigned, 3),
        "cpm_hist": {"edges": CPM_EDGES, "counts": histogram([c for c in cpms if c is not None], CPM_EDGES)},
    }


def productivity_section(prod_tsv: Path):
    """sections.productivity (CONTRACTS section 8) from productivity.tsv: ..., class."""
    classes = {}
    for r in read_tsv(prod_tsv):
        c = (r.get("class") or "").strip() or "unknown"
        classes[c] = classes.get(c, 0) + 1
    return {"class_counts": dict(sorted(classes.items()))}


def copy_if(src: Path, dst: Path) -> bool:
    if src.is_file():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        return True
    return False


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--info", required=True, type=Path)
    ap.add_argument("--env-json", required=True, type=Path)
    ap.add_argument("--steps", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)

    info = json.loads(a.info.read_text())
    steps = a.steps
    asm, ev = steps / "assemble", steps / "eval"
    if not (ev / "bundle.json").is_file():
        die(f"{ev}/bundle.json is missing: structure_metrics.py did not write it")
    for name in ("isoforms.bed", "isoforms.gtf"):
        if not (asm / name).is_file():
            die(f"{asm}/{name} is missing: the adapter did not write it")

    out = a.out
    if out.exists():
        shutil.rmtree(out)
    (out / "tables").mkdir(parents=True)
    (out / "provenance").mkdir()

    # Structure tables and isoforms.
    for t in sorted((ev / "tables").glob("*")):
        if t.is_file():
            shutil.copyfile(t, out / "tables" / t.name)
    for name in ("isoforms.bed", "isoforms.gtf"):
        shutil.copyfile(asm / name, out / name)
    copy_if(steps / "quantify" / "counts.tsv", out / "counts.tsv")
    copy_if(steps / "productivity" / "productivity.tsv", out / "productivity.tsv")

    # Provenance.
    shutil.copyfile(a.env_json, out / "provenance" / "flair_env.json")
    for step in STEP_ORDER:
        d = steps / step
        copy_if(d / "adapter_log.json", out / "provenance" / f"adapter_{step}.json")
        copy_if(d / "skipped.json", out / "provenance" / f"adapter_{step}_skipped.json")
    with open(out / "provenance" / "trace.tsv", "w") as fh:
        fh.write("\t".join(TRACE_COLUMNS) + "\n")
        for row in info.get("trace", []):
            fh.write("\t".join(fmt(row.get(c)) for c in TRACE_COLUMNS) + "\n")

    # bundle.json: keep what structure_metrics wrote, fill what only the pipeline knows.
    b = json.loads((ev / "bundle.json").read_text())
    missing = [k for k in REQUIRED_KEYS if k not in b and k != "pipeline"]
    if missing:
        die(f"eval/bundle.json lacks required key(s) {', '.join(missing)} (assets/schema/bundle.schema.json)")
    for k in ("sample", "run_id"):
        if b.get(k) != info[k]:
            die(f"eval/bundle.json has {k}={b.get(k)!r}, expected {info[k]!r}")
    for k in ("label", "tags"):
        if b.get(k) in (None, "", {}) and info.get(k) not in (None, ""):
            b[k] = info[k]
    # CONTRACTS section 7, A12: core owns pipeline, inputs and compute. Values core knows win;
    # anything else structure_metrics recorded (e.g. its own git commit) is kept.
    pipe = b.get("pipeline") or {}
    for k, v in info["pipeline"].items():
        if v not in (None, ""):
            pipe[k] = v
    if not pipe.get("created"):
        pipe["created"] = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    b["pipeline"] = pipe
    inputs = b.get("inputs") or {}
    for k, v in info["inputs"].items():
        inputs[k] = v           # original absolute paths, not the staged task-local names
    b["inputs"] = inputs
    b["compute"] = info["compute"]
    sections = b.get("sections") or {}
    if (out / "counts.tsv").is_file():
        sections["quant"] = quant_section(out / "counts.tsv")
    if (out / "productivity.tsv").is_file():
        sections["productivity"] = productivity_section(out / "productivity.tsv")
    b["sections"] = sections
    (out / "bundle.json").write_text(json.dumps(b, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
