#!/usr/bin/env python3
"""Build one bundle folder (docs/CONTRACTS.md section 3) from a job's step outputs.

Usage:
  flair_eval_bundle.py --info info.json --env-json flair_env.json --steps DIR --out BUNDLE_DIR \
                       --lib-dir PIPELINE_BIN

sections.quant / sections.productivity are built by `python3 -m flair_eval_lib.sections`
(CONTRACTS section 9, A19), found in --lib-dir.

DIR holds the step folders written by the pipeline: assemble/ and eval/ (required), and
quantify/, productivity/, align/ (own_align runs) when they ran. info.json is written by the
BUNDLE process: sample, run_id, label, tags, inputs, pipeline, compute, trace rows.
Stdlib only.
"""
import argparse
import datetime as dt
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

REQUIRED_KEYS = ["schema", "sample", "run_id", "flair", "pipeline", "summary", "categories",
                 "fsm_subcategories", "junctions", "histograms", "sections"]
# CONTRACTS section 7, A15. Raw units: realtime/duration ms, %cpu percent, peak_rss/rss bytes.
TRACE_COLUMNS = ["task_id", "name", "status", "exit", "realtime", "%cpu", "peak_rss", "rss", "duration"]
STEP_ORDER = ["align", "assemble", "quantify", "productivity", "eval"]


def die(msg):
    sys.exit(f"flair_eval_bundle.py: {msg}")


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    return str(v)


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
    ap.add_argument("--lib-dir", required=True, type=Path,
                    help="directory holding flair_eval_lib (the pipeline's bin/)")
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
    if info.get("run_index") is not None:
        b["run_index"] = int(info["run_index"])          # A17: position in runs.yaml / --flair
    (out / "bundle.json").write_text(json.dumps(b, indent=2) + "\n")

    # A19: sections.quant / sections.productivity come from the one shared implementation.
    env = dict(os.environ)
    env["PYTHONPATH"] = str(a.lib_dir) + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    subprocess.run([sys.executable, "-m", "flair_eval_lib.sections", str(out)], env=env, check=True)

    # A21: how this version quantified (reads | bam), as the adapter recorded it.
    b = json.loads((out / "bundle.json").read_text())
    qlog = out / "provenance" / "adapter_quantify.json"
    quant = (b.get("sections") or {}).get("quant")
    if quant is not None and qlog.is_file():
        qin = json.loads(qlog.read_text()).get("quantify_input")
        if qin:
            quant["input"] = qin
    (out / "bundle.json").write_text(json.dumps(b, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
