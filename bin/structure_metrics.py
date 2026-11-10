#!/usr/bin/env python3
"""Classify one FLAIR run's isoforms against a reference GTF and write its bundle.

Writes BUNDLE_DIR/bundle.json, BUNDLE_DIR/tables/isoforms.tsv and
BUNDLE_DIR/tables/categories.tsv (docs/CONTRACTS.md section 3).

Classification is SQANTI3-compatible (FSM/ISM/NIC/NNC/genic/genic_intron/antisense/
intergenic/fusion; FSM split into CIM/FSM5/FSM3/FSMTE with a 50-bp end window).
See bin/flair_eval_lib/classify.py for the exact rules.

Only main chromosomes and spike-in contigs are evaluated: chrM, unplaced/unlocalised
scaffolds, alts and decoys are dropped from both isoforms and reference.

With --regions, isoforms and reference transcripts that overlap a region are counted;
classification still sees the whole reference. With --genome, novel junctions are split
into canonical (GT-AG, GC-AG, AT-AC) and non-canonical; without it they are counted as
novel only (junctions.motif_checked = false).
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flair_eval_lib import BUNDLE_SCHEMA_ID, __version__  # noqa: E402
from flair_eval_lib.annotation import (  # noqa: E402
    InputError, Regions, assign_isoform_genes, keep_contig, read_bed12, read_gtf_gene_map,
    read_gtf_transcripts)
from flair_eval_lib.classify import FSM_END_WINDOW, Classifier, Reference  # noqa: E402
from flair_eval_lib.fasta import Fasta  # noqa: E402
from flair_eval_lib.metrics import (  # noqa: E402
    normalise_version, ref_in_scope, summarise, write_gene_names, write_tables)
from flair_eval_lib.schema import load_schema, validate  # noqa: E402

_ID_RE = r"^[A-Za-z0-9][A-Za-z0-9._-]*$"


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--isoforms-bed", required=True, help="FLAIR isoforms, BED12")
    p.add_argument("--isoforms-gtf", required=True,
                   help="FLAIR isoforms GTF (gene_id + transcript_id per isoform)")
    p.add_argument("--gtf", required=True, help="reference annotation GTF")
    p.add_argument("--sample", required=True, help="sample name")
    p.add_argument("--run-json", required=True, help="the run map from runs.yaml, as JSON")
    p.add_argument("--env-json", required=True, help="flair_env.json from FLAIR_ENV")
    p.add_argument("--out", required=True, help="bundle directory to write")
    p.add_argument("--regions", help="BED of regions to restrict the evaluation to")
    p.add_argument("--genome", help="genome FASTA (optional; enables junction motif check)")
    p.add_argument("--reads", help="reads path, recorded in bundle inputs only")
    p.add_argument("--bam", help="BAM path, recorded in bundle inputs only")
    p.add_argument("--end-window", type=int, default=FSM_END_WINDOW,
                   help="bp window for FSM end matching (default %(default)s)")
    return p.parse_args(argv)


def _load_json(path: str, what: str) -> dict:
    try:
        with open(path) as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError) as e:
        raise InputError(f"{what} {path}: {e}") from None
    if not isinstance(data, dict):
        raise InputError(f"{what} {path}: expected a JSON object")
    return data


def _pipeline_commit():
    env = os.environ.get("FLAIR_EVAL_COMMIT")
    if env:
        return env
    try:
        out = subprocess.run(["git", "-C", os.path.dirname(os.path.abspath(__file__)),
                              "rev-parse", "HEAD"], capture_output=True, text=True, timeout=10)
        return out.stdout.strip() or None if out.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        return None


def build_bundle(args):
    """Return (bundle dict, isoforms, calls, reference transcripts)."""
    if not re.match(_ID_RE, args.sample):
        raise InputError(f"--sample {args.sample!r} must match {_ID_RE}")
    run = _load_json(args.run_json, "--run-json")
    env = _load_json(args.env_json, "--env-json")
    run_id = run.get("id") or env.get("run_id")
    if not run_id or not re.match(_ID_RE, str(run_id)):
        raise InputError(f"run id {run_id!r} missing or does not match {_ID_RE}")
    if env.get("run_id") and run.get("id") and env["run_id"] != run["id"]:
        raise InputError(f"run id mismatch: run-json {run['id']!r} vs env-json {env['run_id']!r}")
    source = env.get("source")
    if source is None and isinstance(run.get("flair"), dict) and run["flair"]:
        source = next(iter(run["flair"]))
    if source not in ("release", "git", "env", "container"):
        raise InputError(f"flair source {source!r} must be release|git|env|container")

    regions = Regions.from_bed(args.regions) if args.regions else None
    ref_tx = [t for t in read_gtf_transcripts(args.gtf) if keep_contig(t.chrom)]
    if not ref_tx:
        raise InputError(f"{args.gtf}: no reference transcripts on evaluated contigs")
    isoforms = read_bed12(args.isoforms_bed)
    n_read = len(isoforms)
    isoforms = [t for t in isoforms if keep_contig(t.chrom)]
    n_contig = n_read - len(isoforms)
    if regions is not None:
        isoforms = [t for t in isoforms if regions.overlaps(t.chrom, t.start, t.end)]
    assign_isoform_genes(isoforms, read_gtf_gene_map(args.isoforms_gtf))

    ref = Reference(ref_tx)
    clf = Classifier(ref, args.end_window)
    calls = [clf.classify(t) for t in isoforms]
    fasta = Fasta(args.genome) if args.genome else None
    try:
        metrics = summarise(isoforms, calls, ref_in_scope(ref_tx, regions), ref, fasta)
    finally:
        if fasta is not None:
            fasta.close()

    tags = run.get("tags") or {}
    bundle = {
        "schema": BUNDLE_SCHEMA_ID,
        "sample": args.sample,
        "run_id": run_id,
        "label": run.get("label") or run_id,
        "tags": {str(k): str(v) for k, v in tags.items()},
        "flair": {
            "source": source,
            "requested": str(env.get("requested") or ""),
            "version": normalise_version(env.get("flair_version")),
            "commit": env.get("resolved_commit"),
        },
        "pipeline": {
            "version": os.environ.get("FLAIR_EVAL_VERSION", __version__),
            "commit": _pipeline_commit(),
            "created": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        },
        "inputs": {
            "genome": args.genome, "gtf": args.gtf, "reads": args.reads, "bam": args.bam,
            "regions": args.regions,
            "regions_intervals": regions.as_list() if regions else None,
            "isoforms_in_file": n_read,
            "isoforms_dropped_contigs": n_contig,
            "isoforms_outside_regions": n_read - n_contig - len(isoforms),
            "fsm_end_window": args.end_window,
            "transcriptome_args": run.get("transcriptome_args", ""),
        },
        **metrics,
        "compute": {"cpu_hours": None, "peak_rss_gb": None, "wall_minutes": None},
        "sections": {},
    }
    return bundle, isoforms, calls, ref_tx


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        bundle, isoforms, calls, ref_tx = build_bundle(args)
    except InputError as e:
        print(f"structure_metrics.py: error: {e}", file=sys.stderr)
        return 2
    except FileNotFoundError as e:
        print(f"structure_metrics.py: error: {e}", file=sys.stderr)
        return 2
    try:
        errs = validate(bundle, load_schema())
    except OSError as e:
        print(f"structure_metrics.py: warning: schema not found, not validated ({e})",
              file=sys.stderr)
        errs = []
    if errs:
        print("structure_metrics.py: error: bundle fails schema:\n  " + "\n  ".join(errs),
              file=sys.stderr)
        return 3
    os.makedirs(args.out, exist_ok=True)
    write_tables(args.out, isoforms, calls)
    write_gene_names(args.out, calls, ref_tx)
    with open(os.path.join(args.out, "bundle.json"), "w") as fh:
        json.dump(bundle, fh, indent=1)
        fh.write("\n")
    s = bundle["summary"]
    print(f"{bundle['sample']}/{bundle['run_id']}: {s['isoforms']} isoforms, "
          f"FSM {bundle['categories']['FSM']}, refs FSM {s['ref_transcripts_fsm']}/"
          f"{s['ref_transcripts_total']}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
