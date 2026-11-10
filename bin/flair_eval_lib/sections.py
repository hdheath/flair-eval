"""Build the per-bundle `sections.quant` and `sections.productivity` blocks (CONTRACTS section 8)
from the adapter's counts.tsv / productivity.tsv. Core's BUNDLE step can call these; the report
fixtures are made with them too.

    python3 -m flair_eval_lib.sections BUNDLE_DIR [...] [--quant-input reads|bam]

adds both sections when the TSVs exist. `sections.quant.input` (CONTRACTS A21) says whether
FLAIR quantify counted from reads or from the genome BAM: --quant-input if given, else
`quantify_input` from provenance/adapter_quantify.json, else inferred from the quantify command
that log records (`--manifest` = BAM manifest, `-r` = reads manifest), else left out.
"""

from __future__ import annotations

import argparse
import bisect
import json
import os
import sys
from typing import Dict, Optional

from .concordance import read_header_tsv

# Fixed CPM bins so runs compare bin by bin; the last bin is open-ended.
CPM_EDGES = [0, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000]


def quant_section(counts_tsv: str, isoform_ids: Optional[set] = None) -> dict:
    rows = read_header_tsv(counts_tsv, ["isoform_id", "count", "cpm"])
    if isoform_ids is not None:
        rows = [r for r in rows if r["isoform_id"] in isoform_ids]
    counts = [0] * len(CPM_EDGES)
    quantified = 0
    reads = 0.0
    for r in rows:
        c, cpm = float(r["count"] or 0), float(r["cpm"] or 0)
        reads += c
        if c > 0:
            quantified += 1
            counts[max(0, bisect.bisect_right(CPM_EDGES, cpm) - 1)] += 1
    return {"isoforms_quantified": quantified, "reads_assigned": int(round(reads)),
            "cpm_hist": {"edges": list(CPM_EDGES), "counts": counts}}


def productivity_section(productivity_tsv: str, isoform_ids: Optional[set] = None) -> dict:
    rows = read_header_tsv(productivity_tsv, ["isoform_id", "class"])
    classes: Dict[str, int] = {}
    for r in rows:
        if isoform_ids is not None and r["isoform_id"] not in isoform_ids:
            continue
        k = r["class"] or "unknown"
        classes[k] = classes.get(k, 0) + 1
    return {"class_counts": dict(sorted(classes.items()))}


QUANT_INPUTS = ("reads", "bam")


def quant_input_from_log(log: dict) -> Optional[str]:
    """'reads' or 'bam' from an adapter_quantify.json, or None if it cannot tell."""
    v = log.get("quantify_input")
    if v in QUANT_INPUTS:
        return v
    for c in log.get("commands") or []:
        argv = [str(a) for a in c.get("argv") or []]
        if "quantify" not in argv:
            continue
        if "--manifest" in argv:
            return "bam"
        if "-r" in argv or "--reads" in argv:
            return "reads"
    return None


def quant_input_from_dir(bundle_dir: str) -> Optional[str]:
    path = os.path.join(bundle_dir, "provenance", "adapter_quantify.json")
    try:
        with open(path) as fh:
            return quant_input_from_log(json.load(fh))
    except (OSError, ValueError):
        return None


def add_sections(bundle_dir: str, quant_input: Optional[str] = None) -> dict:
    path = os.path.join(bundle_dir, "bundle.json")
    with open(path) as fh:
        bundle = json.load(fh)
    ids = None
    table = os.path.join(bundle_dir, "tables", "isoforms.tsv")
    if os.path.exists(table):
        ids = {r["isoform_id"] for r in read_header_tsv(table, ["isoform_id"])}
    sec = bundle.setdefault("sections", {})
    if os.path.exists(os.path.join(bundle_dir, "counts.tsv")):
        sec["quant"] = quant_section(os.path.join(bundle_dir, "counts.tsv"), ids)
        qi = quant_input or quant_input_from_dir(bundle_dir)
        if qi:
            sec["quant"]["input"] = qi
    if os.path.exists(os.path.join(bundle_dir, "productivity.tsv")):
        sec["productivity"] = productivity_section(os.path.join(bundle_dir, "productivity.tsv"), ids)
    with open(path, "w") as fh:
        json.dump(bundle, fh, indent=1)
        fh.write("\n")
    return bundle


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bundle_dirs", nargs="+")
    ap.add_argument("--quant-input", choices=QUANT_INPUTS)
    a = ap.parse_args()
    for d in a.bundle_dirs:
        b = add_sections(d, a.quant_input)
        q = b.get("sections", {}).get("quant", {})
        print(d, sorted(b.get("sections", {})), "quant input:", q.get("input"), file=sys.stderr)
