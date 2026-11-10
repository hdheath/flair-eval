"""Build the per-bundle `sections.quant` and `sections.productivity` blocks (CONTRACTS section 8)
from the adapter's counts.tsv / productivity.tsv. Core's BUNDLE step can call these; the report
fixtures are made with them too.

    python3 -m flair_eval_lib.sections BUNDLE_DIR   # adds both sections when the TSVs exist
"""

from __future__ import annotations

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


def add_sections(bundle_dir: str) -> dict:
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
    if os.path.exists(os.path.join(bundle_dir, "productivity.tsv")):
        sec["productivity"] = productivity_section(os.path.join(bundle_dir, "productivity.tsv"), ids)
    with open(path, "w") as fh:
        json.dump(bundle, fh, indent=1)
        fh.write("\n")
    return bundle


if __name__ == "__main__":
    for d in sys.argv[1:]:
        b = add_sections(d)
        print(d, sorted(b.get("sections", {})), file=sys.stderr)
