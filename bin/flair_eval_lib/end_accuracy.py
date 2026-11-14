"""Isoform end accuracy against end sets, scored by ted_end_precision.py (CONTRACTS A36-A40).

An *end set* is a BED6 of stranded sites for one end of the transcript:
  kind "truth"      a sample's own orthogonal truth calls from the --truth sheet (CAGE, PRO-cap,
                    RAMPAGE for the 5' end; dRNA for the 3' end). Precision, recall and F1.
  kind "reference"  the end databases (--ref_tss refTSS, --ref_pas PolyA_DB). Precision only:
                    the databases pool many tissues, so recall against them is not reported.

This module does not score anything itself. Every number comes from flair4_ted's end judge,
ted_end_precision.py, copied unchanged into modules/end_accuracy/resources/usr/bin/ (with its
helpers peak_matching.py and signal_utils.py). Here we only plan the judge calls, run them,
read the one-row precision_recall_summary.tsv each call writes, and shape the bundle section.

The judge's rule (documented, not changed): an isoform end is a hit when it lies on the same
strand within `window` bp of any part of a site interval, on chr1-22/X/Y. Single-end precision
and recall are many-to-one (several isoform ends may share one site). Precision = isoform ends
near a site / isoform ends; recall = sites with an isoform end near them / sites in scope. The
paired number counts isoforms whose 5' and 3' ends both hit, deduplicated within a junction chain
by the (5' site, 3' site) pair. A truth BED's interval is the call's 1-bp mode point, so on these
files the rule is "within 50 bp of the mode" (ortho truth project ADR 0001).

Judge calls per job, for each kind separately: every 5' set paired with every 3' set when the
kind has both, otherwise each set alone. A set's single-end numbers are identical in every call
it appears in; the first call's are kept.
"""

from __future__ import annotations

import csv
import json
import os
import re
import subprocess
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from .annotation import InputError
from .ends import md5sum

JUDGE = "ted_end_precision.py"
JUDGE_SOURCE = "flair4_ted bin/evaluation (commit 40fdd21), copied unchanged"
WINDOW = 50
KINDS = ("reference", "truth")
ENDS = ("tss", "tes")
SET_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
REFERENCE_IDS = {"tss": "refTSS", "tes": "PolyA_DB"}

# The judge's column prefix for each end.
_PREFIX = {"tss": "5prime", "tes": "3prime"}

TABLE_COLUMNS = ["kind", "end", "set", "assay", "partner", "isoforms", "isoform_ends",
                 "sites", "sites_hit", "precision", "naive_precision", "dedup_tp",
                 "redundant_calls", "recall", "f1", "unique_pairs", "both_hit",
                 "both_hit_raw", "pair_precision", "pair_naive_precision"]


@dataclass
class EndSet:
    id: str
    kind: str                  # truth | reference
    end: str                   # tss | tes
    assay: str
    path: str                  # staged path the judge reads
    source: Optional[str] = None   # original path (for the record)
    provisional: bool = False

    def meta(self) -> Dict:
        return {"id": self.id, "kind": self.kind, "end": self.end, "assay": self.assay,
                "file": os.path.basename(self.source or self.path), "path": self.source,
                "md5": md5sum(self.path), "provisional": self.provisional}


@dataclass
class Call:
    kind: str
    tss: Optional[EndSet] = None
    tes: Optional[EndSet] = None
    summary: Dict = field(default_factory=dict)

    @property
    def name(self) -> str:
        return "__".join(s.id for s in (self.tss, self.tes) if s is not None)


def read_sets(path: str) -> List[EndSet]:
    """sets.json written by the END_ACCURACY process: a list of EndSet fields."""
    try:
        with open(path) as fh:
            raw = json.load(fh)
    except (OSError, json.JSONDecodeError) as e:
        raise InputError(f"--sets {path}: {e}") from None
    if not isinstance(raw, list):
        raise InputError(f"--sets {path}: expected a JSON list")
    sets, seen = [], set()
    for i, r in enumerate(raw):
        try:
            s = EndSet(id=str(r["id"]), kind=r["kind"], end=r["end"], assay=str(r.get("assay") or r["id"]),
                       path=r["path"], source=r.get("source"), provisional=bool(r.get("provisional", False)))
        except (KeyError, TypeError) as e:
            raise InputError(f"--sets {path}: entry {i}: missing {e}") from None
        if s.kind not in KINDS:
            raise InputError(f"--sets {path}: entry {i}: kind {s.kind!r} must be one of {KINDS}")
        if s.end not in ENDS:
            raise InputError(f"--sets {path}: entry {i}: end {s.end!r} must be one of {ENDS}")
        if not SET_ID.match(s.id):
            raise InputError(f"--sets {path}: entry {i}: id {s.id!r} must match {SET_ID.pattern}")
        if s.id in seen:
            raise InputError(f"--sets {path}: id {s.id!r} is used twice")
        seen.add(s.id)
        if not os.path.isfile(s.path):
            raise InputError(f"--sets {path}: entry {i}: {s.path} does not exist")
        sets.append(s)
    return sets


def plan_calls(sets: Sequence[EndSet]) -> List[Call]:
    """Judge calls, per kind: 5' x 3' when the kind has both ends, otherwise each set alone.
    Order: kinds in KINDS order, then sets in the order given."""
    calls: List[Call] = []
    for kind in KINDS:
        t5 = [s for s in sets if s.kind == kind and s.end == "tss"]
        t3 = [s for s in sets if s.kind == kind and s.end == "tes"]
        if t5 and t3:
            calls += [Call(kind, a, b) for a in t5 for b in t3]
        else:
            calls += [Call(kind, a, None) for a in t5] + [Call(kind, None, b) for b in t3]
    return calls


def regions_args(bed: Optional[str]) -> List[str]:
    """A regions BED as the judge's --region values (chr:start-end, BED coordinates)."""
    if not bed:
        return []
    out = []
    with open(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            p = line.split("\t")
            if len(p) < 3:
                raise InputError(f"{bed}: regions BED needs 3 columns")
            out.append(f"{p[0]}:{int(p[1])}-{int(p[2].strip())}")
    return out


def judge_command(judge: str, isoforms_bed: str, gtf: str, call: Call, outdir: str,
                  regions: Sequence[str] = (), window: int = WINDOW) -> List[str]:
    cmd = [judge, "--isoforms-bed", isoforms_bed, "--gtf", gtf, "--window", str(window),
           "--mode", call.name, "--outdir", outdir]
    if call.tss is not None:
        cmd += ["--peaks-5prime", call.tss.path]
    if call.tes is not None:
        cmd += ["--peaks-3prime", call.tes.path]
    if regions:
        cmd += ["--region", *regions]
    return cmd


def _num(v: str):
    if v is None or v == "":
        return None
    try:
        return int(v)
    except ValueError:
        return float(v)


def read_summary(path: str) -> Dict:
    """The judge's one-row precision_recall_summary.tsv as {column: int | float | None}."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if len(rows) != 1:
        raise InputError(f"{path}: expected one data row, got {len(rows)}")
    return {k: (v if k == "transcriptome_mode" else _num(v)) for k, v in rows[0].items()}


def run_calls(judge: str, isoforms_bed: str, gtf: str, calls: Sequence[Call], workdir: str,
              regions: Sequence[str] = (), window: int = WINDOW) -> None:
    """Run the judge once per call (outputs under workdir/<call name>/) and attach summaries."""
    for c in calls:
        out = os.path.join(workdir, c.name)
        os.makedirs(out, exist_ok=True)
        cmd = judge_command(judge, isoforms_bed, gtf, c, out, regions, window)
        with open(os.path.join(out, "judge.log"), "w") as log:
            res = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
        if res.returncode != 0:
            raise InputError(f"{JUDGE} failed (exit {res.returncode}) for {c.name}; "
                             f"see {os.path.join(out, 'judge.log')}")
        c.summary = read_summary(os.path.join(out, "precision_recall_summary.tsv"))


def _end_row(c: Call, end: str) -> Dict:
    s = c.tss if end == "tss" else c.tes
    p, x = _PREFIX[end], c.summary
    truth = s.kind == "truth"
    return {
        "set": s.id, "kind": s.kind, "end": end, "assay": s.assay,
        "isoforms": x.get("n_isoforms_total"),
        "isoform_ends": x.get(f"{p}_n_isoform_ends"),
        "sites": x.get(f"{p}_n_annot_total"),
        "sites_hit": x.get(f"{p}_n_annot_matched") if truth else None,
        "precision": x.get(f"{p}_precision"),
        "naive_precision": x.get(f"{p}_naive_precision"),
        "dedup_tp": x.get(f"{p}_dedup_tp"),
        "redundant_calls": x.get(f"{p}_redundant_calls"),
        "recall": x.get(f"{p}_recall") if truth else None,
        "f1": x.get(f"{p}_f1") if truth else None,
    }


def _pair_row(c: Call) -> Dict:
    x = c.summary
    return {
        "kind": c.kind, "tss_set": c.tss.id, "tes_set": c.tes.id,
        "isoforms": x.get("paired_n_isoforms"),
        "unique_pairs": x.get("paired_unique_pairs"),
        "both_hit": x.get("paired_both_hit"),
        "both_hit_raw": x.get("paired_both_hit_raw"),
        "precision": x.get("paired_dedup_precision"),
        "naive_precision": x.get("paired_naive_precision"),
    }


def build_section(sets: Sequence[EndSet], calls: Sequence[Call], judge: str,
                  window: int = WINDOW, regions: Sequence[str] = ()) -> Dict:
    """sections.end_accuracy (A38). `ends` has one row per set, `pairs` one per 5' x 3' call."""
    ends: List[Dict] = []
    done = set()
    for c in calls:
        for end in ENDS:
            s = c.tss if end == "tss" else c.tes
            if s is not None and s.id not in done:
                done.add(s.id)
                ends.append(_end_row(c, end))
    order = {s.id: i for i, s in enumerate(sets)}
    ends.sort(key=lambda r: (KINDS.index(r["kind"]), ENDS.index(r["end"]), order[r["set"]]))
    pairs = [_pair_row(c) for c in calls if c.tss is not None and c.tes is not None]
    return {
        "scorer": {"script": JUDGE, "source": JUDGE_SOURCE, "md5": md5sum(judge),
                   "window": window, "regions": list(regions) or None},
        "provisional": any(s.provisional for s in sets if s.kind == "truth"),
        "sets": [s.meta() for s in sets],
        "ends": ends,
        "pairs": pairs,
    }


def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def write_table(path: str, section: Dict) -> None:
    """tables/end_accuracy.tsv (A39): the section's `ends` rows, then its `pairs` rows."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    assay = {s["id"]: s["assay"] for s in section["sets"]}
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(TABLE_COLUMNS)
        for r in section["ends"]:
            row = dict(r, partner=None)
            w.writerow([_fmt(row.get(c)) for c in TABLE_COLUMNS])
        for p in section["pairs"]:
            row = {"kind": p["kind"], "end": "pair", "set": p["tss_set"],
                   "assay": assay.get(p["tss_set"]), "partner": p["tes_set"],
                   "isoforms": p["isoforms"], "unique_pairs": p["unique_pairs"],
                   "both_hit": p["both_hit"], "both_hit_raw": p["both_hit_raw"],
                   "pair_precision": p["precision"], "pair_naive_precision": p["naive_precision"]}
            w.writerow([_fmt(row.get(c)) for c in TABLE_COLUMNS])
