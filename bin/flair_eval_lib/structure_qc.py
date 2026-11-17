"""sections.structure_qc (CONTRACTS section 13): SQANTI3-style quality checks for one run.

Each key follows a figure or column of SQANTI3's QC report, written anew here (no SQANTI3 code
is copied; SQANTI3 is GPL-3):

  subcategories    ISM / NIC / NNC subcategories (classify.Call.detail)
  redundancy       FSM and FSM+ISM isoforms per matched reference transcript
  canonical        multi-exon isoforms whose junctions are all canonical (needs the genome)
  rts              junctions with an RT-switching repeat (needs the genome)
  by_category      transcript length and exon count histograms per structural category
  ref_recovery     reference transcripts in scope, and those with an FSM, by length and exons
  priming          internal-priming risk by category and mono/multi-exon (needs the genome);
                   coding/non-coding is added by sections.add_sections when productivity ran
  junction_support support of junctions in the sample's `junctions` file (short reads when
                   that file comes from short reads; see CONTRACTS section 13)

RT switching (SQANTI3 utilities/rt_switching.py, defaults as sqanti3_qc runs it): for a
junction, take the 12 nt ending one base into the intron at the donor and the 12 nt ending one
base into the downstream exon at the acceptor (10 nt plus 1 nt of wiggle each side), read on
the transcript strand. The junction is RT-switching when the two share an 8-nt stretch, at the
same or different offsets, with at most one mismatch and no indels.
"""

from __future__ import annotations

import os
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from . import CATEGORIES, HISTOGRAM_EDGES
from .annotation import InputError, Transcript, _open
from .classify import Call, Reference
from .fasta import Fasta

SUBCATEGORIES = {
    "ISM": ("3prime_fragment", "5prime_fragment", "internal_fragment", "intron_retention",
            "mono-exon"),
    "NIC": ("combination_of_known_junctions", "combination_of_known_splicesites",
            "intron_retention", "mono-exon_by_intron_retention"),
    "NNC": ("at_least_one_novel_splicesite", "intron_retention"),
}
REDUNDANCY_EDGES = [1, 2, 3, 4, 5, 6, 8, 11]
JUNCTION_TYPES = ("known_canonical", "known_noncanonical", "novel_canonical",
                  "novel_noncanonical")
RTS_PATTERN = 10
RTS_WIGGLE = 1
RTS_MIN_MATCH = 8
RTS_MISMATCHES = 1
RISK_TENTHS = 6          # ends.RISK_A_FRAC in tenths
SPLICED = ("FSM", "ISM", "NIC", "NNC", "fusion")
_COMP = str.maketrans("ACGTN", "TGCAN")


def _hist(values: Iterable[float], edges: Sequence[float]) -> List[int]:
    from .metrics import histogram
    return histogram(list(values), edges)


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


# -- RT switching --------------------------------------------------------------------------
def rts_windows(fasta: Fasta, chrom: str, strand: str, a: int, b: int) -> Tuple[str, str]:
    """(donor-side, acceptor-side) windows of intron [a, b), on the transcript strand."""
    n, w = RTS_PATTERN, RTS_WIGGLE
    if strand == "+":
        return (fasta.fetch(chrom, a - n - w, a + w), fasta.fetch(chrom, b - n - w, b + w))
    return (_revcomp(fasta.fetch(chrom, b - w, b + n + w)),
            _revcomp(fasta.fetch(chrom, a - w, a + n + w)))


def has_repeat(x: str, y: str, k: int = RTS_MIN_MATCH, mismatches: int = RTS_MISMATCHES) -> bool:
    """Do x and y share a k-nt stretch with at most `mismatches` differences?"""
    if len(x) < k or len(y) < k:
        return False
    for i in range(len(x) - k + 1):
        for j in range(len(y) - k + 1):
            diff = 0
            for p in range(k):
                if x[i + p] != y[j + p]:
                    diff += 1
                    if diff > mismatches:
                        break
            if diff <= mismatches:
                return True
    return False


# -- the sample's junctions file ------------------------------------------------------------
def read_junction_file(path: str) -> Dict[Tuple[str, int, int], int]:
    """(chrom, intron start, intron end) -> reads, 0-based half-open, strand ignored.

    STAR SJ.out.tab (9 columns, 1-based first and last intron base, unique reads in column 7)
    or BED (intron [start, end), reads in the score column). Told apart by the column count."""
    out: Dict[Tuple[str, int, int], int] = {}
    n = 0
    with _open(path) as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            try:
                if len(f) == 9 and f[3] in ("0", "1", "2"):
                    key, reads = (f[0], int(f[1]) - 1, int(f[2])), int(f[6])
                elif len(f) >= 5:
                    key, reads = (f[0], int(f[1]), int(f[2])), int(float(f[4]))
                else:
                    raise ValueError("expected STAR SJ.out.tab (9 columns) or BED (5+ columns)")
            except ValueError as e:
                raise InputError(f"{path} line {ln}: {e}") from None
            out[key] = out.get(key, 0) + reads
            n += 1
    if not n:
        raise InputError(f"{path}: no junctions")
    return out


# -- the section ---------------------------------------------------------------------------
def structure_qc_section(isoforms: Sequence[Transcript], calls: Sequence[Call],
                         ref: Reference, ref_scope: Sequence[Transcript],
                         fasta: Optional[Fasta], rows=None,
                         jf: Optional[Dict[Tuple[str, int, int], int]] = None,
                         jf_meta: Optional[Dict] = None) -> Dict:
    sec: Dict = {
        "subcategories": _subcategories(calls),
        "redundancy": _redundancy(calls),
        "by_category": _by_category(isoforms, calls),
        "ref_recovery": _ref_recovery(calls, ref_scope),
    }
    if fasta is not None:
        canon, types = _junction_types(isoforms, ref, fasta)
        sec["canonical"] = _canonical(isoforms, calls, canon)
        sec["rts"] = _rts(isoforms, fasta, types)
        if rows is not None:
            sec["priming"] = _priming(isoforms, rows)
    if jf is not None:
        sec["junction_support"] = _junction_support(isoforms, calls, ref, jf, jf_meta or {})
    return sec


def _subcategories(calls: Sequence[Call]) -> Dict:
    out = {cat: {s: 0 for s in subs} for cat, subs in SUBCATEGORIES.items()}
    for c in calls:
        if c.category in out and c.detail in out[c.category]:
            out[c.category][c.detail] += 1
    return out


def _redundancy(calls: Sequence[Call]) -> Dict:
    fsm: Dict[str, int] = {}
    both: Dict[str, int] = {}
    for c in calls:
        if c.category in ("FSM", "ISM") and c.ref_transcript:
            both[c.ref_transcript] = both.get(c.ref_transcript, 0) + 1
            if c.category == "FSM":
                fsm[c.ref_transcript] = fsm.get(c.ref_transcript, 0) + 1
    return {
        "fsm_isoforms": sum(fsm.values()), "fsm_refs": len(fsm),
        "fsm_refs_multi": sum(1 for v in fsm.values() if v > 1),
        "fsm_ism_isoforms": sum(both.values()), "fsm_ism_refs": len(both),
        "fsm_ism_refs_multi": sum(1 for v in both.values() if v > 1),
        "fsm_hist": {"edges": list(REDUNDANCY_EDGES), "counts": _hist(fsm.values(), REDUNDANCY_EDGES)},
        "fsm_ism_hist": {"edges": list(REDUNDANCY_EDGES),
                         "counts": _hist(both.values(), REDUNDANCY_EDGES)},
    }


def _by_category(isoforms: Sequence[Transcript], calls: Sequence[Call]) -> Dict:
    le, ee = HISTOGRAM_EDGES["transcript_length"], HISTOGRAM_EDGES["exons_per_isoform"]
    lens: Dict[str, List[int]] = {c: [] for c in CATEGORIES}
    exons: Dict[str, List[int]] = {c: [] for c in CATEGORIES}
    for t, c in zip(isoforms, calls):
        lens[c.category].append(t.length)
        exons[c.category].append(t.n_exons)
    return {
        "length": {"edges": list(le), "counts": {c: _hist(v, le) for c, v in lens.items()}},
        "exons": {"edges": list(ee), "counts": {c: _hist(v, ee) for c, v in exons.items()}},
    }


def _ref_recovery(calls: Sequence[Call], ref_scope: Sequence[Transcript]) -> Dict:
    le, ee = HISTOGRAM_EDGES["transcript_length"], HISTOGRAM_EDGES["exons_per_isoform"]
    hit = {c.ref_transcript for c in calls if c.category == "FSM"}
    found = [t for t in ref_scope if t.tid in hit]
    return {
        "length": {"edges": list(le), "all": _hist((t.length for t in ref_scope), le),
                   "fsm": _hist((t.length for t in found), le)},
        "exons": {"edges": list(ee), "all": _hist((t.n_exons for t in ref_scope), ee),
                  "fsm": _hist((t.n_exons for t in found), ee)},
    }


def _junction_types(isoforms: Sequence[Transcript], ref: Reference,
                    fasta: Fasta) -> Tuple[Dict, Dict]:
    """Per distinct stranded junction: canonical? and its known/novel x canonical type."""
    canon: Dict[Tuple[str, str, int, int], bool] = {}
    types: Dict[Tuple[str, str, int, int], str] = {}
    for t in isoforms:
        for a, b in t.introns:
            key = (t.chrom, t.strand, a, b)
            if key in canon:
                continue
            canon[key] = fasta.is_canonical(t.chrom, t.strand, a, b)
            known = ref.has_intron(t.chrom, t.strand, a, b)
            types[key] = ("known_" if known else "novel_") + (
                "canonical" if canon[key] else "noncanonical")
    return canon, types


def _canonical(isoforms: Sequence[Transcript], calls: Sequence[Call], canon: Dict) -> Dict:
    by_cat: Dict[str, List[int]] = {}
    multi = allc = 0
    for t, c in zip(isoforms, calls):
        if not t.introns:
            continue
        ok = all(canon[(t.chrom, t.strand, a, b)] for a, b in t.introns)
        multi += 1
        allc += ok
        n = by_cat.setdefault(c.category, [0, 0])
        n[0] += 1
        n[1] += ok
    return {"multi_exon": multi, "all_canonical": allc,
            "by_category": {k: {"n": v[0], "all_canonical": v[1]} for k, v in sorted(by_cat.items())}}


def _rts(isoforms: Sequence[Transcript], fasta: Fasta, types: Dict) -> Dict:
    hit: Dict[Tuple[str, str, int, int], bool] = {}
    by_type = {k: {"n": 0, "rts": 0} for k in JUNCTION_TYPES}
    for key, typ in sorted(types.items()):
        chrom, strand, a, b = key
        if strand not in ("+", "-") or chrom not in fasta.index:
            continue
        x, y = rts_windows(fasta, chrom, strand, a, b)
        hit[key] = has_repeat(x, y)
        by_type[typ]["n"] += 1
        by_type[typ]["rts"] += hit[key]
    checked = [t for t in isoforms if t.introns and t.strand in ("+", "-") and t.chrom in fasta.index]
    with_rts = sum(1 for t in checked if any(hit.get((t.chrom, t.strand, a, b)) for a, b in t.introns))
    return {"junctions": len(hit), "rts": sum(hit.values()), "by_type": by_type,
            "isoforms": len(checked), "isoforms_with_rts": with_rts,
            "params": {"pattern": RTS_PATTERN, "wiggle": RTS_WIGGLE, "min_match": RTS_MIN_MATCH,
                       "mismatches": RTS_MISMATCHES}}


def _risk(r) -> bool:
    return r.n_a * 10 >= RISK_TENTHS * r.n_down


def _priming(isoforms: Sequence[Transcript], rows) -> Dict:
    exons = {t.tid: t.n_exons for t in isoforms}
    by_cat: Dict[str, List[int]] = {}
    by_exons = {"mono": [0, 0], "multi": [0, 0]}
    for r in rows:
        if not r.n_down:
            continue
        risk = _risk(r)
        v = by_cat.setdefault(r.category, [0, 0])
        v[0] += 1
        v[1] += risk
        g = by_exons["mono" if exons.get(r.isoform_id, 1) == 1 else "multi"]
        g[0] += 1
        g[1] += risk
    as_d = lambda d: {k: {"n": v[0], "risk": v[1]} for k, v in sorted(d.items())}
    return {"by_category": as_d(by_cat), "by_exons": as_d(by_exons)}


def _junction_support(isoforms: Sequence[Transcript], calls: Sequence[Call], ref: Reference,
                jf: Dict[Tuple[str, int, int], int], meta: Dict) -> Dict:
    uniq: Set[Tuple[str, str, int, int]] = {(t.chrom, t.strand, a, b)
                                            for t in isoforms for a, b in t.introns}
    known = {"n": 0, "supported": 0}
    novel = {"n": 0, "supported": 0}
    for chrom, strand, a, b in uniq:
        d = known if ref.has_intron(chrom, strand, a, b) else novel
        d["n"] += 1
        d["supported"] += jf.get((chrom, a, b), 0) > 0
    by_cat: Dict[str, List[int]] = {}
    multi = allsup = 0
    for t, c in zip(isoforms, calls):
        if not t.introns:
            continue
        ok = all(jf.get((t.chrom, a, b), 0) > 0 for a, b in t.introns)
        multi += 1
        allsup += ok
        v = by_cat.setdefault(c.category, [0, 0])
        v[0] += 1
        v[1] += ok
    return {"file": meta.get("file"), "md5": meta.get("md5"), "file_junctions": len(jf),
            "junctions": len(uniq), "supported": known["supported"] + novel["supported"],
            "known": known, "novel": novel,
            "multi_exon": multi, "all_supported": allsup,
            "by_category": {k: {"n": v[0], "all_supported": v[1]} for k, v in sorted(by_cat.items())}}


def coding_priming(ends_tsv: str, productivity_tsv: str) -> Optional[Dict]:
    """priming.by_coding from tables/ends.tsv and productivity.tsv, for sections.add_sections.
    Coding = a FLAIR productivity class with a start codon (PRO, PTC, NST); NGO is non-coding.
    Isoforms without a productivity row are left out."""
    from .concordance import read_header_tsv
    if not (os.path.exists(ends_tsv) and os.path.exists(productivity_tsv)):
        return None
    cls = {r["isoform_id"]: r["class"]
           for r in read_header_tsv(productivity_tsv, ["isoform_id", "class"])}
    out = {"coding": {"n": 0, "risk": 0}, "noncoding": {"n": 0, "risk": 0}}
    for r in read_header_tsv(ends_tsv, ["isoform_id", "down_a_frac"]):
        k = cls.get(r["isoform_id"])
        if not k or r["down_a_frac"] in ("", None):
            continue
        g = out["noncoding" if k == "NGO" else "coding"]
        g["n"] += 1
        g["risk"] += float(r["down_a_frac"]) >= 0.6 - 1e-9
    return out
