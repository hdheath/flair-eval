"""Isoform ends against end databases and the genome (CONTRACTS A27-A29, A31).

Per isoform (tables/ends.tsv, one row per isoform on a main chromosome):

  tss, tes     the isoform's 5' and 3' end positions, 0-based and strand-aware: on '+' the
               first and last base of the model (start, end - 1); on '-' the reverse
               (end - 1, start).
  d_reftss     signed distance from tss to the nearest same-strand refTSS site.
  d_pas        signed distance from tes to the nearest same-strand PolyA_DB site.
  d_ref_tss    signed distance from tss / tes to the 5' / 3' end of the matched reference
  d_ref_tes    transcript (FSM and ISM only; empty otherwise).
  polya_motif  the first hexamer of POLYA_HEXAMERS, in that priority order, found in the 50 nt
               of genome that end at the 3' end (transcript orientation), or empty.
  down_a_frac  the fraction of A in the 20 nt of genome just downstream of the 3' end
               (transcript orientation). >= 0.6 counts as internal-priming risk.

Distance sign (A28): genomic distance in transcript orientation, isoform end minus site, so
negative means the isoform end lies upstream (5') of the site. On '+' that is
`pos - site`; on '-' it is `site - pos`. A site that is an interval [s, e) is at distance 0
when s <= pos < e, otherwise the distance is to its nearer edge base (s or e - 1).
Distances are genomic, not spliced. For the 3' end this is SQANTI3's diff_to_TTS sign
(negative = shorter); for the 5' end it is the opposite of SQANTI3's diff_to_TSS, which is
positive for a 5' extension. A site farther than MAX_SITE_DIST is "no site" (empty).
The distances to the matched reference transcript have no cap.

Isoforms without a strand ('.') get a row with tss/tes from the + convention and every other
column empty; they are left out of the section counts. Only main chromosomes
(chr1-22, X, Y, with or without the 'chr' prefix) are evaluated.

Nearest-site lookups use per (chrom, strand) sorted start arrays, a running maximum of the
ends and bisect, so a genome-wide run costs O(isoforms x log sites).
"""

from __future__ import annotations

import bisect
import csv
import hashlib
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from . import CATEGORIES
from .annotation import InputError, Transcript, _open
from .classify import Call
from .fasta import Fasta

MAIN_CHROM = re.compile(r"^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$")

WINDOWS = [0, 10, 25, 50, 100, 200, 300, 500]
DIST_EDGES = [-1000, -200, -50, -10, 0, 10, 50, 200, 1000]
DOWN_A_EDGES = [i / 10 for i in range(11)]
MAX_SITE_DIST = 10_000
CATEGORY_WINDOW = 50

POLYA_HEXAMERS = ("AATAAA", "ATTAAA", "AGTAAA", "TATAAA", "CATAAA", "GATAAA",
                  "AATATA", "AATACA", "AATAGA", "AAAAAG", "ACTAAA", "AAGAAA")
MOTIF_WINDOW = 50
DOWN_A_WINDOW = 20
RISK_A_FRAC = 0.6
_RISK_TENTHS = 6       # RISK_A_FRAC in tenths, for integer comparison

ENDS_COLUMNS = ["isoform_id", "category", "subcategory", "strand", "tss", "tes",
                "d_reftss", "d_pas", "d_ref_tss", "d_ref_tes", "polya_motif", "down_a_frac"]

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def is_main_chrom(chrom: str) -> bool:
    return bool(MAIN_CHROM.match(chrom))


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def md5sum(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class SiteIndex:
    """Same-strand nearest-site lookup over a BED6 of sites (intervals or 1-bp points).

    Rows are deduplicated by (chrom, strand, start, end): PolyA_DB lists one cleavage site
    once per gene name. Sites on non-main chromosomes are dropped. A site with strand '.'
    is put on both strands.
    """

    def __init__(self, path: str):
        self.path = path
        self.file = os.path.basename(path)
        seen = set()
        by_key: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
        with _open(path) as fh:
            for n, line in enumerate(fh, 1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 6:
                    raise InputError(f"{path}:{n}: site BED needs 6 columns (BED6), got {len(p)}")
                try:
                    s, e = int(p[1]), int(p[2])
                except ValueError:
                    raise InputError(f"{path}:{n}: non-integer coordinates") from None
                if e <= s:
                    raise InputError(f"{path}:{n}: end <= start")
                chrom, strand = p[0], p[5]
                if not is_main_chrom(chrom):
                    continue
                key = (chrom, strand, s, e)
                if key in seen:
                    continue
                seen.add(key)
                for st in (("+", "-") if strand not in ("+", "-") else (strand,)):
                    by_key.setdefault((chrom, st), []).append((s, e))
        if not seen:
            raise InputError(f"{path}: no sites on main chromosomes")
        self.sites = len(seen)
        self._idx: Dict[Tuple[str, str], Tuple[List[int], List[int], List[int]]] = {}
        for key, iv in by_key.items():
            iv.sort()
            starts = [s for s, _ in iv]
            ends = [e for _, e in iv]
            run_max, m = [], -1
            for e in ends:
                m = max(m, e)
                run_max.append(m)
            self._idx[key] = (starts, ends, run_max)

    def meta(self) -> Dict:
        return {"file": self.file, "md5": md5sum(self.path), "sites": self.sites}

    def genomic_offset(self, chrom: str, strand: str, pos: int) -> Optional[int]:
        """pos minus the nearest site base on the genome (0 inside a site), or None."""
        idx = self._idx.get((chrom, strand))
        if idx is None:
            return None
        starts, _ends, run_max = idx
        i = bisect.bisect_right(starts, pos) - 1
        best = None
        if i >= 0:
            if run_max[i] > pos:
                return 0
            best = pos - (run_max[i] - 1)          # > 0: site lies left of pos
        if i + 1 < len(starts):
            right = pos - starts[i + 1]            # < 0: site lies right of pos
            if best is None or -right < best:
                best = right
        return best

    def distance(self, chrom: str, strand: str, pos: int) -> Optional[int]:
        """Signed distance in transcript orientation (negative = pos upstream of the site),
        or None when no same-strand site lies within MAX_SITE_DIST."""
        g = self.genomic_offset(chrom, strand, pos)
        if g is None or abs(g) > MAX_SITE_DIST:
            return None
        return g if strand == "+" else -g


def end_positions(t: Transcript, strand: Optional[str] = None) -> Tuple[int, int]:
    """(5' end, 3' end), 0-based base positions, for the given strand (default t.strand)."""
    s = strand or t.strand
    if s == "-":
        return t.end - 1, t.start
    return t.start, t.end - 1


def oriented(diff_genomic: int, strand: str) -> int:
    return diff_genomic if strand == "+" else -diff_genomic


def _fetch(fasta: Fasta, chrom: str, start: int, end: int) -> str:
    rec = fasta.index.get(chrom)
    if rec is None:
        return ""
    start, end = max(0, start), min(rec[0], end)
    return fasta.fetch(chrom, start, end) if end > start else ""


def upstream_of_3p(fasta: Fasta, t: Transcript, n: int = MOTIF_WINDOW) -> str:
    """The n nt of genome ending at the 3' end (inclusive), in transcript orientation."""
    if t.strand == "+":
        return _fetch(fasta, t.chrom, t.end - n, t.end)
    return revcomp(_fetch(fasta, t.chrom, t.start, t.start + n))


def downstream_of_3p(fasta: Fasta, t: Transcript, n: int = DOWN_A_WINDOW) -> str:
    """The n nt of genome just past the 3' end, in transcript orientation."""
    if t.strand == "+":
        return _fetch(fasta, t.chrom, t.end, t.end + n)
    return revcomp(_fetch(fasta, t.chrom, t.start - n, t.start))


def find_motif(seq: str) -> str:
    for hexamer in POLYA_HEXAMERS:
        if hexamer in seq:
            return hexamer
    return ""


@dataclass
class EndRow:
    isoform_id: str
    category: str
    subcategory: str
    strand: str
    tss: int
    tes: int
    d_reftss: Optional[int] = None
    d_pas: Optional[int] = None
    d_ref_tss: Optional[int] = None
    d_ref_tes: Optional[int] = None
    polya_motif: Optional[str] = None      # None = not checked, "" = none found
    n_a: Optional[int] = None              # A count in the downstream window
    n_down: Optional[int] = None           # bases read in the downstream window

    @property
    def stranded(self) -> bool:
        return self.strand in ("+", "-")

    @property
    def down_a_frac(self) -> Optional[float]:
        if not self.n_down:
            return None
        return self.n_a / self.n_down


def end_rows(isoforms: Sequence[Transcript], calls: Sequence[Call],
             ref_by_id: Dict[str, Transcript], tss_sites: Optional[SiteIndex] = None,
             pas_sites: Optional[SiteIndex] = None,
             fasta: Optional[Fasta] = None) -> List[EndRow]:
    """One EndRow per isoform on a main chromosome, in input order."""
    rows: List[EndRow] = []
    for t, c in zip(isoforms, calls):
        if not is_main_chrom(t.chrom):
            continue
        tss, tes = end_positions(t)
        r = EndRow(t.tid, c.category, c.subcategory, t.strand, tss, tes)
        rows.append(r)
        if c.category in ("FSM", "ISM") and c.ref_transcript in ref_by_id:
            ref = ref_by_id[c.ref_transcript]
            strand = t.strand if r.stranded else ref.strand
            if strand in ("+", "-"):
                r5, r3 = end_positions(ref, strand)
                i5, i3 = end_positions(t, strand)
                r.d_ref_tss = oriented(i5 - r5, strand)
                r.d_ref_tes = oriented(i3 - r3, strand)
        if not r.stranded:
            continue
        if tss_sites is not None:
            r.d_reftss = tss_sites.distance(t.chrom, t.strand, tss)
        if pas_sites is not None:
            r.d_pas = pas_sites.distance(t.chrom, t.strand, tes)
        if fasta is not None and t.chrom in fasta.index:
            r.polya_motif = find_motif(upstream_of_3p(fasta, t))
            down = downstream_of_3p(fasta, t)
            r.n_down = len(down)
            r.n_a = down.count("A")
    return rows


def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.4g}"
    return str(v)


def write_ends_table(out_dir: str, rows: Sequence[EndRow]) -> None:
    tdir = os.path.join(out_dir, "tables")
    os.makedirs(tdir, exist_ok=True)
    with open(os.path.join(tdir, "ends.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(ENDS_COLUMNS)
        for r in rows:
            w.writerow([r.isoform_id, r.category, r.subcategory, r.strand, r.tss, r.tes,
                        _fmt(r.d_reftss), _fmt(r.d_pas), _fmt(r.d_ref_tss), _fmt(r.d_ref_tes),
                        _fmt(r.polya_motif), _fmt(r.down_a_frac)])


# ------------------------------------------------------------------ sections
def dist_hist(values: Sequence[int]) -> Dict:
    """A29 bins: [edge_i, edge_i+1), plus open bins below the first and above the last edge,
    so len(counts) == len(DIST_EDGES) + 1. Exactly 0 falls in [0, 10)."""
    counts = [0] * (len(DIST_EDGES) + 1)
    for v in values:
        counts[bisect.bisect_right(DIST_EDGES, v)] += 1
    return {"edges": list(DIST_EDGES), "counts": counts}


def _half(rows: Sequence[EndRow], attr: str) -> Dict:
    stranded = [r for r in rows if r.stranded]
    vals = [getattr(r, attr) for r in stranded]
    have = [v for v in vals if v is not None]
    return {
        "n": len(stranded),
        "within": [sum(1 for v in have if abs(v) <= w) for w in WINDOWS],
        "dist_hist": dist_hist(have),
        "no_site": len(vals) - len(have),
    }


def ends_section(rows: Sequence[EndRow], tss_sites: Optional[SiteIndex],
                 pas_sites: Optional[SiteIndex], genome_checked: bool) -> Dict:
    """sections.ends (A29). Each half is None without its site file."""
    fsm = [r for r in rows if r.category == "FSM"]
    by_cat = {}
    for cat in CATEGORIES:
        rs = [r for r in rows if r.category == cat and r.stranded]
        by_cat[cat] = {
            "n": len(rs),
            "tss_within_50": (sum(1 for r in rs if r.d_reftss is not None
                                  and abs(r.d_reftss) <= CATEGORY_WINDOW)
                              if tss_sites is not None else None),
            "tes_within_50": (sum(1 for r in rs if r.d_pas is not None
                                  and abs(r.d_pas) <= CATEGORY_WINDOW)
                              if pas_sites is not None else None),
        }
    checked = [r for r in rows if r.polya_motif is not None]
    return {
        "refs": {"tss": tss_sites.meta() if tss_sites is not None else None,
                 "pas": pas_sites.meta() if pas_sites is not None else None},
        "windows": list(WINDOWS),
        "tss": _half(rows, "d_reftss") if tss_sites is not None else None,
        "tes": _half(rows, "d_pas") if pas_sites is not None else None,
        "fsm_tss_hist": dist_hist([r.d_ref_tss for r in fsm if r.d_ref_tss is not None]),
        "fsm_tes_hist": dist_hist([r.d_ref_tes for r in fsm if r.d_ref_tes is not None]),
        "by_category": by_cat,
        "polya_motif": ({"n_with_motif": sum(1 for r in checked if r.polya_motif), "n": len(checked)}
                        if genome_checked else None),
    }


def priming_section(rows: Sequence[EndRow]) -> Dict:
    """sections.priming (A31). Bins follow section 3: [edge_i, edge_i+1), the last bin is
    >= 1.0 (an all-A window), so len(counts) == len(DOWN_A_EDGES). Binning uses integer
    counts, so 0.6 is never lost to floating point."""
    have = [r for r in rows if r.n_down]
    counts = [0] * len(DOWN_A_EDGES)
    risk = 0
    for r in have:
        counts[min(10, (r.n_a * 10) // r.n_down)] += 1
        if r.n_a * 10 >= _RISK_TENTHS * r.n_down:
            risk += 1
    return {"n": len(have), "risk": risk,
            "down_a_hist": {"edges": list(DOWN_A_EDGES), "counts": counts}}
