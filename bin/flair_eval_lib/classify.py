"""SQANTI3-compatible structural classification of isoforms against a reference GTF.

Ported from legacy/bin/evaluation/flair_structural.py (which only had FSM/ISM/NIC/NNC and
mono-exon SEM/SEN) and extended to the full SQANTI3 category set.

Rules, in order. "Same strand" means the isoform's strand; an isoform without a strand
('.') is compared against both strands.

Multi-exon isoforms
  FSM          intron chain equals a reference transcript's chain.
  ISM          intron chain is a contiguous, shorter part of a reference transcript's chain.
  fusion       shares splice sites with two or more reference genes whose spans do not
               overlap each other (separate loci).
  NIC          every splice site (intron start and intron end) is annotated, but the chain
               is new: new combinations of known junctions/sites, or intron retention.
  NNC          uses at least one new splice site, and either shares an annotated site or
               overlaps exons of a same-strand spliced (multi-exon) reference transcript.
  genic        shares no splice site and overlaps same-strand exons only of mono-exon
               reference transcripts.

Mono-exon isoforms
  FSM          overlaps a same-strand mono-exon reference transcript.
  ISM          lies inside one exon of a multi-exon reference transcript.
  NIC          covers a whole same-strand reference intron (mono-exon by intron retention).
  genic        otherwise overlaps same-strand reference exons.

Then, for both, when no same-strand exon is overlapped:
  antisense    overlaps exons of a reference transcript on the opposite strand.
  genic_intron lies within the span of a reference transcript (either strand) without
               touching its exons, i.e. inside an intron.
  intergenic   overlaps no reference transcript.

These follow SQANTI3's priorities, checked against SQANTI3 calls on the fixture runs
(tests/fixtures/bundles/README.md). Known differences: mono-exon FSMs get CIM/FSM5/FSM3/FSMTE
here (SQANTI3 only says "mono-exon"), and SQANTI3 sometimes calls a sub-chain with a
retained intron NIC where this module says ISM.

FSM subcategories (end window 50 bp, the SQANTI3 default), scored against the
best-matching reference transcript with the same chain:
  CIM   both ends within the window   (SQANTI3 reference_match)
  FSM5  5' end novel, 3' end matches  (SQANTI3 alternative_5end)
  FSM3  3' end novel, 5' end matches  (SQANTI3 alternative_3end)
  FSMTE both ends novel               (SQANTI3 alternative_3end5end)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from .annotation import Transcript

FSM_END_WINDOW = 50
_BIN = 100_000


@dataclass
class Call:
    category: str
    subcategory: str = ""
    ref_transcript: str = ""
    ref_gene: str = ""


def _strands(strand: str) -> Tuple[str, ...]:
    return (strand,) if strand in ("+", "-") else ("+", "-")


def _opposite(strand: str) -> Tuple[str, ...]:
    if strand == "+":
        return ("-",)
    if strand == "-":
        return ("+",)
    return ()


class Reference:
    """Indexed reference annotation."""

    def __init__(self, transcripts: Iterable[Transcript]):
        self.tx: List[Transcript] = list(transcripts)
        self._bins: Dict[str, Dict[int, List[int]]] = {}
        self.chains: Dict[Tuple[str, str, Tuple], List[int]] = {}
        self.introns: Set[Tuple[str, str, int, int]] = set()
        self.left_sites: Dict[Tuple[str, str, int], Set[str]] = {}
        self.right_sites: Dict[Tuple[str, str, int], Set[str]] = {}
        self.gene_span: Dict[str, Tuple[str, int, int]] = {}
        self._intron_pos: List[Dict[Tuple[int, int], int]] = []
        for i, t in enumerate(self.tx):
            bins = self._bins.setdefault(t.chrom, {})
            for b in range(t.start // _BIN, (t.end - 1) // _BIN + 1):
                bins.setdefault(b, []).append(i)
            self._intron_pos.append({iv: k for k, iv in enumerate(t.introns)})
            if t.introns:
                self.chains.setdefault((t.chrom, t.strand, t.introns), []).append(i)
            for a, b in t.introns:
                self.introns.add((t.chrom, t.strand, a, b))
                self.left_sites.setdefault((t.chrom, t.strand, a), set()).add(t.gene)
                self.right_sites.setdefault((t.chrom, t.strand, b), set()).add(t.gene)
            c, s, e = self.gene_span.get(t.gene, (t.chrom, t.start, t.end))
            self.gene_span[t.gene] = (c, min(s, t.start), max(e, t.end))

    def overlapping(self, chrom: str, start: int, end: int) -> List[int]:
        """Indices of reference transcripts whose span overlaps [start, end)."""
        bins = self._bins.get(chrom)
        if not bins:
            return []
        hits = set()
        for b in range(start // _BIN, (end - 1) // _BIN + 1):
            for i in bins.get(b, ()):
                t = self.tx[i]
                if t.start < end and start < t.end:
                    hits.add(i)
        return sorted(hits)

    def has_intron(self, chrom: str, strand: str, a: int, b: int) -> bool:
        return any((chrom, s, a, b) in self.introns for s in _strands(strand))

    def is_sub_chain(self, i: int, introns: Tuple) -> bool:
        """True when `introns` is a contiguous run of reference transcript i's introns."""
        k = self._intron_pos[i].get(introns[0])
        if k is None:
            return False
        ref = self.tx[i].introns
        return ref[k:k + len(introns)] == introns


def fsm_subcategory(start: int, end: int, strand: str,
                    candidates: Sequence[Transcript],
                    window: int = FSM_END_WINDOW) -> Tuple[str, Optional[Transcript]]:
    """Best FSM subcategory and the reference transcript it was scored against.

    More matched ends wins, then the smaller total end distance, then the transcript id
    (so the label never depends on GTF order).
    """
    best = None
    for t in sorted(candidates, key=lambda t: (t.start, t.end, t.tid)):
        d_start, d_end = abs(start - t.start), abs(end - t.end)
        d5, d3 = (d_end, d_start) if strand == "-" else (d_start, d_end)
        m5, m3 = d5 <= window, d3 <= window
        rank = (m5 + m3, -(d5 + d3))
        if best is None or rank > best[0]:
            best = (rank, m5, m3, t)
    if best is None:
        return "", None
    _, m5, m3, t = best
    if m5 and m3:
        return "CIM", t
    if m3:
        return "FSM5", t
    if m5:
        return "FSM3", t
    return "FSMTE", t


def _exon_overlap(a: Sequence[Tuple[int, int]], b: Sequence[Tuple[int, int]]) -> int:
    i = j = tot = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if lo < hi:
            tot += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def _pick(cands: Iterable[Transcript], iso: Transcript) -> Transcript:
    """Deterministic choice among candidate references: closest length, then id."""
    return min(cands, key=lambda t: (abs(t.length - iso.length), t.tid))


class Classifier:
    def __init__(self, reference: Reference, end_window: int = FSM_END_WINDOW):
        self.ref = reference
        self.window = end_window

    def classify(self, iso: Transcript) -> Call:
        idx = self.ref.overlapping(iso.chrom, iso.start, iso.end)
        tx = self.ref.tx
        same = [tx[i] for i in idx if tx[i].strand in _strands(iso.strand)]
        opp = [tx[i] for i in idx if tx[i].strand in _opposite(iso.strand)]
        if iso.n_exons > 1:
            call = self._multi(iso, idx, same)
        else:
            call = self._mono(iso, same)
        return call or self._no_shared_site(iso, same, opp)

    # -- multi-exon ------------------------------------------------------------------
    def _multi(self, iso: Transcript, idx: List[int], same: List[Transcript]) -> Optional[Call]:
        ref = self.ref
        fsm: List[Transcript] = []
        for s in _strands(iso.strand):
            fsm += [ref.tx[i] for i in ref.chains.get((iso.chrom, s, iso.introns), ())]
        if fsm:
            sub, best = fsm_subcategory(iso.start, iso.end, iso.strand, fsm, self.window)
            return Call("FSM", sub, best.tid, best.gene)

        ism = [ref.tx[i] for i in idx
               if ref.tx[i].strand in _strands(iso.strand)
               and len(ref.tx[i].introns) > len(iso.introns)
               and ref.is_sub_chain(i, iso.introns)]
        if ism:
            best = _pick(ism, iso)
            return Call("ISM", "", best.tid, best.gene)

        genes: Dict[str, int] = {}
        all_known = True
        for a, b in iso.introns:
            left: Set[str] = set()
            right: Set[str] = set()
            for s in _strands(iso.strand):
                left |= ref.left_sites.get((iso.chrom, s, a), set())
                right |= ref.right_sites.get((iso.chrom, s, b), set())
            if not left or not right:
                all_known = False
            for g in left | right:
                genes[g] = genes.get(g, 0) + 1
        if not genes:
            return None

        loci = self._separate_loci(genes)
        if len(loci) > 1:
            names = sorted(max(locus, key=lambda g: (genes[g], g)) for locus in loci)
            return Call("fusion", "", "", ",".join(names))
        gene = max(genes, key=lambda g: (genes[g], g))
        return Call("NIC" if all_known else "NNC", "", "", gene)

    def _separate_loci(self, genes: Dict[str, int]) -> List[List[str]]:
        spans = sorted((self.ref.gene_span[g][0], self.ref.gene_span[g][1],
                        self.ref.gene_span[g][2], g) for g in genes)
        loci: List[List[str]] = []
        cur_chrom, cur_end = None, -1
        for c, s, e, g in spans:
            if loci and c == cur_chrom and s < cur_end:
                loci[-1].append(g)
                cur_end = max(cur_end, e)
            else:
                loci.append([g])
                cur_chrom, cur_end = c, e
        return loci

    # -- mono-exon -------------------------------------------------------------------
    def _mono(self, iso: Transcript, same: List[Transcript]) -> Optional[Call]:
        s, e = iso.start, iso.end
        mono = [t for t in same if t.n_exons == 1]
        if mono:
            sub, best = fsm_subcategory(s, e, iso.strand, mono, self.window)
            return Call("FSM", sub, best.tid, best.gene)
        multi = [t for t in same if t.n_exons > 1]
        inside = [t for t in multi if any(es <= s and e <= ee for es, ee in t.exons)]
        if inside:
            best = _pick(inside, iso)
            return Call("ISM", "", best.tid, best.gene)
        retained = [t for t in multi if any(s <= a and b <= e for a, b in t.introns)]
        if retained:
            best = _pick(retained, iso)
            return Call("NIC", "", "", best.gene)
        return None

    # -- shared fallback -------------------------------------------------------------
    @staticmethod
    def _no_shared_site(iso: Transcript, same: List[Transcript],
                        opp: List[Transcript]) -> Call:
        def best_by_overlap(cands):
            scored = [(_exon_overlap(iso.exons, t.exons), t) for t in cands]
            scored = [(o, t) for o, t in scored if o > 0]
            if not scored:
                return None
            return max(scored, key=lambda x: (x[0], -abs(x[1].length - iso.length),
                                              x[1].tid))[1]

        if iso.n_exons > 1:
            spliced = best_by_overlap([t for t in same if t.n_exons > 1])
            if spliced is not None:
                return Call("NNC", "", "", spliced.gene)
        hit = best_by_overlap(same)
        if hit is not None:
            return Call("genic", "", "", hit.gene)
        if best_by_overlap(opp) is not None:
            return Call("antisense")
        if same:
            return Call("genic_intron", "", "", _pick(same, iso).gene)
        if opp:
            return Call("genic_intron")
        return Call("intergenic")
