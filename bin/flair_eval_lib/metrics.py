"""Turn classified isoforms into the bundle.json numbers and the bundle tables."""

from __future__ import annotations

import bisect
import csv
import os
import re
from typing import Dict, List, Optional, Sequence

from . import CATEGORIES, FSM_SUBCATEGORIES, HISTOGRAM_EDGES
from .annotation import Regions, Transcript, merge_loci
from .classify import Call, Reference
from .fasta import Fasta

ISOFORM_COLUMNS = ["isoform_id", "gene_id", "chrom", "start", "end", "strand", "n_exons",
                   "length", "category", "subcategory", "ref_transcript", "ref_gene"]

# Categories whose ref_gene counts as "detected": the isoform shares splice sites (or, if
# mono-exon, structure) with that gene. genic/genic_intron only overlap it.
_DETECTING = {"FSM", "ISM", "NIC", "NNC", "fusion"}


def histogram(values: Sequence[float], edges: Sequence[float]) -> List[int]:
    """counts[i] = #values in [edges[i], edges[i+1]); the last bin is >= edges[-1].
    Values below edges[0] go in the first bin."""
    counts = [0] * len(edges)
    for v in values:
        counts[max(0, bisect.bisect_right(edges, v) - 1)] += 1
    return counts


def junction_counts(isoforms: Sequence[Transcript], ref: Reference,
                    fasta: Optional[Fasta]) -> Dict:
    """Distinct junctions (chrom, strand, intron start, intron end) across all isoforms."""
    uniq = sorted({(t.chrom, t.strand, a, b) for t in isoforms for a, b in t.introns})
    known = novel_c = novel_nc = 0
    for chrom, strand, a, b in uniq:
        if ref.has_intron(chrom, strand, a, b):
            known += 1
        elif fasta is not None:
            if fasta.is_canonical(chrom, strand, a, b):
                novel_c += 1
            else:
                novel_nc += 1
    novel = len(uniq) - known
    return {"total": len(uniq), "known": known,
            "novel_canonical": novel_c, "novel_noncanonical": novel_nc,
            "novel": novel, "motif_checked": fasta is not None}


def ref_in_scope(ref_tx: Sequence[Transcript], regions: Optional[Regions]) -> List[Transcript]:
    if regions is None:
        return list(ref_tx)
    return [t for t in ref_tx if regions.overlaps(t.chrom, t.start, t.end)]


def summarise(isoforms: Sequence[Transcript], calls: Sequence[Call],
              ref_scope: Sequence[Transcript], ref: Reference,
              fasta: Optional[Fasta]) -> Dict:
    categories = {c: 0 for c in CATEGORIES}
    subcats = {c: 0 for c in FSM_SUBCATEGORIES}
    fsm_refs, cim_refs, genes_hit = set(), set(), set()
    for c in calls:
        categories[c.category] += 1
        if c.category == "FSM":
            subcats[c.subcategory] += 1
            fsm_refs.add(c.ref_transcript)
            if c.subcategory == "CIM":
                cim_refs.add(c.ref_transcript)
        if c.category in _DETECTING and c.ref_gene:
            genes_hit.update(c.ref_gene.split(","))
    scope_ids = {t.tid for t in ref_scope}
    scope_genes = {t.gene for t in ref_scope}

    per_gene: Dict[str, int] = {}
    for t in isoforms:
        per_gene[t.gene] = per_gene.get(t.gene, 0) + 1
    mono = sum(1 for t in isoforms if t.n_exons == 1)

    hist = {
        "isoforms_per_gene": list(per_gene.values()),
        "transcript_length": [t.length for t in isoforms],
        "exon_length": [e - s for t in isoforms for s, e in t.exons],
        "exons_per_isoform": [t.n_exons for t in isoforms],
    }
    return {
        "summary": {
            "isoforms": len(isoforms),
            "genes": len(per_gene),
            "loci": len(merge_loci(isoforms)),
            "mono_exon": mono,
            "multi_exon": len(isoforms) - mono,
            "ref_transcripts_total": len(scope_ids),
            "ref_transcripts_fsm": len(fsm_refs & scope_ids),
            "ref_transcripts_cim": len(cim_refs & scope_ids),
            "ref_genes_detected": len(genes_hit & scope_genes),
            "ref_genes_total": len(scope_genes),
        },
        "categories": categories,
        "fsm_subcategories": subcats,
        "junctions": junction_counts(isoforms, ref, fasta),
        "histograms": {k: {"edges": list(HISTOGRAM_EDGES[k]),
                           "counts": histogram(v, HISTOGRAM_EDGES[k])}
                       for k, v in hist.items()},
    }


def write_gene_names(out_dir: str, calls: Sequence[Call], ref_tx: Sequence[Transcript]) -> None:
    """tables/genes.tsv: gene_id -> gene_name for every reference gene named in ref_gene.
    Optional extra table; the report uses it to label loci with readable names."""
    used = {g for c in calls for g in c.ref_gene.split(",") if g}
    names = {}
    for t in ref_tx:
        if t.gene in used and t.gene not in names:
            names[t.gene] = t.gene_name
    with open(os.path.join(out_dir, "tables", "genes.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene_id", "gene_name"])
        for g in sorted(names):
            w.writerow([g, names[g]])


def write_tables(out_dir: str, isoforms: Sequence[Transcript], calls: Sequence[Call]) -> None:
    tdir = os.path.join(out_dir, "tables")
    os.makedirs(tdir, exist_ok=True)
    with open(os.path.join(tdir, "isoforms.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(ISOFORM_COLUMNS)
        for t, c in zip(isoforms, calls):
            w.writerow([t.tid, t.gene, t.chrom, t.start, t.end, t.strand, t.n_exons, t.length,
                        c.category, c.subcategory, c.ref_transcript, c.ref_gene])
    counts: Dict = {}
    for c in calls:
        key = (c.category, c.subcategory)
        counts[key] = counts.get(key, 0) + 1
    with open(os.path.join(tdir, "categories.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["category", "subcategory", "n"])
        for cat in CATEGORIES:
            subs = FSM_SUBCATEGORIES if cat == "FSM" else ("",)
            for sub in subs:
                w.writerow([cat, sub, counts.get((cat, sub), 0)])


_VERSION_RE = re.compile(r"(\d+\.\d+(?:\.\d+)?(?:[.+\-][0-9A-Za-z.+\-]*)?)")


def normalise_version(raw: Optional[str]) -> str:
    """'FLAIR 3.0.1' -> '3.0.1'; 'flair 3.0.0+dev' -> '3.0.0+dev'; unknown -> raw or ''."""
    if not raw:
        return ""
    m = _VERSION_RE.search(raw)
    return m.group(1) if m else raw.strip()
