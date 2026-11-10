"""Read isoforms (BED12 / GTF), reference annotations (GTF) and region BEDs.

Coordinates are held in BED convention everywhere: 0-based, half-open.
GTF exons (1-based, inclusive) are converted on read: start - 1, end unchanged.
An intron is (end of the exon before it, start of the exon after it).
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

Exon = Tuple[int, int]

_ATTR_RE = {
    "gene_id": re.compile(r'gene_id "([^"]*)"'),
    "transcript_id": re.compile(r'transcript_id "([^"]*)"'),
}

# Contigs that are never evaluated: mitochondria, unplaced/unlocalised scaffolds, alts,
# decoys and patches. Everything else (chr1..22/X/Y, spike-in contigs such as SIRV1 or
# ERCC-00002) is kept, so the SIRV test set still works.
_DROP_CONTIG = re.compile(
    r"^(chr)?(M|MT)$|_|^(chr)?Un|^(GL|KI|JH|KN|KQ|KV|KZ|ML|MU)\d|^HLA-|^chrEBV$|decoy",
    re.IGNORECASE,
)


def keep_contig(chrom: str) -> bool:
    """True for contigs that are evaluated (main chromosomes and spike-in contigs)."""
    return not _DROP_CONTIG.search(chrom)


def _open(path: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


@dataclass
class Transcript:
    """One transcript model: a reference transcript or a query isoform."""

    tid: str
    gene: str
    chrom: str
    strand: str
    exons: Tuple[Exon, ...]
    introns: Tuple[Exon, ...] = field(init=False)

    def __post_init__(self):
        self.exons = tuple(sorted(self.exons))
        self.introns = tuple((self.exons[i][1], self.exons[i + 1][0])
                             for i in range(len(self.exons) - 1))

    @property
    def start(self) -> int:
        return self.exons[0][0]

    @property
    def end(self) -> int:
        return self.exons[-1][1]

    @property
    def n_exons(self) -> int:
        return len(self.exons)

    @property
    def length(self) -> int:
        return sum(e - s for s, e in self.exons)


class InputError(ValueError):
    """Bad input file content. CLIs turn this into a non-zero exit with a message."""


def _gtf_attr(attrs: str, key: str) -> Optional[str]:
    m = _ATTR_RE[key].search(attrs)
    return m.group(1) if m else None


def read_gtf_transcripts(path: str) -> List[Transcript]:
    """Every transcript in a GTF, from its `exon` lines. Order = first appearance."""
    exons: Dict[str, List[Exon]] = {}
    meta: Dict[str, Tuple[str, str, str]] = {}
    with _open(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                raise InputError(f"{path}:{n}: GTF line has {len(parts)} columns, expected 9")
            if parts[2] != "exon":
                continue
            tid = _gtf_attr(parts[8], "transcript_id")
            if tid is None:
                raise InputError(f"{path}:{n}: exon line has no transcript_id")
            gene = _gtf_attr(parts[8], "gene_id") or ""
            try:
                s, e = int(parts[3]) - 1, int(parts[4])
            except ValueError:
                raise InputError(f"{path}:{n}: non-integer coordinates") from None
            if tid not in meta:
                meta[tid] = (parts[0], parts[6], gene)
                exons[tid] = []
            elif meta[tid][0] != parts[0]:
                # Same transcript_id on two contigs (seen in some FLAIR builds for
                # identical mono-exon models): keep them apart.
                tid = f"{tid}__{parts[0]}"
                if tid not in meta:
                    meta[tid] = (parts[0], parts[6], gene)
                    exons[tid] = []
            exons[tid].append((s, e))
    return [Transcript(t, meta[t][2], meta[t][0], meta[t][1], tuple(_merge_touching(exons[t])))
            for t in meta]


def _merge_touching(exons: List[Exon]) -> List[Exon]:
    """Sort exons and merge overlapping/abutting blocks (0-bp introns are not introns)."""
    out: List[List[int]] = []
    for s, e in sorted(exons):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def read_gtf_gene_map(path: str) -> Dict[str, str]:
    """transcript_id -> gene_id from a GTF (any feature type)."""
    out: Dict[str, str] = {}
    with _open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            tid = _gtf_attr(parts[8], "transcript_id")
            if tid and tid not in out:
                out[tid] = _gtf_attr(parts[8], "gene_id") or ""
    return out


def read_bed12(path: str) -> List[Transcript]:
    """Isoforms from a BED12 (BED6 rows are read as mono-exon). gene is left empty."""
    out: List[Transcript] = []
    seen = set()
    with _open(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                raise InputError(f"{path}:{n}: BED line has {len(p)} columns, need >= 6")
            try:
                start, end = int(p[1]), int(p[2])
                if len(p) >= 12 and p[10].strip(",") and p[11].strip(","):
                    sizes = [int(x) for x in p[10].strip(",").split(",")]
                    starts = [int(x) for x in p[11].strip(",").split(",")]
                    if len(sizes) != len(starts):
                        raise InputError(f"{path}:{n}: blockSizes/blockStarts length differ")
                    exons = [(start + b, start + b + z) for b, z in zip(starts, sizes)]
                else:
                    exons = [(start, end)]
            except ValueError:
                raise InputError(f"{path}:{n}: non-integer BED field") from None
            if end <= start:
                raise InputError(f"{path}:{n}: end <= start")
            name = p[3]
            if name in seen:
                raise InputError(f"{path}:{n}: duplicate isoform name {name!r}")
            seen.add(name)
            strand = p[5] if p[5] in ("+", "-") else "."
            out.append(Transcript(name, "", p[0], strand, tuple(_merge_touching(exons))))
    return out


def assign_isoform_genes(isoforms: List[Transcript], gene_map: Dict[str, str]) -> None:
    """Fill isoform.gene from the isoforms GTF.

    The BED name is the isoform id. Older FLAIR writes BED names `<transcript>_<gene>`
    while its GTF transcript_id is only `<transcript>`, so fall back to that split.
    Isoforms still without a gene get `novel_gene_<n>`, one per overlapping same-strand
    cluster, as CONTRACTS.md section 2 asks of the adapter.
    """
    missing: List[Transcript] = []
    for iso in isoforms:
        g = gene_map.get(iso.tid)
        if not g and "_" in iso.tid:
            head, tail = iso.tid.rsplit("_", 1)
            g = gene_map.get(head) or (tail if tail else None)
        if g:
            iso.gene = g
        else:
            missing.append(iso)
    n = 0
    for cluster in merge_loci(missing):
        n += 1
        for iso in cluster:
            iso.gene = f"novel_gene_{n}"


def merge_loci(models: Iterable[Transcript]) -> List[List[Transcript]]:
    """Group models into loci: same chrom + strand, transitively overlapping spans."""
    by_key: Dict[Tuple[str, str], List[Transcript]] = {}
    for m in models:
        by_key.setdefault((m.chrom, m.strand), []).append(m)
    loci: List[List[Transcript]] = []
    for key in sorted(by_key):
        cur: List[Transcript] = []
        cur_end = -1
        for m in sorted(by_key[key], key=lambda t: (t.start, t.end, t.tid)):
            if cur and m.start < cur_end:
                cur.append(m)
                cur_end = max(cur_end, m.end)
            else:
                if cur:
                    loci.append(cur)
                cur, cur_end = [m], m.end
        if cur:
            loci.append(cur)
    return loci


class Regions:
    """A set of BED intervals. `overlaps` is true when a span touches any interval."""

    def __init__(self, intervals: Dict[str, List[Exon]], source: str = ""):
        self.source = source
        self.by_chrom = {c: _merge_touching(v) for c, v in intervals.items()}

    @classmethod
    def from_bed(cls, path: str) -> "Regions":
        iv: Dict[str, List[Exon]] = {}
        with _open(path) as fh:
            for n, line in enumerate(fh, 1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                p = line.split()
                if len(p) < 3:
                    raise InputError(f"{path}:{n}: regions BED needs chrom, start, end")
                try:
                    iv.setdefault(p[0], []).append((int(p[1]), int(p[2])))
                except ValueError:
                    raise InputError(f"{path}:{n}: non-integer coordinates") from None
        if not iv:
            raise InputError(f"{path}: regions BED is empty")
        return cls(iv, path)

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        for s, e in self.by_chrom.get(chrom, ()):
            if s < end and start < e:
                return True
        return False

    def as_list(self) -> List[List]:
        return [[c, s, e] for c in sorted(self.by_chrom) for s, e in self.by_chrom[c]]


def iter_tsv(path: str) -> Iterator[List[str]]:
    with _open(path) as fh:
        for line in fh:
            if line.strip():
                yield line.rstrip("\n").split("\t")
