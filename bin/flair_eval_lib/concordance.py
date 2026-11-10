"""Run-vs-run concordance and expression agreement for one sample (CONTRACTS.md section 8).

Isoforms from different runs are matched by structure:
- multi-exon: identical intron chain (chrom, strand, introns);
- mono-exon: same chrom and strand, overlapping by at least 50% of *both* lengths. Mono-exon
  models are clustered across all runs of the sample (single linkage), and a cluster counts
  as one structure.

A run's structure set is the set of distinct structures among its evaluated isoforms (the rows
of tables/isoforms.tsv), so two FLAIR isoforms that differ only in their ends count once.
"""

from __future__ import annotations

import math
import os
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .annotation import Transcript, keep_contig, merge_loci, read_bed12

MONO_MIN_OVERLAP = 0.5
TOP_LOCI = 20
MAX_GENE_NAMES = 40
GRID_BINS = 20
_DETECTING = {"FSM", "ISM", "NIC", "NNC", "fusion"}


# ---------------------------------------------------------------- reading
def read_header_tsv(path: str, required: Sequence[str]) -> List[Dict[str, str]]:
    """Rows of a TSV with one header line (CONTRACTS A10). Raises ValueError if columns miss."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        missing = [c for c in required if c not in header]
        if missing:
            raise ValueError(f"{path}: header lacks column(s) {', '.join(missing)}")
        rows = []
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, parts)))
    return rows


def read_gene_names(bundle_dir: str) -> Dict[str, str]:
    path = os.path.join(bundle_dir, "tables", "genes.tsv")
    if not os.path.exists(path):
        return {}
    try:
        return {r["gene_id"]: r.get("gene_name", "") for r in read_header_tsv(path, ["gene_id"])}
    except ValueError:
        return {}


class RunModels:
    """One run's evaluated isoforms plus their table annotations."""

    def __init__(self, run_id: str, isoforms: List[Transcript], info: Dict[str, Tuple[str, str]]):
        self.run_id = run_id
        self.isoforms = isoforms
        self.info = info            # isoform_id -> (category, ref_gene)

    @classmethod
    def from_bundle_dir(cls, run_id: str, bundle_dir: str) -> Optional["RunModels"]:
        bed = os.path.join(bundle_dir, "isoforms.bed")
        if not os.path.exists(bed):
            return None
        isoforms = read_bed12(bed)
        info: Dict[str, Tuple[str, str]] = {}
        table = os.path.join(bundle_dir, "tables", "isoforms.tsv")
        if os.path.exists(table):
            for r in read_header_tsv(table, ["isoform_id", "category", "ref_gene"]):
                info[r["isoform_id"]] = (r["category"], r["ref_gene"])
            isoforms = [t for t in isoforms if t.tid in info]
        else:
            isoforms = [t for t in isoforms if keep_contig(t.chrom)]
        return cls(run_id, isoforms, info)

    def detected_genes(self) -> set:
        out = set()
        for cat, gene in self.info.values():
            if cat in _DETECTING and gene:
                out.update(g for g in gene.split(",") if g)
        return out


# ---------------------------------------------------------------- structures
def assign_structures(runs: Sequence[RunModels]) -> Dict[Tuple[str, str], int]:
    """(run_id, isoform_id) -> structure id, shared across runs."""
    sid: Dict[Tuple[str, str], int] = {}
    chains: Dict[Tuple, int] = {}
    monos: Dict[Tuple[str, str], List[Tuple[int, int, str, str]]] = {}
    for r in runs:
        for t in r.isoforms:
            if t.introns:
                key = (t.chrom, t.strand, t.introns)
                if key not in chains:
                    chains[key] = len(chains)
                sid[(r.run_id, t.tid)] = chains[key]
            else:
                monos.setdefault((t.chrom, t.strand), []).append((t.start, t.end, r.run_id, t.tid))
    next_id = len(chains)
    for key in sorted(monos):
        items = sorted(monos[key])
        parent = list(range(len(items)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        active: List[int] = []
        for i, (s, e, _, _) in enumerate(items):
            active = [j for j in active if items[j][1] > s]
            for j in active:
                js, je = items[j][0], items[j][1]
                ov = min(e, je) - max(s, js)
                if ov >= MONO_MIN_OVERLAP * (e - s) and ov >= MONO_MIN_OVERLAP * (je - js):
                    parent[find(i)] = find(j)
            active.append(i)
        roots: Dict[int, int] = {}
        for i, (_, _, run_id, tid) in enumerate(items):
            root = find(i)
            if root not in roots:
                roots[root] = next_id
                next_id += 1
            sid[(run_id, tid)] = roots[root]
    return sid


# ---------------------------------------------------------------- stats
def _ranks(v: Sequence[float]) -> List[float]:
    order = sorted(range(len(v)), key=lambda i: v[i])
    ranks = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def pearson(x: Sequence[float], y: Sequence[float]) -> Optional[float]:
    n = len(x)
    if n < 3:
        return None
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx == 0 or syy == 0:
        return None
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(sxx * syy)


def spearman(x: Sequence[float], y: Sequence[float]) -> Optional[float]:
    if len(x) < 3:
        return None
    return pearson(_ranks(x), _ranks(y))


def _r(v: Optional[float]) -> Optional[float]:
    return None if v is None else round(v, 4)


# ---------------------------------------------------------------- concordance
def concordance(runs: Sequence[RunModels], baseline: Optional[str],
                gene_names: Dict[str, str]) -> Optional[dict]:
    """Section-8 concordance for one sample. None when fewer than two runs have isoforms.bed."""
    runs = [r for r in runs if r is not None]
    if len(runs) < 2:
        return None
    sid = assign_structures(runs)
    sets = {r.run_id: {sid[(r.run_id, t.tid)] for t in r.isoforms} for r in runs}
    ids = [r.run_id for r in runs]
    union_all = set().union(*sets.values())
    shared_all = set.intersection(*sets.values())

    pairs = []
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            inter = len(sets[a] & sets[b])
            uni = len(sets[a] | sets[b])
            pairs.append({"a": a, "b": b, "shared": inter, "only_a": len(sets[a]) - inter,
                          "only_b": len(sets[b]) - inter,
                          "jaccard": round(inter / uni, 4) if uni else None})

    def label(g):
        return gene_names.get(g) or g

    genes = None
    if baseline in sets:
        base_genes = next(r for r in runs if r.run_id == baseline).detected_genes()
        genes = {}
        for r in runs:
            if r.run_id == baseline:
                continue
            g = r.detected_genes()
            gained, lost = sorted(g - base_genes, key=label), sorted(base_genes - g, key=label)
            genes[r.run_id] = {"gained": len(gained), "lost": len(lost),
                               "gained_names": [label(x) for x in gained[:MAX_GENE_NAMES]],
                               "lost_names": [label(x) for x in lost[:MAX_GENE_NAMES]]}

    return {
        "runs": ids,
        "structures": {k: len(v) for k, v in sets.items()},
        "union": len(union_all),
        "shared_by_all": len(shared_all),
        "pairs": pairs,
        "genes_vs_baseline": genes,
        "baseline": baseline if baseline in sets else None,
        "top_loci": top_loci(runs, sid, gene_names),
        "mono_min_overlap": MONO_MIN_OVERLAP,
    }


def top_loci(runs: Sequence[RunModels], sid: Dict, gene_names: Dict[str, str],
             n: int = TOP_LOCI) -> List[dict]:
    """Loci (same-strand overlapping isoforms, pooled over runs) ranked by disagreement:
    structures not found in every run, then the spread of isoform counts."""
    tagged = []
    owner = {}
    for r in runs:
        for t in r.isoforms:
            u = Transcript(f"{r.run_id}\t{t.tid}", t.gene, t.chrom, t.strand, t.exons)
            owner[u.tid] = (r, t)
            tagged.append(u)
    ids = [r.run_id for r in runs]
    out = []
    for locus in merge_loci(tagged):
        counts = {k: 0 for k in ids}
        per_run: Dict[str, set] = {k: set() for k in ids}
        genes: Dict[str, int] = {}
        for u in locus:
            r, t = owner[u.tid]
            counts[r.run_id] += 1
            per_run[r.run_id].add(sid[(r.run_id, t.tid)])
            cat, g = r.info.get(t.tid, ("", ""))
            if g and cat in _DETECTING | {"genic"}:
                for gg in g.split(","):
                    genes[gg] = genes.get(gg, 0) + 1
        allset = set().union(*per_run.values())
        common = set.intersection(*per_run.values())
        divergent = len(allset) - len(common)
        if divergent == 0:
            continue
        spread = max(counts.values()) - min(counts.values())
        top_genes = sorted(genes, key=lambda g: (-genes[g], g))[:3]
        out.append({
            "chrom": locus[0].chrom, "start": min(u.start for u in locus),
            "end": max(u.end for u in locus), "strand": locus[0].strand,
            "genes": [gene_names.get(g) or g for g in top_genes],
            "structures": len(allset), "shared_by_all": len(common),
            "counts": [counts[k] for k in ids],
            "_rank": (-divergent, -spread, -len(locus), locus[0].chrom, min(u.start for u in locus)),
        })
    out.sort(key=lambda d: d["_rank"])
    for d in out:
        del d["_rank"]
    return out[:n]


# ---------------------------------------------------------------- expression
def read_counts(path: str) -> Dict[str, float]:
    """isoform_id -> CPM from counts.tsv (header: isoform_id, count, cpm)."""
    rows = read_header_tsv(path, ["isoform_id", "cpm"])
    out: Dict[str, float] = {}
    for r in rows:
        try:
            out[r["isoform_id"]] = out.get(r["isoform_id"], 0.0) + float(r["cpm"])
        except ValueError:
            continue
    return out


def quant_agreement(runs: Sequence[RunModels], cpm: Dict[str, Dict[str, float]],
                    baseline: Optional[str]) -> Optional[dict]:
    """CPM agreement between runs, matching isoforms by structure (summed per structure).
    Correlations use log10(CPM + 1) over structures present in both runs."""
    runs = [r for r in runs if r is not None and r.run_id in cpm]
    if len(runs) < 2:
        return None
    sid = assign_structures(runs)
    per: Dict[str, Dict[int, float]] = {}
    for r in runs:
        d: Dict[int, float] = {}
        c = cpm[r.run_id]
        for t in r.isoforms:
            s = sid[(r.run_id, t.tid)]
            d[s] = d.get(s, 0.0) + c.get(t.tid, 0.0)
        per[r.run_id] = d
    ids = [r.run_id for r in runs]
    lg = {k: {s: math.log10(v + 1) for s, v in d.items()} for k, d in per.items()}
    top = max((v for d in lg.values() for v in d.values()), default=1.0) or 1.0
    top = math.ceil(top * 2) / 2
    pairs = []
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            common = sorted(set(per[a]) & set(per[b]))
            x = [lg[a][s] for s in common]
            y = [lg[b][s] for s in common]
            tot_a = sum(per[a].values()) or 1
            tot_b = sum(per[b].values()) or 1
            pair = {"a": a, "b": b, "n": len(common), "pearson": _r(pearson(x, y)),
                    "spearman": _r(spearman(x, y)),
                    "cpm_shared_a": round(sum(per[a][s] for s in common) / tot_a, 4),
                    "cpm_shared_b": round(sum(per[b][s] for s in common) / tot_b, 4)}
            if baseline in (a, b):
                # density grid: x = baseline, y = the other run
                bx, by = (x, y) if a == baseline else (y, x)
                grid = [[0] * GRID_BINS for _ in range(GRID_BINS)]
                for u, v in zip(bx, by):
                    i_ = min(GRID_BINS - 1, int(u / top * GRID_BINS))
                    j_ = min(GRID_BINS - 1, int(v / top * GRID_BINS))
                    grid[j_][i_] += 1
                pair["grid"] = grid
            pairs.append(pair)
    return {"runs": ids, "baseline": baseline if baseline in ids else None,
            "log_max": top, "grid_bins": GRID_BINS, "pairs": pairs}
