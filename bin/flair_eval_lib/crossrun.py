"""Cross-run comparisons computed by build_report.py over one sample's bundles (CONTRACTS A33).

Nothing here is stored per bundle, and nothing depends on which run the reader picks as the
baseline: every result is given for every run pair (or every ordered pair), so the page can
switch the baseline (A22).

- `end_agreement`: Jaccard of 5' (TSS) and 3' (TTS) end clusters between two runs. The ends of
  both runs' evaluated isoforms are pooled per chrom and strand and clustered by single linkage
  (neighbouring ends at most END_WINDOW bp apart join one cluster). A cluster is shared when
  both runs have an end in it; Jaccard = shared / all clusters. Isoforms without a strand are
  left out. Clustering is done per pair, so a third run cannot bridge two clusters.
- `tusco_incremental`: from `sections.tusco.per_gene`. A gene is *found* when its call is TP or
  PTP. For a (baseline, run) pair, `gained` = found in the run and FN in the baseline (A33),
  `lost` = found in the baseline and FN in the run, `both` = found in both.
- `lost_by_decile`: the baseline's detected reference genes (as in concordance) ranked by the
  baseline's gene CPM from counts.tsv (sum of its isoforms' CPM; an isoform assigned to several
  genes splits its CPM evenly), cut into 10 equal-count deciles (ties by gene id), and per other
  run how many genes of each decile it does not detect.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

from .concordance import _DETECTING, RunModels, gene_label

END_WINDOW = 50
DECILES = 10
FOUND = ("TP", "PTP")


# ---------------------------------------------------------------- end clusters
def run_ends(run: RunModels) -> Dict[str, Dict[Tuple[str, str], List[int]]]:
    """{'tss'|'tes': {(chrom, strand): [positions]}} of a run's isoforms (0-based, strand-aware:
    + strand TSS = start, TES = end - 1; - strand TSS = end - 1, TES = start)."""
    out: Dict[str, Dict[Tuple[str, str], List[int]]] = {"tss": {}, "tes": {}}
    for t in run.isoforms:
        if t.strand not in ("+", "-"):
            continue
        first, last = t.start, t.end - 1
        tss, tes = (first, last) if t.strand == "+" else (last, first)
        out["tss"].setdefault((t.chrom, t.strand), []).append(tss)
        out["tes"].setdefault((t.chrom, t.strand), []).append(tes)
    return out


def cluster_jaccard(a: Dict[Tuple[str, str], List[int]], b: Dict[Tuple[str, str], List[int]],
                    window: int = END_WINDOW) -> dict:
    """Shared / only-a / only-b end clusters and their Jaccard (single linkage within `window`)."""
    shared = only_a = only_b = 0
    for key in set(a) | set(b):
        pts = sorted([(p, 0) for p in a.get(key, ())] + [(p, 1) for p in b.get(key, ())])
        i = 0
        while i < len(pts):
            has = [False, False]
            j = i
            while True:
                has[pts[j][1]] = True
                if j + 1 < len(pts) and pts[j + 1][0] - pts[j][0] <= window:
                    j += 1
                else:
                    break
            if has[0] and has[1]:
                shared += 1
            elif has[0]:
                only_a += 1
            else:
                only_b += 1
            i = j + 1
    total = shared + only_a + only_b
    return {"shared": shared, "only_a": only_a, "only_b": only_b,
            "jaccard": round(shared / total, 4) if total else None}


def end_agreement(runs: Sequence[Optional[RunModels]], window: int = END_WINDOW) -> Optional[dict]:
    """TSS/TTS cluster Jaccard for every run pair. None when fewer than two runs have isoforms."""
    runs = [r for r in runs if r is not None]
    if len(runs) < 2:
        return None
    ends = {r.run_id: run_ends(r) for r in runs}
    ids = [r.run_id for r in runs]
    pairs = []
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            pairs.append({"a": a, "b": b,
                          "tss": cluster_jaccard(ends[a]["tss"], ends[b]["tss"], window),
                          "tes": cluster_jaccard(ends[a]["tes"], ends[b]["tes"], window)})
    return {"runs": ids, "window": window, "pairs": pairs}


# ---------------------------------------------------------------- TUSCO
def tusco_incremental(per_gene: Dict[str, Dict[str, str]], gene_names: Dict[str, str]) -> Optional[dict]:
    """per_gene: {run_id: {gene_id: TP|PTP|FP|FN}} for runs with sections.tusco. Every ordered
    (baseline, run) pair, over the genes both runs scored."""
    ids = [r for r, v in per_gene.items() if isinstance(v, dict) and v]
    if len(ids) < 2:
        return None
    pairs = []
    for base in ids:
        for run in ids:
            if run == base:
                continue
            pb, pr = per_gene[base], per_gene[run]
            genes = sorted(set(pb) & set(pr))
            gained = [g for g in genes if pr[g] in FOUND and pb[g] == "FN"]
            lost = [g for g in genes if pb[g] in FOUND and pr[g] == "FN"]
            both = [g for g in genes if pb[g] in FOUND and pr[g] in FOUND]
            pairs.append({"base": base, "run": run, "genes": len(genes),
                          "gained": len(gained), "both": len(both), "lost": len(lost),
                          "gained_genes": [gene_label(g, gene_names) for g in gained],
                          "lost_genes": [gene_label(g, gene_names) for g in lost]})
    return {"runs": ids, "pairs": pairs}


# ---------------------------------------------------------------- lost genes by expression
def gene_cpm(run: RunModels, cpm: Dict[str, float]) -> Dict[str, float]:
    """Reference gene -> summed CPM of the run's isoforms in detecting categories."""
    out: Dict[str, float] = {}
    for t in run.isoforms:
        cat, gene = run.info.get(t.tid, ("", ""))
        if cat not in _DETECTING or not gene:
            continue
        gs = [g for g in gene.split(",") if g]
        share = cpm.get(t.tid, 0.0) / len(gs)
        for g in gs:
            out[g] = out.get(g, 0.0) + share
    return out


def lost_by_decile(runs: Sequence[Optional[RunModels]], cpm: Dict[str, Dict[str, float]],
                   n_bins: int = DECILES) -> Optional[dict]:
    """For every baseline with counts.tsv: its detected genes in CPM deciles (low to high), and
    per other run the number of each decile's genes that run does not detect."""
    runs = [r for r in runs if r is not None]
    if len(runs) < 2:
        return None
    detected = {r.run_id: r.detected_genes() for r in runs}
    by_base = {}
    for base in runs:
        if base.run_id not in cpm:
            continue
        gc = gene_cpm(base, cpm[base.run_id])
        genes = sorted(detected[base.run_id], key=lambda g: (gc.get(g, 0.0), g))
        n = len(genes)
        if n == 0:
            continue
        decile = {g: min(n_bins - 1, i * n_bins // n) for i, g in enumerate(genes)}
        sizes = [0] * n_bins
        top = [None] * n_bins
        for g in genes:
            k = decile[g]
            sizes[k] += 1
            top[k] = round(gc.get(g, 0.0), 3)
        lost = {}
        for r in runs:
            if r.run_id == base.run_id:
                continue
            row = [0] * n_bins
            for g in genes:
                if g not in detected[r.run_id]:
                    row[decile[g]] += 1
            lost[r.run_id] = row
        by_base[base.run_id] = {"genes": n, "n": sizes, "cpm_max": top, "lost": lost}
    if not by_base:
        return None
    return {"deciles": n_bins, "by_base": by_base}
