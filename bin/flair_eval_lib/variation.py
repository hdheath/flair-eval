"""How a gene's isoforms differ: by ends, by splicing, or both (CONTRACTS A32).

Genes are the isoforms' own gene ids (from the FLAIR GTF), on main chromosomes only. For each
gene with >= 2 isoforms:

  chain      an isoform's intron chain together with its strand; every mono-exon isoform
             on a strand has the same (empty) chain.
  ends vary  within a group of isoforms, the 5' ends (or the 3' ends) span more than
             END_DIFF bp: max(position) - min(position) > 50.

  ends_only    one chain only, and its ends vary.
  splice_only  more than one chain, and within every chain the ends do not vary
               (a chain with one isoform never varies).
  both         more than one chain, and the ends vary within at least one chain.
  neither      one chain only and its ends do not vary (redundant models). A32 names three
               classes; this fourth one keeps genes_multi = the sum of all four.
"""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

from .annotation import Transcript
from .ends import end_positions, is_main_chrom

END_DIFF = 50


def _ends_vary(models: Sequence[Transcript]) -> bool:
    if len(models) < 2:
        return False
    fives, threes = zip(*(end_positions(t) for t in models))
    return max(fives) - min(fives) > END_DIFF or max(threes) - min(threes) > END_DIFF


def gene_class(models: Sequence[Transcript]) -> str:
    chains: Dict[Tuple, List[Transcript]] = {}
    for t in models:
        chains.setdefault((t.strand, t.introns), []).append(t)
    vary = any(_ends_vary(g) for g in chains.values())
    if len(chains) == 1:
        return "ends_only" if vary else "neither"
    return "both" if vary else "splice_only"


def variation_section(isoforms: Sequence[Transcript]) -> Dict:
    by_gene: Dict[str, List[Transcript]] = {}
    for t in isoforms:
        if is_main_chrom(t.chrom):
            by_gene.setdefault(t.gene, []).append(t)
    out = {"genes_multi": 0, "ends_only": 0, "splice_only": 0, "both": 0, "neither": 0}
    for models in by_gene.values():
        if len(models) < 2:
            continue
        out["genes_multi"] += 1
        out[gene_class(models)] += 1
    return out
