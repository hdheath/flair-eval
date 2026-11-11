"""TUSCO benchmark (CONTRACTS A30): single-isoform control genes, scored as SQANTI3 does.

Ported from SQANTI3 `src/utilities/report_qc/TUSCO_report.R` at commit
589dfd41e84b6e73ef27cf20d7a823131e02465f ("R:<line>" below cites that file).

Input: SQANTI3's `tusco_human.tsv`: '#' comment lines, then six tab-separated columns
(Ensembl gene, Ensembl transcript, gene symbol, Entrez id, RefSeq mRNA, RefSeq protein)
(R:74-82). The gene set is the distinct (ensembl, refseq, symbol) rows (R:87-89, 94), with
version suffixes stripped (R:125-130).

Rules, as in SQANTI3:
  ID type     The TUSCO column (ensembl, refseq or gene_name) that matches the most genes
              among the isoforms' associated genes; the first of those three wins a tie, and
              gene_name is the fallback when nothing matches (R:132-152). Reference gene ids
              lose their version suffix before matching (R:119, 158-165).
  Associated  Fusion isoforms count for each of their genes (R:111-122). An isoform is a
  transcript  "TUSCO transcript" when one of its genes is a TUSCO gene (R:196-200).
  TP (tx)     reference_match, or a mono-exon FSM against a mono-exon reference with both
              ends within 50 bp (R:202-210).
  TP (gene)   a gene with at least one TP transcript (R:213-214).
  PTP (tx)    FSM or ISM transcripts of genes that are not TP genes (R:216-218).
  FP (tx)     NIC, NNC, genic and fusion transcripts (R:229-231).
  FN (gene)   a TUSCO gene with no TUSCO transcript in any category (R:220-227).
  sensitivity TP genes / TUSCO genes (R:271-273).
  precision   TP transcripts / TUSCO transcripts, SQANTI3's "non-redundant precision"
              (R:254-256). `precision_redundant` = (TP + PTP transcripts) / TUSCO
              transcripts, null when there is no FSM/ISM (R:261-269).
  positive_detection_rate  TP or PTP genes / TUSCO genes (R:275-281).

Gene labels (sections.tusco.per_gene, and the TP/PTP/FP/FN counts, which are gene counts
that sum to genes_in_regions): TP as above; PTP = not TP but has an FSM/ISM transcript; FP =
has TUSCO transcripts, none of them FSM/ISM (any other category, as the R plot's FP group,
R:396-401); FN as above. SQANTI3's own TP/PTP/FP numbers are transcript counts; they are kept
in `transcripts`.

Mapping from our classifier (flair_eval_lib/classify.py) to SQANTI3 categories:
  FSM -> full-splice_match, with CIM/FSM5/FSM3/FSMTE -> reference_match/alternative_5end/
  alternative_3end/alternative_3end5end; ISM -> incomplete-splice_match; NIC -> novel_in_catalog;
  NNC -> novel_not_in_catalog; genic, genic_intron, antisense, intergenic, fusion unchanged.
  The associated gene is `ref_gene` (fusion genes are comma-separated here, '_' in SQANTI3).
Differences that can change a call:
  - Mono-exon FSMs get CIM/FSM5/FSM3/FSMTE here, where SQANTI3 says "mono-exon"; CIM on a
    mono-exon FSM is exactly SQANTI3's mono-exon TP rule (R:207-209), and a mono-exon FSM
    that is not CIM is PTP in both.
  - The FSM subcategory is scored against the reference with the most ends within 50 bp
    (then the smallest total distance); SQANTI3 takes the smallest total distance first. A
    gene can be TP here and PTP in SQANTI3 when two references share the chain.
  - Antisense and intergenic isoforms carry no gene here, so they never count for a TUSCO
    gene (SQANTI3 names them novelGene_*, which do not match TUSCO ids either).
  - genic_intron isoforms carry the gene only when it is on the same strand.

With --regions only TUSCO genes whose reference span overlaps a region count
(`genes_in_regions`); a TUSCO gene that is not in the reference GTF cannot be placed and is
out of scope. Without --regions every TUSCO gene counts (genes_in_regions == genes_total).
"""

from __future__ import annotations

import os
import re
from typing import Dict, List, Optional, Sequence

from .annotation import InputError, Regions, Transcript, _open
from .classify import Call
from .ends import md5sum

_VERSION = re.compile(r"\..*$")
ID_TYPES = ("ensembl", "refseq", "gene_name")
_FP_CATEGORIES = {"NIC", "NNC", "genic", "fusion"}


def _strip(v: str) -> str:
    return _VERSION.sub("", v)


def read_tusco(path: str) -> List[Dict[str, str]]:
    """Distinct TUSCO genes: dicts with ensembl, refseq, gene_name (versions stripped)."""
    out: List[Dict[str, str]] = []
    seen = set()
    with _open(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            p = [x.strip() for x in line.rstrip("\n").split("\t")]
            if len(p) < 5:
                raise InputError(f"{path}:{n}: TUSCO TSV needs 6 columns "
                                 "(ensembl, transcript, name, entrez, refseq, protein)")
            row = {"ensembl": _strip(p[0]), "refseq": _strip(p[4]), "gene_name": p[2]}
            key = (row["ensembl"], row["refseq"], row["gene_name"])
            if key not in seen:
                seen.add(key)
                out.append(row)
    if not out:
        raise InputError(f"{path}: no TUSCO genes")
    return out


def _gene_keys(ref_gene: str, names: Dict[str, str]) -> Dict[str, List[str]]:
    """The associated gene(s) of one isoform, as each ID type would see them."""
    genes = [g for g in ref_gene.split(",") if g]
    return {
        "ensembl": [_strip(g) for g in genes],
        "refseq": [_strip(g) for g in genes],
        "gene_name": [names.get(g, "") for g in genes if names.get(g)],
    }


def tusco_section(tusco: Sequence[Dict[str, str]], calls: Sequence[Call],
                  ref_tx: Sequence[Transcript], regions: Optional[Regions] = None,
                  path: Optional[str] = None) -> Dict:
    names = {t.gene: t.gene_name for t in ref_tx if t.gene_name}
    keyed = [_gene_keys(c.ref_gene, names) for c in calls]

    # ID type (R:132-152)
    assoc = {k: {g for kk in keyed for g in kk[k]} for k in ID_TYPES}
    matches = {k: sum(1 for row in tusco if row[k] and row[k] in assoc[k]) for k in ID_TYPES}
    id_type = max(ID_TYPES, key=lambda k: (matches[k], -ID_TYPES.index(k)))
    if matches[id_type] == 0:
        id_type = "gene_name"

    gid = {row[id_type]: row["ensembl"] or row[id_type] for row in tusco if row[id_type]}

    # scope: TUSCO genes placed in the regions through the reference annotation
    all_genes = [row["ensembl"] or row[id_type] for row in tusco]
    if regions is None:
        scope = set(all_genes)
    else:
        ref_key = {}
        for t in ref_tx:
            k = _strip(t.gene) if id_type != "gene_name" else t.gene_name
            if k in gid and regions.overlaps(t.chrom, t.start, t.end):
                ref_key[gid[k]] = True
        scope = set(ref_key)

    # TUSCO transcripts (R:196-200), each with the TUSCO genes it belongs to
    tx_genes: List[List[str]] = []
    tx_calls: List[Call] = []
    for c, kk in zip(calls, keyed):
        gs = sorted({gid[g] for g in kk[id_type] if g in gid and gid[g] in scope})
        if gs:
            tx_genes.append(gs)
            tx_calls.append(c)

    def is_tp(c: Call) -> bool:
        return c.category == "FSM" and c.subcategory == "CIM"

    tp_genes = {g for gs, c in zip(tx_genes, tx_calls) if is_tp(c) for g in gs}
    fsm_ism_genes = {g for gs, c in zip(tx_genes, tx_calls) if c.category in ("FSM", "ISM")
                     for g in gs}
    any_genes = {g for gs in tx_genes for g in gs}

    # transcript-level counts as SQANTI3 has them (one row per isoform x gene, R:111-122)
    pairs = [(g, c) for gs, c in zip(tx_genes, tx_calls) for g in gs]
    n_tx = len(pairs)
    tp_tx = sum(1 for _g, c in pairs if is_tp(c))
    ptp_tx = sum(1 for g, c in pairs if c.category in ("FSM", "ISM") and g not in tp_genes)
    fp_tx = sum(1 for _g, c in pairs if c.category in _FP_CATEGORIES)
    fsm_ism_tx = sum(1 for _g, c in pairs if c.category in ("FSM", "ISM"))

    per_gene = {}
    for g in sorted(scope):
        if g in tp_genes:
            per_gene[g] = "TP"
        elif g in fsm_ism_genes:
            per_gene[g] = "PTP"
        elif g in any_genes:
            per_gene[g] = "FP"
        else:
            per_gene[g] = "FN"
    counts = {k: sum(1 for v in per_gene.values() if v == k) for k in ("TP", "PTP", "FP", "FN")}
    n_scope = len(scope)

    def ratio(a, b):
        return round(a / b, 6) if b else None

    return {
        "file": os.path.basename(path) if path else None,
        "md5": md5sum(path) if path else None,
        "id_type": id_type,
        "genes_total": len(set(all_genes)),
        "genes_in_regions": n_scope,
        **counts,
        "sensitivity": ratio(counts["TP"], n_scope),
        "precision": ratio(tp_tx, n_tx),
        "precision_redundant": ratio(tp_tx + ptp_tx, n_tx) if fsm_ism_tx else None,
        "positive_detection_rate": ratio(counts["TP"] + counts["PTP"], n_scope),
        "transcripts": {"n": n_tx, "TP": tp_tx, "PTP": ptp_tx, "FP": fp_tx},
        "per_gene": per_gene,
        "gene_names": {(row["ensembl"] or row[id_type]): row["gene_name"] for row in tusco
                       if (row["ensembl"] or row[id_type]) in scope},
    }
