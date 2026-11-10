#!/usr/bin/env python3
"""
end_usage_fraction.py — Per-end usage-fraction landscape + APA/ATSS index
                        + 3'UTR-shortening index (Weber We-F3/F4; LAPA/DaPars PDUI).

Adopts the APA field's *fractional usage* unit instead of binary present/absent.
For each gene (multi-exon isoforms grouped by splice-junction chain, then pooled
to the gene), using per-isoform read counts as molecule weights:

  - relative usage fraction of each alternative TSS and each alternative PAS
    (read-weighted: site usage / total gene reads);
  - an ATSS index and an APA index per gene = 1 - (dominant-site usage fraction)
    i.e. how much usage is spread off the single dominant end (0 = monolithic,
    →1 = highly alternative);
  - a PDUI-style 3'UTR-shortening index = distal-PAS usage / total PAS usage,
    where distal = the PAS giving the longest 3' extension in the txn direction.

TED-confidence stratification (optional, --ted-log): each end's TED reality
(matched by position from the ted_log) is attached, and the summary plot splits
the usage-fraction distribution by high- vs low-reality ends. Demonstrates that
high-confidence TED ends carry the bulk of real alternative usage — a split no
hard-filtered competitor can reproduce.

Outputs (per method):
  <prefix>_end_usage_fraction.gene.tsv   per-gene ATSS/APA/PDUI indices
  <prefix>_end_usage_fraction.site.tsv   per-site usage fractions (+reality)
  <prefix>_end_usage_fraction.png        index distributions + usage-by-reality

EMPTY INPUT: never exits non-zero. Writes header-only TSVs and skips the plot.

Usage:
    python end_usage_fraction.py \\
        --bed ted:isoforms.bed ... \\
        [--read-map ted:flair.isoform.read.map.txt ...] \\
        [--ted-log ted:flair_ted.ted_log.tsv ...] \\
        --output-prefix run --output-dir out/ \\
        [--end-tol 50] [--min-exons 2] [--reality-bin 0.5]
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler
    from signal_utils import (parse_isoforms, gene_from_name,
                              load_read_map, lookup_read_count,
                              group_by_junction_chain)
    from alt_end_enumeration import load_end_reality, reality_at
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler
    from evaluation.signal_utils import (parse_isoforms, gene_from_name,
                                        load_read_map, lookup_read_count,
                                        group_by_junction_chain)
    from evaluation.alt_end_enumeration import load_end_reality, reality_at

apply_rc()

DEFAULT_END_TOL    = 50
DEFAULT_MIN_EXONS  = 2
DEFAULT_REALITY_BIN = 0.5

GENE_TSV_COLUMNS = [
    "method", "gene_id", "strand", "n_isoforms", "total_reads",
    "n_tss", "n_pas", "atss_index", "apa_index", "pdui",
]
SITE_TSV_COLUMNS = [
    "method", "gene_id", "axis", "site_pos", "usage_fraction",
    "reads", "ted_reality",
]


# ── Distinct-site collapse with read-weighting ─────────────────────────────────

def collapse_sites(ends: List[Tuple[int, int]], tol: int) -> List[Tuple[int, int]]:
    """Collapse (pos, reads) into distinct sites within *tol* bp.

    Each cluster's representative position is the read-weighted-most-supported
    member; reads are summed. Returns [(rep_pos, total_reads), ...].
    """
    if not ends:
        return []
    ends = sorted(ends)
    clusters: List[List[Tuple[int, int]]] = [[ends[0]]]
    for pos, reads in ends[1:]:
        if pos - clusters[-1][-1][0] > tol:
            clusters.append([(pos, reads)])
        else:
            clusters[-1].append((pos, reads))
    out = []
    for cl in clusters:
        total = sum(r for _, r in cl)
        rep = max(cl, key=lambda x: x[1])[0]
        out.append((rep, total))
    return out


def _gene_strand(isos: List[dict]) -> str:
    # majority strand (genes are single-stranded in practice)
    plus = sum(1 for i in isos if i["strand"] == "+")
    return "+" if plus >= len(isos) - plus else "-"


def _usage(sites: List[Tuple[int, int]]):
    """[(pos, reads)] → [(pos, reads, usage_fraction)] within this site set."""
    tot = sum(r for _, r in sites) or 1
    return [(pos, reads, reads / tot) for pos, reads in sites]


def _dom_off_fraction(usages) -> float:
    """1 - dominant-site usage fraction (ATSS / APA index for one site set)."""
    if not usages:
        return 0.0
    return 1.0 - max(u for _, _, u in usages)


def _pdui(pas_usage, strand: str) -> Optional[float]:
    """Distal-PAS usage fraction. Distal = furthest downstream in txn dir
    (+: max genomic pos; -: min genomic pos). None when no PAS."""
    if not pas_usage:
        return None
    if strand == "+":
        distal_pos = max(p for p, _, _ in pas_usage)
    else:
        distal_pos = min(p for p, _, _ in pas_usage)
    return next(u for p, _, u in pas_usage if p == distal_pos)


def _chain_usage(isos: List[dict], counts: Dict[str, int], tol: int):
    """Collapse one junction-chain group's ends into read-weighted usages.

    Returns (tss_usage, pas_usage, n_reads) where each usage is a list of
    (pos, reads, usage_fraction) computed WITHIN this chain.
    """
    tss_ends: List[Tuple[int, int]] = []
    pas_ends: List[Tuple[int, int]] = []
    n_reads = 0
    for iso in isos:
        n = lookup_read_count(iso["name"], counts) if counts else 1
        n = max(n, 1)
        n_reads += n
        tss_ends.append((int(iso["tss"]), n))
        pas_ends.append((int(iso["tts"]), n))
    tss_usage = _usage(collapse_sites(tss_ends, tol))
    pas_usage = _usage(collapse_sites(pas_ends, tol))
    return tss_usage, pas_usage, n_reads


def analyze_gene(isos: List[dict], counts: Dict[str, int], tol: int,
                 reality_index: Optional[dict], reality_bin: float):
    """Return (gene_row, site_rows) for one gene.

    FIX 4: alternative-end usage is a WITHIN-junction-chain concept (an
    alternative TSS/PAS for the SAME spliced structure).  The gene's
    multi-exon isoforms are first grouped by splice-junction chain
    (``group_by_junction_chain``); usage fractions and the ATSS/APA/PDUI
    dominant-site indices are computed WITHIN each chain, then POOLED to the
    gene (read-weighted mean of the per-chain indices, weighted by chain
    reads).  This matches the documented "grouped by splice-junction chain,
    then pooled to the gene" contract — previously the code pooled all chains
    of a gene together, mixing structurally distinct isoforms into one
    dominant-end fraction.
    """
    strand = _gene_strand(isos)
    chrom = isos[0]["chrom"]

    # Split the gene into junction-chain groups.  Single-exon isoforms carry no
    # junction chain and are dropped here (the caller already filters by
    # n_exons >= min_exons, but be defensive).  If grouping yields nothing
    # (e.g. malformed/degenerate input) fall back to treating the whole gene as
    # one group so the row is still emitted.
    sjc_groups = group_by_junction_chain(isos, min_exons=1)
    chain_iso_lists = list(sjc_groups.values()) if sjc_groups else [isos]

    site_rows: List[dict] = []
    n_tss = n_pas = 0
    total_reads = 0
    # Read-weighted accumulators for pooled indices.
    atss_num = apa_num = pdui_num = 0.0
    atss_den = apa_den = pdui_den = 0.0

    for chain_isos in chain_iso_lists:
        tss_usage, pas_usage, n_reads = _chain_usage(chain_isos, counts, tol)
        total_reads += n_reads
        n_tss += len(tss_usage)
        n_pas += len(pas_usage)

        atss_num += _dom_off_fraction(tss_usage) * n_reads
        atss_den += n_reads
        apa_num += _dom_off_fraction(pas_usage) * n_reads
        apa_den += n_reads
        pdui_chain = _pdui(pas_usage, strand)
        if pdui_chain is not None:
            pdui_num += pdui_chain * n_reads
            pdui_den += n_reads

        for axis, usages in (("tss", tss_usage), ("pas", pas_usage)):
            ax_real = "tss" if axis == "tss" else "tts"
            for pos, reads, frac in usages:
                rty = None
                if reality_index:
                    rty = reality_at(ax_real, chrom, pos, reality_index, tol)
                site_rows.append({
                    "axis": axis, "site_pos": pos,
                    "usage_fraction": round(frac, 5), "reads": reads,
                    "ted_reality": "" if rty is None else round(rty, 4),
                })

    atss_index = atss_num / atss_den if atss_den else 0.0
    apa_index = apa_num / apa_den if apa_den else 0.0
    pdui = pdui_num / pdui_den if pdui_den else 0.0
    total_reads = total_reads or 1

    gene_row = {
        "strand": strand, "n_isoforms": len(isos), "total_reads": total_reads,
        "n_tss": n_tss, "n_pas": n_pas,
        "atss_index": round(atss_index, 5), "apa_index": round(apa_index, 5),
        "pdui": round(pdui, 5),
    }
    return gene_row, site_rows


def analyze_method(isoforms: List[dict], counts: Dict[str, int], tol: int,
                   min_exons: int, reality_index, reality_bin):
    genes: Dict[str, List[dict]] = defaultdict(list)
    for iso in isoforms:
        if int(iso.get("n_exons", 1)) < min_exons:
            continue
        gid = iso.get("gene_id") or gene_from_name(iso.get("name", ""))
        if not gid:
            continue
        genes[gid].append(iso)

    gene_rows, site_rows = [], []
    for gid, isos in genes.items():
        g, s = analyze_gene(isos, counts, tol, reality_index, reality_bin)
        g["gene_id"] = gid
        gene_rows.append(g)
        for r in s:
            r["gene_id"] = gid
            site_rows.append(r)
    return gene_rows, site_rows


# ── TSVs ───────────────────────────────────────────────────────────────────────

def write_gene_tsv(rows_by_method, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(GENE_TSV_COLUMNS)
        for method, rows in rows_by_method.items():
            for r in sorted(rows, key=lambda x: x["gene_id"]):
                w.writerow([method] + [r[c] for c in GENE_TSV_COLUMNS[1:]])


def write_site_tsv(rows_by_method, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(SITE_TSV_COLUMNS)
        for method, rows in rows_by_method.items():
            for r in sorted(rows, key=lambda x: (x["gene_id"], x["axis"], x["site_pos"])):
                w.writerow([method] + [r[c] for c in SITE_TSV_COLUMNS[1:]])


# ── Plot ────────────────────────────────────────────────────────────────────────

def make_plot(gene_by_method, site_by_method, path, has_reality,
              reality_bin, title_prefix=""):
    methods = [m for m in gene_by_method if gene_by_method[m]]
    if not methods:
        return
    styler = ModeStyler(methods)

    fig, axes = plt.subplots(1, 3, figsize=(W2, W1 * 0.7))

    # Panel 1: ATSS / APA index distributions (ECDF-style step)
    ax = axes[0]
    for m in methods:
        atss = sorted(r["atss_index"] for r in gene_by_method[m])
        if atss:
            y = np.arange(1, len(atss) + 1) / len(atss)
            ax.step(atss, y, where="post", color=styler.color(m), lw=1.2,
                    label=f"{m} ATSS")
    style_ax(ax, xlabel="ATSS index (1 - dominant TSS usage)",
             ylabel="cumulative fraction of genes")
    ax.set_xlim(0, 1)

    # Panel 2: PDUI distribution
    ax = axes[1]
    for m in methods:
        pdui = [r["pdui"] for r in gene_by_method[m] if r["n_pas"] >= 1]
        if pdui:
            ax.hist(pdui, bins=np.linspace(0, 1, 11), histtype="step",
                    lw=1.3, color=styler.color(m), label=m)
    style_ax(ax, xlabel="PDUI (distal PAS usage fraction)", ylabel="genes")
    ax.set_xlim(0, 1)

    # Panel 3: usage fraction by TED reality bin (the TED-unique split)
    ax = axes[2]
    if has_reality:
        for m in methods:
            hi, lo = [], []
            for r in site_by_method[m]:
                rty = r["ted_reality"]
                if rty == "" or rty is None:
                    continue
                (hi if float(rty) >= reality_bin else lo).append(
                    r["usage_fraction"])
            data, labels = [], []
            if lo:
                data.append(lo); labels.append("low")
            if hi:
                data.append(hi); labels.append("high")
            if data:
                bp = ax.boxplot(data, positions=np.arange(len(data)),
                                widths=0.5, showfliers=False,
                                patch_artist=True)
                for patch in bp["boxes"]:
                    patch.set_facecolor(styler.color(m))
                    patch.set_alpha(0.6)
                ax.set_xticks(np.arange(len(labels)))
                ax.set_xticklabels(
                    [f"reality<{reality_bin:g}", f"≥{reality_bin:g}"][:len(labels)],
                    fontsize=5)
                break  # one TED method drives this panel
        style_ax(ax, ylabel="per-site usage fraction")
    else:
        # Fallback: distribution of per-site usage fractions per method
        for m in methods:
            fr = [r["usage_fraction"] for r in site_by_method[m]]
            if fr:
                ax.hist(fr, bins=np.linspace(0, 1, 11), histtype="step",
                        lw=1.3, color=styler.color(m), label=m)
        style_ax(ax, xlabel="per-site usage fraction", ylabel="sites")

    if title_prefix:
        fig.suptitle(f"{title_prefix}per-end usage-fraction landscape", fontsize=8)
    fig.tight_layout(pad=0.4)
    savefig(fig, path)


# ── CLI ────────────────────────────────────────────────────────────────────────

def _parse_pairs(entries):
    out = {}
    for e in entries or []:
        if ":" not in e:
            print(f"WARNING: skipping malformed pair '{e}'", file=sys.stderr)
            continue
        label, path = e.split(":", 1)
        out[label] = path
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed", nargs="+", required=True)
    p.add_argument("--read-map", nargs="+", default=[])
    p.add_argument("--ted-log", nargs="+", default=[])
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--output-dir", default=".")
    p.add_argument("--end-tol", type=int, default=DEFAULT_END_TOL)
    p.add_argument("--min-exons", type=int, default=DEFAULT_MIN_EXONS)
    p.add_argument("--reality-bin", type=float, default=DEFAULT_REALITY_BIN)
    p.add_argument("--title-prefix", default="")
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args(argv)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gene_tsv = out_dir / f"{args.output_prefix}_end_usage_fraction.gene.tsv"
    site_tsv = out_dir / f"{args.output_prefix}_end_usage_fraction.site.tsv"
    png = out_dir / f"{args.output_prefix}_end_usage_fraction.png"

    bed_pairs = _parse_pairs(args.bed)
    map_pairs = _parse_pairs(args.read_map)
    log_pairs = _parse_pairs(args.ted_log)

    gene_by_method, site_by_method = {}, {}
    has_reality = False
    for method, bed_path in bed_pairs.items():
        if not Path(bed_path).exists():
            print(f"WARNING: BED not found for {method}: {bed_path}", file=sys.stderr)
            gene_by_method[method] = []
            site_by_method[method] = []
            continue
        isos = parse_isoforms(bed_path)

        counts = {}
        mp = map_pairs.get(method)
        if mp and mp.lower() != "none" and Path(mp).exists():
            counts = load_read_map(mp)

        reality_index = None
        lp = log_pairs.get(method)
        if lp and lp.lower() != "none" and Path(lp).exists():
            reality_index = load_end_reality(Path(lp))
            if reality_index:
                has_reality = True

        g, s = analyze_method(isos, counts, args.end_tol, args.min_exons,
                              reality_index, args.reality_bin)
        gene_by_method[method] = g
        site_by_method[method] = s
        if args.verbose:
            print(f"  {method}: {len(isos)} isoforms → {len(g)} genes, "
                  f"{len(s)} sites", file=sys.stderr)

    write_gene_tsv(gene_by_method, gene_tsv)
    write_site_tsv(site_by_method, site_tsv)
    print(f"Wrote {gene_tsv}")
    print(f"Wrote {site_tsv}")

    if any(gene_by_method.values()):
        make_plot(gene_by_method, site_by_method, png, has_reality,
                  args.reality_bin, args.title_prefix)
        print(f"Wrote {png}")
    else:
        print("No genes with multi-exon isoforms — skipping plot "
              "(empty input handled gracefully).", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
