#!/usr/bin/env python3
"""
Analyze pairwise TSS/TTS distances between isoforms of the same gene.
Reports how many isoform pairs are near-identical (close at both ends).
"""
import os, sys, argparse
from collections import defaultdict
import numpy as np


def parse_bed_genes(bed_path):
    """Return {gene_id: [(name, chrom, start, end, strand, block_count, block_sizes, block_starts)]}"""
    gene_isoforms = defaultdict(list)
    with open(bed_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.strip().split('\t')
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            name = cols[3]
            strand = cols[5]
            block_count = int(cols[9])
            block_sizes = cols[10].rstrip(',')
            block_starts = cols[11].rstrip(',')
            # Extract gene ID
            parts = name.split('_')
            gene_id = None
            for p in parts:
                if p.startswith('ENSG'):
                    gene_id = p.split('.')[0]
                    break
            if gene_id is None:
                gene_id = parts[-1] if len(parts) > 1 else name
            gene_isoforms[gene_id].append(
                (name, chrom, start, end, strand, block_count, block_sizes, block_starts))
    return gene_isoforms


def parse_gtf_genes(gtf_path):
    """Return {gene_id: [(tx_name, chrom, start, end, strand, n_exons, '', '')]}"""
    tx_exons = defaultdict(list)  # tx_id -> [(start, end)]
    tx_meta = {}  # tx_id -> (chrom, strand, gene_id)
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.strip().split('\t')
            if cols[2] != 'exon':
                continue
            chrom = cols[0]
            start, end = int(cols[3]) - 1, int(cols[4])  # GTF 1-based -> 0-based
            strand = cols[6]
            attrs = cols[8]
            tx_id = gene_id = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    tx_id = attr.split('"')[1]
                elif attr.startswith('gene_id'):
                    gene_id = attr.split('"')[1]
            if tx_id and gene_id:
                tx_exons[tx_id].append((start, end))
                tx_meta[tx_id] = (chrom, strand, gene_id.split('.')[0])

    gene_isoforms = defaultdict(list)
    for tx_id, exons in tx_exons.items():
        exons.sort()
        chrom, strand, gene_id = tx_meta[tx_id]
        tx_start = exons[0][0]
        tx_end = exons[-1][1]
        gene_isoforms[gene_id].append(
            (tx_id, chrom, tx_start, tx_end, strand, len(exons), '', ''))
    return gene_isoforms


def analyze(gene_isoforms, label, verbose=False):
    n_iso = sum(len(v) for v in gene_isoforms.values())
    multi = {g: isos for g, isos in gene_isoforms.items() if len(isos) > 1}

    print(f"\n{'='*80}")
    print(f"{label}")
    print(f"  Total isoforms: {n_iso}, Total genes: {len(gene_isoforms)}, "
          f"Multi-isoform genes: {len(multi)}")

    if verbose:
        for gene in sorted(gene_isoforms.keys()):
            isos = gene_isoforms[gene]
            print(f"\n  Gene {gene} ({len(isos)} isoforms):")
            for rec in sorted(isos, key=lambda x: x[2]):
                name, c, s, e, st = rec[:5]
                n_ex = rec[5]
                print(f"    {st} {s:>12,}-{e:>12,}  ({e-s:>8,} bp, {n_ex} exons)  {name}")

    all_tss_dists = []
    all_tts_dists = []
    close_pairs = []

    for gene, isos in sorted(multi.items()):
        for i in range(len(isos)):
            for j in range(i+1, len(isos)):
                r1, r2 = isos[i], isos[j]
                s1, e1, st1 = r1[2], r1[3], r1[4]
                s2, e2, st2 = r2[2], r2[3], r2[4]
                if st1 == '+':
                    tss_dist = abs(s1 - s2)
                    tts_dist = abs(e1 - e2)
                else:
                    tss_dist = abs(e1 - e2)
                    tts_dist = abs(s1 - s2)
                all_tss_dists.append(tss_dist)
                all_tts_dists.append(tts_dist)
                if tss_dist < 50 and tts_dist < 50:
                    same_exon_count = r1[5] == r2[5]
                    same_blocks = (r1[5] == r2[5] and r1[6] == r2[6])
                    close_pairs.append((gene, r1[0], r2[0], tss_dist, tts_dist,
                                        s1, e1, s2, e2, st1, r1[5], r2[5],
                                        same_blocks))

    if not all_tss_dists:
        print("  No multi-isoform genes to analyze.")
        return

    tss = np.array(all_tss_dists)
    tts = np.array(all_tts_dists)

    print(f"\n  Pairwise TSS distances (n={len(tss)}):")
    print(f"    median={np.median(tss):.0f}  mean={np.mean(tss):.0f}  "
          f"min={np.min(tss)}  max={np.max(tss)}")
    for thresh in [50, 100, 200]:
        n = int(np.sum(tss < thresh))
        print(f"    <{thresh}bp: {n} ({n/len(tss)*100:.1f}%)")

    print(f"\n  Pairwise TTS distances (n={len(tts)}):")
    print(f"    median={np.median(tts):.0f}  mean={np.mean(tts):.0f}  "
          f"min={np.min(tts)}  max={np.max(tts)}")
    for thresh in [50, 100, 200]:
        n = int(np.sum(tts < thresh))
        print(f"    <{thresh}bp: {n} ({n/len(tts)*100:.1f}%)")

    both_50 = sum(1 for t, s in zip(tss, tts) if t < 50 and s < 50)
    both_100 = sum(1 for t, s in zip(tss, tts) if t < 100 and s < 100)
    print(f"\n  BOTH ends <50bp apart:  {both_50}/{len(tss)} ({both_50/len(tss)*100:.1f}%)")
    print(f"  BOTH ends <100bp apart: {both_100}/{len(tss)} ({both_100/len(tss)*100:.1f}%)")

    if close_pairs:
        same_splice = sum(1 for cp in close_pairs if cp[12])
        diff_splice = len(close_pairs) - same_splice
        print(f"\n  Near-identical end pairs (TSS<50 & TTS<50): {len(close_pairs)}")
        print(f"    Same splice structure: {same_splice}")
        print(f"    Different splice structure: {diff_splice}")
        if verbose:
            for cp in close_pairs[:25]:
                gene, n1, n2, td, ed, s1, e1, s2, e2, st, ex1, ex2, same = cp
                tag = "SAME_SPLICE" if same else "DIFF_SPLICE"
                print(f"    {gene}: TSS_d={td}bp TTS_d={ed}bp  "
                      f"exons={ex1}/{ex2}  {tag}  ({st})")
                print(f"      {n1}")
                print(f"      {n2}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed", nargs="*", help="label:path pairs for BED files")
    parser.add_argument("--gtf", nargs="*", help="label:path pairs for GTF files")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    inputs = []
    for entry in (args.bed or []):
        label, path = entry.split(":", 1)
        inputs.append((label, path, "bed"))
    for entry in (args.gtf or []):
        label, path = entry.split(":", 1)
        inputs.append((label, path, "gtf"))

    for label, path, fmt in inputs:
        if not os.path.exists(path):
            print(f"MISSING: {path}")
            continue
        if fmt == "bed":
            gene_isoforms = parse_bed_genes(path)
        else:
            gene_isoforms = parse_gtf_genes(path)
        analyze(gene_isoforms, label, verbose=args.verbose)


if __name__ == "__main__":
    main()
