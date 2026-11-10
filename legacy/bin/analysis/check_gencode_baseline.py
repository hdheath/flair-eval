#!/usr/bin/env python3
"""Check Gencode annotation end distances as baseline."""
from collections import defaultdict
import numpy as np

GTF = "/private/groups/brookslab/reference_annotations/gencode.v38.annotation.gtf"

def parse_region(gtf_path, chrom, start, end):
    tx_exons = defaultdict(list)
    tx_meta = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if cols[0] != chrom:
                continue
            if cols[2] != 'exon':
                continue
            estart, eend = int(cols[3])-1, int(cols[4])
            if eend < start or estart > end:
                continue
            strand = cols[6]
            attrs = cols[8]
            tx_id = gene_id = gene_name = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    tx_id = attr.split('"')[1]
                elif attr.startswith('gene_id'):
                    gene_id = attr.split('"')[1].split('.')[0]
                elif attr.startswith('gene_name'):
                    gene_name = attr.split('"')[1]
            if tx_id and gene_id:
                tx_exons[tx_id].append((estart, eend))
                tx_meta[tx_id] = (chrom, strand, gene_id, gene_name or gene_id)

    gene_isoforms = defaultdict(list)
    for tx_id, exons in tx_exons.items():
        exons.sort()
        chrom, strand, gene_id, gene_name = tx_meta[tx_id]
        gene_isoforms[f"{gene_name}({gene_id})"].append(
            (tx_id, chrom, exons[0][0], exons[-1][1], strand, len(exons)))
    return gene_isoforms

def analyze_region(gene_isoforms, label):
    n_iso = sum(len(v) for v in gene_isoforms.values())
    multi = {g: isos for g, isos in gene_isoforms.items() if len(isos) > 1}
    print(f"\n{'='*80}")
    print(f"{label}")
    print(f"  Total transcripts: {n_iso}, Genes: {len(gene_isoforms)}, Multi-iso: {len(multi)}")
    
    for gene in sorted(gene_isoforms.keys()):
        isos = gene_isoforms[gene]
        print(f"  {gene}: {len(isos)} isoforms")
    
    all_tss = []
    all_tts = []
    for gene, isos in multi.items():
        for i in range(len(isos)):
            for j in range(i+1, len(isos)):
                s1, e1, st1 = isos[i][2], isos[i][3], isos[i][4]
                s2, e2, st2 = isos[j][2], isos[j][3], isos[j][4]
                if st1 == '+':
                    all_tss.append(abs(s1-s2))
                    all_tts.append(abs(e1-e2))
                else:
                    all_tss.append(abs(e1-e2))
                    all_tts.append(abs(s1-s2))
    
    if all_tss:
        tss = np.array(all_tss)
        tts = np.array(all_tts)
        both_50 = sum(1 for t,s in zip(tss,tts) if t<50 and s<50)
        both_100 = sum(1 for t,s in zip(tss,tts) if t<100 and s<100)
        print(f"\n  Pairwise TSS dist (n={len(tss)}): median={np.median(tss):.0f}, "
              f"<50bp={int(np.sum(tss<50))} ({np.sum(tss<50)/len(tss)*100:.1f}%)")
        print(f"  Pairwise TTS dist (n={len(tts)}): median={np.median(tts):.0f}, "
              f"<50bp={int(np.sum(tts<50))} ({np.sum(tts<50)/len(tts)*100:.1f}%)")
        print(f"  BOTH <50bp: {both_50}/{len(tss)} ({both_50/len(tss)*100:.1f}%)")
        print(f"  BOTH <100bp: {both_100}/{len(tss)} ({both_100/len(tss)*100:.1f}%)")

met = parse_region(GTF, "chr7", 116000000, 117000000)
analyze_region(met, "GENCODE v38 - MET locus (chr7:116-117Mb)")

chr22 = parse_region(GTF, "chr22", 0, 52000000)
analyze_region(chr22, "GENCODE v38 - full chr22")
