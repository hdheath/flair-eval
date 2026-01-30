"""
FLAIR structural evaluation functions.

Provides utilities for evaluating isoform predictions including
genic region analysis, splice junction evaluation, and transcript classification.
"""

import subprocess
from typing import Dict, Set, Tuple

SINGLE_EXON_END_WINDOW = 100


def get_chromtoint(file):
    """Parse BED file and extract chromosome intervals"""
    totreads = 0
    chromtoint = {}
    for line in open(file):
        line = line.rstrip().split('\t')
        chrom, strand = line[0], line[5]
        iso, start, end = line[3], int(line[1]), int(line[2])
        if chrom not in chromtoint:
            chromtoint[chrom] = []
        chromtoint[chrom].append((start, end))
        totreads += 1
    return chromtoint, totreads


def get_regions(chromtoint, outfilename):
    """Merge overlapping intervals and write to file"""
    totregions = 0
    out = open(outfilename, 'w')
    for chrom in chromtoint:
        intlist = sorted(chromtoint[chrom])
        newints = []
        laststart, lastend = 0, 0
        for s, e in intlist:
            if s > lastend:
                if lastend != 0:
                    newints.append((laststart, lastend))
                laststart = s
            if e > lastend:
                lastend = e
        newints.append((laststart, lastend))
        for s, e in newints:
            out.write(f'{chrom}\t{s}\t{e}\n')
            totregions += 1
    out.close()
    return totregions


def get_intersect_count(filea, fileb):
    """Use bedtools to count intersecting regions"""
    c = 0
    result = subprocess.run(['bedtools', 'intersect', '-f', '0.5', '-u', '-a', filea, '-b', fileb],
                          capture_output=True, text=True)
    for line in result.stdout.rstrip('\n').split('\n'):
        if line:
            c += 1
    return c


def extract_sj_info(file):
    """Extract splice junction and single exon information from BED file
    Handles both BED6 and BED12 formats.
    """
    read_sjc, read_se_ends = {}, {}
    for line in open(file):
        line = line.rstrip().split('\t')
        chrom, strand = line[0], line[5]
        iso, start, end = line[3], int(line[1]), int(line[2])

        # Check if BED12 format (has exon block info in columns 10-11)
        if len(line) >= 12 and line[10] and line[11]:
            # BED12 format - extract exon blocks
            esizes, estarts = [int(x) for x in line[10].rstrip(',').split(',')], [int(x) for x in line[11].rstrip(',').split(',')]
            exons = [(start+estarts[i], start + estarts[i] + esizes[i]) for i in range(len(esizes))]
        else:
            # BED6 format - treat as single exon
            exons = [(start, end)]

        introns = tuple([(exons[x][1], exons[x+1][0]) for x in range(len(exons)-1)])
        cs = chrom
        if cs not in read_sjc:
            read_sjc[cs] = {}
            read_se_ends[cs] = {}
        if len(introns) > 0:
            if introns not in read_sjc[cs]:
                read_sjc[cs][introns] = 0
            read_sjc[cs][introns] += 1
        else:
            roundedends = (SINGLE_EXON_END_WINDOW * round(start/SINGLE_EXON_END_WINDOW),
                            SINGLE_EXON_END_WINDOW * round(end/SINGLE_EXON_END_WINDOW))
            if roundedends not in read_se_ends[cs]:
                read_se_ends[cs][roundedends] = 0
            read_se_ends[cs][roundedends] += 1
    return read_sjc, read_se_ends


def parse_gtf_transcripts(gtf_file):
    """Parse GTF file to extract transcript exon structures"""
    transcripttoexons = {}
    for line in open(gtf_file):
        if line[0] != '#':
            line = line.split('\t')
            if line[2] == 'exon':
                chrom, strand = line[0], line[6]
                start, end = int(line[3]), int(line[4])
                tname = line[-1].split('transcript_id "')[1].split('"')[0]
                if (chrom, strand) not in transcripttoexons:
                    transcripttoexons[(chrom, strand)] = {}
                if tname not in transcripttoexons[(chrom, strand)]:
                    transcripttoexons[(chrom, strand)][tname] = []
                transcripttoexons[(chrom, strand)][tname].append((start, end))
    return transcripttoexons


def build_reference_structures(transcripttoexons):
    """Build reference junction and single exon structures from GTF"""
    refjuncs, refjuncchains, refseends = {}, {}, {}
    for chrom, strand in transcripttoexons:
        refjuncs[(chrom, strand)] = set()
        refjuncchains[(chrom, strand)] = set()
        refseends[(chrom, strand)] = set()
        for tname in transcripttoexons[(chrom, strand)]:
            exons = sorted(transcripttoexons[(chrom, strand)][tname])
            if len(exons) > 1:
                introns = [(exons[x][1], exons[x+1][0]-1) for x in range(len(exons)-1)]
                refjuncchains[(chrom, strand)].add(tuple(introns))
                refjuncs[(chrom, strand)].update(set(introns))
            else:
                refseends[(chrom, strand)].add(exons[0])
    return refjuncs, refjuncchains, refseends


def classify_transcripts(isoforms_file, refjuncs, refjuncchains, refseends):
    """Classify transcripts as FSM, ISM, NIC, NNC, SEM, SEN"""
    fsm, ism, nic, nnc, sem, sen, tot = 0, 0, 0, 0, 0, 0, 0
    for line in open(isoforms_file):
        line = line.rstrip().split('\t')
        chrom, strand = line[0], line[5]
        iso, start, end = line[3], int(line[1]), int(line[2])
        esizes, estarts = [int(x) for x in line[10].rstrip(',').split(',')], [int(x) for x in line[11].rstrip(',').split(',')]
        exons = [(start+estarts[i], start + estarts[i] + esizes[i]) for i in range(len(esizes))]
        introns = tuple([(exons[x][1], exons[x+1][0]) for x in range(len(exons)-1)])
        tot += 1

        if len(introns) > 0:
            # Check if reference has any transcripts for this (chrom, strand)
            if (chrom, strand) not in refjuncchains:
                # No reference transcripts for this (chrom, strand), classify as NNC
                nnc += 1
            elif introns in refjuncchains[(chrom, strand)]:
                fsm += 1
            else:
                isISM = False
                myjuncstring = str(introns)[1:-1]
                for juncchain in refjuncchains[(chrom, strand)]:
                    if myjuncstring in str(juncchain):
                        isISM = True
                        ism += 1
                        break
                if not isISM:
                    allFound = True
                    if (chrom, strand) in refjuncs:
                        for j in introns:
                            if j not in refjuncs[(chrom, strand)]:
                                allFound = False
                                break
                    else:
                        # No reference junctions for this (chrom, strand)
                        allFound = False
                    if allFound:
                        nic += 1
                    else:
                        nnc += 1
        else:
            endsMatch = False
            if (chrom, strand) in refseends:
                for refstart, refend in refseends[(chrom, strand)]:
                    if abs(start-refstart) < SINGLE_EXON_END_WINDOW and abs(end-refend) < SINGLE_EXON_END_WINDOW:
                        endsMatch = True
                        break
            if endsMatch:
                sem += 1
            else:
                sen += 1

    return fsm, ism, nic, nnc, sem, sen, tot
