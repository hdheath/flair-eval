#!/usr/bin/env python3
"""
FLAIR Evaluation Script

Performs comprehensive evaluation of FLAIR isoform predictions including:
- Genic region analysis
- Splice junction evaluation  
- Transcript classification (FSM, ISM, NIC, NNC, SEM, SEN)

Usage:
    python flair_eval.py --reads-bed <reads.bed> --isoforms-bed <isoforms.bed> --gtf <annotation.gtf> --output <output.txt>
"""

import argparse
import sys
import subprocess

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

def main():
    parser = argparse.ArgumentParser(description='Evaluate FLAIR isoform predictions')
    parser.add_argument('--reads-bed', required=True, help='Input reads BED file')
    parser.add_argument('--isoforms-bed', required=True, help='FLAIR isoforms BED file')
    parser.add_argument('--gtf', required=True, help='Reference annotation GTF file')
    parser.add_argument('--output', required=True, help='Output evaluation summary file (TSV format)')
    parser.add_argument('--verbose', action='store_true', help='Print verbose output')
    
    # Metadata arguments
    parser.add_argument('--test-name', help='Test set name')
    parser.add_argument('--dataset-name', help='Dataset name')
    parser.add_argument('--align-mode', help='Alignment mode')
    parser.add_argument('--partition-mode', help='Partition mode')
    parser.add_argument('--pipeline-mode', help='Pipeline mode (e.g., collapse_with-gtf_default)')
    parser.add_argument('--stage', help='Stage name (e.g., collapse, transcriptome)')
    
    args = parser.parse_args()
    
    if args.verbose:
        print(f"Evaluating: {args.isoforms_bed}")
        print(f"Against reads: {args.reads_bed}")
        print(f"Using annotation: {args.gtf}")
    
    # EVALUATE GENIC REGIONS
    chromtoint, totreads = get_chromtoint(args.reads_bed)
    reads_intervals_file = args.reads_bed.replace('.bed', '.intervals.bed')
    totregions = get_regions(chromtoint, reads_intervals_file)
    
    chromtoint_iso, _ = get_chromtoint(args.isoforms_bed)
    iso_intervals_file = args.isoforms_bed.replace('.bed', '.intervals.bed')
    get_regions(chromtoint_iso, iso_intervals_file)
    
    foundregions = get_intersect_count(reads_intervals_file, iso_intervals_file)
    genicreads = get_intersect_count(args.reads_bed, iso_intervals_file)
    
    # EVALUATE SPLICE JUNCTIONS
    read_sjc, read_se_ends = extract_sj_info(args.reads_bed)
    found_sjc, found_se_ends = extract_sj_info(args.isoforms_bed)
    
    found_subsets = {}
    for cs in found_sjc:
        found_subsets[cs] = set()
        for sjc in found_sjc[cs]:
            for slen in range(len(sjc)-1, 0, -1):
                for i in range(0, len(sjc)-slen+1):
                    found_subsets[cs].add(sjc[i:i+slen])
    
    tot_sjc, sup_sjc, tot_se, sup_se = 0, 0, 0, 0
    subset_sjc = 0
    for cs in read_sjc:
        if cs in found_sjc:
            for sjc in read_sjc[cs]:
                if sjc in found_sjc[cs]:
                    sup_sjc += read_sjc[cs][sjc]
                elif sjc in found_subsets[cs]:
                    subset_sjc += read_sjc[cs][sjc]
                tot_sjc += read_sjc[cs][sjc]
        else:
            for sjc in read_sjc[cs]:
                tot_sjc += read_sjc[cs][sjc]
    
    for cs in read_se_ends:
        if cs in found_se_ends:
            for se in read_se_ends[cs]:
                if se in found_se_ends[cs]:
                    sup_se += read_se_ends[cs][se]
                tot_se += read_se_ends[cs][se]
        else:
            for se in read_se_ends[cs]:
                tot_se += read_se_ends[cs][se]
    
    # EVALUATE TRANSCRIPT CLASSIFICATION
    transcripttoexons = parse_gtf_transcripts(args.gtf)
    refjuncs, refjuncchains, refseends = build_reference_structures(transcripttoexons)
    fsm, ism, nic, nnc, sem, sen, tot = classify_transcripts(args.isoforms_bed, refjuncs, refjuncchains, refseends)
    
    # Write output as simple TSV (one header row, one data row)
    with open(args.output, 'w') as outfile:
        # Build header with metadata columns first
        header = []
        if args.test_name:
            header.append('test_name')
        if args.dataset_name:
            header.append('dataset')
        if args.align_mode:
            header.append('align_mode')
        if args.partition_mode:
            header.append('partition_mode')
        
        # Extract transcriptome_mode from pipeline_mode if provided
        # pipeline_mode format is typically "transcriptome_<mode>" or "collapse_<mode>"
        transcriptome_mode = None
        if args.pipeline_mode:
            if args.pipeline_mode.startswith('transcriptome_'):
                transcriptome_mode = args.pipeline_mode.replace('transcriptome_', '', 1)
            elif args.pipeline_mode.startswith('collapse_'):
                transcriptome_mode = args.pipeline_mode.replace('collapse_', '', 1)
            else:
                transcriptome_mode = args.pipeline_mode
            header.append('transcriptome_mode')
        
        # Add metrics columns
        header.extend(['total_read_regions', 'found_regions', 'genic_reads',
                      'total_sjc', 'supported_sjc', 'subset_sjc', 'total_se', 'supported_se',
                      'FSM', 'ISM', 'NIC', 'NNC', 'SEM', 'SEN'])
        outfile.write('\t'.join(header) + '\n')
        
        # Build data row with metadata values first
        values = []
        if args.test_name:
            values.append(args.test_name)
        if args.dataset_name:
            values.append(args.dataset_name)
        if args.align_mode:
            values.append(args.align_mode)
        if args.partition_mode:
            values.append(args.partition_mode)
        if transcriptome_mode is not None:
            values.append(transcriptome_mode)
        
        # Add metrics values
        values.extend([totregions, foundregions, genicreads,
                      tot_sjc, sup_sjc, subset_sjc, tot_se, sup_se,
                      fsm, ism, nic, nnc, sem, sen])
        outfile.write('\t'.join(str(v) for v in values) + '\n')
    
    if args.verbose:
        print(f"Evaluation complete. Results written to {args.output}")
        print(f"Total transcripts classified: {tot}")

if __name__ == "__main__":
    main()