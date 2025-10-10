import pipettor

SINGLE_EXON_END_WINDOW = 100

BED_READS_FILE = 'WTC11.ENCFF370NFS.chr22.genomealigned.bed'
BED_READS_INTERVALS_FILE = 'WTC11.ENCFF370NFS.chr22.genomealigned.intervals.bed'  # will be created
ANNOT_FILE = 'gencode.v38.annotation.chr22.gtf'

FILE_LIST = [
    'wtc11-chr22-ogflair/ogflair-031725-best.isoforms.bed',
    'wtc11-flair3-093025-collapse-basic.isoforms.bed',
    'wtc11-flair3-093025-collapse-basic-qclip50.isoforms.bed',
    'wtc11-flair3-093025-collapse-basic-v2.isoforms.bed',
    'wtc11-flair3-093025-collapse-basic-v2-qclip50.isoforms.bed',
    'wtc11-flair3-093025-collapse-s1nosubset-v2-qclip50.isoforms.bed',
    'wtc11-flair3-093025-transcriptome-basic.isoforms.bed',
    'wtc11-flair3-093025-transcriptome-basic-qclipplus50.isoforms.bed',
    'wtc11-flair3-093025-transcriptome-basic-cstlessstringent.isoforms.bed',
    'wtc11-flair3-093025-transcriptome-basic-cstlessstringent-qclipplus50.isoforms.bed',
    'wtc11-flair3-093025-transcriptome-basic-cstlessstringent-normends.isoforms.bed',
    'wtc11-flair3-093025-transcriptome-s1nosubset-cstlessstringent-qclipplus50.isoforms.bed',
    'wtc11-flair3-093025-transcriptome-s1nosubset-cstlessstringent-qclipplus50.filt5pct.isoforms.bed',
    'wtc11-flair3-093025-transcriptome-s1nosubset-cstlessstringent.isoforms.bed',
    'wtc11-flair3-093025-transcriptome-s1nosubset-cstlessstringent.filt5pct.isoforms.bed',
    'wtc11-flair3-093025-transcriptome-newdefaults.isoforms.bed',
    ]


def get_chromtoint(file):
    totreads = 0
    chromtoint = {}
    for line in open(file): 
        line = line.rstrip().split('\t')
        chrom, strand = line[0], line[5]
        iso, start, end = line[3], int(line[1]), int(line[2])
        if chrom not in chromtoint: chromtoint[chrom] = []
        chromtoint[chrom].append((start, end))
        totreads += 1
    return chromtoint, totreads

def get_regions(chromtoint, outfilename):
    totregions= 0
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
    c = 0
    dr = pipettor.DataReader()
    pipettor.run([('bedtools', 'intersect', '-f', '0.5', '-u', '-a', filea, '-b', fileb)], stdout=dr)
    for line in dr.data.rstrip('\n').split('\n'):
        c += 1
    return c

def extract_sj_info(file):
    read_sjc, read_se_ends = {}, {}
    for line in open(file): ##bed reads file
        line = line.rstrip().split('\t')
        chrom, strand = line[0], line[5]
        iso, start, end = line[3], int(line[1]), int(line[2])
        esizes, estarts = [int(x) for x in line[10].rstrip(',').split(',')], [int(x) for x in line[11].rstrip(',').split(',')]
        exons = [(start+estarts[i], start + estarts[i] + esizes[i]) for i in range(len(esizes))]
        introns = tuple([(exons[x][1], exons[x+1][0]) for x in range(len(exons)-1)])
        cs = chrom #(chrom, strand) #can't trust read strand
        if cs not in read_sjc:
            read_sjc[cs] = {}
            read_se_ends[cs] = {}
        if len(introns) > 0:
            if introns not in read_sjc[cs]: read_sjc[cs][introns] = 0
            read_sjc[cs][introns] += 1
        else:
            roundedends = (SINGLE_EXON_END_WINDOW * round(start/SINGLE_EXON_END_WINDOW), \
                            SINGLE_EXON_END_WINDOW * round(end/SINGLE_EXON_END_WINDOW))
            if roundedends not in read_se_ends[cs]: read_se_ends[cs][roundedends] = 0
            read_se_ends[cs][roundedends] += 1
    return read_sjc, read_se_ends






# EVALUATE GENIC REGIONS FOUND RELATIVE TO THOSE REPRESENTED IN READS

chromtoint, totreads = get_chromtoint(BED_READS_FILE) ##bed reads file
totregions = get_regions(chromtoint, BED_READS_INTERVALS_FILE)

for file in FILE_LIST:
    chromtoint, _ = get_chromtoint(file)
    get_regions(chromtoint, file.split('.isoforms')[0] + '.intervals.bed')
   
    foundregions = get_intersect_count(BED_READS_INTERVALS_FILE, f'{file.split('.isoforms')[0]}.intervals.bed')
    genicreads = get_intersect_count(BED_READS_FILE, f'{file.split('.isoforms')[0]}.intervals.bed')

    print(file, totregions, foundregions, totreads, genicreads)

print()

# EVALUATE SPLICE JUNCTION CHAINS FOUND RELATIVE TO THOSE REPRESENTED IN READS
# EVALUATE SINGLE EXON READS/ISOFORMS SEPARATELY - round ends to ROUND_ENDS_WINDOW, check match

read_sjc, read_se_ends = extract_sj_info(BED_READS_FILE)

for file in FILE_LIST:
    found_sjc, found_se_ends = extract_sj_info(file)

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


    print(file, tot_sjc, sup_sjc, subset_sjc, tot_se, sup_se)


# EVALUATE TRANSCRIPTS COMPARED TO ANNOTATION
# LABEL EACH TRANSCRIPT AS FSM, ISM, NIC, NNC
# Single exon handled separately again, using SINGLE_EXON_END_WINDOW

transcripttoexons = {}
for line in open(ANNOT_FILE):  # gtf reference file
    if line[0] != '#':
        line = line.split('\t')
        if line[2] == 'exon':
            chrom, strand = line[0], line[6]
            start, end = int(line[3]), int(line[4])
            tname = line[-1].split('transcript_id "')[1].split('"')[0]
            if (chrom, strand) not in transcripttoexons: transcripttoexons[(chrom, strand)] = {}
            if tname not in transcripttoexons[(chrom, strand)]: transcripttoexons[(chrom, strand)][tname] = []
            transcripttoexons[(chrom, strand)][tname].append((start, end))

refjuncs, refjuncchains, refseends = {}, {}, {}
for chrom, strand in transcripttoexons:
    refjuncs[(chrom, strand)] = set()
    refjuncchains[(chrom, strand)] = set()
    refseends[(chrom, strand)] = set()
    for tname in transcripttoexons[(chrom, strand)]:
        exons = sorted(transcripttoexons[(chrom, strand)][tname])
        if len(exons) > 1:
            introns = [(exons[x][1], exons[x+1][0]-1) for x in range(len(exons)-1)]  # remove 1 from exon start coord to match bed
            refjuncchains[(chrom, strand)].add(tuple(introns))
            #print(tuple(introns))
            refjuncs[(chrom, strand)].update(set(introns))
        else:
            refseends[(chrom, strand)].add(exons[0])

for file in FILE_LIST:
    fsm, ism, nic, nnc, sem, sen, tot = 0, 0, 0, 0, 0, 0, 0
    for line in open(file):  # bed isoforms file
        line = line.rstrip().split('\t')
        chrom, strand = line[0], line[5]
        iso, start, end = line[3], int(line[1]), int(line[2])
        esizes, estarts = [int(x) for x in line[10].rstrip(',').split(',')], [int(x) for x in line[11].rstrip(',').split(',')]
        exons = [(start+estarts[i], start + estarts[i] + esizes[i]) for i in range(len(esizes))]
        introns = tuple([(exons[x][1], exons[x+1][0]) for x in range(len(exons)-1)])
        #print('found', introns)
        tot += 1
        if len(introns) > 0:
            if introns in refjuncchains[(chrom, strand)]:
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
                    for j in introns:
                        if j not in refjuncs[(chrom, strand)]:
                            allFound = False
                            break
                    if allFound: nic += 1
                    else: nnc += 1
        else:
            endsMatch = False
            for refstart, refend in refseends[(chrom, strand)]:
                if abs(start-refstart) < SINGLE_EXON_END_WINDOW and abs(end-refend) < SINGLE_EXON_END_WINDOW:
                    endsMatch = True
                    break
            if endsMatch: sem += 1
            else: sen += 1
    print(file, fsm, ism, nic, nnc, sem, sen)

