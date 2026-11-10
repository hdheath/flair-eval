"""
FLAIR structural evaluation functions.

Provides utilities for evaluating isoform predictions including
genic region analysis, splice junction evaluation, and transcript classification.
"""

import subprocess
from typing import Dict, Set, Tuple

SINGLE_EXON_END_WINDOW = 100

# Window (bp) within which an FSM isoform's terminal end counts as matching the
# reference transcript's end.  50 bp is the SQANTI3 default for `reference_match`.
FSM_END_WINDOW = 50

# FSM sub-categories, splitting FSM by which terminal ends are novel.  These are
# the SQANTI3 FSM subcategories under shorter names:
#     CIM   complete isoform match      == SQANTI3 "reference_match"
#     FSM5  5' end novel, 3' matches    == SQANTI3 "alternative_5end"
#     FSM3  3' end novel, 5' matches    == SQANTI3 "alternative_3end"
#     FSMTE both transcript ends novel  == SQANTI3 "alternative_3end5end"
FSM_SUBCATEGORIES = ("CIM", "FSM5", "FSM3", "FSMTE")


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


def get_intersect_count(filea, fileb, timeout=3600):
    """Use bedtools to count intersecting regions.

    Streams bedtools stdout line-by-line so we don't materialise multi-GB
    output as a single Python string (full-genome reads BED with --u can
    easily exceed memory).  Times out after `timeout` seconds with a clear
    exception so an NFS stall doesn't hang the job indefinitely.
    """
    proc = subprocess.Popen(
        ['bedtools', 'intersect', '-f', '0.5', '-u', '-a', filea, '-b', fileb],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    c = 0
    try:
        for line in proc.stdout:
            if line.strip():
                c += 1
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()
        raise RuntimeError(
            f'bedtools intersect timed out after {timeout}s on '
            f'{filea} vs {fileb}; likely NFS / cluster stall.'
        )
    if proc.returncode != 0:
        err = proc.stderr.read() if proc.stderr else ''
        raise RuntimeError(
            f'bedtools intersect failed (rc={proc.returncode}) on '
            f'{filea} vs {fileb}: {err.strip()[:500]}'
        )
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
                tname = line[8].split('transcript_id "')[1].split('"')[0]
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


def build_reference_chain_ends(transcripttoexons):
    """Map each reference intron chain -> the terminal coordinates of the transcripts carrying it.

    Kept SEPARATE from build_reference_structures() on purpose: four call sites unpack that
    function's exact 3-tuple, so widening it would break them all.  This is purely additive.

    `refjuncchains` is a bare set of intron tuples, which throws away which reference transcript
    each chain came from.  That is fine for deciding FSM, but splitting FSM by end novelty
    (CIM/FSM5/FSM3/FSMTE) needs the reference transcript's 5' and 3' coordinates, so we index
    them here by the same intron-chain key.

    COORDINATES: `transcripttoexons` holds raw GTF exons, which are 1-based inclusive, while
    query isoforms come from BED (0-based half-open).  The values returned here are converted
    to BED convention so they can be compared directly against a BED start/end:
        bed_start = gtf_start - 1        bed_end = gtf_end
    The intron keys need no conversion -- (gtf_e1, gtf_s2 - 1) is already identical to the BED
    intron (bed_exon_end, next_bed_exon_start), which is why FSM matching works at all.

    Returns:
        {(chrom, strand): {introns_tuple: set((bed_start, bed_end), ...)}}
        Single-exon reference transcripts are excluded (they have no intron chain; they are
        handled by refseends / the SEM branch).
    """
    refchainends = {}
    for chrom, strand in transcripttoexons:
        chain_ends = {}
        for tname in transcripttoexons[(chrom, strand)]:
            exons = sorted(transcripttoexons[(chrom, strand)][tname])
            if len(exons) <= 1:
                continue
            introns = tuple((exons[x][1], exons[x + 1][0] - 1)
                            for x in range(len(exons) - 1))
            chain_ends.setdefault(introns, set()).add((exons[0][0] - 1, exons[-1][1]))
        refchainends[(chrom, strand)] = chain_ends
    return refchainends


def fsm_subcategory(start, end, strand, ref_ends, window=FSM_END_WINDOW):
    """Split an FSM call into CIM / FSM5 / FSM3 / FSMTE by which terminal ends are novel.

    Args:
        start, end: query isoform coordinates, BED convention (0-based, half-open).
        strand:     '+' or '-'.  On '-' the 5' end is `end` and the 3' end is `start`.
        ref_ends:   iterable of (bed_start, bed_end) for every reference transcript sharing
                    this isoform's intron chain.
        window:     bp tolerance for calling an end "matching".

    One intron chain can be carried by several reference transcripts with different ends, so we
    score every candidate and keep the BEST -- most ends matched first, then smallest total end
    distance.  That mirrors SQANTI3, which classifies against the best-matching reference
    transcript rather than an arbitrary one.

    Returns "FSM" unchanged if `ref_ends` is empty (nothing to compare against).
    """
    best = None
    # sorted(), not raw set iteration: two reference transcripts sharing one chain can tie on both
    # matched-end count AND total end distance (one matches only the 5' end, the other only the
    # 3').  With set order the winner -- and therefore FSM5 vs FSM3 -- would depend on GTF parse
    # order, so simply regenerating the annotation could flip labels with nothing else changed.
    for ref_start, ref_end in sorted(ref_ends):
        d_start, d_end = abs(start - ref_start), abs(end - ref_end)
        if strand == '-':
            d5, d3 = d_end, d_start
        else:
            d5, d3 = d_start, d_end
        match5, match3 = d5 <= window, d3 <= window
        # rank: more matched ends wins; ties broken by the tighter total end distance
        rank = (match5 + match3, -(d5 + d3))
        if best is None or rank > best[0]:
            best = (rank, match5, match3)
    if best is None:
        return "FSM"
    _, match5, match3 = best
    if match5 and match3:
        return "CIM"
    if match3:
        return "FSM5"      # 5' novel, 3' matches
    if match5:
        return "FSM3"      # 3' novel, 5' matches
    return "FSMTE"         # both transcript ends novel


def classify_transcripts(isoforms_file, refjuncs, refjuncchains, refseends):
    """Classify transcripts as FSM, ISM, NIC, NNC, SEM, SEN

    NOTE: counts-only sibling of classify_transcripts_per_isoform(); currently has no callers.
    It deliberately does NOT compute the FSM sub-categories (CIM/FSM5/FSM3/FSMTE) -- those live
    on the per-isoform function, which is what every call site actually uses.
    """
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


def classify_transcripts_per_isoform(isoforms_file, refjuncs, refjuncchains, refseends,
                                     refchainends=None, fsm_end_window=FSM_END_WINDOW):
    """Classify each transcript and return per-isoform labels with coordinates.

    Returns list of dicts:
        [{name, chrom, start, end, strand, category, category_fine, n_exons}, ...]

    `category`      FSM, ISM, NIC, NNC, SEM, SEN  -- unchanged, always present.
    `category_fine` same as `category` for everything EXCEPT FSM, which is split by which
                    terminal ends are novel: CIM / FSM5 / FSM3 / FSMTE (see FSM_SUBCATEGORIES).

    Pass `refchainends` (from build_reference_chain_ends()) to enable the FSM split.  Without it
    `category_fine` is just a copy of `category`, so existing callers keep their exact behaviour
    and nothing that reads `category` changes.
    """
    results = []
    # The ISM test below is a substring search against str(reference_chain).  Stringifying every
    # reference chain again for every query isoform is what made this function O(minutes) per
    # sample genome-wide (profiled: 16.8s of 16.8s on chr22 alone).  Stringify each chain ONCE.
    # Purely a cache: same strings, same first-match break, byte-identical categories out.
    chain_str_cache = {}
    for line in open(isoforms_file):
        parts = line.rstrip().split('\t')
        if len(parts) < 12:
            continue
        chrom, strand = parts[0], parts[5]
        iso_name = parts[3]
        start, end = int(parts[1]), int(parts[2])
        esizes = [int(x) for x in parts[10].rstrip(',').split(',')]
        estarts = [int(x) for x in parts[11].rstrip(',').split(',')]
        exons = [(start + estarts[i], start + estarts[i] + esizes[i])
                 for i in range(len(esizes))]
        introns = tuple([(exons[x][1], exons[x + 1][0])
                         for x in range(len(exons) - 1)])
        n_exons = len(exons)

        fine = None
        if len(introns) > 0:
            if (chrom, strand) not in refjuncchains:
                cat = "NNC"
            elif introns in refjuncchains[(chrom, strand)]:
                cat = "FSM"
                if refchainends is not None:
                    fine = fsm_subcategory(
                        start, end, strand,
                        refchainends.get((chrom, strand), {}).get(introns, ()),
                        fsm_end_window)
            else:
                cat = None
                myjuncstring = str(introns)[1:-1]
                cached = chain_str_cache.get((chrom, strand))
                if cached is None:
                    cached = [str(c) for c in refjuncchains[(chrom, strand)]]
                    chain_str_cache[(chrom, strand)] = cached
                for juncchainstr in cached:
                    if myjuncstring in juncchainstr:
                        cat = "ISM"
                        break
                if cat is None:
                    allFound = True
                    if (chrom, strand) in refjuncs:
                        for j in introns:
                            if j not in refjuncs[(chrom, strand)]:
                                allFound = False
                                break
                    else:
                        allFound = False
                    cat = "NIC" if allFound else "NNC"
        else:
            cat = "SEN"
            if (chrom, strand) in refseends:
                for refstart, refend in refseends[(chrom, strand)]:
                    if (abs(start - refstart) < SINGLE_EXON_END_WINDOW and
                            abs(end - refend) < SINGLE_EXON_END_WINDOW):
                        cat = "SEM"
                        break

        results.append({
            "name": iso_name,
            "chrom": chrom,
            "start": start,
            "end": end,
            "strand": strand,
            "category": cat,
            "category_fine": fine if fine is not None else cat,
            "n_exons": n_exons,
        })

    return results
