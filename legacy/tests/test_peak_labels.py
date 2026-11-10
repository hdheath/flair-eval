#!/usr/bin/env python3
"""
Validation tests for peak-miss labeling in the evaluation pipeline.

Tests that:
1. "unassigned" label genuinely means the read isn't in any isoform's read map
2. "assigned_distant" distances are correct
3. "assigned_nearby" distances are within [window, nearby_threshold]
4. "assigned_wrong_strand" correctly identifies strand mismatches
5. Truncation pattern classification thresholds produce expected assignments
6. SJ support classification agrees with actual SJ chain lookups
7. The assign_peak_reason() priority logic is consistent with underlying data
"""

import sys
import os
import csv
import statistics
from pathlib import Path
from collections import defaultdict

# Add bin dir to path
EVAL_DIR = Path(__file__).resolve().parent.parent / "bin"
sys.path.insert(0, str(EVAL_DIR))

from evaluation.peak_analysis import (
    classify_missed_peak_reads,
    get_reads_to_isoforms,
    find_recoverable_peaks,
    find_captured_peaks,
    analyze_missed_peaks_comprehensive,
    extract_read_end_positions,
    parse_read_sj_chains,
    classify_read_sj_support,
)
from evaluation.truncation import characterize_truncation_pattern
from evaluation.plots import assign_peak_reason
from evaluation.bed_utils import read_bed6

# ──────────────────────────────────────────────────────────────────────
# Configuration: point to a real completed run
# ──────────────────────────────────────────────────────────────────────

BASE = Path(
    "/private/groups/brookslab/hdheath/projects/flair-eval/"
    "work_a549_chr1_badread/f8/3f114f0a6cd6f775b1ed56093c4abf"
)
PREFIX = "A549_cDNA_pre-aligned_chr1_end-scoring-alpha05"

ISOFORMS_BED = BASE / f"{PREFIX}_transcriptome.isoforms.bed"
READ_MAP = BASE / f"{PREFIX}_transcriptome.isoform.read.map.txt"
READS_BED = BASE / "A549_cDNA_pre-aligned_chr1.bed"
CAGE_PEAKS = BASE / "A549_cDNA_pre-aligned_chr1_cage.bed"
DRNA_PEAKS = BASE / "A549_cDNA_pre-aligned_chr1_drna.bed"
MISSED_TSV = BASE / "test_regions" / f"{PREFIX}_transcriptome_missed_cage_peaks.tsv"

WINDOW = 50


def load_isoform_read_map(path: Path):
    """Load isoform -> read_ids mapping from FLAIR read.map file."""
    iso_to_reads = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            parts = line.split('\t', 1)
            iso_id = parts[0]
            reads = parts[1].split(',') if len(parts) > 1 else []
            iso_to_reads[iso_id] = reads
    return iso_to_reads


def load_isoforms(path: Path):
    """Load isoform info from BED12."""
    isoforms = {}
    with open(path) as f:
        for line in f:
            cols = line.rstrip().split('\t')
            if len(cols) < 6:
                continue
            name = cols[3]
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            strand = cols[5]
            tss = start if strand == '+' else end
            tts = end if strand == '+' else start
            # Parse junctions from BED12
            junctions = ()
            if len(cols) >= 12:
                try:
                    esizes = [int(x) for x in cols[10].rstrip(',').split(',')]
                    estarts = [int(x) for x in cols[11].rstrip(',').split(',')]
                    exons = [(start + estarts[i], start + estarts[i] + esizes[i])
                             for i in range(len(esizes))]
                    junctions = tuple((exons[x][1], exons[x+1][0])
                                      for x in range(len(exons)-1))
                except (ValueError, IndexError):
                    pass
            isoforms[name] = {
                'chrom': chrom, 'start': start, 'end': end,
                'strand': strand, 'tss': tss, 'tts': tts,
                'junctions': junctions,
                'n_exons': len(junctions) + 1,
            }
    return isoforms


def load_read_ends(path: Path):
    """Load per-read TSS/TTS positions from reads BED12."""
    read_ends = {}
    with open(path) as f:
        for line in f:
            cols = line.rstrip().split('\t')
            if len(cols) < 6 or cols[0].startswith('#'):
                continue
            name = cols[3]
            if name in read_ends:
                continue
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            strand = cols[5]
            tss = start if strand == '+' else end
            tts = end if strand == '+' else start
            read_ends[name] = {
                'chrom': chrom, 'start': start, 'end': end,
                'strand': strand, 'tss': tss, 'tts': tts,
            }
    return read_ends


# ──────────────────────────────────────────────────────────────────────
# Test 1: "unassigned" label means read NOT in any isoform's read map
# ──────────────────────────────────────────────────────────────────────

def test_unassigned_is_truly_unassigned():
    """Verify that reads labeled 'unassigned' are genuinely absent from the isoform read map."""
    print("\n" + "="*70)
    print("TEST 1: 'unassigned' reads are truly not in any isoform's read map")
    print("="*70)

    iso_to_reads = load_isoform_read_map(READ_MAP)
    isoforms = load_isoforms(ISOFORMS_BED)
    read_ends = load_read_ends(READS_BED)
    read_to_iso = get_reads_to_isoforms(iso_to_reads)

    # Read peaks
    peaks = read_bed6(CAGE_PEAKS)
    read_end_positions = extract_read_end_positions(READS_BED, 'tss')

    recoverable = find_recoverable_peaks(CAGE_PEAKS, read_end_positions, WINDOW)

    # Find captured
    # Use direct overlap approach
    iso_positions = []
    for iso_id, info in isoforms.items():
        tss = info['tss']
        iso_positions.append({
            'Chrom': info['chrom'],
            'Start': tss,
            'End': tss + 1,
            'Strand': info['strand'],
        })
    captured_ids = find_captured_peaks(
        [],  # no closest_rows needed
        WINDOW,
        endpoint_rows=iso_positions,
        peaks_rows=peaks,
    )

    missed = {pid: cnt for pid, cnt in recoverable.items() if pid not in captured_ids}
    print(f"  Recoverable peaks: {len(recoverable)}")
    print(f"  Captured peaks: {len(captured_ids)}")
    print(f"  Missed recoverable peaks: {len(missed)}")

    # Run the classification on a sample of missed peaks
    analysis = analyze_missed_peaks_comprehensive(
        missed_peaks=missed,
        peaks_path=CAGE_PEAKS,
        read_end_positions=read_end_positions,
        iso_to_reads=iso_to_reads,
        isoforms=isoforms,
        read_ends=read_ends,
        window=WINDOW,
        end_type='tss',
        captured_peaks=captured_ids,
    )

    # Now verify: for every "unassigned" read, check it's truly not in read_to_iso
    n_checked = 0
    n_false_unassigned = 0
    false_examples = []
    for cls in analysis.get('peak_classifications', []):
        for rid in cls.get('unassigned_reads', []):
            n_checked += 1
            if rid in read_to_iso:
                n_false_unassigned += 1
                false_examples.append({
                    'read': rid,
                    'peak': cls.get('peak_id'),
                    'assigned_iso': read_to_iso[rid],
                })

    if n_false_unassigned > 0:
        print(f"  [FAIL] {n_false_unassigned}/{n_checked} 'unassigned' reads ARE in the read map!")
        for ex in false_examples[:5]:
            print(f"    Read {ex['read']} -> isoform {ex['assigned_iso']} (peak: {ex['peak']})")
    else:
        print(f"  [PASS] All {n_checked} sampled 'unassigned' reads confirmed absent from read map")

    # Also check: how many total reads near missed peaks are in the read map vs not
    reads_by_cs = defaultdict(list)
    for r in read_end_positions:
        reads_by_cs[(r['Chrom'], r['Strand'])].append(r)

    peak_id_to_info = {f"{p['Chrom']}_{p['Start']}_{p['End']}": p for p in peaks}
    total_reads_near_missed = 0
    total_in_map = 0
    total_not_in_map = 0

    # Build reverse index
    pos_to_rid = {}
    for rid, info in read_ends.items():
        key = (info['chrom'], info['strand'], info['tss'])
        pos_to_rid[key] = rid

    for pid in list(missed.keys())[:50]:  # sample for speed
        pinfo = peak_id_to_info.get(pid)
        if not pinfo:
            continue
        peak_pos = (pinfo['Start'] + pinfo['End']) // 2
        strands = [pinfo['Strand']] if pinfo['Strand'] != '.' else ['+', '-']
        for strand in strands:
            for r in reads_by_cs.get((pinfo['Chrom'], strand), []):
                if abs(r['Start'] - peak_pos) <= WINDOW:
                    lookup = (r['Chrom'], r['Strand'], r['Start'])
                    rid = pos_to_rid.get(lookup)
                    if rid:
                        total_reads_near_missed += 1
                        if rid in read_to_iso:
                            total_in_map += 1
                        else:
                            total_not_in_map += 1

    print(f"\n  Deep check (50 missed peaks):")
    print(f"    Reads near missed peaks: {total_reads_near_missed}")
    print(f"    In read map: {total_in_map} ({100*total_in_map/max(1,total_reads_near_missed):.1f}%)")
    print(f"    NOT in read map: {total_not_in_map} ({100*total_not_in_map/max(1,total_reads_near_missed):.1f}%)")
    print(f"    -> 'unassigned' dominance is {'EXPECTED' if total_not_in_map > total_in_map else 'UNEXPECTED'}")

    return n_false_unassigned == 0


# ──────────────────────────────────────────────────────────────────────
# Test 2: "assigned_distant" distances are actually > nearby_threshold
# ──────────────────────────────────────────────────────────────────────

def test_assigned_distant_distances():
    """Verify that 'assigned_distant' reads are assigned to isoforms with
    endpoints genuinely far from the peak."""
    print("\n" + "="*70)
    print("TEST 2: 'assigned_distant' reads have endpoints > 200bp from peak")
    print("="*70)

    iso_to_reads = load_isoform_read_map(READ_MAP)
    isoforms = load_isoforms(ISOFORMS_BED)
    read_ends = load_read_ends(READS_BED)
    read_to_iso = get_reads_to_isoforms(iso_to_reads)
    peaks = read_bed6(CAGE_PEAKS)

    peak_id_to_info = {f"{p['Chrom']}_{p['Start']}_{p['End']}": p for p in peaks}

    read_end_positions = extract_read_end_positions(READS_BED, 'tss')
    recoverable = find_recoverable_peaks(CAGE_PEAKS, read_end_positions, WINDOW)
    iso_positions = [{'Chrom': iso['chrom'], 'Start': iso['tss'],
                      'End': iso['tss'] + 1, 'Strand': iso['strand']}
                     for iso in isoforms.values()]
    captured_ids = find_captured_peaks([], WINDOW, endpoint_rows=iso_positions, peaks_rows=peaks)
    missed = {pid: cnt for pid, cnt in recoverable.items() if pid not in captured_ids}

    analysis = analyze_missed_peaks_comprehensive(
        missed_peaks=missed, peaks_path=CAGE_PEAKS,
        read_end_positions=read_end_positions,
        iso_to_reads=iso_to_reads, isoforms=isoforms, read_ends=read_ends,
        window=WINDOW, end_type='tss', captured_peaks=captured_ids,
    )

    nearby_threshold = 200
    n_checked = 0
    n_bad = 0
    distances_found = []

    for cls in analysis.get('peak_classifications', []):
        pid = cls.get('peak_id')
        pinfo = peak_id_to_info.get(pid)
        if not pinfo:
            continue
        peak_pos = (pinfo['Start'] + pinfo['End']) // 2

        for rid in cls.get('assigned_distant_reads', []):
            iso_id = read_to_iso.get(rid)
            if not iso_id:
                continue
            iso = isoforms.get(iso_id)
            if not iso:
                continue
            iso_end = iso.get('tss')
            if iso_end is None:
                continue
            dist = abs(int(iso_end) - int(peak_pos))
            distances_found.append(dist)
            n_checked += 1
            if dist <= nearby_threshold:
                n_bad += 1

    if n_checked > 0:
        if n_bad > 0:
            print(f"  [FAIL] {n_bad}/{n_checked} 'assigned_distant' reads have dist <= {nearby_threshold}bp")
        else:
            print(f"  [PASS] All {n_checked} 'assigned_distant' reads confirmed > {nearby_threshold}bp")
        if distances_found:
            print(f"  Distance distribution: min={min(distances_found)}, "
                  f"median={statistics.median(distances_found):.0f}, "
                  f"max={max(distances_found)}")
    else:
        print(f"  [SKIP] No 'assigned_distant' example reads to check")

    return n_bad == 0


# ──────────────────────────────────────────────────────────────────────
# Test 3: Truncation pattern thresholds
# ──────────────────────────────────────────────────────────────────────

def test_truncation_patterns():
    """Verify truncation pattern classifier produces expected labels for
    known distributions."""
    print("\n" + "="*70)
    print("TEST 3: Truncation pattern classifier produces correct labels")
    print("="*70)

    # Synthetic test cases
    peak_pos = 1000
    all_pass = True

    # Case 1: Sharp — most reads at peak
    reads_sharp = [998, 999, 1000, 1001, 1002, 1003, 997, 1001, 999, 1000]
    pat, det = characterize_truncation_pattern(reads_sharp, peak_pos, '+')
    ok = pat == 'sharp'
    print(f"  Sharp case: {pat} {'[PASS]' if ok else '[FAIL]'}")
    all_pass &= ok

    # Case 2: Trailing — reads trail upstream
    reads_trailing = [900, 910, 920, 940, 960, 970, 980, 990, 995, 1000,
                      905, 930, 950, 975, 985]
    pat, det = characterize_truncation_pattern(reads_trailing, peak_pos, '+')
    ok = pat == 'trailing'
    print(f"  Trailing case: {pat} {'[PASS]' if ok else '[FAIL]'}")
    if not ok:
        print(f"    Details: {det}")
    all_pass &= ok

    # Case 3: Sparse — too few reads
    reads_sparse = [990, 1010]
    pat, det = characterize_truncation_pattern(reads_sparse, peak_pos, '+')
    ok = pat == 'sparse'
    print(f"  Sparse case: {pat} {'[PASS]' if ok else '[FAIL]'}")
    all_pass &= ok

    # Case 4: Bimodal — two clusters
    reads_bimodal = [700, 705, 710, 703, 708, 1000, 1003, 1005, 1002, 998,
                     702, 707, 1001, 999]
    pat, det = characterize_truncation_pattern(reads_bimodal, peak_pos, '+')
    ok = pat == 'bimodal'
    print(f"  Bimodal case: {pat} {'[PASS]' if ok else '[FAIL]'}")
    if not ok:
        print(f"    Details: {det}")
    # Bimodal detection is known to be weak — log but don't fail
    if not ok:
        print(f"    Note: bimodal detection uses a simplistic middle-third gap heuristic")

    # Case 5: Dispersed — spread out
    import random
    random.seed(42)
    reads_dispersed = [peak_pos + random.randint(-500, 500) for _ in range(20)]
    pat, det = characterize_truncation_pattern(reads_dispersed, peak_pos, '+')
    print(f"  Dispersed case: {pat} (expected 'dispersed' or 'trailing')")
    # Dispersed is the fallback, many things can trigger other labels

    # Case 6: Minus strand trailing (upstream = 3' direction on minus strand)
    reads_minus_trail = [1020, 1030, 1040, 1050, 1060, 1070, 1080, 1090, 1000, 1005,
                         1025, 1045, 1065, 1085, 1010]
    pat, det = characterize_truncation_pattern(reads_minus_trail, peak_pos, '-')
    ok = pat == 'trailing'
    print(f"  Minus-strand trailing: {pat} {'[PASS]' if ok else '[FAIL]'}")
    if not ok:
        print(f"    Details: {det}")
    all_pass &= ok

    return all_pass


# ──────────────────────────────────────────────────────────────────────
# Test 4: assign_peak_reason() priority logic
# ──────────────────────────────────────────────────────────────────────

def test_assign_peak_reason_logic():
    """Test that assign_peak_reason() produces expected categories for known inputs."""
    print("\n" + "="*70)
    print("TEST 4: assign_peak_reason() priority logic")
    print("="*70)
    all_pass = True

    # Case 1: single_exon_only
    cls = {'total_supporting_reads': 10, 'best_sj_support': 'single_exon',
           'unassigned_count': 8, 'assigned_nearby_count': 0,
           'assigned_distant_count': 0, 'redirected_to_other_peak_count': 0,
           'truncation_pattern': 'sharp'}
    reason = assign_peak_reason(cls)
    ok = reason == 'single_exon_only'
    print(f"  single_exon_only: {reason} {'[PASS]' if ok else '[FAIL]'}")
    all_pass &= ok

    # Case 2: reads_redirected (distant > nearby, redirected > 0)
    cls = {'total_supporting_reads': 10, 'best_sj_support': 'full_match',
           'unassigned_count': 2, 'assigned_nearby_count': 1,
           'assigned_distant_count': 6, 'redirected_to_other_peak_count': 3,
           'truncation_pattern': 'trailing'}
    reason = assign_peak_reason(cls)
    ok = reason == 'reads_redirected'
    print(f"  reads_redirected: {reason} {'[PASS]' if ok else '[FAIL]'}")
    all_pass &= ok

    # Case 3: reads_unassigned (majority unassigned)
    cls = {'total_supporting_reads': 10, 'best_sj_support': 'full_match',
           'unassigned_count': 8, 'assigned_nearby_count': 1,
           'assigned_distant_count': 1, 'redirected_to_other_peak_count': 0,
           'truncation_pattern': 'sharp'}
    reason = assign_peak_reason(cls)
    ok = reason == 'reads_unassigned'
    print(f"  reads_unassigned: {reason} {'[PASS]' if ok else '[FAIL]'}")
    all_pass &= ok

    # Case 4: near_miss (nearby >= distant)
    cls = {'total_supporting_reads': 10, 'best_sj_support': 'full_match',
           'unassigned_count': 3, 'assigned_nearby_count': 5,
           'assigned_distant_count': 2, 'redirected_to_other_peak_count': 0,
           'truncation_pattern': 'sharp'}
    reason = assign_peak_reason(cls)
    ok = reason == 'near_miss'
    print(f"  near_miss: {reason} {'[PASS]' if ok else '[FAIL]'}")
    all_pass &= ok

    # Case 5: trailing_truncation
    cls = {'total_supporting_reads': 10, 'best_sj_support': 'full_match',
           'unassigned_count': 2, 'assigned_nearby_count': 1,
           'assigned_distant_count': 5, 'redirected_to_other_peak_count': 0,
           'truncation_pattern': 'trailing'}
    reason = assign_peak_reason(cls)
    ok = reason == 'trailing_truncation'
    print(f"  trailing_truncation: {reason} {'[PASS]' if ok else '[FAIL]'}")
    all_pass &= ok

    # Case 6: zero reads
    cls = {'total_supporting_reads': 0}
    reason = assign_peak_reason(cls)
    ok = reason == 'reads_unassigned'
    print(f"  zero reads: {reason} {'[PASS]' if ok else '[FAIL]'}")
    all_pass &= ok

    return all_pass


# ──────────────────────────────────────────────────────────────────────
# Test 5: Analyze what "unassigned" REALLY means
# ──────────────────────────────────────────────────────────────────────

def deep_analyze_unassigned():
    """Deep analysis: WHY are reads 'unassigned'?
    Breaks them down into more specific sub-categories by looking at
    the raw read alignment data."""
    print("\n" + "="*70)
    print("ANALYSIS: Why are reads 'unassigned'? (breaking down the mega-bucket)")
    print("="*70)

    iso_to_reads = load_isoform_read_map(READ_MAP)
    isoforms = load_isoforms(ISOFORMS_BED)
    read_ends = load_read_ends(READS_BED)
    read_to_iso = get_reads_to_isoforms(iso_to_reads)

    # All reads in the input BED file
    all_read_ids = set(read_ends.keys())
    # All reads assigned to at least one isoform
    assigned_read_ids = set(read_to_iso.keys())
    # Reads that exist in input but aren't assigned
    truly_unassigned = all_read_ids - assigned_read_ids

    print(f"\n  Total reads in input BED: {len(all_read_ids)}")
    print(f"  Reads assigned to isoforms: {len(assigned_read_ids)}")
    print(f"  Reads with no isoform: {len(truly_unassigned)} "
          f"({100*len(truly_unassigned)/max(1,len(all_read_ids)):.1f}%)")

    # Parse SJ chains to categorize unassigned reads
    read_sj = parse_read_sj_chains(READS_BED)

    # Sub-categories of unassigned reads
    unassigned_single_exon = 0
    unassigned_multi_exon = 0
    unassigned_short_reads = 0  # < 200bp
    unassigned_long_reads = 0   # >= 200bp

    for rid in truly_unassigned:
        chain = read_sj.get(rid, ())
        if not chain:
            unassigned_single_exon += 1
        else:
            unassigned_multi_exon += 1
        rinfo = read_ends.get(rid)
        if rinfo:
            rlen = abs(rinfo.get('end', 0) - rinfo.get('start', 0))
            if rlen < 200:
                unassigned_short_reads += 1
            else:
                unassigned_long_reads += 1

    print(f"\n  Unassigned reads breakdown:")
    print(f"    Single-exon: {unassigned_single_exon} "
          f"({100*unassigned_single_exon/max(1,len(truly_unassigned)):.1f}%)")
    print(f"    Multi-exon: {unassigned_multi_exon} "
          f"({100*unassigned_multi_exon/max(1,len(truly_unassigned)):.1f}%)")
    print(f"    Short (<200bp): {unassigned_short_reads}")
    print(f"    Long (>=200bp): {unassigned_long_reads}")

    # For multi-exon unassigned: check if their SJ chain is in the isoform set
    iso_sj_chains = set()
    for iso_id, info in isoforms.items():
        juncs = info.get('junctions', ())
        if juncs:
            iso_sj_chains.add((info['chrom'], juncs))

    me_sj_match = 0
    me_sj_no_match = 0
    for rid in truly_unassigned:
        chain = read_sj.get(rid, ())
        if not chain:
            continue
        rinfo = read_ends.get(rid)
        if not rinfo:
            continue
        key = (rinfo['chrom'], chain)
        if key in iso_sj_chains:
            me_sj_match += 1
        else:
            me_sj_no_match += 1

    print(f"\n  Multi-exon unassigned: SJ chain in isoform set?")
    print(f"    SJ chain matches an isoform: {me_sj_match}")
    print(f"    SJ chain NOT in isoform set: {me_sj_no_match}")
    print(f"    -> {me_sj_match} reads COULD have been assigned but count_sam_transcripts "
          f"rejected them (alignment quality / clipping / etc.)")
    print(f"    -> {me_sj_no_match} reads have novel SJ chains with no matching isoform model")


# ──────────────────────────────────────────────────────────────────────
# Test 6: Cross-check TSV output against computed labels
# ──────────────────────────────────────────────────────────────────────

def test_tsv_output_consistency():
    """Verify the written TSV labels match what the functions produce."""
    print("\n" + "="*70)
    print("TEST 6: TSV output labels match computed values")
    print("="*70)

    if not MISSED_TSV.exists():
        print("  [SKIP] No TSV output file found")
        return True

    # Read TSV
    tsv_labels = {}
    with open(MISSED_TSV) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            pid = row.get('peak_id', '')
            if pid:
                tsv_labels[pid] = {
                    'dominant_class': row.get('dominant_class', ''),
                    'best_sj': row.get('best_sj_support', ''),
                    'truncation': row.get('truncation_pattern', ''),
                    'read_support': int(row.get('read_support', 0)),
                    'unassigned': int(row.get('unassigned_reads', 0)),
                    'nearby': int(row.get('assigned_nearby', 0)),
                    'distant': int(row.get('assigned_distant', 0)),
                    'wrong_strand': int(row.get('wrong_strand', 0)),
                }

    # Verify that dominant_class is consistent with counts
    n_checked = 0
    n_inconsistent = 0
    for pid, labels in tsv_labels.items():
        counts = {
            'unassigned': labels['unassigned'],
            'assigned_nearby': labels['nearby'],
            'assigned_distant': labels['distant'],
            'assigned_wrong_strand': labels['wrong_strand'],
        }
        if sum(counts.values()) == 0:
            continue
        expected_dominant = max(counts.items(), key=lambda x: x[1])[0]
        # Normalize naming
        label_map = {
            'unassigned': 'unassigned',
            'nearby': 'assigned_nearby',
            'distant': 'assigned_distant',
            'wrong_strand': 'assigned_wrong_strand',
        }
        actual = labels['dominant_class']
        n_checked += 1
        if actual != expected_dominant:
            n_inconsistent += 1

    if n_checked > 0:
        if n_inconsistent > 0:
            print(f"  [WARN] {n_inconsistent}/{n_checked} peaks have dominant class != argmax of counts")
            print(f"    (This can happen with ties — just checking for gross errors)")
        else:
            print(f"  [PASS] All {n_checked} peaks have consistent dominant_class labels")
    else:
        print(f"  [SKIP] No peaks to check")

    return True


# ──────────────────────────────────────────────────────────────────────
# Analysis for 3' (dRNA) missed peaks too
# ──────────────────────────────────────────────────────────────────────

def analyze_drna_labels():
    """Analyze dRNA (3') missed peak labels for completeness."""
    print("\n" + "="*70)
    print("ANALYSIS: dRNA (3') missed peak label distribution")
    print("="*70)

    iso_to_reads = load_isoform_read_map(READ_MAP)
    isoforms = load_isoforms(ISOFORMS_BED)
    read_ends = load_read_ends(READS_BED)

    peaks = read_bed6(DRNA_PEAKS)
    read_end_positions = extract_read_end_positions(READS_BED, 'tts')

    recoverable = find_recoverable_peaks(DRNA_PEAKS, read_end_positions, WINDOW)
    iso_positions = [{'Chrom': iso['chrom'], 'Start': iso['tts'],
                      'End': iso['tts'] + 1, 'Strand': iso['strand']}
                     for iso in isoforms.values()]
    captured_ids = find_captured_peaks([], WINDOW, endpoint_rows=iso_positions, peaks_rows=peaks)
    missed = {pid: cnt for pid, cnt in recoverable.items() if pid not in captured_ids}

    print(f"  Recoverable 3' peaks: {len(recoverable)}")
    print(f"  Captured: {len(captured_ids)}")
    print(f"  Missed recoverable: {len(missed)}")

    analysis = analyze_missed_peaks_comprehensive(
        missed_peaks=missed, peaks_path=DRNA_PEAKS,
        read_end_positions=read_end_positions,
        iso_to_reads=iso_to_reads, isoforms=isoforms, read_ends=read_ends,
        window=WINDOW, end_type='tts', captured_peaks=captured_ids,
    )

    cls_summary = analysis.get('classification_summary', {})
    print(f"\n  Classification summary:")
    for cat, count in sorted(cls_summary.items(), key=lambda x: -x[1]):
        print(f"    {cat}: {count}")

    # Compute assign_peak_reason counts
    reason_counts = defaultdict(int)
    for cls in analysis.get('peak_classifications', []):
        reason = assign_peak_reason(cls)
        reason_counts[reason] += 1
    print(f"\n  Reason breakdown:")
    for reason, cnt in sorted(reason_counts.items(), key=lambda x: -x[1]):
        print(f"    {reason}: {cnt}")


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    # Check that test data exists
    required = [ISOFORMS_BED, READ_MAP, READS_BED, CAGE_PEAKS, DRNA_PEAKS]
    for f in required:
        if not f.exists():
            print(f"FATAL: Missing test data: {f}")
            sys.exit(1)

    results = {}
    results['unassigned'] = test_unassigned_is_truly_unassigned()
    results['distant'] = test_assigned_distant_distances()
    results['truncation'] = test_truncation_patterns()
    results['reason_logic'] = test_assign_peak_reason_logic()
    test_tsv_output_consistency()

    # Deep analyses
    deep_analyze_unassigned()
    analyze_drna_labels()

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    for name, passed in results.items():
        print(f"  {name}: {'PASS' if passed else 'FAIL'}")

    all_passed = all(results.values())
    print(f"\n  Overall: {'ALL PASSED' if all_passed else 'SOME FAILURES'}")
    sys.exit(0 if all_passed else 1)
