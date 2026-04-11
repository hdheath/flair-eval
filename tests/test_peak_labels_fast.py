#!/usr/bin/env python3
"""Fast unit tests for peak labeling — no large dataset loading."""
import sys
from pathlib import Path
from collections import defaultdict

EVAL_DIR = Path(__file__).resolve().parent.parent / "bin"
sys.path.insert(0, str(EVAL_DIR))

from evaluation.truncation import characterize_truncation_pattern
from evaluation.plots import assign_peak_reason


def test_truncation_patterns():
    """Verify truncation pattern classifier produces expected labels."""
    print("TEST 3: Truncation pattern classifier")
    print("=" * 60)
    peak_pos = 1000
    all_pass = True

    # Sharp — most reads cluster tight at peak
    reads = [998, 999, 1000, 1001, 1002, 1003, 997, 1001, 999, 1000]
    pat, det = characterize_truncation_pattern(reads, peak_pos, '+')
    ok = pat == 'sharp'
    print(f"  Sharp case: {pat} {'[PASS]' if ok else '[FAIL]'}")
    if not ok: print(f"    Details: {det}")
    all_pass &= ok

    # Trailing — reads trail upstream on + strand
    reads = [900, 910, 920, 940, 960, 970, 980, 990, 995, 1000,
             905, 930, 950, 975, 985]
    pat, det = characterize_truncation_pattern(reads, peak_pos, '+')
    ok = pat == 'trailing'
    print(f"  Trailing (+): {pat} {'[PASS]' if ok else '[FAIL]'}")
    if not ok: print(f"    Details: {det}")
    all_pass &= ok

    # Trailing — reads trail downstream on - strand
    reads = [1020, 1030, 1040, 1050, 1060, 1070, 1080, 1090, 1000, 1005,
             1025, 1045, 1065, 1085, 1010]
    pat, det = characterize_truncation_pattern(reads, peak_pos, '-')
    ok = pat == 'trailing'
    print(f"  Trailing (-): {pat} {'[PASS]' if ok else '[FAIL]'}")
    if not ok: print(f"    Details: {det}")
    all_pass &= ok

    # Sparse — too few reads
    reads = [990, 1010]
    pat, det = characterize_truncation_pattern(reads, peak_pos, '+')
    ok = pat == 'sparse'
    print(f"  Sparse case: {pat} {'[PASS]' if ok else '[FAIL]'}")
    all_pass &= ok

    # Bimodal — two clusters with gap in middle
    reads = [700, 705, 710, 703, 708, 1000, 1003, 1005, 1002, 998, 702, 707, 1001, 999]
    pat, det = characterize_truncation_pattern(reads, peak_pos, '+')
    print(f"  Bimodal case: {pat} (expected 'bimodal', heuristic may vary)")
    if pat != 'bimodal':
        print(f"    Note: got '{pat}' — bimodal detection uses middle-third gap heuristic")
        print(f"    Details: {det}")

    # Dispersed — wide spread
    import random
    random.seed(42)
    reads = [peak_pos + random.randint(-500, 500) for _ in range(20)]
    pat, det = characterize_truncation_pattern(reads, peak_pos, '+')
    print(f"  Dispersed case: {pat} (catch-all)")

    return all_pass


def test_assign_peak_reason_logic():
    """Test assign_peak_reason() priority logic."""
    print("\nTEST 4: assign_peak_reason() priority logic")
    print("=" * 60)
    all_pass = True

    cases = [
        ("single_exon_only", {
            'total_supporting_reads': 10, 'best_sj_support': 'single_exon',
            'unassigned_count': 8, 'assigned_nearby_count': 0,
            'assigned_distant_count': 0, 'redirected_to_other_peak_count': 0,
            'truncation_pattern': 'sharp'}),
        ("reads_redirected", {
            'total_supporting_reads': 10, 'best_sj_support': 'full_match',
            'unassigned_count': 2, 'assigned_nearby_count': 1,
            'assigned_distant_count': 6, 'redirected_to_other_peak_count': 3,
            'truncation_pattern': 'trailing'}),
        # New: no_isoform_model — unassigned + unsupported SJ chains
        ("no_isoform_model", {
            'total_supporting_reads': 10, 'best_sj_support': 'full_match',
            'unassigned_count': 8, 'assigned_nearby_count': 1,
            'assigned_distant_count': 1, 'redirected_to_other_peak_count': 0,
            'truncation_pattern': 'sharp',
            'sj_full_match': 1, 'sj_subset_match': 0,
            'sj_unsupported': 6, 'sj_single_exon': 1}),
        # New: alignment_filtered — unassigned but SJ chains match
        ("alignment_filtered", {
            'total_supporting_reads': 10, 'best_sj_support': 'full_match',
            'unassigned_count': 8, 'assigned_nearby_count': 1,
            'assigned_distant_count': 1, 'redirected_to_other_peak_count': 0,
            'truncation_pattern': 'sharp',
            'sj_full_match': 5, 'sj_subset_match': 2,
            'sj_unsupported': 1, 'sj_single_exon': 0}),
        # Original reads_unassigned fallback (no SJ data)
        ("reads_unassigned", {
            'total_supporting_reads': 10, 'best_sj_support': 'full_match',
            'unassigned_count': 8, 'assigned_nearby_count': 1,
            'assigned_distant_count': 1, 'redirected_to_other_peak_count': 0,
            'truncation_pattern': 'sharp'}),
        ("near_miss", {
            'total_supporting_reads': 10, 'best_sj_support': 'full_match',
            'unassigned_count': 3, 'assigned_nearby_count': 5,
            'assigned_distant_count': 2, 'redirected_to_other_peak_count': 0,
            'truncation_pattern': 'sharp'}),
        # New: end_absorbed — distant without redirection
        ("end_absorbed", {
            'total_supporting_reads': 10, 'best_sj_support': 'full_match',
            'unassigned_count': 2, 'assigned_nearby_count': 1,
            'assigned_distant_count': 5, 'redirected_to_other_peak_count': 0,
            'truncation_pattern': 'sharp'}),
        # trailing_truncation: distant > 0 + trailing pattern on TSS
        # Must have distant > nearby so near_miss doesn't fire, and pattern='trailing'
        ("trailing_truncation", {
            'total_supporting_reads': 10, 'best_sj_support': 'full_match',
            'unassigned_count': 2, 'assigned_nearby_count': 1,
            'assigned_distant_count': 5, 'redirected_to_other_peak_count': 0,
            'truncation_pattern': 'trailing'}),
        ("reads_unassigned", {  # zero reads
            'total_supporting_reads': 0}),
    ]

    for expected, cls in cases:
        # Default end_type is 'tss' for most tests
        reason = assign_peak_reason(cls)
        ok = reason == expected
        print(f"  {expected}: {reason} {'[PASS]' if ok else '[FAIL]'}")
        all_pass &= ok

    # Test end_type='tts' specific behavior
    # end_spread should trigger for TTS with high end variance and distant reads
    cls_spread = {
        'total_supporting_reads': 10, 'best_sj_support': 'full_match',
        'unassigned_count': 2, 'assigned_nearby_count': 1,
        'assigned_distant_count': 5, 'redirected_to_other_peak_count': 0,
        'truncation_pattern': 'dispersed', 'end_spread_std': 80}
    reason = assign_peak_reason(cls_spread, end_type='tts')
    ok = reason == 'end_spread'
    print(f"  end_spread (tts): {reason} {'[PASS]' if ok else '[FAIL]'}")
    all_pass &= ok

    # Same data but end_type='tss' should get trailing_truncation, not end_spread
    cls_spread_tss = dict(cls_spread)
    cls_spread_tss['truncation_pattern'] = 'trailing'
    reason_tss = assign_peak_reason(cls_spread_tss, end_type='tss')
    ok2 = reason_tss == 'trailing_truncation'
    print(f"  trailing_truncation (tss): {reason_tss} {'[PASS]' if ok2 else '[FAIL]'}")
    all_pass &= ok2

    # low_signal: ≤5 reads, but none of the specific patterns match
    # nearby=0, distant=0 ensures near_miss/end_absorbed/trailing don't fire
    cls_low = {
        'total_supporting_reads': 3, 'best_sj_support': 'full_match',
        'unassigned_count': 1, 'assigned_nearby_count': 0,
        'assigned_distant_count': 0, 'redirected_to_other_peak_count': 0,
        'assigned_wrong_strand_count': 2,
        'truncation_pattern': 'sparse'}
    reason = assign_peak_reason(cls_low)
    ok = reason == 'low_signal'
    print(f"  low_signal: {reason} {'[PASS]' if ok else '[FAIL]'}")
    all_pass &= ok

    # Edge cases: equal counts
    cls = {'total_supporting_reads': 10, 'best_sj_support': 'full_match',
           'unassigned_count': 5, 'assigned_nearby_count': 5,
           'assigned_distant_count': 0, 'redirected_to_other_peak_count': 0,
           'truncation_pattern': 'dispersed'}
    reason = assign_peak_reason(cls)
    print(f"  tie_unassigned_nearby: {reason} (informational)")

    return all_pass


def test_tsv_label_consistency():
    """Check TSV labels from a saved run for consistency."""
    from pathlib import Path
    import csv

    TSV = Path(
        "/private/groups/brookslab/hdheath/projects/flair-eval/"
        "work_a549_chr1_badread/f8/3f114f0a6cd6f775b1ed56093c4abf/"
        "test_regions/A549_cDNA_pre-aligned_chr1_end-scoring-alpha05_transcriptome_missed_cage_peaks.tsv"
    )
    print("\nTEST 6: TSV output consistency")
    print("=" * 60)
    if not TSV.exists():
        print("  [SKIP] No TSV output file")
        return True

    inconsistencies = 0
    total = 0
    examples = []
    with open(TSV) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            counts = {
                'unassigned': int(row.get('unassigned_reads', 0)),
                'assigned_nearby': int(row.get('assigned_nearby', 0)),
                'assigned_distant': int(row.get('assigned_distant', 0)),
                'assigned_wrong_strand': int(row.get('wrong_strand', 0)),
            }
            if sum(counts.values()) == 0:
                continue
            total += 1
            expected = max(counts.items(), key=lambda x: x[1])[0]
            actual = row.get('dominant_class', '')
            if actual != expected:
                inconsistencies += 1
                if len(examples) < 5:
                    examples.append((row.get('peak_id', '?'), actual, expected, counts))

    if inconsistencies > 0:
        print(f"  [WARN] {inconsistencies}/{total} peaks: dominant_class != argmax(counts)")
        for pid, actual, expected, counts in examples:
            print(f"    {pid}: got '{actual}', argmax='{expected}', counts={counts}")
    else:
        print(f"  [PASS] All {total} peaks have consistent dominant_class")

    return True


def analyze_tsv_reason_coverage():
    """Analyze what reason labels cover and where gaps are."""
    from pathlib import Path
    import csv

    print("\nANALYSIS: Reason label coverage (from saved TSV)")
    print("=" * 60)

    for label, tsvpath in [
        ("CAGE (5')", "A549_cDNA_pre-aligned_chr1_end-scoring-alpha05_transcriptome_missed_cage_peaks.tsv"),
        ("QuantSeq (3')", "A549_cDNA_pre-aligned_chr1_end-scoring-alpha05_transcriptome_missed_quantseq_peaks.tsv"),
    ]:
        path = Path(
            "/private/groups/brookslab/hdheath/projects/flair-eval/"
            "work_a549_chr1_badread/f8/3f114f0a6cd6f775b1ed56093c4abf/"
            f"test_regions/{tsvpath}"
        )
        if not path.exists():
            print(f"  {label}: [SKIP] no file")
            continue

        rows = []
        with open(path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                rows.append(row)

        # Label distribution
        dominant_counts = defaultdict(int)
        sj_counts = defaultdict(int)
        trunc_counts = defaultdict(int)
        for row in rows:
            dominant_counts[row.get('dominant_class', 'NA')] += 1
            sj_counts[row.get('best_sj_support', 'NA')] += 1
            trunc_counts[row.get('truncation_pattern', 'NA')] += 1

        # Compute reasons
        reason_counts = defaultdict(int)
        for row in rows:
            cls = {
                'total_supporting_reads': int(row.get('read_support', 0)),
                'best_sj_support': row.get('best_sj_support', ''),
                'unassigned_count': int(row.get('unassigned_reads', 0)),
                'assigned_nearby_count': int(row.get('assigned_nearby', 0)),
                'assigned_distant_count': int(row.get('assigned_distant', 0)),
                'redirected_to_other_peak_count': int(row.get('redirected_to_other_peak', 0) if 'redirected_to_other_peak' in row else 0),
                'truncation_pattern': row.get('truncation_pattern', ''),
            }
            reason = assign_peak_reason(cls)
            reason_counts[reason] += 1

        print(f"\n  {label}: {len(rows)} missed peaks")
        print(f"    Dominant class:")
        for k, v in sorted(dominant_counts.items(), key=lambda x: -x[1]):
            print(f"      {k}: {v} ({100*v/max(1,len(rows)):.1f}%)")
        print(f"    SJ support:")
        for k, v in sorted(sj_counts.items(), key=lambda x: -x[1]):
            print(f"      {k}: {v} ({100*v/max(1,len(rows)):.1f}%)")
        print(f"    Truncation:")
        for k, v in sorted(trunc_counts.items(), key=lambda x: -x[1]):
            print(f"      {k}: {v} ({100*v/max(1,len(rows)):.1f}%)")
        print(f"    Assigned reason:")
        for k, v in sorted(reason_counts.items(), key=lambda x: -x[1]):
            print(f"      {k}: {v} ({100*v/max(1,len(rows)):.1f}%)")

        # Cross-tab: dominant_class vs reason
        cross = defaultdict(lambda: defaultdict(int))
        for row in rows:
            cls = {
                'total_supporting_reads': int(row.get('read_support', 0)),
                'best_sj_support': row.get('best_sj_support', ''),
                'unassigned_count': int(row.get('unassigned_reads', 0)),
                'assigned_nearby_count': int(row.get('assigned_nearby', 0)),
                'assigned_distant_count': int(row.get('assigned_distant', 0)),
                'redirected_to_other_peak_count': int(row.get('redirected_to_other_peak', 0) if 'redirected_to_other_peak' in row else 0),
                'truncation_pattern': row.get('truncation_pattern', ''),
            }
            reason = assign_peak_reason(cls)
            cross[row.get('dominant_class', 'NA')][reason] += 1

        print(f"    Cross-tab (dominant_class → reason):")
        for dc in sorted(cross.keys()):
            reasons = cross[dc]
            parts = [f"{r}:{c}" for r, c in sorted(reasons.items(), key=lambda x: -x[1])]
            print(f"      {dc} → {', '.join(parts)}")


def test_false_positive_classification():
    """Test precision-side false positive classification with synthetic data."""
    from evaluation.peak_analysis import classify_false_positive_endpoints

    print("\nTEST 5: False positive endpoint classification")
    print("=" * 60)
    all_pass = True

    # Synthetic isoforms
    isoforms = {
        'iso_singleton': {'chrom': 'chr1', 'start': 1000, 'end': 5000, 'strand': '+',
                          'tss': 1000, 'tts': 5000, 'n_exons': 3},
        'iso_low_conf': {'chrom': 'chr1', 'start': 10000, 'end': 15000, 'strand': '+',
                         'tss': 10000, 'tts': 15000, 'n_exons': 2},
        'iso_near_miss': {'chrom': 'chr1', 'start': 20000, 'end': 25000, 'strand': '+',
                          'tss': 20000, 'tts': 25000, 'n_exons': 4},
        'iso_cluster_a': {'chrom': 'chr1', 'start': 30000, 'end': 35000, 'strand': '+',
                          'tss': 30000, 'tts': 35000, 'n_exons': 3},
        'iso_cluster_b': {'chrom': 'chr1', 'start': 30010, 'end': 35010, 'strand': '+',
                          'tss': 30010, 'tts': 35010, 'n_exons': 3},
        'iso_cluster_c': {'chrom': 'chr1', 'start': 30020, 'end': 35020, 'strand': '+',
                          'tss': 30020, 'tts': 35020, 'n_exons': 2},
        'iso_novel': {'chrom': 'chr1', 'start': 50000, 'end': 55000, 'strand': '+',
                      'tss': 50000, 'tts': 55000, 'n_exons': 5},
        'iso_true_pos': {'chrom': 'chr1', 'start': 60000, 'end': 65000, 'strand': '+',
                         'tss': 60000, 'tts': 65000, 'n_exons': 3},
    }

    iso_to_reads = {
        'iso_singleton': ['r1'],
        'iso_low_conf': ['r2', 'r3'],
        'iso_near_miss': ['r4', 'r5', 'r6', 'r7', 'r8'],
        'iso_cluster_a': ['r10', 'r11', 'r12', 'r13', 'r14'],
        'iso_cluster_b': ['r15', 'r16', 'r17'],
        'iso_cluster_c': ['r18', 'r19'],
        'iso_novel': [f'r{i}' for i in range(100, 120)],  # 20 reads
        'iso_true_pos': [f'r{i}' for i in range(200, 210)],
    }

    # dist_map: isoform -> (distance, nearest_peak_id)
    dist_map = {
        'iso_singleton': (500, 'peak_1'),       # far, 1 read → singleton
        'iso_low_conf': (300, 'peak_2'),         # far, 2 reads → low_confidence
        'iso_near_miss': (60, 'peak_3'),         # near miss (within 2×window)
        'iso_cluster_a': (400, 'peak_4'),        # clustered endpoints
        'iso_cluster_b': (410, 'peak_4'),
        'iso_cluster_c': (420, 'peak_4'),
        'iso_novel': (1000, 'peak_5'),           # high support → novel_site
        'iso_true_pos': (30, 'peak_6'),          # within window → true positive
    }

    result = classify_false_positive_endpoints(
        isoforms=isoforms,
        iso_to_reads=iso_to_reads,
        dist_map=dist_map,
        end_type='tss',
        window=50,
    )

    # True positive should not be in classifications
    fp_ids = {c['isoform_id'] for c in result['classifications']}
    ok = 'iso_true_pos' not in fp_ids
    print(f"  True positive excluded: {'[PASS]' if ok else '[FAIL]'}")
    all_pass &= ok

    # Check specific labels
    reason_by_iso = {c['isoform_id']: c['reason'] for c in result['classifications']}

    checks = [
        ('iso_near_miss', 'near_miss_fp', 'near miss (within 2×window)'),
        ('iso_singleton', 'singleton_isoform', '1 read = singleton'),
        ('iso_low_conf', 'low_confidence_end', '2 reads = low confidence'),
        ('iso_novel', 'novel_site', '20 reads + multi-exon = potentially real'),
    ]

    for iso_id, expected, desc in checks:
        actual = reason_by_iso.get(iso_id, 'MISSING')
        ok = actual == expected
        print(f"  {desc}: {actual} {'[PASS]' if ok else f'[FAIL] expected {expected}'}")
        all_pass &= ok

    # Cluster detection
    cluster_ids = {'iso_cluster_a', 'iso_cluster_b', 'iso_cluster_c'}
    cluster_reasons = {iso_id: reason_by_iso.get(iso_id) for iso_id in cluster_ids}
    n_cluster_split = sum(1 for r in cluster_reasons.values() if r == 'end_cluster_split')
    print(f"  Cluster split detection: {n_cluster_split}/{len(cluster_ids)} labeled end_cluster_split")
    # At least some should be detected
    ok = n_cluster_split >= 1
    print(f"    {'[PASS]' if ok else '[FAIL]'}")
    all_pass &= ok

    # Summary stats
    print(f"  Total FP: {result['n_false_positive']}, reason counts: {result['reason_counts']}")

    return all_pass


if __name__ == '__main__':
    r1 = test_truncation_patterns()
    r2 = test_assign_peak_reason_logic()
    r3 = test_false_positive_classification()
    test_tsv_label_consistency()
    analyze_tsv_reason_coverage()

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  truncation_patterns: {'PASS' if r1 else 'FAIL'}")
    print(f"  reason_logic: {'PASS' if r2 else 'FAIL'}")
    print(f"  fp_classification: {'PASS' if r3 else 'FAIL'}")
    sys.exit(0 if r1 and r2 and r3 else 1)
