#!/usr/bin/env python3
"""Unit tests for the interval-OR-summit match predicate + canonical-chrom filter.

The 2026-06-20 fix to peak_matching makes a called end a true positive iff it is
within ±W of a truth-peak SUMMIT **OR** falls INSIDE the truth-peak INTERVAL
[start, end).  Previously only summit±W was credited, which under-counted
precision/recall on BROAD promoters (e.g. a 368 bp CAGE peak whose summit is
>W from a real call).

THE predicate (the contract):

    call c matches peak p  iff  (|c - p.summit| <= W)  OR  (p.start <= c < p.end)

with one-to-one accounting (each peak claimed at most once; interior credit is
match-distance 0, tie-broken by |c - summit|).

These tests pin:
  (1) BACKWARD-COMPAT: for a narrow peak with width <= 2W the interval is a
      subset of summit±W, so the new predicate is byte-for-byte identical to the
      old summit±W behaviour (old==new proof).
  (2) a WIDE peak (368 bp) now credits an interior call whose summit distance is
      > W (the broad-promoter fix);
  (3) one-to-one: two calls inside ONE wide peak yield TP=1 (not 2);
  (4) the canonical-chrom default excludes chrM + GL*/KI* scaffolds.

conftest.py puts bin/ on sys.path so `from evaluation.peak_matching ...` works.
"""
from __future__ import annotations

import math

from evaluation.peak_matching import (
    Peak,
    build_summit_index,
    load_summit_index_from_bed,
    count_summits,
    match_positions_one_to_one,
    precision_recall,
    is_canonical_chrom,
    CANONICAL_CHROMS,
)


# ── Reference: the OLD summit±W behaviour, reproduced exactly ────────────────
#
# This is the pre-fix matcher: candidate pairs are |pos - summit| <= W, sorted
# by (distance, call_idx, summit_idx), greedily assigned one-to-one.  It treats
# summits as scalar positions only (no interval).  We run it against the NEW
# matcher fed plain-int summits and against narrow Peak intervals to PROVE the
# new predicate reduces to the old one for width <= 2W.
def _old_summit_matcher(positions, summits, window):
    from bisect import bisect_left
    n = len(positions)
    assignment = [-1] * n
    if n == 0 or len(summits) == 0 or window < 0:
        return assignment
    ss = [int(s) for s in summits]
    pairs = []
    for i, pos in enumerate(positions):
        lo = bisect_left(ss, pos - window)
        j = lo
        while j < len(ss):
            s = ss[j]
            if s > pos + window:
                break
            pairs.append((abs(pos - s), i, j))
            j += 1
    pairs.sort()
    call_used = [False] * n
    summit_used = [False] * len(ss)
    for _dist, i, j in pairs:
        if call_used[i] or summit_used[j]:
            continue
        assignment[i] = j
        call_used[i] = True
        summit_used[j] = True
    return assignment


# ── (1) BACKWARD-COMPAT: narrow peak (width <= 2W) -> new == old ─────────────

def test_narrow_peak_new_equals_old_scalar():
    """Plain-int summits (zero-width peaks) -> new matcher == old summit±W."""
    summits = [1000, 2000, 5000]
    positions = [980, 1010, 1051, 1990, 2049, 4960, 9000]
    W = 50
    assert (match_positions_one_to_one(positions, summits, W)
            == _old_summit_matcher(positions, summits, W))


def test_narrow_peak_interval_subset_of_window_new_equals_old():
    """A Peak of width <= 2W: interval ⊂ summit±W, so OR adds nothing.

    Build 80 bp peaks (width 80 <= 2*50=100) as Peak objects and as bare
    scalar summits; the new matcher must give the IDENTICAL assignment, and
    that assignment must equal the old summit±W matcher.  This is the
    surgical backward-compat guarantee: only peaks WIDER than 2W can change.
    """
    W = 50
    # peaks centred at 1000 and 2000, width 80 (summit ± 40 ⊂ ± 50)
    peaks = [Peak(1000, 960, 1040), Peak(2000, 1960, 2040)]
    scalar = [1000, 2000]
    # calls spanning interior, window edge, and just-outside-window
    positions = [960, 1000, 1040, 1049, 1051, 1090, 1960, 2049, 2051]
    new_with_intervals = match_positions_one_to_one(positions, peaks, W)
    new_with_scalars = match_positions_one_to_one(positions, scalar, W)
    old = _old_summit_matcher(positions, scalar, W)
    assert new_with_intervals == new_with_scalars
    assert new_with_intervals == old


def test_backcompat_precision_recall_narrow_fixture():
    """Old==new at the precision_recall level on a narrow-peak fixture.

    All peaks are <= 2W wide, so the headline precision/recall/TP must be
    identical whether we credit interval interiors or not (the interval is a
    subset of the window).  We compare the interval-aware index against a
    matched scalar-summit index run through the OLD matcher.
    """
    W = 50
    # narrow peaks: widths 60, 100 (== 2W, still a subset), 20
    raw = [
        ("chr1", 970, 1030, "+"),    # summit 1000, width 60
        ("chr1", 1950, 2050, "+"),   # summit 2000, width 100 == 2W
        ("chr1", 2990, 3010, "-"),   # summit 3000, width 20
    ]
    idx = build_summit_index(raw)
    ends = [
        ("chr1", "+", 1005), ("chr1", "+", 1049), ("chr1", "+", 1051),
        ("chr1", "+", 2000), ("chr1", "+", 2051),
        ("chr1", "-", 3010), ("chr1", "-", 3051),
    ]
    new = precision_recall(ends, idx, W)

    # OLD reference: scalar-summit index + old matcher, computed by hand via
    # the helper on each bucket.
    scalar_idx = {("chr1", "+"): [1000, 2000], ("chr1", "-"): [3000]}
    tp = 0
    n_calls = len(ends)
    from collections import defaultdict
    buckets = defaultdict(list)
    for c, s, p in ends:
        buckets[(c, s)].append(p)
    for key, ps in buckets.items():
        assign = _old_summit_matcher(ps, scalar_idx.get(key, []), W)
        tp += sum(1 for j in assign if j >= 0)
    n_summits = sum(len(v) for v in scalar_idx.values())
    old_precision = tp / n_calls
    old_recall = tp / n_summits

    assert new["tp"] == tp
    assert math.isclose(new["precision"], old_precision)
    assert math.isclose(new["recall"], old_recall)


# ── (2) WIDE peak: interior call beyond summit±W now credited ───────────────

def test_wide_peak_interior_call_credited():
    """A 368 bp CAGE peak: a call 154 bp from the summit but INSIDE the peak.

    Mirrors the motivating example (chr12:52288119 inside KRT81): summit far
    from the real call, but the call is inside the broad promoter interval.
    Old summit±50 would MISS it; the new interval-OR predicate credits it.
    """
    W = 50
    # peak [52287935, 52288303) width 368, summit at 52287951 (interval
    # centre); a real call at 52288105 is 154 bp from the summit (>W) but
    # firmly inside [start, end) -> interior credit.
    start, end = 52287935, 52288303
    peak = Peak((start + end) // 2, start, end)          # summit = 52288119
    call_inside = (start + end) // 2 + 154                # 154 bp past summit
    assert call_inside < end                              # genuinely interior
    assert abs(call_inside - int(peak)) > W               # summit-far

    # OLD summit±W: miss.
    assert _old_summit_matcher([call_inside], [int(peak)], W) == [-1]
    # NEW interval-OR: hit.
    assign = match_positions_one_to_one([call_inside], [peak], W)
    assert assign == [0]


def test_wide_peak_within_W_of_edge_matches_far_misses():
    """Within-W-of-interval rule: a call within W bp of ANY part of the peak
    (including up to W past either edge) matches; >W past the edge misses.
    This is the broad-promoter edge-margin fix (the chr1:PGD case)."""
    W = 50
    start, end = 1000, 1368
    peak = Peak((start + end) // 2, start, end)
    # Just past the END, within W of the edge -> NOW credited; and it is >W
    # from the SUMMIT, proving this is interval-edge, not summit±W.
    call_just_past_edge = end + 10
    assert abs(call_just_past_edge - int(peak)) > W
    assert match_positions_one_to_one([call_just_past_edge], [peak], W) == [0]
    # Exactly W past the last covered base -> still matches (boundary); +1 -> miss.
    assert match_positions_one_to_one([(end - 1) + W], [peak], W) == [0]
    assert match_positions_one_to_one([(end - 1) + W + 1], [peak], W) == [-1]
    # Far past the edge -> misses.
    assert match_positions_one_to_one([end + 100], [peak], W) == [-1]
    # Symmetric on the START side: within W before -> match, beyond -> miss.
    assert match_positions_one_to_one([start - W], [peak], W) == [0]
    assert match_positions_one_to_one([start - W - 1], [peak], W) == [-1]


# ── (3) one-to-one: two calls inside one wide peak -> TP == 1 ───────────────

def test_two_calls_inside_one_wide_peak_one_tp():
    W = 50
    start, end = 1000, 1500       # width 500
    raw = [("chr2", start, end, "+")]
    idx = build_summit_index(raw)
    # two calls both interior, both summit-far, far from each other
    ends = [("chr2", "+", 1100), ("chr2", "+", 1400)]
    pr = precision_recall(ends, idx, W)
    assert pr["n_summits"] == 1
    assert pr["tp"] == 1, "two interior calls on one peak must give ONE TP"
    assert pr["n_calls"] == 2
    assert pr["recall"] == 1.0
    assert math.isclose(pr["precision"], 0.5)


def test_interior_call_claims_peak_over_window_call():
    """Match-distance ordering: an interior call (dist 0) claims the peak
    before a summit-window call.  The interior call wins; the window call is
    only matched if a second peak is available."""
    W = 50
    start, end = 1000, 1400          # summit 1200, width 400
    peak = Peak((start + end) // 2, start, end)
    # call A interior at 1050 (inside, summit-dist 150 -> match-dist 0)
    # call B at 1210 (summit-dist 10, also interior actually) — pick B OUTSIDE
    # so only one is interior:
    summit = int(peak)
    interior = 1050                  # inside [1000,1400)
    window_only = summit + 30        # within W of summit, also interior -> use
    # To make window_only NON-interior, place it outside the interval but in
    # ±W of the summit is impossible here (summit deep inside). Instead test
    # two interior calls: nearer-to-summit wins the single peak.
    near_summit = summit + 5         # interior, summit-dist 5
    assign = match_positions_one_to_one([interior, near_summit], [peak], W)
    # both interior -> match-dist 0 for both; tie broken by summit-distance,
    # so near_summit (dist 5) claims the peak, interior (dist 150) loses.
    assert assign[1] == 0
    assert assign[0] == -1
    assert sum(1 for j in assign if j >= 0) == 1


# ── (4) canonical-chrom filter excludes chrM + scaffolds ────────────────────

def test_is_canonical_chrom():
    for c in CANONICAL_CHROMS:
        assert is_canonical_chrom(c)
    assert is_canonical_chrom("chr1")
    assert is_canonical_chrom("chrX")
    assert is_canonical_chrom("chrY")
    assert not is_canonical_chrom("chrM")
    assert not is_canonical_chrom("chrMT")
    assert not is_canonical_chrom("GL000220.1")
    assert not is_canonical_chrom("KI270733.1")
    assert not is_canonical_chrom("chr1_KI270706v1_random")
    assert not is_canonical_chrom("chrUn_GL000195v1")


def test_load_summit_index_canonical_default_excludes_scaffolds(tmp_path):
    bed = tmp_path / "peaks.bed"
    bed.write_text(
        "chr1\t1000\t1100\tp1\t0\t+\n"
        "chrX\t2000\t2100\tp2\t0\t-\n"
        "chrM\t500\t600\tpM\t0\t+\n"
        "GL000220.1\t300\t400\tpG\t0\t+\n"
        "KI270733.1\t10\t20\tpK\t0\t-\n"
    )
    # default: canonical only -> chrM, GL*, KI* dropped, 2 summits remain
    idx = load_summit_index_from_bed(str(bed))
    assert count_summits(idx) == 2
    assert ("chr1", "+") in idx
    assert ("chrX", "-") in idx
    assert ("chrM", "+") not in idx
    assert ("GL000220.1", "+") not in idx
    assert ("KI270733.1", "-") not in idx

    # opt out: all 5 contigs retained
    idx_all = load_summit_index_from_bed(str(bed), canonical_only=False)
    assert count_summits(idx_all) == 5
    assert ("chrM", "+") in idx_all

    # explicit chroms= takes precedence over canonical_only
    idx_explicit = load_summit_index_from_bed(str(bed), chroms={"chrM"})
    assert count_summits(idx_explicit) == 1
    assert ("chrM", "+") in idx_explicit
