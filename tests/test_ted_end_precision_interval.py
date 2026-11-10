#!/usr/bin/env python3
"""Interval-predicate tests for ted_end_precision's own helpers.

ted_end_precision.py routes the dedup/paired metric through peak_matching
(auto-fixed) but keeps its OWN single-end helper `_many_to_one_match_counts`
(the single-end headline) and the paired-dedup `_nearest_annot`.

The 2026-06-21 fix makes the SINGLE-END many-to-one path use the same
within-W-of-INTERVAL predicate as peak_matching._match_distance / the shared
matcher / stratified_pr (so the single-end 5'/3' headline reconciles across all
eval scripts):

    single-end:  call c matches peak p  iff  _match_distance(c, p) <= W
                 (supported zone [lo-W, (hi-1)+W]; credits a call within W of an
                  EDGE of a wide peak, not only within W of the summit).

The PAIRED-dedup `_nearest_annot` retains its own predicate
    (|c - summit| <= W) OR (start <= c < end)
because Fig C / the paired headline are built on it deliberately.

These tests pin:
  * single-end many-to-one: a wide-peak interior call is credited, and an
    edge-near call (within W of an interval edge but >W from summit) is now
    credited under the corrected predicate;
  * paired-dedup `_nearest_annot`: interior credit on a wide peak, summit±W on
    narrow peaks (unchanged), and a stable (summit,start,end) identity so two
    interior calls collapse to the SAME peak.

conftest.py puts bin/ on sys.path.
"""
from __future__ import annotations

import math
from collections import defaultdict

from evaluation.peak_matching import build_summit_index
from evaluation.ted_end_precision import (
    _many_to_one_match_counts,
    _nearest_annot,
    _summit_index_from_intervals,
    parse_peaks_bed,
)


def _iso(chrom, strand, tss, tts):
    return {"chrom": chrom, "strand": strand, "tss": tss, "tts": tts}


# ── single-end many-to-one: wide-peak interior credit ──────────────────────

def test_many_to_one_wide_peak_interior_credit():
    """A wide peak: a summit-far interior call is now BOTH a precision hit and
    recovers the peak (recall), where old summit±W would credit neither."""
    W = 50
    # peak [1000, 1500) width 500, summit 1250
    peaks = {("chr1", "+"): [(1000, 1500)]}
    idx = _summit_index_from_intervals(peaks)
    # call at 1100: 150 bp from summit (>W) but inside -> hit
    isos = [_iso("chr1", "+", 1100, 0)]
    ends_near, n_calls, peaks_rec, n_summits = _many_to_one_match_counts(
        isos, idx, "tss", W)
    assert ends_near == 1
    assert n_calls == 1
    assert peaks_rec == 1
    assert n_summits == 1


def test_many_to_one_within_W_of_edge_credited():
    """CORRECTED predicate: a call within W of an interval EDGE is credited even
    when it is >W from the summit (the documented [lo-W, (hi-1)+W] zone)."""
    W = 50
    # peaks: summit 1000 [970,1030) (w 60), summit 2000 [1950,2050) (w 100)
    peaks = {("chr1", "+"): [(970, 1030), (1950, 2050)]}
    idx = _summit_index_from_intervals(peaks)
    isos = [
        _iso("chr1", "+", 1005, 0),   # inside -> hit
        _iso("chr1", "+", 1049, 0),   # within W of summit AND edge -> hit
        _iso("chr1", "+", 1051, 0),   # 51 from summit but 22 from right edge -> HIT (corrected)
        _iso("chr1", "+", 2000, 0),   # on summit -> hit
        _iso("chr1", "+", 2051, 0),   # 51 from summit but 2 from right edge -> HIT (corrected)
        _iso("chr1", "+",  919, 0),   # 81 from summit, 51 from left edge 970 -> MISS
    ]
    ends_near, n_calls, peaks_rec, n_summits = _many_to_one_match_counts(
        isos, idx, "tss", W)
    # within-W-of-interval reference, computed directly via _match_distance:
    from evaluation.peak_matching import _match_distance
    ref_ends_near = sum(
        1 for i in isos
        if any(_match_distance(i["tss"], p) <= W for p in idx[("chr1", "+")]))
    assert ends_near == ref_ends_near == 5   # all but the 919 left-edge miss
    assert peaks_rec == 2
    assert n_calls == 6
    assert n_summits == 2


def test_many_to_one_two_interior_calls_recall_one_peak():
    """Many-to-one ALLOWS duplicate precision (single-end policy) but a peak is
    recovered ONCE: two interior calls -> ends_near=2, peaks_recovered=1."""
    W = 50
    peaks = {("chr1", "+"): [(1000, 1500)]}
    idx = _summit_index_from_intervals(peaks)
    isos = [_iso("chr1", "+", 1100, 0), _iso("chr1", "+", 1400, 0)]
    ends_near, n_calls, peaks_rec, n_summits = _many_to_one_match_counts(
        isos, idx, "tss", W)
    assert ends_near == 2          # both interior calls credited (duplication ok)
    assert peaks_rec == 1          # the single peak is recovered once
    assert n_summits == 1


# ── paired-dedup _nearest_annot: identity + interval credit ─────────────────

def test_nearest_annot_wide_peak_interior():
    W = 50
    intervals = [(1000, 1500)]     # summit 1250
    # interior, summit-far -> matches, returns (summit, start, end)
    key = _nearest_annot(1100, intervals, W)
    assert key == (1250, 1000, 1500)


def test_nearest_annot_narrow_peak_backcompat():
    """Narrow peak: summit±W behaviour, miss just outside the window."""
    W = 50
    intervals = [(970, 1030)]      # summit 1000, width 60
    assert _nearest_annot(1049, intervals, W) == (1000, 970, 1030)
    assert _nearest_annot(1051, intervals, W) is None    # 51 from summit -> miss


def test_nearest_annot_stable_identity_two_interior_calls():
    """Two interior calls on the SAME wide peak return the SAME identity, so
    paired-dedup collapses them (one unique pair)."""
    W = 50
    intervals = [(1000, 1500)]
    k1 = _nearest_annot(1100, intervals, W)
    k2 = _nearest_annot(1400, intervals, W)
    assert k1 == k2
    assert len({k1, k2}) == 1


def test_nearest_annot_picks_interior_peak_over_window_peak():
    """When a call is interior to one peak and within ±W of another summit, the
    interior peak (match-distance 0) wins."""
    W = 50
    # peak A interior [1000, 1400) summit 1200; peak B narrow [1430, 1450)
    # summit 1440.  Call at 1420 is 20 bp from B's summit (window hit) but is
    # NOT interior to A (1420 < 1400 is false).  So B should win here.
    intervals = [(1000, 1400), (1430, 1450)]
    # call at 1390 IS interior to A and 50 from B-summit edge:
    key = _nearest_annot(1390, intervals, W)
    assert key == (1200, 1000, 1400)     # interior to A wins


# ── canonical-chrom filtering on parse_peaks_bed ────────────────────────────

def test_parse_peaks_bed_canonical_default(tmp_path):
    bed = tmp_path / "p.bed"
    bed.write_text(
        "chr1\t100\t200\tp\t0\t+\n"
        "chrM\t10\t20\tpm\t0\t+\n"
        "GL000220.1\t5\t9\tpg\t0\t-\n"
    )
    peaks = parse_peaks_bed(str(bed))            # canonical default
    assert ("chr1", "+") in peaks
    assert ("chrM", "+") not in peaks
    assert ("GL000220.1", "-") not in peaks

    peaks_all = parse_peaks_bed(str(bed), canonical_only=False)
    assert ("chrM", "+") in peaks_all
    assert ("GL000220.1", "-") in peaks_all
