#!/usr/bin/env python3
"""Unit tests for the shared one-to-one peak matcher (the R3 fix).

These pin the eval-contract definition of a true positive (2026-05-28
direction memo §5 item 2):

  * each CAGE/dRNA peak cluster collapses to ONE summit (peak centre when
    no signal bedGraph is supplied),
  * a called end matches the NEAREST summit within ±W,
  * each summit is matched AT MOST ONCE (greedy nearest-first),
  * the window sweep ±10/±25/±50/±100 behaves monotonically.

The headline regression: two called ends near a single summit must yield
ONE true positive, not two — the recall inflator the no-filter TED was
most exposed to.

conftest.py puts bin/ on sys.path, so `from evaluation.peak_matching ...`
resolves.
"""
from __future__ import annotations

import math

from evaluation.peak_matching import (
    build_summit_index,
    load_summit_index_from_bed,
    count_summits,
    match_positions_one_to_one,
    match_ends,
    precision_recall,
    window_sweep_precision_recall,
    SWEEP_WINDOWS,
    _interval_summit,
)


# ── Summit collapse ───────────────────────────────────────────────────────

def test_summit_is_interval_centre_without_signal():
    # [950, 1050) -> centre 1000
    assert _interval_summit(950, 1050) == 1000
    assert _interval_summit(1000, 1001) == 1000


def test_summit_uses_signal_peak_when_provided():
    # Dominant signal at 1020 inside [950, 1050) overrides the centre.
    signal = [(960, 1.0), (1020, 9.0), (1040, 2.0)]
    assert _interval_summit(950, 1050, signal) == 1020


def test_summit_signal_tie_breaks_to_centre():
    # Two equal maxima; the one closer to the centre (1000) wins.
    signal = [(970, 5.0), (1030, 5.0)]
    # |970-1000| == |1030-1000| == 30 -> deterministic: first seen (970)
    # since center-dist ties keep the earlier; assert it is one of them and
    # equidistant (the contract only requires determinism).
    s = _interval_summit(950, 1050, signal)
    assert s in (970, 1030)
    assert abs(s - 1000) == 30


def test_build_summit_index_one_summit_per_peak():
    peaks = [
        ("chr1", 950, 1050, "+"),
        ("chr1", 2000, 2100, "+"),
        ("chr1", 3000, 3100, "-"),
    ]
    idx = build_summit_index(peaks)
    assert count_summits(idx) == 3
    assert idx[("chr1", "+")] == [1000, 2050]
    assert idx[("chr1", "-")] == [3050]


# ── One-to-one greedy matching ────────────────────────────────────────────

def test_two_calls_one_summit_yields_one_tp():
    """THE headline R3 test: two called ends near one summit -> ONE TP."""
    summits = [1000]
    # Two calls within +-50 of the single summit.
    positions = [990, 1010]
    assign = match_positions_one_to_one(positions, summits, window=50)
    matched = [j for j in assign if j >= 0]
    assert len(matched) == 1, "only one call may claim the single summit"
    # Greedy-nearest: 990 (dist 10) and 1010 (dist 10) tie; deterministic
    # tie-break (call index) gives the summit to call 0.
    assert assign[0] == 0
    assert assign[1] == -1


def test_greedy_assigns_nearest_first():
    """When two summits and two calls compete, nearest pair wins first."""
    summits = [1000, 1100]
    # call A=1005 (dist 5 to s0, 95 to s1); call B=1098 (dist 2 to s1)
    positions = [1005, 1098]
    assign = match_positions_one_to_one(positions, summits, window=50)
    assert assign[0] == 0  # A -> summit 1000
    assert assign[1] == 1  # B -> summit 1100
    assert len([j for j in assign if j >= 0]) == 2


def test_greedy_steals_nearer_call_for_shared_summit():
    """Two calls near the same summit: the NEARER call claims it, the
    farther call is left unmatched (no second summit in range)."""
    summits = [1000]
    positions = [1040, 1003]  # 1003 is nearer (dist 3) than 1040 (dist 40)
    assign = match_positions_one_to_one(positions, summits, window=50)
    assert assign[1] == 0          # nearer call (1003) wins
    assert assign[0] == -1         # farther call (1040) unmatched
    assert sum(1 for j in assign if j >= 0) == 1


def test_unmatched_outside_window():
    summits = [1000]
    positions = [1200]
    assign = match_positions_one_to_one(positions, summits, window=50)
    assert assign == [-1]


def test_empty_inputs():
    assert match_positions_one_to_one([], [1000], 50) == []
    assert match_positions_one_to_one([1000], [], 50) == [-1]


# ── Bucketed match_ends (chrom/strand aware) ──────────────────────────────

def test_match_ends_strand_aware():
    idx = build_summit_index([("chr1", 950, 1050, "+")])  # summit 1000 on +
    ends = [
        ("chr1", "+", 1005),   # matches
        ("chr1", "-", 1005),   # wrong strand -> no summit there
    ]
    matched = match_ends(ends, idx, window=50)
    assert matched[0] == "chr1:+:1000"
    assert matched[1] is None


def test_match_ends_one_summit_one_call_across_bucket():
    """A cluster of 5 calls near one summit in the same bucket -> 1 match."""
    idx = build_summit_index([("chr1", 950, 1050, "+")])
    ends = [("chr1", "+", p) for p in (980, 990, 1000, 1010, 1020)]
    matched = match_ends(ends, idx, window=50)
    assert sum(1 for m in matched if m is not None) == 1


# ── precision / recall accounting ─────────────────────────────────────────

def test_precision_recall_cluster_not_inflated():
    """5 calls on 1 summit: TP=1, recall=1/1, precision=1/5."""
    idx = build_summit_index([("chr1", 950, 1050, "+")])
    ends = [("chr1", "+", p) for p in (980, 990, 1000, 1010, 1020)]
    pr = precision_recall(ends, idx, window=50)
    assert pr["tp"] == 1
    assert pr["n_calls"] == 5
    assert pr["n_summits"] == 1
    assert pr["recall"] == 1.0
    assert math.isclose(pr["precision"], 1 / 5)


def test_precision_recall_two_summits_two_calls():
    idx = build_summit_index([
        ("chr1", 950, 1050, "+"),    # summit 1000
        ("chr1", 1950, 2050, "+"),   # summit 2000
    ])
    ends = [("chr1", "+", 1005), ("chr1", "+", 2003)]
    pr = precision_recall(ends, idx, window=50)
    assert pr["tp"] == 2
    assert pr["recall"] == 1.0
    assert pr["precision"] == 1.0


# ── window sweep ──────────────────────────────────────────────────────────

def test_window_sweep_is_monotonic_in_tp():
    """Wider windows can only match the same or more calls/summits."""
    idx = build_summit_index([("chr1", 1000, 1001, "+")])  # summit 1000
    # calls at increasing distance: 5, 20, 40, 80 bp away.
    ends = [("chr1", "+", 1000 + d) for d in (5, 20, 40, 80)]
    sweep = window_sweep_precision_recall(ends, idx)
    tps = [sweep[w]["tp"] for w in SWEEP_WINDOWS]  # windows 10/25/50/100
    # Still one-to-one (only one summit) so tp caps at 1 in every window,
    # but recall must be non-decreasing and reach 1 by +-100.
    assert tps == [1, 1, 1, 1]
    # Now with 4 separate summits, one per call distance, tp grows with W.
    idx2 = build_summit_index([
        ("chr1", 1005, 1006, "+"),
        ("chr1", 1020, 1021, "+"),
        ("chr1", 1040, 1041, "+"),
        ("chr1", 1080, 1081, "+"),
    ])
    ends2 = [("chr1", "+", 1000)]  # single call at 1000
    sweep2 = window_sweep_precision_recall(ends2, idx2)
    tp2 = [sweep2[w]["tp"] for w in SWEEP_WINDOWS]
    # one call can match at most one summit -> tp capped at 1, but it only
    # finds a summit once W >= 5.
    assert sweep2[10]["tp"] == 1
    # recall denominators are 4 summits, so recall is 1/4 at most here.
    assert math.isclose(sweep2[10]["recall"], 1 / 4)


def test_window_recall_nondecreasing_with_more_calls():
    idx = build_summit_index([
        ("chr1", 1005, 1006, "+"),
        ("chr1", 1020, 1021, "+"),
        ("chr1", 1040, 1041, "+"),
        ("chr1", 1080, 1081, "+"),
    ])
    # 4 calls, each exactly on a summit -> all match at every window >= 0.
    ends = [("chr1", "+", p) for p in (1005, 1020, 1040, 1080)]
    sweep = window_sweep_precision_recall(ends, idx)
    recalls = [sweep[w]["recall"] for w in SWEEP_WINDOWS]
    assert recalls == [1.0, 1.0, 1.0, 1.0]


# ── BED loader on the tiny fixtures ───────────────────────────────────────

def test_load_summit_index_from_tiny_cage(data_dir):
    idx = load_summit_index_from_bed(str(data_dir / "tiny_cage_peaks.bed"))
    # 6 CAGE peaks in the fixture -> 6 summits.
    assert count_summits(idx) == 6
    # cage_peak_1 at [950,1050) on + -> summit 1000.
    assert 1000 in idx[("chr_test", "+")]
