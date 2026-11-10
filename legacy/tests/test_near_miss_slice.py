#!/usr/bin/env python3
"""Unit tests for the near-miss eval slice (FIX 3 / PI-synthesis A5).

Covers BOTH near-miss falloff views in bin/evaluation/near_miss_slice.py:

  1. The pure falloff math over model scores (summarize_falloff /
     localization_score) — a sharp localizer vs a flat region detector.
  2. The isoform-end-driven slice (bin_ends_by_signed_distance) — taking the
     ACTUAL emitted ends, binning by signed distance to the nearest truth
     summit (one-to-one matched), and reporting mean TED reality + per-bin
     precision (the end-side near-miss falloff the review asks for).

All synthetic; no genome, model pickle, or chr22 run required. conftest puts
bin/ on sys.path so `evaluation.*` resolves.
"""
from __future__ import annotations

import math

from evaluation.near_miss_slice import (
    signed_offsets,
    summarize_falloff,
    localization_score,
    signed_distance_to_nearest_summit,
    bin_ends_by_signed_distance,
    _bin_label,
    DEFAULT_DIST_BIN_EDGES,
)


# ── 1. Pure model-score falloff math ───────────────────────────────────────

def test_signed_offsets_symmetric_and_includes_zero():
    assert signed_offsets([5, 10]) == [-10, -5, 0, 5, 10]
    assert signed_offsets([]) == [0]
    # dedupe + sort
    assert signed_offsets([10, 10, 5]) == [-10, -5, 0, 5, 10]


def test_sharp_localizer_scores_high():
    # True site scores 0.9, decoys fall off hard -> good localization (~1.0).
    sbo = {
        0:   [0.90, 0.92, 0.88],
        10:  [0.40, 0.35, 0.42],
        -10: [0.38, 0.41, 0.39],
    }
    summ = summarize_falloff(sbo)
    # decoys score well BELOW true -> negative delta, frac_ge_true ~ 0
    assert summ[10]["mean_delta_from_true"] < 0
    assert summ[10]["frac_ge_true"] == 0.0
    assert localization_score(summ) == 1.0


def test_flat_region_detector_scores_near_half():
    # Decoys score ~ as high as the true site -> flat curve, ~0.5 localization.
    sbo = {
        0:   [0.80, 0.82, 0.78],
        10:  [0.81, 0.79, 0.80],
        -10: [0.79, 0.83, 0.77],
    }
    summ = summarize_falloff(sbo)
    assert abs(summ[10]["mean_delta_from_true"]) < 0.05
    loc = localization_score(summ)
    assert 0.2 <= loc <= 0.8  # not a clean localizer


def test_summarize_requires_true_offset():
    try:
        summarize_falloff({10: [0.5]})
    except ValueError:
        return
    raise AssertionError("expected ValueError when offset 0 missing")


# ── 2. Isoform-end-driven near-miss falloff ────────────────────────────────

def test_signed_distance_sign_convention():
    summits = [1000]
    # + strand: call downstream of summit -> positive
    assert signed_distance_to_nearest_summit([1010], summits, "+") == [10]
    assert signed_distance_to_nearest_summit([990], summits, "+") == [-10]
    # - strand: sign flips (downstream in transcription dir = lower coord)
    assert signed_distance_to_nearest_summit([1010], summits, "-") == [-10]
    assert signed_distance_to_nearest_summit([990], summits, "-") == [10]


def test_bin_label_edges():
    edges = DEFAULT_DIST_BIN_EDGES
    assert _bin_label(edges, 0) == "[0,5)"
    assert _bin_label(edges, 7) == "[5,10)"
    assert _bin_label(edges, -3) == "[-5,0)"
    assert _bin_label(edges, 999) == ">=100"
    assert _bin_label(edges, -999) == "<-100"


def test_bin_ends_precision_falls_off_with_distance():
    """A near-miss falloff: ends ON summits are TPs with high reality; ends far
    from summits are not TPs and carry low reality. Per-bin precision must drop
    as |distance| grows."""
    # Three truth summits on chr1 +; well separated so one-to-one is clean.
    summit_index = {("chr1", "+"): [1000, 2000, 3000]}
    # Called ends:
    #  - exactly on each summit (dist 0)            -> TP, high reality
    #  - one 8 bp off summit 1000 (within W=50)     -> TP (one-to-one already
    #       took 1000 for the exact call, so this 1008 falls to the next summit?
    #       no: separate summits; greedy gives exact calls their own summit and
    #       1008 has no remaining nearby summit -> not a TP) low reality
    #  - one 200 bp off any summit                  -> not a TP, low reality
    ends = [
        ("chr1", "+", 1000),   # on summit
        ("chr1", "+", 2000),   # on summit
        ("chr1", "+", 3000),   # on summit
        ("chr1", "+", 1008),   # near-miss 8bp (summit already claimed)
        ("chr1", "+", 1200),   # far miss 200bp
    ]
    reality = [0.9, 0.85, 0.88, 0.30, 0.10]
    bins = bin_ends_by_signed_distance(
        ends, summit_index, window=50, reality=reality)

    # The on-summit bin [0,5) holds the 3 exact calls, all one-to-one TPs.
    on = bins["[0,5)"]
    assert on["n"] == 3
    assert on["precision"] == 1.0
    assert on["mean_reality"] > 0.8

    # The 8bp near-miss lands in [5,10): not a TP (its summit was claimed),
    # low reality.
    near = bins["[5,10)"]
    assert near["n"] == 1
    assert near["precision"] == 0.0
    assert near["mean_reality"] < 0.5

    # The 200bp far miss lands in >=100: not a TP, lowest reality.
    far = bins[">=100"]
    assert far["n"] == 1
    assert far["precision"] == 0.0
    assert far["mean_reality"] < 0.2

    # Falloff: precision monotone non-increasing from on-summit outward.
    assert on["precision"] >= near["precision"] >= far["precision"]


def test_bin_ends_reality_optional():
    summit_index = {("chr1", "+"): [1000]}
    ends = [("chr1", "+", 1000), ("chr1", "+", 1050)]
    bins = bin_ends_by_signed_distance(ends, summit_index, window=10)
    # no reality vector -> mean_reality NaN, precision/dist still computed
    for b in bins.values():
        assert math.isnan(b["mean_reality"])
        assert math.isfinite(b["precision"])


def test_bin_ends_no_summit_bucket_ignored():
    # Ends on a (chrom, strand) with no summits are not near-miss observations.
    summit_index = {("chr1", "+"): [1000]}
    ends = [("chr2", "+", 5000), ("chr1", "+", 1000)]
    bins = bin_ends_by_signed_distance(ends, summit_index, window=50)
    total = sum(int(b["n"]) for b in bins.values())
    assert total == 1  # only the chr1 end counted


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-v", "--tb=short"]))
