#!/usr/bin/env python3
"""Unit tests for analysis/novel_end_pr/novel_end_pr.py — the novel-end carve-out.

Pins the carve-out semantics:
  * an orthogonal summit (or a called end) is NOVEL iff it is strictly MORE than
    `dist` bp from the nearest annotated end on its (chrom, strand);
  * filter_index_novel drops non-novel summits and keeps the index sorted;
  * the PR it feeds (evaluation.peak_matching.precision_recall) is unchanged, so
    a novel call on a novel summit scores TP=1 / P=R=1, and an annotated call
    that is filtered out cannot inflate novel precision.

conftest.py puts bin/ on sys.path; we additionally load the standalone analysis
module by path (it lives outside the package tree).
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest

_NEP_PATH = (Path(__file__).resolve().parent.parent
             / "analysis" / "novel_end_pr" / "novel_end_pr.py")


@pytest.fixture(scope="module")
def nep():
    spec = importlib.util.spec_from_file_location("novel_end_pr", _NEP_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── distance + novelty primitives ──────────────────────────────────────────────

def test_min_dist_basic(nep):
    arr = [100, 500, 5000]
    assert nep._min_dist(100, arr) == 0
    assert nep._min_dist(120, arr) == 20      # nearest is 100
    assert nep._min_dist(480, arr) == 20      # nearest is 500
    assert nep._min_dist(10_000, arr) == 5000
    assert nep._min_dist(0, []) == math.inf   # no annotation → infinitely far


def test_is_novel_threshold_is_strict(nep):
    ann = {("chr1", "+"): [1000]}
    # exactly 100 bp away is NOT novel (strictly greater than dist required)
    assert nep._is_novel("chr1", "+", 1100, ann, 100) is False
    assert nep._is_novel("chr1", "+", 1101, ann, 100) is True
    # other strand / chrom with no annotation is novel
    assert nep._is_novel("chr1", "-", 1000, ann, 100) is True
    assert nep._is_novel("chrX", "+", 1000, ann, 100) is True


def test_filter_index_novel_keeps_only_far_summits_sorted(nep):
    ann = {("chr1", "+"): [1000, 5000]}
    idx = {("chr1", "+"): [1010, 3000, 4950], ("chr2", "+"): [10]}
    out = nep.filter_index_novel(idx, ann, 100)
    # 1010 (10bp) and 4950 (50bp) drop; 3000 (1990bp) stays; chr2 has no annot
    assert out[("chr1", "+")] == [3000]
    assert out[("chr2", "+")] == [10]
    # sortedness preserved for downstream bisect-based matching
    assert out[("chr1", "+")] == sorted(out[("chr1", "+")])


# ── carve-out feeds the unchanged PR correctly ─────────────────────────────────

def test_precision_recall_on_novel_subset(nep):
    ann = {("chr1", "+"): [1000]}
    full_index = nep.load_summit_index_from_bed  # noqa: F841  (sanity: importable)
    # one novel summit at 3000, one annotated-adjacent summit at 1010
    idx = {("chr1", "+"): [1010, 3000]}
    novel_index = nep.filter_index_novel(idx, ann, 100)
    assert nep.count_summits(novel_index) == 1

    # calls: one novel (3000, matches), one near-annotation (1010, filtered out)
    all_calls = [("chr1", "+", 1010), ("chr1", "+", 3000)]
    novel_calls = [c for c in all_calls
                   if nep._is_novel(c[0], c[1], c[2], ann, 100)]
    assert novel_calls == [("chr1", "+", 3000)]

    pr = nep.precision_recall(novel_calls, novel_index, 50)
    assert pr["tp"] == 1
    assert pr["precision"] == 1.0      # the one novel call matched
    assert pr["recall"] == 1.0         # the one novel summit recovered

    # If we (wrongly) kept all calls against the novel index, the annotated call
    # would be an unmatched FP — precision must fall. This is why the default
    # carve-out filters calls too.
    pr_calls_unfiltered = nep.precision_recall(all_calls, novel_index, 50)
    assert pr_calls_unfiltered["tp"] == 1
    assert pr_calls_unfiltered["precision"] == 0.5
