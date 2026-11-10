#!/usr/bin/env python3
"""Regression tests for sjc_alt_end_analysis._analyse_assembler classification.

The bug we fixed (May 8, 2026): the original `else: n_fully_redundant += 1`
branch caught EVERY case where neither end was unique within the SJC group,
without checking whether the shared (p5, p3) tuple was actually
ortho-validated. Two isoforms sharing (CAGE_X, None) were being classified
as `fully_redundant` even though their `None` 3'-end signature meant
neither was duplicating an ortho-validated coordinate. Those rows now
route into `fully_unsup`.

These tests pin the corrected classification using a deliberately small
fixture with known input → known output:

  tests/data/tiny_redundant_isoforms.bed (4 isoforms, 1 SJC group):
    REDA / REDB: tss=1000, tts ∈ {3500, 3510}
                 → p5 = cage_peak_1, p3 = drna_peak_2  (both supported)
    REDC / REDD: tss=1000, tts ∈ {6000, 6100}
                 → p5 = cage_peak_1, p3 = None         (3' missing dRNA)

Pre-fix: the classifier would report n_fully_redundant=4, n_fully_unsup=0.
Post-fix: the classifier correctly reports n_fully_redundant=2 (REDA+REDB),
n_fully_unsup=2 (REDC+REDD).

If anyone reverts the classifier or breaks the (no_5p or no_3p) gate
in the bottom branch, this test catches it.
"""
from __future__ import annotations

import pytest

# conftest.py adds bin/ to sys.path
from evaluation.sjc_alt_end_analysis import _analyse_assembler
from evaluation.signal_utils import parse_isoforms
from evaluation.ted_end_precision import parse_peaks_bed


@pytest.fixture(scope="module")
def fixture_data(tiny_redundant_isoforms_bed, tiny_cage_peaks, tiny_drna_peaks):
    """Parse the fixture once for the whole test module."""
    isoforms = parse_isoforms(str(tiny_redundant_isoforms_bed))
    # Fixtures live on the synthetic 'chr_test' contig, which the new
    # canonical-chrom default (chr1-22,X,Y) would drop — opt out here.
    cage_peaks = parse_peaks_bed(str(tiny_cage_peaks), canonical_only=False)
    qs_peaks = parse_peaks_bed(str(tiny_drna_peaks), canonical_only=False)
    return isoforms, cage_peaks, qs_peaks


@pytest.fixture(scope="module")
def analysis_result(fixture_data):
    """Run _analyse_assembler once; reuse the result across assertions."""
    isoforms, cage_peaks, qs_peaks = fixture_data
    return _analyse_assembler("test_label", isoforms, cage_peaks, qs_peaks)


# ── Sanity checks ────────────────────────────────────────────────────


def test_fixture_parses_4_isoforms(fixture_data):
    """The fixture must contain exactly 4 isoforms."""
    isoforms, _, _ = fixture_data
    assert len(isoforms) == 4, f"expected 4 isoforms, got {len(isoforms)}"


def test_fixture_groups_into_one_sjc(fixture_data):
    """All 4 isoforms must share one (chrom, strand, junctions) key."""
    isoforms, _, _ = fixture_data
    from evaluation.signal_utils import group_by_junction_chain
    groups = group_by_junction_chain(isoforms)
    assert len(groups) == 1, (
        f"expected 1 SJC group (all 4 isos share intron 1200-3000), "
        f"got {len(groups)} groups: {list(groups.keys())}"
    )
    members = next(iter(groups.values()))
    assert len(members) == 4


# ── The core regression assertions ──────────────────────────────────


def test_no_isoform_is_fully_unique(analysis_result):
    """Every isoform shares a (p5, p3) signature with at least one
    sibling, so n_fully_unique must be 0."""
    assert analysis_result["n_fully_unique"] == 0


def test_redundant_count_is_2_post_fix(analysis_result):
    """REGRESSION: n_fully_redundant must equal 2 (the REDA+REDB pair
    that genuinely shares (cage_peak_1, drna_peak_2)). Pre-fix this
    would have been 4 — anyone re-introducing the old `else: redundant`
    branch will trip this assertion."""
    assert analysis_result["n_fully_redundant"] == 2, (
        f"expected n_fully_redundant=2 (only REDA+REDB share an "
        f"ortho-validated peak pair), got "
        f"n_fully_redundant={analysis_result['n_fully_redundant']}. "
        f"This usually means the (no_5p or no_3p) gate in "
        f"_analyse_assembler's bottom else branch was removed, "
        f"re-routing redund-with-None cases back to fully_redundant."
    )


def test_unsup_count_is_2_post_fix(analysis_result):
    """REGRESSION: n_fully_unsup must absorb the 2 isos sharing
    (cage_peak_1, None) — REDC + REDD. Pre-fix this would have been 0
    (those isoforms were mis-routed to fully_redundant)."""
    assert analysis_result["n_fully_unsup"] == 2, (
        f"expected n_fully_unsup=2 (the REDC+REDD pair has p3=None "
        f"on both isos and should be classified as unsupported), got "
        f"n_fully_unsup={analysis_result['n_fully_unsup']}. "
        f"This means the (no_5p or no_3p) gate was removed and "
        f"redund-with-None isoforms are being mis-classified again."
    )


def test_one_axis_only_unique_counts_are_zero(analysis_result):
    """In this fixture, p5 is identical (cage_peak_1) for all 4 isoforms,
    so no isoform's 5p is unique. p3 is shared in pairs (drna_peak_2 for
    REDA+REDB; None for REDC+REDD), so no isoform's 3p is unique either.
    Both n_unique_5p_only and n_unique_3p_only must be 0."""
    assert analysis_result["n_unique_5p_only"] == 0
    assert analysis_result["n_unique_3p_only"] == 0


def test_classification_is_total(analysis_result):
    """The 5 categories must sum to exactly 4 (every isoform
    classified into exactly one bucket)."""
    total = (
        analysis_result["n_fully_unique"]
        + analysis_result["n_unique_5p_only"]
        + analysis_result["n_unique_3p_only"]
        + analysis_result["n_fully_redundant"]
        + analysis_result["n_fully_unsup"]
    )
    assert total == 4, (
        f"all 4 isoforms in the fixture must be classified into "
        f"exactly one bucket; got total={total}. Categories: "
        f"unique={analysis_result['n_fully_unique']}, "
        f"5p_only={analysis_result['n_unique_5p_only']}, "
        f"3p_only={analysis_result['n_unique_3p_only']}, "
        f"redundant={analysis_result['n_fully_redundant']}, "
        f"unsup={analysis_result['n_fully_unsup']}"
    )


# ── Panel C semantics: pairwise distance lists ──────────────────────


def test_panel_c_redundant_pairs_only_from_supported(analysis_result):
    """Panel C's `pairs_fully_redundant` list must include ONLY
    pairs where both p5 AND p3 are non-None and shared. The REDA-REDB
    pair (CAGE+dRNA shared) goes here; the REDC-REDD pair (CAGE shared,
    None p3) goes to `pairs_partial` because Panel C's same_5p / same_3p
    check requires non-None matches."""
    pairs = analysis_result["pairs_fully_redundant"]
    assert len(pairs) == 1, (
        f"expected 1 pair (REDA-REDB) in pairs_fully_redundant, "
        f"got {len(pairs)} pairs: {pairs}"
    )
    # The (5p_dist, 3p_dist) for REDA(tss=1000,tts=3500)
    # vs REDB(tss=1000,tts=3510): d5=0, d3=10
    d5, d3 = pairs[0]
    assert d5 == 0
    assert d3 == 10


def test_panel_c_partial_includes_none_signature(analysis_result):
    """The REDC-REDD pair shares (CAGE, None) — their 5p is shared (so
    same_5p=True) but their 3p is None on both, so same_3p=False (the
    Panel C check requires non-None equality). They route into
    pairs_partial. Plus the cross pairs (REDA<->REDC etc.) where 5p is
    shared but 3p differs end up there too. Total: 5 cross-pair entries
    (one for REDC-REDD, four for REDA/REDB × REDC/REDD = 4 cross-pairs)."""
    partial = analysis_result["pairs_partial"]
    # REDA-REDC, REDA-REDD, REDB-REDC, REDB-REDD, plus REDC-REDD itself
    # = 5 pairs total in pairs_partial.
    assert len(partial) == 5, (
        f"expected 5 partial pairs, got {len(partial)}: {partial}"
    )
