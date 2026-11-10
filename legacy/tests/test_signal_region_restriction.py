"""Tests for Fix B: region-restricted signal bedGraph loading.

Two surfaces:
  1. ``chroms_from_regions`` — derives the chromosome set from CLI --region
     tokens (None/[] => None = whole genome).
  2. ``load_signal_tracks`` / ``BedGraphTrack.from_file`` — when given a
     ``chroms`` set, must NOT load intervals on any other chromosome.

These pin the contract that prevents genome-wide signal from leaking into
region-restricted evaluation metrics.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from evaluation.signal_utils import (
    chroms_from_regions,
    load_signal_tracks,
    BedGraphTrack,
)


# ─── chroms_from_regions ────────────────────────────────────────────────

def test_chroms_from_regions_multi():
    assert chroms_from_regions(["chr22:1-100", "chr3:5-6"]) == {"chr22", "chr3"}


def test_chroms_from_regions_none():
    assert chroms_from_regions(None) is None


def test_chroms_from_regions_empty_list():
    assert chroms_from_regions([]) is None


def test_chroms_from_regions_dedup_same_chrom():
    # Two windows on the same chromosome collapse to one entry.
    assert chroms_from_regions(["chr22:1-100", "chr22:500-600"]) == {"chr22"}


def test_chroms_from_regions_bare_chrom_no_colon():
    # A token without ':' is treated as a whole chromosome name.
    assert chroms_from_regions(["chrX"]) == {"chrX"}


def test_chroms_from_regions_skips_empty_tokens():
    assert chroms_from_regions(["", "chr1:1-2", None]) == {"chr1"}


# ─── loader restriction ─────────────────────────────────────────────────

@pytest.fixture
def two_chrom_bedgraph(tmp_path: Path) -> Path:
    """A 4-column bedGraph with chr22 AND chrX intervals."""
    bg = tmp_path / "sig.bedgraph"
    bg.write_text(
        "chr22\t100\t200\t5.0\n"
        "chr22\t300\t400\t7.0\n"
        "chrX\t100\t200\t9.0\n"
        "chrX\t500\t600\t3.0\n"
    )
    return bg


def test_bedgraphtrack_from_file_restricts_chroms(two_chrom_bedgraph: Path):
    track = BedGraphTrack.from_file(two_chrom_bedgraph, chroms={"chr22"})
    assert "chr22" in track.data
    assert "chrX" not in track.data, "chrX intervals leaked despite chroms={'chr22'}"
    # chrX query returns 0 (no data) while chr22 returns the loaded value.
    assert track.query("chrX", 100, 200) == 0.0
    assert track.query("chr22", 100, 200) == pytest.approx(5.0)


def test_bedgraphtrack_from_file_no_restriction_loads_all(two_chrom_bedgraph: Path):
    track = BedGraphTrack.from_file(two_chrom_bedgraph, chroms=None)
    assert "chr22" in track.data
    assert "chrX" in track.data


def test_load_signal_tracks_restricts_chroms(tmp_path: Path):
    """load_signal_tracks must thread chroms through to every track so no
    out-of-region chromosome is present in any of the four returned tracks."""
    files = {}
    for name in ("cage_p", "cage_m", "qs_p", "qs_m"):
        p = tmp_path / f"{name}.bedgraph"
        p.write_text(
            "chr22\t100\t200\t4.0\n"
            "chrX\t100\t200\t8.0\n"
        )
        files[name] = str(p)

    cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
        files["cage_p"], files["cage_m"], files["qs_p"], files["qs_m"],
        chroms={"chr22"},
    )
    for track in (cage_p, cage_m, qs_p, qs_m):
        assert "chr22" in track.data
        assert "chrX" not in track.data, "chrX leaked through load_signal_tracks"
        assert track.query("chrX", 100, 200) == 0.0
