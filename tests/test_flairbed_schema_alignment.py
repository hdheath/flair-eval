"""Tests that eval scripts read the NEW 30-column FlairBed schema correctly.

The new FLAIR/TED ``*.isoforms.bed`` is a 30-column FlairBed (see
``flair-ted-52826/src/flair/flair_bed.as``):

  col 7/8  thickStart/thickEnd  -> CDS genomic span (INLINE; the separate
                                   *.isoforms.CDS.bed no longer exists)
  col 13   gene_id (FLG...)     -> the real gene id (col 4 is a BARE FLT id,
                                   NOT a gene id)
  col 16   read_support
  col 17   frac_support
  col 18   productivity         -> PRO/PTC/NGO/NST (DEDICATED column; NOT
                                   embedded in the name)
  col 23-30 TED confidence      -> may be empty strings when TED not run

These tests pin down that:
  * signal_utils.parse_bed12 reads gene_id from col 13, productivity from
    col 18, CDS from 7/8, read_support from 16, and tolerates empty TED cols;
  * alt_end_cds_analysis / cds_coverage_recovery now produce non-empty data
    from the regular isoforms.bed (they were silently empty under the old
    name-regex productivity parse);
  * plain (non-FlairBed) BED12 still falls back to gene_from_name.

Fixture: tests/data/tiny_flairbed_30col.bed
  FLG0001: FLT0001 PRO, FLT0002 PTC, FLT0003 NGO  (chr1 +, shared SJC)
  FLG0002: FLT0004 PRO, FLT0005 NST               (chr2 -, shared SJC)
  TED cols populated on FLT0001 & FLT0004, empty on the rest.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parent.parent / "bin"
for p in (BIN, BIN / "evaluation"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

DATA = Path(__file__).resolve().parent / "data"
FLAIRBED = DATA / "tiny_flairbed_30col.bed"
PLAIN_BED12 = DATA / "tiny_isoforms.bed"  # 12-col, non-FlairBed

import signal_utils as su  # noqa: E402
import alt_end_cds_analysis as cds  # noqa: E402
import cds_coverage_recovery as ccr  # noqa: E402


# ── signal_utils.parse_bed12 reads the inline FlairBed columns ──────────────

def test_parse_bed12_reads_inline_flairbed_fields():
    isos = su.parse_bed12(FLAIRBED)
    assert len(isos) == 5
    by_name = {i["name"]: i for i in isos}

    f1 = by_name["FLT0001"]
    # gene_id from col 13 (FLG...), NOT the bare FLT name
    assert f1["gene_id"] == "FLG0001"
    assert f1["is_flair_bed"] is True
    # productivity from col 18
    assert f1["productivity"] == "PRO"
    # CDS from cols 7/8
    assert f1["thick_start"] == 1100
    assert f1["thick_end"] == 4800
    # read_support col 16, frac_support col 17
    assert f1["read_support"] == 19
    assert abs(f1["frac_support"] - 0.61) < 1e-9
    # ref ids cols 14/15 (trailing comma stripped)
    assert f1["ref_transcript_id"] == "ENST001.1"
    assert f1["ref_gene_mapping"] == "ENSG001.1"
    # TED cols 23-30 populated
    assert f1["ted"]["TSS_reality"] == "0.85"
    assert f1["ted"]["TTS_reality"] == "0.55"


def test_parse_bed12_distinct_genes_not_per_isoform():
    """The bug we fixed: gene_from_name on bare FLT names made every
    isoform its own gene. Col-13 grouping yields the real 2 genes."""
    isos = su.parse_bed12(FLAIRBED)
    assert {i["gene_id"] for i in isos} == {"FLG0001", "FLG0002"}


def test_parse_bed12_productivity_classes_present():
    isos = su.parse_bed12(FLAIRBED)
    prods = sorted(i["productivity"] for i in isos)
    assert prods == ["NGO", "NST", "PRO", "PRO", "PTC"]


def test_parse_bed12_empty_ted_cols_are_blank_strings():
    isos = su.parse_bed12(FLAIRBED)
    f2 = next(i for i in isos if i["name"] == "FLT0002")
    # TED not run for this isoform -> all 8 TED cols empty strings
    assert set(f2["ted"].values()) == {""}
    # productivity is still read (col 18 present even without TED)
    assert f2["productivity"] == "PTC"


def test_parse_bed12_plain_bed12_falls_back_to_name():
    """Non-FlairBed 12-col input must NOT be treated as FlairBed; gene_id
    falls back to gene_from_name and productivity is None."""
    if not PLAIN_BED12.exists():
        pytest.skip("tiny_isoforms.bed fixture missing")
    isos = su.parse_bed12(PLAIN_BED12)
    assert isos, "plain BED12 should still parse"
    for i in isos:
        assert i["is_flair_bed"] is False
        assert i["productivity"] is None


# ── alt_end_cds_analysis reads CDS from 7/8 and productivity from col 18 ─────

def test_alt_end_cds_parser_inline_productivity_and_cds():
    isos = cds.parse_bed12_cds(FLAIRBED)
    by_name = {i["name"]: i for i in isos}
    f1 = by_name["FLT0001"]
    assert f1["productivity"] == "PRO"          # col 18, not name regex
    assert f1["gene_id"] == "FLG0001"           # col 13
    assert (f1["thick_start"], f1["thick_end"]) == (1100, 4800)
    # NGO: productivity present but thick == span -> has_cds True (predicted),
    # but has_coding_span False (no ORF span -> no UTR partition)
    f3 = by_name["FLT0003"]
    assert f3["productivity"] == "NGO"
    assert cds.has_cds(f3) is True
    assert cds.has_coding_span(f3) is False
    # PRO/PTC with real ORF span -> both True
    assert cds.has_coding_span(f1) is True


def test_alt_end_cds_analysis_produces_nonempty(tmp_path):
    """Previously silently empty (name-regex found no productivity).
    Now per-iso rows are emitted and the SJC sibling pair (PRO/PTC on the
    same chain) yields a productivity-transition pair row."""
    isos = cds.parse_bed12_cds(FLAIRBED)
    per_iso, per_pair = cds.analyse_one_label("FLAIR", isos)
    assert not per_iso.empty
    assert len(per_iso) == 5
    assert set(per_iso["productivity"]) >= {"PRO", "PTC", "NGO", "NST"}
    # FLG0001 has FLT0001(PRO) & FLT0002(PTC) on the same junction chain,
    # both with a coding span -> at least one sibling pair with a
    # PRO->PTC productivity transition.
    assert not per_pair.empty
    assert (per_pair["productivity_transition"] == "PRO->PTC").any()


# ── cds_coverage_recovery now emits data instead of skipping ─────────────────

def test_cds_coverage_recovery_nonempty():
    isos = cds.parse_bed12_cds(FLAIRBED)
    # CONSERVATIVE gate (hdheath 2026-05-29): only isoforms with a real
    # coding span (thickStart != thickEnd -> PRO/PTC) count as "CDS
    # predicted".  NGO (FLT0003) and NST (FLT0005) have thick == span and
    # are NO-CDS, so 3 of the 5 isoforms qualify (PRO,PTC on FLG0001; PRO on
    # FLG0002).
    assert sum(1 for i in isos if ccr.has_cds(i)) == 3
    by_gene = ccr.gene_productivity_classes(isos)
    assert set(by_gene) == {"FLG0001", "FLG0002"}
    # FLG0001 keeps PRO + PTC (both real coding spans); NGO is dropped.
    assert by_gene["FLG0001"]["PRO"] == 1
    assert by_gene["FLG0001"]["PTC"] == 1
    assert by_gene["FLG0001"]["NGO"] == 0
    # FLG0002 keeps PRO; NST is dropped (no coding span).
    assert by_gene["FLG0002"]["PRO"] == 1
    assert by_gene["FLG0002"]["NST"] == 0


def test_cds_coverage_recovery_gate_is_conservative_but_altend_keeps_all_four():
    """The two analyses use DIFFERENT CDS gates on purpose:

      * cds_coverage_recovery.has_cds == has_coding_span -> NGO/NST are
        NO-CDS (dropped from the coverage-recovery diff);
      * alt_end_cds_analysis.has_cds == productivity-present -> NGO/NST
        still count, so productivity-transition classification sees all
        four classes.
    """
    isos = cds.parse_bed12_cds(FLAIRBED)
    by_name = {i["name"]: i for i in isos}
    ngo = by_name["FLT0003"]  # NGO, thick == span
    nst = by_name["FLT0005"]  # NST, thick == span

    # coverage-recovery's gate (conservative): NGO/NST -> NO-CDS
    assert ccr.has_cds(ngo) is False
    assert ccr.has_cds(nst) is False

    # alt_end_cds's own has_cds (productivity present): NGO/NST -> still CDS
    assert cds.has_cds(ngo) is True
    assert cds.has_cds(nst) is True

    # alt_end_cds still classifies all 4 productivity classes across isoforms
    per_iso, _ = cds.analyse_one_label("FLAIR", isos)
    assert set(per_iso["productivity"]) == {"PRO", "PTC", "NGO", "NST"}
