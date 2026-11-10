#!/usr/bin/env python3
"""Byte-identical equivalence: building TED's annotated end/intron sets from the
region-subset GtfData (the in-memory, per-partition path) must equal the
whole-GTF file parse restricted to that chrom.

PERFORMANCE FIX (claim 1): _build_ted_scorer runs once per forked partition
worker; the legacy parse_gtf_ends / parse_gtf_introns re-open + re-scan the
ENTIRE annotation GTF in every worker. The new helpers
ted.annotated_ends_from_gtf_data / annotated_introns_from_gtf_data derive the
same sets from the region-subset GtfData each worker already holds. This module
proves the in-memory build is byte-identical to the file parse so the fix is a
pure redundancy removal, not a behavior change.

The reference path is exactly what the host runs:
  flair.gtf_io.gtf_data_parser(...)  (== flair_transcriptome.py:2070)
  GtfData.subset_for_region(chrom, start, end)  (== partition_runner.py:161)
  ted.annotated_ends_from_gtf_data(subset, chrom)  (== _build_ted_scorer)

Coordinate note: parse_gtf_ends reads the *transcript* feature line as
start=int(cols[3])-1, end=int(cols[4]); GtfData._parse_coordinates uses the
SAME convention, so GtfTranscript.start/.end already match with no ±1 fixup.
The fixture deliberately uses the GENCODE `transcript` feature label (which is
what parse_gtf_ends keys on, cols[2]=='transcript', and which the production
annotation GTF uses) so both parsers see the identical transcript set.

Loads ted.py standalone (mirrors test_ted_two_pass.py); skips if the clone is
absent. Imports flair.gtf_io for the reference parser (the installed flair pkg).
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_CLONE = Path("/private/home/hdheath/projects/flair4_ted/flair-ted-52826")
_TED_PY = _CLONE / "src/flair/ted.py"

pytestmark = pytest.mark.skipif(
    not _TED_PY.exists(), reason=f"TED clone not found: {_TED_PY}")


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ted = _load(_TED_PY, "ted_under_test_gtf_data_equiv")

# Reference parser from the installed flair package (same as the host).
gtf_io = pytest.importorskip(
    "flair.gtf_io", reason="installed flair.gtf_io needed for reference parse")
gtf_data_parser = gtf_io.gtf_data_parser
GtfAttrsSet = gtf_io.GtfAttrsSet
TRANSCRIPT_EXON_FEATURES = gtf_io.TRANSCRIPT_EXON_FEATURES


# A couple of genes, >=2 transcripts each, a '-' strand gene, two chroms (so a
# subset can be exercised). GENCODE-style `transcript` feature lines (1-based,
# end-inclusive GTF) + their exons.
_GTF = (
    # ---- chr1, gene g1 (+ strand), two transcripts -------------------------
    'chr1\ttest\ttranscript\t101\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "g1.t1";\n'
    'chr1\ttest\texon\t101\t300\t.\t+\t.\tgene_id "g1"; transcript_id "g1.t1";\n'
    'chr1\ttest\texon\t801\t1000\t.\t+\t.\tgene_id "g1"; transcript_id "g1.t1";\n'
    'chr1\ttest\texon\t1601\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "g1.t1";\n'
    'chr1\ttest\ttranscript\t501\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "g1.t2";\n'
    'chr1\ttest\texon\t501\t1000\t.\t+\t.\tgene_id "g1"; transcript_id "g1.t2";\n'
    'chr1\ttest\texon\t1601\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "g1.t2";\n'
    # ---- chr1, gene g2 (- strand), two transcripts -------------------------
    'chr1\ttest\ttranscript\t5001\t7000\t.\t-\t.\tgene_id "g2"; transcript_id "g2.t1";\n'
    'chr1\ttest\texon\t5001\t5500\t.\t-\t.\tgene_id "g2"; transcript_id "g2.t1";\n'
    'chr1\ttest\texon\t6001\t6500\t.\t-\t.\tgene_id "g2"; transcript_id "g2.t1";\n'
    'chr1\ttest\texon\t6801\t7000\t.\t-\t.\tgene_id "g2"; transcript_id "g2.t1";\n'
    'chr1\ttest\ttranscript\t5001\t6500\t.\t-\t.\tgene_id "g2"; transcript_id "g2.t2";\n'
    'chr1\ttest\texon\t5001\t5500\t.\t-\t.\tgene_id "g2"; transcript_id "g2.t2";\n'
    'chr1\ttest\texon\t6001\t6500\t.\t-\t.\tgene_id "g2"; transcript_id "g2.t2";\n'
    # ---- chr2, gene g3 (+ strand), single + spliced ------------------------
    'chr2\ttest\ttranscript\t201\t1500\t.\t+\t.\tgene_id "g3"; transcript_id "g3.t1";\n'
    'chr2\ttest\texon\t201\t600\t.\t+\t.\tgene_id "g3"; transcript_id "g3.t1";\n'
    'chr2\ttest\texon\t1101\t1500\t.\t+\t.\tgene_id "g3"; transcript_id "g3.t1";\n'
    'chr2\ttest\ttranscript\t301\t900\t.\t-\t.\tgene_id "g3"; transcript_id "g3.t2";\n'   # opposite strand same locus
    'chr2\ttest\texon\t301\t900\t.\t-\t.\tgene_id "g3"; transcript_id "g3.t2";\n'         # single-exon, no intron
)


def _write_gtf(tmp_path):
    p = tmp_path / "mini.gtf"
    p.write_text(_GTF)
    return str(p)


def _parse(gtf_file):
    return gtf_data_parser(gtf_file, attrs=GtfAttrsSet.FLAIR,
                           include_features=TRANSCRIPT_EXON_FEATURES)


# ════════════════════════════════════════════════════════════════════════
# ENDS — annotated_ends_from_gtf_data == parse_gtf_ends restricted to chrom
# ════════════════════════════════════════════════════════════════════════
@pytest.mark.parametrize("chrom", ["chr1", "chr2"])
def test_ends_equivalence_whole_chrom(tmp_path, chrom):
    """Full-chrom subset: in-memory ends == file parse for that chrom, exactly
    (same tss/tts sets, same sort order)."""
    gtf_file = _write_gtf(tmp_path)

    # Reference: whole-GTF file parse, restricted to chrom.
    file_ends = ted.parse_gtf_ends(gtf_file)
    ref = file_ends.get(chrom, {"tss": [], "tts": []})

    # New path: parse -> subset_for_region(whole chrom) -> in-memory ends.
    gtf_data = _parse(gtf_file)
    sub = gtf_data.subset_for_region(chrom, 0, 10_000_000)
    got = ted.annotated_ends_from_gtf_data(sub, chrom).get(
        chrom, {"tss": [], "tts": []})

    assert got["tss"] == ref["tss"], f"{chrom} TSS mismatch"
    assert got["tts"] == ref["tts"], f"{chrom} TTS mismatch"
    # sanity: the '-' gene g2 on chr1 puts a TSS at its (0-based) end and a TTS
    # at its start, so the strand-aware assignment is genuinely exercised.
    if chrom == "chr1":
        assert 7000 in got["tss"] and 5000 in got["tts"]


def test_ends_strand_aware_minus_gene(tmp_path):
    """Explicit '-' strand check mirroring parse_gtf_ends lines 672-677:
    tss=end, tts=start for '-'; tss=start, tts=end for '+'."""
    gtf_file = _write_gtf(tmp_path)
    gtf_data = _parse(gtf_file)
    sub = gtf_data.subset_for_region("chr1", 0, 10_000_000)
    ends = ted.annotated_ends_from_gtf_data(sub, "chr1")["chr1"]
    # g2.t1: transcript 5001..7000 '-' => start0=5000 (TTS), end=7000 (TSS)
    # g2.t2: transcript 5001..6500 '-' => start0=5000 (TTS), end=6500 (TSS)
    assert 7000 in ends["tss"] and 6500 in ends["tss"]
    assert 5000 in ends["tts"]
    # g1 (+): transcript 101..2000 => start0=100 (TSS), end=2000 (TTS)
    #         transcript 501..2000 => start0=500 (TSS), end=2000 (TTS)
    assert 100 in ends["tss"] and 500 in ends["tss"]
    assert 2000 in ends["tts"]


def test_ends_subregion_only_contains_subregion_transcripts(tmp_path):
    """A narrower region-subset yields ONLY that sub-region's transcripts' ends.
    Restrict chr1 to the g1 locus (101..2000): g2 (5001..7000) must be absent."""
    gtf_file = _write_gtf(tmp_path)
    gtf_data = _parse(gtf_file)

    sub = gtf_data.subset_for_region("chr1", 100, 2100)  # g1 only
    ends = ted.annotated_ends_from_gtf_data(sub, "chr1")["chr1"]
    # g1 ends present:
    assert 100 in ends["tss"] and 500 in ends["tss"]
    assert 2000 in ends["tts"]
    # g2 ends ABSENT (its transcripts do not overlap [100, 2100)):
    assert 7000 not in ends["tss"] and 6500 not in ends["tss"]
    assert 5000 not in ends["tts"]
    # and the in-memory subset must equal building from the subset's
    # transcripts directly (no leakage from outside the window):
    spans = [(t.chrom, t.start, t.end, t.strand) for t in sub.transcripts]
    assert ted.build_annotated_ends(spans)["chr1"]["tss"] == ends["tss"]
    assert ted.build_annotated_ends(spans)["chr1"]["tts"] == ends["tts"]


def test_ends_chrom_filter_excludes_other_chroms(tmp_path):
    """A subset is single-chrom in the partition model, but the helper also
    hard-filters on chrom: asking for chr2 ends from a chr1 subset is empty."""
    gtf_file = _write_gtf(tmp_path)
    gtf_data = _parse(gtf_file)
    sub = gtf_data.subset_for_region("chr1", 0, 10_000_000)
    assert ted.annotated_ends_from_gtf_data(sub, "chr2") == {}


# ════════════════════════════════════════════════════════════════════════
# INTRONS — annotated_introns_from_gtf_data == parse_gtf_introns / chrom
# ════════════════════════════════════════════════════════════════════════
@pytest.mark.parametrize("chrom", ["chr1", "chr2"])
def test_introns_equivalence_whole_chrom(tmp_path, chrom):
    """In-memory intron set == whole-GTF parse_gtf_introns restricted to the
    (chrom, strand) keys for that chrom, exactly (TOL-binned donor/acceptor)."""
    gtf_file = _write_gtf(tmp_path)

    file_introns = ted.parse_gtf_introns(gtf_file)
    ref = {k: v for k, v in file_introns.items() if k[0] == chrom}

    gtf_data = _parse(gtf_file)
    sub = gtf_data.subset_for_region(chrom, 0, 10_000_000)
    got = ted.annotated_introns_from_gtf_data(sub, chrom)

    assert got == ref, f"{chrom} intron set mismatch"


def test_introns_values_match_known_junctions(tmp_path):
    """Spot-check concrete TOL-binned introns for g1.t1 on chr1 (+).
    exons 101..300, 801..1000, 1601..2000 (1-based, inclusive) => 0-based ends
    300/1000 donors, next-exon 0-based starts 800/1600 acceptors."""
    gtf_file = _write_gtf(tmp_path)
    gtf_data = _parse(gtf_file)
    sub = gtf_data.subset_for_region("chr1", 0, 10_000_000)
    got = ted.annotated_introns_from_gtf_data(sub, "chr1")
    TOL = ted._JUNC_TOL
    plus = got[("chr1", "+")]
    assert (300 // TOL, 800 // TOL) in plus      # g1.t1 intron 1
    assert (1000 // TOL, 1600 // TOL) in plus     # g1.t1 / g1.t2 shared intron 2
    # '-' gene g2 introns present on the '-' key:
    minus = got[("chr1", "-")]
    assert (5500 // TOL, 6000 // TOL) in minus
    assert (6500 // TOL, 6800 // TOL) in minus


def test_introns_subregion_excludes_outside_transcripts(tmp_path):
    """Narrow region-subset -> only its transcripts' introns. Restrict chr1 to
    the g1 locus: the g2 (- strand) intron key must be absent entirely."""
    gtf_file = _write_gtf(tmp_path)
    gtf_data = _parse(gtf_file)
    sub = gtf_data.subset_for_region("chr1", 100, 2100)  # g1 only
    got = ted.annotated_introns_from_gtf_data(sub, "chr1")
    assert ("chr1", "+") in got
    assert ("chr1", "-") not in got   # g2 excluded by the window


def test_introns_single_exon_transcript_contributes_nothing(tmp_path):
    """g3.t2 on chr2 is single-exon (no consecutive exon pair) -> contributes no
    intron; the only chr2 intron key is the '+' one from g3.t1."""
    gtf_file = _write_gtf(tmp_path)
    gtf_data = _parse(gtf_file)
    sub = gtf_data.subset_for_region("chr2", 0, 10_000_000)
    got = ted.annotated_introns_from_gtf_data(sub, "chr2")
    assert set(got.keys()) == {("chr2", "+")}
    TOL = ted._JUNC_TOL
    assert (600 // TOL, 1100 // TOL) in got[("chr2", "+")]   # g3.t1 intron
