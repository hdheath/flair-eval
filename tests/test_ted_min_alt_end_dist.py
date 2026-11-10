#!/usr/bin/env python3
"""Tests for --ted_min_alt_end_dist (opt-in near-duplicate end suppression).

The flag adds a minimum strand-aware distance rule to FLAIR-TED's multi-end
emission (the only place >1 end per junction chain is emitted). When
--ted_min_alt_end_dist D > 0 and --max_ends > 1, alternative ends are kept
GREEDILY in descending-support order: the highest-support end at a peak is
kept, and any candidate within D bp of an already-kept end on BOTH axes (TSS
AND TTS) is suppressed (its reads merged into the nearest kept end, conserving
support). D=0 (default / absent) is a strict no-op => BYTE-IDENTICAL output.

Two layers:
  1. Unit (no subprocess) on flair_transcriptome.apply_min_alt_end_dist:
       * 2bp apart, D=50 -> only the higher-support end emitted
       * 200bp apart, D=50 -> both emitted
       * D=0 -> both always (no-op, no read merge)
       * strand-aware (- strand uses end for TSS / start for TTS)
       * support conserved (skipped end's reads merged into nearest kept)
  2. End-to-end byte-identity: on the tiny fixture, the full `flair
     transcriptome --ted` output is byte-identical between (flag ABSENT) and
     (--ted_min_alt_end_dist 0). This is the default-off byte-identical proof.
"""
from __future__ import annotations

import hashlib
import os
import subprocess
import sys
import tempfile
from collections import namedtuple
from pathlib import Path

import pytest

# ─── Locations (skip cleanly if anything is missing) ─────────────────────
ENV_PREFIX = Path("/private/home/hdheath/miniforge3/envs/flair-ted-52826")
FLAIR_BIN = ENV_PREFIX / "bin" / "flair"
CLONE = Path("/private/home/hdheath/projects/flair4_ted/flair-ted-52826")
DATA = Path(__file__).resolve().parent / "data"
MODELS = Path("/private/home/hdheath/projects/artifact_classifier/model_weights")
TSS_MODEL = MODELS / "tss" / "xgboost.pkl"
TTS_MODEL = MODELS / "tts" / "xgboost.pkl"

_CLONE_SRC = CLONE / "src"
if str(_CLONE_SRC) not in sys.path:
    sys.path.insert(0, str(_CLONE_SRC))

try:
    from flair import flair_transcriptome as _ft
    _FT_IMPORTABLE = True
except Exception:  # pragma: no cover - clone/deps unavailable
    _FT_IMPORTABLE = False


# ─── Unit layer: apply_min_alt_end_dist ──────────────────────────────────
class _StubIso:
    """Lightweight Isoform stand-in: only start/end/strand/reads/num_reads/
    genomic_length are read by apply_min_alt_end_dist."""
    def __init__(self, start, end, strand, nreads):
        self.start = start
        self.end = end
        self.strand = strand
        self.reads = [object() for _ in range(nreads)]
        self.juncs = ((1200, 1500),)

    @property
    def num_reads(self):
        return len(self.reads)

    @property
    def genomic_length(self):
        return self.end - self.start


@pytest.mark.skipif(not _FT_IMPORTABLE, reason="flair_transcriptome not importable")
def test_two_ends_2bp_apart_D50_keeps_higher_support_only():
    """2bp apart, D=50 -> only the higher-support end emitted; the skipped
    end's reads are merged into it (support conserved)."""
    hi = _StubIso(1000, 5000, "+", 10)
    lo = _StubIso(1002, 5000, "+", 3)   # 2bp upstream TSS, identical TTS
    out = _ft.apply_min_alt_end_dist([hi, lo], 50)
    assert len(out) == 1, "near-dup (2bp) end not suppressed at D=50"
    assert out[0] is hi, "kept end is not the higher-support one"
    assert out[0].num_reads == 13, (
        f"reads not merged (support not conserved): {out[0].num_reads} != 13")


@pytest.mark.skipif(not _FT_IMPORTABLE, reason="flair_transcriptome not importable")
def test_two_ends_200bp_apart_D50_keeps_both():
    """200bp apart, D=50 -> both emitted (distinct on the TSS axis)."""
    hi = _StubIso(1000, 5000, "+", 10)
    lo = _StubIso(1200, 5000, "+", 3)   # 200bp upstream TSS
    out = _ft.apply_min_alt_end_dist([hi, lo], 50)
    assert len(out) == 2, "distinct (200bp) alt end was wrongly suppressed"
    assert {o.start for o in out} == {1000, 1200}
    # No merge happened: support split stays 10 / 3.
    assert sorted(o.num_reads for o in out) == [3, 10]


@pytest.mark.skipif(not _FT_IMPORTABLE, reason="flair_transcriptome not importable")
def test_D0_is_noop_both_always_kept_no_merge():
    """D=0 -> both ends always emitted, input returned unchanged, no read
    merge (the byte-identical default contract at the unit level)."""
    hi = _StubIso(1000, 5000, "+", 10)
    lo = _StubIso(1002, 5000, "+", 3)   # would be a near-dup at D>0
    inp = [hi, lo]
    out = _ft.apply_min_alt_end_dist(inp, 0)
    assert out is inp, "D=0 should return the exact input list (no-op)"
    assert len(out) == 2
    assert hi.num_reads == 10 and lo.num_reads == 3, (
        "D=0 must not merge reads (no-op)")


@pytest.mark.skipif(not _FT_IMPORTABLE, reason="flair_transcriptome not importable")
def test_negative_D_is_noop():
    """Defensive: D<0 behaves like OFF."""
    hi = _StubIso(1000, 5000, "+", 10)
    lo = _StubIso(1002, 5000, "+", 3)
    out = _ft.apply_min_alt_end_dist([hi, lo], -1)
    assert len(out) == 2


@pytest.mark.skipif(not _FT_IMPORTABLE, reason="flair_transcriptome not importable")
def test_single_end_is_noop():
    """One isoform -> nothing to suppress, returned unchanged even at D>0."""
    only = _StubIso(1000, 5000, "+", 10)
    out = _ft.apply_min_alt_end_dist([only], 100)
    assert out == [only]


@pytest.mark.skipif(not _FT_IMPORTABLE, reason="flair_transcriptome not importable")
def test_distinct_on_tts_axis_only_is_kept():
    """OR-on-axis: same TSS but a distinct (>=D) TTS -> kept (introduces a new
    distinct 3' end / peak pair)."""
    hi = _StubIso(1000, 5000, "+", 10)
    lo = _StubIso(1000, 5300, "+", 3)   # identical TSS, 300bp distinct TTS
    out = _ft.apply_min_alt_end_dist([hi, lo], 50)
    assert len(out) == 2, "distinct-TTS alt end wrongly suppressed (OR-on-axis)"


@pytest.mark.skipif(not _FT_IMPORTABLE, reason="flair_transcriptome not importable")
def test_minus_strand_is_strand_aware():
    """- strand: TSS == genomic end, TTS == genomic start. Two ends whose
    genomic END differs by 2bp (TSS near-dup) and genomic START is identical
    (TTS identical) are near-dups -> suppressed at D=50; mirror of + strand."""
    hi = _StubIso(1000, 5000, "-", 10)
    lo = _StubIso(1000, 5002, "-", 3)   # 2bp on the - strand TSS (genomic end)
    out = _ft.apply_min_alt_end_dist([hi, lo], 50)
    assert len(out) == 1, "minus-strand near-dup not suppressed (not strand-aware?)"
    assert out[0].num_reads == 13
    # A 200bp - strand TSS difference is distinct -> kept.
    hi2 = _StubIso(1000, 5000, "-", 10)
    lo2 = _StubIso(1000, 5200, "-", 3)
    out2 = _ft.apply_min_alt_end_dist([hi2, lo2], 50)
    assert len(out2) == 2, "distinct minus-strand alt end wrongly suppressed"


@pytest.mark.skipif(not _FT_IMPORTABLE, reason="flair_transcriptome not importable")
def test_filter_ends_allow_multiple_default_byte_identical():
    """filter_ends_allow_multiple with min_alt_end_dist=0 (default) must return
    exactly the legacy [:max_ends] slice (object identity preserved)."""
    a = _StubIso(1000, 5000, "+", 10)
    b = _StubIso(1002, 5000, "+", 8)
    c = _StubIso(1004, 5000, "+", 6)
    legacy = _ft.filter_ends_allow_multiple([a, b, c], sjc_support=1, max_ends=3)
    assert legacy == [a, b, c]
    explicit0 = _ft.filter_ends_allow_multiple(
        [a, b, c], sjc_support=1, max_ends=3, min_alt_end_dist=0)
    assert explicit0 == [a, b, c]
    # And with the flag ON the near-dups collapse to one.
    a2 = _StubIso(1000, 5000, "+", 10)
    b2 = _StubIso(1002, 5000, "+", 8)
    c2 = _StubIso(1004, 5000, "+", 6)
    on = _ft.filter_ends_allow_multiple(
        [a2, b2, c2], sjc_support=1, max_ends=3, min_alt_end_dist=50)
    assert len(on) == 1 and on[0].num_reads == 24


# ─── End-to-end byte-identity: flag absent == --ted_min_alt_end_dist 0 ────
_missing = [str(p) for p in (FLAIR_BIN, CLONE, DATA, TSS_MODEL, TTS_MODEL)
            if not p.exists()]


def _run_ted(outdir: Path, extra_args=()):
    env = dict(os.environ)
    env["PATH"] = f"{ENV_PREFIX / 'bin'}:" + env.get("PATH", "")
    out = outdir / "tiny"
    cmd = [
        str(FLAIR_BIN), "transcriptome",
        "-b", str(DATA / "tiny_reads.bam"),
        "-g", str(DATA / "tiny_genome.fa"),
        "-f", str(DATA / "tiny_annotation.gtf"),
        "--junction_tab", str(DATA / "tiny_junctions.tab"),
        "-t", "1", "-o", str(out), "--sample_name", "ted_test",
        "--ted", "--ted_legacy_reality",
        "--tss_model", str(TSS_MODEL), "--tts_model", str(TTS_MODEL),
        "--ted_log", "--keep_intermediate",
        *extra_args,
    ]
    res = subprocess.run(cmd, env=env, cwd=str(outdir),
                         capture_output=True, text=True)
    assert res.returncode == 0, (
        f"flair transcriptome --ted exited {res.returncode}\n"
        f"STDERR tail:\n{res.stderr[-2000:]}")
    return out


def _md5(path: Path):
    return hashlib.md5(path.read_bytes()).hexdigest()


@pytest.mark.skipif(bool(_missing),
                    reason=f"TED graft prerequisites missing: {_missing}")
def test_end_to_end_flag_absent_equals_D0_byte_identical():
    """The full --ted output is BYTE-IDENTICAL between (flag ABSENT) and
    (--ted_min_alt_end_dist 0). This is the default-off byte-identical proof
    over the real pipeline (md5 of .isoforms.bed / .gtf / .ted_log.tsv /
    .firstpass.unfiltered.bed). Also runs --max_ends 3 to exercise the
    multi-end branch where the flag actually engages — at D=0 it must STILL be
    byte-identical to flag-absent at the same --max_ends."""
    suffixes = [".isoforms.bed", ".isoforms.gtf", ".ted_log.tsv",
                ".firstpass.unfiltered.bed", ".firstpass.bed"]
    for me in ("1", "3"):
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            (tdp / "absent").mkdir(); (tdp / "d0").mkdir()
            absent = _run_ted(tdp / "absent", extra_args=("--max_ends", me))
            d0 = _run_ted(tdp / "d0",
                          extra_args=("--max_ends", me, "--ted_min_alt_end_dist", "0"))
            for sfx in suffixes:
                pa = Path(str(absent) + sfx)
                pd = Path(str(d0) + sfx)
                if not pa.exists():
                    continue
                assert _md5(pa) == _md5(pd), (
                    f"[max_ends={me}] {sfx} differs between flag-absent and "
                    f"--ted_min_alt_end_dist 0 — default is NOT byte-identical")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v", "--tb=short"]))
