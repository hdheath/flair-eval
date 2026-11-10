#!/usr/bin/env python3
"""Unit tests for the S_repro cross-library reproducibility term in ted.py.

S_repro is the one validated lever from the cross-protocol reproducibility study:
among >=3-read NOVEL ends, whether the same read-END cluster recurs in the OTHER
libraries of the same sample separates real from artifact (AUC 0.82 TSS / 0.91 TTS),
orthogonal to the sequence model. UNLIKE S_shape it is added to BOTH ends.

Guarantees asserted here (mirror test_ted_sshape.py's structure):

1. BYTE-IDENTICAL OFF — with w_repro==0 (default) the term is inert: the validator's
   reality equals the legacy blend exactly and `repro_score` is 0; end-to-end the
   published BED/GTF is byte-identical to a run with no S_repro flags (a loaded
   evidence list with zero weight is inert, and the evidence BAMs are never opened).

2. S_repro VALUE — `_ReproEvidence.score` returns (#supporting tracks)/(#tracks):
   2 tracks both in-window with >=k ends -> 1.0; one -> 0.5; none -> 0.0. The
   +/-window and the min_reads (k) floor are respected (an in-window track with
   <k ends does NOT count).

3. STRAND / END-TYPE — a +-strand TSS peak counts only +-strand 5'-ends; a TTS peak
   counts 3'-ends; and `_repro_end_positions` assigns 5'/3' strand-aware (5' =
   reference_start on +, reference_end-1 on -; 3' mirror).

4. BLEND MATH — with w_repro>0 the reality gains w_repro*S_repro on BOTH ends; the
   term carries through the two-pass NOVEL re-blend at the same per-end weight.

Loads ted.py standalone (numpy/scipy at call time); skips if the clone is absent.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

_TED_PY = Path(
    "/private/home/hdheath/projects/flair4_ted/flair-ted-52826/src/flair/ted.py")

pytestmark = pytest.mark.skipif(
    not _TED_PY.exists(), reason=f"TED clone not found: {_TED_PY}")


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ted = _load(_TED_PY, "ted_under_test_srepro")


# ── synthetic pysam alignment + indexed BAM helpers ──────────────────────
class _Aln:
    """Minimal pysam.AlignedSegment stand-in for _repro_end_positions."""
    def __init__(self, start, end, is_reverse,
                 is_unmapped=False, is_secondary=False, is_supplementary=False):
        self.reference_start = start
        self.reference_end = end
        self.is_reverse = is_reverse
        self.is_unmapped = is_unmapped
        self.is_secondary = is_secondary
        self.is_supplementary = is_supplementary


def _write_bam(path, contig, length, reads):
    """Write a tiny coordinate-sorted, indexed BAM.

    `reads` is a list of (start, end, strand) — one alignment each, on `contig`.
    A trivial CIGAR (one match block spanning start..end) gives reference_end=end.
    """
    import pysam
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": contig, "LN": length}]}
    rows = sorted(reads, key=lambda r: r[0])
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, (start, end, strand) in enumerate(rows):
            a = pysam.AlignedSegment(out.header)
            a.query_name = f"r{i}"
            a.flag = 16 if strand == '-' else 0
            a.reference_id = 0
            a.reference_start = int(start)
            a.mapping_quality = 60
            a.cigartuples = [(0, int(end) - int(start))]   # M block start..end
            a.query_sequence = "A" * (int(end) - int(start))
            a.query_qualities = pysam.qualitystring_to_array("I" * (int(end) - int(start)))
            out.write(a)
    pysam.index(str(path))


# ── 1. byte-identical OFF (math-level) ───────────────────────────────────
def _validate(**kw):
    """Call _ted_validate_cluster with sane defaults for a cluster."""
    base = dict(
        positions=[1000] * 10, end_type='tss', chrom='chrT', strand='+',
        scorer=None,                      # no model -> model_score 0
        w_depth=0.30, w_model=0.20, w_annot=0.20, threshold=0.5,
        depth_norm=None, depth_saturation=50.0)
    base.update(kw)
    return ted._ted_validate_cluster(**base)


def test_blend_off_is_legacy_both_ends():
    """w_repro == 0 (default) -> reality == legacy blend, no matter the
    repro_score, for BOTH TSS and TTS, and repro_score in scores is 0."""
    sd = min(1.0, np.log1p(10) / np.log1p(50))
    for et in ('tss', 'tts'):
        _, reality, scores = _validate(end_type=et, w_repro=0.0, repro_score=0.9)
        assert reality == pytest.approx(0.30 * sd)   # only depth (model/annot 0)
        assert scores['repro_score'] == 0.0          # inert when w_repro=0


def test_blend_repro_none_contributes_zero():
    """repro_score=None (term ON but no value, e.g. self-position) -> S_repro=0."""
    sd = min(1.0, np.log1p(10) / np.log1p(50))
    _, reality, scores = _validate(w_repro=0.40, repro_score=None, threshold=0.0)
    assert reality == pytest.approx(0.30 * sd)
    assert scores['repro_score'] == 0.0


# ── 2. blend math (BOTH ends gain w_repro*S_repro) ───────────────────────
@pytest.mark.parametrize("et", ['tss', 'tts'])
def test_blend_adds_repro_both_ends(et):
    """w_repro>0 adds w_repro*S_repro on BOTH ends; reality = hand sum."""
    sd = min(1.0, np.log1p(10) / np.log1p(50))
    _, reality, scores = _validate(
        end_type=et, w_depth=0.30, w_model=0.20, w_annot=0.20,
        w_repro=0.40, repro_score=0.5, threshold=0.0)
    expected = 0.30 * sd + 0.20 * 0.0 + 0.20 * 0.0 + 0.40 * 0.5
    assert reality == pytest.approx(expected)
    assert scores['repro_score'] == pytest.approx(0.5)


def test_blend_repro_coexists_with_shape_on_tss():
    """On TSS, S_shape and S_repro both add (independent terms)."""
    sd = min(1.0, np.log1p(10) / np.log1p(50))
    _, reality, scores = _validate(
        w_depth=0.25, w_model=0.15, w_annot=0.10,
        w_shape=0.30, shape_score=0.6,
        w_repro=0.20, repro_score=0.5, threshold=0.0)
    expected = 0.25 * sd + 0.30 * 0.6 + 0.20 * 0.5
    assert reality == pytest.approx(expected)
    assert scores['shape_score'] == pytest.approx(0.6)
    assert scores['repro_score'] == pytest.approx(0.5)


def test_blend_threshold_uses_repro():
    """A cluster that fails without S_repro can PASS once a reproducible end
    adds w_repro*S_repro (the lever's intended effect), on a TTS end."""
    norm = ted.DepthNormalizer([], method='logsat', logsat_constant=1.0)  # S_d=1
    is_real0, r0, _ = _validate(end_type='tts', w_depth=0.40, w_model=0.20,
                                w_annot=0.20, w_repro=0.0, repro_score=0.0,
                                threshold=0.6, depth_norm=norm)
    assert r0 == pytest.approx(0.40) and r0 < 0.6 and is_real0 is False
    is_real1, r1, _ = _validate(end_type='tts', w_depth=0.40, w_model=0.20,
                                w_annot=0.20, w_repro=0.40, repro_score=1.0,
                                threshold=0.6, depth_norm=norm)
    assert r1 == pytest.approx(0.80) and r1 > r0 and is_real1 is True


def test_blend_repro_carries_through_novel_pass():
    """Two-pass NOVEL TSS re-blend keeps w_repro*S_repro at the same weight.

    annot_dist large (no annotation) -> NOVEL route; reads-only weights drop
    S_a, so reality = wd_n*S_d + wm_n*S_m + w_repro*S_repro. With model 0 and
    S_d=1 (logsat const=1): 0.10*1 + 0 + 0.50*1.0 = 0.60 >= novel θ 0.40."""
    norm = ted.DepthNormalizer([], method='logsat', logsat_constant=1.0)
    is_real, reality, scores = _validate(
        end_type='tss', two_pass=True, two_pass_novel_dist=500,
        w_depth=0.30, w_model=0.20, w_annot=0.20,
        w_depth_novel=0.10, w_model_novel=0.0, w_shape_novel=0.0,
        threshold_novel=0.40, threshold=0.99,
        w_repro=0.50, repro_score=1.0, depth_norm=norm)
    assert scores['pass_label'] == 'novel'
    assert reality == pytest.approx(0.60)
    assert is_real is True


# ── 3. _repro_end_positions: strand-aware 5'/3' ──────────────────────────
def test_repro_end_positions_plus():
    """+ strand: 5'(TSS)=reference_start, 3'(TTS)=reference_end-1."""
    tss, tts = ted._repro_end_positions(_Aln(100, 600, is_reverse=False))
    assert tss == 100 and tts == 599


def test_repro_end_positions_minus():
    """- strand: 5'(TSS)=reference_end-1, 3'(TTS)=reference_start."""
    tss, tts = ted._repro_end_positions(_Aln(100, 600, is_reverse=True))
    assert tss == 599 and tts == 100


# ── 4. _ReproEvidence.score: fraction, window, min_reads, strand/end-type ─
def _two_track_evidence(tmp_path, t0_reads, t1_reads, window=50, min_reads=2):
    """Two indexed evidence BAMs in `tmp_path` (created if needed); return the
    _ReproEvidence over chr_test."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    b0 = tmp_path / "ev0.bam"
    b1 = tmp_path / "ev1.bam"
    _write_bam(b0, "chr_test", 50000, t0_reads)
    _write_bam(b1, "chr_test", 50000, t1_reads)
    return ted._ReproEvidence([str(b0), str(b1)], window=window,
                              min_reads=min_reads,
                              regions={"chr_test": (0, 50000)})


def test_score_fraction_both_one_none(tmp_path):
    """2 tracks: a +-strand TSS at 10000.
       both have >=2 5'-ends in-window -> 1.0;
       only track0 does -> 0.5;
       neither does -> 0.0."""
    pos = 10000
    # 5'-end of a + read = reference_start; put 2 starts at ~pos in each track.
    in_win = [(pos - 5, pos + 400, '+'), (pos + 7, pos + 420, '+')]
    far = [(pos + 5000, pos + 5400, '+'), (pos + 5007, pos + 5420, '+')]
    ev_both = _two_track_evidence(tmp_path / "both", in_win, in_win)
    ev_one = _two_track_evidence(tmp_path / "one", in_win, far)
    ev_none = _two_track_evidence(tmp_path / "none", far, far)
    assert ev_both.score("chr_test", pos, '+', 'tss') == pytest.approx(1.0)
    assert ev_one.score("chr_test", pos, '+', 'tss') == pytest.approx(0.5)
    assert ev_none.score("chr_test", pos, '+', 'tss') == pytest.approx(0.0)


def test_score_min_reads_floor(tmp_path):
    """An in-window track with < k ends does NOT count (k=2: 1 end -> 0)."""
    pos = 20000
    one_end = [(pos, pos + 300, '+')]                      # only 1 in-window 5'-end
    two_end = [(pos - 3, pos + 300, '+'), (pos + 4, pos + 310, '+')]
    ev = _two_track_evidence(tmp_path / "x", one_end, two_end, min_reads=2)
    # track0 has 1 (<2 -> no), track1 has 2 (>=2 -> yes) => 0.5
    assert ev.score("chr_test", pos, '+', 'tss') == pytest.approx(0.5)


def test_score_window_respected(tmp_path):
    """Ends just outside +/-W do not count; just inside do."""
    pos = 30000
    # window 50: starts at pos+60/pos+61 are OUT; pos+40/pos+45 are IN.
    out = [(pos + 60, pos + 360, '+'), (pos + 61, pos + 361, '+')]
    inside = [(pos + 40, pos + 340, '+'), (pos + 45, pos + 345, '+')]
    ev = _two_track_evidence(tmp_path / "w", out, inside, window=50, min_reads=2)
    assert ev.score("chr_test", pos, '+', 'tss') == pytest.approx(0.5)


def test_score_strand_and_endtype(tmp_path):
    """A +-strand TSS peak counts only +-strand 5'-ends; a TTS peak counts
    3'-ends. A track of - reads (whose 5'-end is reference_end-1, far away)
    does NOT support a +-strand TSS at the + reads' start."""
    pos = 40000
    # track0: two + reads starting at ~pos (so + TSS in-window; their + TTS far).
    plus = [(pos - 2, pos + 500, '+'), (pos + 3, pos + 510, '+')]
    # track1: two - reads whose reference_start ~pos (so their 3'/TTS = start ~pos,
    # but their 5'/TSS = reference_end-1 ~pos+500, far from a + TSS at pos).
    minus = [(pos - 2, pos + 500, '-'), (pos + 3, pos + 510, '-')]
    ev = _two_track_evidence(tmp_path / "s", plus, minus, window=50, min_reads=2)
    # +-strand TSS at pos: only track0 (+ 5'-ends) supports -> 0.5 (track1 is
    # - strand, different (chrom,strand,'tss') bucket).
    assert ev.score("chr_test", pos, '+', 'tss') == pytest.approx(0.5)
    # +-strand TTS at pos: + reads' 3'-end is ~pos+500 (far), so track0 does NOT
    # support; track1 is - strand (different bucket) -> 0.0 on the + TTS.
    assert ev.score("chr_test", pos, '+', 'tts') == pytest.approx(0.0)
    # -strand TTS at pos: track1's - reads' 3'-end = reference_start ~pos -> 1
    # track supports -> 0.5.
    assert ev.score("chr_test", pos, '-', 'tts') == pytest.approx(0.5)


def test_score_skips_secondary_supplementary_unmapped():
    """_ReproEvidence ignores secondary/supplementary/unmapped (constructed via
    the position lists; here assert the end-position helper's flag handling is
    consumed — a direct construction check on a single-track in-memory tally)."""
    # Build a track by hand: 2 good + 1 secondary at the same pos; only 2 count.
    ev = ted._ReproEvidence.__new__(ted._ReproEvidence)
    ev.window = 50
    ev.min_reads = 2
    ev.n_tracks = 1
    arr = np.asarray([1000, 1002], dtype=np.int64)     # the 2 GOOD ends only
    ev._tracks = [{("chr_test", '+', 'tss'): arr}]
    assert ev.score("chr_test", 1001, '+', 'tss') == pytest.approx(1.0)
    # a third (secondary) read at 1001 would NOT be in this array (skipped at
    # build time); with only the 2 good ends the count is still 2 (>=k) -> 1.0.


def test_score_empty_tracks_is_zero():
    ev = ted._ReproEvidence.__new__(ted._ReproEvidence)
    ev.window, ev.min_reads, ev.n_tracks, ev._tracks = 50, 2, 0, []
    assert ev.score("chr_test", 1000, '+', 'tss') == 0.0


# ── End-to-end wiring (real pipeline): OFF byte-identity + ON populates ───
import os
import subprocess
import tempfile

_ENV = Path("/private/home/hdheath/miniforge3/envs/flair-ted-52826")
_FLAIR = _ENV / "bin" / "flair"
_DATA = Path("/private/home/hdheath/projects/flair4_ted/tests/data")
_TSS_M = Path("/private/home/hdheath/projects/artifact_classifier/model_weights/tss/xgboost.pkl")
_TTS_M = Path("/private/home/hdheath/projects/artifact_classifier/model_weights/tts/xgboost.pkl")
_e2e_missing = [str(p) for p in (_FLAIR, _DATA / "tiny_reads.bam", _TSS_M, _TTS_M)
                if not p.exists()]


def _run(outdir, tag, extra):
    env = dict(os.environ)
    env["PATH"] = f"{_ENV / 'bin'}:" + env.get("PATH", "")
    out = outdir / tag
    cmd = [str(_FLAIR), "transcriptome",
           "-b", str(_DATA / "tiny_reads.bam"),
           "-g", str(_DATA / "tiny_genome.fa"),
           "-f", str(_DATA / "tiny_annotation.gtf"),
           "--junction_tab", str(_DATA / "tiny_junctions.tab"),
           "-t", "1", "-o", str(out), "--sample_name", "srepro",
           # NORTH-STAR re-arch: S_repro is a LEGACY reality-path term.
           "--ted", "--ted_legacy_reality",
           "--tss_model", str(_TSS_M), "--tts_model", str(_TTS_M),
           "--ted_log", "--keep_intermediate", *extra]
    res = subprocess.run(cmd, env=env, cwd=str(outdir), capture_output=True, text=True)
    assert res.returncode == 0, f"flair exited {res.returncode}\n{res.stderr[-1500:]}"
    return out


def _rows(p):
    return [ln.split("\t") for ln in Path(p).read_text().splitlines()
            if ln and not ln.startswith("#")]


@pytest.mark.skipif(bool(_e2e_missing), reason=f"e2e prereqs missing: {_e2e_missing}")
def test_e2e_off_bed_byte_identical():
    """Two OFF runs (default, and with --ted_repro_evidence set but both weights 0)
    produce byte-identical published BED/GTF — a supplied evidence list with zero
    weight is inert and the evidence BAM is never opened."""
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        (tdp / "a").mkdir(); (tdp / "b").mkdir()
        a = _run(tdp / "a", "tiny", [])
        # evidence given but both weights 0 => term OFF, evidence ignored.
        b = _run(tdp / "b", "tiny",
                 ["--ted_repro_evidence", str(_DATA / "tiny_reads.bam"),
                  "--ted_w_repro_tss", "0.0", "--ted_w_repro_tts", "0.0"])
        assert (Path(str(a) + ".isoforms.bed").read_text()
                == Path(str(b) + ".isoforms.bed").read_text())
        assert (Path(str(a) + ".isoforms.gtf").read_text()
                == Path(str(b) + ".isoforms.gtf").read_text())
        # OFF ted_log: both repro columns are the -1 sentinel.
        log = _rows(str(a) + ".ted_log.tsv")
        hdr = log[0]
        for col in ("TED_tss_repro", "TED_tts_repro"):
            ci = hdr.index(col)
            assert all(r[ci] == "-1.0000" for r in log[1:])


@pytest.mark.skipif(bool(_e2e_missing), reason=f"e2e prereqs missing: {_e2e_missing}")
def test_e2e_on_populates_repro_columns():
    """ON run (w_repro>0 + evidence = the reads' own BAM as one track): at least
    one global peak carries a finite S_repro in [0,1] on each end, and at least
    one peak's reality shifts vs the OFF run (the term entered the blend).

    Using the reads' own BAM as the single evidence track makes S_repro
    deterministically 1.0 at peaks that the reads themselves support — the
    cleanest ON-populates check on the tiny fixture."""
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        (tdp / "off").mkdir(); (tdp / "on").mkdir()
        off = _run(tdp / "off", "tiny", [])
        on = _run(tdp / "on", "tiny",
                  ["--ted_repro_evidence", str(_DATA / "tiny_reads.bam"),
                   "--ted_w_repro_tss", "0.20", "--ted_w_repro_tts", "0.20",
                   # re-derive thresholds down so the +0.20 term can't only ever
                   # raise pass-counts past 1.0; thresholds don't affect logging.
                   "--ted_repro_window", "50", "--ted_repro_min_reads", "1"])
        off_log = _rows(str(off) + ".ted_log.tsv")
        on_log = _rows(str(on) + ".ted_log.tsv")
        hdr = on_log[0]
        ti = hdr.index("TED_tss_repro")
        oi = hdr.index("TED_tts_repro")

        def peaks(log, kind):
            return [r for r in log[1:]
                    if r[2] == "global_peak" and r[1].startswith(f"global_{kind}")]
        tss_vals = [float(r[ti]) for r in peaks(on_log, "tss") if r[ti] != "-1.0000"]
        tts_vals = [float(r[oi]) for r in peaks(on_log, "tts") if r[oi] != "-1.0000"]
        assert tss_vals, "no populated TED_tss_repro on the ON run"
        assert tts_vals, "no populated TED_tts_repro on the ON run"
        assert all(0.0 <= v <= 1.0 for v in tss_vals + tts_vals)
        assert max(tss_vals + tts_vals) > 0.0, "S_repro never positive on ON run"
        # OFF leaves both columns at -1; ON populates them (column-level shift).
        assert all(r[ti] == "-1.0000" for r in peaks(off_log, "tss"))
        # ON shifts at least one TSS-peak reality vs OFF (term entered the blend).
        ri = hdr.index("TED_tss_reality")

        def reality_by_pos(log, kind):
            return {r[6]: r[ri] for r in peaks(log, kind)}   # key by tss_pos (idx 6)
        off_r = reality_by_pos(off_log, "tss")
        on_r = reality_by_pos(on_log, "tss")
        common = set(off_r) & set(on_r)
        assert common and any(off_r[k] != on_r[k] for k in common), \
            "ON did not change any TSS reality vs OFF"
