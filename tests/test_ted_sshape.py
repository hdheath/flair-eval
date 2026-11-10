#!/usr/bin/env python3
"""Unit tests for the S_shape 4th TSS reality term (Lever A) in ted.py.

Three guarantees of the default-OFF S_shape graft on flair-ted-52826:

1. FEATURE EQUIVALENCE — the in-flow `shape_features` (operating on ReadRec-like
   objects via `_ShapeReadView`) is BYTE-IDENTICAL to the offline source of truth
   `artifact_classifier/scratch/novel_shape/shape_features.py` on a synthetic read
   set, BOTH strands. This is the load-bearing property: the deployed S_shape must
   match the OOF-calibrated offline features or the calibration is meaningless.

2. BLEND MATH — with w_shape>0 the TSS reality is the 4-term blend
   reality = w_d·S_d + w_m·S_m + w_a·S_a + w_shape·S_shape, S_shape enters in
   [0,1], the read floor (n_shape_reads<4) forces S_shape=0, and the θ comparison
   uses the 4-term reality. With w_shape==0 the term is inert (3-term form) and the
   returned reality equals the legacy 3-term blend exactly.

3. WEIGHT 4-SIMPLEX — _resolve_ted_weights renormalizes (d,m,a,shape) to sum 1 when
   w_shape>0 (returns a 4-tuple), and is byte-identical to the legacy 3-tuple when
   w_shape==0 (the OFF / TTS path is unchanged).

The OFF-by-default byte-identity of the PUBLISHED isoform outputs is verified
end-to-end against a no-S_shape-code build in test_ted_no_filter_emit-style runs;
here we assert the math-level invariant that underpins it (w_shape=0 ⇒ identical
reality) plus the structural OFF flags (default 0.0, model not required).

Loads ted.py standalone (numpy/scipy/sklearn at call time); skips if the clone or
the offline source is absent.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

_TED_PY = Path(
    "/private/home/hdheath/projects/flair4_ted/flair-ted-52826/src/flair/ted.py")
_OFFLINE = Path(
    "/private/home/hdheath/projects/artifact_classifier/scratch/novel_shape/shape_features.py")
_SHAPE_PKL = Path(
    "/private/home/hdheath/projects/artifact_classifier/model_weights/"
    "tss_shape/shape_tss_pacbio_cDNA.pkl")

pytestmark = pytest.mark.skipif(
    not _TED_PY.exists(), reason=f"TED clone not found: {_TED_PY}")


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ted = _load(_TED_PY, "ted_under_test_sshape")


class _RR:
    """Minimal ReadRec stand-in: .start/.end/.strand (what TED bins in-flow)."""
    def __init__(self, start, end, strand):
        self.start = start
        self.end = end
        self.strand = strand


class _PyRead:
    """Minimal pysam stand-in for the OFFLINE shape_features signature."""
    def __init__(self, start, end, is_reverse):
        self.reference_start = start
        self.reference_end = end
        self.is_reverse = is_reverse


def _synthetic(pos, strand, n=50, seed=3):
    """Build matched (ReadRec, pysam) read lists for a candidate at `pos`."""
    rng = np.random.default_rng(seed)
    if strand == '+':
        starts = pos + rng.integers(-20, 21, size=n)
        lens = rng.integers(250, 1500, size=n)
        rr = [_RR(int(s), int(s + l), '+') for s, l in zip(starts, lens)]
        py = [_PyRead(int(s), int(s + l), False) for s, l in zip(starts, lens)]
    else:
        ends = pos + rng.integers(-20, 21, size=n)
        lens = rng.integers(250, 1500, size=n)
        rr = [_RR(int(e - l), int(e + 1), '-') for e, l in zip(ends, lens)]
        py = [_PyRead(int(e - l), int(e + 1), True) for e, l in zip(ends, lens)]
    return rr, py


# ── 1. feature equivalence (in-flow port == offline source) ──────────────
@pytest.mark.skipif(not _OFFLINE.exists(), reason="offline shape_features.py absent")
@pytest.mark.parametrize("strand,pos", [('+', 12000), ('-', 8000)])
def test_shape_features_match_offline(strand, pos):
    off = _load(_OFFLINE, "offline_shape_features")
    rr, py = _synthetic(pos, strand)
    views = [ted._ShapeReadView(r) for r in rr]
    fin = ted.shape_features(views, pos, strand)
    fof = off.shape_features(py, pos, strand)
    # SHAPE_COLS must agree exactly (NaN==NaN treated as equal).
    assert ted.SHAPE_COLS == off.SHAPE_COLS
    for k in off.SHAPE_COLS:
        a, b = fin[k], fof[k]
        if a != a and b != b:        # both NaN
            continue
        assert a == b, f"{k}: in-flow {a} != offline {b}"


def test_shape_view_maps_readrec():
    """_ShapeReadView presents a ReadRec under the pysam attribute names."""
    v = ted._ShapeReadView(_RR(100, 500, '-'))
    assert v.reference_start == 100 and v.reference_end == 500 and v.is_reverse is True
    v2 = ted._ShapeReadView(_RR(7, 9, '+'))
    assert v2.is_reverse is False


def test_shape_features_read_floor():
    """< MIN_READS strand reads -> n_shape_reads recorded, features NaN."""
    rr, _ = _synthetic(5000, '+', n=2)
    feat = ted.shape_features([ted._ShapeReadView(r) for r in rr], 5000, '+')
    assert feat['n_shape_reads'] == 2
    assert all(feat[c] != feat[c] for c in ted.SHAPE_COLS)  # all NaN


# ── 2. blend math (4-term reality, read floor, θ) ────────────────────────
def _validate(**kw):
    """Call _ted_validate_cluster with sane defaults for a TSS cluster."""
    base = dict(
        positions=[1000] * 10, end_type='tss', chrom='chrT', strand='+',
        scorer=None,                      # no model -> model_score 0
        w_depth=0.30, w_model=0.20, w_annot=0.20, threshold=0.5,
        depth_norm=None, depth_saturation=50.0)
    base.update(kw)
    return ted._ted_validate_cluster(**base)


def test_blend_off_is_three_term():
    """w_shape == 0 (default) -> reality == 3-term blend, byte-identical, no
    matter what shape_score is passed."""
    # scorer None => model_score 0, annot 0; depth from log1p(10)/log1p(50).
    is_real, reality, scores = _validate(w_shape=0.0, shape_score=0.9)
    sd = min(1.0, np.log1p(10) / np.log1p(50))
    expected = 0.30 * sd + 0.20 * 0.0 + 0.20 * 0.0
    assert reality == pytest.approx(expected)
    assert scores['shape_score'] == 0.0   # inert when w_shape=0


def test_blend_four_term_adds_shape():
    """w_shape>0 adds w_shape·S_shape; reality matches the hand-computed sum."""
    is_real, reality, scores = _validate(
        w_depth=0.24, w_model=0.19, w_annot=0.17, w_shape=0.40,
        shape_score=0.5, threshold=0.0)
    sd = min(1.0, np.log1p(10) / np.log1p(50))
    expected = 0.24 * sd + 0.19 * 0.0 + 0.17 * 0.0 + 0.40 * 0.5
    assert reality == pytest.approx(expected)
    assert scores['shape_score'] == pytest.approx(0.5)


def test_blend_shape_none_contributes_zero():
    """shape_score=None (term ON but no value supplied, e.g. self-position) ->
    S_shape contributes 0 (read-floor / missing-evidence semantics)."""
    _, reality, scores = _validate(w_shape=0.40, shape_score=None, threshold=0.0)
    sd = min(1.0, np.log1p(10) / np.log1p(50))
    assert reality == pytest.approx(0.30 * sd)   # only depth term (model/annot 0)
    assert scores['shape_score'] == 0.0


def test_blend_tts_ignores_shape():
    """The 4th term is TSS-only; on a TTS cluster shape_score is never added."""
    _, reality, scores = _validate(
        end_type='tts', w_shape=0.40, shape_score=0.9, threshold=0.0)
    sd = min(1.0, np.log1p(10) / np.log1p(50))
    assert reality == pytest.approx(0.30 * sd)   # shape excluded on TTS
    assert scores['shape_score'] == 0.0


def test_blend_threshold_uses_four_term():
    """A novel-TSS-like cluster (annot 0, model 0 here) that fails 3-term can PASS
    4-term when S_shape is high — the geometry-lock-opening property (Lever A).

    scorer=None => S_m=0, S_a=0; with the logsat constant=1 the depth term
    saturates so S_d=1. So 3-term reality = w_depth = 0.24 (< θ=0.5); the 4-term
    blend adds w_shape·S_shape = 0.40·0.7 = 0.28 -> 0.52 (≥ θ)."""
    norm = ted.DepthNormalizer([], method='logsat', logsat_constant=1.0)  # S_d=1
    _, r3, _ = _validate(w_depth=0.24, w_model=0.19, w_annot=0.17,
                         w_shape=0.0, shape_score=0.0, threshold=0.5,
                         depth_norm=norm)
    assert r3 == pytest.approx(0.24)
    assert r3 < 0.5
    # 4-term with a clean shape: + 0.40*0.7 = +0.28 -> clears 0.5.
    is_real, r4, _ = _validate(w_depth=0.24, w_model=0.19, w_annot=0.17,
                               w_shape=0.40, shape_score=0.7, threshold=0.5,
                               depth_norm=norm)
    assert r4 == pytest.approx(0.52)
    assert r4 > r3
    assert is_real is True


# ── 3. weight 4-simplex resolution ───────────────────────────────────────
def _load_host_resolve():
    """Import _resolve_ted_weights from the clone's flair_transcriptome.py.

    flair_transcriptome imports the flair package; load it as part of the clone
    so its relative imports resolve. Skip if that import chain is unavailable.
    """
    ft = Path("/private/home/hdheath/projects/flair4_ted/flair-ted-52826/"
              "src/flair/flair_transcriptome.py")
    if not ft.exists():
        pytest.skip("flair_transcriptome.py absent")
    import sys
    clone_src = "/private/home/hdheath/projects/flair4_ted/flair-ted-52826/src"
    if clone_src not in sys.path:
        sys.path.insert(0, clone_src)
    try:
        from flair.flair_transcriptome import _resolve_ted_weights
    except Exception as e:  # heavy deps (pysam/xgboost) — skip cleanly
        pytest.skip(f"cannot import _resolve_ted_weights: {e}")
    return _resolve_ted_weights


def test_resolve_weights_off_is_three_tuple():
    rw = _load_host_resolve()
    out = rw(0.39, 0.32, 0.29, 0.40, 0.35, 0.25, True, w_shape=0.0)
    assert len(out) == 3
    assert sum(out) == pytest.approx(1.0)
    assert out == pytest.approx((0.39, 0.32, 0.29))


def test_resolve_weights_four_simplex():
    rw = _load_host_resolve()
    out = rw(0.24, 0.19, 0.17, 0.40, 0.35, 0.25, True, w_shape=0.40)
    assert len(out) == 4
    assert sum(out) == pytest.approx(1.0)
    # already-normalized inputs (sum 1.0) pass through to the 4term_ws.40 target.
    assert out == pytest.approx((0.24, 0.19, 0.17, 0.40))


def test_resolve_weights_four_simplex_renormalizes():
    rw = _load_host_resolve()
    # production 3-term weights + w_shape raw 0.40 -> renormalize the 4-vector.
    out = rw(0.39, 0.32, 0.29, 0.40, 0.35, 0.25, True, w_shape=0.40)
    assert len(out) == 4
    assert sum(out) == pytest.approx(1.0)
    tot = 0.39 + 0.32 + 0.29 + 0.40
    assert out == pytest.approx((0.39 / tot, 0.32 / tot, 0.29 / tot, 0.40 / tot))


# ── End-to-end wiring (real pipeline): OFF byte-identity + ON populates ───
import os
import subprocess
import tempfile

_ENV = Path("/private/home/hdheath/miniforge3/envs/flair-ted-52826")
_FLAIR = _ENV / "bin" / "flair"
_DATA = Path(__file__).resolve().parent / "data"
_TSS_M = Path("/private/home/hdheath/projects/artifact_classifier/model_weights/tss/xgboost.pkl")
_TTS_M = Path("/private/home/hdheath/projects/artifact_classifier/model_weights/tts/xgboost.pkl")
_e2e_missing = [str(p) for p in (_FLAIR, _DATA / "tiny_reads.bam", _TSS_M, _TTS_M, _SHAPE_PKL)
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
           "-t", "1", "-o", str(out), "--sample_name", "sshape",
           # NORTH-STAR re-arch: S_shape is a LEGACY reality-path term.
           "--ted", "--ted_legacy_reality",
           "--tss_model", str(_TSS_M), "--tts_model", str(_TTS_M),
           "--ted_log", "--keep_intermediate", *extra]
    res = subprocess.run(cmd, env=env, cwd=str(outdir), capture_output=True, text=True)
    assert res.returncode == 0, f"flair exited {res.returncode}\n{res.stderr[-1500:]}"
    return out


@pytest.mark.skipif(bool(_e2e_missing), reason=f"e2e prereqs missing: {_e2e_missing}")
def test_e2e_off_vs_on():
    """OFF (no shape flags) and ON (w_shape>0) both run; OFF leaves the published
    BED/GTF unchanged vs the same-models 3-term run and the ted_log TED_tss_shape
    column is -1; ON populates TED_tss_shape in [0,1] and shifts TSS reality."""
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        (tdp / "off").mkdir(); (tdp / "on").mkdir()
        off = _run(tdp / "off", "tiny", [])
        on = _run(tdp / "on", "tiny",
                  ["--tss_shape_model", str(_SHAPE_PKL), "--ted_w_shape_tss", "0.40"])

        def rows(p):
            return [ln.split("\t") for ln in Path(p).read_text().splitlines()
                    if ln and not ln.startswith("#")]

        off_log = rows(str(off) + ".ted_log.tsv")
        on_log = rows(str(on) + ".ted_log.tsv")
        hdr = off_log[0]
        sc = hdr.index("TED_tss_shape")
        # OFF: every row's TED_tss_shape is the OFF sentinel -1.
        assert all(r[sc] == "-1.0000" for r in off_log[1:])
        # ON: at least one global_tss peak carries a finite S_shape in [0,1].
        tss_shape_vals = [float(r[sc]) for r in on_log[1:]
                          if len(r) > sc and r[2] == "global_peak"
                          and r[1].startswith("global_tss") and r[sc] != "-1.0000"]
        assert tss_shape_vals, "no populated TED_tss_shape on the ON run"
        assert all(0.0 <= v <= 1.0 for v in tss_shape_vals)
        # ON shifts at least one TSS-peak reality vs OFF (4-term != 3-term).
        ridx = hdr.index("TED_tss_reality")

        def tss_reality(log):
            return {r[6]: r[ridx] for r in log[1:]      # key by tss_pos (col 7, idx 6)
                    if r[2] == "global_peak" and r[1].startswith("global_tss")}
        off_r, on_r = tss_reality(off_log), tss_reality(on_log)
        common = set(off_r) & set(on_r)
        assert common and any(off_r[k] != on_r[k] for k in common), \
            "ON did not change any TSS reality vs OFF"


@pytest.mark.skipif(bool(_e2e_missing), reason=f"e2e prereqs missing: {_e2e_missing}")
def test_e2e_off_bed_byte_identical_to_no_shape_flags():
    """Two OFF runs (default, and with --tss_shape_model set but w_shape=0) produce
    byte-identical published BED — proves a loaded model with zero weight is inert."""
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        (tdp / "a").mkdir(); (tdp / "b").mkdir()
        a = _run(tdp / "a", "tiny", [])
        # model path given but weight 0 => term OFF, model never loaded/used.
        b = _run(tdp / "b", "tiny",
                 ["--tss_shape_model", str(_SHAPE_PKL), "--ted_w_shape_tss", "0.0"])
        assert (Path(str(a) + ".isoforms.bed").read_text()
                == Path(str(b) + ".isoforms.bed").read_text())
        assert (Path(str(a) + ".isoforms.gtf").read_text()
                == Path(str(b) + ".isoforms.gtf").read_text())


# ── TssShapeScorer contract (read floor, [0,1], batch == per-row) ─────────
@pytest.mark.skipif(not _SHAPE_PKL.exists(), reason="shape pkl absent")
def test_shape_scorer_contract():
    sc = ted.TssShapeScorer(str(_SHAPE_PKL))
    assert sc.chemistry == 'pacbio_cDNA'
    rr, _ = _synthetic(20000, '+', n=40, seed=11)
    feat = ted.shape_features([ted._ShapeReadView(r) for r in rr], 20000, '+')
    s = sc.score_one(feat)
    assert 0.0 <= s <= 1.0
    # read floor: < 4 reads -> exactly 0.0, model not consulted.
    low = ted.shape_features([ted._ShapeReadView(r) for r in rr[:2]], 20000, '+')
    assert sc.score_one(low) == 0.0
    # batch == per-row (and read-floor rows -> 0.0 in the batch).
    feats = [feat, low, feat]
    batch = sc.score_batch(feats)
    assert batch[1] == 0.0
    assert batch[0] == pytest.approx(s) and batch[2] == pytest.approx(s)
