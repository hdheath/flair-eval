#!/usr/bin/env python3
"""Unit tests for the R2 dataset-unique depth score (DepthNormalizer).

These assert the four guarantees of the per-library, per-end depth
normalization grafted into flair-ted-52826/src/flair/ted.py:

1. PERCENTILE PARSING  — 'pNN' -> percentile, plain int -> absolute count,
   'auto' -> p75, and out-of-range percentiles raise.
2. PER-END INDEPENDENCE — TSS and TTS saturate at independently chosen
   points (different saturation specs -> different saturation counts/scores).
3. SCALE-PARITY (the load-bearing R2 property) — for the default 'pctrank'
   method with a PERCENTILE saturation, multiplying every read count in a
   library by a constant k leaves every peak's depth_score UNCHANGED. This is
   what lets WTC11-fit logreg weights apply to A549 (same transfer-breaking
   risk class as uncalibrated model_score). Simulated by two libraries that
   differ only by a depth scale factor.
4. OLD-BEHAVIOR REPRODUCTION — 'logsat' with a constant reproduces the legacy
   min(1, log1p(n)/log1p(S)) numbers exactly (A/B comparison flag).

Imports DepthNormalizer / parse_depth_saturation directly from the clone's
ted.py (no FLAIR pipeline run needed). Skips cleanly if the clone is absent.
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

_TED_PY = Path(
    "/private/home/hdheath/projects/flair4_ted/flair-ted-52826/src/flair/ted.py"
)
pytestmark = pytest.mark.skipif(
    not _TED_PY.exists(), reason=f"TED clone not found: {_TED_PY}")


def _load_ted():
    """Import ted.py as a standalone module (numpy/scipy only at call time)."""
    spec = importlib.util.spec_from_file_location("ted_under_test", _TED_PY)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ted = _load_ted()
DepthNormalizer = ted.DepthNormalizer
parse_depth_saturation = ted.parse_depth_saturation


# ── 1. percentile / absolute parsing ────────────────────────────────────
def test_parse_percentile():
    assert parse_depth_saturation("p50") == ("pctile", 50.0)
    assert parse_depth_saturation("P90") == ("pctile", 90.0)
    assert parse_depth_saturation("p99.5") == ("pctile", 99.5)


def test_parse_absolute():
    assert parse_depth_saturation("50") == ("abs", 50.0)
    assert parse_depth_saturation("300") == ("abs", 300.0)


def test_parse_auto_is_p75():
    assert parse_depth_saturation("auto") == ("pctile", 75.0)


def test_parse_rejects_bad_percentile():
    for bad in ("p0", "p150", "p-3"):
        with pytest.raises(ValueError):
            parse_depth_saturation(bad)


def test_parse_rejects_nonpositive_absolute():
    with pytest.raises(ValueError):
        parse_depth_saturation("0")


# ── 2. per-end independence ──────────────────────────────────────────────
def test_per_end_independent_saturation():
    """TSS and TTS use independent saturation points -> different scores."""
    rng = np.random.default_rng(0)
    dist = rng.integers(3, 800, size=300).tolist()
    n_tss = DepthNormalizer(dist, method="pctrank", saturation_spec="p50")
    n_tts = DepthNormalizer(dist, method="pctrank", saturation_spec="p90")
    assert n_tss.sat_count != n_tts.sat_count, (
        "TSS p50 and TTS p90 must resolve to different saturation counts")
    # A mid-depth peak scores higher under the lower (p50) saturation.
    mid = int(np.median(dist))
    assert n_tss.score(mid) >= n_tts.score(mid)


def test_per_end_distributions_are_separate():
    """Each end normalizes against its OWN distribution (not a shared one)."""
    tss_dist = [5, 5, 5, 6, 7]           # shallow end
    tts_dist = [500, 600, 700, 800, 900]  # deep end
    n_tss = DepthNormalizer(tss_dist, method="pctrank", saturation_spec="p50")
    n_tts = DepthNormalizer(tts_dist, method="pctrank", saturation_spec="p50")
    # A 600-read peak is saturating for the shallow TSS distribution but
    # only median for the deep TTS distribution.
    assert n_tss.score(600) == pytest.approx(1.0)
    assert n_tts.score(600) < 1.0


# ── 3. SCALE-PARITY (the R2 guarantee) ───────────────────────────────────
@pytest.mark.parametrize("spec", ["p50", "p90", "p25"])
@pytest.mark.parametrize("k", [2, 7, 100])
def test_scale_parity_pctrank_percentile(spec, k):
    """Two libraries differing ONLY by a depth scale factor k must yield
    IDENTICAL depth scores under pctrank + percentile saturation.

    This is the WTC11->A549 transfer guarantee: the depth feature is a
    within-dataset relative quantity (rank), so weights transfer.
    """
    rng = np.random.default_rng(42)
    lib_a = rng.integers(3, 800, size=400).tolist()
    lib_b = [x * k for x in lib_a]   # B = A scaled by k, nothing else changed
    norm_a = DepthNormalizer(lib_a, method="pctrank", saturation_spec=spec)
    norm_b = DepthNormalizer(lib_b, method="pctrank", saturation_spec=spec)
    # Saturation count must scale with the library (tracks its distribution).
    assert norm_b.sat_count == pytest.approx(norm_a.sat_count * k, rel=0.02)
    # Every peak's depth score must be EXACTLY equal across the two libs.
    for a, b in zip(lib_a, lib_b):
        assert norm_a.score(a) == pytest.approx(norm_b.score(b), abs=1e-12)


def test_absolute_saturation_is_not_scale_invariant():
    """Documented contrast: an ABSOLUTE count saturation is a fixed raw
    threshold and is NOT scale-invariant (so percentile is the parity mode)."""
    rng = np.random.default_rng(3)
    lib_a = rng.integers(3, 800, size=300).tolist()
    lib_b = [x * 5 for x in lib_a]
    norm_a = DepthNormalizer(lib_a, method="pctrank", saturation_spec="50")
    norm_b = DepthNormalizer(lib_b, method="pctrank", saturation_spec="50")
    diffs = [abs(norm_a.score(a) - norm_b.score(b)) for a, b in zip(lib_a, lib_b)]
    assert max(diffs) > 0.1, "absolute-count saturation should break parity"


def test_pctrank_scores_in_unit_interval_and_finite():
    rng = np.random.default_rng(9)
    dist = rng.integers(1, 1000, size=200).tolist()
    norm = DepthNormalizer(dist, method="pctrank", saturation_spec="p50")
    for n in (1, 2, 50, 500, 999, 100000):
        s = norm.score(n)
        assert math.isfinite(s) and 0.0 <= s <= 1.0


# ── 4. old-behavior (logsat) reproduction ────────────────────────────────
def test_logsat_reproduces_legacy_formula():
    dist = [10, 20, 30, 40, 50]
    norm = DepthNormalizer(dist, method="logsat", logsat_constant=50)
    for n in (1, 5, 50, 500, 5000):
        expected = min(1.0, math.log1p(n) / math.log1p(50))
        assert norm.score(n) == pytest.approx(expected, abs=1e-12)


def test_logsat_constant_overrides_per_end_spec():
    """For logsat, an explicit constant (the legacy --ted_depth_saturation)
    governs both ends -> exact old-behavior reproduction regardless of spec."""
    dist = [10, 20, 30]
    a = DepthNormalizer(dist, method="logsat", saturation_spec="p50",
                        logsat_constant=50)
    b = DepthNormalizer(dist, method="logsat", saturation_spec="p90",
                        logsat_constant=50)
    assert a.sat_count == b.sat_count == 50.0
    assert a.score(100) == b.score(100)


# ── 5. end-to-end fixture: new flags drive the real pipeline ─────────────
import os  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402

_ENV = Path("/private/home/hdheath/miniforge3/envs/flair-ted-52826")
_FLAIR = _ENV / "bin" / "flair"
_DATA = Path(__file__).resolve().parent / "data"
_MODELS = Path("/private/home/hdheath/projects/artifact_classifier/model_weights")
_TSS_M = _MODELS / "tss" / "xgboost_calibrated_2026-05-28.pkl"
_TTS_M = _MODELS / "tts" / "xgboost_calibrated_2026-05-28.pkl"
_e2e_missing = [str(p) for p in (_FLAIR, _DATA, _TSS_M, _TTS_M) if not p.exists()]

# TSS_depth / TTS_depth columns in the 30-col FlairBed.
_TSS_DEPTH, _TTS_DEPTH = 22, 26


def _run(outdir, extra):
    env = dict(os.environ)
    env["PATH"] = f"{_ENV / 'bin'}:" + env.get("PATH", "")
    out = outdir / "tiny"
    cmd = [
        str(_FLAIR), "transcriptome",
        "-b", str(_DATA / "tiny_reads.bam"), "-g", str(_DATA / "tiny_genome.fa"),
        "-f", str(_DATA / "tiny_annotation.gtf"),
        "--junction_tab", str(_DATA / "tiny_junctions.tab"),
        "-t", "1", "-o", str(out), "--sample_name", "dn_test", "--ted", "--ted_log",
        *extra,
    ]
    res = subprocess.run(cmd, env=env, cwd=str(outdir), capture_output=True, text=True)
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stderr[-1500:]}"
    rows = [ln.split("\t") for ln in (Path(str(out) + ".isoforms.bed")
            ).read_text().splitlines() if ln and not ln.startswith("#")]
    return rows


@pytest.mark.skipif(bool(_e2e_missing),
                    reason=f"e2e prerequisites missing: {_e2e_missing}")
@pytest.mark.parametrize("extra", [
    [],                                                       # default pctrank p50
    ["--ted_depth_saturation_tss", "p90", "--ted_depth_saturation_tts", "p90"],  # percentile
    ["--ted_depth_saturation_tss", "50", "--ted_depth_saturation_tts", "25"],    # absolute, per-end
    ["--ted_depth_norm", "logsat"],                          # old behavior flag
    ["--ted_depth_norm", "logsat", "--ted_depth_saturation", "50"],
])
def test_e2e_depth_flags_finite(extra):
    with tempfile.TemporaryDirectory() as td:
        rows = _run(Path(td), extra)
    assert rows, f"no isoforms emitted for {extra}"
    for r in rows:
        for idx in (_TSS_DEPTH, _TTS_DEPTH):
            v = float(r[idx])
            assert math.isfinite(v) and 0.0 <= v <= 1.0, (
                f"depth col {idx}={v} out of [0,1] for {r[3]} ({extra})")


@pytest.mark.skipif(bool(_e2e_missing),
                    reason=f"e2e prerequisites missing: {_e2e_missing}")
def test_e2e_logsat_reproduces_prior_numbers():
    """--ted_depth_norm logsat --ted_depth_saturation 50 reproduces the
    legacy log1p(n)/log1p(50) depth numbers (recorded in the lab notebook:
    0.9445 / 0.7743 on the tiny fixture)."""
    with tempfile.TemporaryDirectory() as td:
        rows = _run(Path(td), ["--ted_depth_norm", "logsat",
                               "--ted_depth_saturation", "50"])
    depth_vals = {round(float(r[_TSS_DEPTH]), 4) for r in rows} | \
                 {round(float(r[_TTS_DEPTH]), 4) for r in rows}
    # The legacy fixture depth scores are a subset of {0.7743, 0.9445}.
    assert depth_vals <= {0.7743, 0.9445}, (
        f"logsat depth numbers drifted from legacy: {depth_vals}")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v", "--tb=short"]))
