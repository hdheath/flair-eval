#!/usr/bin/env python3
"""5'-end (TSS) trust / scoring behavior for the GLOBAL-only TED.

WHY THIS FILE WAS REWRITTEN
---------------------------
The original version of this file tested the *per-SJC* TED implementation:
HDBSCAN-on-starts + window-on-ends collapse (`ted_collapse_end_groups`), a
`_build_scorer_config` continuous-scoring container, and a
`--trust_ends_5prime` count_sam_transcripts flag that special-cased the 5'
edge. The GLOBAL-only port (flair-ted-52826, branch ted-port-52826) DELETED
all of that: there is no per-SJC clustering, no `ted_collapse_end_groups`, no
`_build_scorer_config`, and count_sam_transcripts is Colette's dev version
(no `trust_ends_5prime`). So every import the old file relied on is gone and
the file failed on a clean checkout (14 failed + 29 errors).

WHAT THE OLD FILE WAS REALLY ABOUT — and where it lives now
-----------------------------------------------------------
The biological intent was 5'/TSS end *trust* and *scoring*:
  * the 5' end should be model-scored / validated,
  * the 3' end is handled separately (and, for dRNA, trusted),
  * the per-end reality score depends on depth + model + annotation.

In the global-only TED that intent is realized by:
  * ``ted.add_ted_args`` — registers ``--ted_model_ends`` (which ends use the
    XGBoost boundary model) and ``--ted_trust_3prime`` (dRNA 3'-trust),
  * ``ted.TedScorer`` — holds the per-end (tss/tts) models + annotated ends,
  * ``ted._ted_validate_cluster`` — THE reality blend:
        reality = w_depth*depth + w_model*model + w_annot*annot
    with TSS=5' and TTS=3' handled strand-aware, the TTS model short-circuited
    to 1.0 under ``trust_3prime``.

This file now tests that SAME intent against those current symbols. It imports
``ted.py`` standalone (numpy/scipy only) exactly like test_ted_depth_norm.py,
so no FLAIR pipeline / pysam / xgboost is needed. Models and the genome
sequence extractor are duck-typed mocks; annotated ends come from the real
tiny GTF fixture via the real ``parse_gtf_ends``.

COVERAGE DELTA (intent that NO LONGER EXISTS in global-only TED, dropped here):
  * count_sam_transcripts ``--trust_ends_5prime`` / ``check_firstlastexon`` /
    ``check_exonenddist`` 5'-relaxed-3'-strict edge logic — removed with the
    old count_sam edits (dev count_sam is used as-is).
  * ``ted_collapse_end_groups`` HDBSCAN-starts/window-ends collapse and its
    [score,start,end,strand,name,[names]] firstpass/final output contracts.
  * ``_build_scorer_config`` / ``_build_scorer_from_config`` picklable config.
  * the SQANTI junction_match_type (FSM/NIC/NNC) boost in the confidence blend
    — that blend was the per-SJC collapse score, not the global reality score.
  * ``--ted`` preset auto-setting ``trust_ends_5prime`` in add_preset_args.
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pytest

# ─────────────────────────────────────────────────────────────────────────────
# Load the GLOBAL-only ted.py standalone (numpy/scipy only at call time).
# ─────────────────────────────────────────────────────────────────────────────
_TED_PY = Path(
    "/private/home/hdheath/projects/flair4_ted/flair-ted-52826/src/flair/ted.py"
)
pytestmark = pytest.mark.skipif(
    not _TED_PY.exists(), reason=f"TED clone not found: {_TED_PY}")

DATA_DIR = Path(__file__).resolve().parent / "data"
_TINY_GTF = DATA_DIR / "tiny_annotation.gtf"


def _load_ted():
    spec = importlib.util.spec_from_file_location("ted_under_test_5prime", _TED_PY)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ted = _load_ted()
add_ted_args = ted.add_ted_args
TedScorer = ted.TedScorer
parse_gtf_ends = ted.parse_gtf_ends
_ted_validate_cluster = ted._ted_validate_cluster


# ─────────────────────────────────────────────────────────────────────────────
# Test doubles: a duck-typed XGBoost model + sequence-feature extractor.
# _ted_validate_cluster only ever calls extractor.extract_batch([(chrom,pos,strand)])
# and model.predict_proba(X)[0, 1], so these stand in without xgboost/pysam.
# ─────────────────────────────────────────────────────────────────────────────
class _FixedModel:
    """predict_proba returns a fixed P(positive) for every row."""

    def __init__(self, p):
        self._p = float(p)

    def predict_proba(self, X):
        n = len(X)
        return np.array([[1.0 - self._p, self._p]] * n)


class _RecordingExtractor:
    """Returns one finite feature row per query and records the (chrom,pos,strand)
    tuples it was asked to extract — so a test can assert WHICH end was scored."""

    def __init__(self):
        self.calls = []

    def extract_batch(self, queries):
        self.calls.extend(queries)
        return np.zeros((len(queries), 4), dtype=float)


def _make_scorer(*, tss_p=None, tts_p=None, annotated_ends=None):
    """Build a TedScorer that mirrors the --ted_model_ends gate: a model is
    present for an end ONLY when its probability is given (None = that end's
    model was not loaded, e.g. --ted_model_ends=tss leaves tts_model=None)."""
    tss_model = _FixedModel(tss_p) if tss_p is not None else None
    tts_model = _FixedModel(tts_p) if tts_p is not None else None
    return TedScorer(
        tss_model=tss_model,
        tts_model=tts_model,
        flair_extractor_tss=_RecordingExtractor() if tss_model is not None else None,
        flair_extractor_tts=_RecordingExtractor() if tts_model is not None else None,
        annotated_ends=annotated_ends or {},
    )


# A junction-chain-sized cluster of read 5'/3' ends near a position.
def _cluster(center, n=10, jitter=2):
    return [center + (i % (jitter + 1)) for i in range(n)]


# ═════════════════════════════════════════════════════════════════════════════
# 1. The --ted_model_ends / --ted_trust_3prime CLI surface (add_ted_args)
#    Replaces the old "TestAddPresetArgs / --ted sets trust_ends_5prime".
# ═════════════════════════════════════════════════════════════════════════════
class TestTedEndControlArgs:
    def _parse(self, argv):
        p = argparse.ArgumentParser()
        add_ted_args(p)
        return p.parse_args(argv)

    def test_model_ends_default_is_both(self):
        args = self._parse(['--ted'])
        assert args.ted_model_ends == 'both'

    def test_model_ends_accepts_per_end_choices(self):
        for choice in ('both', 'tss', 'tts', 'none'):
            assert self._parse(['--ted', '--ted_model_ends', choice]).ted_model_ends == choice

    def test_model_ends_rejects_garbage(self):
        with pytest.raises(SystemExit):
            self._parse(['--ted', '--ted_model_ends', 'fiveprime'])

    def test_trust_3prime_default_off(self):
        assert self._parse(['--ted']).ted_trust_3prime is False

    def test_trust_3prime_flag_sets_true(self):
        assert self._parse(['--ted', '--ted_trust_3prime']).ted_trust_3prime is True

    def test_model_ends_and_trust_3prime_independent(self):
        args = self._parse(['--ted', '--ted_model_ends', 'tss', '--ted_trust_3prime'])
        assert args.ted_model_ends == 'tss'
        assert args.ted_trust_3prime is True


# ═════════════════════════════════════════════════════════════════════════════
# 2. The reality blend (_ted_validate_cluster) — formula + bounds.
#    Replaces the old "test_confidence_weight_formula" composite-score tests,
#    retargeted to the global reality score: w_depth*depth + w_model*model
#    + w_annot*annot.
# ═════════════════════════════════════════════════════════════════════════════
class TestRealityBlend:
    def test_reality_formula_matches_weights(self):
        scorer = _make_scorer(tss_p=0.8)
        positions = _cluster(1100, n=10)
        ok, reality, scores = _ted_validate_cluster(
            positions, 'tss', 'chr_test', '+', scorer,
            w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.50,
            depth_saturation=50.0)
        expected = (0.40 * scores['depth_score']
                    + 0.35 * scores['model_score']
                    + 0.25 * scores['annotation_score'])
        assert scores['model_score'] == pytest.approx(0.8)
        assert reality == pytest.approx(expected)

    def test_reality_is_bounded_0_1(self):
        scorer = _make_scorer(tss_p=1.0)
        positions = _cluster(1100, n=40)
        _, reality, _ = _ted_validate_cluster(
            positions, 'tss', 'chr_test', '+', scorer,
            w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.50,
            depth_saturation=50.0)
        assert 0.0 <= reality <= 1.0

    def test_threshold_pass_fail(self):
        scorer = _make_scorer(tss_p=0.9)
        positions = _cluster(1100, n=30)
        ok_lo, reality, _ = _ted_validate_cluster(
            positions, 'tss', 'chr_test', '+', scorer,
            w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.01,
            depth_saturation=50.0)
        ok_hi, _, _ = _ted_validate_cluster(
            positions, 'tss', 'chr_test', '+', scorer,
            w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.99,
            depth_saturation=50.0)
        assert bool(ok_lo) is True and bool(ok_hi) is False
        assert 0.01 <= reality <= 0.99


# ═════════════════════════════════════════════════════════════════════════════
# 3. 5' (TSS) model score AFFECTS 5' reality — the core surviving intent.
#    Replaces "test_scorer_boosts_confidence / test_scorer_penalizes_bad_ends".
# ═════════════════════════════════════════════════════════════════════════════
class TestFivePrimeModelScoring:
    def _tss_reality(self, tss_p):
        scorer = _make_scorer(tss_p=tss_p)
        positions = _cluster(1100, n=10)
        _, reality, scores = _ted_validate_cluster(
            positions, 'tss', 'chr_test', '+', scorer,
            w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.50,
            depth_saturation=50.0)
        return reality, scores

    def test_high_5prime_model_boosts_reality(self):
        low, _ = self._tss_reality(0.1)
        high, _ = self._tss_reality(0.9)
        assert high > low, f"high 5' model score ({high:.3f}) should beat low ({low:.3f})"

    def test_5prime_model_score_threaded_through(self):
        _, scores = self._tss_reality(0.73)
        assert scores['model_score'] == pytest.approx(0.73)

    def test_no_5prime_model_zeroes_model_term(self):
        """--ted_model_ends excluding TSS -> no tss_model -> model_score 0."""
        scorer = _make_scorer(tss_p=None)  # TSS model not loaded
        positions = _cluster(1100, n=10)
        _, _, scores = _ted_validate_cluster(
            positions, 'tss', 'chr_test', '+', scorer,
            w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.50,
            depth_saturation=50.0)
        assert scores['model_score'] == 0.0


# ═════════════════════════════════════════════════════════════════════════════
# 4. --ted_model_ends=tss : ONLY the TSS is model-scored; the TTS is not.
#    Direct replacement for the old "--ted sets 5'-only trust" intent.
# ═════════════════════════════════════════════════════════════════════════════
class TestModelEndsTssOnly:
    def test_tss_scored_tts_not_when_model_ends_tss(self):
        # model_ends=tss => tss_model present, tts_model None (mirrors _build_ted_scorer).
        scorer = _make_scorer(tss_p=0.85, tts_p=None)

        _, _, tss_scores = _ted_validate_cluster(
            _cluster(1100), 'tss', 'chr_test', '+', scorer,
            w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.50,
            depth_saturation=50.0)
        _, _, tts_scores = _ted_validate_cluster(
            _cluster(3500), 'tts', 'chr_test', '+', scorer,
            w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.50,
            depth_saturation=50.0)

        assert tss_scores['model_score'] == pytest.approx(0.85), "TSS must be model-scored"
        assert tts_scores['model_score'] == 0.0, "TTS must NOT be model-scored under model_ends=tss"

    def test_extractor_called_only_for_tss(self):
        scorer = _make_scorer(tss_p=0.85, tts_p=None)
        _ted_validate_cluster(
            _cluster(1100), 'tss', 'chr_test', '+', scorer,
            w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.50,
            depth_saturation=50.0)
        assert len(scorer._flair_extractor_tss.calls) == 1
        assert scorer._flair_extractor_tts is None


# ═════════════════════════════════════════════════════════════════════════════
# 5. Strand-aware TSS=5' / TTS=3' coordinate convention.
#    Replaces "test_tts_position_correct_plus/minus_strand".
#    The classifier is trained on the inclusive boundary base, so a closing
#    (exclusive) edge is model-scored at centroid-1; the open edge at centroid.
# ═════════════════════════════════════════════════════════════════════════════
class TestStrandAwareEndConvention:
    def _scored_pos(self, end_type, strand):
        scorer = _make_scorer(tss_p=0.5, tts_p=0.5)
        centroid = 2000
        _ted_validate_cluster(
            [centroid] * 6, end_type, 'chr_test', strand, scorer,
            w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.50,
            depth_saturation=50.0)
        extractor = (scorer._flair_extractor_tts if end_type == 'tts'
                     else scorer._flair_extractor_tss)
        return extractor.calls[0][1], centroid  # model_pos, centroid

    def test_plus_strand_tss_is_open_edge(self):
        # + strand TSS = read start (inclusive) -> model_pos == centroid
        pos, centroid = self._scored_pos('tss', '+')
        assert pos == centroid

    def test_plus_strand_tts_is_closing_edge(self):
        # + strand TTS = read end (exclusive) -> model_pos == centroid - 1
        pos, centroid = self._scored_pos('tts', '+')
        assert pos == centroid - 1

    def test_minus_strand_tss_is_closing_edge(self):
        # - strand TSS = read end (exclusive) -> model_pos == centroid - 1
        pos, centroid = self._scored_pos('tss', '-')
        assert pos == centroid - 1

    def test_minus_strand_tts_is_open_edge(self):
        # - strand TTS = read start (inclusive) -> model_pos == centroid
        pos, centroid = self._scored_pos('tts', '-')
        assert pos == centroid


# ═════════════════════════════════════════════════════════════════════════════
# 6. dRNA 3'-trust complement (--ted_trust_3prime).
#    The biological complement of 5' trust: for dRNA the 3' end is real, so the
#    TTS model is short-circuited to 1.0 while the 5' (TSS) keeps full scoring.
#    Replaces the old per-SJC asymmetric-trust scenario tests.
# ═════════════════════════════════════════════════════════════════════════════
class TestDrnaThreePrimeTrust:
    def test_trust_3prime_forces_tts_model_high(self):
        # tts_p deliberately LOW (0.05) but trust must override it to 1.0.
        scorer = _make_scorer(tss_p=0.5, tts_p=0.05)
        _, _, scores = _ted_validate_cluster(
            _cluster(3500), 'tts', 'chr_test', '+', scorer,
            w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.50,
            depth_saturation=50.0, trust_3prime=True)
        assert scores['model_score'] == 1.0, "trust_3prime must force TTS model_score to 1.0"

    def test_trust_3prime_does_not_touch_tss(self):
        # The 5' end keeps full model scoring even with trust_3prime on.
        scorer = _make_scorer(tss_p=0.42, tts_p=0.05)
        _, _, scores = _ted_validate_cluster(
            _cluster(1100), 'tss', 'chr_test', '+', scorer,
            w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.50,
            depth_saturation=50.0, trust_3prime=True)
        assert scores['model_score'] == pytest.approx(0.42), \
            "trust_3prime must NOT alter the 5' (TSS) model score"

    def test_trust_3prime_raises_tts_reality_vs_low_model(self):
        positions = _cluster(3500, n=10)
        scorer = _make_scorer(tss_p=0.5, tts_p=0.05)
        kw = dict(w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.50,
                  depth_saturation=50.0)
        _, r_untrusted, _ = _ted_validate_cluster(
            positions, 'tts', 'chr_test', '+', scorer, **kw)
        _, r_trusted, _ = _ted_validate_cluster(
            positions, 'tts', 'chr_test', '+', scorer, trust_3prime=True, **kw)
        assert r_trusted > r_untrusted, \
            "dRNA 3'-trust should raise TTS reality (model term 0.05 -> 1.0)"


# ═════════════════════════════════════════════════════════════════════════════
# 7. Annotation term uses real annotated TSS/TTS from the tiny GTF fixture.
#    Confirms 5' (TSS) reality rises when the cluster sits ON an annotated TSS,
#    and that the annotation distance is end-type / strand specific.
# ═════════════════════════════════════════════════════════════════════════════
class TestAnnotationProximity:
    @pytest.fixture(scope="class")
    def annotated_ends(self):
        assert _TINY_GTF.exists(), f"missing fixture {_TINY_GTF}"
        return parse_gtf_ends(str(_TINY_GTF))

    def test_gtf_parses_per_end_annotations(self, annotated_ends):
        assert 'chr_test' in annotated_ends
        assert annotated_ends['chr_test']['tss'], "expected annotated TSS positions"
        assert annotated_ends['chr_test']['tts'], "expected annotated TTS positions"

    def test_on_annotated_tss_scores_higher_than_far(self, annotated_ends):
        scorer = _make_scorer(tss_p=0.5, annotated_ends=annotated_ends)
        # GENE1 (+) starts at GTF 1001 -> 0-based annotated TSS == 1000.
        on_tss = annotated_ends['chr_test']['tss'][0]
        kw = dict(w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.50,
                  depth_saturation=50.0, annot_half_life=100.0, annot_max_dist=500)
        _, _, near = _ted_validate_cluster(
            [on_tss] * 8, 'tss', 'chr_test', '+', scorer, **kw)
        _, _, far = _ted_validate_cluster(
            [on_tss + 2000] * 8, 'tss', 'chr_test', '+', scorer, **kw)
        assert near['annotation_score'] > far['annotation_score']
        assert near['annotation_score'] == pytest.approx(1.0)  # exp(-0/half_life)
        assert far['annotation_score'] == 0.0                  # beyond annot_max_dist

    def test_annotation_distance_reported(self, annotated_ends):
        scorer = _make_scorer(tss_p=0.5, annotated_ends=annotated_ends)
        on_tss = annotated_ends['chr_test']['tss'][0]
        _, _, scores = _ted_validate_cluster(
            [on_tss + 30] * 8, 'tss', 'chr_test', '+', scorer,
            w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.50,
            depth_saturation=50.0, annot_half_life=100.0, annot_max_dist=500)
        assert scores['annot_dist_bp'] == 30


# ═════════════════════════════════════════════════════════════════════════════
# 8. Depth term: 5' reality rises with read support (more reads = higher depth).
#    Replaces "test_density_is_primary_weight" (density-dominant) intent.
# ═════════════════════════════════════════════════════════════════════════════
class TestDepthSupport:
    def test_more_reads_higher_depth_score(self):
        scorer = _make_scorer(tss_p=0.5)
        kw = dict(w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.50,
                  depth_saturation=50.0)
        _, _, few = _ted_validate_cluster(
            _cluster(1100, n=3), 'tss', 'chr_test', '+', scorer, **kw)
        _, _, many = _ted_validate_cluster(
            _cluster(1100, n=40), 'tss', 'chr_test', '+', scorer, **kw)
        assert many['depth_score'] > few['depth_score']

    def test_depth_dominant_weight_drives_reality(self):
        """With w_depth the largest weight, a high-depth/low-model cluster can
        out-score a low-depth/high-model one — depth is the primary signal."""
        scorer_lowmodel = _make_scorer(tss_p=0.1)
        scorer_himodel = _make_scorer(tss_p=0.9)
        kw = dict(w_depth=0.40, w_model=0.35, w_annot=0.25, threshold=0.50,
                  depth_saturation=50.0)
        _, deep_lowmodel, _ = _ted_validate_cluster(
            _cluster(1100, n=50), 'tss', 'chr_test', '+', scorer_lowmodel, **kw)
        _, shallow_himodel, _ = _ted_validate_cluster(
            _cluster(1100, n=2), 'tss', 'chr_test', '+', scorer_himodel, **kw)
        assert deep_lowmodel > shallow_himodel


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
