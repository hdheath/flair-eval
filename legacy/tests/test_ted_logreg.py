#!/usr/bin/env python3
"""Regression tests for analysis/ted_logreg/ted_logreg.py.

Pins three bugs that were fixed in the L1-clip / class-weight / Youden's-J
patch:

  1. Silent clip-to-zero of negative regression coefficients (lines 241-246
     in the pre-fix code). Caused the model coefficient to vanish on real
     data because TP and FP both have high model_score, making the raw
     coefficient slightly negative under class_weight=balanced. The fix
     keeps both signed weights (for diagnostics) and clipped non-negative
     deploy weights (for TED's JSON), and surfaces both via the pipeline
     object and the coefficients TSV.

  2. class_weight="balanced" over-correcting the mild 1:2 imbalance and
     biasing the model coefficient negative. Removed in the fix.

  3. Youden's-J threshold selection over-penalising FPs at high TPR.
     Replaced with F1-on-held-out-folds (5-fold KFold) to match the
     metric the confusion matrix evaluation uses.

Strategy: build a synthetic 1000-row fixture where
  - TED_*_model is CONSTANT across rows (deterministically non-informative),
  - TED_*_annot perfectly separates joint_tp=1 from joint_tp=0,
  - TED_depth is mildly informative.

Assertions verify that
  - The deploy weight for model is zero (because the raw coefficient is
    ~0 — there's no variance) or near-zero,
  - |w_annot_deploy| > |w_depth_deploy| > |w_model_deploy|,
  - signed weights are reported on the pipeline,
  - The F1-on-held-out threshold lands near the analytic boundary
    (annot threshold separates the classes exactly).

The pre-fix code would have produced different (broken) numbers — the
exact failure mode depending on which bug we re-introduced. This test
catches all three.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Add analysis/ted_logreg/ to sys.path so the test can import the
# module. ted_logreg.py lives outside bin/, so conftest's path setup
# doesn't reach it.
REPO_ROOT = Path(__file__).resolve().parent.parent
LOGREG_DIR = REPO_ROOT / "analysis" / "ted_logreg"
if str(LOGREG_DIR) not in sys.path:
    sys.path.insert(0, str(LOGREG_DIR))

import ted_logreg  # noqa: E402


# ── Synthetic fixture ─────────────────────────────────────────────────


@pytest.fixture(scope="module")
def synth_df() -> pd.DataFrame:
    """1000-row DataFrame mimicking the schema fit_end() expects.

    Layout:
      - joint_tp ∈ {0, 1} with 400 positives, 600 negatives (mild imbalance
        matching real PacBio data: ~1752 TP / 3461 FP).
      - TED_depth: positives drawn from N(0.7, 0.15), negatives from
        N(0.4, 0.15) — mildly separating (AUC ~ 0.85).
      - TED_tss_model and TED_tts_model: constant = 0.85 for ALL rows
        (zero variance). The logreg's coefficient on this feature must
        come out at exactly zero because there's no signal to fit.
      - TED_tss_annot and TED_tts_annot: positives drawn from N(0.9, 0.05),
        negatives from N(0.1, 0.05) — perfectly separating (AUC ~ 1.0).
    """
    rng = np.random.default_rng(seed=42)
    n_pos, n_neg = 400, 600
    n = n_pos + n_neg

    joint_tp = np.concatenate([np.ones(n_pos), np.zeros(n_neg)]).astype(np.int8)

    depth = np.concatenate([
        rng.normal(0.7, 0.15, n_pos),
        rng.normal(0.4, 0.15, n_neg),
    ]).clip(0.0, 1.0)

    model = np.full(n, 0.85)  # non-informative — constant across rows

    annot = np.concatenate([
        rng.normal(0.9, 0.05, n_pos),
        rng.normal(0.1, 0.05, n_neg),
    ]).clip(0.0, 1.0)

    # Shuffle so the train/test KFold splits don't trivially partition
    # the data into pure-positive / pure-negative folds.
    perm = rng.permutation(n)

    df = pd.DataFrame({
        "junc_id": [f"chr_test:{i}:+" for i in range(n)],
        "status": ["pass"] * n,
        "tss_pos": rng.integers(1000, 100000, size=n),
        "tts_pos": rng.integers(1000, 100000, size=n),
        "TED_depth": depth,
        "TED_tss_model": model,
        "TED_tts_model": model,
        "TED_tss_annot": annot,
        "TED_tts_annot": annot,
        "joint_tp": joint_tp,
    }).iloc[perm].reset_index(drop=True)
    return df


# ── Tests ─────────────────────────────────────────────────────────────


def test_fit_end_returns_expected_tuple_shape(synth_df):
    """fit_end must return (pipeline, weights, threshold, roc, pr_curve)."""
    out = ted_logreg.fit_end(synth_df, "tss")
    assert len(out) == 5
    pipe, deploy_weights, threshold, roc, pr_curve = out
    assert deploy_weights.shape == (3,)
    assert isinstance(threshold, float)
    assert 0.0 <= threshold <= 1.0
    assert len(roc) == 3  # (fpr, tpr, roc_thresh)
    assert len(pr_curve) == 3  # (prec_arr, rec_arr, thresh_arr)


def test_signed_weights_attached_to_pipeline(synth_df):
    """REGRESSION: the pre-fix code zeroed any negative coefficient. The
    fix exposes signed weights (|·|-normalised, sum |w|=1) on the
    pipeline so the user can see anti-predictive features as negative."""
    pipe, _, _, _, _ = ted_logreg.fit_end(synth_df, "tss")
    assert hasattr(pipe, "signed_weights_")
    sw = pipe.signed_weights_
    assert sw.shape == (3,)
    # Sum of magnitudes ~ 1 (since model_coef is exactly 0, only depth
    # and annot contribute, but the |·|-normalisation should still sum
    # to 1 over the surviving features).
    assert abs(np.abs(sw).sum() - 1.0) < 1e-6


def test_deploy_weights_sum_to_one_and_nonneg(synth_df):
    """Deploy weights must be non-negative and sum to 1 — TED's JSON
    config requires this. The clip-then-renormalise step is the
    deploy-time guard the fix retains."""
    _, w, _, _, _ = ted_logreg.fit_end(synth_df, "tss")
    assert (w >= 0).all()
    assert abs(w.sum() - 1.0) < 1e-6


def test_model_weight_is_zero_when_feature_is_constant(synth_df):
    """REGRESSION: a constant (non-informative) feature must produce a
    zero deploy weight. This is the *correct* behaviour of clip-to-zero
    when the raw coefficient is genuinely ~0. The pre-fix bug was that
    real-data features with marginally-negative coefficients ALSO got
    clipped to zero; that's tested separately via real-data validation,
    not here. This test pins the constant-feature edge case."""
    _, w_tss, _, _, _ = ted_logreg.fit_end(synth_df, "tss")
    _, w_tts, _, _, _ = ted_logreg.fit_end(synth_df, "tts")
    # Constant model feature → coefficient ~ 0 → deploy weight ~ 0
    assert w_tss[1] < 0.05, f"TSS model weight should be ~0, got {w_tss[1]}"
    assert w_tts[1] < 0.05, f"TTS model weight should be ~0, got {w_tts[1]}"


def test_weight_ordering_matches_feature_informativeness(synth_df):
    """Annot perfectly separates (AUC ~1), depth mildly separates
    (AUC ~0.85), model is non-informative (constant). After clipping
    negatives to zero, the deploy weights should reflect that ordering:
        w_annot > w_depth > w_model
    """
    for end in ("tss", "tts"):
        _, w, _, _, _ = ted_logreg.fit_end(synth_df, end)
        assert w[2] > w[0] > w[1], (
            f"{end}: expected w_annot ({w[2]:.3f}) > w_depth ({w[0]:.3f}) "
            f"> w_model ({w[1]:.3f}); not the case"
        )


def test_f1_threshold_in_reasonable_range(synth_df):
    """The F1-on-held-out threshold should land somewhere on the
    composite-score axis where the perfectly-separating annot feature
    cleanly partitions the labels. Because annot dominates the composite,
    the threshold ends up near the midpoint of the composite range
    (0.3 < t < 0.7). Pre-fix Youden's-J on the same data routinely
    pushed the threshold above 0.85, collapsing recall."""
    for end in ("tss", "tts"):
        _, _, t, _, _ = ted_logreg.fit_end(synth_df, end)
        assert 0.20 < t < 0.80, (
            f"{end}: F1 threshold {t:.2f} outside reasonable range "
            "[0.20, 0.80] — possible regression to a Youden's-J-style "
            "selector that over-penalises FPs at high TPR"
        )


def test_no_class_weight_balanced(synth_df):
    """REGRESSION: class_weight='balanced' over-corrected the imbalance
    on real data and pushed the model coefficient slightly negative.
    Confirm the LogisticRegression inside the Pipeline has no
    class_weight set."""
    pipe, _, _, _, _ = ted_logreg.fit_end(synth_df, "tss")
    lr = pipe.named_steps["lr"]
    # sklearn's default class_weight is None.
    assert lr.class_weight is None, (
        "LogisticRegression has class_weight set — the fix removed this "
        "because 'balanced' over-corrected the mild 1:2 imbalance and "
        "biased the model coefficient negative on real data."
    )


def test_fit_end_runs_on_minimal_imbalance(synth_df):
    """Sanity: fit_end accepts the mild 400/600 imbalance without
    raising. The pre-fix code did too — but with class_weight=balanced
    the result was biased. With class_weight=None the result is
    interpretable."""
    # Should not raise.
    pipe, w, t, roc, pr = ted_logreg.fit_end(synth_df, "tss")
    assert (w >= 0).all()
    assert abs(w.sum() - 1.0) < 1e-6


def test_fit_end_raises_on_no_variance():
    """If every row has the same joint_tp, the function must refuse to
    fit. Pinned because a future refactor that silently zeroes
    coefficients on this edge case would mask data-pipeline bugs."""
    df = pd.DataFrame({
        "junc_id": ["chr_test:0:+"] * 10,
        "TED_depth": [0.5] * 10,
        "TED_tss_model": [0.5] * 10,
        "TED_tss_annot": [0.5] * 10,
        "TED_tts_model": [0.5] * 10,
        "TED_tts_annot": [0.5] * 10,
        "joint_tp": [1] * 10,  # all positive — no negatives to learn from
    })
    with pytest.raises(ValueError, match="No variance"):
        ted_logreg.fit_end(df, "tss")


# ── R3 fix: one-to-one summit matching in ted_logreg ──────────────────────


def _summit_peaks(positions, strand="+", chrom="chr_test"):
    """Build a summit index (one summit per 1-bp peak) for ted_logreg."""
    rows = [(chrom, p, p + 1, strand) for p in positions]
    from evaluation.peak_matching import build_summit_index
    return build_summit_index(rows)


def test_assign_joint_tp_one_to_one_not_inflated():
    """Three candidate ends near ONE CAGE summit and ONE dRNA summit must
    yield exactly ONE joint TP, not three (the R3 recall inflator)."""
    cage = _summit_peaks([1000])   # single TSS summit
    drna = _summit_peaks([5000])   # single TTS summit
    df = pd.DataFrame({
        # all three rows fall within +-50 of the single summit pair
        "junc_id": ["chr_test:a:+", "chr_test:a:+", "chr_test:a:+"],
        "tss_pos": [990, 1000, 1010],
        "tts_pos": [4990, 5000, 5010],
    })
    out = ted_logreg.assign_joint_tp(df, cage, drna)
    assert out["joint_tp"].sum() == 1, (
        "three calls on one summit pair must collapse to a single joint TP"
    )


def test_assign_joint_tp_greedy_nearest_wins():
    """The nearest candidate claims the summit; farther ones do not."""
    cage = _summit_peaks([1000])
    drna = _summit_peaks([5000])
    df = pd.DataFrame({
        "junc_id": ["chr_test:a:+", "chr_test:a:+"],
        "tss_pos": [1002, 1040],   # 1002 nearer than 1040
        "tts_pos": [5001, 5045],
    })
    out = ted_logreg.assign_joint_tp(df, cage, drna)
    assert out.iloc[0]["joint_tp"] == 1
    assert out.iloc[1]["joint_tp"] == 0


def test_assign_joint_tp_separate_summits_both_count():
    """Two calls each near their OWN summit pair both count as TPs."""
    cage = _summit_peaks([1000, 2000])
    drna = _summit_peaks([5000, 6000])
    df = pd.DataFrame({
        "junc_id": ["chr_test:a:+", "chr_test:a:+"],
        "tss_pos": [1005, 2005],
        "tts_pos": [5005, 6005],
    })
    out = ted_logreg.assign_joint_tp(df, cage, drna)
    assert out["joint_tp"].sum() == 2


def test_dedup_pr_recall_capped_by_summits():
    """Recall numerator cannot exceed the number of distinct summits even
    when many survivors cluster on a single summit (the inflation guard)."""
    cage = _summit_peaks([1000])
    drna = _summit_peaks([5000])
    # 4 distinct (junc, tss, tts) survivors, all on the one summit pair.
    df = pd.DataFrame({
        "junc_id": ["chr_test:a:+"] * 4,
        "tss_pos": [985, 995, 1005, 1015],
        "tts_pos": [4985, 4995, 5005, 5015],
    })
    out = ted_logreg._dedup_pr_from_survivors(df, cage, drna)
    (prec_5, rec_5, f1_5, prec_3, rec_3, f1_3,
     paired_prec, paired_rec, paired_f1) = out
    # one summit -> recall 1.0 (1/1); precision 1/4 (one of four pairs hits)
    assert rec_5 == 1.0 and rec_3 == 1.0
    assert abs(prec_5 - 0.25) < 1e-9
    assert abs(prec_3 - 0.25) < 1e-9
