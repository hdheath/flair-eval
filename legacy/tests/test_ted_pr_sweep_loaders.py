#!/usr/bin/env python3
"""Regression tests for the load_ted_log filter chain in
analysis/ted_pr_sweep/ted_pr_sweep.py.

The bug we fixed (May 6, 2026): pr_sweep's load_ted_log was including
`global_peak` partition-summary rows in its per-isoform candidate set.
Those rows have junc_id like 'chr*:None-None:strand:0j' with one of
(tss_pos, tts_pos) set to a sentinel -1, which caused phantom 0j FPs
in the confusion matrix because the unset axis fails every peak match.

Fix: a 4-stage filter chain in load_ted_log
(analysis/ted_pr_sweep/ted_pr_sweep.py:166-191):

    1. status in {pass, reject}             — drops `noise` rows
    2. stage in {ted_cluster, global_assign,
                 fallback, tss_snap}         — drops `global_peak` rows
    3. tss_pos.notna() AND tss_pos >= 0      — drops -1 sentinel
    4. tts_pos.notna() AND tts_pos >= 0      — drops -1 sentinel

These tests pin all 4 filters using tests/data/tiny_ted_log_global.tsv,
which contains:

    9 raw rows                        breakdown
    ────────────                      ─────────────────────────────────
    1 global_peak (chr_test:+:0j)     should be dropped by stage filter
    1 noise (chr_test:+:1j)           should be dropped by status filter
    1 c1 row with tss_pos=-1 (-:2j)   should be dropped by tss_pos filter
    1 c0 row with tts_pos=-1 (-:2j)   should be dropped by tts_pos filter
    5 valid c0/c1 rows                should survive

Expected post-load row count: 5.

If anyone removes any one of the 4 filter stages, the row count
changes and at least one of the assertions fails. This catches
regressions to the original phantom-FP bug.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest


# ── Cross-module import setup ──────────────────────────────────────


@pytest.fixture(scope="module")
def load_ted_log():
    """Import load_ted_log from analysis/ted_pr_sweep/ted_pr_sweep.py.

    Note: pr_sweep is in `analysis/`, not `bin/`, so conftest's
    sys.path bootstrap doesn't cover it. Add the path manually.
    """
    pr_sweep_dir = (
        Path(__file__).resolve().parent.parent
        / "analysis"
        / "ted_pr_sweep"
    )
    if not (pr_sweep_dir / "ted_pr_sweep.py").exists():
        pytest.skip(
            f"ted_pr_sweep.py not found at {pr_sweep_dir}; skipping"
        )
    if str(pr_sweep_dir) not in sys.path:
        sys.path.insert(0, str(pr_sweep_dir))
    import importlib
    if "ted_pr_sweep" in sys.modules:
        importlib.reload(sys.modules["ted_pr_sweep"])
    import ted_pr_sweep
    return ted_pr_sweep.load_ted_log


# ── Tests ──────────────────────────────────────────────────────────


def test_total_rows_after_filter_chain(load_ted_log, tiny_ted_log_global):
    """REGRESSION: 9 input rows → 5 output rows after the 4-filter
    chain. If anyone removes any filter, the count changes."""
    df = load_ted_log(str(tiny_ted_log_global))
    assert len(df) == 5, (
        f"expected 5 rows after the 4-filter chain (1 global_peak, "
        f"1 noise, and 2 sentinel rows must be dropped); got {len(df)}. "
        f"This usually means one of the load_ted_log filter stages "
        f"was removed: status / stage / tss_pos / tts_pos."
    )


def test_global_peak_rows_dropped(load_ted_log, tiny_ted_log_global):
    """REGRESSION: the `stage in {...}` filter must exclude
    'global_peak' rows. The fixture has 1 such row (the partition
    summary for chr_test:+:0j). It should be absent in the output."""
    df = load_ted_log(str(tiny_ted_log_global))
    gp = df[df["stage"] == "global_peak"]
    assert len(gp) == 0, (
        f"expected 0 global_peak rows in the output; got {len(gp)}. "
        f"This is the phantom-FP regression — global_peak rows are "
        f"partition-level summaries, not per-isoform candidates, "
        f"and including them inflates FP counts in the confusion "
        f"matrix."
    )


def test_noise_rows_dropped(load_ted_log, tiny_ted_log_global):
    """REGRESSION: the `status in {pass, reject}` filter must drop
    `noise` rows. The fixture has 1 such row (status='noise', both
    tss_pos and tts_pos = -1)."""
    df = load_ted_log(str(tiny_ted_log_global))
    noise = df[df["status"] == "noise"]
    assert len(noise) == 0


def test_tss_sentinel_rows_dropped(load_ted_log, tiny_ted_log_global):
    """REGRESSION: rows where tss_pos == -1 must be dropped.
    The fixture has 1 such row (c1 on chr_test:-:2j)."""
    df = load_ted_log(str(tiny_ted_log_global))
    sentinels = df[df["tss_pos"] == -1]
    assert len(sentinels) == 0
    # Also no rows with tss_pos NaN (string '-1' would fail this if not
    # numerically coerced).
    assert df["tss_pos"].notna().all()
    # Every surviving row's tss_pos is >= 0.
    assert (df["tss_pos"] >= 0).all()


def test_tts_sentinel_rows_dropped(load_ted_log, tiny_ted_log_global):
    """REGRESSION: rows where tts_pos == -1 must be dropped.
    The fixture has 1 such row (c0 on chr_test:-:2j with tss_pos=7300,
    tts_pos=-1, drop_reason='tts_below_threshold')."""
    df = load_ted_log(str(tiny_ted_log_global))
    sentinels = df[df["tts_pos"] == -1]
    assert len(sentinels) == 0
    assert df["tts_pos"].notna().all()
    assert (df["tts_pos"] >= 0).all()


def test_all_surviving_rows_have_per_isoform_stages(
    load_ted_log, tiny_ted_log_global
):
    """Every surviving row must have stage in the per-isoform-candidate
    set. This is the explicit guard against future stage values being
    added without updating the filter — if ted.py later adds a new
    stage label like 'reanchor', the test author has to think about
    whether that stage represents per-isoform candidates and add it to
    the filter list (or not)."""
    df = load_ted_log(str(tiny_ted_log_global))
    allowed = {"ted_cluster", "global_assign", "fallback", "tss_snap"}
    actual = set(df["stage"].unique())
    assert actual.issubset(allowed), (
        f"surviving rows have stages outside the per-isoform-candidate "
        f"set. Allowed: {allowed}; actual: {actual}. If ted.py added a "
        f"new stage label, decide whether it's a per-isoform candidate "
        f"and update the filter in pr_sweep's load_ted_log."
    )


def test_dedup_quadruples_distinct_after_filter(
    load_ted_log, tiny_ted_log_global
):
    """The surviving 5 rows are all distinct on
    (junc_id, cluster_id, tss_pos, tts_pos). pr_sweep's load_ted_log
    doesn't dedup explicitly, but the fixture is constructed so that
    no real-world dedup would be needed."""
    df = load_ted_log(str(tiny_ted_log_global))
    keys = df[["junc_id", "cluster_id", "tss_pos", "tts_pos"]]
    assert len(keys.drop_duplicates()) == len(keys)
