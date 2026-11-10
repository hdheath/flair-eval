#!/usr/bin/env python3
"""Regression tests for the load_ted_log dedup behavior in
ted/bin/plot_per_locus.py (sister repo).

The bug we fixed (May 8, 2026): the original load_ted_log dedup was
`drop_duplicates(subset=["junc_id", "cluster_id"], keep="last")`. In
per-SJC mode each cluster_id is 1:1 with a chain so this is correct,
but in --ted_global mode the same per-stage cluster_id (`c0`, `c1`)
gets reused across stages with DIFFERENT (tss_pos, tts_pos) coords —
each row is a stage-specific routing decision. Deduping on
(junc_id, cluster_id) alone destroys all but the last stage's
end-pair, silently dropping real isoforms from the plot.

Fix: include (tss_pos, tts_pos) in the dedup key so each unique
end-pair survives even when its c-label collides with another
end-pair's c-label in a different stage.

Fixture: tests/data/tiny_ted_log_global.tsv contains 4 rows on the
+:3j chain where cluster_id `c0` appears at TWO distinct
(tss_pos, tts_pos) values, and `c1` likewise. Pre-fix dedup
would collapse these 4 rows to 2 (one c0 + one c1). Post-fix
dedup keeps all 4.

This test imports load_ted_log from the SISTER ted/ repo (not
flair4_ted). The conftest fixture `ted_repo_bin` resolves the path
once. If the cross-repo layout changes, fix conftest first.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest


# ── Cross-repo import setup ────────────────────────────────────────


@pytest.fixture(scope="module")
def load_ted_log(ted_repo_bin):
    """Import load_ted_log from the sister ted/ repo's bin/.

    Raises pytest.skip if the sister repo isn't present (e.g. when
    running tests in CI without the ted/ checkout). Local development
    should always have both repos side-by-side.
    """
    if not (ted_repo_bin / "plot_per_locus.py").exists():
        pytest.skip(
            f"sister ted/ repo not found at {ted_repo_bin}; skipping "
            f"cross-repo loader tests"
        )
    if str(ted_repo_bin) not in sys.path:
        sys.path.insert(0, str(ted_repo_bin))
    # Reload-safe: importlib invalidates caches
    import importlib
    if "plot_per_locus" in sys.modules:
        importlib.reload(sys.modules["plot_per_locus"])
    import plot_per_locus
    return plot_per_locus.load_ted_log


# ── Tests ──────────────────────────────────────────────────────────


def test_load_returns_all_unique_quadruples(load_ted_log, tiny_ted_log_global):
    """REGRESSION: dedup key must be (junc_id, cluster_id, tss_pos,
    tts_pos). In the fixture, the +:3j chain has 4 rows with c0/c1
    each appearing twice at distinct (tss_pos, tts_pos). All 4 must
    survive. If anyone reverts the dedup key to (junc_id, cluster_id)
    alone, only 2 of those 4 rows survive — this assertion catches
    that."""
    df = load_ted_log(str(tiny_ted_log_global))

    # Restrict to the +:3j chain to isolate the c-label collision case.
    sub = df[df["junc_id"] == "chr_test:None-None:+:3j"]
    assert len(sub) == 4, (
        f"expected 4 rows on +:3j chain (c0 at 2 distinct end-pairs, "
        f"c1 at 2 distinct end-pairs), got {len(sub)}. "
        f"This usually means load_ted_log's dedup key was reverted to "
        f"only (junc_id, cluster_id), which destroys distinct end-pairs "
        f"that share a c-label across stages."
    )

    # Both c0 entries must be present with their distinct tss_pos values.
    c0 = sub[sub["cluster_id"] == "c0"]
    assert len(c0) == 2
    assert set(c0["tss_pos"].astype(int)) == {2000, 2400}

    c1 = sub[sub["cluster_id"] == "c1"]
    assert len(c1) == 2
    assert set(c1["tss_pos"].astype(int)) == {2000, 2400}


def test_load_preserves_global_peak_row(load_ted_log, tiny_ted_log_global):
    """load_ted_log itself must NOT filter by stage. The `global_peak`
    rows (partition-summary; tss_pos>=0 but tts_pos=-1 or vice versa)
    are filtered downstream (in the plotter, in pr_sweep), but at the
    loader layer they should remain. This separation lets each consumer
    apply its own filter — the plotter ignores them, pr_sweep drops
    them, and a future analysis could use them for partition-level
    summaries."""
    df = load_ted_log(str(tiny_ted_log_global))
    gp = df[df["stage"] == "global_peak"]
    assert len(gp) == 1, (
        f"expected 1 global_peak row to survive load_ted_log, got "
        f"{len(gp)}. load_ted_log shouldn't filter by stage; that's "
        f"each consumer's responsibility."
    )


def test_load_preserves_sentinel_rows(load_ted_log, tiny_ted_log_global):
    """Rows with tss_pos=-1 or tts_pos=-1 (sentinels for "TSS-only
    peak" / "TTS-only peak" / "noise") must survive load_ted_log
    intact. Filtering of these is downstream (sjc_alt_end_analysis,
    ted_pr_sweep, etc.). The fixture has 4 such rows."""
    df = load_ted_log(str(tiny_ted_log_global))
    sentinels = df[(df["tss_pos"] == -1) | (df["tts_pos"] == -1)]
    assert len(sentinels) == 4, (
        f"expected 4 sentinel rows (1 global_peak with tts=-1, 1 c0 "
        f"with tts=-1, 1 c1 with tss=-1, 1 noise with both -1), got "
        f"{len(sentinels)}"
    )


def test_load_total_row_count(load_ted_log, tiny_ted_log_global):
    """The fixture has 9 rows; all should be distinct on the
    (junc_id, cluster_id, tss_pos, tts_pos) dedup key, so the loader
    returns 9. If anyone tightens the dedup key further (e.g. drops
    tss_pos/tts_pos), this drops below 9 and the test catches it."""
    df = load_ted_log(str(tiny_ted_log_global))
    assert len(df) == 9, (
        f"expected 9 rows after dedup (all 9 are distinct on the "
        f"(junc_id, cluster_id, tss_pos, tts_pos) key); got {len(df)}"
    )


def test_load_casts_numeric_columns(load_ted_log, tiny_ted_log_global):
    """The loader casts tss_pos / tts_pos / n_reads / etc. to numeric
    dtype. Without this, the dedup key uses string comparison and
    rows like '2000' vs '2000.0' would not match. Sanity-check the
    types are numeric."""
    df = load_ted_log(str(tiny_ted_log_global))
    for col in ("n_reads", "jc_n_reads_total", "tss_pos", "tts_pos",
                "TED_tss_reality", "TED_tts_reality", "TED_confidence",
                "TED_depth"):
        if col in df.columns:
            # Either int or float dtype is fine; both numeric.
            assert df[col].dtype.kind in "if", (
                f"column {col} should be numeric (int/float); got "
                f"{df[col].dtype}"
            )
