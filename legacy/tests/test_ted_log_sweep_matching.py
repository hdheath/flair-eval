#!/usr/bin/env python3
"""W2 regression: the weight-selection sweep uses one-to-one summit matching.

`ted_log_sweep.py` computes the F1 that SELECTS the per-tech logreg weights
and thresholds.  Before this pass it matched candidate ends to truth peaks by
interval membership (`_hit_peak`) and counted recall by unique interval id —
many-to-one, so a fat cluster of surviving candidates near one CAGE/dRNA peak
inflated the selection F1 and biased the fitted weights toward a leaky
objective (W2 in the 2026-05-28 graft/eval review).

After the fix the sweep collapses each peak to ONE summit and assigns
survivors one-to-one (greedy nearest), so a cluster of N survivors on one
summit yields exactly ONE true positive — the same TP the headline F1 uses.

conftest.py puts bin/ on sys.path so `from evaluation.ted_log_sweep ...`
resolves.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from evaluation.ted_log_sweep import (
    _parse_peaks,
    _assign_peak_ids,
    _dedup_pr_one_to_one,
    _one_to_one_tp,
    _make_arrays,
)
from evaluation.peak_matching import build_summit_index, count_summits


def _row(tss_pos, tts_pos, jid):
    return {"junc_id": jid, "tss_pos": tss_pos, "tts_pos": tts_pos,
            "TED_depth": 1.0, "TED_tss_model": 1.0, "TED_tss_annot": 1.0,
            "TED_tts_model": 1.0, "TED_tts_annot": 1.0, "status": "pass"}


def test_parse_peaks_collapses_to_summits(tmp_path):
    bed = tmp_path / "cage.bed"
    bed.write_text(
        "chr_test\t950\t1050\tc1\t100\t+\n"
        "chr_test\t1950\t2050\tc2\t100\t+\n"
    )
    idx = _parse_peaks(str(bed))
    assert count_summits(idx) == 2
    assert idx[("chr_test", "+")] == [1000, 2000]


def test_one_to_one_tp_cluster_yields_one():
    # 5 calls all near the SAME single summit (1000) -> exactly one TP.
    summit_index = {("chr_test", "+"): [1000]}
    chrom = np.array(["chr_test"] * 5)
    strand = np.array(["+"] * 5)
    pos = np.array([980, 990, 1000, 1010, 1020], dtype=np.int64)
    assert _one_to_one_tp(chrom, strand, pos, summit_index) == 1


def test_one_to_one_tp_distinct_summits():
    # 3 calls each within W of a distinct summit; a 4th is 51 bp out (miss).
    summit_index = {("chr_test", "+"): [1000, 2000, 3000]}
    chrom = np.array(["chr_test"] * 4)
    strand = np.array(["+"] * 4)
    pos = np.array([1005, 2005, 3005, 4051], dtype=np.int64)
    assert _one_to_one_tp(chrom, strand, pos, summit_index) == 3


def test_one_to_one_tp_fallback_to_second_nearest():
    """FIX 3 regression: a call whose nearest summit is stolen by a closer
    call falls back to its next in-window summit.

    summits=[1000,1040], calls=[1030,1050], W=50.
      call 1050 is 10 bp from 1040 (nearest); call 1030 is 30 bp from 1040
      and 30 bp from 1000.  The old nearest-only precompute recorded BOTH
      calls' nearest as summit 1040, so deduping uids gave TP=1.  The shared
      matcher assigns 1050->1040 and 1030->1000, giving the correct TP=2.
    """
    summit_index = {("chr_test", "+"): [1000, 1040]}
    chrom = np.array(["chr_test"] * 2)
    strand = np.array(["+"] * 2)
    pos = np.array([1030, 1050], dtype=np.int64)
    assert _one_to_one_tp(chrom, strand, pos, summit_index) == 2


def test_sweep_precision_not_inflated_by_cluster(tmp_path):
    """THE W2 headline: 5 survivors near ONE CAGE summit -> TP=1, not 5."""
    cage = tmp_path / "cage.bed"
    cage.write_text("chr_test\t950\t1050\tc1\t100\t+\n")     # summit 1000
    drna = tmp_path / "drna.bed"
    drna.write_text("chr_test\t3450\t3550\td1\t100\t+\n")    # summit 3500

    cage_idx = _parse_peaks(str(cage))
    drna_idx = _parse_peaks(str(drna))

    # 5 candidate clusters, each a distinct junction chain, all with a TSS
    # within +-50 of summit 1000 and a (shared) far-away TTS.
    df = pd.DataFrame([
        _row(980,  9000, "chr_test:None-None:+:1j"),
        _row(990,  9000, "chr_test:None-None:+:2j"),
        _row(1000, 9000, "chr_test:None-None:+:3j"),
        _row(1010, 9000, "chr_test:None-None:+:4j"),
        _row(1020, 9000, "chr_test:None-None:+:5j"),
    ])
    df = _assign_peak_ids(df, cage_idx, drna_idx)

    (chrom_arr, strand_arr, tss_pos, tts_pos,
     cage_peaks, drna_peaks, total_cage, total_drna) = _make_arrays(df)
    assert total_cage == 1
    mask = np.ones(len(df), dtype=bool)   # all survive
    p5, r5, f5, p3, r3, f3 = _dedup_pr_one_to_one(
        mask, chrom_arr, strand_arr, tss_pos, tts_pos,
        cage_peaks, drna_peaks, total_cage, total_drna)

    # One-to-one: TP_5 = 1 over 5 survivors -> precision 1/5, recall 1/1.
    assert p5 == 1 / 5, "precision must NOT be the inflated 5/5"
    assert r5 == 1.0


def test_sweep_signal_summit(tmp_path):
    """Signal bedGraph moves the CAGE summit; matching follows it."""
    cage = tmp_path / "cage.bed"
    cage.write_text("chr_test\t900\t1100\tc1\t100\t+\n")    # centre 1000
    # Signal dominant at 1080.
    signal = {("chr_test", "+"): [(910, 1.0), (1080, 9.0)]}
    idx_sig = _parse_peaks(str(cage), signal=signal)
    idx_centre = _parse_peaks(str(cage))
    assert idx_sig[("chr_test", "+")] == [1080]
    assert idx_centre[("chr_test", "+")] == [1000]

    drna = tmp_path / "drna.bed"
    drna.write_text("chr_test\t3450\t3550\td1\t100\t+\n")
    drna_idx = _parse_peaks(str(drna))

    # Call at 1075: 75 bp from centre (miss) but 5 bp from signal summit (hit).
    df = pd.DataFrame([_row(1075, 9000, "chr_test:None-None:+:1j")])

    df_sig = _assign_peak_ids(df.copy(), idx_sig, drna_idx)
    df_centre = _assign_peak_ids(df.copy(), idx_centre, drna_idx)
    mask = np.ones(1, dtype=bool)

    a = _make_arrays(df_sig)
    b = _make_arrays(df_centre)
    # _make_arrays → (chrom, strand, tss_pos, tts_pos, cage_peaks, drna_peaks,
    #                 total_cage, total_drna)
    p5_sig = _dedup_pr_one_to_one(mask, *a)[0]
    p5_centre = _dedup_pr_one_to_one(mask, *b)[0]
    assert p5_sig == 1.0
    assert p5_centre == 0.0
