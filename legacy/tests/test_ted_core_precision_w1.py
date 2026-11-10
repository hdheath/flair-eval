#!/usr/bin/env python3
"""W1 regression: ted_core.tss_tts_metrics precision is no longer inflated.

Before this pass, `tss_tts_metrics` computed PRECISION as a per-isoform
bedtools-closest fraction (m / n_tx — the fraction of isoforms whose end is
within the window of ANY peak, many-to-one) while RECALL went through the
shared one-to-one summit matcher.  The two halves of the headline F1
therefore disagreed on what a true positive was, and the no-filter TED
emitter (densest end set) got a structurally inflated precision: a fat
cluster of called ends near one CAGE/dRNA summit counted every isoform as a
hit.

After the W1 fix both precision and recall come from the SAME one-to-one TP
(bed_utils.one_to_one_match_counts), so a cluster of N ends near one summit
contributes exactly ONE TP on both sides:
    precision = TP / n_calls     recall = TP / n_summits

This test builds a tiny isoforms BED with 5 isoforms whose 5' ends all sit
within +-50 bp of ONE CAGE summit and asserts:
  * OLD many-to-one precision would have been 5/5 = 1.0,
  * NEW one-to-one precision is 1/5 = 0.2 (TP=1),
  * precision and recall share the same numerator.

Requires bedtools on PATH (the eval harness dependency).
"""
from __future__ import annotations

import shutil

import pytest

from evaluation.ted_core import tss_tts_metrics

pytestmark = pytest.mark.skipif(
    shutil.which("bedtools") is None, reason="bedtools required for tss_tts_metrics"
)


def _write_bed(path, rows):
    with open(path, "w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def _make_clustered_isoforms(path):
    """5 single-exon + isoforms whose 5' (start) ends cluster near 1000.

    BED12.  + strand, so TSS = start.  Starts at 980/990/1000/1010/1020 all
    fall within +-50 bp of the CAGE summit at 1000.  3' ends are spread far
    apart so they do not collide on the (single, distant) dRNA summit.
    """
    rows = []
    starts = [980, 990, 1000, 1010, 1020]
    for i, s in enumerate(starts):
        end = 4000 + i * 500   # distinct, far 3' ends
        rows.append([
            "chr_test", s, end, f"ISO{i}_GENE", 100, "+",
            s, end, "0,0,0", 1, end - s, 0,
        ])
    _write_bed(path, rows)


def test_w1_precision_not_inflated_by_cluster(tmp_path):
    iso_bed = tmp_path / "clustered.isoforms.bed"
    _make_clustered_isoforms(iso_bed)

    cage = tmp_path / "cage.bed"
    _write_bed(cage, [["chr_test", 950, 1050, "cage1", 100, "+"]])  # summit 1000

    peaks = {"prime5": cage, "prime3": None,
             "ref_prime5": None, "ref_prime3": None}

    metrics, _ = tss_tts_metrics(iso_bed, peaks, window=50, n_isoforms=5)

    # One-to-one TP = 1 across the 5-call cluster.
    # precision = TP / n_calls = 1/5;  recall = TP / n_summits = 1/1.
    assert metrics["5prime_precision"] == pytest.approx(1 / 5), (
        "precision must be the one-to-one 1/5, NOT the inflated per-isoform 5/5"
    )
    assert metrics["5prime_recall"] == pytest.approx(1.0)
    # The OLD many-to-one precision would have been 1.0 (all 5 isoforms near
    # the peak).  Guard against a regression back to it.
    assert metrics["5prime_precision"] < 0.99


def test_w1_precision_recall_share_tp(tmp_path):
    """Two isoforms on two distinct summits -> precision == recall == 1.0,
    and three isoforms (one extra on a shared summit) drop precision but not
    recall — proving the shared-TP accounting."""
    iso_bed = tmp_path / "two.isoforms.bed"
    # iso A start 1000 (summit1), iso B start 2000 (summit2),
    # iso C start 1005 (also near summit1 -> competes, loses).
    rows = []
    for i, s in enumerate((1000, 2000, 1005)):
        rows.append(["chr_test", s, 3000 + i * 400, f"ISO{i}", 100, "+",
                     s, 3000 + i * 400, "0,0,0", 1, 3000 + i * 400 - s, 0])
    _write_bed(iso_bed, rows)

    cage = tmp_path / "cage.bed"
    _write_bed(cage, [
        ["chr_test", 950, 1050, "c1", 100, "+"],     # summit 1000
        ["chr_test", 1950, 2050, "c2", 100, "+"],    # summit 2000
    ])
    peaks = {"prime5": cage, "prime3": None,
             "ref_prime5": None, "ref_prime3": None}
    metrics, _ = tss_tts_metrics(iso_bed, peaks, window=50, n_isoforms=3)

    # TP = 2 (both summits), n_calls = 3, n_summits = 2.
    assert metrics["5prime_precision"] == pytest.approx(2 / 3)
    assert metrics["5prime_recall"] == pytest.approx(2 / 2)


def test_w1_signal_summit_threading(tmp_path):
    """signal_tracks routes a dominant-signal summit through tss_tts_metrics."""
    iso_bed = tmp_path / "sig.isoforms.bed"
    # Single isoform whose 5' end (1075) is 75 bp from the peak CENTRE (1000,
    # a miss at W=50) but 5 bp from the SIGNAL summit (1080, a hit).
    _write_bed(iso_bed, [
        ["chr_test", 1075, 3000, "ISO0", 100, "+",
         1075, 3000, "0,0,0", 1, 1925, 0],
    ])
    cage = tmp_path / "cage.bed"
    _write_bed(cage, [["chr_test", 900, 1100, "c1", 100, "+"]])  # centre 1000
    peaks = {"prime5": cage, "prime3": None,
             "ref_prime5": None, "ref_prime3": None}
    signal = {"prime5": {("chr_test", "+"): [(910, 1.0), (1080, 9.0)]}}

    m_sig, _ = tss_tts_metrics(iso_bed, peaks, window=50, n_isoforms=1,
                               signal_tracks=signal)
    m_nosig, _ = tss_tts_metrics(iso_bed, peaks, window=50, n_isoforms=1)

    assert m_sig["5prime_precision"] == pytest.approx(1.0)    # signal summit hit
    assert m_nosig["5prime_precision"] == pytest.approx(0.0)  # centre miss
