#!/usr/bin/env python3
"""Integration tests for the GLOBAL-only TED graft on flair-ted-52826.

Asserts the two locked guarantees of the port:

1. NO-FILTER EMIT (hdheath locked decision #4): every clustered isoform is
   emitted with its confidences; nothing is dropped on the reality/confidence
   threshold. Concretely: the number of candidate isoforms written to the
   firstpass intermediates equals the number of TED clusters (global_assign
   passes) in the ted_log — i.e. the no-filter redundant-end bypass does not
   collapse ends.

2. FINITE PER-END CONFIDENCES: the 8 strand-aware FlairBed columns
   (TSS/TTS × depth/model/annot/reality) are populated with finite values in
   [0, 1] (annot/reality) for every emitted isoform on the --ted path.

The test drives the real `flair transcriptome --ted` pipeline through the
dedicated conda env (flair-ted-52826) so it exercises the actual graft, not a
re-implementation. It is skipped if the env / models / clone are unavailable.
"""
from __future__ import annotations

import os
import subprocess
import tempfile
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

_missing = [str(p) for p in (FLAIR_BIN, CLONE, DATA, TSS_MODEL, TTS_MODEL) if not p.exists()]
pytestmark = pytest.mark.skipif(
    bool(_missing), reason=f"TED graft prerequisites missing: {_missing}")


def _run_ted(outdir: Path, extra_args=()):
    """Run `flair transcriptome --ted` on the tiny fixture; return (outprefix)."""
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
        "--ted", "--tss_model", str(TSS_MODEL), "--tts_model", str(TTS_MODEL),
        "--ted_log", "--keep_intermediate",
        *extra_args,
    ]
    res = subprocess.run(cmd, env=env, cwd=str(outdir),
                         capture_output=True, text=True)
    assert res.returncode == 0, (
        f"flair transcriptome --ted exited {res.returncode}\n"
        f"STDERR tail:\n{res.stderr[-2000:]}")
    return out


def _read_rows(path: Path):
    return [ln.rstrip("\n").split("\t")
            for ln in path.read_text().splitlines() if ln and not ln.startswith("#")]


@pytest.fixture(scope="module")
def ted_run():
    with tempfile.TemporaryDirectory() as td:
        out = _run_ted(Path(td))
        # Snapshot everything we need before the tmpdir is cleaned up.
        bed = _read_rows(Path(str(out) + ".isoforms.bed"))
        ted_log = _read_rows(Path(str(out) + ".ted_log.tsv"))
        reallyunfilt = _read_rows(Path(str(out) + ".firstpass.reallyunfiltered.bed"))
        unfilt = _read_rows(Path(str(out) + ".firstpass.unfiltered.bed"))
        gtf = Path(str(out) + ".isoforms.gtf").read_text()
    return {"bed": bed, "ted_log": ted_log,
            "reallyunfilt": reallyunfilt, "unfilt": unfilt, "gtf": gtf}


def test_outputs_produced(ted_run):
    assert len(ted_run["bed"]) > 0, "no isoforms emitted on --ted path"


def test_no_filter_emit_clusters_not_dropped(ted_run):
    """No-filter: the redundant-end bypass must not drop any cluster.

    #candidate isoforms (firstpass.unfiltered) == #raw clusters
    (firstpass.reallyunfiltered) == #global_assign passes in the ted_log.
    """
    n_clusters_log = sum(1 for r in ted_run["ted_log"]
                         if len(r) > 3 and r[2] == "global_assign" and r[3] == "pass")
    n_reallyunfilt = len(ted_run["reallyunfilt"])
    n_unfilt = len(ted_run["unfilt"])
    assert n_clusters_log > 0, "expected at least one TED cluster in ted_log"
    assert n_reallyunfilt == n_clusters_log, (
        f"raw clusters ({n_reallyunfilt}) != ted_log global_assign passes ({n_clusters_log})")
    # The no-filter redundant/confidence bypass must not collapse ends:
    assert n_unfilt == n_reallyunfilt, (
        f"no-filter dropped ends: unfiltered ({n_unfilt}) != raw clusters ({n_reallyunfilt})")


def test_global_peaks_not_threshold_dropped_silently(ted_run):
    """Any peak/cluster that fails a threshold is *logged*, never silently gone.

    With no-filter emit there should be no reject rows that vanish without a
    diagnostic; reality is recorded for every global_peak row.
    """
    peak_rows = [r for r in ted_run["ted_log"]
                 if len(r) > 3 and r[2] == "global_peak"]
    assert peak_rows, "no global_peak rows in ted_log"
    # TED_confidence column (index 10 in the ted_log header) is always written.
    for r in peak_rows:
        assert len(r) >= 11, f"ted_log row too short: {r}"


# Column layout of the FlairBed: bed12 (0-11) + 10 flair cols (12-21) + 8 TED
# cols (22-29): TSS_depth, TSS_model, TSS_annot, TSS_reality, TTS_depth,
# TTS_model, TTS_annot, TTS_reality.
_TED_COL_IDX = {
    "TSS_depth": 22, "TSS_model": 23, "TSS_annot": 24, "TSS_reality": 25,
    "TTS_depth": 26, "TTS_model": 27, "TTS_annot": 28, "TTS_reality": 29,
}


def test_eight_ted_columns_present_and_finite(ted_run):
    import math
    rows = ted_run["bed"]
    assert rows, "no isoforms to check columns on"
    for r in rows:
        assert len(r) == 30, f"expected 30 BED columns, got {len(r)}: {r[:6]}"
        for name, idx in _TED_COL_IDX.items():
            val = r[idx]
            assert val != "", f"{name} empty for isoform {r[3]} on --ted path"
            f = float(val)
            assert math.isfinite(f), f"{name}={val} not finite for {r[3]}"
            # depth/annot/reality are scores in [0,1]; model is a probability.
            assert 0.0 <= f <= 1.0, f"{name}={f} out of [0,1] for {r[3]}"


def test_ted_columns_flow_to_gtf(ted_run):
    gtf = ted_run["gtf"]
    for name in _TED_COL_IDX:
        assert name in gtf, f"{name} not found as a GTF attribute"


def test_snap_drop_fix_threshold_does_not_drop_reads():
    """FIX 2 (snap-drop / no-filter): a reality threshold of 1.0 rejects EVERY
    global peak, but under no-filter NOTHING is dropped on that confidence
    threshold — reads still snap to their nearest peak and emit, carrying the
    (low) reality forward as a column.

    Pre-fix, only threshold-approved peaks were snap targets, so an all-reject
    threshold dropped every read -> 0 isoforms. Post-fix the # emitted isoforms
    still equals the # clusters and the per-end confidences are all present.
    """
    with tempfile.TemporaryDirectory() as td:
        out = _run_ted(Path(td), extra_args=("--ted_threshold", "1.0"))
        bed = _read_rows(Path(str(out) + ".isoforms.bed"))
        ted_log = _read_rows(Path(str(out) + ".ted_log.tsv"))
        reallyunfilt = _read_rows(Path(str(out) + ".firstpass.reallyunfiltered.bed"))
        unfilt = _read_rows(Path(str(out) + ".firstpass.unfiltered.bed"))

    # Every global peak should be rejected at threshold=1.0 ...
    peak_status = [r[3] for r in ted_log
                   if len(r) > 3 and r[2] == "global_peak"]
    assert peak_status, "no global_peak rows in ted_log"
    assert all(s == "reject" for s in peak_status), (
        f"expected ALL peaks rejected at threshold=1.0, got {set(peak_status)}")

    # ... yet reads still snap and clusters still emit (no confidence drop).
    n_clusters_log = sum(1 for r in ted_log
                         if len(r) > 3 and r[2] == "global_assign" and r[3] == "pass")
    assert n_clusters_log > 0, (
        "FIX 2 regression: all reads dropped at snap time under all-reject "
        "threshold (confidence-threshold drop survived the no-filter fix)")
    # #emitted isoforms == #clusters: no-filter parity holds even all-reject.
    assert len(reallyunfilt) == n_clusters_log, (
        f"raw clusters ({len(reallyunfilt)}) != global_assign passes ({n_clusters_log})")
    assert len(unfilt) == len(reallyunfilt), (
        f"no-filter dropped ends: unfiltered ({len(unfilt)}) != raw ({len(reallyunfilt)})")
    assert len(bed) > 0, "no isoforms emitted under all-reject threshold (snap-drop fix failed)"

    # Per-end confidences all present and finite for every emitted isoform —
    # carrying the (low) reality as a column is the whole point of the fix.
    import math
    for r in bed:
        assert len(r) == 30, f"expected 30 BED cols, got {len(r)}"
        for name, idx in _TED_COL_IDX.items():
            v = r[idx]
            assert v != "", f"{name} empty for {r[3]} under all-reject threshold"
            assert math.isfinite(float(v)), f"{name}={v} not finite for {r[3]}"


def test_chemistry_preset_runs_clean():
    """FIX 1: the --ted_chemistry preset bundle runs end-to-end (cDNA -> the
    sigma_tts=15 default fills in) and preserves no-filter parity."""
    with tempfile.TemporaryDirectory() as td:
        out = _run_ted(Path(td), extra_args=("--ted_chemistry", "ont_cdna"))
        bed = _read_rows(Path(str(out) + ".isoforms.bed"))
        ted_log = _read_rows(Path(str(out) + ".ted_log.tsv"))
    assert bed, "no isoforms emitted under --ted_chemistry ont_cdna"
    n_clusters = sum(1 for r in ted_log
                     if len(r) > 3 and r[2] == "global_assign" and r[3] == "pass")
    assert n_clusters > 0
    # cDNA preset must NOT trust the 3' (only dRNA does); TTS_model stays real.
    tts_models = {float(r[_TED_COL_IDX["TTS_model"]]) for r in bed}
    assert any(v != 1.0 for v in tts_models), (
        "cDNA --ted_chemistry should NOT short-circuit TTS_model to 1.0")


def test_chemistry_drna_preset_trusts_3prime():
    """FIX 1: --ted_chemistry ont_drna enables 3'-trust via the bundle (no
    explicit --ted_trust_3prime needed) -> TTS_model forced to 1.0."""
    with tempfile.TemporaryDirectory() as td:
        out = _run_ted(Path(td), extra_args=("--ted_chemistry", "ont_drna"))
        bed = _read_rows(Path(str(out) + ".isoforms.bed"))
    assert bed, "no isoforms emitted under --ted_chemistry ont_drna"
    for r in bed:
        assert float(r[_TED_COL_IDX["TTS_model"]]) == 1.0, (
            f"ont_drna preset should trust 3' (TTS_model==1.0): {r[3]}")


def test_dRNA_trust_short_circuits_tts_model():
    """--ted_trust_3prime forces TTS_model to 1.0 while keeping TSS_model real."""
    with tempfile.TemporaryDirectory() as td:
        out = _run_ted(Path(td), extra_args=("--ted_trust_3prime",))
        rows = _read_rows(Path(str(out) + ".isoforms.bed"))
    assert rows, "no isoforms emitted under --ted_trust_3prime"
    for r in rows:
        assert float(r[_TED_COL_IDX["TTS_model"]]) == 1.0, (
            f"TTS_model not trusted (==1.0) under --ted_trust_3prime: {r[3]}")
    # At least one isoform should still have a non-trivial (real) TSS model score.
    tss_models = {float(r[_TED_COL_IDX["TSS_model"]]) for r in rows}
    assert any(v != 1.0 for v in tss_models), (
        "TSS_model should remain real (not all 1.0) under --ted_trust_3prime")


# ─── FIX 2b: distance-miss reads are KEPT at their own end position ───────
# Unit-level (no subprocess): drives ted_global_assign_reads directly so it
# runs even when the full FLAIR env / models are unavailable. Asserts that a
# read whose end is FAR from any approved peak is still emitted (not dropped)
# at its own end position, carrying a finite TED reality/confidence column.
import sys as _sys
from collections import namedtuple as _namedtuple

_CLONE_SRC = CLONE / "src"
if str(_CLONE_SRC) not in _sys.path:
    _sys.path.insert(0, str(_CLONE_SRC))

try:
    import numpy as _np
    from flair import ted as _ted
    _TED_IMPORTABLE = True
except Exception:  # pragma: no cover - clone/deps unavailable
    _TED_IMPORTABLE = False

_RR = _namedtuple("_RR", "chrom strand start end name polyA")
_JUNC = _namedtuple("_JUNC", "start end")   # duck-types isoform_data.Junc
_VK = dict(w_depth=0.4, w_model=0.35, w_annot=0.25, threshold=0.5,
           depth_saturation=50.0, annot_half_life=100.0, annot_max_dist=500,
           trust_3prime=False)


def _emitted_exons_valid(rep_start, rep_end, juncs):
    """Mirror isoform_data.get_exons geometry and return True iff every exon
    has start <= end (i.e. would NOT trip the Exon() assert). rep_start/rep_end
    are the snapped genomic ends; juncs is the chain (genomic order)."""
    if not juncs:
        return rep_start <= rep_end
    exons = [(rep_start, juncs[0].start)]
    for i in range(len(juncs) - 1):
        exons.append((juncs[i].end, juncs[i + 1].start))
    exons.append((juncs[-1].end, rep_end))
    return all(s <= e for s, e in exons)


@pytest.mark.skipif(not _TED_IMPORTABLE, reason="flair.ted not importable")
def test_distance_miss_read_kept_at_own_position_not_dropped():
    """FIX 2b: a read whose 5'/3' end is far from ANY approved peak is NOT
    dropped at snap time. It is emitted at its own end position carrying a
    finite TED reality/confidence column (a logged diagnostic, not a removal).
    Normal snapping still applies when a peak IS within the window."""
    import math
    approved = {
        ("chr1", "+"): {
            "tss": _np.array([1000], dtype=_np.int64),
            "tts": _np.array([5000], dtype=_np.int64),
            "tss_scores": {1000: {"depth_score": 0.5, "model_score": 0.6,
                                  "annotation_score": 0.4, "annot_dist_bp": 10}},
            "tts_scores": {5000: {"depth_score": 0.5, "model_score": 0.6,
                                  "annotation_score": 0.4, "annot_dist_bp": 10}},
            "tss_reality": {1000: 0.5}, "tts_reality": {5000: 0.5},
            "tss_approved": _np.array([1000], dtype=_np.int64),
            "tts_approved": _np.array([5000], dtype=_np.int64),
        }
    }
    reads = [
        _RR("chr1", "+", 1005, 5002, "near", (False, False)),    # snaps to peaks
        _RR("chr1", "+", 20000, 25000, "far", (False, False)),   # far from any peak
    ]
    out = _ted.ted_global_assign_reads(
        reads, "chr1", approved, tss_snap=50, tts_snap=25, scorer=None,
        val_kwargs_tss=_VK, val_kwargs_tts=_VK,
        threshold_tss=0.5, threshold_tts=0.5)

    # No drop: #emitted clusters == #distinct (snapped/self) end pairs == 2.
    assert len(out) == 2, (
        f"distance-miss read was DROPPED at snap time: got {len(out)} groups, "
        "expected 2 (snap-drop fix regression)")

    by_read = {}
    for conf, rs, re_, strand, repname, names, scores, tag in out:
        for nm in names:
            by_read[nm] = (conf, rs, re_, scores)
    assert "far" in by_read, "far read not emitted"
    assert "near" in by_read, "near read not emitted"

    # The near read snapped to the approved peaks; the far read kept its own
    # end position (no rounding to the distant peak).
    near_conf, near_s, near_e, _ = by_read["near"]
    far_conf, far_s, far_e, far_scores = by_read["far"]
    assert (near_s, near_e) == (1000, 5000), "near read did not snap to peaks"
    assert (far_s, far_e) == (20000, 25000), (
        "far read was snapped away from its own position instead of kept")

    # Finite confidence carried as a column for the kept self-position read.
    assert math.isfinite(far_conf), "far read confidence not finite"
    for k in ("TED_tss_reality", "TED_tts_reality", "TED_confidence"):
        assert math.isfinite(float(far_scores[k])), f"{k} not finite for far read"


# ─── SPLICE-BOUNDS CLAMP: snapped/kept ends never cross splice structure ──
# Regression for the chr22 WTC11 AssertionError (isoform_data.py:68
# Exon(start<=end) via get_exons -> Exon(start, juncs[0].start)). A snapped or
# kept biological 5'/3' end that lands on the wrong side of a splice site must
# be CLAMPED back to a valid coordinate — NOT dropped — so get_exons always
# receives valid geometry. Covers BOTH strands and BOTH ends.

@pytest.mark.skipif(not _TED_IMPORTABLE, reason="flair.ted not importable")
def test_plus_strand_tss_peak_past_first_donor_is_clamped_not_dropped():
    """+ strand: the nearest approved TSS peak sits DOWNSTREAM of the first
    splice donor. Snapping there would make rep_start >= juncs[0].start and
    blow up get_exons. The read must be emitted (no drop) with rep_start
    clamped strictly inside the first exon, valid geometry, finite confidence."""
    import math
    # Chain: exon1 [800, 1200), intron, exon2 [1500, 5200). First donor=1200.
    juncs = [_JUNC(1200, 1500)]
    # TSS peak at 1400 is PAST the donor (1200); TTS peak at 5000 is fine.
    approved = {
        ("chr1", "+"): {
            "tss": _np.array([1400], dtype=_np.int64),
            "tts": _np.array([5000], dtype=_np.int64),
            "tss_scores": {1400: {"depth_score": 0.5, "model_score": 0.6,
                                  "annotation_score": 0.4, "annot_dist_bp": 10}},
            "tts_scores": {5000: {"depth_score": 0.5, "model_score": 0.6,
                                  "annotation_score": 0.4, "annot_dist_bp": 10}},
            "tss_reality": {1400: 0.5}, "tts_reality": {5000: 0.5},
            "tss_approved": _np.array([1400], dtype=_np.int64),
            "tts_approved": _np.array([5000], dtype=_np.int64),
        }
    }
    # Read starts at 1390 (its own TSS already past the donor), end 5002.
    reads = [_RR("chr1", "+", 1390, 5002, "overshoot", (False, False))]
    out = _ted.ted_global_assign_reads(
        reads, "chr1", approved, tss_snap=50, tts_snap=25, scorer=None,
        val_kwargs_tss=_VK, val_kwargs_tts=_VK,
        threshold_tss=0.5, threshold_tts=0.5, chain_juncs=juncs)

    assert len(out) == 1, f"read was DROPPED instead of clamped: {len(out)} groups"
    conf, rep_start, rep_end, strand, repname, names, scores, tag = out[0]
    # rep_start must be strictly before the first donor (valid first exon).
    assert rep_start < juncs[0].start, (
        f"rep_start {rep_start} not clamped before first donor {juncs[0].start}")
    assert rep_end >= juncs[-1].end, f"rep_end {rep_end} crosses last acceptor"
    assert _emitted_exons_valid(rep_start, rep_end, juncs), (
        f"emitted exon geometry invalid: start={rep_start} end={rep_end} juncs={juncs}")
    assert math.isfinite(conf) and math.isfinite(float(scores["TED_confidence"]))


@pytest.mark.skipif(not _TED_IMPORTABLE, reason="flair.ted not importable")
def test_plus_strand_tts_peak_before_last_acceptor_is_clamped_not_dropped():
    """+ strand symmetric case: nearest approved TTS peak sits UPSTREAM of the
    last splice acceptor. Must clamp rep_end strictly after the last acceptor."""
    # Chain: exon1 [800,1200) intron exon2 [1500,5200). Last acceptor=1500.
    juncs = [_JUNC(1200, 1500)]
    # TTS peak at 1450 is BEFORE the last acceptor (1500); TSS peak 900 ok.
    approved = {
        ("chr1", "+"): {
            "tss": _np.array([900], dtype=_np.int64),
            "tts": _np.array([1450], dtype=_np.int64),
            "tss_scores": {900: {"depth_score": 0.5, "model_score": 0.6,
                                 "annotation_score": 0.4, "annot_dist_bp": 10}},
            "tts_scores": {1450: {"depth_score": 0.5, "model_score": 0.6,
                                  "annotation_score": 0.4, "annot_dist_bp": 10}},
            "tss_reality": {900: 0.5}, "tts_reality": {1450: 0.5},
            "tss_approved": _np.array([900], dtype=_np.int64),
            "tts_approved": _np.array([1450], dtype=_np.int64),
        }
    }
    reads = [_RR("chr1", "+", 905, 1460, "shorttts", (False, False))]
    out = _ted.ted_global_assign_reads(
        reads, "chr1", approved, tss_snap=50, tts_snap=25, scorer=None,
        val_kwargs_tss=_VK, val_kwargs_tts=_VK,
        threshold_tss=0.5, threshold_tts=0.5, chain_juncs=juncs)
    assert len(out) == 1, "read dropped instead of clamped (TTS case)"
    _c, rep_start, rep_end, _s, _rn, _ns, _sc, _t = out[0]
    assert rep_end > juncs[-1].end, (
        f"rep_end {rep_end} not clamped after last acceptor {juncs[-1].end}")
    assert _emitted_exons_valid(rep_start, rep_end, juncs)


@pytest.mark.skipif(not _TED_IMPORTABLE, reason="flair.ted not importable")
def test_minus_strand_both_ends_clamped_not_dropped():
    """- strand: biological TSS = genomic END, biological TTS = genomic START.
    A TSS peak BEFORE the last acceptor (genomic end too small) and a TTS peak
    PAST the first donor (genomic start too large) must both be clamped, on the
    correct genomic sides, with valid emitted geometry and no drop."""
    # Chain genomic: exon1 [800,1200) intron exon2 [1500,5200).
    # first donor=1200 (genomic start ceiling), last acceptor=1500 (end floor).
    juncs = [_JUNC(1200, 1500)]
    # - strand: TSS peak (genomic end) at 1450 < last acceptor 1500 -> invalid.
    #           TTS peak (genomic start) at 1400 > first donor 1200 -> invalid.
    approved = {
        ("chr1", "-"): {
            "tss": _np.array([1450], dtype=_np.int64),
            "tts": _np.array([1400], dtype=_np.int64),
            "tss_scores": {1450: {"depth_score": 0.5, "model_score": 0.6,
                                  "annotation_score": 0.4, "annot_dist_bp": 10}},
            "tts_scores": {1400: {"depth_score": 0.5, "model_score": 0.6,
                                  "annotation_score": 0.4, "annot_dist_bp": 10}},
            "tss_reality": {1450: 0.5}, "tts_reality": {1400: 0.5},
            "tss_approved": _np.array([1450], dtype=_np.int64),
            "tts_approved": _np.array([1400], dtype=_np.int64),
        }
    }
    # - strand read: genomic start=1410 (TTS side), genomic end=1460 (TSS side).
    reads = [_RR("chr1", "-", 1410, 1460, "minus_overshoot", (False, False))]
    out = _ted.ted_global_assign_reads(
        reads, "chr1", approved, tss_snap=50, tts_snap=25, scorer=None,
        val_kwargs_tss=_VK, val_kwargs_tts=_VK,
        threshold_tss=0.5, threshold_tts=0.5, chain_juncs=juncs)
    assert len(out) == 1, "minus-strand read dropped instead of clamped"
    _c, rep_start, rep_end, _s, _rn, _ns, _sc, _t = out[0]
    # Genomic invariants hold regardless of strand.
    assert rep_start < juncs[0].start, (
        f"genomic rep_start {rep_start} not < first donor {juncs[0].start}")
    assert rep_end > juncs[-1].end, (
        f"genomic rep_end {rep_end} not > last acceptor {juncs[-1].end}")
    assert _emitted_exons_valid(rep_start, rep_end, juncs), (
        "minus-strand emitted geometry invalid")


@pytest.mark.skipif(not _TED_IMPORTABLE, reason="flair.ted not importable")
def test_valid_ends_unchanged_and_single_exon_untouched():
    """Clamp must be a no-op when ends are already valid, and single-exon
    chains (no juncs) must be left entirely alone (chain_juncs empty/None)."""
    juncs = [_JUNC(1200, 1500)]
    approved = {
        ("chr1", "+"): {
            "tss": _np.array([1000], dtype=_np.int64),
            "tts": _np.array([5000], dtype=_np.int64),
            "tss_scores": {1000: {"depth_score": 0.5, "model_score": 0.6,
                                  "annotation_score": 0.4, "annot_dist_bp": 10}},
            "tts_scores": {5000: {"depth_score": 0.5, "model_score": 0.6,
                                  "annotation_score": 0.4, "annot_dist_bp": 10}},
            "tss_reality": {1000: 0.5}, "tts_reality": {5000: 0.5},
            "tss_approved": _np.array([1000], dtype=_np.int64),
            "tts_approved": _np.array([5000], dtype=_np.int64),
        }
    }
    reads = [_RR("chr1", "+", 1002, 5002, "valid", (False, False))]
    out = _ted.ted_global_assign_reads(
        reads, "chr1", approved, tss_snap=50, tts_snap=25, scorer=None,
        val_kwargs_tss=_VK, val_kwargs_tts=_VK,
        threshold_tss=0.5, threshold_tts=0.5, chain_juncs=juncs)
    assert len(out) == 1
    _c, rep_start, rep_end, *_ = out[0]
    # Unchanged snap to valid peaks.
    assert (rep_start, rep_end) == (1000, 5000), "valid ends were perturbed"

    # Single-exon: no juncs -> no clamp, snap as usual.
    out_se = _ted.ted_global_assign_reads(
        reads, "chr1", approved, tss_snap=50, tts_snap=25, scorer=None,
        val_kwargs_tss=_VK, val_kwargs_tts=_VK,
        threshold_tss=0.5, threshold_tts=0.5, chain_juncs=None)
    assert len(out_se) == 1
    _c2, rs2, re2, *_ = out_se[0]
    assert _emitted_exons_valid(rs2, re2, [])


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v", "--tb=short"]))
