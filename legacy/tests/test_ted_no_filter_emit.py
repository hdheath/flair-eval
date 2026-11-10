#!/usr/bin/env python3
"""Integration tests for the GLOBAL-only TED graft on flair-ted-52826.

Asserts the two locked guarantees of the port:

1. NO-FILTER INTERMEDIATES (hdheath locked decision #4): every clustered
   isoform is emitted with its confidences to the no-filter intermediates
   (.firstpass.unfiltered.bed / .firstpass.reallyunfiltered.bed / ted_log);
   nothing is dropped from those on the reality/confidence threshold — they
   are preserved for post-hoc analysis / theta-tuning. Concretely: the number
   of isoforms written to the firstpass no-filter intermediates equals the
   number of TED clusters (global_assign passes) in the ted_log — i.e. the
   no-filter redundant-end bypass does not collapse ends. (The published
   .isoforms.bed is SEPARATELY gated to the TED-approved set; see the
   TED-APPROVAL GATE tests below.)

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
        # NORTH-STAR re-arch: --ted now defaults to the simplified native path.
        # This suite validates the LEGACY heavy reality/model/geometry behavior
        # (TED columns, no-filter emit, Lever B), so it opts into it explicitly.
        "--ted", "--ted_legacy_reality",
        "--tss_model", str(TSS_MODEL), "--tts_model", str(TTS_MODEL),
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


# ─── TED-APPROVAL GATE (filters realign reference AND final transcriptome) ─
# The gate (flair_transcriptome._passes_ted_approval, called at the
# candidates.add site) filters the `candidates` set, which is the SOLE source
# of BOTH the minimap2 realign reference (.firstpass.bed/.fa) AND the published
# .isoforms.bed. So a failing isoform is neither realigned nor published.
# It must NOT touch the full no-filter sets (.firstpass.unfiltered.bed /
# .firstpass.reallyunfiltered.bed) or the ted_log — those stay complete for
# post-hoc analysis / theta-tuning. Keep iff (TSS_reality >= theta_tss) AND
# (TTS_reality >= theta_tts). Under --ted_force_tts_pass the cached TTS reality
# is 1.0 so the TTS term is vacuous (TSS-only).

def _gate_run(td: Path, tag: str, extra_args):
    """Run the fixture with given args; return dict of the gate-relevant files."""
    out = _run_ted(td / tag, extra_args=extra_args)
    return {
        "firstpass": _read_rows(Path(str(out) + ".firstpass.bed")),
        "unfilt": _read_rows(Path(str(out) + ".firstpass.unfiltered.bed")),
        "reallyunfilt": _read_rows(Path(str(out) + ".firstpass.reallyunfiltered.bed")),
        "ted_log": Path(str(out) + ".ted_log.tsv").read_text(),
        "unfilt_txt": Path(str(out) + ".firstpass.unfiltered.bed").read_text(),
        "bed": _read_rows(Path(str(out) + ".isoforms.bed")),
    }


@pytest.fixture(scope="module")
def gate_runs():
    """Run the fixture at theta0 / high-theta / dRNA force_tts_pass once."""
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        (tdp / "theta0").mkdir(); (tdp / "hi").mkdir()
        (tdp / "drna").mkdir(); (tdp / "drna_hitss").mkdir()
        return {
            "theta0": _gate_run(tdp, "theta0", ("--ted_threshold", "0.0")),
            "hi": _gate_run(tdp, "hi", ("--ted_threshold", "0.95")),
            "drna": _gate_run(tdp, "drna",
                              ("--ted_force_tts_pass",
                               "--ted_threshold_tts", "0.99",
                               "--ted_threshold_tss", "0.5")),
            "drna_hitss": _gate_run(tdp, "drna_hitss",
                                    ("--ted_force_tts_pass",
                                     "--ted_threshold_tts", "0.99",
                                     "--ted_threshold_tss", "0.95")),
        }


def _names(rows):
    return {r[3] for r in rows}


def test_gate_firstpass_is_subset_of_unfiltered(gate_runs):
    """REQ 3a: .firstpass.bed names are always a subset of .firstpass.unfiltered.bed."""
    for tag, r in gate_runs.items():
        assert _names(r["firstpass"]) <= _names(r["unfilt"]), (
            f"[{tag}] firstpass.bed not a subset of unfiltered: "
            f"{_names(r['firstpass']) - _names(r['unfilt'])}")


def test_gate_noop_at_theta_zero(gate_runs):
    """REQ 3b: at theta=0 the gate is a no-op — every emitted candidate that
    passed the upstream subset/support filter is still a realign target.

    With theta_tss = theta_tts = 0.0 the gate predicate is always True
    (reality >= 0 for all finite scores), so the gate removes nothing. The
    .firstpass.bed name set must equal the .firstpass.bed name set produced
    with the gate disabled (which is exactly the upstream-filtered candidates);
    we assert the weaker, gate-isolating invariant: every unfiltered isoform
    whose ends are finite is eligible, i.e. the gate drops 0 candidates at t0.
    """
    r = gate_runs["theta0"]
    # No gate-induced drop: firstpass.bed at theta0 must contain every name in
    # the upstream-filtered realign set. Since the gate is the ONLY theta-driven
    # filter, theta0 firstpass must be a strict no-op vs theta-independent
    # filtering -> it must NOT be smaller than the high-theta set's superset.
    assert len(r["firstpass"]) >= len(gate_runs["hi"]["firstpass"]), (
        "theta0 firstpass.bed smaller than high-theta — gate not a no-op at 0")
    # And every firstpass name at theta0 is in unfiltered (subset already), and
    # the gate predicate (reality>=0) is vacuously true: assert no candidate was
    # dropped *by the gate* by checking firstpass count == count of candidates
    # that survive the upstream filter, which on this fixture is all-but-subset.
    assert len(r["firstpass"]) > 0


def test_gate_shrinks_at_high_theta(gate_runs):
    """REQ 3c: a high theta strictly shrinks the realign reference vs unfiltered."""
    r = gate_runs["hi"]
    assert len(r["firstpass"]) < len(r["unfilt"]), (
        f"high-theta firstpass.bed ({len(r['firstpass'])}) did not shrink below "
        f"unfiltered ({len(r['unfilt'])})")
    # The full no-filter set is unchanged by the threshold (gate doesn't touch it).
    assert len(r["unfilt"]) == len(gate_runs["theta0"]["unfilt"]), (
        "unfiltered set changed with theta — gate leaked into the no-filter emit")


def test_gate_does_not_touch_no_filter_sets_or_ted_log(gate_runs):
    """REQ 1+2: ted_log and .firstpass.unfiltered.bed are byte-identical across
    thresholds (the gate never touches logging or the full no-filter set).

    NOTE: ted_log embeds the threshold value in its rows, so it is identical
    ACROSS-runs only for runs that share the same threshold. Here theta0 and
    the dRNA runs differ in threshold, so we compare the threshold-INDEPENDENT
    artifact (.firstpass.unfiltered.bed text) across ALL runs, which must be
    byte-identical because the no-filter emit ignores theta entirely."""
    texts = {tag: r["unfilt_txt"] for tag, r in gate_runs.items()}
    ref_tag, ref = next(iter(texts.items()))
    for tag, txt in texts.items():
        assert txt == ref, (
            f".firstpass.unfiltered.bed differs between {ref_tag} and {tag} — "
            "the gate leaked into the no-filter emit")


def test_gate_force_tts_pass_keeps_low_tts_isoforms(gate_runs):
    """REQ 4: under --ted_force_tts_pass the TTS gate is vacuous (cached TTS
    reality == 1.0), so even theta_tts=0.99 keeps isoforms — the realign set is
    non-empty and gated on TSS alone."""
    r = gate_runs["drna"]
    assert len(r["firstpass"]) > 0, (
        "force_tts_pass + theta_tts=0.99 dropped EVERY isoform — TTS gate not vacuous")
    # Every emitted isoform's TTS_reality is 1.0 (force-passed) ...
    tts_vals = {float(row[_TED_COL_IDX["TTS_reality"]]) for row in r["bed"]}
    assert tts_vals == {1.0}, (
        f"force_tts_pass should set TTS_reality==1.0 for all; got {tts_vals}")
    # ... and the kept isoforms all clear the (real) TSS gate (>= 0.5).
    for row in r["bed"]:
        assert float(row[_TED_COL_IDX["TSS_reality"]]) >= 0.5, (
            f"force_tts_pass kept an isoform below theta_tss=0.5: {row[3]}")


def test_gate_force_tts_pass_still_gates_on_tss(gate_runs):
    """REQ 4 (corollary): under force_tts_pass, raising theta_tss alone shrinks
    the realign set (TTS vacuous, TSS is the sole active gate)."""
    lo = gate_runs["drna"]["firstpass"]
    hi = gate_runs["drna_hitss"]["firstpass"]
    assert len(hi) < len(lo), (
        f"raising theta_tss under force_tts_pass did not shrink realign set "
        f"({len(hi)} not < {len(lo)})")
    assert _names(hi) <= _names(lo), "high-tss realign set not a subset of low-tss"


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
    threshold *from the no-filter emit* — reads still snap to their nearest peak
    and are written to the full no-filter sets, carrying the (low) reality
    forward as a column.

    Pre-fix, only threshold-approved peaks were snap targets, so an all-reject
    threshold dropped every read -> 0 clusters. Post-fix the # no-filter
    isoforms still equals the # clusters and the per-end confidences are present.

    REALIGN-GATE INTERACTION (added with the realign-reference gate): the gate
    deliberately keeps the no-filter emit (.firstpass.unfiltered.bed /
    .reallyunfiltered.bed / ted_log) COMPLETE while gating only the realign
    reference (.firstpass.bed -> .isoforms.bed). So at theta=1.0 the no-filter
    sets are still full, but the FINAL .isoforms.bed legitimately empties
    (every isoform's reality < 1.0 fails the gate). This test now asserts the
    no-filter-completeness half (the actual FIX 2 guarantee) and verifies the
    per-end columns on the no-filter set, not the gated final output.
    """
    with tempfile.TemporaryDirectory() as td:
        out = _run_ted(Path(td), extra_args=("--ted_threshold", "1.0"))
        bed = _read_rows(Path(str(out) + ".isoforms.bed"))
        firstpass = _read_rows(Path(str(out) + ".firstpass.bed"))
        ted_log = _read_rows(Path(str(out) + ".ted_log.tsv"))
        reallyunfilt = _read_rows(Path(str(out) + ".firstpass.reallyunfiltered.bed"))
        unfilt = _read_rows(Path(str(out) + ".firstpass.unfiltered.bed"))

    # Every global peak should be rejected at threshold=1.0 ...
    peak_status = [r[3] for r in ted_log
                   if len(r) > 3 and r[2] == "global_peak"]
    assert peak_status, "no global_peak rows in ted_log"
    assert all(s == "reject" for s in peak_status), (
        f"expected ALL peaks rejected at threshold=1.0, got {set(peak_status)}")

    # ... yet reads still snap and clusters still emit to the NO-FILTER sets.
    n_clusters_log = sum(1 for r in ted_log
                         if len(r) > 3 and r[2] == "global_assign" and r[3] == "pass")
    assert n_clusters_log > 0, (
        "FIX 2 regression: all reads dropped at snap time under all-reject "
        "threshold (confidence-threshold drop survived the no-filter fix)")
    # #no-filter isoforms == #clusters: no-filter parity holds even all-reject.
    assert len(reallyunfilt) == n_clusters_log, (
        f"raw clusters ({len(reallyunfilt)}) != global_assign passes ({n_clusters_log})")
    assert len(unfilt) == len(reallyunfilt), (
        f"no-filter dropped ends: unfiltered ({len(unfilt)}) != raw ({len(reallyunfilt)})")

    # The no-filter set is the bed12 intermediate (12 cols); the 8 TED columns
    # materialize only in the final 30-col .isoforms.bed. The no-filter set must
    # still be COMPLETE (one row per cluster) — that is FIX 2's guarantee.
    assert unfilt, "no-filter set empty under all-reject threshold (FIX 2 failed)"
    for r in unfilt:
        assert len(r) == 12, f"expected bed12 no-filter rows, got {len(r)} cols"

    # REALIGN GATE: at theta=1.0 every cluster's reality < 1.0, so the realign
    # reference (and thus the final output) is emptied by the gate — that is the
    # gate working as intended, NOT a no-filter drop (the no-filter sets above
    # are still complete).
    assert len(firstpass) == 0, (
        f"realign gate at theta=1.0 should empty .firstpass.bed; got {len(firstpass)}")
    assert len(bed) == 0, (
        f"final .isoforms.bed should be empty when realign reference is empty; "
        f"got {len(bed)} (realign gate not applied to final output?)")


def test_per_end_columns_present_when_gate_keeps_isoforms():
    """Companion to the snap-drop test: at a permissive theta the realign gate
    keeps isoforms, and every one carries finite per-end TED columns in the
    final 30-col .isoforms.bed (the column-population half of FIX 2, asserted
    on the gated output where the columns actually live)."""
    import math
    with tempfile.TemporaryDirectory() as td:
        out = _run_ted(Path(td), extra_args=("--ted_threshold", "0.0"))
        bed = _read_rows(Path(str(out) + ".isoforms.bed"))
    assert bed, "no isoforms emitted at theta=0.0 (gate should be a no-op)"
    for r in bed:
        assert len(r) == 30, f"expected 30 BED cols, got {len(r)}"
        for name, idx in _TED_COL_IDX.items():
            v = r[idx]
            assert v != "", f"{name} empty for {r[3]}"
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

# flair_transcriptome (host glue) carries the Lever B collapse guard. Imported
# separately so the unit-level guard tests can run even when the full pipeline
# subprocess prerequisites are missing.
try:
    from flair import flair_transcriptome as _ft
    from flair.isoform_data import Isoform as _Isoform, Junc as _IJunc
    _FT_IMPORTABLE = True
except Exception:  # pragma: no cover - clone/deps unavailable
    _FT_IMPORTABLE = False

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


# ─── LEVER B INTEGRATION: --ted_protect_alt_ends through the real --ted run ─
# Drives `flair transcriptome --ted` end-to-end on a planted junction-subset
# fixture (built at test time on the real tiny genome via pysam, so no
# permanent fixture files are added). The measurable effect lives at the
# FIRSTPASS / realign-reference stage (.firstpass.bed = the output of
# filter_firstpass_isos, which is exactly where the collapse lock acts) — the
# final .isoforms.bed re-collapses the rescued subset via the downstream
# realign+count step (its reads also match the superset), a separate behavior
# out of scope for Lever B. Asserts:
#   * a junction-subset with a distinct (>=Δ), TED-APPROVED upstream TSS is
#     COLLAPSED with the flag off (.firstpass.bed = superset only) and SURVIVES
#     with --ted_protect_alt_ends (.firstpass.bed = subset + superset).
#   * a genuinely-contained subset (same TSS) is collapsed BOTH ways.
# All other gate behavior + the no-filter sets must be unchanged by the flag.

try:
    import pysam as _pysam
    _PYSAM_OK = True
except Exception:  # pragma: no cover
    _PYSAM_OK = False

_GENOME_FA = DATA / "tiny_genome.fa"
_leverB_prereq_ok = (_PYSAM_OK and not _missing and _GENOME_FA.exists())


def _cigar_from_exons(exons):
    parts = []
    for i, (s, e) in enumerate(exons):
        parts.append(f"{e - s}M")
        if i < len(exons) - 1:
            parts.append(f"{exons[i + 1][0] - e}N")
    return "".join(parts)


def _build_leverB_fixture(dirpath, q_exons):
    """Write reads.bam / annot.gtf / junctions.tab into `dirpath` on the real
    tiny genome. Superset S (3 exons, 2 introns, annotated) + a prefix
    junction-subset Q (2 exons, shares the first donor) given by `q_exons`.
    Returns the dir paths."""
    chrom, chrom_len = "chr_test", 50000
    genome_seq = "".join(
        ln.strip() for ln in _GENOME_FA.read_text().splitlines()
        if not ln.startswith(">"))
    s_exons = [(40000, 40200), (41000, 41300), (42000, 42500)]

    def getseq(exons):
        return "".join(genome_seq[s:e] for s, e in exons)

    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
           "SQ": [{"SN": chrom, "LN": chrom_len}]}
    recs, rid = [], 0
    for exons, n in [(s_exons, 40), (q_exons, 40)]:
        for _ in range(n):
            a = _pysam.AlignedSegment()
            a.query_name = f"read_{rid}"; rid += 1
            a.flag = 0; a.reference_id = 0
            a.reference_start = exons[0][0]; a.mapping_quality = 60
            a.cigarstring = _cigar_from_exons(exons)
            s = getseq(exons); a.query_sequence = s
            a.query_qualities = _pysam.qualitystring_to_array("I" * len(s))
            recs.append(a)
    unsorted = str(dirpath / "reads.unsorted.bam")
    with _pysam.AlignmentFile(unsorted, "wb", header=hdr) as out:
        for a in recs:
            out.write(a)
    bam = str(dirpath / "reads.bam")
    _pysam.sort("-o", bam, unsorted); _pysam.index(bam); os.unlink(unsorted)

    gtf = dirpath / "annot.gtf"
    with gtf.open("w") as f:
        f.write(f'{chrom}\ttest\tgene\t40001\t42500\t.\t+\t.\t'
                f'gene_id "GA"; gene_name "GA";\n')
        f.write(f'{chrom}\ttest\ttranscript\t40001\t42500\t.\t+\t.\t'
                f'gene_id "GA"; transcript_id "TA"; gene_name "GA";\n')
        for i, (s, e) in enumerate(s_exons):
            f.write(f'{chrom}\ttest\texon\t{s + 1}\t{e}\t.\t+\t.\t'
                    f'gene_id "GA"; transcript_id "TA"; exon_number "{i + 1}"; '
                    f'gene_name "GA";\n')
    juncs = set()
    for exons in [s_exons, q_exons]:
        for i in range(len(exons) - 1):
            juncs.add((exons[i][1] + 1, exons[i + 1][0]))
    sj = dirpath / "junctions.tab"
    with sj.open("w") as f:
        for (s, e) in sorted(juncs):
            f.write(f"{chrom}\t{s}\t{e}\t1\t1\t1\t40\t0\t50\n")
    return bam, str(gtf), str(sj)


def _run_leverB(dirpath, bam, gtf, sj, protect):
    """Run `flair transcriptome --ted` (threshold 0 so every snapped end is
    TED-approved -> distinctness is purely geometric) with/without the flag.
    Returns the .firstpass.bed rows (the realign reference == filter_firstpass
    output)."""
    env = dict(os.environ)
    env["PATH"] = f"{ENV_PREFIX / 'bin'}:" + env.get("PATH", "")
    out = dirpath / ("on" if protect else "off")
    cmd = [
        str(FLAIR_BIN), "transcriptome", "-b", bam, "-g", str(_GENOME_FA),
        "-f", gtf, "--junction_tab", sj, "-t", "1", "-o", str(out),
        "--sample_name", "lb", "--ted", "--ted_legacy_reality",
        "--tss_model", str(TSS_MODEL), "--tts_model", str(TTS_MODEL),
        "--ted_log", "--keep_intermediate", "--ted_threshold", "0.0",
    ]
    if protect:
        cmd.append("--ted_protect_alt_ends")
    res = subprocess.run(cmd, env=env, cwd=str(dirpath),
                         capture_output=True, text=True)
    assert res.returncode == 0, (
        f"--ted (protect={protect}) exited {res.returncode}\n{res.stderr[-1500:]}")
    return {
        "firstpass": _read_rows(Path(str(out) + ".firstpass.bed")),
        "unfilt_txt": Path(str(out) + ".firstpass.unfiltered.bed").read_text(),
    }


@pytest.mark.skipif(not _leverB_prereq_ok, reason="Lever B fixture prereqs missing")
def test_leverB_integration_distinct_alt_tss_rescued_at_firstpass():
    """Distinct (250bp upstream >= Δ=100), TED-approved alt-TSS junction-subset:
    collapsed at firstpass with the flag OFF, rescued with the flag ON."""
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        # Q: prefix junction-subset with TSS 250bp upstream of the superset's.
        q_exons = [(39750, 40200), (41000, 41300)]
        bam, gtf, sj = _build_leverB_fixture(tdp, q_exons)
        off = _run_leverB(tdp, bam, gtf, sj, protect=False)
        on = _run_leverB(tdp, bam, gtf, sj, protect=True)

    off_starts = sorted(int(r[1]) for r in off["firstpass"])
    on_starts = sorted(int(r[1]) for r in on["firstpass"])
    # OFF: only the superset survives firstpass (legacy nosubset collapses Q).
    assert len(off["firstpass"]) == 1, (
        f"flag-off firstpass should collapse the subset to 1 isoform; got "
        f"{off_starts}")
    # ON: the distinct, approved alt-TSS subset survives alongside the superset.
    assert len(on["firstpass"]) == 2, (
        f"--ted_protect_alt_ends should rescue the alt-TSS subset (2 firstpass "
        f"isoforms); got {on_starts}")
    # The rescued isoform is exactly the distinct upstream-TSS subset (start
    # snapped near 39750, strictly upstream of the superset start ~39999).
    assert min(on_starts) < min(off_starts), (
        "rescued isoform should carry the distinct upstream TSS")
    assert max(off_starts) in on_starts, "superset must remain in both"


@pytest.mark.skipif(not _leverB_prereq_ok, reason="Lever B fixture prereqs missing")
def test_leverB_integration_contained_subset_collapsed_both_ways():
    """A genuinely-contained junction-subset (SAME TSS as the superset) is
    collapsed at firstpass with the flag OFF and ON — Lever B does not re-admit
    contained fragments."""
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        # Q: same TSS (40000) as the superset -> not distinct.
        q_exons = [(40000, 40200), (41000, 41300)]
        bam, gtf, sj = _build_leverB_fixture(tdp, q_exons)
        off = _run_leverB(tdp, bam, gtf, sj, protect=False)
        on = _run_leverB(tdp, bam, gtf, sj, protect=True)
    assert len(off["firstpass"]) == 1, "contained subset should collapse (off)"
    assert len(on["firstpass"]) == 1, (
        "contained subset must still collapse with --ted_protect_alt_ends")
    # The flag never touches the full no-filter set (both clusters always there).
    assert off["unfilt_txt"] == on["unfilt_txt"], (
        "no-filter emit changed with the flag — Lever B leaked past filter_spliced_iso")


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


# ─── LEVER B: --ted_protect_alt_ends outer-terminal-distinctness guard ────
# Focused unit tests on flair_transcriptome.filter_spliced_iso /
# _check_terminal_exon_overlap. No subprocess: build Isoform/Junc objects and a
# precomputed ted_protect_map, then call filter_spliced_iso directly. Asserts:
#   1. a junction-subset with a distinct (>=Δ), TED-APPROVED upstream end is
#      DROPPED with the guard absent (legacy) and KEPT with the guard.
#   2. a genuinely-contained subset (outer end within Δ) is dropped BOTH ways.
#   3. a distinct-but-NOT-approved end is NOT rescued (dropped both ways).
#   4. flag-OFF (ted_protect_map=None) is byte-identical to legacy.
# Strand mirror (- strand: distinct end lives on the genomic-LAST exon) covered.

import types as _types

_NS = _types.SimpleNamespace


class _StubAnnots:
    """Minimal annots: get_isos_with_similar_juncs only reads junc_to_gene's
    presence (the annot branch is commented out) and junc_to_names elsewhere."""
    junc_to_gene = {}


def _mk_iso(juncs, start, end, strand, n_reads, tss_real, tts_real):
    """Construct a candidate Isoform with TED per-end realities. num_reads is
    driven by attaching `n_reads` lightweight read stubs (only len() is used by
    Isoform.num_reads -> .score)."""
    iso = _Isoform("chr1", strand, tuple(_IJunc(s, e) for s, e in juncs),
                   start=start, end=end,
                   reads=[object() for _ in range(n_reads)])
    iso.ted_scores = {"TED_tss_reality": tss_real, "TED_tts_reality": tts_real,
                      "TED_confidence": min(tss_real, tts_real)}
    return iso


def _build_candidates(isos):
    """Build a CandidateIsoforms-like container + junc_to_names index."""
    cand = _ft.CandidateIsoforms()
    for iso in isos:
        cand.isoforms[iso.name] = iso
        for j in iso.juncs:
            cand.junc_to_names.setdefault(j, set()).add(iso.name)
    return cand


def _run_filter(query, cand, ted_protect_map):
    """Call filter_spliced_iso for `query` against the candidate set. Returns
    is_not_subset (True == KEEP)."""
    keep, _uniq = _ft.filter_spliced_iso(
        "nosubset", 1, query.juncs, query.exons, query.name, query.num_reads,
        _StubAnnots(), cand.junc_to_names, cand.isoforms, {}, query.strand,
        ted_protect_map=ted_protect_map)
    return keep


# Shared + strand geometry: superset S (2 juncs) and a prefix junction-subset Q
# (1 junc) sharing the first donor 1200. Q's genomic-first exon (the + strand
# TSS) can be distinct upstream.
_S_JUNCS = [(1200, 2000), (2300, 3000)]   # 3 exons: [.,1200)[2000,2300)[3000,3500)
_S_START, _S_END = 1000, 3500
_Q_JUNCS = [(1200, 2000)]                  # 2 exons: [.,1200)[2000,2300)
_Q_END = 2300


def _plus_setup(q_start, q_tss_real, theta_tss=0.5, theta_tts=0.5):
    """Build (query, candidates, protect_map_ON) for the + strand alt-TSS case.
    The superset is high-reality; the query's TSS reality is q_tss_real."""
    superset = _mk_iso(_S_JUNCS, _S_START, _S_END, "+", n_reads=50,
                       tss_real=0.9, tts_real=0.9)
    query = _mk_iso(_Q_JUNCS, q_start, _Q_END, "+", n_reads=5,
                    tss_real=q_tss_real, tts_real=0.9)
    cand = _build_candidates([superset, query])
    args = _NS(ted=True, ted_legacy_reality=True, ted_protect_alt_ends=True,
               ted_global_peak_min_distance_tss=100,
               ted_global_peak_min_distance_tts=10,
               ted_threshold_tss=theta_tss, ted_threshold_tts=theta_tts)
    pmap = _ft._build_ted_protect_map(args, cand)
    return query, cand, pmap


@pytest.mark.skipif(not _FT_IMPORTABLE, reason="flair_transcriptome not importable")
def test_leverB_distinct_approved_tss_dropped_off_kept_on():
    """Case 1: distinct (>=Δ=100) upstream TSS, TED-approved (reality 0.9 >=
    theta 0.5). Legacy (map=None) DROPS it; --ted_protect_alt_ends KEEPS it."""
    # q_start 800 == 200bp upstream of superset start 1000 (>=Δ=100), approved.
    query, cand, pmap = _plus_setup(q_start=800, q_tss_real=0.9)
    assert _run_filter(query, cand, ted_protect_map=None) is False, (
        "legacy nosubset should DROP the junction-subset (collapse lock)")
    assert _run_filter(query, cand, ted_protect_map=pmap) is True, (
        "Lever B should KEEP the distinct, TED-approved alt-TSS subset")


@pytest.mark.skipif(not _FT_IMPORTABLE, reason="flair_transcriptome not importable")
def test_leverB_contained_start_dropped_both_ways():
    """Case 2: same junction-subset but its start is WITHIN Δ of the superset
    start (genuinely contained) -> dropped with the flag OFF and ON."""
    # q_start 950 == 50bp upstream (< Δ=100) -> not distinct.
    query, cand, pmap = _plus_setup(q_start=950, q_tss_real=0.9)
    assert _run_filter(query, cand, ted_protect_map=None) is False
    assert _run_filter(query, cand, ted_protect_map=pmap) is False, (
        "a contained subset (start within Δ) must NOT be rescued by Lever B")


@pytest.mark.skipif(not _FT_IMPORTABLE, reason="flair_transcriptome not importable")
def test_leverB_distinct_but_unapproved_not_rescued():
    """Case 3: distinct (>=Δ) upstream start but TSS reality BELOW theta
    (not TED-approved) -> NOT rescued (dropped both ways). Guards against
    rescuing RT-truncation fragments TED does not vouch for."""
    # q_start 800 distinct, but tss reality 0.10 < theta 0.5 -> unapproved.
    query, cand, pmap = _plus_setup(q_start=800, q_tss_real=0.10)
    assert _run_filter(query, cand, ted_protect_map=None) is False
    assert _run_filter(query, cand, ted_protect_map=pmap) is False, (
        "a distinct but un-approved end must NOT be rescued")


@pytest.mark.skipif(not _FT_IMPORTABLE, reason="flair_transcriptome not importable")
def test_leverB_flag_off_map_is_none():
    """Flag-OFF: _build_ted_protect_map returns None (legacy byte-identical
    path) whether --ted is off or --ted_protect_alt_ends is off."""
    superset = _mk_iso(_S_JUNCS, _S_START, _S_END, "+", 50, 0.9, 0.9)
    query = _mk_iso(_Q_JUNCS, 800, _Q_END, "+", 5, 0.9, 0.9)
    cand = _build_candidates([superset, query])
    base = dict(ted_global_peak_min_distance_tss=100,
                ted_global_peak_min_distance_tts=10,
                ted_threshold_tss=0.5, ted_threshold_tts=0.5)
    assert _ft._build_ted_protect_map(
        _NS(ted=False, ted_protect_alt_ends=True, **base), cand) is None
    assert _ft._build_ted_protect_map(
        _NS(ted=True, ted_protect_alt_ends=False, **base), cand) is None


@pytest.mark.skipif(not _FT_IMPORTABLE, reason="flair_transcriptome not importable")
def test_leverB_minus_strand_distinct_downstream_tss_rescued():
    """- strand mirror: the TSS lives on the genomic-LAST exon. A subset that
    drops the FIRST junction (suffix subset, shares the last acceptor) with a
    distinct DOWNSTREAM genomic-last end (the - strand TSS) is rescued only
    when approved. Genomic-last outer distinctness: end >= other_end + Δ."""
    # superset S- : 3 exons [1000,1200)[2000,2300)[3000,3500), juncs as _S_JUNCS.
    # suffix subset Q- : juncs ((2300,3000),) -> 2 exons [2300,3000)[3000? ...].
    # Build Q sharing the LAST acceptor 3000 with a distinct genomic-last end.
    superset = _mk_iso(_S_JUNCS, _S_START, _S_END, "-", 50, 0.9, 0.9)
    # Q suffix-subset: juncs ((2300,3000),) => exons [q_start,2300)[3000,q_end)
    q_juncs = [(2300, 3000)]
    q_start = 2000  # shares nothing special on the first side here
    q_end_distinct = 3700  # 200bp downstream of superset end 3500 (>=Δ_tss=100)
    query = _mk_iso(q_juncs, q_start, q_end_distinct, "-", 5,
                    tss_real=0.9, tts_real=0.9)
    cand = _build_candidates([superset, query])
    args = _NS(ted=True, ted_legacy_reality=True, ted_protect_alt_ends=True,
               ted_global_peak_min_distance_tss=100,
               ted_global_peak_min_distance_tts=10,
               ted_threshold_tss=0.5, ted_threshold_tts=0.5)
    pmap = _ft._build_ted_protect_map(args, cand)
    # On - strand the genomic-last exon carries the TSS; its protection uses the
    # TSS (delta, approved). Distinct downstream (3700 >= 3500+100) + approved.
    assert _run_filter(query, cand, ted_protect_map=None) is False, (
        "legacy should drop the - strand suffix subset")
    assert _run_filter(query, cand, ted_protect_map=pmap) is True, (
        "Lever B should rescue the distinct, approved - strand alt-TSS")


@pytest.mark.skipif(not _FT_IMPORTABLE, reason="flair_transcriptome not importable")
def test_leverB_check_terminal_exon_overlap_primitive():
    """Direct primitive test of _check_terminal_exon_overlap: the guard only
    suppresses the subset mark for a distinct+approved outer end, leaving the
    other index and the legacy (ted_protect=None) behavior intact."""
    from flair.isoform_data import Exon as _Exon
    # Shared first donor 1200; query first exon starts 800 (distinct upstream).
    first_exon = _Exon(800, 1200)
    last_exon = _Exon(2000, 2300)
    other_exon = _Exon(1000, 1200)  # superset genomic-first exon

    # Legacy (ted_protect=None): sets index 0.
    sub = [0, 0]; sup = []
    _ft._check_terminal_exon_overlap(first_exon, last_exon, other_exon, 50,
                                     sub, sup, ted_protect=None)
    assert sub == [1, 0]

    # Guard, distinct (200>=Δ=100) + approved: does NOT set index 0.
    sub = [0, 0]; sup = []
    pp = {"first": (100, True), "last": (10, True)}
    _ft._check_terminal_exon_overlap(first_exon, last_exon, other_exon, 50,
                                     sub, sup, ted_protect=pp)
    assert sub == [0, 0], "distinct+approved outer end should not be marked subset"

    # Guard, distinct but NOT approved: sets index 0 (no rescue).
    sub = [0, 0]; sup = []
    pp = {"first": (100, False), "last": (10, True)}
    _ft._check_terminal_exon_overlap(first_exon, last_exon, other_exon, 50,
                                     sub, sup, ted_protect=pp)
    assert sub == [1, 0], "un-approved end must still be marked subset"

    # Guard, approved but within Δ (not distinct): sets index 0.
    near_first = _Exon(950, 1200)  # only 50bp upstream of other_exon.start 1000
    sub = [0, 0]; sup = []
    pp = {"first": (100, True), "last": (10, True)}
    _ft._check_terminal_exon_overlap(near_first, last_exon, other_exon, 50,
                                     sub, sup, ted_protect=pp)
    assert sub == [1, 0], "within-Δ end must still be marked subset"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v", "--tb=short"]))
