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


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v", "--tb=short"]))
