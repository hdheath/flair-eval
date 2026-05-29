"""Regression tests for the EMPTY/0-isoform eval-chain landmines.

Single-chromosome / single-sample runs (e.g. WTC11 chr22) routinely produce
tools that emit ZERO isoforms in the region (IsoSeq found 0 on chr22). An
empty-but-valid result must NOT abort the pipeline: every FATAL eval process
(no errorStrategy 'ignore') must exit 0 and emit a valid (possibly
header-only) output file.

These tests lock in that contract for the scripts that run inside the FATAL
Nextflow processes the chr22 run executes:
  * gtf_to_bed12.py           (run under ted.py subprocess.run(check=True))
  * ted.py                    (Evaluation; simplified GTF path + BED path)
  * flair_eval.py             (Evaluation)
  * synthesize_evaluations.py (Evaluation)
  * ted_end_precision.py      (TedEndPrecision; peaks path + gtf-only path)
  * firstpass_vs_final.py     (FirstpassComparison)

They invoke each script exactly as the Nextflow module does (subprocess,
checking the exit code) on a 0-row BED / GTF, then assert the named output
the module will `mv`/consume actually exists with a header.

Run env: must have pandas + matplotlib + bedtools + samtools on PATH
(nextflow_env mirrors the pipeline's analysis env).
"""
from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
BIN = REPO_ROOT / "bin"
EVAL = BIN / "evaluation"
DATA = REPO_ROOT / "tests" / "data"

# bedtools is required by ted.py's TSS/TTS metrics; skip the whole module
# if it isn't on PATH so the suite stays green in lean envs.
pytestmark = pytest.mark.skipif(
    shutil.which("bedtools") is None,
    reason="bedtools required for the empty-input eval-chain end-to-end test",
)


def _resolve_script(name: str) -> str:
    """Map a script name to its absolute path under bin/ or bin/evaluation/."""
    for base in (BIN, EVAL):
        cand = base / name
        if cand.exists():
            return str(cand)
    raise FileNotFoundError(name)


def _run(cmd, cwd, extra_pythonpath):
    """Run a script, returning the CompletedProcess (capture for diagnostics).

    cmd[0] is a bare script name; it is resolved to an absolute path so the
    test is independent of cwd (the Nextflow modules invoke by absolute path).
    """
    import os
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(p) for p in extra_pythonpath] + [env.get("PYTHONPATH", "")]
    )
    script, *rest = cmd
    return subprocess.run(
        [sys.executable, _resolve_script(script), *rest],
        cwd=str(cwd),
        env=env,
        capture_output=True,
        text=True,
    )


@pytest.fixture()
def empty_bed(tmp_path) -> Path:
    p = tmp_path / "empty.bed"
    p.write_text("")
    return p


@pytest.fixture()
def empty_gtf(tmp_path) -> Path:
    p = tmp_path / "empty.gtf"
    p.write_text("")
    return p


def _assert_has_header(path: Path):
    assert path.exists(), f"expected output not written: {path}"
    with open(path) as fh:
        first = fh.readline()
    assert first.strip(), f"output is completely empty (no header): {path}"
    # A real header row has at least one column name.
    assert len(next(csv.reader([first], delimiter="\t"))) >= 1


def test_gtf_to_bed12_empty_exits_zero(empty_gtf, tmp_path):
    """0-transcript GTF must NOT exit non-zero (ted.py runs it check=True)."""
    out = tmp_path / "out.bed"
    cp = _run(["gtf_to_bed12.py", "--gtf", str(empty_gtf), "--output", str(out)],
              cwd=tmp_path, extra_pythonpath=[BIN])
    assert cp.returncode == 0, cp.stderr
    assert out.exists()
    assert out.read_text() == ""  # valid empty BED12


def test_ted_simplified_gtf_path_empty(empty_gtf, tmp_path):
    """ted.py --gtf-input (Bambu/IsoQuant/IsoSeq path) on 0 isoforms."""
    plots = tmp_path / "ted_plots"
    plots.mkdir()
    out = tmp_path / "ted.tsv"
    cp = _run([
        "ted.py", "--gtf-input", str(empty_gtf), "--skip-read-metrics",
        "--genome", str(DATA / "tiny_genome.fa"),
        "--gtf", str(DATA / "tiny_annotation.gtf"),
        "--prime5-peaks", str(DATA / "tiny_cage_peaks.bed"),
        "--prime3-peaks", str(DATA / "tiny_drna_peaks.bed"),
        "--window", "50", "--stage", "transcriptome",
        "--pipeline-mode", "transcriptome_isoseq",
        "--plot-output-dir", str(plots), "--plot-dpi", "72",
        "--output", str(out),
    ], cwd=tmp_path, extra_pythonpath=[BIN])
    assert cp.returncode == 0, cp.stderr
    _assert_has_header(out)


def test_ted_bed_path_empty(empty_bed, tmp_path):
    """ted.py --isoforms-bed (FLAIR path) on a 0-isoform BED, full plotting."""
    plots = tmp_path / "ted_plots"
    plots.mkdir()
    out = tmp_path / "ted.tsv"
    cp = _run([
        "ted.py", "--isoforms-bed", str(empty_bed),
        "--map-file", str(DATA / "tiny_read_map.txt"),
        "--bam", str(DATA / "tiny_reads.bam"),
        "--reads-bed", str(DATA / "tiny_reads.bed"),
        "--genome", str(DATA / "tiny_genome.fa"),
        "--gtf", str(DATA / "tiny_annotation.gtf"),
        "--prime5-peaks", str(DATA / "tiny_cage_peaks.bed"),
        "--prime3-peaks", str(DATA / "tiny_drna_peaks.bed"),
        "--window", "50", "--stage", "transcriptome",
        "--pipeline-mode", "transcriptome_flair_ted",
        "--plot-output-dir", str(plots), "--plot-dpi", "72",
        "--output", str(out),
    ], cwd=tmp_path, extra_pythonpath=[BIN])
    assert cp.returncode == 0, cp.stderr
    _assert_has_header(out)


def test_flair_eval_and_synthesize_empty(empty_gtf, tmp_path):
    """flair_eval.py + synthesize_evaluations.py on 0 isoforms."""
    plots = tmp_path / "ted_plots"
    plots.mkdir()
    ted_out = tmp_path / "ted.tsv"
    # produce a ted.tsv first (simplified path) so synthesize has both inputs
    cp = _run([
        "ted.py", "--gtf-input", str(empty_gtf), "--skip-read-metrics",
        "--genome", str(DATA / "tiny_genome.fa"),
        "--gtf", str(DATA / "tiny_annotation.gtf"),
        "--window", "50", "--stage", "transcriptome",
        "--pipeline-mode", "transcriptome_isoseq",
        "--plot-output-dir", str(plots), "--plot-dpi", "72",
        "--output", str(ted_out),
    ], cwd=tmp_path, extra_pythonpath=[BIN])
    assert cp.returncode == 0, cp.stderr

    flair_out = tmp_path / "flair_eval.tsv"
    cats = tmp_path / "categories.tsv"
    cp = _run([
        "flair_eval.py",
        "--reads-bed", str(DATA / "tiny_reads.bed"),
        "--gtf-input", str(empty_gtf),
        "--gtf", str(DATA / "tiny_annotation.gtf"),
        "--pipeline-mode", "transcriptome_isoseq", "--stage", "transcriptome",
        "--plot-output-dir", str(plots), "--plot-prefix", "WTC11_isoseq",
        "--categories-output", str(cats), "--output", str(flair_out),
    ], cwd=tmp_path, extra_pythonpath=[BIN])
    assert cp.returncode == 0, cp.stderr
    _assert_has_header(flair_out)

    synth_out = tmp_path / "evaluation.tsv"
    cp = _run([
        "synthesize_evaluations.py",
        "--ted-files", str(ted_out), "--flair-files", str(flair_out),
        "--output", str(synth_out), "--test-name", "T",
    ], cwd=tmp_path, extra_pythonpath=[BIN])
    assert cp.returncode == 0, cp.stderr
    _assert_has_header(synth_out)


def test_ted_end_precision_empty_writes_both_tsvs(empty_bed, tmp_path):
    """TedEndPrecision module does an unconditional `mv per_junction_chain.tsv`
    and `mv precision_recall_summary.tsv` — both must always be written."""
    for label, peaks in (("peaks", True), ("gtf_only", False)):
        outdir = tmp_path / label
        cmd = [
            "ted_end_precision.py", "--isoforms-bed", str(empty_bed),
            "--gtf", str(DATA / "tiny_annotation.gtf"),
            "--mode", "isoseq", "--window", "50", "--outdir", str(outdir),
        ]
        if peaks:
            cmd += ["--peaks-5prime", str(DATA / "tiny_cage_peaks.bed"),
                    "--peaks-3prime", str(DATA / "tiny_drna_peaks.bed")]
        cp = _run(cmd, cwd=tmp_path, extra_pythonpath=[EVAL])
        assert cp.returncode == 0, cp.stderr
        _assert_has_header(outdir / "precision_recall_summary.tsv")
        _assert_has_header(outdir / "per_junction_chain.tsv")


def test_firstpass_vs_final_empty_writes_tsv_and_png(empty_bed, tmp_path):
    """FirstpassComparison requires BOTH firstpass_vs_final.tsv and .png as
    declared outputs — they must exist even when both BEDs are empty."""
    outdir = tmp_path / "fvf"
    cp = _run([
        "firstpass_vs_final.py",
        "--firstpass-bed", str(empty_bed), "--final-bed", str(empty_bed),
        "--gtf", str(DATA / "tiny_annotation.gtf"),
        "--peaks-5prime", str(DATA / "tiny_cage_peaks.bed"),
        "--peaks-3prime", str(DATA / "tiny_drna_peaks.bed"),
        "--mode", "isoseq", "--window", "50", "--outdir", str(outdir),
    ], cwd=tmp_path, extra_pythonpath=[EVAL])
    assert cp.returncode == 0, cp.stderr
    _assert_has_header(outdir / "firstpass_vs_final.tsv")
    assert (outdir / "firstpass_vs_final.png").exists(), \
        "PNG output missing — Nextflow declares it a required output"
