"""phoenix (UCSC PRISM) helpers: paths, the driver job script, sbatch submission, presets."""
import subprocess
from pathlib import Path

import pytest

from flair_eval import phoenix
from flair_eval import validate as V


def test_absolute_paths(tmp_path):
    base = str(tmp_path)
    assert phoenix.absolute("r.fq.gz", base) == f"{base}/r.fq.gz"
    assert phoenix.absolute("../x/r.fq.gz", base) == f"{tmp_path.parent}/x/r.fq.gz"
    assert phoenix.absolute("/abs/r.fq.gz", base) == "/abs/r.fq.gz"
    assert phoenix.absolute("hg38", base) == "hg38"  # a genome preset is not a path
    assert phoenix.absolute("s3://bucket/r.fq.gz", base) == "s3://bucket/r.fq.gz"
    assert phoenix.absolute("", base) == "" and phoenix.absolute(None, base) is None
    assert phoenix.absolute_list("a,/b, c/*", base) == f"{base}/a,/b,{base}/c/*"


def test_absolute_expands_user_and_vars(monkeypatch):
    monkeypatch.setenv("USER", "someone")
    monkeypatch.setenv("HOME", "/home/someone")
    assert phoenix.absolute("/groups/$USER/out") == "/groups/someone/out"
    assert phoenix.absolute("~/x") == "/home/someone/x"


def test_in_home_and_node_local(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert phoenix.in_home("out", home=str(tmp_path))
    assert phoenix.in_home(str(tmp_path), home=str(tmp_path))
    assert not phoenix.in_home("/elsewhere/out", home=str(tmp_path))
    assert not phoenix.in_home(str(tmp_path) + "-other/out", home=str(tmp_path))
    assert phoenix.node_local("/data/tmp/me/out") and phoenix.node_local("/tmp")
    assert not phoenix.node_local("/private/groups/brookslab/me/out")


def test_sbatch_script():
    import datetime
    text = phoenix.sbatch_script("nextflow run hdheath/flair-eval -r v2 \\\n  --outdir /g/o",
                                 "/g/my out", version="9.9",
                                 now=datetime.datetime(2026, 11, 10, 9, 30))
    lines = text.splitlines()
    assert lines[0] == "#!/bin/bash"
    assert lines[1:7] == ["#SBATCH --job-name=flair-eval", "#SBATCH --partition=long",
                          "#SBATCH --cpus-per-task=2", "#SBATCH --mem=8G",
                          "#SBATCH --time=7-00:00:00",
                          "#SBATCH --output=/g/my out/pipeline_info/nextflow-%j.log"]
    assert "# Written by flair-eval 9.9 on 2026-11-10 09:30." in lines
    assert "cd '/g/my out'" in lines  # quoted for the shell
    # shared lab caches (env_cache, conda.cacheDir) need group-writable lock files
    assert "umask 002" in lines and lines.index("umask 002") < lines.index("cd '/g/my out'")
    assert text.endswith("  --outdir /g/o\n")
    report = phoenix.sbatch_script("x", "/o", phoenix.REPORT_JOB)
    assert "#SBATCH --partition=medium" in report and "#SBATCH --time=12:00:00" in report


def test_sbatch_script_is_valid_bash(tmp_path):
    p = tmp_path / "launch.sbatch"
    p.write_text(phoenix.sbatch_script("nextflow run x \\\n  --outdir /o", "/o"))
    assert subprocess.run(["bash", "-n", str(p)]).returncode == 0


def test_submit_parses_job_id():
    seen = []

    def run(argv, **k):
        seen.append(argv)
        return subprocess.CompletedProcess(argv, 0, stdout="812345;phoenix\n", stderr="")
    assert phoenix.submit(Path("/o/launch.sbatch"), run=run) == "812345"
    assert seen == [["sbatch", "--parsable", "/o/launch.sbatch"]]


def test_submit_errors():
    def failing(argv, **k):
        return subprocess.CompletedProcess(argv, 1, stdout="", stderr="invalid partition\n")
    with pytest.raises(RuntimeError, match="invalid partition"):
        phoenix.submit(Path("x"), run=failing)

    def missing(argv, **k):
        raise FileNotFoundError("sbatch")
    with pytest.raises(RuntimeError, match="sbatch failed"):
        phoenix.submit(Path("x"), run=missing)


# ---------------------------------------------------------------- validation with presets
def test_validate_accepts_preset_without_gtf():
    assert V.validate_samples([{"sample": "s", "reads": "r", "genome": "hg38"}]) == []
    got = [(i.path, i.code) for i in V.validate_samples([{"sample": "s", "reads": "r",
                                                           "genome": "/ref/g.fa"}])]
    assert got == [(("samples", 0, "gtf"), "required")]


def test_check_paths_preset(monkeypatch, tmp_path):
    fa, gtf = tmp_path / "g.fa", tmp_path / "a.gtf"
    fa.write_text("")
    gtf.write_text("")
    monkeypatch.setitem(phoenix.GENOME_PRESETS, "hg38",
                        {"fasta": str(fa), "gtf": str(gtf), "note": "test"})
    rows = [{"sample": "s", "reads": str(fa), "genome": "hg38"}]
    assert V.check_paths(rows) == []  # no "file not found: hg38"
    gtf.unlink()
    got = V.check_paths(rows)
    assert [(i.path, i.code, i.level) for i in got] == [
        (("samples", 0, "genome"), "preset_missing", "warning")]
    assert "not on this machine" in got[0].message
    # an own GTF replaces the preset's, so its absence no longer matters
    assert V.check_paths([dict(rows[0], gtf=str(fa))]) == []


def test_cwd_only_paths(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = tmp_path / "cfg"
    cfg.mkdir()
    (tmp_path / "r.fq.gz").write_text("")
    (cfg / "a.gtf").write_text("")
    rows = [{"sample": "s", "reads": "r.fq.gz", "genome": "hg38", "gtf": "a.gtf"}]
    got = V.cwd_only_paths(rows, cfg)
    assert [(i.path, i.code) for i in got] == [(("samples", 0, "reads"), "cwd_relative")]
