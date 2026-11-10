import json

import pytest
import yaml

from flair_eval import cli, env_check
from flair_eval.config import runs_yaml, samples_csv

from .test_validate import BUNDLE


@pytest.fixture
def here(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    for f in ("g.fa", "a.gtf", "r.fq.gz"):
        (tmp_path / f).write_text("")
    return tmp_path


def lines_of(text):
    return [x.strip().rstrip("\\").strip() for x in text.splitlines()]


def test_help_and_version(capsys):
    assert cli.main([]) == 0
    assert "init" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        cli.main(["--version"])


def test_init_from_flags_then_check(here, capsys):
    assert cli.main(["init", "--reads", "r.fq.gz", "--genome", "g.fa", "--gtf", "a.gtf",
                     "--flair", "3.0.1,dev", "--flair_args", "--stringent"]) == 0
    out = capsys.readouterr().out
    assert "Wrote samples.csv (1 sample) and runs.yaml (2 runs)" in out
    assert "--flair_args=--stringent" in out
    runs = yaml.safe_load((here / "runs.yaml").read_text())["runs"]
    assert [r["transcriptome_args"] for r in runs] == ["--stringent"] * 2
    assert cli.main(["check"]) == 0
    assert "All good." in capsys.readouterr().out
    # refuses to overwrite
    assert cli.main(["init", "--reads", "r.fq.gz", "--genome", "g.fa", "--gtf", "a.gtf",
                     "--flair", "3.0.1"]) == 1
    assert "--force" in capsys.readouterr().err


def test_init_errors(here, capsys):
    assert cli.main(["init", "--reads", "r.fq.gz", "--flair", "3.0.1"]) == 1  # no genome/gtf
    assert "missing required column 'genome'" in capsys.readouterr().err
    assert cli.main(["init", "--reads", "r.fq.gz", "--genome", "g", "--gtf", "a",
                     "--flair", "master"]) == 2
    assert cli.main(["init", "--flair", "3.0.1"]) == 2


def test_init_interactive(here, monkeypatch, capsys):
    answers = iter(["r.fq.gz", "", "g.fa", "a.gtf", "", "", "n", "3.0.1", "", "n"])
    monkeypatch.setattr("builtins.input", lambda p="": next(answers))
    assert cli.main(["init", "--dir", "cfg"]) == 0
    assert (here / "cfg" / "samples.csv").exists()
    out = capsys.readouterr().out
    assert "flair-eval check --samples cfg/samples.csv --runs cfg/runs.yaml" in out


def test_check_reports_errors_and_strict(here, capsys):
    (here / "samples.csv").write_text("sample,reads,genome,gtf\ns 1,missing.fq,g.fa,a.gtf\n")
    (here / "runs.yaml").write_text("runs:\n  - id: a\n    flair: {release: 3.0}\n")
    assert cli.main(["check"]) == 1
    out = capsys.readouterr().out
    assert "samples.csv row 1: sample: 's 1' is not a valid name" in out
    assert 'put it in quotes, e.g. release: "3.0"' in out
    assert "warning: samples.csv row 1: reads: file not found: missing.fq" in out
    (here / "samples.csv").write_text("sample,reads,genome,gtf\ns1,missing.fq,g.fa,a.gtf\n")
    (here / "runs.yaml").write_text('runs:\n  - id: a\n    flair: {release: "3.0"}\n')
    assert cli.main(["check"]) == 0
    assert cli.main(["check", "--strict"]) == 1
    assert cli.main(["check", "--no-files", "--strict"]) == 0


def test_check_nothing(here, capsys):
    assert cli.main(["check"]) == 2


def test_check_bundles(here, capsys):
    d = here / "out" / "bundles" / "s" / "r"
    d.mkdir(parents=True)
    (d / "bundle.json").write_text(json.dumps(BUNDLE))
    assert cli.main(["check", "--bundles", "out"]) == 0
    bad = dict(BUNDLE, schema="nope")
    (d / "bundle.json").write_text(json.dumps(bad))
    assert cli.main(["check", "--bundles", "out"]) == 1
    assert "should be 'flair-eval/bundle/1'" in capsys.readouterr().out
    assert cli.main(["check", "--bundles", "empty-dir"]) == 1


def test_run_quick_dry_run(here, capsys):
    rc = cli.main(["run", "--reads", "r.fq.gz", "--genome", "g.fa", "--gtf", "a.gtf",
                   "--flair", "3.0.1,dev", "--flair_args", "--stringent --check_splice",
                   "--profile", "slurm", "--outdir", "/big/out", "--baseline", "flair-dev",
                   "--dry-run", "--", "-with-tower"])
    assert rc == 0
    got = lines_of(capsys.readouterr().out)
    assert "nextflow run hdheath/flair-eval -r v2" in got
    assert "--flair_args='--stringent --check_splice'" in got
    assert "--sample r" in got and "-profile slurm" in got and "-with-tower" in got


def test_run_files_dry_run_picks_up_cwd_files(here, capsys, sample, two_runs):
    (here / "samples.csv").write_text(samples_csv([dict(sample, reads="r.fq.gz",
                                                        genome="g.fa", gtf="a.gtf")]))
    (here / "runs.yaml").write_text(runs_yaml(two_runs))
    assert cli.main(["run", "--dry-run", "-r", "1.2.3", "--resume"]) == 0
    got = lines_of(capsys.readouterr().out)
    assert "nextflow run hdheath/flair-eval -r 1.2.3" in got
    assert "--input samples.csv" in got and "--runs runs.yaml" in got and "-resume" in got


def test_run_mixed_sources(here, capsys, two_runs):
    (here / "runs.yaml").write_text(runs_yaml(two_runs))
    assert cli.main(["run", "--reads", "r.fq.gz", "--genome", "g.fa", "--gtf", "a.gtf",
                     "--runs", "runs.yaml", "--dry-run"]) == 0
    got = lines_of(capsys.readouterr().out)
    assert "--reads r.fq.gz" in got and "--runs runs.yaml" in got
    assert not any(x.startswith("--input") for x in got)


def test_run_refuses_bad_input(here, capsys, two_runs):
    assert cli.main(["run", "--reads", "r.fq.gz", "--genome", "g.fa", "--gtf", "a.gtf",
                     "--flair", "3.0.1", "--baseline", "nope", "--dry-run"]) == 1
    assert "--baseline 'nope' is not a run id (flair-3.0.1)" in capsys.readouterr().err
    assert cli.main(["run", "--reads", "r.fq.gz", "--genome", "g.fa", "--gtf", "a.gtf",
                     "--flair", "master", "--dry-run"]) == 1
    assert cli.main(["run", "--dry-run", "--flair", "3.0.1"]) == 2  # no samples
    assert cli.main(["run", "--dry-run", "--reads", "r.fq.gz", "--genome", "g.fa",
                     "--gtf", "a.gtf"]) == 2  # no runs
    (here / "s.csv").write_text("sample,reads,genome,gtf\ns,r.fq.gz,g.fa,a.gtf\n")
    assert cli.main(["run", "--samples", "s.csv", "--reads", "r.fq.gz", "--flair", "3.0.1",
                     "--dry-run"]) == 2
    (here / "r.yaml").write_text(runs_yaml(two_runs))
    assert cli.main(["run", "--samples", "s.csv", "--runs", "r.yaml", "--flair", "3.0.1",
                     "--dry-run"]) == 2


def test_run_checks_tools_and_launches(here, monkeypatch, capsys):
    args = ["run", "--reads", "r.fq.gz", "--genome", "g.fa", "--gtf", "a.gtf", "--flair", "3.0.1"]
    monkeypatch.setattr(env_check, "check_tools", lambda: (["nextflow is not on PATH"], []))
    assert cli.main(args) == 2
    assert "nextflow is not on PATH" in capsys.readouterr().err
    calls = []
    monkeypatch.setattr(env_check, "check_tools", lambda: ([], ["old java"]))
    monkeypatch.setattr(cli.subprocess, "call", lambda argv: calls.append(argv) or 0)
    assert cli.main(args) == 0
    assert calls[0][:5] == ["nextflow", "run", "hdheath/flair-eval", "-r", "v2"]
    assert cli.main(args + ["--skip-checks"]) == 0 and len(calls) == 2


def test_report(here, capsys):
    assert cli.main(["report", "out1", "out2", "--title", "Nov runs", "--baseline", "flair-3.0.1",
                     "--outdir", "rep", "--dry-run"]) == 0
    got = lines_of(capsys.readouterr().out)
    assert "--report_only true" in got and "--bundles out1,out2" in got
    assert "--report_title 'Nov runs'" in got and "--baseline flair-3.0.1" in got
    assert "--outdir rep" in got


def test_extra_args_only_for_run_and_report(here):
    with pytest.raises(SystemExit):
        cli.main(["check", "--", "-x"])
