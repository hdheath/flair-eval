import json

import pytest
import yaml

from flair_eval import cli, env_check
from flair_eval.config import runs_yaml, samples_csv

from .test_validate import BUNDLE


@pytest.fixture(autouse=True)
def fixed_profile(monkeypatch):
    """Profile and phoenix detection depend on the machine; pin them (not phoenix)."""
    monkeypatch.setattr(cli.doctor, "suggest_profile", lambda which=None: "local,conda")
    monkeypatch.setattr(cli.phoenix, "detect", lambda **k: False)
    monkeypatch.setattr(cli.phoenix, "submit", _no_sbatch)


def _no_sbatch(*a, **k):  # pragma: no cover - a test reaching this is a bug
    raise AssertionError("tests must not submit SLURM jobs")


@pytest.fixture
def on_phoenix(monkeypatch):
    monkeypatch.setattr(cli.phoenix, "detect", lambda **k: True)
    monkeypatch.setattr(cli.doctor, "suggest_profile", lambda which=None: "")


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
    answers = iter(["n", "r.fq.gz", "", "g.fa", "a.gtf", "", "", "n", "3.0.1", "", "n",
                    "/groups/out"])
    monkeypatch.setattr("builtins.input", lambda p="": next(answers))
    assert cli.main(["init", "--dir", "cfg"]) == 0
    assert (here / "cfg" / "samples.csv").exists()
    out = capsys.readouterr().out
    assert "flair-eval check --samples cfg/samples.csv --runs cfg/runs.yaml" in out
    assert "flair-eval run --samples cfg/samples.csv --runs cfg/runs.yaml --outdir /groups/out" \
        in out


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
    assert "nextflow run hdheath/flair-eval -r main" in got
    assert "--flair_args='--stringent --check_splice'" in got
    assert "--sample r" in got and "-profile slurm" in got and "-with-tower" in got


def test_run_files_dry_run_picks_up_cwd_files(here, capsys, sample, two_runs):
    (here / "samples.csv").write_text(samples_csv([dict(sample, reads="r.fq.gz",
                                                        genome="g.fa", gtf="a.gtf")]))
    (here / "runs.yaml").write_text(runs_yaml(two_runs))
    assert cli.main(["run", "--dry-run", "-r", "1.2.3", "--resume", "--outdir", "out"]) == 0
    got = lines_of(capsys.readouterr().out)
    assert "nextflow run hdheath/flair-eval -r 1.2.3" in got
    assert "--input samples.csv" in got and "--runs runs.yaml" in got and "-resume" in got


def test_run_mixed_sources(here, capsys, two_runs):
    (here / "runs.yaml").write_text(runs_yaml(two_runs))
    assert cli.main(["run", "--reads", "r.fq.gz", "--genome", "g.fa", "--gtf", "a.gtf",
                     "--runs", "runs.yaml", "--dry-run", "--outdir", "out"]) == 0
    got = lines_of(capsys.readouterr().out)
    assert "--reads r.fq.gz" in got and "--runs runs.yaml" in got
    assert not any(x.startswith("--input") for x in got)


def test_run_refuses_bad_input(here, capsys, two_runs):
    assert cli.main(["run", "--reads", "r.fq.gz", "--genome", "g.fa", "--gtf", "a.gtf",
                     "--flair", "3.0.1", "--baseline", "nope", "--dry-run", "--outdir", "o"]) == 1
    assert "--baseline 'nope' is not a run id (flair-3.0.1)" in capsys.readouterr().err
    assert cli.main(["run", "--reads", "r.fq.gz", "--genome", "g.fa", "--gtf", "a.gtf",
                     "--flair", "master", "--dry-run", "--outdir", "o"]) == 1
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
    args = ["run", "--reads", "r.fq.gz", "--genome", "g.fa", "--gtf", "a.gtf", "--flair", "3.0.1",
            "--outdir", "out"]
    monkeypatch.setattr(env_check, "check_tools", lambda: (["nextflow is not on PATH"], []))
    assert cli.main(args) == 2
    assert "nextflow is not on PATH" in capsys.readouterr().err
    calls = []
    monkeypatch.setattr(env_check, "check_tools", lambda: ([], ["old java"]))
    monkeypatch.setattr(cli.subprocess, "call", lambda argv: calls.append(argv) or 0)
    assert cli.main(args) == 0
    assert calls[0][:5] == ["nextflow", "run", "hdheath/flair-eval", "-r", "main"]
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


def test_init_test_data_choice(here, monkeypatch, capsys):
    monkeypatch.setattr("builtins.input", lambda p="": "y")
    assert cli.main(["init"]) == 0
    assert "flair-eval run --test" in capsys.readouterr().out
    assert not (here / "samples.csv").exists()


def test_auto_profile(here, capsys, sample, two_runs):
    (here / "samples.csv").write_text(samples_csv([dict(sample, reads="r.fq.gz",
                                                        genome="g.fa", gtf="a.gtf")]))
    (here / "runs.yaml").write_text(runs_yaml(two_runs))
    assert cli.main(["run", "--dry-run", "--outdir", "/big/out"]) == 0
    out = capsys.readouterr().out
    assert "Using -profile local,conda (picked for this machine" in out
    assert "-profile local,conda" in lines_of(out)
    assert "--outdir /big/out" in lines_of(out)


@pytest.mark.parametrize("argv", [
    ["run", "--reads", "r.fq.gz", "--genome", "g.fa", "--gtf", "a.gtf", "--flair", "3.0.1"],
    ["run", "--test"],
    ["report", "out1"],
])
def test_outdir_is_required(here, capsys, argv):
    from flair_eval import phoenix
    assert cli.main(argv + ["--dry-run"]) == 2
    err = capsys.readouterr().err
    assert "--outdir is required" in err and phoenix.OUTDIR_EXAMPLE in err


def test_outdir_in_home_warns(here, monkeypatch, capsys):
    monkeypatch.setenv("HOME", str(here))
    assert cli.main(["run", "--test", "--dry-run", "--outdir", "results"]) == 0
    assert "is in your home folder" in capsys.readouterr().err
    assert cli.main(["run", "--test", "--dry-run", "--outdir", "/elsewhere/results"]) == 0
    assert "home folder" not in capsys.readouterr().err


def test_run_test_data(here, capsys):
    assert cli.main(["run", "--test", "--dry-run", "--outdir", "t"]) == 0
    got = lines_of(capsys.readouterr().out)
    assert "-profile test,local,conda" in got and "--outdir t" in got
    assert not any(x.startswith("--input") or x.startswith("--runs") for x in got)
    assert cli.main(["run", "--test", "--profile", "slurm", "--dry-run", "--outdir", "t"]) == 0
    assert "-profile test,slurm" in lines_of(capsys.readouterr().out)
    assert cli.main(["run", "--test", "--flair", "3.0.1", "--dry-run", "--outdir", "t"]) == 2


def test_doctor_command(here, monkeypatch, capsys):
    from flair_eval import doctor
    monkeypatch.setattr(doctor, "run_checks", lambda outdir=None: [
        doctor.Check("ok", "Java", "version 17"),
        doctor.Check("fail", "Nextflow", "not found", "conda install -c bioconda nextflow")])
    assert cli.main(["doctor"]) == 1
    out = capsys.readouterr().out
    assert "FAIL  Nextflow: not found" in out and "conda install -c bioconda nextflow" in out
    assert "-profile local,conda" in out


def test_run_local_checkout(here, capsys):
    assert cli.main(["run", "--pipeline", "./flair-eval", "--reads", "r.fq.gz", "--genome", "g.fa",
                     "--gtf", "a.gtf", "--flair", "3.0.1", "--dry-run", "--outdir", "o",
                     "--", "-stub"]) == 0
    got = lines_of(capsys.readouterr().out)
    assert "nextflow run ./flair-eval" in got and "-stub" in got


# ---------------------------------------------------------------- genome presets
def test_genome_preset_quick_run(here, capsys):
    assert cli.main(["run", "--reads", "r.fq.gz", "--genome", "hg38", "--flair", "2.2.0,3.0.1,dev",
                     "--outdir", "/big/out", "--dry-run"]) == 0
    cap = capsys.readouterr()
    got = lines_of(cap.out)
    assert "--genome hg38" in got and not any(x.startswith("--gtf") for x in got)
    assert "file not found" not in cap.err


def test_genome_preset_in_samples_csv(here, capsys, two_runs):
    (here / "samples.csv").write_text("sample,reads,genome,gtf\ns1,r.fq.gz,hg38,\n")
    (here / "runs.yaml").write_text(runs_yaml(two_runs))
    assert cli.main(["check"]) == 0
    out = capsys.readouterr().out
    assert "missing required column 'gtf'" not in out and "file not found: hg38" not in out
    # without a preset the gtf is still required
    (here / "samples.csv").write_text("sample,reads,genome,gtf\ns1,r.fq.gz,g.fa,\n")
    assert cli.main(["check"]) == 1
    assert "missing required column 'gtf'" in capsys.readouterr().out


def test_init_with_genome_preset(here, capsys):
    assert cli.main(["init", "--reads", "r.fq.gz", "--genome", "hg38", "--flair", "3.0.1",
                     "--outdir", "/groups/x"]) == 0
    out = capsys.readouterr().out
    assert "--outdir /groups/x" in out
    assert (here / "samples.csv").read_text().splitlines()[1] == "r,r.fq.gz,,hg38,,,"


def test_cli_genome_presets_match_nextflow_config(repo):
    """flair_eval.phoenix.GENOME_PRESETS must equal params.genomes in nextflow.config."""
    import re
    from flair_eval import phoenix
    p = repo / "nextflow.config"
    if not p.exists():
        pytest.skip("not in the pipeline repository")
    text = p.read_text()
    m = re.search(r"^\s*genomes\s*=\s*\[", text, re.M)
    assert m, "params.genomes not found in nextflow.config"
    depth, i = 1, m.end()
    while depth:  # the matching ] of `genomes = [`
        depth += {"[": 1, "]": -1}.get(text[i], 0)
        i += 1
    block = text[m.end():i - 1]
    nf = {}
    for key, body in re.findall(r"^\s*([A-Za-z0-9_.-]+)\s*:\s*\[(.*?)\]", block, re.S | re.M):
        nf[key] = dict(re.findall(r"(\w+)\s*:\s*'([^']*)'", body))
    assert sorted(nf) == sorted(phoenix.GENOME_PRESETS)
    for key, preset in phoenix.GENOME_PRESETS.items():
        assert nf[key]["fasta"] == preset["fasta"], key
        assert nf[key]["gtf"] == preset["gtf"], key


# ---------------------------------------------------------------- phoenix launch
def test_phoenix_run_dry_run_prints_sbatch_script(here, on_phoenix, capsys):
    out_abs = str(here / "res")
    assert cli.main(["run", "--reads", "r.fq.gz", "--genome", "hg38", "--flair", "3.0.1,dev",
                     "--outdir", "res", "--regions", "x.bed", "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "phoenix detected" in out
    assert f"would be {out_abs}/pipeline_info/launch.sbatch" in out
    for line in ("#SBATCH --job-name=flair-eval", "#SBATCH --partition=long",
                 "#SBATCH --cpus-per-task=2", "#SBATCH --mem=8G", "#SBATCH --time=7-00:00:00",
                 f"#SBATCH --output={out_abs}/pipeline_info/nextflow-%j.log", f"cd {out_abs}"):
        assert line in out.splitlines(), line
    got = lines_of(out)
    assert f"--reads {here}/r.fq.gz" in got and "--genome hg38" in got  # preset kept as is
    assert f"--outdir {out_abs}" in got and f"--regions {here}/x.bed" in got
    assert not any(x.startswith("-profile") or " -profile " in x for x in got)
    assert not (here / "res").exists()  # dry run writes nothing


def test_phoenix_run_submits(here, on_phoenix, monkeypatch, capsys, two_runs):
    (here / "samples.csv").write_text("sample,reads,genome,gtf\ns1,r.fq.gz,g.fa,a.gtf\n")
    (here / "runs.yaml").write_text(runs_yaml(two_runs))
    monkeypatch.setattr(env_check, "check_tools", lambda: ([], []))
    submitted = []
    monkeypatch.setattr(cli.phoenix, "submit", lambda path: submitted.append(path) or "4242")
    monkeypatch.setattr(cli.subprocess, "call", _no_sbatch)
    assert cli.main(["run", "--outdir", "res", "--resume"]) == 0
    script = here / "res" / "pipeline_info" / "launch.sbatch"
    assert submitted == [script]
    text = script.read_text()
    assert text.startswith("#!/bin/bash\n#SBATCH")
    assert f"umask 002\ncd {here}/res\n" in text
    got = lines_of(text)
    assert f"--input {here}/samples.csv" in got and f"--runs {here}/runs.yaml" in got
    assert "-resume" in got
    out = capsys.readouterr().out
    assert "SLURM job 4242" in out and "squeue -j 4242" in out
    assert f"{here}/res/pipeline_info/nextflow-4242.log" in out
    assert f"{here}/res/report/report.html" in out


def test_phoenix_sbatch_error(here, on_phoenix, monkeypatch, capsys):
    monkeypatch.setattr(env_check, "check_tools", lambda: ([], []))

    def fail(path):
        raise RuntimeError("sbatch failed: invalid partition")
    monkeypatch.setattr(cli.phoenix, "submit", fail)
    assert cli.main(["run", "--test", "--outdir", "res"]) == 1
    err = capsys.readouterr().err
    assert "invalid partition" in err and "launch.sbatch" in err


def test_phoenix_foreground(here, on_phoenix, monkeypatch, capsys):
    calls = []
    monkeypatch.setattr(env_check, "check_tools", lambda: ([], []))
    monkeypatch.setattr(cli.subprocess, "call", lambda argv: calls.append(argv) or 0)
    assert cli.main(["run", "--reads", "r.fq.gz", "--genome", "hg38", "--flair", "3.0.1",
                     "--outdir", "res", "--foreground"]) == 0
    assert calls and "-profile" not in calls[0]
    assert calls[0][calls[0].index("--reads") + 1] == "r.fq.gz"  # as typed, no job script
    assert not (here / "res").exists()


def test_phoenix_test_profile_and_report(here, on_phoenix, capsys):
    assert cli.main(["run", "--test", "--outdir", "t", "--dry-run"]) == 0
    assert "-profile test,phoenix" in lines_of(capsys.readouterr().out)
    assert cli.main(["report", "b1", "--outdir", "rep", "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "#SBATCH --partition=medium" in out
    assert f"--bundles {here}/b1" in lines_of(out)


def test_phoenix_cwd_relative_paths_in_samples_csv(here, on_phoenix, capsys, two_runs):
    (here / "cfg").mkdir()
    (here / "cfg" / "samples.csv").write_text("sample,reads,genome,gtf\ns1,r.fq.gz,g.fa,a.gtf\n")
    (here / "runs.yaml").write_text(runs_yaml(two_runs))
    assert cli.main(["run", "--samples", "cfg/samples.csv", "--outdir", "o", "--dry-run"]) == 1
    assert "the SLURM job launches in --outdir" in capsys.readouterr().err
    assert cli.main(["run", "--samples", "cfg/samples.csv", "--outdir", "o", "--dry-run",
                     "--foreground"]) == 0


def test_phoenix_doctor(here, on_phoenix, monkeypatch, capsys):
    from flair_eval import doctor
    monkeypatch.setattr(doctor, "run_checks", lambda outdir=None: [
        doctor.Check("ok", "Java", "version 17")])
    assert cli.main(["doctor"]) == 0
    out = capsys.readouterr().out
    assert "phoenix (UCSC PRISM) detected: no -profile needed" in out
    assert "shared lab cache" in out and "-profile local" not in out


def test_phoenix_init_wizard_defaults(here, on_phoenix, monkeypatch, capsys):
    answers = iter(["n", "r.fq.gz", "", "", "", "", "", "n", "3.0.1", "", "n", ""])
    monkeypatch.setattr("builtins.input", lambda p="": next(answers))
    monkeypatch.setattr(cli, "suggested_outdir", lambda: "/groups/u/flair-eval/10_11_26_x")
    assert cli.main(["init"]) == 0
    assert "r,r.fq.gz,,hg38,,," in (here / "samples.csv").read_text()
    out = capsys.readouterr().out
    assert "flair-eval run --outdir /groups/u/flair-eval/10_11_26_x" in out
    assert "in a SLURM job, not on the login node" in out
