import shutil
import subprocess

import pytest

from flair_eval.command import Options, build_command, format_command, opt, quote, \
    report_command


def test_quick_form(sample, two_runs):
    c = build_command([sample], two_runs)
    assert c.mode == "quick"
    assert c.argv == ["nextflow", "run", "hdheath/flair-eval", "-r", "main",
                      "--reads", "/data/HG002.fq.gz", "--genome", "/ref/hg38.fa",
                      "--gtf", "/ref/gencode.gtf", "--sample", "HG002",
                      "--flair", "3.0.1,dev", "--outdir", "<outdir>"]


def test_quick_form_genome_preset_without_gtf(sample, two_runs):
    s = {"sample": "HG002", "reads": "/data/HG002.fq.gz", "genome": "hg38"}
    c = build_command([s], two_runs, Options(outdir="/big/out"))
    assert c.mode == "quick"
    assert c.argv[5:11] == ["--reads", "/data/HG002.fq.gz", "--genome", "hg38",
                            "--sample", "HG002"]
    assert "--gtf" not in c.argv


def test_quick_form_with_args_bam_platform_and_options(sample, two_runs):
    s = {"sample": "s1", "bam": "/b/s1.bam", "genome": "g.fa", "gtf": "a.gtf",
         "platform": "pacbio", "library": "cDNA"}
    runs = [dict(r, transcriptome_args="--stringent --check_splice") for r in two_runs]
    opts = Options(ref="2.0.0", profile="slurm", outdir="/big/out", resume=True,
                   work_dir="/scratch/work",
                   params={"baseline": "flair-dev", "report_title": "My runs",
                           "align_flair": "3.0.1", "regions": None, "bundles": ""})
    c = build_command([s], runs, opts)
    assert c.mode == "quick"
    assert c.argv[3:7] == ["-r", "2.0.0", "-profile", "slurm"]
    assert "--flair_args=--stringent --check_splice" in c.argv
    assert c.argv[-9:] == ["--baseline", "flair-dev", "--report_title", "My runs",
                           "--outdir", "/big/out", "-work-dir", "/scratch/work", "-resume"]
    assert "--align_flair" not in c.argv  # equals the pipeline default
    assert c.argv[-1] == "-resume"
    i = c.argv.index("--platform")
    assert c.argv[i - 2:i + 4] == ["--sample", "s1", "--platform", "pacbio", "--library", "cDNA"]


@pytest.mark.parametrize("change,reason", [
    (lambda s, r: ([s, dict(s, sample="s2")], r), "2 samples"),
    (lambda s, r: ([dict(s, cage="c.bed")], r), "extra tracks (cage)"),
    (lambda s, r: ([dict(s, bam="x.bam")], r), "both reads and bam"),
    (lambda s, r: ([s], [dict(r[0], id="mine"), r[1]]), "run mine: its id"),
    (lambda s, r: ([s], [dict(r[0], label="FLAIR 3"), r[1]]), "labels need runs.yaml"),
    (lambda s, r: ([s], [dict(r[0], tags={"v": "3"}), r[1]]), "tags need runs.yaml"),
    (lambda s, r: ([s], [dict(r[0], quantify=False), r[1]]), "quantify is not the default"),
    (lambda s, r: ([s], [dict(r[0], own_align=True), r[1]]), "own_align is not the default"),
    (lambda s, r: ([s], [dict(r[0], transcriptome_args="--stringent"), r[1]]),
     "different transcriptome_args"),
    (lambda s, r: ([s], []), "no runs"),
])
def test_files_form_reasons(sample, two_runs, change, reason):
    samples, runs = change(sample, two_runs)
    c = build_command(samples, runs)
    assert c.mode == "files"
    assert any(reason in x for x in c.reasons), c.reasons
    assert c.argv[5:9] == ["--input", "samples.csv", "--runs", "runs.yaml"]


def test_label_equal_to_id_still_quick(sample, two_runs):
    runs = [dict(two_runs[0], label="flair-3.0.1"), two_runs[1]]
    assert build_command([sample], runs).mode == "quick"


def test_prefer_files(sample, two_runs):
    c = build_command([sample], two_runs, Options(prefer_files=True, samples_path="c/s.csv",
                                                  runs_path="c/r.yaml"))
    assert c.mode == "files" and c.reasons == ["files requested"]
    assert c.argv[5:9] == ["--input", "c/s.csv", "--runs", "c/r.yaml"]


def test_report_command():
    c = report_command(["out1", "out2/bundles"], Options(params={"report_title": "T"}))
    assert c.argv == ["nextflow", "run", "hdheath/flair-eval", "-r", "main", "--report_only",
                      "true", "--bundles", "out1,out2/bundles", "--report_title", "T",
                      "--outdir", "<outdir>"]


@pytest.mark.parametrize("arg,quoted", [
    ("/a/b.fq.gz", "/a/b.fq.gz"), ("a b", "'a b'"), ("it's", "'it'\"'\"'s'"), ("", "''"),
    ("--flair_args=--x --y", "--flair_args='--x --y'"), ("$HOME", "'$HOME'"),
    ("a\n", "'a\n'"),
])
def test_quote(arg, quoted):
    assert quote(arg) == quoted


def test_opt_uses_equals_for_dash_values():
    assert opt("--flair_args", "--stringent") == ["--flair_args=--stringent"]
    assert opt("--report_title", "x") == ["--report_title", "x"]
    assert opt("--baseline", 5) == ["--baseline", "5"]


def test_format_command(sample, two_runs):
    runs = [dict(r, transcriptome_args="--stringent") for r in two_runs]
    text = format_command(build_command([sample], runs, Options(resume=True)).argv)
    lines = text.split(" \\\n")
    assert lines[0] == "nextflow run hdheath/flair-eval -r main"
    assert "  --flair_args=--stringent" in lines
    assert "  --outdir '<outdir>'" in lines  # no default: the pipeline requires --outdir
    assert lines[-1] == "  -resume"


def test_shell_round_trip(sample, two_runs):
    """The one-line and multi-line forms parse back to the same argv in a real shell."""
    import shlex
    runs = [dict(r, transcriptome_args="--support 5 --stringent") for r in two_runs]
    s = dict(sample, reads="/data/my reads/it's.fq.gz", sample="s1")
    c = build_command([s], runs)
    assert shlex.split(str(c)) == c.argv
    assert shlex.split(c.pretty().replace(" \\\n", " ")) == c.argv
    if shutil.which("bash"):
        out = subprocess.run(["bash", "-c", 'printf "%s\\n" ' + c.pretty()[len("nextflow"):]],
                             capture_output=True, text=True).stdout.splitlines()
        assert out == c.argv[1:]


@pytest.mark.skipif(shutil.which("nextflow") is None, reason="nextflow not on PATH")
def test_nextflow_reads_flair_args_value(tmp_path):
    """Guard for the `--flair_args --stringent` trap: Nextflow must see the value."""
    nf = tmp_path / "t.nf"
    nf.write_text('params.flair_args = ""\nworkflow { println "ARGS=[${params.flair_args}]" }\n')
    argv = ["nextflow", "-q", "run", str(nf)] + opt("--flair_args", "--stringent --check_splice")
    out = subprocess.run(argv, capture_output=True, text=True, cwd=tmp_path, timeout=300).stdout
    assert "ARGS=[--stringent --check_splice]" in out


def test_local_pipeline_has_no_revision(sample, two_runs):
    c = build_command([sample], two_runs, Options(pipeline="/home/me/flair-eval", ref="9.9"))
    assert c.argv[:4] == ["nextflow", "run", "/home/me/flair-eval", "--reads"]
    assert format_command(c.argv).split(" \\\n")[0] == "nextflow run /home/me/flair-eval"
    c = build_command([sample], two_runs, Options(pipeline="me/fork", ref="dev"))
    assert c.argv[:5] == ["nextflow", "run", "me/fork", "-r", "dev"]
