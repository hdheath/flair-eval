import subprocess
from collections import namedtuple

from flair_eval import doctor

Usage = namedtuple("Usage", "total used free")


def runner(outputs):
    def run(cmd, **k):
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr=outputs.get(cmd[0], ""))
    return run


def statuses(checks):
    return {c.name: c.status for c in checks}


def test_suggest_profile():
    have = lambda *names: (lambda x: "/bin/" + x if x in names else None)  # noqa: E731
    assert doctor.suggest_profile(have("sbatch", "conda"), on_phoenix=False) == "slurm,conda"
    assert doctor.suggest_profile(have("apptainer"), on_phoenix=False) == "local,singularity"
    assert doctor.suggest_profile(have("docker"), on_phoenix=False) == "local,docker"
    assert doctor.suggest_profile(have(), on_phoenix=False) == "local"
    assert doctor.suggest_profile(have(), on_phoenix=None) == "local"  # no sbatch: not phoenix
    # phoenix: the pipeline's default profile is phoenix, so no -profile at all
    assert doctor.suggest_profile(have("sbatch", "conda"), on_phoenix=True) == ""


def test_phoenix_detection():
    from flair_eval import phoenix
    sb = lambda x: "/usr/bin/sbatch" if x == "sbatch" else None  # noqa: E731
    lab = lambda p: p == phoenix.LAB_ROOT  # noqa: E731
    assert phoenix.detect(which=sb, is_dir=lab) is True
    assert phoenix.detect(which=lambda x: None, is_dir=lab) is False   # no SLURM
    assert phoenix.detect(which=sb, is_dir=lambda p: False) is False   # another cluster


def test_phoenix_doctor_report(tmp_path):
    from flair_eval import phoenix
    which = lambda x: "/bin/" + x  # noqa: E731
    out = {"java": 'openjdk version "17.0.9"', "nextflow": "nextflow version 25.10.0"}
    checks = doctor.run_checks(str(tmp_path), which=which, run=runner(out),
                               disk_usage=lambda p: Usage(0, 0, 500e9), home="/nonexistent",
                               on_phoenix=True, is_dir=lambda p: True)
    cluster = [c for c in checks if c.name == "Cluster"][0]
    assert "phoenix" in cluster.detail and cluster.status == "ok"
    envs = [c for c in checks if c.name == "FLAIR envs"][0]
    assert phoenix.ENV_CACHE in envs.detail
    lines = []
    assert doctor.report(checks, "", out=lines.append, on_phoenix=True) == 0
    text = "\n".join(lines)
    assert "phoenix (UCSC PRISM) detected: no -profile needed" in text
    assert phoenix.ENV_CACHE in text and "Suggested profile" not in text
    # home outdir on phoenix points at groups storage
    checks = doctor.run_checks(str(tmp_path / "out"), which=which, run=runner(out),
                               disk_usage=lambda p: Usage(0, 0, 500e9), home=str(tmp_path),
                               on_phoenix=True, is_dir=lambda p: False)
    st = statuses(checks)
    assert st["Output folder"] == "warn" and st["FLAIR envs"] == "warn"
    assert phoenix.OUTDIR_EXAMPLE in [c for c in checks if c.name == "Output folder"][0].fix


def test_nextflow_version():
    assert doctor.nextflow_version("nextflow version 25.10.0.10289") == (25, 10)
    assert doctor.nextflow_version("junk") is None


def test_all_good(tmp_path):
    which = lambda x: "/bin/" + x  # noqa: E731
    out = {"java": 'openjdk version "17.0.9"', "nextflow": "nextflow version 25.10.0"}
    checks = doctor.run_checks(str(tmp_path), which=which, run=runner(out),
                               disk_usage=lambda p: Usage(0, 0, 500e9), home="/nonexistent",
                               on_phoenix=False)
    assert set(statuses(checks).values()) == {"ok"}
    lines = []
    assert doctor.report(checks, "slurm,conda", out=lines.append) == 0
    assert any("-profile slurm,conda" in x for x in lines)


def test_problems_and_fixes(tmp_path):
    which = lambda x: "/bin/" + x if x in ("java", "nextflow") else None  # noqa: E731
    out = {"java": 'java version "1.8.0_292"', "nextflow": "nextflow version 23.10.1"}
    home = str(tmp_path)
    checks = doctor.run_checks(str(tmp_path / "out"), which=which, run=runner(out),
                               disk_usage=lambda p: Usage(0, 0, 10e9), home=home)
    st = statuses(checks)
    assert st["Java"] == "fail" and st["Nextflow"] == "fail"
    assert st["conda"] == "warn" and st["Containers"] == "warn" and st["git"] == "warn"
    assert st["Disk"] == "warn" and st["Output folder"] == "warn"
    assert "self-update" in [c for c in checks if c.name == "Nextflow"][0].fix
    assert doctor.report(checks, "local", out=lambda m: None) == 1


def test_missing_tools():
    checks = doctor.run_checks(".", which=lambda x: None, run=runner({}),
                               disk_usage=lambda p: Usage(0, 0, 100e9), home="/nonexistent")
    st = statuses(checks)
    assert st["Java"] == "fail" and st["Nextflow"] == "fail"
