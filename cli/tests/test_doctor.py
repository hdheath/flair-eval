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
    assert doctor.suggest_profile(have("sbatch", "conda")) == "slurm,conda"
    assert doctor.suggest_profile(have("apptainer")) == "local,singularity"
    assert doctor.suggest_profile(have("docker")) == "local,docker"
    assert doctor.suggest_profile(have()) == "local"


def test_nextflow_version():
    assert doctor.nextflow_version("nextflow version 25.10.0.10289") == (25, 10)
    assert doctor.nextflow_version("junk") is None


def test_all_good(tmp_path):
    which = lambda x: "/bin/" + x  # noqa: E731
    out = {"java": 'openjdk version "17.0.9"', "nextflow": "nextflow version 25.10.0"}
    checks = doctor.run_checks(str(tmp_path), which=which, run=runner(out),
                               disk_usage=lambda p: Usage(0, 0, 500e9), home="/nonexistent")
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
