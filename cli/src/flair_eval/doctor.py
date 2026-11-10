"""`flair-eval doctor`: check this machine before a first launch, with a fix for each problem.

Also picks a Nextflow profile that suits the machine (`suggest_profile`), which
`flair-eval run` uses when no --profile is given.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from .env_check import JAVA_HINT, MIN_JAVA, NEXTFLOW_HINT, java_major

MIN_NEXTFLOW = (25, 4)
LOW_DISK_GB = 50


@dataclass
class Check:
    status: str  # "ok", "warn" or "fail"
    name: str
    detail: str
    fix: str = ""


def nextflow_version(text: str) -> tuple | None:
    m = re.search(r"version\s+(\d+)\.(\d+)", text)
    return (int(m.group(1)), int(m.group(2))) if m else None


def suggest_profile(which=shutil.which) -> str:
    """Executor + software profile for this machine, e.g. 'slurm,conda' or 'local,singularity'."""
    executor = "slurm" if which("sbatch") else "local"
    if which("conda"):
        software = "conda"
    elif which("singularity") or which("apptainer"):
        software = "singularity"
    elif which("docker"):
        software = "docker"
    else:
        software = ""
    return ",".join(p for p in (executor, software) if p)


def _version_of(cmd: list, run) -> str:
    try:
        p = run(cmd, capture_output=True, text=True, timeout=120)
        return (p.stdout or "") + (p.stderr or "")
    except (OSError, subprocess.SubprocessError):
        return ""


def run_checks(outdir: str | None = None, which=shutil.which, run=subprocess.run,
               disk_usage=shutil.disk_usage, home: str | None = None) -> list[Check]:
    checks: list[Check] = []
    v = sys.version_info
    checks.append(Check("ok", "Python", f"{v.major}.{v.minor}.{v.micro}"))

    if which("java") is None:
        checks.append(Check("fail", "Java", "not found", JAVA_HINT))
    else:
        major = java_major(_version_of(["java", "-version"], run))
        if major is None:
            checks.append(Check("warn", "Java", "found, but `java -version` gave no version"))
        elif major < MIN_JAVA:
            checks.append(Check("fail", "Java", f"version {major}; Nextflow needs {MIN_JAVA}+",
                                JAVA_HINT))
        else:
            checks.append(Check("ok", "Java", f"version {major}"))

    if which("nextflow") is None:
        checks.append(Check("fail", "Nextflow", "not found", NEXTFLOW_HINT))
    else:
        nv = nextflow_version(_version_of(["nextflow", "-v"], run))
        want = ".".join(map(str, MIN_NEXTFLOW))
        if nv is None:
            checks.append(Check("warn", "Nextflow", "found, but `nextflow -v` gave no version"))
        elif nv < MIN_NEXTFLOW:
            checks.append(Check("fail", "Nextflow", f"version {nv[0]}.{nv[1]}; flair-eval needs "
                                f"{want}+", "Update it: nextflow self-update"))
        else:
            checks.append(Check("ok", "Nextflow", f"version {nv[0]}.{nv[1]}"))

    conda = which("conda")
    sing = which("singularity") or which("apptainer")
    docker = which("docker")
    if conda:
        checks.append(Check("ok", "conda", "found: FLAIR releases and git branches can be built"))
    else:
        checks.append(Check("warn", "conda", "not found: release and git runs need it to build FLAIR "
                            "(env: and container: runs do not)",
                            "Install Miniforge: https://conda-forge.org/download/ "
                            "(or `module load conda` on a cluster)"))
    if sing or docker:
        names = [n for n, ok in (("singularity/apptainer", sing), ("docker", docker)) if ok]
        checks.append(Check("ok", "Containers", " and ".join(names) + " found"))
    else:
        checks.append(Check("warn", "Containers", "no singularity, apptainer or docker: "
                            "container: runs will not work", "Only needed for container: runs."))
    if not which("git"):
        checks.append(Check("warn", "git", "not found: git runs (dev, main, git:...) need it",
                            "conda install -c conda-forge git"))
    if which("sbatch"):
        checks.append(Check("ok", "Cluster", "SLURM found (sbatch): jobs can run on the cluster"))
    else:
        checks.append(Check("ok", "Cluster", "no SLURM: jobs run on this machine"))

    target = Path(outdir or ".").expanduser().resolve()
    probe = target
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    try:
        free_gb = disk_usage(str(probe)).free / 1e9
        msg = f"{free_gb:,.0f} GB free at {probe}"
        if free_gb < LOW_DISK_GB:
            checks.append(Check("warn", "Disk", msg, f"Use --outdir and --work-dir on storage with "
                                f"more than {LOW_DISK_GB} GB free."))
        else:
            checks.append(Check("ok", "Disk", msg))
    except OSError:
        checks.append(Check("warn", "Disk", f"could not read free space at {probe}"))
    home = home if home is not None else os.path.expanduser("~")
    if home and (str(target) + "/").startswith(home.rstrip("/") + "/"):
        checks.append(Check("warn", "Output folder", f"{target} is in your home folder, which often "
                            "has a small quota", "Pass --outdir and --work-dir on project or "
                            "scratch storage."))
    return checks


def report(checks: list[Check], profile: str, out=print) -> int:
    mark = {"ok": "ok  ", "warn": "WARN", "fail": "FAIL"}
    out("flair-eval doctor\n")
    for c in checks:
        out(f"  {mark[c.status]}  {c.name}: {c.detail}")
        for line in c.fix.splitlines():
            out(f"        {line}")
    fails = [c for c in checks if c.status == "fail"]
    out("")
    out(f"Suggested profile for this machine: -profile {profile}" if profile else
        "No container or conda found; install conda to build FLAIR envs.")
    if fails:
        out(f"\n{len(fails)} problem(s) to fix before launching.")
        return 1
    out("\nReady. Next: flair-eval init (or try the test data: flair-eval run --test)")
    return 0
