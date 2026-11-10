"""The UCSC PRISM "phoenix" SLURM cluster, which flair-eval assumes by default.

On phoenix the pipeline needs no -profile (its `standard` profile is conf/phoenix.config:
SLURM + conda + the shared lab FLAIR env cache), `--genome hg38` names a lab reference, and
the Nextflow driver runs as its own SLURM job instead of on the login node.

Detection, the genome presets and the driver job script live here so the rest of the
package has one place to ask. Keep GENOME_PRESETS in step with `params.genomes` in
nextflow.config and GENOME_PRESETS in frontend/builder/logic.js (tests check both).
"""

from __future__ import annotations

import datetime
import os
import re
import shutil
import subprocess
from pathlib import Path

#: Lab storage on phoenix; its presence (with sbatch) is how phoenix is recognised.
LAB_ROOT = "/private/groups/brookslab"
#: Shared FLAIR env cache set by conf/phoenix.config.
ENV_CACHE = f"{LAB_ROOT}/flair-eval/env_cache"
#: Where a comparison's --outdir should go (home has a ~100 GB quota).
OUTDIR_EXAMPLE = f"{LAB_ROOT}/$USER/flair-eval/dd_mm_yy_name"

#: `--genome <key>` presets: lab references on phoenix (same as params.genomes).
GENOME_PRESETS = {
    "hg38": {
        "fasta": f"{LAB_ROOT}/reference_sequence/GRCh38.primary_assembly.genome.fa",
        "gtf": f"{LAB_ROOT}/reference_annotations/gencode.v49.basic.annotation.gtf",
        "note": "GRCh38 primary assembly + GENCODE v49 basic",
    },
}

#: SLURM request for the Nextflow driver. It only schedules tasks, so it is small but long:
#: every pipeline step is its own SLURM job (conf/slurm.config).
RUN_JOB = {"partition": "long", "cpus": 2, "mem": "8G", "time": "7-00:00:00"}
#: Report-only launches build one HTML file; the driver finishes within hours.
REPORT_JOB = {"partition": "medium", "cpus": 2, "mem": "8G", "time": "12:00:00"}
JOB_NAME = "flair-eval"

_URI = re.compile(r"^[a-z][a-z0-9+.-]*://")


def detect(which=shutil.which, is_dir=os.path.isdir) -> bool:
    """True on phoenix: SLURM's sbatch is on PATH and the lab storage is mounted."""
    return bool(which("sbatch")) and bool(is_dir(LAB_ROOT))


def is_genome_preset(value) -> bool:
    return isinstance(value, str) and value in GENOME_PRESETS


def preset_names() -> str:
    return ", ".join(GENOME_PRESETS)


# ---------------------------------------------------------------- paths
def in_home(path: str, home: str | None = None) -> bool:
    """True when `path` (relative to the current folder) is inside the home folder."""
    home = home if home is not None else os.path.expanduser("~")
    if not home or home == "/":
        return False
    p = os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
    h = os.path.abspath(home)
    return p == h or p.startswith(h.rstrip("/") + "/")


#: Folders that exist only on the machine you are logged in to, not on the compute nodes.
NODE_LOCAL = ("/tmp", "/data/tmp", "/scratch", "/var/tmp")


def node_local(path: str) -> bool:
    """True when `path` is on login-node-local storage that SLURM jobs cannot see."""
    p = os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
    return any(p == d or p.startswith(d + "/") for d in NODE_LOCAL)


def absolute(path: str | None, base: str | None = None) -> str | None:
    """`path` as an absolute local path (from `base`, default the current folder).
    Empty values, remote URIs and genome preset keys are returned unchanged."""
    if not path or _URI.match(path) or is_genome_preset(path):
        return path
    p = os.path.expanduser(os.path.expandvars(path))
    return os.path.normpath(os.path.join(base or os.getcwd(), p))


def absolute_list(value: str | None, base: str | None = None) -> str | None:
    """Comma-separated paths or globs (e.g. --bundles), each made absolute."""
    if not value:
        return value
    return ",".join(absolute(v.strip(), base) for v in value.split(",") if v.strip())


# ---------------------------------------------------------------- driver job
def info_dir(outdir: str) -> Path:
    return Path(outdir) / "pipeline_info"


def script_path(outdir: str) -> Path:
    return info_dir(outdir) / "launch.sbatch"


def log_pattern(outdir: str) -> str:
    return str(info_dir(outdir) / "nextflow-%j.log")


def sbatch_script(command_text: str, outdir: str, job: dict | None = None,
                  version: str = "", now: datetime.datetime | None = None) -> str:
    """Job script that runs the Nextflow driver from inside `outdir` (absolute).

    Launching there keeps `.nextflow/` (the -resume cache) and `.nextflow.log` with the
    results instead of in the home folder."""
    from .command import quote
    job = job or RUN_JOB
    stamp = (now or datetime.datetime.now()).strftime("%Y-%m-%d %H:%M")
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={JOB_NAME}",
        f"#SBATCH --partition={job['partition']}",
        f"#SBATCH --cpus-per-task={job['cpus']}",
        f"#SBATCH --mem={job['mem']}",
        f"#SBATCH --time={job['time']}",
        f"#SBATCH --output={log_pattern(outdir)}",
        f"# Written by flair-eval{' ' + version if version else ''} on {stamp}.",
        "# Runs the Nextflow driver on a compute node; every pipeline step is its own SLURM job.",
        "# The driver launches inside --outdir, so .nextflow/ (for -resume) lives there too.",
        f"# Submit again with: sbatch {quote(str(script_path(outdir)))}",
        "set -euo pipefail",
        "# Group-writable files: lab members share the env and conda caches (and their locks).",
        "umask 002",
        f"cd {quote(outdir)}",
        command_text,
        "",
    ]
    return "\n".join(lines)


def submit(script: Path, run=subprocess.run) -> str:
    """sbatch the script; return the job id. Raises RuntimeError with sbatch's message."""
    try:
        p = run(["sbatch", "--parsable", str(script)], capture_output=True, text=True,
                timeout=120)
    except (OSError, subprocess.SubprocessError) as e:
        raise RuntimeError(f"sbatch failed: {e}") from e
    if p.returncode != 0:
        raise RuntimeError(f"sbatch failed: {(p.stderr or p.stdout or '').strip()}")
    job = (p.stdout or "").strip().split(";")[0]
    if not job:
        raise RuntimeError("sbatch printed no job id")
    return job
