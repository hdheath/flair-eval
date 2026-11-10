"""Check that Java and Nextflow are installed before launching, with install hints."""

from __future__ import annotations

import re
import shutil
import subprocess

MIN_JAVA = 17

JAVA_HINT = f"""Java {MIN_JAVA} or newer is needed to run Nextflow. Install one of:
    conda install -c conda-forge openjdk=17
    curl -s https://get.sdkman.io | bash && sdk install java 17.0.13-tem
or load it on your cluster, e.g. `module load java/17`."""

NEXTFLOW_HINT = """Nextflow is not on PATH. Install one of:
    conda create -y -n nextflow -c conda-forge -c bioconda nextflow openjdk=17
    conda activate nextflow
    curl -s https://get.nextflow.io | bash && mv nextflow ~/bin/   (any directory on PATH)
Docs: https://www.nextflow.io/docs/latest/install.html"""


def java_major(version_text: str) -> int | None:
    """Parse `java -version` output: '1.8.0_292' -> 8, '17.0.2' -> 17, '21' -> 21."""
    m = re.search(r'version "?(\d+)(?:\.(\d+))?', version_text)
    if not m:
        return None
    major = int(m.group(1))
    if major == 1 and m.group(2):
        major = int(m.group(2))
    return major


def check_tools(which=shutil.which, run=subprocess.run) -> tuple[list[str], list[str]]:
    """Return (errors, warnings). `which` and `run` are injectable for tests."""
    errors, warnings = [], []
    if which("java") is None:
        errors.append("java not found.\n" + JAVA_HINT)
    else:
        try:
            p = run(["java", "-version"], capture_output=True, text=True, timeout=30)
            major = java_major((p.stderr or "") + (p.stdout or ""))
        except (OSError, subprocess.SubprocessError):
            major = None
        if major is None:
            warnings.append("could not read the Java version (`java -version`).")
        elif major < MIN_JAVA:
            warnings.append(f"Java {major} found; current Nextflow needs Java {MIN_JAVA}+.\n"
                            + JAVA_HINT)
    if which("nextflow") is None:
        errors.append(NEXTFLOW_HINT)
    return errors, warnings
