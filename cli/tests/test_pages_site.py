"""The GitHub Pages site (.github/workflows/pages.yml) must still assemble."""
import re
import subprocess
import sys

import pytest


@pytest.fixture
def workflow(repo):
    p = repo / ".github" / "workflows" / "pages.yml"
    if not p.exists():
        pytest.skip("not in the pipeline repository")
    return p.read_text()


def test_workflow_copies_every_file_the_builder_loads(repo, workflow):
    html = (repo / "frontend" / "builder" / "index.html").read_text()
    for src in re.findall(r'<script src="([^"]+)"', html):
        assert f"frontend/builder/{src}" in workflow, src
    assert "frontend/builder/welcome.html _site/index.html" in workflow


def test_welcome_links_point_into_the_site_layout(repo, workflow):
    html = (repo / "frontend" / "builder" / "welcome.html").read_text()
    local = {h for h in re.findall(r'href="([^"#]+)"', html) if not h.startswith("http")}
    assert local == {"builder/", "demo/report.html"}
    assert "_site/builder/" in workflow and "--out _site/demo/report.html" in workflow


def test_demo_report_builds_from_fixtures(repo, workflow, tmp_path):
    script, bundles = repo / "bin" / "build_report.py", repo / "tests" / "fixtures" / "bundles"
    if not script.exists() or not bundles.is_dir():
        pytest.skip("report code or fixtures not here")
    out = tmp_path / "report.html"
    p = subprocess.run([sys.executable, str(script), "--bundles", str(bundles), "--out", str(out),
                        "--baseline", "flair-3.0.1", "--title", "Demo"],
                       capture_output=True, text=True, timeout=300)
    assert p.returncode == 0, p.stderr
    assert out.stat().st_size > 10_000
