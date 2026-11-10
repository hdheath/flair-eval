"""Tests for bin/build_report.py and the report assets."""

import copy
import json
import os
import re
import shutil
import subprocess
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
BIN = os.path.join(ROOT, "bin")
FIXTURES = os.path.join(ROOT, "tests", "fixtures", "bundles")
sys.path.insert(0, BIN)

from flair_eval_lib.schema import load_schema, validate  # noqa: E402

SCRIPT = os.path.join(BIN, "build_report.py")


def build(*args):
    return subprocess.run([sys.executable, SCRIPT, *map(str, args)], capture_output=True, text=True)


def embedded(path):
    html = open(path, encoding="utf-8").read()
    m = re.search(r'<script type="application/json" id="report-data">(.*?)</script>', html, re.S)
    assert m, "no embedded data"
    return html, json.loads(m.group(1))


def fixture_bundles():
    out = []
    for root, _, files in os.walk(FIXTURES):
        if "bundle.json" in files:
            out.append(os.path.join(root, "bundle.json"))
    return sorted(out)


def write_bundle(d, b):
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "bundle.json"), "w") as fh:
        json.dump(b, fh)


def test_fixtures_are_valid_bundles():
    files = fixture_bundles()
    assert len(files) >= 6
    runs, samples = set(), set()
    for f in files:
        b = json.load(open(f))
        assert validate(b, load_schema()) == [], f
        runs.add(b["run_id"])
        samples.add(b["sample"])
        tdir = os.path.join(os.path.dirname(f), "tables")
        iso = open(os.path.join(tdir, "isoforms.tsv")).read().splitlines()
        assert len(iso) - 1 == b["summary"]["isoforms"]
        cats = [r.split("\t") for r in open(os.path.join(tdir, "categories.tsv")).read().splitlines()[1:]]
        assert sum(int(r[2]) for r in cats) == b["summary"]["isoforms"]
        assert sum(b["categories"].values()) == b["summary"]["isoforms"]
        assert sum(b["fsm_subcategories"].values()) == b["categories"]["FSM"]
        assert b["inputs"]["gtf"] and "/" not in b["inputs"]["gtf"]   # no site paths committed
    assert len(runs) >= 3 and len(samples) >= 1


def test_build_from_fixtures(tmp_path):
    out = tmp_path / "r" / "report.html"
    res = build("--bundles", FIXTURES, "--out", out, "--title", "Demo <b>&</b>")
    assert res.returncode == 0, res.stderr
    html, data = embedded(out)
    n = len(fixture_bundles())
    assert len(data["bundles"]) == n
    assert data["baseline"] == data["runs"][0]
    assert set(data["samples"]) == {json.load(open(f))["sample"] for f in fixture_bundles()}
    # self-contained: no external scripts, styles, fonts or images
    assert not re.search(r'<script[^>]+src=', html)
    assert "<link" not in html and "@import" not in html
    assert not re.search(r'url\(\s*["\']?https?:', html)
    assert not re.search(r'fetch\(|XMLHttpRequest', html)
    # title escaped
    assert "Demo &lt;b&gt;&amp;&lt;/b&gt;" in html
    # no per-isoform tables embedded
    first_iso = open(os.path.join(os.path.dirname(fixture_bundles()[0]), "tables", "isoforms.tsv")).read().splitlines()[1].split("\t")[0]
    assert first_iso not in html


def test_baseline_option(tmp_path):
    data_runs = embedded_runs(tmp_path)
    out = tmp_path / "b.html"
    res = build("--bundles", FIXTURES, "--out", out, "--baseline", data_runs[-1])
    assert res.returncode == 0, res.stderr
    assert embedded(out)[1]["baseline"] == data_runs[-1]
    res = build("--bundles", FIXTURES, "--out", out, "--baseline", "no-such-run")
    assert res.returncode != 0 and "baseline" in res.stderr


def embedded_runs(tmp_path):
    out = tmp_path / "x.html"
    assert build("--bundles", FIXTURES, "--out", out).returncode == 0
    return embedded(out)[1]["runs"]


def test_newest_bundle_wins_and_glob(tmp_path):
    b = json.load(open(fixture_bundles()[0]))
    old = copy.deepcopy(b)
    old["pipeline"]["created"] = "2020-01-01T00:00:00+00:00"
    old["summary"]["isoforms"] = 1
    new = copy.deepcopy(b)
    new["pipeline"]["created"] = "2030-01-01T00:00:00+00:00"
    new["summary"]["isoforms"] = 99999
    write_bundle(tmp_path / "launch1" / "bundles" / b["sample"] / b["run_id"], old)
    write_bundle(tmp_path / "launch2" / "bundles" / b["sample"] / b["run_id"], new)
    out = tmp_path / "r.html"
    res = build("--bundles", str(tmp_path / "launch*" / "bundles"), "--out", out)
    assert res.returncode == 0, res.stderr
    data = embedded(out)[1]
    assert len(data["bundles"]) == 1
    assert data["bundles"][0]["summary"]["isoforms"] == 99999


def test_unknown_sections_and_missing_optional_parts(tmp_path):
    b = json.load(open(fixture_bundles()[0]))
    b["sections"] = {"some_future_section": {"x": [1, 2, 3]}}
    b.pop("compute", None)
    b.pop("inputs", None)
    b["label"] = "</script><script>alert(1)</script>"
    write_bundle(tmp_path / "one", b)
    out = tmp_path / "r.html"
    res = build("--bundles", tmp_path / "one", "--out", out)
    assert res.returncode == 0, res.stderr
    html, data = embedded(out)
    assert data["bundles"][0]["label"] == b["label"]
    assert html.count("</script>") == 2          # the data block and the code block only
    assert data["bundles"][0]["sections"]["some_future_section"] == {"x": [1, 2, 3]}


def test_rejects_bad_input(tmp_path):
    out = tmp_path / "r.html"
    res = build("--bundles", tmp_path / "missing", "--out", out)
    assert res.returncode != 0 and "error" in res.stderr
    (tmp_path / "empty").mkdir()
    res = build("--bundles", tmp_path / "empty", "--out", out)
    assert res.returncode != 0
    b = json.load(open(fixture_bundles()[0]))
    del b["summary"]
    write_bundle(tmp_path / "bad", b)
    res = build("--bundles", tmp_path / "bad", "--out", out)
    assert res.returncode != 0 and "summary" in res.stderr
    assert not out.exists()


def test_size_under_5mb_for_50_bundles(tmp_path):
    src = [json.load(open(f)) for f in fixture_bundles()]
    for i in range(50):
        b = copy.deepcopy(src[i % len(src)])
        b["sample"] = f"sample{i // 10}"
        b["run_id"] = f"run-{i % 10}"
        write_bundle(tmp_path / "many" / b["sample"] / b["run_id"], b)
    out = tmp_path / "r.html"
    res = build("--bundles", tmp_path / "many", "--out", out)
    assert res.returncode == 0, res.stderr
    assert len(embedded(out)[1]["bundles"]) == 50
    assert os.path.getsize(out) < 5 * 1024 * 1024


def test_help():
    res = build("--help")
    assert res.returncode == 0 and "--bundles" in res.stdout


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_report_js_parses():
    res = subprocess.run(["node", "--check", os.path.join(ROOT, "assets", "report", "report.js")],
                         capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
