"""Unit tests for bin/flair_env.py with fake conda, pip, git, PyPI and apptainer."""
import json
import os
import stat
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FLAIR_ENV = os.path.join(ROOT, "bin", "flair_env.py")
sys.path.insert(0, os.path.join(ROOT, "bin"))
import flair_env  # noqa: E402

FIELDS = ["run_id", "source", "requested", "resolved_commit", "flair_version", "conda_prefix",
          "container", "python", "pip_freeze_sha256", "created"]
SHA = "d033600f4a60c7929d1aaa6b5b6c9d4b1e4a2cc4"

# A fake env python: handles `-m pip install X`, `-m pip freeze` and the platform one-liner.
FAKE_PY = r'''
import os, sys
prefix = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
a = sys.argv[1:]
if a[:3] == ["-m", "pip", "install"]:
    target = a[-1]
    if "==" in target:
        ver = target.split("==")[1]
    else:
        ver = "3.0.0+dev" if os.path.exists(os.path.join(target, "pyproject.toml")) else "?"
    flair = os.path.join(prefix, "bin", "flair")
    with open(flair, "w") as fh:
        fh.write("#!/bin/sh\necho 'FLAIR %s'\n" % ver)
    os.chmod(flair, 0o755)
    open(os.path.join(prefix, "installed.txt"), "w").write(target)
elif a[:3] == ["-m", "pip", "freeze"]:
    print("flair-brookslab==%s" % open(os.path.join(prefix, "installed.txt")).read())
    print("pysam==0.23.3")
elif a and a[0] == "-c":
    print("3.12.7")
'''

FAKE_CONDA = r'''
import os, sys, time
a = sys.argv[1:]
log = os.environ["FAKE_CONDA_LOG"]
with open(log, "a") as fh:
    fh.write(" ".join(a) + "\n")
if os.environ.get("FAKE_CONDA_FAIL"):
    sys.exit(1)
time.sleep(float(os.environ.get("FAKE_CONDA_SLEEP", "0")))
prefix = a[a.index("-p") + 1]
os.makedirs(os.path.join(prefix, "bin"))
py = os.path.join(prefix, "bin", "python")
with open(py, "w") as fh:
    fh.write("#!%s\n" % sys.executable + open(os.environ["FAKE_PY_SRC"]).read())
os.chmod(py, 0o755)
for spec in a:
    if spec.startswith("flair="):          # bioconda route: flair comes with the env
        ver = spec.split("=")[1]
        open(os.path.join(prefix, "installed.txt"), "w").write("bioconda " + ver)
        f = os.path.join(prefix, "bin", "flair")
        open(f, "w").write("#!/bin/sh\necho 'FLAIR %s'\n" % ver)
        os.chmod(f, 0o755)
'''

FAKE_GIT = r'''
import os, sys
SHA = "%s"
a = sys.argv[1:]
open(os.environ["FAKE_GIT_LOG"], "a").write(" ".join(a) + "\n")
src = None
if a[0] == "-C":
    src, a = a[1], a[2:]
if a[0] == "ls-remote":
    print("%%s\trefs/heads/dev" %% SHA)
    print("1111111111111111111111111111111111111111\trefs/tags/v9")
    print("2222222222222222222222222222222222222222\trefs/tags/v9^{}")
elif a[0] == "init":
    os.makedirs(a[-1], exist_ok=True)
elif a[0] == "checkout":
    open(os.path.join(src, "pyproject.toml"), "w").write('requires-python = ">=3.12,<3.13"\n')
    os.makedirs(os.path.join(src, "misc"))
    open(os.path.join(src, "misc", "flair_conda_env.yaml"), "w").write(
        "dependencies:\n  - python >=3.12\n  - minimap2 =2.30\n  - ucsc-bedtobigbed\n  - r-ggplot2\n"
        "  - pip\n  - pip:\n      - flair-brookslab\n")
    open(os.path.join(src, ".head"), "w").write(SHA)
elif a[0] == "rev-parse":
    p = os.path.join(src, ".head")
    if not os.path.exists(p):
        sys.exit(128)
    print(open(p).read())
''' % SHA

FAKE_APPTAINER = r'''
import sys
a = sys.argv[1:]
assert a[0] == "exec"
cmd = a[2:]
if cmd == ["flair", "--version"]:
    print("FLAIR 3.0.1")
elif cmd[:2] == ["python", "-c"]:
    print("3.13.1")
elif cmd[:4] == ["python", "-m", "pip", "freeze"]:
    print("flair-brookslab==3.0.1")
'''


def write_exe(path, body, shebang=True):
    with open(path, "w") as fh:
        fh.write(("#!%s\n" % sys.executable if shebang else "") + body)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


@pytest.fixture
def world(tmp_path):
    b = tmp_path / "fakebin"
    b.mkdir()
    write_exe(str(b / "mamba"), FAKE_CONDA)
    write_exe(str(b / "git"), FAKE_GIT)
    write_exe(str(b / "apptainer"), FAKE_APPTAINER)
    pysrc = tmp_path / "fake_py.py"
    pysrc.write_text(FAKE_PY)
    pypi = tmp_path / "pypi" / "flair-brookslab"
    for v, req in (("3.0.1", ">=3.12,<3.15"), ("2.2.0", ">=3.12,<3.13")):
        (pypi / v).mkdir(parents=True)
        (pypi / v / "json").write_text(json.dumps({"info": {"requires_python": req, "version": v}}))
    rel = tmp_path / "relyaml"
    rel.mkdir()
    (rel / "v3.0.1.yaml").write_text("dependencies:\n  - python >=3.12\n  - samtools =1.22.1\n"
                                     "  - intron-prospector >=1.5  # junctions\n  - pip\n")
    env = dict(os.environ)
    env.update({"FLAIR_EVAL_RELEASE_YAML": "file://" + str(rel) + "/%s.yaml"})
    env.update({"PATH": str(b) + os.pathsep + env["PATH"], "FAKE_PY_SRC": str(pysrc),
                "FAKE_CONDA_LOG": str(tmp_path / "conda.log"), "FAKE_GIT_LOG": str(tmp_path / "git.log"),
                "FLAIR_EVAL_PYPI": "file://" + str(tmp_path / "pypi"), "FLAIR_EVAL_CONDA": str(b / "mamba"),
                "CONDA_PKGS_DIRS": str(tmp_path / "pkgs")})
    for k in ("FLAIR_EVAL_TOOL_PINS", "FLAIR_EVAL_CONTAINER_RUNTIME"):
        env.pop(k, None)
    cache = tmp_path / "cache"

    class W:
        pass
    w = W()
    w.tmp, w.env, w.cache, w.bin = tmp_path, env, cache, b

    def resolve(run, out_name="flair_env.json", wait=True, extra=()):
        rj = tmp_path / ("run_%s.json" % run["id"])
        rj.write_text(json.dumps(run))
        out = tmp_path / out_name
        argv = [sys.executable, FLAIR_ENV, "resolve", "--run-json", str(rj), "--cache", str(cache),
                "--out", str(out)] + list(extra)
        if not wait:
            return subprocess.Popen(argv, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True), out
        p = subprocess.run(argv, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           universal_newlines=True)
        rec = json.load(open(out)) if p.returncode == 0 else None
        return p, rec

    w.resolve = resolve
    w.conda_calls = lambda: [l for l in open(env["FAKE_CONDA_LOG"]).read().splitlines()] \
        if os.path.exists(env["FAKE_CONDA_LOG"]) else []
    return w


def check_record(rec):
    assert list(rec) == FIELDS
    assert rec["created"].endswith("+00:00")


def test_help():
    p = subprocess.run([sys.executable, FLAIR_ENV, "resolve", "--help"], stdout=subprocess.PIPE,
                       universal_newlines=True)
    assert p.returncode == 0 and "--run-json" in p.stdout
    assert subprocess.run([sys.executable, FLAIR_ENV], stdout=subprocess.PIPE).returncode == 2


def test_release_pip_route_and_cache(world):
    p, rec = world.resolve({"id": "flair-3.0.1", "flair": {"release": "3.0.1"}})
    assert p.returncode == 0, p.stderr
    check_record(rec)
    assert rec["source"] == "release" and rec["requested"] == "3.0.1"
    assert rec["flair_version"] == "FLAIR 3.0.1" and rec["resolved_commit"] is None
    assert rec["conda_prefix"] == str(world.cache / "envs" / "release-3.0.1")
    assert rec["python"] == "3.12.7" and len(rec["pip_freeze_sha256"]) == 64
    calls = world.conda_calls()
    assert len(calls) == 1
    assert "python=3.12" in calls[0] and "minimap2=2.30" in calls[0] and "samtools=1.22.1" in calls[0]
    assert "-c conda-forge -c bioconda" in calls[0]
    assert "intron-prospector>=1.5" in calls[0] and "samtools=1.22.1" in calls[0]
    assert open(os.path.join(rec["conda_prefix"], "installed.txt")).read() == "flair-brookslab==3.0.1"
    marker = json.load(open(os.path.join(rec["conda_prefix"], ".flair-eval-env.json")))
    assert marker["status"] == "complete" and marker["route"] == "pip"
    # second resolve reuses the cache
    p, rec2 = world.resolve({"id": "again", "flair": {"release": "3.0.1"}})
    assert p.returncode == 0 and len(world.conda_calls()) == 1
    assert rec2["conda_prefix"] == rec["conda_prefix"]
    assert rec2["pip_freeze_sha256"] == rec["pip_freeze_sha256"]


def test_release_bioconda_fallback(world):
    p, rec = world.resolve({"id": "old", "flair": {"release": "2.1.1"}})   # not on (fake) PyPI
    assert p.returncode == 0, p.stderr
    assert rec["flair_version"] == "FLAIR 2.1.1"
    call = world.conda_calls()[0]
    assert "flair=2.1.1" in call and "python=" not in call
    marker = json.load(open(os.path.join(rec["conda_prefix"], ".flair-eval-env.json")))
    assert marker["route"] == "bioconda"


def test_tool_pins_override(world):
    world.env["FLAIR_EVAL_TOOL_PINS"] = "minimap2=2.28"
    p, rec = world.resolve({"id": "r", "flair": {"release": "2.2.0"}})
    assert p.returncode == 0, p.stderr
    call = world.conda_calls()[0]
    assert "minimap2=2.28" in call and "samtools" not in call


def test_concurrent_resolves_build_once(world):
    world.env["FAKE_CONDA_SLEEP"] = "1.5"
    procs = [world.resolve({"id": "r%d" % i, "flair": {"release": "3.0.1"}}, out_name="o%d.json" % i,
                           wait=False) for i in range(3)]
    recs = []
    for proc, out in procs:
        _, err = proc.communicate(timeout=120)
        assert proc.returncode == 0, err
        recs.append(json.load(open(out)))
    assert len(world.conda_calls()) == 1
    assert len({r["conda_prefix"] for r in recs}) == 1


def test_incomplete_env_is_rebuilt(world):
    prefix = world.cache / "envs" / "release-3.0.1"
    (prefix / "bin").mkdir(parents=True)          # a crashed build: no marker
    p, rec = world.resolve({"id": "r", "flair": {"release": "3.0.1"}})
    assert p.returncode == 0, p.stderr
    assert len(world.conda_calls()) == 1
    assert os.path.exists(prefix / ".flair-eval-env.json")


def test_build_failure_leaves_no_marker(world):
    world.env["FAKE_CONDA_FAIL"] = "1"
    p, rec = world.resolve({"id": "r", "flair": {"release": "3.0.1"}})
    assert p.returncode == 1
    err = p.stderr
    assert "could not set up FLAIR for run 'r' (release 3.0.1)" in err
    assert "What failed: `conda create (base env with Python and tools)` exited with code 1" in err
    assert "Command:     " in err and " create -y -p " in err
    assert "Log:         %s" % (world.cache / "logs" / "release-3.0.1.build.log") in err
    assert "What to try:" in err and "rebuilds it" in err
    assert not os.path.exists(world.cache / "envs" / "release-3.0.1" / ".flair-eval-env.json")
    assert os.path.exists(world.cache / "logs" / "release-3.0.1.build.log")


def test_git_branch_resolves_to_sha(world):
    p, rec = world.resolve({"id": "flair-dev", "flair": {"git": "https://example.org/flair", "ref": "dev"}})
    assert p.returncode == 0, p.stderr
    check_record(rec)
    assert rec["source"] == "git" and rec["resolved_commit"] == SHA
    assert rec["requested"] == "https://example.org/flair@dev"
    assert rec["conda_prefix"] == str(world.cache / "envs" / ("git-" + SHA[:12]))
    assert open(os.path.join(rec["conda_prefix"], "installed.txt")).read() == str(world.cache / "src" / SHA)
    call = world.conda_calls()[0]
    assert "python=3.12" in call
    assert "ucsc-bedtobigbed" in call and "r-ggplot2" not in call and call.count("minimap2") == 1
    gitlog = open(world.env["FAKE_GIT_LOG"]).read()
    assert "fetch -q --depth 1 origin " + SHA in gitlog
    # same commit through a SHA ref: cache hit
    p, rec2 = world.resolve({"id": "by-sha", "flair": {"git": "https://example.org/flair", "ref": SHA}})
    assert p.returncode == 0 and rec2["conda_prefix"] == rec["conda_prefix"]
    assert len(world.conda_calls()) == 1


def test_git_annotated_tag_is_dereferenced():
    out = "1111111111111111111111111111111111111111\trefs/tags/v9\n" \
          "2222222222222222222222222222222222222222\trefs/tags/v9^{}\n"
    orig = flair_env.capture
    try:
        flair_env.capture = lambda argv, env=None, timeout=900: (0, out, "")
        assert flair_env.resolve_git_ref("u", "v9") == "2" * 40
        assert flair_env.resolve_git_ref("u", "abc1234") is None     # short SHA needs a clone
        with pytest.raises(flair_env.EnvError):
            flair_env.resolve_git_ref("u", "nope")
    finally:
        flair_env.capture = orig


def test_env_passthrough(world, tmp_path):
    prefix = tmp_path / "myenv"
    (prefix / "bin").mkdir(parents=True)
    write_exe(str(prefix / "bin" / "flair"), "echo 'FLAIR 3.0.0+dev'\n", shebang=False)
    with open(prefix / "bin" / "flair") as fh:
        body = fh.read()
    (prefix / "bin" / "flair").write_text("#!/bin/sh\n" + body)
    os.symlink(sys.executable, str(prefix / "bin" / "python"))
    sp = prefix / "lib" / "python3.12" / "site-packages" / "flair_brookslab-3.0.0.dist-info"
    sp.mkdir(parents=True)
    (sp / "direct_url.json").write_text(json.dumps({"url": "https://github.com/x/flair",
                                                    "vcs_info": {"vcs": "git", "commit_id": SHA}}))
    p, rec = world.resolve({"id": "env-myenv", "flair": {"env": str(prefix)}})
    assert p.returncode == 0, p.stderr
    check_record(rec)
    assert rec["conda_prefix"] == str(prefix) and rec["source"] == "env"
    assert rec["flair_version"] == "FLAIR 3.0.0+dev" and rec["resolved_commit"] == SHA
    assert rec["python"] and world.conda_calls() == []


def test_env_missing_or_broken(world, tmp_path):
    p, _ = world.resolve({"id": "x", "flair": {"env": str(tmp_path / "nope")}})
    assert p.returncode == 1 and "does not exist" in p.stderr and "What to try" in p.stderr
    broken = tmp_path / "broken"
    (broken / "bin").mkdir(parents=True)
    (broken / "bin" / "flair").write_text("#!/bin/sh\necho boom >&2; exit 1\n")
    os.chmod(str(broken / "bin" / "flair"), 0o755)
    p, _ = world.resolve({"id": "x", "flair": {"env": str(broken)}})
    assert p.returncode == 1 and "| boom" in p.stderr and "FLAIR does not start" in p.stderr


def test_container_passthrough(world, tmp_path):
    img = tmp_path / "flair.sif"
    img.write_text("")
    p, rec = world.resolve({"id": "ctr-flair", "flair": {"container": str(img)}})
    assert p.returncode == 0, p.stderr
    check_record(rec)
    assert rec["container"] == str(img) and rec["conda_prefix"] is None
    assert rec["flair_version"] == "FLAIR 3.0.1" and rec["python"] == "3.13.1"
    p, rec = world.resolve({"id": "ctr2", "flair": {"container": "docker://example/flair:3.0.1"}})
    assert p.returncode == 0 and rec["container"] == "docker://example/flair:3.0.1"


@pytest.mark.parametrize("flair", [
    {}, {"release": "3.0.1", "env": "/x"}, {"git": "u"}, {"release": "3.0.1", "ref": "dev"},
    {"release": ""}, {"bogus": "1"}, "3.0.1",
])
def test_bad_run_maps(world, flair):
    p, _ = world.resolve({"id": "bad", "flair": flair})
    assert p.returncode != 0


def test_bad_run_json(world, tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text("{not json")
    p = subprocess.run([sys.executable, FLAIR_ENV, "resolve", "--run-json", str(bad), "--cache",
                        str(tmp_path / "c"), "--out", str(tmp_path / "o.json")], env=world.env,
                       stderr=subprocess.PIPE, universal_newlines=True)
    assert p.returncode == 2


@pytest.mark.parametrize("spec,expected", [
    (">=3.12,<3.13", "3.12"), (">=3.12,<3.15", "3.12"), (">=3.13", "3.13"), (None, "3.12"),
    ("~=3.10", "3.12"), (">=3.8,<3.11", "3.10"), ("==3.11.*", "3.11"),
])
def test_pick_python(spec, expected):
    assert flair_env.pick_python(spec) == expected


def test_pick_python_impossible():
    with pytest.raises(flair_env.EnvError):
        flair_env.pick_python(">=4")


def test_requires_python_of(tmp_path):
    (tmp_path / "pyproject.toml").write_text('[project]\nrequires-python = ">=3.12,<3.13"\n')
    assert flair_env.requires_python_of(str(tmp_path)) == ">=3.12,<3.13"


def test_shipped_conda_extras():
    y = """name: x
channels:
  - conda-forge
dependencies:
  - python >=3.12,<3.15
  - minimap2 =2.30
  - pysam >=0.24.1,<0.25.0   # flair's range
  - intron-prospector >=1.5
  - ucsc-bedtobigbed
  - longshot >=1.0
  - r-argparse
  - bioconductor-deseq2
  - pip
  - pip:
      - flair-brookslab==v3.0.0+dev
"""
    pins = flair_env.tool_pins("minimap2=2.30 samtools=1.22.1")
    assert flair_env.shipped_conda_extras(y, pins) == ["pysam>=0.24.1,<0.25.0", "intron-prospector>=1.5",
                                                      "ucsc-bedtobigbed"]
    assert flair_env.shipped_conda_extras(None, pins) == []


def test_git_bad_ref_message(world):
    world.env["PATH"] = world.env["PATH"]
    p, _ = world.resolve({"id": "g", "flair": {"git": "https://example.org/flair", "ref": "no-such-branch"}})
    assert p.returncode == 1
    assert "'no-such-branch' is not a branch, tag or commit of https://example.org/flair" in p.stderr
    assert "Command:     git ls-remote" in p.stderr


def test_no_conda_message(world):
    world.env.pop("FLAIR_EVAL_CONDA")
    world.env["PATH"] = os.path.dirname(sys.executable) + os.pathsep + "/usr/bin:/bin"
    p, _ = world.resolve({"id": "r", "flair": {"release": "3.0.1"}})
    if p.returncode == 0:
        pytest.skip("a conda executable is on the minimal PATH")
    assert "needs mamba, micromamba or conda" in p.stderr and "FLAIR_EVAL_CONDA" in p.stderr
