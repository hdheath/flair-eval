#!/usr/bin/env python3
"""flair_env.py: turn a run's `flair:` source into a FLAIR environment (docs/CONTRACTS.md section 2).

Usage (contract):
  flair_env.py resolve --run-json RUN.json --cache DIR --out flair_env.json

Sources (runs.yaml `flair:`, exactly one):
  release: "3.0.1"            conda env + `pip install flair-brookslab==3.0.1` (PyPI); bioconda
                              `flair=3.0.1` only when PyPI does not have that version
  git: URL, ref: REF          REF (branch, tag or SHA) -> commit SHA, then that commit is
                              pip-installed into a conda env cached by SHA
  env: /abs/conda/prefix      used as is, after `flair --version` works
  container: URI | /x.sif     used as is, after `apptainer exec URI flair --version` works

Built envs share one recipe so runs differ only in FLAIR: Python picked from the package's
`requires-python` (3.12 when allowed), the external tools pinned in TOOL_PINS, and the runtime-only
conda packages that FLAIR's own misc/flair_conda_env.yaml lists at that commit/tag (e.g. dev needs
ucsc-bedtobigbed; the R stack and variant tools are left out). Envs are
cached as DIR/envs/<key> (key = release-<v> | git-<sha12>), and a file lock (DIR/locks/<key>.lock)
makes concurrent tasks build a key only once. A prefix without its completion marker is a broken
build and is rebuilt.

Writes flair_env.json:
  {"run_id","source","requested","resolved_commit","flair_version","conda_prefix","container",
   "python","pip_freeze_sha256","created"}

Stdlib only. Environment overrides: FLAIR_EVAL_CONDA (conda/mamba executable),
FLAIR_EVAL_TOOL_PINS ("minimap2=2.30 samtools=1.22.1 bedtools=2.31.1"),
FLAIR_EVAL_PYPI (PyPI JSON base URL), FLAIR_EVAL_CONTAINER_RUNTIME (apptainer|singularity).
"""
import argparse
import contextlib
import datetime
import fcntl
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request

ENV_TOOL_VERSION = "1.0.0"
PYPI_PACKAGE = "flair-brookslab"
DEFAULT_TOOL_PINS = "minimap2=2.30 samtools=1.22.1 bedtools=2.31.1"
CHANNELS = ["conda-forge", "bioconda"]
PYTHON_CHOICES = ["3.12", "3.13", "3.11", "3.10"]   # first one requires-python allows wins
MARKER = ".flair-eval-env.json"
THREAD_ENV = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
              "XGBOOST_NTHREAD": "1"}
SHA_RE = re.compile(r"^[0-9a-f]{40}$")


class EnvError(Exception):
    pass


def now_iso():
    return datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat()


def log(msg):
    sys.stderr.write("flair_env.py: %s\n" % msg)
    sys.stderr.flush()


def child_env(prefix=None):
    env = dict(os.environ)
    env.update(THREAD_ENV)
    if prefix:
        env["PATH"] = os.path.join(prefix, "bin") + os.pathsep + env.get("PATH", "")
        env["CONDA_PREFIX"] = prefix
    for k in ("PYTHONPATH", "PYTHONHOME", "PYTHONUSERBASE"):
        env.pop(k, None)
    env["PYTHONNOUSERSITE"] = "1"
    return env


def run(argv, log_fh=None, env=None, cwd=None, timeout=None):
    """Run argv; stream output to log_fh; raise EnvError on failure."""
    cmd = " ".join(argv)
    log("+ " + cmd)
    if log_fh:
        log_fh.write("\n$ %s\n" % cmd)
        log_fh.flush()
    try:
        p = subprocess.run(argv, stdout=log_fh or subprocess.DEVNULL, stderr=subprocess.STDOUT,
                           env=env or child_env(), cwd=cwd, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as e:
        raise EnvError("could not run %s: %s" % (cmd, e))
    if p.returncode != 0:
        raise EnvError("%s exited %d%s" % (cmd, p.returncode,
                                            " (log %s)" % log_fh.name if log_fh else ""))


def capture(argv, env=None, timeout=900):
    try:
        p = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           universal_newlines=True, env=env or child_env(), timeout=timeout)
        return p.returncode, p.stdout, p.stderr
    except (OSError, subprocess.TimeoutExpired) as e:
        return 127, "", str(e)


# --------------------------------------------------------------------------------------------
# Locking
# --------------------------------------------------------------------------------------------
@contextlib.contextmanager
def key_lock(cache, key):
    """Exclusive POSIX lock (fcntl.lockf: honoured across nodes on CephFS/NFS) for one cache key."""
    d = os.path.join(cache, "locks")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, key + ".lock")
    fh = open(path, "a+")
    t = time.time()
    try:
        try:
            fcntl.lockf(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            log("waiting for lock %s (another task is building %s)" % (path, key))
            fcntl.lockf(fh, fcntl.LOCK_EX)
            log("got lock after %.0f s" % (time.time() - t))
        fh.seek(0)
        fh.truncate()
        fh.write("%s pid %d %s\n" % (os.uname()[1], os.getpid(), now_iso()))
        fh.flush()
        yield
    finally:
        try:
            fcntl.lockf(fh, fcntl.LOCK_UN)
        finally:
            fh.close()


# --------------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------------
def conda_exe(explicit=None):
    for c in (explicit, os.environ.get("FLAIR_EVAL_CONDA")):
        if c:
            if shutil.which(c) or os.path.exists(c):
                return shutil.which(c) or c
            raise EnvError("conda executable not found: %s" % c)
    for name in ("mamba", "micromamba", "conda"):
        p = shutil.which(name)
        if p:
            return p
    raise EnvError("no mamba/micromamba/conda on PATH (set FLAIR_EVAL_CONDA)")


def tool_pins(explicit=None):
    s = explicit if explicit is not None else os.environ.get("FLAIR_EVAL_TOOL_PINS", DEFAULT_TOOL_PINS)
    return [t for t in s.replace(",", " ").split() if t]


def version_key(v):
    return tuple(int(x) for x in re.findall(r"\d+", v)[:3])


def python_ok(spec, py):
    """Does python version `py` (e.g. '3.12') satisfy a requires-python spec?"""
    if not spec:
        return True
    pv = version_key(py)
    for clause in [c.strip() for c in spec.split(",") if c.strip()]:
        m = re.match(r"(~=|==|!=|<=|>=|<|>)\s*([\d.*]+)", clause)
        if not m:
            continue
        op, ver = m.groups()
        if "*" in ver:
            want = version_key(ver.replace("*", ""))
            hit = pv[:len(want)] == want
            if (op == "==" and not hit) or (op == "!=" and hit):
                return False
            continue
        vv = version_key(ver)
        n = max(len(vv), 2)
        a, b = (pv + (0, 0, 0))[:n], (vv + (0, 0, 0))[:n]
        ok = {"<": a < b, "<=": a <= b, ">": a > b, ">=": a >= b, "==": a == b, "!=": a != b,
              "~=": a >= b and a[:len(vv) - 1] == b[:len(vv) - 1]}[op]
        if not ok:
            return False
    return True


def pick_python(spec):
    for py in PYTHON_CHOICES:
        if python_ok(spec, py):
            return py
    raise EnvError("no supported Python (%s) satisfies requires-python %r" % (PYTHON_CHOICES, spec))


def pypi_info(version):
    """(requires_python, True) if flair-brookslab==version is on PyPI, (None, False) if not.
    Raises EnvError when PyPI cannot be reached."""
    base = os.environ.get("FLAIR_EVAL_PYPI", "https://pypi.org/pypi").rstrip("/")
    url = "%s/%s/%s/json" % (base, PYPI_PACKAGE, version)
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            d = json.load(r)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None, False
        raise EnvError("PyPI query failed (%s): %s" % (url, e))
    except (urllib.error.URLError, OSError, ValueError) as e:
        raise EnvError("PyPI query failed (%s): %s" % (url, e))
    return (d.get("info") or {}).get("requires_python"), True


def flair_version(prefix=None, container=None, runtime=None):
    if container:
        argv = [runtime, "exec", container, "flair", "--version"]
        env = child_env()
    else:
        argv = [os.path.join(prefix, "bin", "flair"), "--version"]
        env = child_env(prefix)
    code, out, err = capture(argv, env=env)
    text = (out + "\n" + err).strip()
    m = None
    for line in text.splitlines()[::-1]:
        m = re.search(r"FLAIR\s+v?(\S+)", line)
        if m:
            break
    if code != 0 or not m:
        raise EnvError("`%s` failed (exit %d): %s" % (" ".join(argv), code, text[-800:]))
    return "FLAIR " + m.group(1)


def python_and_freeze(prefix=None, container=None, runtime=None):
    def call(args):
        if container:
            return capture([runtime, "exec", container] + args, env=child_env())
        py = os.path.join(prefix, "bin", "python")
        if not os.path.exists(py):
            return 127, "", "no python in env"
        return capture([py] + args[1:], env=child_env(prefix))

    code, out, _ = call(["python", "-c", "import platform; print(platform.python_version())"])
    pyver = out.strip() if code == 0 and out.strip() else None
    code, out, _ = call(["python", "-m", "pip", "freeze"])
    freeze = out if code == 0 else None
    sha = hashlib.sha256(freeze.encode()).hexdigest() if freeze is not None else None
    return pyver, sha, freeze


def detect_commit(prefix):
    """FLAIR's commit for a pass-through env: pip's direct_url.json (vcs or local dir), or the
    git checkout of an editable install. None when unknown (e.g. a PyPI wheel)."""
    sp = []
    for root, dirs, _ in os.walk(os.path.join(prefix, "lib")):
        if root.count(os.sep) - prefix.count(os.sep) > 3:
            dirs[:] = []
            continue
        if os.path.basename(root) == "site-packages":
            sp.append(root)
            dirs[:] = []
    for d in sp:
        for name in sorted(os.listdir(d)):
            if not (name.lower().startswith(("flair_brookslab", "flair-")) and name.endswith(".dist-info")):
                continue
            du = os.path.join(d, name, "direct_url.json")
            if not os.path.exists(du):
                continue
            try:
                info = json.load(open(du))
            except (OSError, ValueError):
                continue
            commit = (info.get("vcs_info") or {}).get("commit_id")
            if commit:
                return commit
            url = info.get("url", "")
            if url.startswith("file://"):
                src = url[len("file://"):]
                code, out, _ = capture(["git", "-C", src, "rev-parse", "HEAD"])
                if code == 0 and SHA_RE.match(out.strip()):
                    c2, st, _ = capture(["git", "-C", src, "status", "--porcelain", "--untracked-files=no"])
                    if c2 == 0 and st.strip():
                        log("note: %s has uncommitted changes; resolved_commit is its HEAD" % src)
                    return out.strip()
    return None


def read_marker(prefix):
    p = os.path.join(prefix, MARKER)
    if not os.path.exists(p):
        return None
    try:
        return json.load(open(p))
    except (OSError, ValueError):
        return None


def write_marker(prefix, data):
    tmp = os.path.join(prefix, MARKER + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=2)
        fh.write("\n")
    os.replace(tmp, os.path.join(prefix, MARKER))


def conda_create(conda, prefix, packages, log_fh):
    argv = [conda, "create", "-y", "-p", prefix, "--override-channels"]
    for c in CHANNELS:
        argv += ["-c", c]
    run(argv + packages, log_fh=log_fh, env=build_env())


def build_env():
    env = child_env()
    env.setdefault("PIP_NO_CACHE_DIR", "1")
    env["PIP_DISABLE_PIP_VERSION_CHECK"] = "1"
    return env


def pip_install(prefix, target, log_fh):
    run([os.path.join(prefix, "bin", "python"), "-m", "pip", "install", "--no-cache-dir", target],
        log_fh=log_fh, env=child_env(prefix))


# --------------------------------------------------------------------------------------------
# Git
# --------------------------------------------------------------------------------------------
def resolve_git_ref(url, ref):
    """Commit SHA for a branch/tag/SHA without cloning when possible."""
    if SHA_RE.match(ref.lower()):
        return ref.lower()
    code, out, err = capture(["git", "ls-remote", url, ref, "refs/heads/" + ref, "refs/tags/" + ref,
                              "refs/tags/%s^{}" % ref], timeout=300)
    if code != 0:
        raise EnvError("git ls-remote %s failed: %s" % (url, err.strip()[-500:]))
    refs = {}
    for line in out.splitlines():
        parts = line.split()
        if len(parts) == 2:
            refs[parts[1]] = parts[0]
    for name in ("refs/tags/%s^{}" % ref, "refs/tags/" + ref, "refs/heads/" + ref, ref):
        if name in refs:
            return refs[name]
    if re.match(r"^[0-9a-fA-F]{7,39}$", ref):
        return None          # short SHA: needs a clone (resolve_short_sha)
    raise EnvError("ref %r not found in %s" % (ref, url))


def fetch_source(url, sha, src, log_fh):
    """Check out `sha` of `url` into `src` (shallow fetch by SHA, full clone as fallback)."""
    if os.path.isdir(src):
        shutil.rmtree(src)
    os.makedirs(src)
    try:
        run(["git", "init", "-q", src], log_fh=log_fh)
        run(["git", "-C", src, "remote", "add", "origin", url], log_fh=log_fh)
        run(["git", "-C", src, "fetch", "-q", "--depth", "1", "origin", sha], log_fh=log_fh, timeout=1800)
        run(["git", "-C", src, "checkout", "-q", "FETCH_HEAD"], log_fh=log_fh)
    except EnvError:
        log("shallow fetch by SHA failed; cloning")
        shutil.rmtree(src, ignore_errors=True)
        run(["git", "clone", "-q", url, src], log_fh=log_fh, timeout=3600)
        run(["git", "-C", src, "checkout", "-q", sha], log_fh=log_fh)
    code, out, _ = capture(["git", "-C", src, "rev-parse", "HEAD"])
    if code != 0 or out.strip() != sha:
        raise EnvError("checkout of %s gave %s, wanted %s" % (url, out.strip(), sha))


def resolve_short_sha(url, ref, cache, log_fh):
    tmp = os.path.join(cache, "src", "_resolve-%d" % os.getpid())
    shutil.rmtree(tmp, ignore_errors=True)
    try:
        run(["git", "clone", "-q", "--filter=blob:none", "--no-checkout", url, tmp], log_fh=log_fh, timeout=3600)
        code, out, err = capture(["git", "-C", tmp, "rev-parse", "--verify", ref + "^{commit}"])
        if code != 0:
            raise EnvError("ref %r not found in %s" % (ref, url))
        return out.strip()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# conda packages in FLAIR's misc/flair_conda_env.yaml that flair-eval does not need: the R stack
# (diffexp/diffsplice), variant tools, build/doc tooling, and what TOOL_PINS / pip already give.
SHIPPED_SKIP = re.compile(r"^(python|pip|r-.*|bioconductor-.*|longshot|whatshap|poetry|sphinx.*|"
                          r"ucsc-hubcheck|pytest.*|flake8)$")
GITHUB_RAW = "https://raw.githubusercontent.com/BrooksLabUCSC/flair/%s/misc/flair_conda_env.yaml"


def shipped_conda_extras(yaml_text, pins):
    """Runtime conda packages a FLAIR checkout declares beyond our base recipe (e.g. dev needs
    ucsc-bedtobigbed, which pip cannot install). Parses the top-level `dependencies:` list only."""
    pinned = {re.split(r"[=<>! ]", p, maxsplit=1)[0] for p in pins}
    extras, in_deps = [], False
    for raw in (yaml_text or "").splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        if re.match(r"^dependencies\s*:", line):
            in_deps = True
            continue
        if in_deps and not line.startswith(" "):
            break
        m = re.match(r"^\s{0,4}-\s+(.+)$", line) if in_deps else None
        if not m or m.group(1).rstrip().endswith(":"):
            continue
        spec = m.group(1).strip().strip("'\"")
        name = re.split(r"[=<>! ]", spec, maxsplit=1)[0]
        if SHIPPED_SKIP.match(name) or name in pinned:
            continue
        extras.append(spec.replace(" ", ""))
    return extras


def shipped_yaml_for_release(version):
    """misc/flair_conda_env.yaml at tag v<version>, or None when it cannot be fetched."""
    url = os.environ.get("FLAIR_EVAL_RELEASE_YAML", GITHUB_RAW) % ("v" + version)
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            return r.read().decode()
    except (urllib.error.URLError, OSError, ValueError):
        return None


def requires_python_of(src):
    for name in ("pyproject.toml", "setup.cfg"):
        p = os.path.join(src, name)
        if os.path.exists(p):
            m = re.search(r"""requires[-_]python\s*[=:]\s*["']?([^"'\n]+)["']?""", open(p).read())
            if m:
                return m.group(1).strip()
    p = os.path.join(src, "setup.py")
    if os.path.exists(p):
        m = re.search(r"""python_requires\s*=\s*["']([^"']+)["']""", open(p).read())
        if m:
            return m.group(1)
    return None


# --------------------------------------------------------------------------------------------
# Resolution per source
# --------------------------------------------------------------------------------------------
def ensure_built(cache, key, builder, recipe):
    """Return prefix of a complete env for `key`, building it once under the key's lock."""
    prefix = os.path.join(cache, "envs", key)
    marker = read_marker(prefix)
    if marker and marker.get("status") == "complete":
        return prefix, marker, False
    with key_lock(cache, key):
        marker = read_marker(prefix)
        if marker and marker.get("status") == "complete":
            log("env %s was built by another task meanwhile" % key)
            return prefix, marker, False
        if os.path.exists(prefix):
            log("removing incomplete env %s" % prefix)
            shutil.rmtree(prefix)
        os.makedirs(os.path.join(cache, "logs"), exist_ok=True)
        log_path = os.path.join(cache, "logs", "%s.build.log" % key)
        t0 = time.time()
        with open(log_path, "w") as log_fh:
            log_fh.write("# flair_env.py %s build of %s at %s\n# recipe %s\n"
                         % (ENV_TOOL_VERSION, key, now_iso(), json.dumps(recipe)))
            try:
                extra = builder(prefix, log_fh) or {}
            except Exception:
                log("build failed; see %s" % log_path)
                raise
        marker = dict(recipe)
        marker.update(extra)
        marker.update({"status": "complete", "key": key, "built": now_iso(),
                       "build_seconds": round(time.time() - t0, 1), "build_log": log_path,
                       "flair_env_tool": ENV_TOOL_VERSION})
        write_marker(prefix, marker)
        return prefix, marker, True


def resolve_release(version, cache, opts):
    key = "release-%s" % re.sub(r"[^A-Za-z0-9._+-]", "_", version)
    pins = tool_pins(opts.tool_pins)
    prefix = os.path.join(cache, "envs", key)
    marker = read_marker(prefix)
    if not (marker and marker.get("status") == "complete"):
        # decide the route before locking (network only; cheap)
        try:
            req_py, on_pypi = pypi_info(version)
        except EnvError as e:
            if opts.route == "pip":
                raise
            log("%s; falling back to bioconda" % e)
            req_py, on_pypi = None, False
        route = opts.route if opts.route != "auto" else ("pip" if on_pypi else "bioconda")
        if route == "pip" and not on_pypi:
            raise EnvError("%s==%s is not on PyPI" % (PYPI_PACKAGE, version))
        py = pick_python(req_py) if route == "pip" else None
        extras = shipped_conda_extras(shipped_yaml_for_release(version), pins) if route == "pip" else []
        recipe = {"source": "release", "requested": version, "route": route, "python_spec": py,
                  "tool_pins": pins, "conda_extras": extras, "channels": CHANNELS}
    else:
        recipe = {}

    conda = conda_exe(opts.conda) if not (marker and marker.get("status") == "complete") else None

    def builder(pfx, log_fh):
        if recipe["route"] == "pip":
            conda_create(conda, pfx, ["python=%s" % recipe["python_spec"]] + pins + recipe["conda_extras"]
                         + ["pip"], log_fh)
            pip_install(pfx, "%s==%s" % (PYPI_PACKAGE, version), log_fh)
        else:
            conda_create(conda, pfx, ["flair=%s" % version] + pins, log_fh)
        v = flair_version(prefix=pfx)
        if version not in v:
            raise EnvError("built env reports %r, wanted %s" % (v, version))
        return {}

    prefix, marker, built = ensure_built(cache, key, builder, recipe)
    return {"source": "release", "requested": version, "resolved_commit": None,
            "conda_prefix": prefix, "container": None, "_built": built, "_marker": marker}


def resolve_git(url, ref, cache, opts):
    sha = resolve_git_ref(url, ref)
    if sha is None:
        os.makedirs(os.path.join(cache, "logs"), exist_ok=True)
        with open(os.path.join(cache, "logs", "resolve-%d.log" % os.getpid()), "w") as fh:
            sha = resolve_short_sha(url, ref, cache, fh)
    key = "git-%s" % sha[:12]
    pins = tool_pins(opts.tool_pins)
    recipe = {"source": "git", "requested": "%s@%s" % (url, ref), "git_url": url, "git_ref": ref,
              "commit": sha, "tool_pins": pins, "channels": CHANNELS}
    prefix = os.path.join(cache, "envs", key)
    marker = read_marker(prefix)
    conda = conda_exe(opts.conda) if not (marker and marker.get("status") == "complete") else None

    def builder(pfx, log_fh):
        src = os.path.join(cache, "src", sha)
        with key_lock(cache, "src-" + sha[:12]):
            if not (os.path.isdir(src) and capture(["git", "-C", src, "rev-parse", "HEAD"])[1].strip() == sha):
                fetch_source(url, sha, src, log_fh)
        spec = requires_python_of(src)
        py = pick_python(spec)
        shipped = os.path.join(src, "misc", "flair_conda_env.yaml")
        extras = shipped_conda_extras(open(shipped).read(), pins) if os.path.exists(shipped) else []
        conda_create(conda, pfx, ["python=%s" % py] + pins + extras + ["pip", "git"], log_fh)
        pip_install(pfx, src, log_fh)
        flair_version(prefix=pfx)
        return {"python_spec": py, "requires_python": spec, "src": src, "conda_extras": extras}

    prefix, marker, built = ensure_built(cache, key, builder, recipe)
    return {"source": "git", "requested": "%s@%s" % (url, ref), "resolved_commit": sha,
            "conda_prefix": prefix, "container": None, "_built": built, "_marker": marker}


def resolve_env(path, cache, opts):
    prefix = os.path.abspath(os.path.expanduser(path))
    if not os.path.isdir(prefix):
        raise EnvError("env not found: %s" % prefix)
    if not os.path.exists(os.path.join(prefix, "bin", "flair")):
        raise EnvError("no bin/flair in env %s" % prefix)
    return {"source": "env", "requested": path, "resolved_commit": detect_commit(prefix),
            "conda_prefix": prefix, "container": None}


def container_runtime(explicit=None):
    for c in (explicit, os.environ.get("FLAIR_EVAL_CONTAINER_RUNTIME"), "apptainer", "singularity"):
        if c and shutil.which(c):
            return shutil.which(c)
    raise EnvError("no apptainer/singularity on PATH for a container source")


def resolve_container(uri, cache, opts):
    if not re.match(r"^[a-z]+://", uri):
        p = os.path.abspath(os.path.expanduser(uri))
        if not os.path.exists(p):
            raise EnvError("container image not found: %s" % p)
        uri = p
    return {"source": "container", "requested": uri, "resolved_commit": None,
            "conda_prefix": None, "container": uri}


RESOLVERS = {"release": resolve_release, "git": resolve_git, "env": resolve_env,
             "container": resolve_container}


def resolve(run_map, cache, opts):
    run_id = run_map.get("id")
    if not run_id:
        raise EnvError("run has no id")
    flair = run_map.get("flair")
    if not isinstance(flair, dict):
        raise EnvError("run %s: `flair` must be a map with one of %s" % (run_id, sorted(RESOLVERS)))
    sources = [k for k in RESOLVERS if k in flair]
    extra = sorted(set(flair) - set(RESOLVERS) - {"ref"})
    if len(sources) != 1 or extra:
        raise EnvError("run %s: `flair` needs exactly one of release|git|env|container, got %s"
                       % (run_id, sorted(flair)))
    src = sources[0]
    val = flair[src]
    if not isinstance(val, str) or not val:
        raise EnvError("run %s: flair.%s must be a non-empty string" % (run_id, src))
    if src == "git":
        ref = flair.get("ref")
        if not isinstance(ref, str) or not ref:
            raise EnvError("run %s: flair.git needs `ref`" % run_id)
        res = resolve_git(val, ref, cache, opts)
    elif "ref" in flair:
        raise EnvError("run %s: `ref` only goes with `git`" % run_id)
    else:
        res = RESOLVERS[src](val, cache, opts)

    runtime = container_runtime(opts.container_runtime) if res["container"] else None
    fv = flair_version(prefix=res["conda_prefix"], container=res["container"], runtime=runtime)
    pyver, freeze_sha, freeze = python_and_freeze(prefix=res["conda_prefix"], container=res["container"],
                                                  runtime=runtime)
    if res.get("_built") and freeze is not None and res["conda_prefix"]:
        with open(os.path.join(res["conda_prefix"], ".flair-eval-pip-freeze.txt"), "w") as fh:
            fh.write(freeze)
    return {
        "run_id": run_id,
        "source": res["source"],
        "requested": res["requested"],
        "resolved_commit": res["resolved_commit"],
        "flair_version": fv,
        "conda_prefix": res["conda_prefix"],
        "container": res["container"],
        "python": pyver,
        "pip_freeze_sha256": freeze_sha,
        "created": now_iso(),
    }


def cmd_resolve(a):
    try:
        with open(a.run_json) as fh:
            run_map = json.load(fh)
    except (OSError, ValueError) as e:
        log("error: cannot read --run-json %s: %s" % (a.run_json, e))
        return 2
    cache = os.path.abspath(os.path.expanduser(a.cache))
    try:
        os.makedirs(cache, exist_ok=True)
        if "CONDA_PKGS_DIRS" not in os.environ:
            # keep conda's package cache off $HOME (quota) and next to the envs it serves
            os.environ["CONDA_PKGS_DIRS"] = os.path.join(cache, "conda-pkgs")
        rec = resolve(run_map, cache, a)
    except EnvError as e:
        log("error: %s" % e)
        return 1
    tmp = a.out + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(rec, fh, indent=2)
        fh.write("\n")
    os.replace(tmp, a.out)
    log("%s -> %s (%s)" % (rec["run_id"], rec["conda_prefix"] or rec["container"], rec["flair_version"]))
    return 0


def build_parser():
    p = argparse.ArgumentParser(prog="flair_env.py", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd")
    r = sub.add_parser("resolve", help="resolve a run's flair source to an env; write flair_env.json")
    r.add_argument("--run-json", required=True, help="the run map (from runs.yaml) as JSON")
    r.add_argument("--cache", required=True, help="env cache directory (params.env_cache)")
    r.add_argument("--out", required=True, help="flair_env.json to write")
    r.add_argument("--conda", default=None, help="conda/mamba executable (default: mamba, micromamba, conda)")
    r.add_argument("--tool-pins", default=None, help='conda specs added to built envs (default "%s")'
                   % DEFAULT_TOOL_PINS)
    r.add_argument("--route", choices=["auto", "pip", "bioconda"], default="auto",
                   help="release build route: auto = PyPI if it has the version, else bioconda")
    r.add_argument("--container-runtime", default=None)
    return p


def main(argv=None):
    p = build_parser()
    a = p.parse_args(argv)
    if a.cmd != "resolve":
        p.print_help()
        return 2
    return cmd_resolve(a)


if __name__ == "__main__":
    sys.exit(main())
