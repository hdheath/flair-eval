"""The `--flair` token mapping (docs/CONTRACTS.md section 1, "One-command mode").

Each comma-separated token becomes one run:

    3.0.1 | release:3.0.1     -> {release: "3.0.1"}                       id flair-3.0.1
    dev | main | git:<ref>    -> {git: DEFAULT_GIT_URL, ref: <ref>}        id flair-<ref>
    git:<url>@<ref>           -> {git: <url>, ref: <ref>}                  id flair-<ref>
    env:/path                 -> {env: "/path"}                            id env-<basename>
    container:<uri>           -> {container: <uri>}                        id ctr-<basename-without-tag>

Details the contract leaves open, decided here (and mirrored in frontend/builder/logic.js):
- A bare token is a release when it starts with a digit; `dev` and `main` are git refs;
  any other bare word is an error that suggests `git:<word>`.
- In `git:<x>`, x is `<url>@<ref>` when the part before the last `@` looks like a URL
  (contains `://`, starts with `git@`, `/` or `file:`), otherwise x is a ref on the default URL.
- Characters outside `[A-Za-z0-9._-]` in a derived id are replaced by `-`
  (so `git:feature/x` gives id `flair-feature-x`).
- basename = text after the last `/` (trailing `/` ignored). For containers the tag is
  `:<tag>` and/or `@<digest>` on the basename; `/x/flair.sif` keeps `.sif` (id `ctr-flair.sif`).
"""

from __future__ import annotations

import re
from dataclasses import dataclass

DEFAULT_GIT_URL = "https://github.com/BrooksLabUCSC/flair"
GIT_SHORTCUTS = ("dev", "main")
ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class TokenError(ValueError):
    """A --flair token that cannot be mapped to a run."""


@dataclass(frozen=True)
class TokenRun:
    token: str
    id: str
    flair: dict


def _sanitize(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "-", s)


def _basename(path: str) -> str:
    p = path.rstrip("/")
    return p.rsplit("/", 1)[-1]


def _looks_like_url(s: str) -> bool:
    return "://" in s or s.startswith(("git@", "/", "file:"))


def _container_name(uri: str) -> str:
    base = _basename(uri.split("://", 1)[-1])
    base = base.split("@", 1)[0]  # digest
    if ":" in base:  # tag
        base = base.rsplit(":", 1)[0]
    return base


def _make_id(prefix: str, rest: str, token: str) -> str:
    rid = prefix + _sanitize(rest)
    if not ID_RE.fullmatch(rid) or not rest:
        raise TokenError(f"--flair token {token!r} gives an unusable run id {rid!r}; "
                         "write this run in runs.yaml instead")
    return rid


def parse_token(token: str) -> TokenRun:
    """Map one --flair token to (id, flair source)."""
    t = token.strip()
    if not t:
        raise TokenError("empty --flair token (check for a stray comma)")

    if t.startswith("release:"):
        v = t[len("release:"):].strip()
        if not v:
            raise TokenError(f"--flair token {t!r}: give a version, e.g. release:3.0.1")
        return TokenRun(t, _make_id("flair-", v, t), {"release": v})

    if t.startswith("git:"):
        rest = t[len("git:"):].strip()
        if not rest:
            raise TokenError(f"--flair token {t!r}: give a branch, tag or SHA, e.g. git:dev")
        url, ref = DEFAULT_GIT_URL, rest
        if _looks_like_url(rest):
            at = rest.rfind("@")
            url, ref = (rest[:at], rest[at + 1:]) if at > 0 else ("", "")
            if not _looks_like_url(url) or not ref or "/" in ref and ":" in ref:
                raise TokenError(f"--flair token {t!r}: a git URL needs @<ref>, "
                                 "e.g. git:https://github.com/me/flair@my-branch")
        return TokenRun(t, _make_id("flair-", ref, t), {"git": url, "ref": ref})

    if t.startswith("env:"):
        p = t[len("env:"):].strip()
        if not p:
            raise TokenError(f"--flair token {t!r}: give a conda env path, e.g. env:/path/to/env")
        return TokenRun(t, _make_id("env-", _basename(p), t), {"env": p})

    if t.startswith("container:"):
        uri = t[len("container:"):].strip()
        if not uri:
            raise TokenError(f"--flair token {t!r}: give an image, e.g. "
                             "container:docker://org/flair:3.0.1 or container:/path/flair.sif")
        return TokenRun(t, _make_id("ctr-", _container_name(uri), t), {"container": uri})

    if t in GIT_SHORTCUTS:
        return TokenRun(t, "flair-" + t, {"git": DEFAULT_GIT_URL, "ref": t})

    if re.match(r"[0-9]", t):
        return TokenRun(t, _make_id("flair-", t, t), {"release": t})

    raise TokenError(f"unknown --flair token {t!r}. Use a release (3.0.1), dev, main, "
                     f"git:{t} for another branch, env:/path or container:<uri>")


def parse_tokens(value: str) -> list[TokenRun]:
    """Parse a comma-separated --flair value. Run ids must come out unique."""
    if value is None or not str(value).strip():
        raise TokenError("--flair is empty; give at least one token, e.g. --flair 3.0.1")
    runs = [parse_token(t) for t in str(value).split(",")]
    seen: dict[str, str] = {}
    for r in runs:
        if r.id in seen:
            raise TokenError(f"--flair tokens {seen[r.id]!r} and {r.token!r} both give run id "
                             f"{r.id!r}; write these runs in runs.yaml with distinct ids")
        seen[r.id] = r.token
    return runs


def runs_from_tokens(value: str, flair_args: str = "") -> list[dict]:
    """The runs.yaml-style run maps that a --flair value (+ --flair_args) stands for."""
    return [{"id": r.id, "flair": dict(r.flair), "transcriptome_args": flair_args or ""}
            for r in parse_tokens(value)]


def token_for_flair(flair: dict) -> str | None:
    """Shortest token whose flair source equals `flair`, or None if no token can express it."""
    if not isinstance(flair, dict) or len(flair) == 0:
        return None
    if set(flair) == {"release"}:
        v = str(flair["release"])
        tok = v if re.match(r"[0-9]", v) else "release:" + v
    elif set(flair) == {"git", "ref"}:
        url, ref = str(flair["git"]), str(flair["ref"])
        if url == DEFAULT_GIT_URL:
            tok = ref if ref in GIT_SHORTCUTS else "git:" + ref
        else:
            tok = f"git:{url}@{ref}"
    elif set(flair) == {"env"}:
        tok = "env:" + str(flair["env"])
    elif set(flair) == {"container"}:
        tok = "container:" + str(flair["container"])
    else:
        return None
    if "," in tok or tok != tok.strip():
        return None
    try:
        back = parse_token(tok)
    except TokenError:
        return None
    return tok if back.flair == flair else None


def token_for_run(run: dict) -> str | None:
    """Token for a run when --flair would recreate exactly this run id and source."""
    tok = token_for_flair(run.get("flair"))
    if tok is None:
        return None
    return tok if parse_token(tok).id == run.get("id") else None
