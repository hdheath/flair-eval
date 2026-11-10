import pytest

from flair_eval.tokens import DEFAULT_GIT_URL, TokenError, parse_token, parse_tokens, \
    runs_from_tokens, token_for_flair, token_for_run

# Every row of the CONTRACTS.md section 1 token table.
CONTRACT_TABLE = [
    ("3.0.1", {"release": "3.0.1"}, "flair-3.0.1"),
    ("release:3.0.1", {"release": "3.0.1"}, "flair-3.0.1"),
    ("dev", {"git": DEFAULT_GIT_URL, "ref": "dev"}, "flair-dev"),
    ("main", {"git": DEFAULT_GIT_URL, "ref": "main"}, "flair-main"),
    ("git:ted-port", {"git": DEFAULT_GIT_URL, "ref": "ted-port"}, "flair-ted-port"),
    ("git:https://github.com/me/flair@my-branch",
     {"git": "https://github.com/me/flair", "ref": "my-branch"}, "flair-my-branch"),
    ("env:/path/to/conda/envs/flair-3.0.1", {"env": "/path/to/conda/envs/flair-3.0.1"},
     "env-flair-3.0.1"),
    ("container:docker://quay.io/biocontainers/flair:3.0.1--pyhdfd78af_0",
     {"container": "docker://quay.io/biocontainers/flair:3.0.1--pyhdfd78af_0"}, "ctr-flair"),
    ("container:/images/flair.sif", {"container": "/images/flair.sif"}, "ctr-flair.sif"),
]


@pytest.mark.parametrize("token,flair,rid", CONTRACT_TABLE)
def test_contract_table(token, flair, rid):
    r = parse_token(token)
    assert (r.flair, r.id) == (flair, rid)


@pytest.mark.parametrize("token,flair,rid", [
    ("git:git@github.com:me/flair.git@abc123",
     {"git": "git@github.com:me/flair.git", "ref": "abc123"}, "flair-abc123"),
    ("git:/home/me/flair@dev", {"git": "/home/me/flair", "ref": "dev"}, "flair-dev"),
    ("git:feature/ends", {"git": DEFAULT_GIT_URL, "ref": "feature/ends"}, "flair-feature-ends"),
    ("env:/x/myenv/", {"env": "/x/myenv/"}, "env-myenv"),
    ("container:docker://org/flair@sha256:abcd", {"container": "docker://org/flair@sha256:abcd"},
     "ctr-flair"),
    ("container:docker://flair:3.0.1", {"container": "docker://flair:3.0.1"}, "ctr-flair"),
    ("  2.2.0 ", {"release": "2.2.0"}, "flair-2.2.0"),
    ("release:v3", {"release": "v3"}, "flair-v3"),
])
def test_open_details(token, flair, rid):
    r = parse_token(token)
    assert (r.flair, r.id) == (flair, rid)


@pytest.mark.parametrize("bad", ["", "  ", "master", "release:", "git:", "env:", "container:",
                                 "git:https://github.com/me/flair",
                                 "git:git@github.com:me/flair.git", "env:/"])
def test_bad_tokens(bad):
    with pytest.raises(TokenError):
        parse_token(bad)


def test_unknown_word_suggests_git():
    with pytest.raises(TokenError, match="git:master"):
        parse_token("master")


def test_parse_tokens_list_and_duplicates():
    assert [r.id for r in parse_tokens("3.0.1,dev")] == ["flair-3.0.1", "flair-dev"]
    with pytest.raises(TokenError, match="both give run id"):
        parse_tokens("3.0.1,release:3.0.1")
    with pytest.raises(TokenError, match="stray comma"):
        parse_tokens("3.0.1,,dev")
    with pytest.raises(TokenError):
        parse_tokens("")


def test_runs_from_tokens_applies_flair_args():
    runs = runs_from_tokens("3.0.1,dev", "--stringent")
    assert [r["transcriptome_args"] for r in runs] == ["--stringent", "--stringent"]
    assert runs[1]["flair"] == {"git": DEFAULT_GIT_URL, "ref": "dev"}


@pytest.mark.parametrize("token,flair,rid", CONTRACT_TABLE)
def test_round_trip(token, flair, rid):
    tok = token_for_flair(flair)
    assert parse_token(tok).flair == flair
    assert token_for_run({"id": rid, "flair": flair}) == tok


def test_token_for_run_rejects_other_ids_and_sources():
    assert token_for_run({"id": "my-name", "flair": {"release": "3.0.1"}}) is None
    assert token_for_flair({"release": 3.0}) is None  # unquoted YAML number
    assert token_for_flair({"release": "3.0.1", "env": "/x"}) is None
    assert token_for_flair({"env": "/a,b"}) is None
    assert token_for_flair({}) is None
    assert token_for_flair({"git": DEFAULT_GIT_URL, "ref": "dev"}) == "dev"
    assert token_for_flair({"git": DEFAULT_GIT_URL, "ref": "x"}) == "git:x"
