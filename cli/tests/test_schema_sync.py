"""The bundled schema copies (pip package, builder page) must match assets/schema."""
import importlib.util
import json

import pytest

from flair_eval.validate import schema


def _sync(repo):
    p = repo / "cli" / "tools" / "sync_schemas.py"
    spec = importlib.util.spec_from_file_location("sync_schemas", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_copies_are_current(repo):
    if not (repo / "assets" / "schema").is_dir():
        pytest.skip("not in the pipeline repository")
    assert _sync(repo).main(["--check"]) == 0, "run: python cli/tools/sync_schemas.py"


def test_package_schemas_load(repo):
    for name in ("samples", "runs", "bundle"):
        s = schema(name)
        assert s["$id"].startswith("flair-eval/")
        src = repo / "assets" / "schema" / f"{name}.schema.json"
        if src.exists():
            assert s == json.loads(src.read_text())
