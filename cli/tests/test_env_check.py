import subprocess

import pytest

from flair_eval.env_check import check_tools, java_major


@pytest.mark.parametrize("text,major", [
    ('openjdk version "17.0.2" 2022-01-18', 17), ('java version "1.8.0_292"', 8),
    ('openjdk version "21" 2023-09-19', 21), ("garbage", None),
])
def test_java_major(text, major):
    assert java_major(text) == major


def fake_run(stderr):
    def run(*a, **k):
        return subprocess.CompletedProcess(a, 0, stdout="", stderr=stderr)
    return run


def test_all_present():
    errs, warns = check_tools(which=lambda x: "/usr/bin/" + x,
                              run=fake_run('openjdk version "17.0.9"'))
    assert errs == [] and warns == []


def test_missing_tools_give_install_hints():
    errs, _ = check_tools(which=lambda x: None)
    assert len(errs) == 2
    assert "openjdk" in errs[0] and "get.nextflow.io" in errs[1]


def test_old_java_warns():
    errs, warns = check_tools(which=lambda x: "/usr/bin/" + x,
                              run=fake_run('java version "1.8.0_292"'))
    assert errs == [] and "Java 8 found" in warns[0]
