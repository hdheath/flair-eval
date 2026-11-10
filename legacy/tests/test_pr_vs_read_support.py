#!/usr/bin/env python3
"""Unit tests for analysis/pr_vs_read_support/pr_vs_read_support.py.

Pins the read-support parsing — "first numeric column after the id" — across the
counts formats in this project (FLAIR two-count, isoseq/isoquant/lraa #id/count,
bambu TXNAME/GENEID/sample with a non-numeric GENEID column).
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_PATH = (Path(__file__).resolve().parent.parent
         / "analysis" / "pr_vs_read_support" / "pr_vs_read_support.py")


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("pr_vs_read_support", _PATH)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_flair_two_count_uses_first_numeric(mod, tmp_path):
    # FLAIR: id \t count \t count2 (no header) — support = first count
    p = _write(tmp_path, "flair.counts.txt",
               "FLT00000001\t489\t494\nFLT00000002\t1\t1\n")
    c = mod.load_counts(p)
    assert c["FLT00000001"] == 489.0
    assert c["FLT00000002"] == 1.0


def test_isoseq_isoquant_header_and_single_count(mod, tmp_path):
    p = _write(tmp_path, "iso.tsv", "#feature_id\tcount\nPB.1.13\t6\nENST9\t0.00\n")
    c = mod.load_counts(p)
    assert c["PB.1.13"] == 6.0
    assert c["ENST9"] == 0.0


def test_bambu_skips_nonnumeric_geneid(mod, tmp_path):
    # bambu: TXNAME \t GENEID \t sample — GENEID is non-numeric, skip to sample
    p = _write(tmp_path, "bambu.txt",
               "TXNAME\tGENEID\tA549\nBambuTx1\tENSG00000228794.13\t27.33567\n")
    c = mod.load_counts(p)
    assert c["BambuTx1"] == pytest.approx(27.33567)


def test_threshold_filtering_via_support(mod, tmp_path):
    p = _write(tmp_path, "c.tsv", "#feature_id\tcount\na\t1\nb\t3\nc\t10\n")
    c = mod.load_counts(p)
    for T, expect in [(1, {"a", "b", "c"}), (3, {"b", "c"}), (5, {"c"})]:
        kept = {k for k, v in c.items() if v >= T}
        assert kept == expect, (T, kept)
