"""Tests for build_truth_union.py: a line's 5' truth union as one truth BED (ted-clustering map ticket 05).

Every call of the higher-tier set is kept at its reference point. A call of the next set pairs
one-to-one with a kept call within 50 bp (the sweep of the ortho-truth agreement baselines) and
is then the same site; an unpaired call is added. Tiers are given in order (CAGE, PRO-cap,
RAMPAGE). Coordinates are BED-style; a truth call is a 1-bp reference point."""

import os
import subprocess
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
VENDORED = os.path.join(ROOT, "modules", "end_accuracy", "resources", "usr", "bin")
BUILDER = os.path.join(VENDORED, "build_truth_union.py")
sys.path.insert(0, VENDORED)

import build_truth_union as union  # noqa: E402

HEADER = "#chrom\tstart\tend\tcall_id\tscore\tstrand\text_start\text_end"


def truth(tmp_path, name, calls):
    """A truth BED of (chrom, strand, point) calls, with a header like the ortho-truth files."""
    path = tmp_path / f"{name}.truth.bed"
    lines = [HEADER] + [f"{c}\t{p}\t{p + 1}\t{name}_{i}\t1\t{s}\t{p}\t{p + 1}"
                        for i, (c, s, p) in enumerate(calls, 1)]
    path.write_text("".join(line + "\n" for line in lines))
    return str(path)


def merged(tmp_path, *tiers):
    """The union of (name, calls) tiers as (set, chrom, strand, point) tuples."""
    built = union.build_union([(name, union.read_calls(truth(tmp_path, name, calls)))
                               for name, calls in tiers], window=50)
    return [(c.set, c.chrom, c.strand, c.point) for c in built]


def test_a_site_seen_by_both_sets_is_one_call(tmp_path):
    assert merged(tmp_path, ("CAGE", [("chr1", "+", 1000)]), ("PROcap", [("chr1", "+", 1030)])) == [
        ("CAGE", "chr1", "+", 1000)]


def test_a_next_tier_call_beyond_50_bp_is_added(tmp_path):
    assert merged(tmp_path, ("CAGE", [("chr1", "+", 1000)]), ("PROcap", [("chr1", "+", 1051)])) == [
        ("CAGE", "chr1", "+", 1000), ("PROcap", "chr1", "+", 1051)]


def test_pairing_is_one_to_one(tmp_path):
    """Two PRO-cap calls near one CAGE call: one pairs with it, the other is a site of its own.
    The sweep pairs the leftmost free call, so 980 pairs and 1040 is added."""
    assert merged(tmp_path, ("CAGE", [("chr1", "+", 1000)]),
                  ("PROcap", [("chr1", "+", 980), ("chr1", "+", 1040)])) == [
        ("CAGE", "chr1", "+", 1000), ("PROcap", "chr1", "+", 1040)]


def test_strands_never_pair(tmp_path):
    assert merged(tmp_path, ("CAGE", [("chr1", "+", 1000)]), ("PROcap", [("chr1", "-", 1000)])) == [
        ("CAGE", "chr1", "+", 1000), ("PROcap", "chr1", "-", 1000)]


def test_third_tier_pairs_against_calls_the_second_added(tmp_path):
    got = merged(tmp_path, ("CAGE", [("chr1", "+", 1000)]), ("PROcap", [("chr1", "+", 2000)]),
                 ("RAMPAGE", [("chr1", "+", 2010), ("chr1", "+", 1020), ("chr2", "+", 5)]))
    assert got == [("CAGE", "chr1", "+", 1000), ("PROcap", "chr1", "+", 2000), ("RAMPAGE", "chr2", "+", 5)]


def test_higher_tier_calls_are_all_kept_even_when_close(tmp_path):
    assert merged(tmp_path, ("CAGE", [("chr1", "+", 1000), ("chr1", "+", 1020)]),
                  ("PROcap", [("chr1", "+", 1010)])) == [
        ("CAGE", "chr1", "+", 1000), ("CAGE", "chr1", "+", 1020)]


def test_read_calls_rejects_a_wider_interval(tmp_path):
    path = tmp_path / "wide.bed"
    path.write_text("chr1\t1000\t1050\tc\t1\t+\n")
    with pytest.raises(ValueError, match="1-bp"):
        union.read_calls(str(path))


def test_cli_writes_a_truth_bed_the_judge_reads(tmp_path):
    cage = truth(tmp_path, "CAGE", [("chr1", "+", 1000)])
    procap = truth(tmp_path, "PROcap", [("chr1", "+", 1030), ("chr1", "-", 400)])
    out = tmp_path / "A549.union5.truth.bed"
    r = subprocess.run([sys.executable, BUILDER, "--out", str(out), f"CAGE={cage}", f"PROcap={procap}"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    lines = out.read_text().splitlines()
    assert lines[0].startswith("#chrom\tstart\tend\tcall_id\tscore\tstrand\tsource_set")
    assert lines[1:] == ["chr1\t400\t401\tPROcap_2\t1\t-\tPROcap", "chr1\t1000\t1001\tCAGE_1\t1\t+\tCAGE"]
    sys.path.insert(0, VENDORED)
    pytest.importorskip("scipy", reason="the judge's signal_utils imports scipy")
    import ted_end_precision as judge
    peaks = judge.parse_peaks_bed(str(out))
    assert peaks == {("chr1", "+"): [(1000, 1001)], ("chr1", "-"): [(400, 401)]}


def test_cli_needs_two_sets(tmp_path):
    cage = truth(tmp_path, "CAGE", [("chr1", "+", 1000)])
    r = subprocess.run([sys.executable, BUILDER, "--out", str(tmp_path / "u.bed"), f"CAGE={cage}"],
                       capture_output=True, text=True)
    assert r.returncode != 0 and "two" in r.stderr
