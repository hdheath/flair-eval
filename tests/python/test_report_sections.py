"""Tests for the round-2 report sections (CONTRACTS section 8): concordance, CPM agreement,
the quant/productivity section builders, and how build_report.py embeds them."""

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

from flair_eval_lib.annotation import Transcript  # noqa: E402
from flair_eval_lib.concordance import (  # noqa: E402
    RunModels, assign_structures, concordance, pearson, quant_agreement, read_counts,
    read_header_tsv, spearman)
from flair_eval_lib.sections import (  # noqa: E402
    CPM_EDGES, add_sections, productivity_section, quant_section)


def T(name, exons, strand="+", chrom="chr1"):
    return Transcript(name, "", chrom, strand, tuple(exons))


def run(rid, isoforms, info=None):
    return RunModels(rid, isoforms, info or {t.tid: ("FSM", "") for t in isoforms})


# --------------------------------------------------------------------- matching
def test_multi_exon_match_by_intron_chain_ignores_ends():
    a = run("a", [T("a1", [(100, 200), (300, 400)])])
    b = run("b", [T("b1", [(50, 200), (300, 480)]), T("b2", [(100, 210), (300, 400)])])
    sid = assign_structures([a, b])
    assert sid[("a", "a1")] == sid[("b", "b1")]
    assert sid[("a", "a1")] != sid[("b", "b2")]


def test_multi_exon_needs_same_strand_and_chrom():
    a = run("a", [T("a1", [(100, 200), (300, 400)])])
    b = run("b", [T("b1", [(100, 200), (300, 400)], "-"), T("b2", [(100, 200), (300, 400)], "+", "chr2")])
    sid = assign_structures([a, b])
    assert len({sid[("a", "a1")], sid[("b", "b1")], sid[("b", "b2")]}) == 3


@pytest.mark.parametrize("other,same", [
    ((100, 200), True),          # identical
    ((140, 240), True),          # 60% of both
    ((160, 260), False),         # 40%
    ((100, 400), False),         # 100% of a, 33% of b: not reciprocal
])
def test_mono_exon_reciprocal_half_overlap(other, same):
    a = run("a", [T("a1", [(100, 200)])])
    b = run("b", [T("b1", [other])])
    sid = assign_structures([a, b])
    assert (sid[("a", "a1")] == sid[("b", "b1")]) is same


def test_mono_exon_strand_matters():
    sid = assign_structures([run("a", [T("a1", [(100, 200)], "+")]), run("b", [T("b1", [(100, 200)], "-")])])
    assert sid[("a", "a1")] != sid[("b", "b1")]


# --------------------------------------------------------------------- concordance
def test_concordance_pairs_genes_and_loci():
    a = run("a", [T("a1", [(100, 200), (300, 400)]), T("a2", [(1000, 1100)]), T("a3", [(5000, 5100), (5200, 5300)])],
            {"a1": ("FSM", "G1"), "a2": ("FSM", "G2"), "a3": ("NNC", "G3")})
    b = run("b", [T("b1", [(90, 200), (300, 410)]), T("b2", [(1020, 1120)]), T("b3", [(7000, 7100), (7200, 7300)])],
            {"b1": ("FSM", "G1"), "b2": ("ISM", "G2"), "b3": ("NIC", "G4")})
    c = concordance([a, b], "a", {"G1": "GENE1", "G4": "GENE4"})
    assert c["structures"] == {"a": 3, "b": 3}
    assert c["union"] == 4 and c["shared_by_all"] == 2
    p = c["pairs"][0]
    assert (p["shared"], p["only_a"], p["only_b"], p["jaccard"]) == (2, 1, 1, 0.5)
    g = c["genes_vs_baseline"]["b"]
    assert (g["gained"], g["lost"]) == (1, 1)
    assert g["gained_names"] == ["GENE4"] and g["lost_names"] == ["G3"]
    loci = c["top_loci"]
    assert {(l["start"], l["end"]) for l in loci} == {(5000, 5300), (7000, 7300)}
    for l in loci:
        assert l["structures"] == 1 and l["shared_by_all"] == 0 and sorted(l["counts"]) == [0, 1]
    assert {tuple(l["genes"]) for l in loci} == {("G3",), ("GENE4",)}


def test_concordance_needs_two_runs():
    assert concordance([run("a", [T("a1", [(1, 5)])]), None], "a", {}) is None


def test_top_loci_ranked_by_divergence():
    shared = [T("x", [(100, 200), (300, 400)])]
    a = run("a", shared + [T("d1", [(5000, 5100), (5200, 5300)]), T("d2", [(5000, 5100), (5250, 5300)]),
                           T("e1", [(9000, 9100), (9200, 9300)])])
    b = run("b", [T("y", [(100, 200), (300, 400)])])
    loci = concordance([a, b], "a", {})["top_loci"]
    assert [(l["start"], l["structures"]) for l in loci] == [(5000, 2), (9000, 1)]


# --------------------------------------------------------------------- expression
def test_correlations():
    assert pearson([1, 2, 3], [2, 4, 6]) == pytest.approx(1.0)
    assert spearman([1, 2, 3, 4], [10, 20, 30, 1000]) == pytest.approx(1.0)
    assert pearson([1, 1, 1], [1, 2, 3]) is None
    assert pearson([1, 2], [1, 2]) is None


def test_quant_agreement_matches_by_structure_and_sums():
    a = run("a", [T("a1", [(100, 200), (300, 400)]), T("a2", [(110, 200), (300, 400)]),
                  T("a3", [(1000, 1100), (1200, 1300)]), T("a4", [(2000, 2100), (2200, 2300)])])
    b = run("b", [T("b1", [(100, 200), (300, 400)]), T("b3", [(1000, 1100), (1200, 1300)]),
                  T("b4", [(2000, 2100), (2200, 2300)]), T("b5", [(8000, 8100), (8200, 8300)])])
    cpm = {"a": {"a1": 50, "a2": 50, "a3": 10, "a4": 1000}, "b": {"b1": 99, "b3": 11, "b4": 900, "b5": 40}}
    q = quant_agreement([a, b], cpm, "a")
    p = q["pairs"][0]
    assert p["n"] == 3 and p["pearson"] > 0.99
    assert p["cpm_shared_b"] == pytest.approx(1010 / 1050, abs=1e-3)
    assert sum(map(sum, p["grid"])) == 3
    assert q["baseline"] == "a"


def test_read_header_tsv_and_counts(tmp_path):
    f = tmp_path / "counts.tsv"
    f.write_text("isoform_id\tcount\tcpm\ni1\t3\t300.5\ni2\t1\t100\n")
    assert read_counts(str(f)) == {"i1": 300.5, "i2": 100.0}
    f.write_text("i1\t3\t300.5\n")
    with pytest.raises(ValueError, match="header"):
        read_header_tsv(str(f), ["isoform_id", "cpm"])


def test_section_builders(tmp_path):
    c = tmp_path / "counts.tsv"
    c.write_text("isoform_id\tcount\tcpm\ni1\t10\t0.5\ni2\t0\t0\ni3\t5\t50\ni4\t99\t40000\nother\t7\t9\n")
    q = quant_section(str(c), {"i1", "i2", "i3", "i4"})
    assert q["isoforms_quantified"] == 3 and q["reads_assigned"] == 114
    assert q["cpm_hist"]["edges"] == CPM_EDGES and len(q["cpm_hist"]["counts"]) == len(CPM_EDGES)
    assert q["cpm_hist"]["counts"][0] == 1 and q["cpm_hist"]["counts"][-1] == 1
    p = tmp_path / "productivity.tsv"
    p.write_text("isoform_id\tcds_start\tcds_end\tclass\ni1\t1\t9\tPRO\ni2\t\t\tNGO\ni3\t1\t9\tPRO\nx\t1\t2\tPTC\n")
    assert productivity_section(str(p), {"i1", "i2", "i3"}) == {"class_counts": {"NGO": 1, "PRO": 2}}


def test_add_sections_on_a_fixture_copy(tmp_path):
    src = os.path.join(FIXTURES, "wtc11_pacbio", "flair-3.0.1")
    d = tmp_path / "b"
    shutil.copytree(src, d)
    b = json.load(open(d / "bundle.json"))
    b["sections"] = {}
    json.dump(b, open(d / "bundle.json", "w"))
    out = add_sections(str(d))
    assert out["sections"]["quant"]["reads_assigned"] > 0
    assert out["sections"]["productivity"]["class_counts"]["PRO"] > 0


# --------------------------------------------------------------------- report wiring
def build(*args):
    return subprocess.run([sys.executable, os.path.join(BIN, "build_report.py"), *map(str, args)],
                          capture_output=True, text=True)


def embedded(path):
    html = open(path, encoding="utf-8").read()
    return json.loads(re.search(r'id="report-data">(.*?)</script>', html, re.S).group(1))


def test_fixtures_carry_round2_inputs():
    for sample in os.listdir(FIXTURES):
        sd = os.path.join(FIXTURES, sample)
        if not os.path.isdir(sd):
            continue
        for rid in os.listdir(sd):
            d = os.path.join(sd, rid)
            for f in ("isoforms.bed", "counts.tsv", "productivity.tsv", "tables/genes.tsv"):
                assert os.path.exists(os.path.join(d, f)), (d, f)
            b = json.load(open(os.path.join(d, "bundle.json")))
            assert set(b["sections"]) >= {"quant", "productivity"}
            assert b["compute"]["cpu_hours"] is not None
            assert open(os.path.join(d, "counts.tsv")).readline().split() == ["isoform_id", "count", "cpm"]


def test_report_embeds_concordance_and_quant(tmp_path):
    out = tmp_path / "r.html"
    res = build("--bundles", FIXTURES, "--out", out, "--baseline", "flair-3.0.1")
    assert res.returncode == 0, res.stderr
    d = embedded(out)
    assert d["runs"][0] == "flair-3.0.1"            # baseline first: it keeps colour slot 1
    for sample in d["samples"]:
        c = d["by_sample"][sample]["concordance"]
        n = len(c["runs"])
        assert len(c["pairs"]) == n * (n - 1) // 2
        assert c["baseline"] == "flair-3.0.1"
        assert 0 < len(c["top_loci"]) <= 20
        assert all(len(l["counts"]) == n for l in c["top_loci"])
        assert set(c["genes_vs_baseline"]) == set(c["runs"]) - {"flair-3.0.1"}
        q = d["by_sample"][sample]["quant"]
        assert all(p["pearson"] is None or -1 <= p["pearson"] <= 1 for p in q["pairs"])
        assert sum(1 for p in q["pairs"] if "grid" in p) == n - 1


def test_default_baseline_is_name_order(tmp_path):
    out = tmp_path / "r.html"
    assert build("--bundles", FIXTURES, "--out", out).returncode == 0
    d = embedded(out)
    assert d["baseline"] == sorted(d["runs"])[0] == d["runs"][0]


def test_no_round2_inputs_means_no_sections(tmp_path):
    src = os.path.join(FIXTURES, "wtc11_ont")
    for rid in os.listdir(src):
        d = tmp_path / "b" / rid
        os.makedirs(d)
        shutil.copy(os.path.join(src, rid, "bundle.json"), d / "bundle.json")
    out = tmp_path / "r.html"
    assert build("--bundles", tmp_path / "b", "--out", out).returncode == 0
    assert embedded(out)["by_sample"] == {}


def test_size_budget_with_full_bundles(tmp_path):
    n = 0
    for copy in range(5):
        for sample in ("wtc11_pacbio", "wtc11_ont"):
            for rid in sorted(os.listdir(os.path.join(FIXTURES, sample))):
                dst = tmp_path / "many" / f"{sample}_{copy}" / rid
                shutil.copytree(os.path.join(FIXTURES, sample, rid), dst)
                b = json.load(open(dst / "bundle.json"))
                b["sample"] = f"{sample}_{copy}"
                json.dump(b, open(dst / "bundle.json", "w"))
                n += 1
    assert n == 50
    out = tmp_path / "r.html"
    res = build("--bundles", tmp_path / "many", "--out", out)
    assert res.returncode == 0, res.stderr
    assert len(embedded(out)["by_sample"]) == 10
    assert os.path.getsize(out) < 5 * 1024 * 1024


def test_structure_metrics_writes_gene_names():
    rows = read_header_tsv(os.path.join(FIXTURES, "wtc11_pacbio", "flair-3.0.1", "tables", "genes.tsv"),
                           ["gene_id", "gene_name"])
    names = {r["gene_name"] for r in rows}
    assert "TP53" in names
