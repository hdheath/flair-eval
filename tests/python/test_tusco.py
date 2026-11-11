"""Tests for bin/flair_eval_lib/tusco.py: SQANTI3 TUSCO rules (CONTRACTS A30)."""

import os
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "bin"))

from flair_eval_lib.annotation import InputError, Regions, Transcript  # noqa: E402
from flair_eval_lib.classify import Call  # noqa: E402
from flair_eval_lib.tusco import read_tusco, tusco_section  # noqa: E402

TSV = ["# Generated from x", "# All genes: 6",
       "ENSG01\tENST01\tGA\t1\tNM_01\tNP_01",
       "ENSG02\tENST02\tGB\t2\tNM_02\tNP_02",
       "ENSG03\tENST03\tGC\t3\tNM_03\tNP_03",
       "ENSG04\tENST04\tGD\t4\tNM_04\tNP_04",
       "ENSG05\tENST05\tGE\t5\tNM_05\tNP_05",
       "ENSG06\tENST06\tGF\t6\tNM_06\tNP_06",
       "ENSG06\tENST06\tGF\t6\tNM_06\tNP_06"]          # duplicate row counts once


@pytest.fixture
def tusco(tmp_path):
    p = tmp_path / "tusco.tsv"
    p.write_text("\n".join(TSV) + "\n")
    return read_tusco(str(p)), str(p)


def ref(gene, start, name=None, chrom="chr1"):
    return Transcript(f"T_{gene}", f"{gene}.7", chrom, "+", ((start, start + 100),),
                      name or "")


REF = [ref("ENSG01", 0, "GA"), ref("ENSG02", 1000, "GB"), ref("ENSG03", 2000, "GC"),
       ref("ENSG04", 3000, "GD"), ref("ENSG05", 4000, "GE"), ref("ENSG06", 5000, "GF"),
       ref("ENSG99", 6000, "OTHER")]

CALLS = [
    Call("FSM", "CIM", "T_ENSG01", "ENSG01.7"),       # GA: TP
    Call("FSM", "FSM5", "T_ENSG01", "ENSG01.7"),      # GA: FSM in a TP gene: not PTP (R:217-218)
    Call("NIC", "", "", "ENSG01.7"),                  # GA: FP transcript, gene stays TP
    Call("ISM", "", "T_ENSG02", "ENSG02.7"),          # GB: PTP
    Call("FSM", "FSM3", "T_ENSG02", "ENSG02.7"),      # GB: PTP
    Call("genic_intron", "", "", "ENSG03.7"),         # GC: TUSCO transcript, not FP-category
    Call("fusion", "", "", "ENSG04.7,ENSG99.7"),      # GD: FP via fusion split
    Call("antisense"),                                # no gene
    Call("FSM", "CIM", "T_ENSG99", "ENSG99.7"),       # not a TUSCO gene
]                                                     # GE, GF: FN


def test_read_tusco(tusco):
    rows, _ = tusco
    assert len(rows) == 6
    assert rows[0] == {"ensembl": "ENSG01", "refseq": "NM_01", "gene_name": "GA"}


def test_read_tusco_errors(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("# x\nENSG1\tENST1\n")
    with pytest.raises(InputError, match="6 columns"):
        read_tusco(str(p))


def test_sqanti3_rules(tusco):
    rows, path = tusco
    s = tusco_section(rows, CALLS, REF, None, path)
    assert s["id_type"] == "ensembl"
    assert s["genes_total"] == 6 and s["genes_in_regions"] == 6
    assert s["per_gene"] == {"ENSG01": "TP", "ENSG02": "PTP", "ENSG03": "FP", "ENSG04": "FP",
                             "ENSG05": "FN", "ENSG06": "FN"}
    assert (s["TP"], s["PTP"], s["FP"], s["FN"]) == (1, 1, 2, 2)
    # transcripts: 7 TUSCO transcripts; TP 1, PTP 2 (GB), FP 2 (NIC + fusion part)
    assert s["transcripts"] == {"n": 7, "TP": 1, "PTP": 2, "FP": 2}
    assert s["sensitivity"] == pytest.approx(1 / 6, abs=1e-6)
    assert s["precision"] == pytest.approx(1 / 7, abs=1e-6)
    assert s["precision_redundant"] == pytest.approx(3 / 7, abs=1e-6)
    assert s["positive_detection_rate"] == pytest.approx(2 / 6, abs=1e-6)
    assert s["gene_names"]["ENSG01"] == "GA"
    assert s["file"] == "tusco.tsv" and len(s["md5"]) == 32


def test_gene_name_ids_when_reference_has_no_ensembl(tusco):
    rows, _ = tusco
    refs = [Transcript("t1", "sym1", "chr1", "+", ((0, 100),), "GA")]
    s = tusco_section(rows, [Call("FSM", "CIM", "t1", "sym1")], refs)
    assert s["id_type"] == "gene_name"
    assert s["per_gene"]["ENSG01"] == "TP" and s["TP"] == 1 and s["FN"] == 5


def test_regions_limit_scope(tusco):
    rows, _ = tusco
    regions = Regions({"chr1": [(0, 1500)]})                 # GA and GB only
    s = tusco_section(rows, CALLS, REF, regions)
    assert s["genes_total"] == 6 and s["genes_in_regions"] == 2
    assert s["per_gene"] == {"ENSG01": "TP", "ENSG02": "PTP"}
    assert s["sensitivity"] == 0.5


def test_no_tusco_gene_in_regions(tusco):
    rows, _ = tusco
    s = tusco_section(rows, CALLS, REF, Regions({"chr2": [(0, 10)]}))
    assert s["genes_in_regions"] == 0 and s["per_gene"] == {}
    assert s["sensitivity"] is None and s["precision"] is None
    assert s["precision_redundant"] is None


def test_mono_exon_cim_is_tp(tusco):
    rows, _ = tusco
    s = tusco_section(rows, [Call("FSM", "CIM", "T_ENSG05", "ENSG05.7")], REF)
    assert s["per_gene"]["ENSG05"] == "TP"
