"""Tests for bin/flair_eval_lib/ends.py (CONTRACTS A28, A29, A31) and the round-3 CLI flags of
bin/structure_metrics.py. Coordinates are BED-style (0-based, half-open)."""

import json
import os
import subprocess
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)

from flair_eval_lib.annotation import InputError, Transcript  # noqa: E402
from flair_eval_lib.classify import Call  # noqa: E402
from flair_eval_lib.ends import (  # noqa: E402
    DIST_EDGES, DOWN_A_EDGES, ENDS_COLUMNS, WINDOWS, SiteIndex, dist_hist, end_positions,
    end_rows, ends_section, find_motif, is_main_chrom, priming_section, write_ends_table)
from flair_eval_lib.fasta import Fasta  # noqa: E402
from flair_eval_lib.schema import load_schema, validate  # noqa: E402


def tx(exons, strand="+", chrom="chr1", name="q", gene="g"):
    return Transcript(name, gene, chrom, strand, tuple(exons))


def write(path, lines):
    path.write_text("".join(line + "\n" for line in lines))
    return str(path)


@pytest.fixture
def sites(tmp_path):
    # + strand: interval [1000,1010) and point 2000; - strand: point 5000; '.' at 9000
    # the duplicate 2000 row (other gene name, as in PolyA_DB) counts once
    return SiteIndex(write(tmp_path / "s.bed", [
        "chr1\t1000\t1010\tt1\t0\t+",
        "chr1\t2000\t2001\tp1|GENE1\t0\t+",
        "chr1\t2000\t2001\tp1|GENE2\t0\t+",
        "chr1\t5000\t5001\tp2\t0\t-",
        "chr1\t9000\t9001\tp3\t0\t.",
        "chrUn_x\t10\t20\tdrop\t0\t+",
    ]))


def test_site_index_dedupes_and_drops_non_main(sites):
    assert sites.sites == 4
    assert sites.meta()["file"] == "s.bed" and len(sites.meta()["md5"]) == 32


@pytest.mark.parametrize("strand,pos,want", [
    ("+", 1005, 0),          # inside the interval
    ("+", 1000, 0),
    ("+", 1009, 0),
    ("+", 1010, 1),          # one past the last base: downstream of the site -> positive
    ("+", 990, -10),         # upstream of the site on + -> negative
    ("+", 1600, -400),       # 2000 (-400) is nearer than the interval's last base (+591)
    ("+", 1999, -1),
    ("+", 2003, 3),
    ("-", 5010, -10),        # on -, a higher coordinate is upstream -> negative
    ("-", 4990, 10),
    ("+", 9000, 0),          # '.' site is on both strands
    ("-", 9000, 0),
    ("-", 30000, None),      # no site within 10 kb
])
def test_signed_distance(sites, strand, pos, want):
    assert sites.distance("chr1", strand, pos) == want


def test_distance_other_chrom_and_cap(sites):
    assert sites.distance("chr2", "+", 100) is None
    assert sites.distance("chr1", "+", 12000) == 3000      # from the '.' site at 9000
    assert sites.distance("chr1", "+", 19001) is None       # 10,001 bp from 9000


def test_bad_site_bed(tmp_path):
    with pytest.raises(InputError, match="BED6"):
        SiteIndex(write(tmp_path / "b.bed", ["chr1\t1\t2"]))
    with pytest.raises(InputError, match="no sites"):
        SiteIndex(write(tmp_path / "e.bed", ["chrM\t1\t2\tx\t0\t+"]))


def test_end_positions_are_strand_aware():
    t = tx([(100, 200), (300, 400)], "+")
    assert end_positions(t) == (100, 399)
    assert end_positions(tx([(100, 200), (300, 400)], "-")) == (399, 100)


def test_main_chrom():
    assert is_main_chrom("chr1") and is_main_chrom("chrX") and is_main_chrom("22")
    assert not is_main_chrom("chrM") and not is_main_chrom("chr1_KI270706v1_random")
    assert not is_main_chrom("SIRV1")


def test_find_motif_uses_priority_order():
    assert find_motif("CCATTAAACCAATAAACC") == "AATAAA"     # AATAAA ranks first
    assert find_motif("GGGAAGAAAGG") == "AAGAAA"
    assert find_motif("CCCCCC") == ""


def test_dist_hist_bins():
    h = dist_hist([-5000, -1000, -201, -200, -11, -10, -1, 0, 9, 10, 999, 1000, 5000])
    assert h["edges"] == DIST_EDGES and len(h["counts"]) == len(DIST_EDGES) + 1
    #               <-1000 [-1000,-200) [-200,-50) [-50,-10) [-10,0) [0,10) [10,50) [50,200) [200,1000) >=1000
    assert h["counts"] == [1, 2, 1, 1, 2, 2, 1, 0, 1, 2]


# genome for end_rows: chr1 is 'C' everywhere except the planted motif and A-runs below.
GENOME_LEN = 12000


def genome(tmp_path):
    seq = ["C"] * GENOME_LEN
    # + isoform A ends at 1999 (end 2000): AATAAA at 1970-1976, then 20 x A downstream
    seq[1970:1976] = "AATAAA"
    seq[2000:2020] = "A" * 20
    # - isoform B starts at 5000 (3' end on -): downstream is [4980,5000) on +; T on + = A on -
    seq[4980:4992] = "T" * 12                   # 12/20 = 0.6 exactly
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "\n".join("".join(seq[i:i + 60]) for i in range(0, GENOME_LEN, 60))
                  + "\n>chrM\nACGT\n")
    return str(fa)


def test_end_rows_and_sections(tmp_path, sites):
    ref = {"R1": tx([(1000, 1200), (1500, 2010)], "+", name="R1"),
           "R2": tx([(5000, 5300), (5600, 6000)], "-", name="R2")}
    isos = [
        tx([(1003, 1200), (1500, 2000)], "+", name="A"),     # FSM vs R1
        tx([(5000, 5300), (5600, 5900)], "-", name="B"),     # FSM vs R2 (5' end 100 bp short)
        tx([(7000, 7100)], ".", name="C"),                   # unstranded
        tx([(10, 50)], "+", chrom="chrM", name="M"),         # not a main chromosome
    ]
    calls = [Call("FSM", "CIM", "R1", "g"), Call("FSM", "FSM5", "R2", "g"),
             Call("intergenic"), Call("intergenic")]
    fasta = Fasta(genome(tmp_path))
    try:
        rows = end_rows(isos, calls, ref, sites, sites, fasta)
    finally:
        fasta.close()
    assert [r.isoform_id for r in rows] == ["A", "B", "C"]
    a, b, c = rows
    assert (a.tss, a.tes) == (1003, 1999)
    assert a.d_reftss == 0                  # inside [1000,1010)
    assert a.d_pas == -1                    # 3' end one base upstream of site 2000
    assert (a.d_ref_tss, a.d_ref_tes) == (3, -10)
    assert a.polya_motif == "AATAAA" and a.down_a_frac == 1.0
    assert (b.tss, b.tes) == (5899, 5000)
    assert (b.d_ref_tss, b.d_ref_tes) == (100, 0)    # 5' end 100 bp downstream of R2's
    assert b.d_pas == 0
    assert b.polya_motif == "" and b.down_a_frac == pytest.approx(0.6)
    assert c.d_reftss is None and c.polya_motif is None and c.down_a_frac is None

    sec = ends_section(rows, sites, None, True)
    assert sec["refs"]["pas"] is None and sec["tes"] is None
    assert sec["windows"] == WINDOWS
    assert sec["tss"]["n"] == 2                      # the unstranded isoform is left out
    assert len(sec["tss"]["within"]) == len(WINDOWS)
    assert sec["by_category"]["FSM"] == {"n": 2, "tss_within_50": 1, "tes_within_50": None}
    assert sec["fsm_tss_hist"]["counts"][DIST_EDGES.index(50) + 1] == 1     # +100 in [50,200)
    assert sum(sec["fsm_tes_hist"]["counts"]) == 2
    assert sec["polya_motif"] == {"n_with_motif": 1, "n": 2}

    pri = priming_section(rows)
    assert pri["n"] == 2 and pri["risk"] == 2        # 1.0 and exactly 0.6
    assert pri["down_a_hist"]["edges"] == DOWN_A_EDGES
    assert len(pri["down_a_hist"]["counts"]) == len(DOWN_A_EDGES)
    assert pri["down_a_hist"]["counts"][6] == 1 and pri["down_a_hist"]["counts"][10] == 1

    write_ends_table(str(tmp_path / "out"), rows)
    lines = (tmp_path / "out" / "tables" / "ends.tsv").read_text().splitlines()
    assert lines[0].split("\t") == ENDS_COLUMNS
    assert lines[1].split("\t") == ["A", "FSM", "CIM", "+", "1003", "1999", "0", "-1", "3",
                                    "-10", "AATAAA", "1"]
    assert lines[3].split("\t")[6:] == ["", "", "", "", "", ""]


def test_genome_edges_are_clipped(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\nAAAAACCCCC\n")
    f = Fasta(str(fa))
    try:
        rows = end_rows([tx([(0, 10)], "+"), tx([(0, 10)], "-", name="m")],
                        [Call("intergenic")] * 2, {}, fasta=f)
    finally:
        f.close()
    assert rows[0].n_down == 0 and rows[0].down_a_frac is None       # nothing past the end
    assert rows[1].n_down == 0


# --------------------------------------------------------------- CLI
from test_structure_metrics import _run_cli  # noqa: E402


def _cli(tmp_path, extra):
    return _run_cli(tmp_path, extra)


def test_cli_without_new_inputs_writes_only_variation(tmp_path):
    res, out = _cli(tmp_path, [])
    assert res.returncode == 0, res.stderr
    b = json.loads((out / "bundle.json").read_text())
    assert set(b["sections"]) == {"variation"}
    lines = (out / "tables" / "ends.tsv").read_text().splitlines()
    assert lines[0].split("\t") == ENDS_COLUMNS and len(lines) == 1 + 6     # chrM dropped


def test_cli_with_refs_genome_and_tusco(tmp_path):
    tss = write(tmp_path / "tss.bed", ["chr1\t95\t105\tt\t0\t+", "chr1\t5590\t5600\tt\t0\t-"])
    pas = write(tmp_path / "pas.bed", ["chr1\t599\t600\tp\t0\t+"])
    tus = write(tmp_path / "tusco.tsv", ["# comment", "# All genes: 2",
                                         "A\tA1\tGENEA\t1\tNM_1\tNP_1",
                                         "Z\tZ1\tGENEZ\t2\tNM_2\tNP_2"])
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "C" * 30000 + "\n")
    res, out = _cli(tmp_path, ["--ref-tss", tss, "--ref-pas", pas, "--tusco", tus,
                               "--genome", str(fa)])
    assert res.returncode == 0, res.stderr
    b = json.loads((out / "bundle.json").read_text())
    assert validate(b, load_schema()) == []
    sec = b["sections"]
    assert set(sec) == {"ends", "tusco", "priming", "variation"}
    assert sec["ends"]["refs"]["tss"]["sites"] == 2 and sec["ends"]["refs"]["pas"]["sites"] == 1
    assert sec["ends"]["tss"]["n"] == 6 and sec["ends"]["tes"]["n"] == 6
    assert sec["ends"]["polya_motif"]["n"] == 6
    assert sec["priming"]["n"] == 6 and sec["priming"]["risk"] == 0
    assert sec["tusco"]["genes_total"] == 2 and sec["tusco"]["per_gene"] == {"A": "TP", "Z": "FN"}
    assert sec["variation"]["genes_multi"] == 0
    assert b["inputs"]["ref_tss"] == tss


def test_cli_bad_tusco_fails(tmp_path):
    bad = write(tmp_path / "t.tsv", ["# only comments"])
    res, _ = _cli(tmp_path, ["--tusco", bad])
    assert res.returncode == 2 and "no TUSCO genes" in res.stderr


def test_cli_help_lists_new_flags():
    res = subprocess.run([sys.executable, os.path.join(BIN, "structure_metrics.py"), "--help"],
                         capture_output=True, text=True)
    assert res.returncode == 0
    for flag in ("--ref-tss", "--ref-pas", "--tusco"):
        assert flag in res.stdout
