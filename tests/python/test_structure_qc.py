"""Tests for sections.structure_qc (CONTRACTS section 13): SQANTI3 subcategories in classify.py,
bin/flair_eval_lib/structure_qc.py, the coding/non-coding priming split in sections.py, and the
structure_metrics.py --junction-file flag."""

import json
import os
import random
import subprocess
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)

from flair_eval_lib import CATEGORIES  # noqa: E402
from flair_eval_lib.annotation import InputError, Transcript, read_gtf_transcripts  # noqa: E402
from flair_eval_lib.classify import Classifier, Reference  # noqa: E402
from flair_eval_lib.ends import end_rows  # noqa: E402
from flair_eval_lib.fasta import Fasta  # noqa: E402
from flair_eval_lib.schema import load_schema, validate  # noqa: E402
from flair_eval_lib.sections import add_sections  # noqa: E402
from flair_eval_lib.structure_qc import (  # noqa: E402
    SUBCATEGORIES, coding_priming, has_repeat, read_junction_file, rts_windows,
    structure_qc_section)

# Reference, BED-style exons.
#   gene A (+): A1 = 100-200, 300-400, 500-600       introns (200,300) (400,500)
#               A2 = 100-200, 500-600, 700-800       introns (200,500) (600,700)
#   gene D (+): D1 = 1100-1200, 1300-1400, 1500-1600, 1700-1800
#   gene C (-): C1 = 5000-5100, 5300-5400, 5500-5600
REF = [
    ("A1", "A", "chr1", "+", [(100, 200), (300, 400), (500, 600)]),
    ("A2", "A", "chr1", "+", [(100, 200), (500, 600), (700, 800)]),
    ("D1", "D", "chr1", "+", [(1100, 1200), (1300, 1400), (1500, 1600), (1700, 1800)]),
    ("C1", "C", "chr1", "-", [(5000, 5100), (5300, 5400), (5500, 5600)]),
]


def write_gtf(path, models):
    with open(path, "w") as fh:
        for tid, gene, chrom, strand, exons in models:
            attrs = f'gene_id "{gene}"; transcript_id "{tid}";'
            fh.write(f"{chrom}\tsyn\ttranscript\t{exons[0][0] + 1}\t{exons[-1][1]}\t.\t{strand}\t.\t{attrs}\n")
            for s, e in exons:
                fh.write(f"{chrom}\tsyn\texon\t{s + 1}\t{e}\t.\t{strand}\t.\t{attrs}\n")


def write_bed(path, models):
    with open(path, "w") as fh:
        for name, chrom, strand, exons in models:
            s, e = exons[0][0], exons[-1][1]
            sizes = ",".join(str(b - a) for a, b in exons)
            starts = ",".join(str(a - s) for a, _ in exons)
            fh.write(f"{chrom}\t{s}\t{e}\t{name}\t1000\t{strand}\t{s}\t{s}\t0,0,0\t"
                     f"{len(exons)}\t{sizes},\t{starts},\n")


@pytest.fixture(scope="module")
def ref_tx(tmp_path_factory):
    gtf = tmp_path_factory.mktemp("ref") / "ref.gtf"
    write_gtf(gtf, REF)
    return read_gtf_transcripts(str(gtf))


@pytest.fixture(scope="module")
def clf(ref_tx):
    return Classifier(Reference(ref_tx))


def iso(exons, strand="+", name="q", gene="G"):
    return Transcript(name, gene, "chr1", strand, tuple(exons))


# --------------------------------------------------------------------- subcategories (A42)
DETAILS = [
    # name, exons, strand, category, detail
    ("ism_3prime_plus", [(300, 400), (500, 600)], "+", "ISM", "3prime_fragment"),
    ("ism_5prime_plus", [(100, 200), (300, 400)], "+", "ISM", "5prime_fragment"),
    ("ism_internal", [(1300, 1400), (1500, 1600)], "+", "ISM", "internal_fragment"),
    ("ism_3prime_d", [(1300, 1400), (1500, 1600), (1700, 1800)], "+", "ISM", "3prime_fragment"),
    # minus strand: the genomic-left junction is the 3' end
    ("ism_3prime_minus", [(5000, 5100), (5300, 5400)], "-", "ISM", "3prime_fragment"),
    ("ism_5prime_minus", [(5300, 5400), (5500, 5600)], "-", "ISM", "5prime_fragment"),
    ("ism_intron_retention", [(1100, 1200), (1300, 1600)], "+", "ISM", "intron_retention"),
    ("ism_mono", [(310, 390)], "+", "ISM", "mono-exon"),
    ("nic_known_junctions", [(100, 200), (300, 400), (500, 600), (700, 800)], "+",
     "NIC", "combination_of_known_junctions"),
    ("nic_known_sites", [(100, 200), (300, 400), (700, 800)], "+",
     "NIC", "combination_of_known_splicesites"),
    ("nic_intron_retention", [(100, 200), (300, 600), (700, 800)], "+", "NIC", "intron_retention"),
    ("nic_mono_ir", [(150, 350)], "+", "NIC", "mono-exon_by_intron_retention"),
    ("nnc_novel_site", [(100, 200), (320, 400), (500, 600)], "+",
     "NNC", "at_least_one_novel_splicesite"),
    ("nnc_intron_retention", [(100, 200), (320, 600)], "+", "NNC", "intron_retention"),
]


@pytest.mark.parametrize("name,exons,strand,cat,detail", DETAILS, ids=[d[0] for d in DETAILS])
def test_subcategory(clf, name, exons, strand, cat, detail):
    c = clf.classify(iso(exons, strand))
    assert (c.category, c.detail) == (cat, detail)


def test_every_subcategory_is_covered():
    got = {(d[3], d[4]) for d in DETAILS}
    want = {(c, s) for c, subs in SUBCATEGORIES.items() for s in subs}
    assert got == want


def test_fsm_and_other_classes_have_no_detail(clf):
    assert clf.classify(iso([(100, 200), (300, 400), (500, 600)])).detail == ""
    assert clf.classify(iso([(40000, 40100), (40300, 40400)])).detail == ""


# --------------------------------------------------------------------- RT switching
def test_has_repeat():
    assert has_repeat("CCAAAAGGGGTT", "TTTAAAAGGGGA")          # shared 8-mer AAAAGGGG
    assert has_repeat("CCAAAAGGGGTT", "TTTAAACGGGGA")          # one mismatch
    assert not has_repeat("CCAAAAGGGGTT", "TTTAACCGGGGA")      # two mismatches
    assert not has_repeat("AAAAAAA", "AAAAAAA")                # shorter than 8
    assert not has_repeat("ACGTACGTACGT", "TTTTTTTTTTTT")


def _genome(tmp_path, seq):
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    return Fasta(str(fa))


def test_rts_windows_strands(tmp_path):
    seq = "".join("ACGT"[(i * 7 + i // 5) % 4] for i in range(200))
    f = _genome(tmp_path, seq)
    rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    # + strand: 12 nt ending one base into the intron, and ending one base into the next exon
    assert rts_windows(f, "chr1", "+", 50, 120) == (seq[39:51], seq[109:121])
    # - strand: the same, read on the minus strand (donor at the right-hand end)
    assert rts_windows(f, "chr1", "-", 50, 120) == (rc(seq[119:131]), rc(seq[49:61]))
    f.close()


def _planted_genome():
    """GT-AG intron (40, 100) whose donor and acceptor windows share AAAAGGGG; GT-AG intron
    (120, 180) with no repeat."""
    rnd = random.Random(7)
    seq = [rnd.choice("ACGT") for _ in range(260)]
    seq[29:37] = "AAAAGGGG"           # in the donor window [29, 41)
    seq[89:97] = "AAAAGGGG"           # in the acceptor window [89, 101)
    seq[40:42], seq[98:100] = "GT", "AG"
    seq[120:122], seq[178:180] = "GT", "AG"
    return "".join(seq)


# --------------------------------------------------------------------- junctions file
def test_read_junction_file_star_and_bed(tmp_path):
    star = tmp_path / "SJ.out.tab"
    star.write_text("chr1\t201\t300\t1\t1\t1\t5\t2\t30\nchr1\t401\t500\t1\t1\t1\t0\t4\t30\n")
    assert read_junction_file(str(star)) == {("chr1", 200, 300): 5, ("chr1", 400, 500): 0}
    bed = tmp_path / "j.bed"
    bed.write_text("track name=x\nchr1\t200\t300\tsj1\t3\t+\nchr1\t200\t300\tsj1\t4\t-\n")
    assert read_junction_file(str(bed)) == {("chr1", 200, 300): 7}     # strand ignored


@pytest.mark.parametrize("text", ["", "chr1\t1\t2\n", "chr1\tx\t300\tj\t3\t+\n"])
def test_read_junction_file_rejects_bad_input(tmp_path, text):
    p = tmp_path / "bad.bed"
    p.write_text(text)
    with pytest.raises(InputError):
        read_junction_file(str(p))


# --------------------------------------------------------------------- the section
def test_section_counts(tmp_path, ref_tx, clf):
    isoforms = [
        iso([(100, 200), (300, 400), (500, 600)], name="fsm_a"),
        iso([(120, 200), (300, 400), (500, 590)], name="fsm_b"),           # same reference
        iso([(300, 400), (500, 600)], name="ism"),
        iso([(100, 200), (320, 400), (500, 600)], name="nnc"),
        iso([(310, 390)], name="mono"),
    ]
    calls = [clf.classify(t) for t in isoforms]
    jf = {("chr1", 200, 300): 4, ("chr1", 400, 500): 0}
    q = structure_qc_section(isoforms, calls, Reference(ref_tx), ref_tx, None, None, jf,
                             {"file": "SJ.out.tab", "md5": None})
    assert set(q) == {"subcategories", "redundancy", "by_category", "ref_recovery", "junction_support"}
    assert q["subcategories"]["ISM"]["3prime_fragment"] == 1
    assert q["subcategories"]["ISM"]["mono-exon"] == 1
    assert q["subcategories"]["NNC"]["at_least_one_novel_splicesite"] == 1
    r = q["redundancy"]
    assert (r["fsm_isoforms"], r["fsm_refs"], r["fsm_refs_multi"]) == (2, 1, 1)
    assert (r["fsm_ism_isoforms"], r["fsm_ism_refs"], r["fsm_ism_refs_multi"]) == (4, 1, 1)
    assert r["fsm_hist"]["counts"][:2] == [0, 1]
    for k in ("length", "exons"):
        counts = q["by_category"][k]["counts"]
        assert set(counts) == set(CATEGORIES)
        assert sum(counts["FSM"]) == 2 and sum(counts["ISM"]) == 2 and sum(counts["NNC"]) == 1
    rr = q["ref_recovery"]["exons"]
    assert sum(rr["all"]) == len(ref_tx) and sum(rr["fsm"]) == 1
    j = q["junction_support"]
    # distinct junctions: (200,300) (400,500) (200,320); only (200,300) has reads
    assert (j["junctions"], j["supported"], j["file_junctions"]) == (3, 1, 2)
    assert j["known"] == {"n": 2, "supported": 1} and j["novel"] == {"n": 1, "supported": 0}
    assert (j["multi_exon"], j["all_supported"]) == (4, 0)
    assert j["by_category"]["FSM"] == {"n": 2, "all_supported": 0}


def test_section_with_genome(tmp_path):
    seq = _planted_genome()
    f = _genome(tmp_path, seq)
    ref = [Transcript("R1", "R", "chr1", "+", ((10, 40), (100, 120), (180, 240)))]
    isoforms = [Transcript("i1", "R", "chr1", "+", ((10, 40), (100, 120), (180, 240))),
                Transcript("i2", "R", "chr1", "+", ((10, 40), (100, 150))),
                Transcript("i3", "R", "chr1", "+", ((60, 80),))]
    reference = Reference(ref)
    clf = Classifier(reference)
    calls = [clf.classify(t) for t in isoforms]
    rows = end_rows(isoforms, calls, {t.tid: t for t in ref}, None, None, f)
    q = structure_qc_section(isoforms, calls, reference, ref, f, rows)
    assert q["canonical"] == {"multi_exon": 2, "all_canonical": 2,
                              "by_category": {"FSM": {"n": 1, "all_canonical": 1},
                                              "ISM": {"n": 1, "all_canonical": 1}}}
    t = q["rts"]
    assert (t["junctions"], t["rts"], t["isoforms"], t["isoforms_with_rts"]) == (2, 1, 2, 2)
    assert t["by_type"]["known_canonical"] == {"n": 2, "rts": 1}
    assert t["params"] == {"pattern": 10, "wiggle": 1, "min_match": 8, "mismatches": 1}
    p = q["priming"]
    assert set(p) == {"by_category", "by_exons"}
    assert p["by_exons"]["mono"]["n"] == 1 and p["by_exons"]["multi"]["n"] == 2
    assert sum(v["n"] for v in p["by_category"].values()) == 3
    f.close()


# --------------------------------------------------------------------- coding vs non-coding
ENDS_HEAD = ("isoform_id category subcategory strand tss tes d_reftss d_pas d_ref_tss "
             "d_ref_tes polya_motif down_a_frac").split()


def _write_ends(path, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(ENDS_HEAD) + "\n")
        for iid, frac in rows:
            fh.write("\t".join([iid, "FSM", "CIM", "+", "1", "9", "", "", "", "", "", frac]) + "\n")


def test_coding_priming(tmp_path):
    ends, prod = tmp_path / "ends.tsv", tmp_path / "productivity.tsv"
    _write_ends(ends, [("a", "0.6"), ("b", "0.1"), ("c", "0.9"), ("d", ""), ("e", "1")])
    prod.write_text("isoform_id\tclass\na\tPRO\nb\tPTC\nc\tNGO\nd\tPRO\n")
    assert coding_priming(str(ends), str(prod)) == {"coding": {"n": 2, "risk": 1},
                                                    "noncoding": {"n": 1, "risk": 1}}
    assert coding_priming(str(ends), str(tmp_path / "missing.tsv")) is None


def test_add_sections_adds_by_coding(tmp_path):
    os.makedirs(tmp_path / "tables")
    _write_ends(tmp_path / "tables" / "ends.tsv", [("a", "0.7")])
    (tmp_path / "productivity.tsv").write_text("isoform_id\tclass\na\tNST\n")
    pri = {"by_category": {}, "by_exons": {"mono": {"n": 0, "risk": 0}, "multi": {"n": 1, "risk": 1}}}
    (tmp_path / "bundle.json").write_text(json.dumps({"sections": {"structure_qc": {"priming": pri}}}))
    b = add_sections(str(tmp_path))
    assert b["sections"]["structure_qc"]["priming"]["by_coding"]["coding"] == {"n": 1, "risk": 1}
    assert b["sections"]["productivity"] == {"class_counts": {"NST": 1}}


# --------------------------------------------------------------------- CLI
def test_cli_junction_file(tmp_path):
    seq = _planted_genome()
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + seq + "\n")
    gtf = tmp_path / "ref.gtf"
    write_gtf(gtf, [("R1", "R", "chr1", "+", [(10, 40), (100, 120), (180, 240)])])
    bed = tmp_path / "iso.bed"
    models = [("i1", "chr1", "+", [(10, 40), (100, 120), (180, 240)]),
              ("i2", "chr1", "+", [(10, 40), (100, 150)])]
    write_bed(bed, models)
    igtf = tmp_path / "iso.gtf"
    write_gtf(igtf, [(n, "G1", c, s, e) for n, c, s, e in models])
    sj = tmp_path / "SJ.out.tab"
    sj.write_text("chr1\t41\t100\t1\t1\t1\t9\t0\t40\n")
    rj, ej = tmp_path / "run.json", tmp_path / "env.json"
    rj.write_text(json.dumps({"id": "r1"}))
    ej.write_text(json.dumps({"source": "release", "flair_version": "3.0.1"}))
    out = tmp_path / "b"
    res = subprocess.run([sys.executable, os.path.join(BIN, "structure_metrics.py"),
                          "--isoforms-bed", str(bed), "--isoforms-gtf", str(igtf), "--gtf", str(gtf),
                          "--sample", "s", "--run-json", str(rj), "--env-json", str(ej),
                          "--genome", str(fa), "--junction-file", str(sj), "--out", str(out)],
                         capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    b = json.loads((out / "bundle.json").read_text())
    assert validate(b, load_schema()) == []
    q = b["sections"]["structure_qc"]
    assert {"subcategories", "redundancy", "by_category", "ref_recovery", "canonical", "rts",
            "priming", "junction_support"} == set(q)
    j = q["junction_support"]
    assert j["file"] == "SJ.out.tab" and len(j["md5"]) == 32
    assert (j["junctions"], j["supported"], j["all_supported"], j["multi_exon"]) == (2, 1, 1, 2)
    assert b["inputs"]["junction_file"] == str(sj)
    assert q["subcategories"]["ISM"]["5prime_fragment"] == 1


def test_schema_rejects_unknown_structure_qc_key():
    s = load_schema()
    sec = {"subcategories": {"ISM": {k: 0 for k in SUBCATEGORIES["ISM"]},
                             "NIC": {k: 0 for k in SUBCATEGORIES["NIC"]},
                             "NNC": {k: 0 for k in SUBCATEGORIES["NNC"]}}}
    qc_schema = s["properties"]["sections"]["properties"]["structure_qc"]
    assert any("redundancy" in e for e in validate(sec, qc_schema))
    sec["subcategories"]["ISM"]["bogus"] = 1
    assert any("bogus" in e for e in validate(sec, qc_schema))
