"""Unit tests for bin/structure_metrics.py and bin/flair_eval_lib (classifier + metrics).

Every SQANTI category and FSM subcategory gets a tiny synthetic case. Coordinates below are
BED-style (0-based, half-open); the helper writes them to GTF as 1-based inclusive.
"""

import json
import os
import subprocess
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)

from flair_eval_lib import CATEGORIES, HISTOGRAM_EDGES  # noqa: E402
from flair_eval_lib.annotation import (  # noqa: E402
    InputError, Regions, Transcript, assign_isoform_genes, keep_contig, read_bed12,
    read_gtf_transcripts)
from flair_eval_lib.classify import Classifier, Reference, fsm_subcategory  # noqa: E402
from flair_eval_lib.fasta import Fasta  # noqa: E402
from flair_eval_lib.metrics import histogram, normalise_version  # noqa: E402
from flair_eval_lib.schema import load_schema, validate  # noqa: E402

# Reference, BED-style exons.
#   gene A (+): A1 = 100-200, 300-400, 500-600       introns (200,300) (400,500)
#               A2 = 100-200, 500-600, 700-800       introns (200,500) (600,700)
#   gene M (+): M1 = 2000-3000 (mono-exon)
#   gene C (-): C1 = 5000-5100, 5300-5400, 5500-5600
#   gene B (+): B1 = 10000-10100, 10300-10400, 10500-10600
#   gene T (+): T1 = 20000-21000 (mono-exon only gene)
REF = [
    ("A1", "A", "chr1", "+", [(100, 200), (300, 400), (500, 600)]),
    ("A2", "A", "chr1", "+", [(100, 200), (500, 600), (700, 800)]),
    ("M1", "M", "chr1", "+", [(2000, 3000)]),
    ("C1", "C", "chr1", "-", [(5000, 5100), (5300, 5400), (5500, 5600)]),
    ("B1", "B", "chr1", "+", [(10000, 10100), (10300, 10400), (10500, 10600)]),
    ("T1", "T", "chr1", "+", [(20000, 21000)]),
    ("X1", "X", "chrM", "+", [(10, 500)]),
]


def write_gtf(path, models):
    with open(path, "w") as fh:
        fh.write("##description: synthetic\n")
        for tid, gene, chrom, strand, exons in models:
            s, e = exons[0][0], exons[-1][1]
            attrs = f'gene_id "{gene}"; transcript_id "{tid}";'
            fh.write(f"{chrom}\tsyn\ttranscript\t{s + 1}\t{e}\t.\t{strand}\t.\t{attrs}\n")
            for es, ee in exons:
                fh.write(f"{chrom}\tsyn\texon\t{es + 1}\t{ee}\t.\t{strand}\t.\t{attrs}\n")


def write_bed(path, models):
    with open(path, "w") as fh:
        for name, chrom, strand, exons in models:
            s, e = exons[0][0], exons[-1][1]
            sizes = ",".join(str(b - a) for a, b in exons)
            starts = ",".join(str(a - s) for a, _ in exons)
            fh.write(f"{chrom}\t{s}\t{e}\t{name}\t1000\t{strand}\t{s}\t{s}\t0,0,0\t"
                     f"{len(exons)}\t{sizes},\t{starts},\n")


@pytest.fixture(scope="module")
def clf(tmp_path_factory):
    d = tmp_path_factory.mktemp("ref")
    gtf = d / "ref.gtf"
    write_gtf(gtf, REF)
    return Classifier(Reference(read_gtf_transcripts(str(gtf))))


def iso(exons, strand="+", chrom="chr1", name="q"):
    return Transcript(name, "", chrom, strand, tuple(exons))


CASES = [
    # name, exons, strand, category, subcategory, ref_transcript, ref_gene
    ("fsm_cim", [(100, 200), (300, 400), (500, 600)], "+", "FSM", "CIM", "A1", "A"),
    ("fsm_cim_within_window", [(130, 200), (300, 400), (500, 640)], "+", "FSM", "CIM", "A1", "A"),
    ("fsm5_plus", [(20, 200), (300, 400), (500, 600)], "+", "FSM", "FSM5", "A1", "A"),
    ("fsm3_plus", [(100, 200), (300, 400), (500, 700)], "+", "FSM", "FSM3", "A1", "A"),
    ("fsmte_plus", [(20, 200), (300, 400), (500, 700)], "+", "FSM", "FSMTE", "A1", "A"),
    # minus strand: the 5' end is the right-hand end
    ("fsm5_minus", [(5000, 5100), (5300, 5400), (5500, 5700)], "-", "FSM", "FSM5", "C1", "C"),
    ("fsm3_minus", [(4900, 5100), (5300, 5400), (5500, 5600)], "-", "FSM", "FSM3", "C1", "C"),
    ("ism", [(300, 400), (500, 600)], "+", "ISM", "", "A1", "A"),
    ("ism_other_ref", [(550, 600), (700, 800)], "+", "ISM", "", "A2", "A"),
    ("nic_new_combination", [(100, 200), (300, 400), (500, 600), (700, 800)], "+",
     "NIC", "", "", "A"),
    ("nnc_novel_acceptor", [(100, 200), (320, 400), (500, 600)], "+", "NNC", "", "", "A"),
    # no shared splice site: NNC when it overlaps exons of a spliced gene (SQANTI3)
    ("nnc_no_shared_site", [(150, 250), (450, 550)], "+", "NNC", "", "", "A"),
    ("genic_multi_mono_gene", [(20100, 20300), (20500, 20900)], "+", "genic", "", "", "T"),
    ("genic_intron_multi", [(410, 430), (450, 490)], "+", "genic_intron", "", "", "A"),
    ("antisense_multi", [(5020, 5080), (5320, 5380)], "+", "antisense", "", "", ""),
    ("intergenic_multi", [(50000, 50100), (50300, 50400)], "+", "intergenic", "", "", ""),
    ("intergenic_no_ref_contig", [(100, 200), (300, 400)], "+", "intergenic", "", "", ""),
    ("fusion", [(100, 200), (300, 400), (10300, 10400), (10500, 10600)], "+",
     "fusion", "", "", "A,B"),
    # mono-exon
    ("mono_fsm_cim", [(2010, 2990)], "+", "FSM", "CIM", "M1", "M"),
    ("mono_fsm3", [(2000, 3100)], "+", "FSM", "FSM3", "M1", "M"),
    ("mono_ism", [(310, 390)], "+", "ISM", "", "A1", "A"),
    ("mono_nic_intron_retention", [(150, 350)], "+", "NIC", "", "", "A"),
    ("mono_nic_ir_ends_in_intron", [(150, 450)], "+", "NIC", "", "", "A"),
    ("mono_genic", [(180, 250)], "+", "genic", "", "", "A"),
    ("mono_genic_intron", [(210, 290)], "+", "genic_intron", "", "", "A"),
    ("mono_antisense", [(120, 180)], "-", "antisense", "", "", ""),
    # inside an opposite-strand intron, no exon touched: genic_intron (SQANTI3)
    ("mono_genic_intron_opposite", [(410, 490)], "-", "genic_intron", "", "", ""),
    ("mono_intergenic", [(40000, 41000)], "+", "intergenic", "", "", ""),
]


@pytest.mark.parametrize("name,exons,strand,cat,sub,rtx,rgene", CASES, ids=[c[0] for c in CASES])
def test_category(clf, name, exons, strand, cat, sub, rtx, rgene):
    chrom = "chr2" if name == "intergenic_no_ref_contig" else "chr1"
    call = clf.classify(iso(exons, strand, chrom))
    assert (call.category, call.subcategory, call.ref_transcript, call.ref_gene) == \
        (cat, sub, rtx, rgene)


def test_every_category_and_subcategory_is_covered():
    cats = {c[3] for c in CASES}
    subs = {c[4] for c in CASES if c[4]}
    assert cats == set(CATEGORIES)
    assert subs == {"CIM", "FSM5", "FSM3", "FSMTE"}


def test_unstranded_isoform_matches_either_strand(clf):
    assert clf.classify(iso([(5000, 5100), (5300, 5400), (5500, 5600)], ".")).category == "FSM"


def test_fsm_subcategory_prefers_best_reference_deterministically():
    a = Transcript("r1", "g", "c", "+", ((100, 200), (300, 400)))
    b = Transcript("r2", "g", "c", "+", ((40, 200), (300, 460)))
    sub, best = fsm_subcategory(45, 455, "+", [a, b])
    assert (sub, best.tid) == ("CIM", "r2")
    sub2, best2 = fsm_subcategory(45, 455, "+", [b, a])
    assert (sub2, best2.tid) == (sub, best.tid)


def test_histogram_bins_last_open_ended():
    edges = HISTOGRAM_EDGES["isoforms_per_gene"]
    counts = histogram([1, 1, 2, 7, 21, 500], edges)
    assert len(counts) == len(edges)
    assert counts[0] == 2 and counts[1] == 1 and counts[edges.index(6)] == 1
    assert counts[-1] == 2


def test_keep_contig():
    assert keep_contig("chr1") and keep_contig("chrX") and keep_contig("SIRV1")
    assert keep_contig("ERCC-00002") and keep_contig("17")
    for bad in ("chrM", "MT", "chr1_KI270706v1_random", "chrUn_GL000220v1", "GL000195.1",
                "KI270728.1", "chrEBV"):
        assert not keep_contig(bad), bad


def test_normalise_version():
    assert normalise_version("FLAIR 3.0.1") == "3.0.1"
    assert normalise_version("flair 3.0.0+dev") == "3.0.0+dev"
    assert normalise_version("2.2.0") == "2.2.0"
    assert normalise_version(None) == ""


def test_gene_assignment_fallbacks(tmp_path):
    gtf = tmp_path / "iso.gtf"
    write_gtf(gtf, [("ENST1", "ENSG1", "chr1", "+", [(0, 10)])])
    from flair_eval_lib.annotation import read_gtf_gene_map
    isos = [Transcript("ENST1_ENSG1", "", "chr1", "+", ((0, 10),)),
            Transcript("FLT9", "", "chr1", "+", ((100, 200),)),
            Transcript("FLT10", "", "chr1", "+", ((150, 250),)),
            Transcript("FLT11", "", "chr1", "+", ((900, 950),))]
    assign_isoform_genes(isos, read_gtf_gene_map(str(gtf)))
    assert isos[0].gene == "ENSG1"
    assert isos[1].gene == isos[2].gene != isos[3].gene
    assert isos[3].gene.startswith("novel_gene_")


def test_bed_errors(tmp_path):
    p = tmp_path / "bad.bed"
    p.write_text("chr1\t10\t5\tx\t0\t+\n")
    with pytest.raises(InputError):
        read_bed12(str(p))
    p.write_text("chr1\t1\t5\tx\t0\t+\nchr1\t1\t5\tx\t0\t+\n")
    with pytest.raises(InputError, match="duplicate"):
        read_bed12(str(p))


def test_fasta_motifs(tmp_path):
    #           0         1         2
    #           0123456789012345678901234
    seq = "AAAAAGTAAAAAAAGAAAAACTAAAAAACAAAAA"
    fa = tmp_path / "g.fa"
    fa.write_text(">chrT\n" + "\n".join(seq[i:i + 10] for i in range(0, len(seq), 10)) + "\n")
    f = Fasta(str(fa))
    assert f.fetch("chrT", 5, 7) == "GT"
    assert f.fetch("chrT", 8, 12) == seq[8:12]
    assert f.is_canonical("chrT", "+", 5, 15)       # GT ... AG
    assert not f.is_canonical("chrT", "-", 5, 15)
    assert f.is_canonical("chrT", "-", 20, 29)      # CT ... AC on + = GT-AG on -
    assert f.is_canonical("chrT", ".", 20, 29)
    assert f.fetch("chrT", 30, 99) == ""
    f.close()


# --------------------------------------------------------------------- CLI end to end
def _run_cli(tmp_path, extra=(), isoforms=None, run=None, env=None, sample="s1"):
    gtf = tmp_path / "ref.gtf"
    write_gtf(gtf, REF)
    isoforms = isoforms or [
        ("i_fsm", "chr1", "+", [(100, 200), (300, 400), (500, 600)]),
        ("i_fsm5", "chr1", "+", [(20, 200), (300, 400), (500, 600)]),
        ("i_ism", "chr1", "+", [(300, 400), (500, 600)]),
        ("i_nnc", "chr1", "+", [(100, 200), (320, 400), (500, 600)]),
        ("i_mono", "chr1", "+", [(2010, 2990)]),
        ("i_far", "chr1", "+", [(10000, 10100), (10300, 10400), (10500, 10600)]),
        ("i_mt", "chrM", "+", [(10, 500)]),
    ]
    bed = tmp_path / "iso.bed"
    write_bed(bed, isoforms)
    igtf = tmp_path / "iso.gtf"
    write_gtf(igtf, [(n, "G_" + n.split("_")[1], c, s, e) for n, c, s, e in isoforms])
    rj = tmp_path / "run.json"
    rj.write_text(json.dumps(run or {"id": "flair-3.0.1", "flair": {"release": "3.0.1"},
                                     "label": "FLAIR 3.0.1", "tags": {"version": "3.0.1"}}))
    ej = tmp_path / "env.json"
    ej.write_text(json.dumps(env or {"run_id": "flair-3.0.1", "source": "release",
                                     "requested": "3.0.1", "resolved_commit": None,
                                     "flair_version": "FLAIR 3.0.1"}))
    out = tmp_path / "bundle"
    cmd = [sys.executable, os.path.join(BIN, "structure_metrics.py"),
           "--isoforms-bed", str(bed), "--isoforms-gtf", str(igtf), "--gtf", str(gtf),
           "--sample", sample, "--run-json", str(rj), "--env-json", str(ej),
           "--out", str(out), *extra]
    return subprocess.run(cmd, capture_output=True, text=True), out


def test_cli_writes_valid_bundle(tmp_path):
    res, out = _run_cli(tmp_path)
    assert res.returncode == 0, res.stderr
    b = json.loads((out / "bundle.json").read_text())
    assert validate(b, load_schema()) == []
    try:
        import jsonschema
        jsonschema.validate(b, load_schema())
    except ImportError:
        pass
    assert b["sample"] == "s1" and b["run_id"] == "flair-3.0.1"
    assert b["flair"] == {"source": "release", "requested": "3.0.1", "version": "3.0.1",
                          "commit": None}
    s = b["summary"]
    assert s["isoforms"] == 6                       # chrM dropped
    assert b["inputs"]["isoforms_dropped_contigs"] == 1
    assert s["mono_exon"] == 1 and s["multi_exon"] == 5
    assert s["ref_transcripts_total"] == 6          # chrM reference excluded
    assert s["ref_transcripts_fsm"] == 3            # A1, M1, B1
    assert s["ref_transcripts_cim"] == 3
    assert s["ref_genes_detected"] == 3             # A, M, B
    assert b["categories"]["FSM"] == 4 and b["categories"]["ISM"] == 1
    assert b["categories"]["NNC"] == 1
    assert b["fsm_subcategories"] == {"CIM": 3, "FSM5": 1, "FSM3": 0, "FSMTE": 0}
    j = b["junctions"]
    # distinct junctions: (200,300) (400,500) (200,320) (10100,10300) (10400,10500)
    assert (j["total"], j["known"], j["novel"], j["motif_checked"]) == (5, 4, 1, False)
    for k, h in b["histograms"].items():
        assert h["edges"] == HISTOGRAM_EDGES[k] and len(h["counts"]) == len(h["edges"])
    assert sum(b["histograms"]["exons_per_isoform"]["counts"]) == 6
    assert sum(b["histograms"]["isoforms_per_gene"]["counts"]) == s["genes"]

    rows = (out / "tables" / "isoforms.tsv").read_text().splitlines()
    assert rows[0].split("\t") == ["isoform_id", "gene_id", "chrom", "start", "end", "strand",
                                   "n_exons", "length", "category", "subcategory",
                                   "ref_transcript", "ref_gene"]
    assert len(rows) == 7
    first = rows[1].split("\t")
    assert first[0] == "i_fsm" and first[8:] == ["FSM", "CIM", "A1", "A"]
    cats = [r.split("\t") for r in (out / "tables" / "categories.tsv").read_text().splitlines()]
    assert cats[0] == ["category", "subcategory", "n"]
    assert sum(int(r[2]) for r in cats[1:]) == 6
    assert ["FSM", "FSM5", "1"] in cats


def test_cli_regions_and_genome(tmp_path):
    reg = tmp_path / "r.bed"
    reg.write_text("chr1\t0\t1000\n")
    seq = ["A"] * 1200
    seq[200:202], seq[318:320] = "GT", "AG"          # novel intron (200,320) is GT-AG
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "".join(seq) + "\n")
    res, out = _run_cli(tmp_path, ["--regions", str(reg), "--genome", str(fa)])
    assert res.returncode == 0, res.stderr
    b = json.loads((out / "bundle.json").read_text())
    assert b["summary"]["isoforms"] == 4            # mono at 2000 and B at 10000 are outside
    assert b["summary"]["ref_transcripts_total"] == 2   # A1, A2
    assert b["inputs"]["isoforms_outside_regions"] == 2
    j = b["junctions"]
    assert j["motif_checked"] is True
    assert (j["known"], j["novel_canonical"], j["novel_noncanonical"]) == (2, 1, 0)
    assert j["known"] + j["novel_canonical"] + j["novel_noncanonical"] == j["total"]


@pytest.mark.parametrize("bad", ["sample", "run_id_mismatch", "source"])
def test_cli_rejects_bad_input(tmp_path, bad):
    if bad == "sample":
        res, _ = _run_cli(tmp_path, sample="bad sample!")
    elif bad == "run_id_mismatch":
        res, _ = _run_cli(tmp_path, env={"run_id": "other", "source": "release",
                                         "requested": "3.0.1", "flair_version": "3.0.1"})
    else:
        res, _ = _run_cli(tmp_path, env={"run_id": "flair-3.0.1", "source": "pip",
                                         "requested": "3.0.1", "flair_version": "3.0.1"})
    assert res.returncode != 0
    assert "error" in res.stderr


def test_cli_help():
    res = subprocess.run([sys.executable, os.path.join(BIN, "structure_metrics.py"), "--help"],
                         capture_output=True, text=True)
    assert res.returncode == 0 and "--isoforms-bed" in res.stdout


def test_regions_overlap():
    r = Regions({"chr1": [(100, 200), (150, 300)]})
    assert r.by_chrom["chr1"] == [(100, 300)]
    assert r.overlaps("chr1", 299, 400) and not r.overlaps("chr1", 300, 400)
    assert not r.overlaps("chr2", 0, 10**9)


def test_schema_validator_catches_errors():
    schema = load_schema()
    assert validate({"schema": "x"}, schema)
    errs = validate({"schema": "flair-eval/bundle/1", "sample": "bad name"}, schema)
    assert any("does not match" in e for e in errs)
    assert any("missing required key 'summary'" in e for e in errs)
