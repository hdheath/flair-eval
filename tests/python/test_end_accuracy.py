"""Tests for bin/end_accuracy.py and bin/flair_eval_lib/end_accuracy.py (CONTRACTS section 12).

The numbers come from the real judge, ted_end_precision.py, run on a hand-built case where every
hit and miss is worked out below. Coordinates are BED-style (0-based, half-open)."""

import filecmp
import json
import os
import subprocess
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
BIN = os.path.join(ROOT, "bin")
VENDORED = os.path.join(ROOT, "modules", "end_accuracy", "resources", "usr", "bin")
JUDGE = os.path.join(VENDORED, "ted_end_precision.py")
sys.path.insert(0, BIN)

from flair_eval_lib.annotation import InputError  # noqa: E402
from flair_eval_lib.end_accuracy import (  # noqa: E402
    TABLE_COLUMNS, EndSet, plan_calls, read_sets, read_summary, regions_args)
from flair_eval_lib.schema import load_schema, validate  # noqa: E402

pytest.importorskip("scipy", reason="the judge's signal_utils imports scipy")


def write(path, lines):
    path.write_text("".join(line + "\n" for line in lines))
    return str(path)


# ------------------------------------------------------------------ the vendored judge
@pytest.mark.parametrize("name", ["peak_matching.py", "signal_utils.py"])
def test_vendored_judge_helpers_are_unchanged(name):
    """The helpers are used as written: byte-identical to the legacy copy (= flair4_ted 40fdd21)."""
    assert filecmp.cmp(os.path.join(VENDORED, name),
                       os.path.join(ROOT, "legacy", "bin", "evaluation", name), shallow=False)


@pytest.mark.parametrize("peaks", [("CAGE", "dRNA"), ("CAGE", None), (None, "dRNA"), (None, None)])
def test_extended_judge_keeps_the_legacy_columns(case, tmp_path, peaks):
    """The judge gained per-chain columns (ted-clustering ticket 12). Its legacy columns, in order,
    and per_junction_chain.tsv are exactly the legacy judge's (flair4_ted 40fdd21)."""
    files = {"CAGE": os.path.join(case["tmp"], "CAGE.truth.bed"),
             "dRNA": os.path.join(case["tmp"], "dRNA.truth.bed")}
    args = ["--isoforms-bed", case["iso"], "--gtf", case["gtf"], "--mode", "m"]
    for flag, name in zip(("--peaks-5prime", "--peaks-3prime"), peaks):
        if name:
            args += [flag, files[name]]
    outs = {}
    for label, script in (("legacy", os.path.join(ROOT, "legacy", "bin", "evaluation",
                                                  "ted_end_precision.py")),
                          ("vendored", JUDGE)):
        out = tmp_path / label
        r = subprocess.run([sys.executable, script, *args, "--outdir", str(out)],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        outs[label] = out
    legacy = (outs["legacy"] / "precision_recall_summary.tsv").read_text().splitlines()
    new = (outs["vendored"] / "precision_recall_summary.tsv").read_text().splitlines()
    width = len(legacy[0].split("\t"))
    assert [line.split("\t")[:width] for line in new] == [line.split("\t") for line in legacy]
    assert len(new[0].split("\t")) > width
    assert ((outs["vendored"] / "per_junction_chain.tsv").read_text()
            == (outs["legacy"] / "per_junction_chain.tsv").read_text())


def test_judge_is_executable():
    assert os.access(JUDGE, os.X_OK)


# ------------------------------------------------------------------ planning
def es(i, kind, end):
    return EndSet(id=i, kind=kind, end=end, assay=i, path="/dev/null")


def test_plan_pairs_every_tss_set_with_every_tes_set_per_kind():
    sets = [es("CAGE", "truth", "tss"), es("PROcap", "truth", "tss"), es("dRNA", "truth", "tes"),
            es("refTSS", "reference", "tss"), es("PolyA_DB", "reference", "tes")]
    names = [(c.kind, c.name) for c in plan_calls(sets)]
    assert names == [("reference", "refTSS__PolyA_DB"),
                     ("truth", "CAGE__dRNA"), ("truth", "PROcap__dRNA")]


def test_plan_single_end_sets_run_alone_and_kinds_never_mix():
    sets = [es("CAGE", "truth", "tss"), es("RAMPAGE", "truth", "tss"),
            es("PolyA_DB", "reference", "tes")]
    names = [c.name for c in plan_calls(sets)]
    assert names == ["PolyA_DB", "CAGE", "RAMPAGE"]


def test_regions_args(tmp_path):
    bed = write(tmp_path / "r.bed", ["# comment", "chr19\t0\t2000000", "chr1\t10\t20\tname"])
    assert regions_args(bed) == ["chr19:0-2000000", "chr1:10-20"]
    assert regions_args(None) == []


def test_read_sets_rejects_bad_entries(tmp_path):
    f = tmp_path / "x.bed"
    f.write_text("")
    good = {"id": "CAGE", "kind": "truth", "end": "tss", "path": str(f)}
    for bad, msg in [({"kind": "other"}, "kind"), ({"end": "5p"}, "end"),
                     ({"id": "bad id"}, "id"), ({"path": str(tmp_path / "no.bed")}, "does not exist")]:
        p = tmp_path / "s.json"
        p.write_text(json.dumps([dict(good, **bad)]))
        with pytest.raises(InputError, match=msg):
            read_sets(str(p))
    p.write_text(json.dumps([good, good]))
    with pytest.raises(InputError, match="twice"):
        read_sets(str(p))


# ------------------------------------------------------------------ end to end with the judge
@pytest.fixture
def case(tmp_path):
    """Four scored isoforms (+ one on chrM that the judge drops):
         iso1 chr1 + [1000,2000) 2 exons   5' 1000   3' 2000
         iso2 chr1 + [5000,6000) 2 exons   5' 5000   3' 6000
         iso3 chr1 - [10000,11000) 2 exons 5' 11000  3' 10000
         iso4 chr2 + [1000,1500) mono      5' 1000   3' 1500
       (the judge takes BED start/end as the ends: tss_tts() in signal_utils.py)"""
    iso = write(tmp_path / "isoforms.bed", [
        "chr1\t1000\t2000\tiso1\t0\t+\t1000\t2000\t0\t2\t200,200,\t0,800,",
        "chr1\t5000\t6000\tiso2\t0\t+\t5000\t6000\t0\t2\t100,100,\t0,900,",
        "chr1\t10000\t11000\tiso3\t0\t-\t10000\t11000\t0\t2\t100,100,\t0,900,",
        "chr2\t1000\t1500\tiso4\t0\t+\t1000\t1500\t0\t1\t500,\t0,",
        "chrM\t100\t900\tiso5\t0\t+\t100\t900\t0\t1\t800,\t0,",
    ])
    gtf = write(tmp_path / "ann.gtf", [
        'chr1\tt\ttranscript\t1001\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
        'chr1\tt\texon\t1001\t1200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
        'chr1\tt\texon\t1801\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
    ])
    hdr = "#chrom\tstart\tend\tcall_id\tscore\tstrand\text_start\text_end"
    # 5' truth: iso1 hit (20 bp), iso2 miss (100 bp), iso3 hit (10 bp), chr3 call unmatched
    cage = write(tmp_path / "CAGE.truth.bed", [
        hdr, "chr1\t1020\t1021\tc1\t1\t+\t1015\t1030", "chr1\t5100\t5101\tc2\t1\t+\t5100\t5101",
        "chr1\t10990\t10991\tc3\t1\t-\t10980\t10995", "chr3\t500\t501\tc4\t1\t+\t500\t501"])
    # 3' truth: iso1 hit (10 bp), iso2 hit (0 bp), iso4 miss (100 bp), iso3's call on the wrong strand
    drna = write(tmp_path / "dRNA.truth.bed", [
        hdr, "chr1\t2010\t2011\td1\t1\t+\t2010\t2011", "chr1\t6000\t6001\td2\t1\t+\t6000\t6001",
        "chr2\t1600\t1601\td3\t1\t+\t1600\t1601", "chr1\t10000\t10001\td4\t1\t+\t10000\t10001"])
    # references: refTSS hits iso1 (inside) and iso4 (41 bp from the interval edge);
    # PolyA_DB hits iso2 (20 bp), listed twice as for two genes
    reftss = write(tmp_path / "refTSS.bed", ["chr1\t990\t1010\tr1\t0\t+", "chr2\t900\t960\tr2\t0\t+"])
    pas = write(tmp_path / "pas.bed", ["chr1\t5980\t5981\tp1|A\t0\t+", "chr1\t5980\t5981\tp1|B\t0\t+"])
    sets = [
        {"id": "CAGE", "kind": "truth", "end": "tss", "assay": "CAGE", "path": cage,
         "source": "/truth/W.CAGE.truth.bed", "provisional": True},
        {"id": "dRNA", "kind": "truth", "end": "tes", "assay": "dRNA", "path": drna, "provisional": True},
        {"id": "refTSS", "kind": "reference", "end": "tss", "assay": "refTSS", "path": reftss},
        {"id": "PolyA_DB", "kind": "reference", "end": "tes", "assay": "PolyA_DB", "path": pas},
    ]
    sj = tmp_path / "sets.json"
    sj.write_text(json.dumps(sets))
    return {"iso": iso, "gtf": gtf, "sets": str(sj), "out": str(tmp_path / "out"), "tmp": tmp_path}


def run_cli(c, *extra):
    return subprocess.run([sys.executable, os.path.join(BIN, "end_accuracy.py"),
                           "--isoforms-bed", c["iso"], "--gtf", c["gtf"], "--sets", c["sets"],
                           "--out", c["out"], "--judge", JUDGE, *extra],
                          capture_output=True, text=True)


@pytest.fixture
def result(case):
    r = run_cli(case)
    assert r.returncode == 0, r.stderr
    with open(os.path.join(case["out"], "end_accuracy.json")) as fh:
        return json.load(fh)


def row(sec, set_id):
    return next(r for r in sec["ends"] if r["set"] == set_id)


def test_truth_precision_recall_f1(result):
    cage, drna = row(result, "CAGE"), row(result, "dRNA")
    assert cage["isoforms"] == 4                      # chrM isoform dropped
    assert cage["precision"] == pytest.approx(0.5)    # iso1, iso3 of 4
    assert (cage["sites_hit"], cage["sites"]) == (2, 4)
    assert cage["recall"] == pytest.approx(0.5)
    assert cage["f1"] == pytest.approx(0.5)
    assert drna["precision"] == pytest.approx(0.5)    # iso1, iso2
    assert (drna["sites_hit"], drna["sites"]) == (2, 4)
    assert drna["recall"] == pytest.approx(0.5)


def test_reference_rows_are_precision_only(result):
    ref, pas = row(result, "refTSS"), row(result, "PolyA_DB")
    assert ref["precision"] == pytest.approx(0.5)     # iso1, iso4
    assert pas["precision"] == pytest.approx(0.25)    # iso2
    for r in (ref, pas):
        assert r["recall"] is None and r["f1"] is None and r["sites_hit"] is None


def test_pairs(result):
    pairs = {(p["tss_set"], p["tes_set"]): p for p in result["pairs"]}
    assert set(pairs) == {("CAGE", "dRNA"), ("refTSS", "PolyA_DB")}
    t = pairs[("CAGE", "dRNA")]
    assert t["both_hit"] == 1 and t["isoforms"] == 4  # only iso1 has both ends on truth
    assert t["precision"] == pytest.approx(0.25)
    assert pairs[("refTSS", "PolyA_DB")]["both_hit"] == 0


def test_section_meta_and_order(result):
    assert [(r["kind"], r["end"], r["set"]) for r in result["ends"]] == [
        ("reference", "tss", "refTSS"), ("reference", "tes", "PolyA_DB"),
        ("truth", "tss", "CAGE"), ("truth", "tes", "dRNA")]
    assert result["provisional"] is True
    assert result["scorer"]["script"] == "ted_end_precision.py" and result["scorer"]["window"] == 50
    meta = {s["id"]: s for s in result["sets"]}
    assert meta["CAGE"]["file"] == "W.CAGE.truth.bed" and len(meta["CAGE"]["md5"]) == 32


def test_section_validates_against_bundle_schema(result):
    schema = load_schema()
    sec_schema = schema["properties"]["sections"]["properties"]["end_accuracy"]
    assert validate(result, sec_schema) == []


def test_single_end_numbers_do_not_depend_on_the_partner(case, tmp_path):
    """A 5' set alone gives the same 5' numbers as when paired with a 3' set."""
    full = run_cli(case)
    assert full.returncode == 0, full.stderr
    sets = json.loads(open(case["sets"]).read())
    alone = tmp_path / "alone.json"
    alone.write_text(json.dumps([s for s in sets if s["id"] == "CAGE"]))
    out2 = str(tmp_path / "out2")
    r = subprocess.run([sys.executable, os.path.join(BIN, "end_accuracy.py"), "--isoforms-bed",
                        case["iso"], "--gtf", case["gtf"], "--sets", str(alone), "--out", out2,
                        "--judge", JUDGE], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    a = json.load(open(os.path.join(case["out"], "end_accuracy.json")))
    b = json.load(open(os.path.join(out2, "end_accuracy.json")))
    assert row(a, "CAGE") == row(b, "CAGE")
    assert b["pairs"] == []


def test_table(case, result):
    with open(os.path.join(case["out"], "tables", "end_accuracy.tsv")) as fh:
        lines = [line.rstrip("\n").split("\t") for line in fh]
    assert lines[0] == TABLE_COLUMNS
    body = [dict(zip(TABLE_COLUMNS, x)) for x in lines[1:]]
    assert [(b["kind"], b["end"], b["set"], b["partner"]) for b in body] == [
        ("reference", "tss", "refTSS", ""), ("reference", "tes", "PolyA_DB", ""),
        ("truth", "tss", "CAGE", ""), ("truth", "tes", "dRNA", ""),
        ("reference", "pair", "refTSS", "PolyA_DB"), ("truth", "pair", "CAGE", "dRNA")]
    assert body[0]["recall"] == "" and float(body[2]["recall"]) == pytest.approx(0.5)
    # the judge's own outputs are kept per call
    assert os.path.isfile(os.path.join(case["out"], "judge", "CAGE__dRNA",
                                       "precision_recall_summary.tsv"))


def test_regions_restrict_sites_and_isoforms(case):
    reg = write(case["tmp"] / "reg.bed", ["chr1\t0\t3000"])
    r = run_cli(case, "--regions", reg)
    assert r.returncode == 0, r.stderr
    sec = json.load(open(os.path.join(case["out"], "end_accuracy.json")))
    cage = row(sec, "CAGE")
    assert cage["isoforms"] == 1 and (cage["sites_hit"], cage["sites"]) == (1, 1)
    assert sec["scorer"]["regions"] == ["chr1:0-3000"]


def test_minus_strand_5prime_end_is_the_bed_end(case, tmp_path):
    """Documents the judge's convention (not changed): on '-' the 5' end is BED end, one past the
    first transcribed base (end - 1). A call 50 bp from end - 1 but 51 bp from end is a miss."""
    iso = write(tmp_path / "m.bed", ["chr1\t10000\t11000\tm\t0\t-\t10000\t11000\t0\t1\t1000,\t0,"])
    t = write(tmp_path / "t.bed", ["chr1\t10949\t10950\tc\t1\t-"])        # 10999 - 50
    sj = tmp_path / "m.json"
    sj.write_text(json.dumps([{"id": "CAGE", "kind": "truth", "end": "tss", "path": t}]))
    out = str(tmp_path / "mo")
    r = subprocess.run([sys.executable, os.path.join(BIN, "end_accuracy.py"), "--isoforms-bed", iso,
                        "--gtf", case["gtf"], "--sets", str(sj), "--out", out, "--judge", JUDGE],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert row(json.load(open(os.path.join(out, "end_accuracy.json"))), "CAGE")["precision"] == 0


def test_missing_judge_is_a_clear_error(case):
    r = subprocess.run([sys.executable, os.path.join(BIN, "end_accuracy.py"), "--isoforms-bed",
                        case["iso"], "--gtf", case["gtf"], "--sets", case["sets"], "--out",
                        case["out"], "--judge", "/no/such/judge.py"], capture_output=True, text=True)
    assert r.returncode == 2 and "not found" in r.stderr


def test_read_summary_types(tmp_path):
    p = write(tmp_path / "s.tsv", ["transcriptome_mode\tn_isoforms_total\t5prime_precision\t3prime_recall",
                                   "x\t3\t0.500000\t"])
    assert read_summary(p) == {"transcriptome_mode": "x", "n_isoforms_total": 3,
                               "5prime_precision": 0.5, "3prime_recall": None}
