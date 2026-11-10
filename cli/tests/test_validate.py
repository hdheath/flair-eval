import json

import pytest

from flair_eval import validate as V


def codes(issues):
    return sorted(((i.path, i.code) for i in issues), key=lambda x: ([str(p) for p in x[0]], x[1]))


def test_valid_sample_and_runs(sample, two_runs):
    assert V.validate_samples([sample]) == []
    assert V.validate_runs({"runs": two_runs}) == []
    bam = {"sample": "s", "bam": "x.bam", "genome": "g", "gtf": "a", "platform": "ont",
           "library": "dRNA", "cage": "c.bed", "pas": "p.bed", "junctions": "SJ.out.tab",
           "spikein_gtf": "sirv.gtf"}
    assert V.validate_samples([bam]) == []


@pytest.mark.parametrize("row,expected", [
    ({"sample": "s", "genome": "g", "gtf": "a"}, [(("samples", 0, "reads"), "either_reads_bam")]),
    ({"sample": "a b", "reads": "r", "genome": "g", "gtf": "a"},
     [(("samples", 0, "sample"), "pattern")]),
    ({"sample": "s", "reads": "r", "gtf": "a"}, [(("samples", 0, "genome"), "required")]),
    ({"sample": "s", "reads": "r", "genome": "g", "gtf": "a", "platform": "nanopore"},
     [(("samples", 0, "platform"), "enum")]),
    ({"sample": "s", "reads": "r", "genome": "g", "gtf": "a", "library": "cdna"},
     [(("samples", 0, "library"), "enum")]),
    ({"sample": "s", "reads": "r", "genome": "g", "gtf": "a", "colour": "x"},
     [(("samples", 0, "colour"), "unknown_key")]),
])
def test_sample_errors(row, expected):
    assert codes(V.validate_samples([row])) == expected


def test_duplicate_samples_and_header(sample):
    issues = V.validate_samples([sample, dict(sample)], header=["sample", "reads", "genome",
                                                                "gtf", "extra"])
    assert codes(issues) == [(("samples", 1, "sample"), "duplicate"),
                             (("samples", "header", "extra"), "unknown_key")]
    assert "already used by sample 1" in [i for i in issues if i.code == "duplicate"][0].message
    assert codes(V.validate_samples([])) == [(("samples",), "min_items")]


@pytest.mark.parametrize("flair,needle", [
    ({"git": "x"}, "git needs a ref"),
    ({"release": "3.0.1", "env": "/x"}, "only one FLAIR source"),
    ({"ref": "dev"}, "ref only goes with git"),
    ({"release": 3.0}, 'quotes, e.g. release: "3.0"'),
    ({"release": ""}, "release should be non-empty"),
    ({"git": "u", "ref": ""}, "ref should be non-empty"),
    ({"version": "3"}, "unknown key 'version' in flair"),
    ({}, "exactly one FLAIR source"),
    ("3.0.1", "exactly one FLAIR source"),
])
def test_flair_source_messages(flair, needle):
    issues = V.validate_runs({"runs": [{"id": "a", "flair": flair}]})
    assert codes(issues) == [(("runs", 0, "flair"), "flair_source")]
    assert needle in issues[0].message


def test_run_errors():
    doc = {"runs": [{"id": "bad id", "flair": {"release": "3"}, "quantify": "yes",
                     "tags": {"n": 3}, "colour": 1},
                    {"flair": {"release": "3"}},
                    {"id": "x", "flair": {"release": "3"}}, {"id": "x", "flair": {"env": "/e"}}]}
    assert codes(V.validate_runs(doc)) == [
        (("runs", 0, "colour"), "unknown_key"), (("runs", 0, "id"), "pattern"),
        (("runs", 0, "quantify"), "type"), (("runs", 0, "tags", "n"), "type"),
        (("runs", 1, "id"), "required"), (("runs", 3, "id"), "duplicate")]
    msgs = {i.path: i.message for i in V.validate_runs(doc)}
    assert msgs[("runs", 0, "quantify")] == "should be true or false, not text"
    assert 'put it in quotes, e.g. "3"' in msgs[("runs", 0, "tags", "n")]


@pytest.mark.parametrize("doc,expected", [
    ({}, [(("runs",), "required")]),
    ({"runs": []}, [(("runs",), "min_items")]),
    ([{"id": "a"}], [(("runs",), "type")]),
    ("x", [(("runs",), "type")]),
    ({"runs": [], "extra": 1}, [(("extra",), "unknown_key"), (("runs",), "min_items")]),
])
def test_runs_document_shape(doc, expected):
    assert codes(V.validate_runs(doc)) == expected


def test_issue_formatting():
    i = V.Issue(("samples", 2, "gtf"), "required", "missing required column 'gtf'")
    assert str(i) == "error: samples.csv row 3: gtf: missing required column 'gtf'"
    assert V.format_path(("runs", 0, "flair")) == "runs.yaml run 1: flair"
    assert V.format_path(()) == ""


def test_check_paths(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "g.fa").write_text(">c\nA\n")
    (tmp_path / "x.bam").write_text("")
    rows = [{"sample": "s", "bam": "x.bam", "genome": "g.fa", "gtf": "missing.gtf",
             "cage": "s3://bucket/cage.bed"}]
    got = codes(V.check_paths(rows))
    assert got == [(("samples", 0, "bam"), "missing_index"), (("samples", 0, "gtf"), "missing_file")]
    (tmp_path / "x.bam.bai").write_text("")
    assert codes(V.check_paths(rows)) == [(("samples", 0, "gtf"), "missing_file")]
    assert all(i.level == "warning" for i in V.check_paths(rows))


def test_check_files(tmp_path):
    assert codes(V.check_files(tmp_path / "no.csv", tmp_path / "no.yaml")) == [
        (("runs",), "missing_file"), (("samples",), "missing_file")]
    bad = tmp_path / "bad.yaml"
    bad.write_text("runs: [")
    assert [i.code for i in V.check_files(None, bad)] == ["yaml"]


BUNDLE = {
    "schema": "flair-eval/bundle/1", "sample": "s", "run_id": "r", "label": "r", "tags": {},
    "flair": {"source": "release", "requested": "3.0.1", "version": "3.0.1", "commit": None},
    "pipeline": {"version": "2.0.0", "commit": None, "created": "2026-11-09T00:00:00Z"},
    "inputs": {}, "summary": {k: 0 for k in ["isoforms", "genes", "loci", "mono_exon",
                                              "multi_exon", "ref_transcripts_total",
                                              "ref_transcripts_fsm", "ref_transcripts_cim",
                                              "ref_genes_detected"]},
    "categories": {k: 0 for k in ["FSM", "ISM", "NIC", "NNC", "genic", "genic_intron",
                                  "antisense", "intergenic", "fusion"]},
    "fsm_subcategories": {"CIM": 0, "FSM5": 0, "FSM3": 0, "FSMTE": 0},
    "junctions": {"total": 0, "known": 0, "novel_canonical": 0, "novel_noncanonical": 0},
    "histograms": {k: {"edges": [0, 1], "counts": [0, 0]} for k in
                   ["isoforms_per_gene", "transcript_length", "exon_length", "exons_per_isoform"]},
    "compute": {"cpu_hours": None, "peak_rss_gb": None, "wall_minutes": None},
    "sections": {},
}


def test_bundle_validation():
    assert V.validate_bundle(json.loads(json.dumps(BUNDLE))) == []
    bad = json.loads(json.dumps(BUNDLE))
    bad["schema"] = "flair-eval/bundle/0"
    bad["summary"]["isoforms"] = -1
    del bad["junctions"]
    assert codes(V.validate_bundle(bad)) == [(("bundle", "junctions"), "required"),
                                             (("bundle", "schema"), "const"),
                                             (("bundle", "summary", "isoforms"), "minimum")]
