import csv
import io

import pytest
import yaml

from flair_eval import config
from flair_eval.validate import validate_runs, validate_samples


@pytest.mark.parametrize("path,name", [
    ("/x/HG002.fq.gz", "HG002"), ("HG002.rep1.fastq", "HG002.rep1"), ("/a/b/s1.bam", "s1"),
    ("reads.FA.GZ", "reads"), ("s.fasta", "s"), ("weird", "weird"), (".fq.gz", ".fq.gz"),
])
def test_default_sample_name(path, name):
    assert config.default_sample_name(path) == name


def test_samples_csv_columns_and_quoting(sample):
    other = dict(sample, sample="s2", reads="/d/a,b.fq", cage="/t/cage.bed")
    text = config.samples_csv([sample, other])
    rows = list(csv.reader(io.StringIO(text)))
    assert rows[0] == config.SAMPLE_COLUMNS + ["cage"]
    assert rows[2][1] == "/d/a,b.fq"
    assert '"/d/a,b.fq"' in text
    assert config.samples_csv([sample]).splitlines()[0] == ",".join(config.SAMPLE_COLUMNS)


def test_samples_round_trip(tmp_path, sample):
    p = tmp_path / "samples.csv"
    p.write_text(config.samples_csv([sample]))
    header, rows = config.read_samples(p)
    assert header == config.SAMPLE_COLUMNS
    assert rows == [sample]  # empty cells dropped
    assert validate_samples(rows, header) == []


def test_read_samples_bom_blank_lines_and_spaces(tmp_path):
    p = tmp_path / "s.csv"
    p.write_text("﻿sample, reads ,genome,gtf\n s1 ,r.fq,g.fa,a.gtf\n,,,\n", encoding="utf-8")
    header, rows = config.read_samples(p)
    assert header == ["sample", "reads", "genome", "gtf"]
    assert rows == [{"sample": "s1", "reads": "r.fq", "genome": "g.fa", "gtf": "a.gtf"}]


def test_runs_yaml_parses_back_and_validates(two_runs):
    runs = two_runs + [{"id": "ted", "flair": {"env": "/envs/ted"},
                        "transcriptome_args": "--ted --ted_chemistry pacbio_cdna",
                        "label": 'TED "local"', "tags": {"branch": "ted-port", "n": "3"},
                        "quantify": False}]
    text = config.runs_yaml(runs)
    doc = yaml.safe_load(text)
    assert validate_runs(doc) == []
    got = doc["runs"][2]
    assert got["label"] == 'TED "local"'
    assert got["tags"] == {"branch": "ted-port", "n": "3"}
    assert got["quantify"] is False and got["productivity"] is True and got["own_align"] is False
    assert doc["runs"][0]["flair"] == {"release": "3.0.1"}
    assert isinstance(doc["runs"][0]["flair"]["release"], str)


def test_runs_yaml_quotes_number_like_versions():
    doc = yaml.safe_load(config.runs_yaml([{"id": "a", "flair": {"release": "3.0"}}]))
    assert doc["runs"][0]["flair"]["release"] == "3.0"


def test_read_runs_errors(tmp_path):
    p = tmp_path / "r.yaml"
    p.write_text("runs: [\n")
    with pytest.raises(config.YamlError):
        config.read_runs(p)
    p.write_text("")
    assert config.read_runs(p) == {}
