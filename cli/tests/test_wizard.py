import pytest
import yaml

from flair_eval import config
from flair_eval.validate import validate_runs, validate_samples
from flair_eval.wizard import Abort, Wizard, write_configs


def scripted(answers):
    it = iter(answers)
    prompts, said = [], []

    def ask(p):
        prompts.append(p)
        try:
            return next(it)
        except StopIteration:
            raise EOFError
    return Wizard(ask=ask, say=said.append), prompts, said


def test_two_samples_two_versions_plus_settings_run():
    w, prompts, said = scripted([
        "n",                                                                         # test data?
        "/data/HG002.fq.gz", "", "/ref/hg38.fa", "/ref/gencode.gtf", "pacbio", "",  # sample 1
        "y",
        "/data/s2.bam", "", "", "", "ont", "dRNA",                                  # sample 2
        "n",
        "3.0.1,dev", "", "2",                                                       # runs
        "y", "", "4", "", "Strict 3.0.1",                                           # extra run
        "n",
    ])
    samples, runs = w.run()
    assert samples == [
        {"sample": "HG002", "reads": "/data/HG002.fq.gz", "genome": "/ref/hg38.fa",
         "gtf": "/ref/gencode.gtf", "platform": "pacbio"},
        {"sample": "s2", "bam": "/data/s2.bam", "genome": "/ref/hg38.fa",
         "gtf": "/ref/gencode.gtf", "platform": "ont", "library": "dRNA"},
    ]
    assert [r["id"] for r in runs] == ["flair-3.0.1", "flair-dev", "flair-dev-stringent_splice"]
    assert runs[0]["transcriptome_args"] == ""
    assert runs[1]["transcriptome_args"] == "--stringent"
    assert runs[2]["flair"]["ref"] == "dev"  # defaults to the previous run's FLAIR
    assert runs[2]["label"] == "Strict 3.0.1"
    assert validate_samples(samples) == [] and validate_runs({"runs": runs}) == []


def test_reasks_on_bad_answers():
    w, prompts, said = scripted([
        "n", "/d/r.fq", "bad name", "ok", "", "g.fa", "a.gtf", "nanopore", "", "", "n",
        "master", "3.0.1", "--support 10", "n",
    ])
    samples, runs = w.run()
    assert samples[0]["sample"] == "ok"
    assert "platform" not in samples[0]
    assert runs == [{"id": "flair-3.0.1", "flair": {"release": "3.0.1"},
                     "transcriptome_args": "--support 10"}]
    text = "\n".join(said)
    assert "Not a valid name" in text and "This one is required" in text
    assert "Use one of pacbio/ont" in text and "git:master" in text


def test_duplicate_run_id_rejected():
    w, _, said = scripted(["n", "r.fq", "", "g", "a", "", "", "n", "3.0.1", "", "y", "3.0.1", "1",
                           "flair-3.0.1", "second", "", "n"])
    _, runs = w.run()
    assert [r["id"] for r in runs] == ["flair-3.0.1", "second"]
    assert any("already used" in s for s in said)


def test_eof_aborts():
    w, _, _ = scripted(["n", "r.fq"])
    with pytest.raises(Abort):
        w.run()


def test_write_configs(tmp_path, sample, two_runs):
    sp, rp = write_configs(tmp_path / "cfg", [sample], two_runs)
    assert config.read_samples(sp)[1] == [sample]
    assert yaml.safe_load(rp.read_text())["runs"][1]["id"] == "flair-dev"
    with pytest.raises(FileExistsError):
        write_configs(tmp_path / "cfg", [sample], two_runs)
    write_configs(tmp_path / "cfg", [sample], two_runs, force=True)


def test_test_data_shortcut():
    w, _, _ = scripted(["y"])
    assert w.run() is None


def test_platform_and_library_defaults_from_file_name():
    from flair_eval.wizard import guess_library, guess_platform
    assert guess_platform("/d/WTC11_ONT_cDNA.fq.gz") == "ont"
    assert guess_platform("/d/sample.hifi_reads.bam") == "pacbio"
    assert guess_platform("/d/montana.fq") == ""
    assert guess_library("/d/WTC11_ONT_cDNA.fq.gz") == "cDNA"
    assert guess_library("/d/x_dRNA.fastq") == "dRNA"
    w, _, _ = scripted(["n", "/d/a549_nanopore_drna.fq.gz", "", "g", "a", "", "", "n",
                        "3.0.1", "", "n"])
    samples, _ = w.run()
    assert samples[0]["platform"] == "ont" and samples[0]["library"] == "dRNA"
    w, _, _ = scripted(["n", "/d/a549_nanopore.fq.gz", "", "g", "a", "-", "", "n", "3.0.1", "",
                        "n"])
    assert "platform" not in w.run()[0][0]
