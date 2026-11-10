"""Unit tests for bin/flair_adapter.py against fake `flair` executables, one per CLI generation.

Generations (see docs/flair_versions.md):
  v2    FLAIR 2.2.0    align, correct + collapse, reads quantify, predictProductivity
  v301  FLAIR 3.0.1    align, transcriptome (--parallelmode), reads quantify, predictProductivity
  ted   3.0.0+dev TED  align, transcriptome --sample_name, reads quantify + --isoform_bed, 31-col BED
  dev   3.0.0+dev      no align, transcriptome --sample_name, BAM quantify, inline CDS, 25-col BED
"""
import json
import os
import stat
import subprocess
import sys
import textwrap

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ADAPTER = os.path.join(ROOT, "bin", "flair_adapter.py")

FAKE_FLAIR = r'''
import json, os, sys
GEN = os.environ["FAKE_FLAIR_GEN"]
VERSION = {"v2": "2.2.0", "v301": "3.0.1", "ted": "3.0.0+dev", "dev": "3.0.0+dev"}[GEN]
argv = sys.argv[1:]
calls = os.environ.get("FAKE_FLAIR_CALLS")
if calls:
    with open(calls, "a") as fh:
        fh.write(json.dumps(argv) + "\n")

def opt(*names, default=None):
    for i, a in enumerate(argv):
        if a in names and i + 1 < len(argv):
            return argv[i + 1]
    return default

if argv == ["--version"]:
    print("FLAIR " + VERSION); sys.exit(0)
if argv == ["--help"]:
    if GEN == "v2":
        print("usage: flair module --help\nmodule: align, correct, collapse, quantify\n\nVersion 2.2.0")
    elif GEN == "dev":
        print("usage: flair [-h] [--version] <subcommand> ...\n\noptions:\n  -h, --help  help\n\n"
              "subcommands:\n  <subcommand>\n    transcriptome       Build a transcriptome\n"
              "    quantify            Quantify isoform expression\n    combine             Combine")
    elif GEN == "ted":
        print("usage: flair [-h] [--version]\n             {align,transcriptome,quantify,combine} ...")
    else:
        print("usage: flair [-h] [--version]\n             {align,correct,transcriptome,collapse,quantify,combine} ...")
    sys.exit(0)

SUBS = {"v2": ["align", "correct", "collapse", "quantify"],
        "v301": ["align", "correct", "transcriptome", "collapse", "quantify"],
        "ted": ["align", "transcriptome", "quantify"],
        "dev": ["transcriptome", "quantify"]}[GEN]
sub = argv[0]
if sub not in SUBS:
    sys.stderr.write("flair: invalid choice: %r\n" % sub); sys.exit(2)

HELP = {
 ("v2", "collapse"): "usage: flair collapse -g genome.fa -q <query.bed> -r <reads.fq> [options]\n\n"
     "  -f GTF, --gtf GTF\n  --generate_map\n  --temp_dir TEMP_DIR\n  -t T, --threads T\n  --stringent\n\n"
     "required named arguments:\n  -q QUERY, --query QUERY\n  -g GENOME, --genome GENOME\n  -r READS, --reads READS\n",
 ("v2", "correct"): "usage: flair correct -q query.bed12 -g genome.fa [options]\n\n  -f GTF, --gtf GTF\n  --nvrna\n"
     "  -t T, --threads T\n\nrequired named arguments:\n  -q QUERY, --query QUERY\n  -g GENOME, --genome GENOME\n",
 ("v301", "transcriptome"): "usage: transcriptome [-h] [-b GENOMEALIGNEDBAM] [-g GENOME] [-o OUTPUT]\n"
     "                     [-t THREADS] [-f GTF] [--stringent] [--parallelmode PARALLELMODE]\n\n"
     "options:\n  -b, --genomealignedbam X\n  --stringent\n  --parallelmode PARALLELMODE\n",
 ("ted", "transcriptome"): "usage: transcriptome [-h] -b GENOME_ALIGNED_BAM -g GENOME --sample_name\n"
     "                     SAMPLE_NAME [-o OUTPUT] [-f ANNOT_GTF] [--parallel_mode PARALLEL_MODE]\n"
     "                     [--generate_map] [--ted]\n\n  --stringent_x\n\nrequired named arguments:\n"
     "  -b GENOME_ALIGNED_BAM, --genome_aligned_bam GENOME_ALIGNED_BAM\n  -g GENOME, --genome GENOME\n"
     "  --sample_name SAMPLE_NAME\n",
 ("dev", "transcriptome"): "usage: flair transcriptome [-h] -b GENOME_ALIGNED_BAM -g GENOME --sample_name\n"
     "                           SAMPLE_NAME [-o OUTPUT] [-f ANNOT_GTF] [--parallel_mode PARALLEL_MODE]\n"
     "                           [--temp_dir TEMP_DIR] [--generate_map] [--no_stringent] [--directRNA]\n\n"
     "required named arguments:\n  -b GENOME_ALIGNED_BAM, --genome_aligned_bam GENOME_ALIGNED_BAM\n"
     "  -g GENOME, --genome GENOME\n  --sample_name SAMPLE_NAME\n",
 ("dev", "quantify"): "usage: flair quantify [-h] --manifest MANIFEST -g GENOME --isoform_bed\n"
     "                      ISOFORM_BED [-o OUTPUT] [-t THREADS] [--temp_dir TEMP_DIR]\n",
 ("dev-long", "quantify"): "usage: flair quantify [-h] --manifest MANIFEST --genome GENOME --isoform_bed\n"
     "                      ISOFORM_BED [--output OUTPUT] [--threads THREADS] [--temp_dir TEMP_DIR]\n",
 ("ted", "quantify"): "usage: quantify [-h] -r R -i I [-o O] [-t T] [--temp_dir TEMP_DIR]\n"
     "                [--isoform_bed ISOFORMS]\n\nrequired named arguments:\n  -r R, --reads_manifest R\n"
     "  -i I, --isoforms I\n  --isoform_bed ISOFORMS, --isoformbed ISOFORMS\n",
}
if "--help" in argv:
    text = None
    if GEN == "dev" and os.environ.get("FAKE_DEV_LONG_ONLY"):
        text = HELP.get(("dev-long", sub))
    if text is None:
        text = HELP.get((GEN, sub))
    if text and GEN == "dev" and sub == "transcriptome" and os.environ.get("FAKE_DEV_NO_PARALLEL"):
        text = text.replace("[--parallel_mode PARALLEL_MODE]", "[--region_bed REGION_BED]")
    if text is None:
        if sub == "quantify":
            text = ("usage: flair quantify -r reads_manifest.tsv -i isoforms.fa [options]\n\n"
                    "  --temp_dir TEMP_DIR\n  --isoform_bed ISOFORMS\n\nrequired named arguments:\n"
                    "  -r R, --reads_manifest R\n  -i I, --isoforms I\n")
        elif sub == "align":
            text = "usage: flair_align -g genome.fa -r reads [options]\n\n  --nvrna\n  -t T, --threads T\n"
        else:
            text = "usage: flair %s [options]\n" % sub
    print(text); sys.exit(0)

fail = os.environ.get("FAKE_FLAIR_FAIL")
if fail == sub:
    sys.stderr.write("Traceback: fake failure in %s\n" % sub); sys.exit(3)

def w(path, text):
    with open(path, "w") as fh:
        fh.write(text)

# three isoforms: an annotated one, a novel one with a coordinate gene, and (2.x/3.0.1) one
# whose gene the GTF cannot give.
ISOS = [("chr1", 100, 900, "+", [(100, 300), (500, 900)]),
        ("chr1", 2000, 2600, "-", [(2000, 2200), (2400, 2600)]),
        ("chr2", 50, 450, "+", [(50, 450)])]

def bed12(c, s, e, st, ex, name, score):
    sizes = ",".join(str(b - a) for a, b in ex) + ","
    starts = ",".join(str(a - s) for a, b in ex) + ","
    return [c, str(s), str(e), name, str(score), st, str(s), str(e), "0,0,0", str(len(ex)), sizes, starts]

if sub == "align":
    o = opt("-o", "--output")
    w(o + ".bam", "readA\nreadB\nreadC\n"); w(o + ".bam.bai", ""); w(o + ".bed", "")
    sys.exit(0)

if sub == "correct":
    o = opt("-o", "--output")
    w(o + "_all_corrected.bed", open(opt("-q", "--query")).read())
    sys.exit(0)

if sub in ("transcriptome", "collapse"):
    o = opt("-o", "--output")
    if GEN in ("v2", "v301"):
        names = ["ENST1.1_ENSG1.1", "readB_chr1:2000", "readC_chr2:50"]
        genes = ["ENSG1.1", "chr1:2000", None]
        tids = ["ENST1.1", "readB", "readC"]
    else:
        names = ["FLT00000001", "FLT00000002", "FLT00000003"]
        genes = ["FLG00000001", "FLG00000002", ""]
        tids = names
    rows, gtf, fa, mp = [], [], [], []
    for (c, s, e, st, ex), n, g, t in zip(ISOS, names, genes, tids):
        r = bed12(c, s, e, st, ex, n, 1500)
        if GEN == "dev":
            r[6], r[7] = str(s + 10), str(e - 10)
            r += [g, "ENST1.1", "ENSG1.1,", "5", "0.5", "PRO" if n.endswith("1") else "NGO", "basic",
                  "", "", "s1,", "", "", "FLP1"]
        elif GEN == "ted":
            r += [g, "ENST1.1", "ENSG1.1,", "5", "0.5", "PTC", "basic", "", "", "s1,"] + [""] * 9
        rows.append("\t".join(r))
        if g is not None:
            gtf.append('%s\tFLAIR\ttranscript\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s";' % (c, s + 1, e, st, g, t))
        fa.append(">%s some description\nACGTACGT\nACGT" % n)
        mp.append("%s\tr%s1,r%s2" % (n, n[-1], n[-1]))
    rows.append(rows[0])   # FLAIR occasionally repeats a name; the adapter keeps the first
    w(o + ".isoforms.bed", "\n".join(rows) + "\n")
    w(o + ".isoforms.gtf", "\n".join(gtf) + "\n")
    w(o + ".isoforms.fa", "\n".join(fa) + "\n>orphan\nAAAA\n")
    if GEN == "dev":
        w(o + ".isoform.read.map.tsv", "name\treads\n" + "\n".join(mp) + "\n")
        w(o + ".isoform.counts.tsv", "isoform_id\tfull_length_reads\ttotal_reads\n")
    elif GEN != "v2" or "--generate_map" in argv:
        w(o + ".isoform.read.map.txt", "\n".join(mp) + "\n")
    sys.exit(0)

if sub == "quantify":
    o = opt("-o", "--output")
    if GEN == "dev":
        if os.environ.get("FAKE_DEV_LONG_ONLY") and "-g" in argv:
            sys.stderr.write("flair quantify: error: unrecognized arguments: -g\n"); sys.exit(2)
        assert "-g" in argv or "--genome" in argv
        for need in ("--manifest", "--isoform_bed"):
            assert need in argv, need
        nb = [l.split("\t") for l in open(opt("--isoform_bed")).read().splitlines()]
        assert all(len(f) >= 25 for f in nb), "dev quantify needs FLAIR's 25-column BED"
        w(o + ".counts.tsv", "gene_id\tisoform_id\ts1\nFLG00000001\tFLT00000001\t10\nFLG00000002\tFLT00000002\t30\n")
    else:
        man = open(opt("-r", "--reads_manifest")).read().split("\t")
        assert "_" not in man[0], "FLAIR 3.0.1 rejects '_' in sample ids"
        if GEN == "ted":
            assert "--isoform_bed" in argv
        names = [l[1:].split()[0] for l in open(opt("-i", "--isoforms")) if l.startswith(">")]
        # omits the last isoform (zero reads); 2.x/3.0.1 write floats
        body = "".join("%s\t%s\n" % (n, "3.0" if i == 0 else "1.0") for i, n in enumerate(names[:-1]))
        w(o + ".counts.tsv", "ids\ts1_c1_b1\n" + body)
    sys.exit(0)
sys.exit(0)
'''

FAKE_PP = r'''
import os, sys
GEN = os.environ["FAKE_FLAIR_GEN"]
argv = sys.argv[1:]
def opt(n):
    return argv[argv.index(n) + 1]
if "--help" in argv:
    if GEN == "ted":
        print("usage: predictProductivity [-h] -i INPUT_ISOFORMS -g GTF -f GENOME_FASTA -o OUTPUT [--is_flair_bed]")
    else:
        print("usage: predictProductivity -i isoforms.bed -f genome.fa -g annotations.gtf -o outputprefix\n"
              "  --append_column\n  --firstTIS\n  --longestORF")
    sys.exit(0)
if "--longestORF" not in argv and "--firstTIS" not in argv:
    sys.stderr.write("error: one of the arguments --firstTIS --longestORF is required\n"); sys.exit(2)
assert ".bed" not in os.getcwd()
with open(opt("-o") + ".bed", "w") as out:
    for i, line in enumerate(open(opt("-i"))):
        f = line.rstrip("\n").split("\t")
        f[6], f[7] = str(int(f[1]) + 5), str(int(f[2]) - 5)
        out.write("\t".join(f + [["PRO", "PTC", "NST"][i % 3]]) + "\n")
open(opt("-o") + ".info.tsv", "w").write("#isoname\n")
'''

FAKE_SAMTOOLS = r'''
import shutil, sys
a = sys.argv[1:]
if a[0] == "faidx":
    open(a[1] + ".fai", "w").write("chr1\t4\t6\t4\t5\n")
elif a[0] == "index":
    open(a[1] + ".bai", "w").close()
elif a[0] == "fastq":
    for name in open(a[-1]).read().split():
        print("@%s\nACGT\n+\nIIII" % name)
elif a[0] == "view" and "-b" in a:
    shutil.copyfile(a[-1], a[a.index("-o") + 1])
elif a[0] == "view":
    for name in open(a[-1]).read().split():
        if name != "readC":            # readC lies outside the regions
            print("%s\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGT\t*" % name)
'''

FAKE_BAM2BED = r'''
import sys
bam = sys.argv[sys.argv.index("-i") + 1]
for i, name in enumerate(open(bam).read().split()):
    print("chr1\t%d\t%d\t%s\t60\t+\t%d\t%d\t0\t1\t100,\t0," % (i * 1000, i * 1000 + 100, name, i * 1000, i * 1000 + 100))
'''


def write_exe(path, body):
    with open(path, "w") as fh:
        fh.write("#!%s\n%s" % (sys.executable, body))
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


@pytest.fixture
def fake(tmp_path):
    """Make a fake FLAIR env for a generation; returns a runner for the adapter CLI."""
    inputs = tmp_path / "in"
    inputs.mkdir()
    (inputs / "reads.fa").write_text(">readA x\nACGT\n>readB\nACGT\n>readC\nACGT\n")
    (inputs / "reads.fq").write_text("@readA\nACGT\n+\nIIII\n@readC\nACGT\n+\nIIII\n")
    (inputs / "genome.fa").write_text(">chr1\nACGT\n")
    (inputs / "genome.fa.fai").write_text("chr1\t4\t6\t4\t5\n")    # core stages it (A23)
    (inputs / "annot.gtf").write_text('chr1\tx\texon\t1\t10\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
    (inputs / "shared.bam").write_text("readA\nreadB\nreadC\n")
    (inputs / "shared.bam.bai").write_text("")
    (inputs / "regions.bed").write_text("chr1\t0\t5000\n")

    def make(gen, with_pp=True, fail=None, extra_env=None):
        b = tmp_path / ("env_" + gen) / "bin"
        b.mkdir(parents=True, exist_ok=True)
        write_exe(str(b / "flair"), FAKE_FLAIR)
        write_exe(str(b / "samtools"), FAKE_SAMTOOLS)
        if with_pp and gen != "dev":
            write_exe(str(b / "predictProductivity"), FAKE_PP)
        os.symlink(sys.executable, str(b / "python")) if not (b / "python").exists() else None
        pkg = tmp_path / ("env_" + gen) / "pylib" / "flair"
        pkg.mkdir(parents=True, exist_ok=True)
        (pkg / "__init__.py").write_text("")
        (pkg / "bam2Bed12.py").write_text(FAKE_BAM2BED)
        calls = tmp_path / ("calls_" + gen + ".jsonl")
        env = dict(os.environ)
        env.update({"FAKE_FLAIR_GEN": gen, "FAKE_FLAIR_CALLS": str(calls),
                    "PATH": str(b) + os.pathsep + env["PATH"],
                    "PYTHONPATH": str(tmp_path / ("env_" + gen) / "pylib")})
        env.pop("FLAIR_EVAL_FLAIR", None)
        env.pop("FAKE_DEV_LONG_ONLY", None)
        env.pop("FAKE_DEV_NO_PARALLEL", None)
        env.update(extra_env or {})
        if fail:
            env["FAKE_FLAIR_FAIL"] = fail
        else:
            env.pop("FAKE_FLAIR_FAIL", None)

        def run(*args, check=None):
            p = subprocess.run([sys.executable, ADAPTER] + [str(a) for a in args], env=env,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            if check is not None:
                assert p.returncode == check, p.stderr
            return p

        def flair_calls():
            if not calls.exists():
                return []
            out = [json.loads(l) for l in calls.read_text().splitlines()]
            return [c for c in out if "--help" not in c and c != ["--version"]]

        run.calls = flair_calls
        run.inputs = inputs
        run.tmp = tmp_path
        return run
    return make


def assemble(run, out, *extra, check=None):
    i = run.inputs
    return run("assemble", "--bam", i / "shared.bam", "--reads", i / "reads.fa", "--genome", i / "genome.fa",
               "--gtf", i / "annot.gtf", "--args", "", "--threads", "2", "--out", out, *extra, check=check)


def argval(call, *names):
    """Value after whichever of `names` (short or long form) the command used."""
    for n in names:
        if n in call:
            return call[call.index(n) + 1]
    raise AssertionError("none of %s in %s" % (names, call))


def read_tsv(path):
    lines = open(path).read().splitlines()
    return [l.split("\t") for l in lines]


def check_assembly(out):
    bed = read_tsv(out / "isoforms.bed")
    assert all(len(f) == 12 for f in bed)
    names = [f[3] for f in bed]
    assert len(names) == len(set(names)) == 3
    assert all(0 <= int(f[4]) <= 1000 for f in bed)
    assert not any(f[10].endswith(",") for f in bed)
    gtf = read_tsv(out / "isoforms.gtf")
    tx = [f for f in gtf if f[2] == "transcript"]
    assert len(tx) == 3
    for f in gtf:
        assert 'gene_id "' in f[8] and 'transcript_id "' in f[8]
    gtf_names = {f[8].split('transcript_id "')[1].split('"')[0] for f in tx}
    assert gtf_names == set(names)
    fa = open(out / "isoforms.fa").read()
    assert [l[1:] for l in fa.splitlines() if l.startswith(">")] == names
    assert "orphan" not in fa
    log = json.load(open(out / "adapter_log.json"))
    assert log["status"] == "ok" and log["step"] == "assemble"
    assert log["commands"] and all(c["exit_code"] == 0 for c in log["commands"])
    assert log["outputs"]["isoforms"] == 3
    assert (out / "isoforms.flair.bed").exists()
    assert not (out / ".work").exists()
    return bed, gtf, log


def genes_of(gtf):
    return {f[8].split('transcript_id "')[1].split('"')[0]: f[8].split('gene_id "')[1].split('"')[0]
            for f in gtf if f[2] == "transcript"}


# --------------------------------------------------------------------------------------------
def test_help_and_bad_input(fake):
    run = fake("v301")
    assert run("--help", check=0).stdout.count("assemble") >= 1
    for step in ("align", "assemble", "quantify", "productivity"):
        assert run(step, "--help", check=0).stdout.startswith("usage")
    assert run(check=None).returncode != 0
    p = run("assemble", "--bam", "/nope.bam", "--genome", "/nope.fa", "--gtf", "/nope.gtf",
            "--out", run.tmp / "bad")
    assert p.returncode == 2 and "not found" in p.stderr
    assert json.load(open(run.tmp / "bad" / "adapter_log.json"))["status"] == "failed"


@pytest.mark.parametrize("gen,cli,align,qmode", [
    ("v2", "correct+collapse", True, "reads_manifest"),
    ("v301", "transcriptome", True, "reads_manifest"),
    ("ted", "transcriptome", True, "reads_manifest"),
    ("dev", "transcriptome", False, "bam_manifest"),
])
def test_detect(fake, gen, cli, align, qmode):
    prof = json.loads(fake(gen)("detect", "--json", check=0).stdout)
    assert prof["assemble_cli"] == cli
    assert prof["align"] is align
    assert prof["quantify_mode"] == qmode
    if gen in ("ted", "dev"):
        assert prof["transcriptome"]["needs_sample_name"]
        assert prof["transcriptome"]["parallel_flag"] == "--parallel_mode"
        assert prof["quantify"]["needs_isoform_bed"]
    if gen == "v301":
        assert prof["transcriptome"]["parallel_flag"] == "--parallelmode"
        assert not prof["transcriptome"]["needs_sample_name"]
        assert not prof["quantify"]["needs_isoform_bed"]
    assert bool(prof["predict_productivity"]) is (gen != "dev")


def test_align(fake):
    run = fake("v301")
    out = run.tmp / "align"
    run("align", "--reads", run.inputs / "reads.fa", "--genome", run.inputs / "genome.fa",
        "--threads", "3", "--out", out, check=0)
    assert (out / "aligned.bam").exists() and (out / "aligned.bam.bai").exists()
    call = run.calls()[-1]
    assert call[0] == "align" and argval(call, "-t", "--threads") == "3"
    log = json.load(open(out / "adapter_log.json"))
    assert log["flair_version"] == "3.0.1" and log["commands"][0]["argv"][1] == "align"


def test_align_skipped_on_dev(fake):
    run = fake("dev")
    out = run.tmp / "align"
    run("align", "--reads", run.inputs / "reads.fa", "--genome", run.inputs / "genome.fa",
        "--out", out, check=0)
    sk = json.load(open(out / "skipped.json"))
    assert sk["step"] == "align" and "align" in sk["reason"]
    assert json.load(open(out / "adapter_log.json"))["status"] == "skipped"
    assert not (out / "aligned.bam").exists()


def test_assemble_v301(fake):
    run = fake("v301")
    out = run.tmp / "asm"
    assemble(run, out, check=0)
    bed, gtf, log = check_assembly(out)
    g = genes_of(gtf)
    assert g["ENST1.1_ENSG1.1"] == "ENSG1.1"
    assert g["readB_chr1:2000"] == "chr1:2000"
    assert g["readC_chr2:50"] == "chr2:50"          # from the name suffix when the GTF has none
    tx = [c for c in run.calls() if c[0] == "transcriptome"][0]
    assert tx[tx.index("--parallelmode") + 1] == "byregion"
    assert "--sample_name" not in tx
    rm = read_tsv(out / "read_map.tsv")
    assert rm[0] == ["read_id", "isoform_id"] and len(rm) == 7
    assert any("duplicate" in w for w in log["warnings"])


def test_assemble_v2_uses_correct_collapse(fake):
    run = fake("v2")
    out = run.tmp / "asm"
    assemble(run, out, check=0)
    _, gtf, log = check_assembly(out)
    subs = [c[0] for c in run.calls()]
    assert subs.index("correct") < subs.index("collapse")
    assert "transcriptome" not in subs
    col = [c for c in run.calls() if c[0] == "collapse"][0]
    assert "--generate_map" in col and argval(col, "-r", "--reads").endswith("reads.fa")
    labels = [c["label"] for c in log["commands"]]
    assert labels[:3] == ["bam2Bed12", "flair_correct", "flair_collapse"]


def test_assemble_v2_without_reads_extracts_them_from_bam(fake):
    """A6: --reads is optional; 2.x collapse then gets the BAM's primary reads."""
    run = fake("v2")
    i = run.inputs
    out = run.tmp / "asm"
    run("assemble", "--bam", i / "shared.bam", "--genome", i / "genome.fa", "--gtf", i / "annot.gtf",
        "--out", out, check=0)
    check_assembly(out)
    col = [c for c in run.calls() if c[0] == "collapse"][0]
    assert argval(col, "-r", "--reads").endswith("reads_from_bam.fq")
    log = json.load(open(out / "adapter_log.json"))
    fq = [c for c in log["commands"] if c["label"] == "samtools_fastq"][0]
    assert fq["argv"][1:4] == ["fastq", "-F", "0x904"]
    assert any("extracted from the BAM" in n for n in log["notes"])


def test_assemble_v301_without_reads(fake):
    run = fake("v301")
    i = run.inputs
    out = run.tmp / "asm"
    run("assemble", "--bam", i / "shared.bam", "--genome", i / "genome.fa", "--gtf", i / "annot.gtf",
        "--out", out, check=0)
    check_assembly(out)


def test_assemble_v2_regions_subsets_reads(fake):
    run = fake("v2")
    out = run.tmp / "asm"
    assemble(run, out, "--regions", run.inputs / "regions.bed", "--keep-work", check=0)
    col = [c for c in run.calls() if c[0] == "collapse"][0]
    sub_reads = argval(col, "-r", "--reads")
    assert "regions.reads" in sub_reads
    assert open(sub_reads).read().count(">") == 2     # readC is outside the regions
    labels = [c["label"] for c in json.load(open(out / "adapter_log.json"))["commands"]]
    assert labels[0] == "samtools_view_regions"


@pytest.mark.parametrize("gen", ["ted", "dev"])
def test_assemble_sample_name_and_genes(fake, gen):
    run = fake(gen)
    out = run.tmp / "asm"
    assemble(run, out, "--regions", run.inputs / "regions.bed", check=0)
    _, gtf, log = check_assembly(out)
    g = genes_of(gtf)
    assert g["FLT00000001"] == "FLG00000001"
    assert g["FLT00000003"] == "novel_gene_1"        # FLAIR gave no gene
    tx = [c for c in run.calls() if c[0] == "transcriptome"][0]
    assert tx[tx.index("--sample_name") + 1] == "flaireval"
    assert tx[tx.index("--parallel_mode") + 1] == "byregion"
    assert "--generate_map" in tx
    assert argval(tx, "-b", "--genome_aligned_bam", "--genomealignedbam").endswith("regions.bam")
    if gen == "dev":
        assert "--temp_dir" in tx
    ncols = len(read_tsv(out / "isoforms.flair.bed")[0])
    assert ncols == (25 if gen == "dev" else 31)
    assert len(read_tsv(out / "read_map.tsv")) == 7


def test_assemble_user_args(fake):
    run = fake("v301")
    out = run.tmp / "asm"
    i = run.inputs
    run("assemble", "--bam", i / "shared.bam", "--genome", i / "genome.fa", "--gtf", i / "annot.gtf",
        "--args", "--stringent --parallelmode bychrom --no_such_flag 3", "--out", out, check=0)
    tx = [c for c in run.calls() if c[0] == "transcriptome"][0]
    assert tx[-5:] == ["--stringent", "--parallelmode", "bychrom", "--no_such_flag", "3"]
    assert tx.count("--parallelmode") == 1          # the user's choice wins
    log = json.load(open(out / "adapter_log.json"))
    assert any("--no_such_flag" in w for w in log["warnings"])


def test_assemble_failure_exit_code(fake):
    run = fake("v301", fail="transcriptome")
    out = run.tmp / "asm"
    p = assemble(run, out)
    assert p.returncode == 3
    log = json.load(open(out / "adapter_log.json"))
    assert log["status"] == "failed" and log["commands"][-1]["exit_code"] == 3
    assert (out / ".work").exists()                  # kept for debugging
    assert not (out / "skipped.json").exists()


@pytest.mark.parametrize("gen", ["v2", "v301", "ted"])
def test_quantify_reads(fake, gen):
    run = fake(gen)
    asm = run.tmp / "asm"
    assemble(run, asm, check=0)
    out = run.tmp / "q"
    run("quantify", "--isoforms-fa", asm / "isoforms.fa", "--reads", run.inputs / "reads.fa",
        "--sample", "my_sample", "--threads", "2", "--out", out, check=0)
    rows = read_tsv(out / "counts.tsv")
    assert rows[0] == ["isoform_id", "count", "cpm"]
    names = [l[1:] for l in open(asm / "isoforms.fa").read().splitlines() if l.startswith(">")]
    assert [r[0] for r in rows[1:]] == names
    counts = [float(r[1]) for r in rows[1:]]
    assert counts[-1] == 0                           # zero-read row filled in
    assert abs(sum(float(r[2]) for r in rows[1:]) - 1e6) < 1
    assert float(rows[1][2]) == pytest.approx(3 / 4 * 1e6, rel=1e-6)
    q = [c for c in run.calls() if c[0] == "quantify"][0]
    if gen == "ted":
        assert q[q.index("--isoform_bed") + 1].endswith("isoforms.flair.bed")
    log = json.load(open(out / "adapter_log.json"))
    assert log["status"] == "ok" and "my_sample" in json.dumps(log)


def test_quantify_dev_needs_bam(fake):
    run = fake("dev")
    asm = run.tmp / "asm"
    assemble(run, asm, check=0)
    out = run.tmp / "q"
    run("quantify", "--isoforms-fa", asm / "isoforms.fa", "--reads", run.inputs / "reads.fa",
        "--sample", "s", "--out", out, check=0)
    sk = json.load(open(out / "skipped.json"))
    assert "--bam" in sk["reason"]
    out2 = run.tmp / "q2"
    run("quantify", "--isoforms-fa", asm / "isoforms.fa", "--reads", run.inputs / "reads.fa",
        "--sample", "s", "--bam", run.inputs / "shared.bam", "--genome", run.inputs / "genome.fa",
        "--isoforms-bed", asm / "isoforms.bed", "--out", out2, check=0)
    rows = read_tsv(out2 / "counts.tsv")
    assert rows[1:] == [["FLT00000001", "10", "250000.0000"], ["FLT00000002", "30", "750000.0000"],
                        ["FLT00000003", "0", "0.0000"]]
    assert not (out2 / "skipped.json").exists()


def test_quantify_flair_bed_flag(fake):
    """A8: FLAIR's own BED comes through --flair-bed; staging symlinks are not followed."""
    run = fake("ted")
    asm = run.tmp / "asm"
    assemble(run, asm, check=0)
    stage = run.tmp / "stage"
    stage.mkdir()
    os.symlink(str(asm / "isoforms.fa"), str(stage / "isoforms.fa"))
    os.symlink(str(asm / "isoforms.flair.bed"), str(stage / "native.bed"))
    out = run.tmp / "q"
    run("quantify", "--isoforms-fa", stage / "isoforms.fa", "--reads", run.inputs / "reads.fa",
        "--sample", "s", "--out", out, check=0)
    assert "--flair-bed" in json.load(open(out / "skipped.json"))["reason"]
    out2 = run.tmp / "q2"
    run("quantify", "--isoforms-fa", stage / "isoforms.fa", "--reads", run.inputs / "reads.fa",
        "--sample", "s", "--flair-bed", stage / "native.bed", "--out", out2, check=0)
    assert (out2 / "counts.tsv").exists()
    q = [c for c in run.calls() if c[0] == "quantify"][-1]
    assert q[q.index("--isoform_bed") + 1].endswith("native.bed")
    p = run("quantify", "--isoforms-fa", stage / "isoforms.fa", "--reads", run.inputs / "reads.fa",
            "--sample", "s", "--flair-bed", run.tmp / "nope.bed", "--out", run.tmp / "q3")
    assert p.returncode == 2 and "--flair-bed" in p.stderr and "not found" in p.stderr


def test_quantify_bam_only_sample_skips_reads_versions(fake):
    """A6/A7: without reads, a version that quantifies from reads writes skipped.json."""
    run = fake("v301")
    asm = run.tmp / "asm"
    assemble(run, asm, check=0)
    out = run.tmp / "q"
    run("quantify", "--isoforms-fa", asm / "isoforms.fa", "--sample", "s", "--bam", run.inputs / "shared.bam",
        "--genome", run.inputs / "genome.fa", "--isoforms-bed", asm / "isoforms.bed", "--out", out, check=0)
    assert "BAM-only" in json.load(open(out / "skipped.json"))["reason"]


def test_quantify_dev_from_bam_without_reads(fake):
    """A7: dev counts from the BAM, so a BAM-only sample still gets counts."""
    run = fake("dev")
    asm = run.tmp / "asm"
    assemble(run, asm, check=0)
    out = run.tmp / "q"
    run("quantify", "--isoforms-fa", asm / "isoforms.fa", "--sample", "s", "--bam", run.inputs / "shared.bam",
        "--genome", run.inputs / "genome.fa", "--isoforms-bed", asm / "isoforms.bed",
        "--flair-bed", asm / "isoforms.flair.bed", "--out", out, "--keep-work", check=0)
    assert len(read_tsv(out / "counts.tsv")) == 4
    q = [c for c in run.calls() if c[0] == "quantify"][-1]
    man = open(q[q.index("--manifest") + 1]).read().split("\t")
    assert man[3].strip().endswith("shared.bam")


@pytest.mark.parametrize("gen,command,flag", [("v301", "align", "--nvrna"), ("v2", "correct", "--nvrna"),
                                              ("dev", "transcriptome", "--directRNA")])
def test_library_drna_flags(fake, gen, command, flag):
    run = fake(gen)
    out = run.tmp / "o"
    i = run.inputs
    if command == "align":
        run("align", "--reads", i / "reads.fa", "--genome", i / "genome.fa", "--library", "drna",
            "--platform", "ONT", "--out", out, check=0)
    else:
        assemble(run, out, "--library", "dRNA", "--platform", "ont", check=0)
    call = [c for c in run.calls() if c[0] == command][0]
    assert flag in call
    log = json.load(open(out / "adapter_log.json"))
    assert log["inputs"]["library"] == "dRNA" and log["inputs"]["platform"] == "ont"
    assert log["library_flags"] == [{"command": command, "flag": flag, "library": "dRNA"}]


def test_library_cdna_and_unmapped_versions_add_nothing(fake):
    run = fake("v301")
    out = run.tmp / "o"
    assemble(run, out, "--library", "dRNA", check=0)      # 3.0.1 transcriptome has no dRNA flag
    tx = [c for c in run.calls() if c[0] == "transcriptome"][0]
    assert "--directRNA" not in tx and "--nvrna" not in tx
    log = json.load(open(out / "adapter_log.json"))
    assert "library_flags" not in log and any("--directRNA" in n for n in log["notes"])
    run2 = fake("dev")
    assemble(run2, run.tmp / "o2", "--library", "cDNA", check=0)
    assert "--directRNA" not in [c for c in run2.calls() if c[0] == "transcriptome"][0]


def test_bad_library_value(fake):
    run = fake("v301")
    p = assemble(run, run.tmp / "o", "--library", "polyA")
    assert p.returncode == 2 and "cDNA|dRNA" in p.stderr


def test_failure_message_names_command_and_log(fake):
    run = fake("v301", fail="transcriptome")
    out = run.tmp / "asm"
    p = assemble(run, out)
    assert p.returncode == 3
    err = p.stderr
    assert "the assemble step failed (FLAIR 3.0.1)" in err
    assert "What failed: `flair transcriptome` exited with code 3" in err
    assert "Command:" in err and " transcriptome --genomealignedbam " in err
    assert "Log:         %s" % (out / "logs" / "01_flair_transcriptome.log") in err
    assert "| Traceback: fake failure in transcriptome" in err
    log = json.load(open(out / "adapter_log.json"))
    assert log["failed_command"]["exit_code"] == 3


def test_missing_output_message(fake, tmp_path):
    run = fake("v301")
    b = tmp_path / "env_v301" / "bin"
    write_exe(str(b / "flair"), FAKE_FLAIR.replace('if sub in ("transcriptome", "collapse"):',
                                                   'if sub == "transcriptome": sys.exit(0)\nif sub == "collapse":', 1))
    p = assemble(run, tmp_path / "asm")
    assert p.returncode == 1
    assert "`flair transcriptome` finished without error but FLAIR wrote no isoform BED" in p.stderr
    assert "Log:" in p.stderr


@pytest.mark.parametrize("gen,method", [("v2", "predictProductivity"), ("v301", "predictProductivity"),
                                        ("ted", "inline"), ("dev", "inline")])
def test_productivity(fake, gen, method):
    run = fake(gen)
    asm = run.tmp / "asm"
    assemble(run, asm, check=0)
    out = run.tmp / "p"
    run("productivity", "--isoforms-bed", asm / "isoforms.bed", "--genome", run.inputs / "genome.fa",
        "--gtf", run.inputs / "annot.gtf", "--out", out, check=0)
    rows = read_tsv(out / "productivity.tsv")
    assert rows[0] == ["isoform_id", "cds_start", "cds_end", "class"]
    assert len(rows) == 4
    assert {r[3] for r in rows[1:]} <= {"PRO", "PTC", "NST", "NGO"}
    log = json.load(open(out / "adapter_log.json"))
    assert log["outputs"]["method"] == method
    if method == "predictProductivity":
        pp = [c for c in log["commands"] if c["label"] == "predictProductivity"][0]
        assert "--append_column" in pp["argv"] and "--longestORF" in pp["argv"]


def test_productivity_dev_flair_bed_flag(fake):
    run = fake("dev")
    asm = run.tmp / "asm"
    assemble(run, asm, check=0)
    other = run.tmp / "elsewhere"
    other.mkdir()
    (other / "isoforms.bed").write_text(open(asm / "isoforms.bed").read())
    out = run.tmp / "p"
    run("productivity", "--isoforms-bed", other / "isoforms.bed", "--genome", run.inputs / "genome.fa",
        "--gtf", run.inputs / "annot.gtf", "--out", out, check=0)
    assert (out / "skipped.json").exists()        # dev has no predictProductivity and no native BED given
    out2 = run.tmp / "p2"
    run("productivity", "--isoforms-bed", other / "isoforms.bed", "--genome", run.inputs / "genome.fa",
        "--gtf", run.inputs / "annot.gtf", "--flair-bed", asm / "isoforms.flair.bed", "--out", out2, check=0)
    assert json.load(open(out2 / "adapter_log.json"))["outputs"]["method"] == "inline"


def test_productivity_skipped_without_any_source(fake):
    run = fake("v301", with_pp=False)
    asm = run.tmp / "asm"
    assemble(run, asm, check=0)
    out = run.tmp / "p"
    run("productivity", "--isoforms-bed", asm / "isoforms.bed", "--genome", run.inputs / "genome.fa",
        "--gtf", run.inputs / "annot.gtf", "--out", out, check=0)
    assert "predictProductivity" in json.load(open(out / "skipped.json"))["reason"]


def test_thread_env_exported(fake, tmp_path):
    run = fake("v301")
    # flair is a fake; check the adapter sets the OMP family for FLAIR's children
    b = tmp_path / "env_v301" / "bin"
    write_exe(str(b / "flair"), FAKE_FLAIR.replace(
        'if argv == ["--version"]:',
        'open(os.environ["FAKE_FLAIR_CALLS"] + ".env", "w").write(os.environ.get("OMP_NUM_THREADS", ""))\n'
        'if argv == ["--version"]:', 1))
    out = tmp_path / "asm"
    assemble(run, out, check=0)
    assert open(str(tmp_path / "calls_v301.jsonl") + ".env").read() == "1"


def test_stdlib_only():
    src = open(ADAPTER).read()
    imports = set()
    for line in src.splitlines():
        line = line.strip()
        if line.startswith("import "):
            imports.update(x.strip().split(".")[0] for x in line[7:].split(","))
        elif line.startswith("from ") and " import " in line:
            imports.add(line.split()[1].split(".")[0])
    stdlib = getattr(sys, "stdlib_module_names", None)
    if stdlib is None:
        pytest.skip("needs Python 3.10+")
    assert imports <= set(stdlib) | {"__future__"}, imports - set(stdlib)


def test_dev_quantify_long_genome_flag(fake):
    """A20: dev 57040d2/425b260 quantify accept only --genome; the flag comes from --help."""
    run = fake("dev", extra_env={"FAKE_DEV_LONG_ONLY": "1"})
    asm = run.tmp / "asm"
    assemble(run, asm, check=0)
    out = run.tmp / "q"
    run("quantify", "--isoforms-fa", asm / "isoforms.fa", "--sample", "s", "--bam", run.inputs / "shared.bam",
        "--genome", run.inputs / "genome.fa", "--flair-bed", asm / "isoforms.flair.bed",
        "--threads", "2", "--out", out, check=0)
    q = [c for c in run.calls() if c[0] == "quantify"][-1]
    assert "-g" not in q and argval(q, "--genome").endswith("genome.fa")
    assert argval(q, "--threads") == "2" and "-t" not in q
    assert len(read_tsv(out / "counts.tsv")) == 4
    # the older help (short -g listed only) keeps the short form
    run2 = fake("dev")
    asm2 = run2.tmp / "asm2"
    assemble(run2, asm2, check=0)
    run2("quantify", "--isoforms-fa", asm2 / "isoforms.fa", "--sample", "s", "--bam", run2.inputs / "shared.bam",
         "--genome", run2.inputs / "genome.fa", "--flair-bed", asm2 / "isoforms.flair.bed",
         "--out", run2.tmp / "q2", check=0)
    assert "-g" in [c for c in run2.calls() if c[0] == "quantify"][-1]


@pytest.mark.parametrize("no_parallel", [False, True])
def test_dev_parallel_mode_byregion_when_supported(fake, no_parallel):
    """A20: --parallel_mode byregion whenever transcriptome --help lists it; nothing otherwise."""
    run = fake("dev", extra_env={"FAKE_DEV_NO_PARALLEL": "1"} if no_parallel else None)
    out = run.tmp / "asm"
    assemble(run, out, check=0)
    tx = [c for c in run.calls() if c[0] == "transcriptome"][0]
    log = json.load(open(out / "adapter_log.json"))
    if no_parallel:
        assert "--parallel_mode" not in tx and "bychrom" not in tx
        assert any("no --parallel_mode" in n for n in log["notes"])
    else:
        assert argval(tx, "--parallel_mode") == "byregion"


@pytest.mark.parametrize("gen,expected", [("v2", "reads"), ("v301", "reads"), ("ted", "reads"), ("dev", "bam")])
def test_quantify_input_recorded(fake, gen, expected):
    """A21: adapter_log.json says whether this version counts reads or the genome BAM."""
    run = fake(gen)
    asm = run.tmp / "asm"
    assemble(run, asm, check=0)
    out = run.tmp / "q"
    run("quantify", "--isoforms-fa", asm / "isoforms.fa", "--reads", run.inputs / "reads.fa", "--sample", "s",
        "--bam", run.inputs / "shared.bam", "--genome", run.inputs / "genome.fa",
        "--flair-bed", asm / "isoforms.flair.bed", "--out", out, check=0)
    assert json.load(open(out / "adapter_log.json"))["quantify_input"] == expected
    assert (out / "counts.tsv").exists()
    # recorded on a skip too, and null for other steps
    out2 = run.tmp / "q2"
    run("quantify", "--isoforms-fa", asm / "isoforms.fa", "--sample", "s", "--out", out2,
        *(["--bam", run.inputs / "shared.bam"]), check=0)
    assert json.load(open(out2 / "adapter_log.json"))["quantify_input"] == expected
    assert json.load(open(asm / "adapter_log.json"))["quantify_input"] is None


def commands(out):
    return json.load(open(out / "adapter_log.json"))["commands"]


@pytest.mark.parametrize("gen", ["v2", "v301", "dev"])
def test_assemble_builds_missing_fai_once(fake, gen):
    """A23: a genome without .fai is indexed once with samtools faidx before FLAIR runs."""
    run = fake(gen)
    genome = run.inputs / "genome.fa"
    os.remove(str(run.inputs / "genome.fa.fai"))
    out = run.tmp / "asm"
    assemble(run, out, check=0)
    assert (run.inputs / "genome.fa.fai").exists()
    labels = [c["label"] for c in commands(out)]
    assert labels.count("samtools_faidx") == 1
    first_flair = min(i for i, l in enumerate(labels) if l.startswith(("flair_", "bam2Bed12")))
    assert labels.index("samtools_faidx") < first_flair
    log = json.load(open(out / "adapter_log.json"))
    assert log["genome_index"] == {"genome": str(genome), "fai": str(genome) + ".fai", "action": "built"}
    assert any("A23" in n for n in log["notes"])
    # a second step finds the index and does not rebuild it
    out2 = run.tmp / "asm2"
    assemble(run, out2, check=0)
    assert "samtools_faidx" not in [c["label"] for c in commands(out2)]
    assert json.load(open(out2 / "adapter_log.json"))["genome_index"]["action"] == "existing"


def test_dev_quantify_indexes_genome_before_flair(fake):
    run = fake("dev")
    asm = run.tmp / "asm"
    assemble(run, asm, check=0)
    os.remove(str(run.inputs / "genome.fa.fai"))
    out = run.tmp / "q"
    run("quantify", "--isoforms-fa", asm / "isoforms.fa", "--sample", "s", "--bam", run.inputs / "shared.bam",
        "--genome", run.inputs / "genome.fa", "--flair-bed", asm / "isoforms.flair.bed", "--out", out, check=0)
    labels = [c["label"] for c in commands(out)]
    assert labels == ["samtools_faidx", "flair_quantify"]
    assert (run.inputs / "genome.fa.fai").exists()
    assert json.load(open(out / "adapter_log.json"))["genome_index"]["action"] == "built"


def test_reads_quantify_does_not_touch_genome(fake):
    run = fake("v301")
    asm = run.tmp / "asm"
    assemble(run, asm, check=0)
    os.remove(str(run.inputs / "genome.fa.fai"))
    out = run.tmp / "q"
    run("quantify", "--isoforms-fa", asm / "isoforms.fa", "--reads", run.inputs / "reads.fa", "--sample", "s",
        "--genome", run.inputs / "genome.fa", "--out", out, check=0)
    assert "samtools_faidx" not in [c["label"] for c in commands(out)]
    assert json.load(open(out / "adapter_log.json"))["genome_index"] is None


def test_fai_built_in_work_dir_when_genome_dir_read_only(fake, tmp_path):
    run = fake("v301")
    ro = tmp_path / "ro"
    ro.mkdir()
    (ro / "genome.fa").write_text(">chr1\nACGT\n")
    os.chmod(str(ro), 0o555)
    try:
        if os.access(str(ro), os.W_OK):
            pytest.skip("running as a user who can write anywhere")
        i = run.inputs
        out = tmp_path / "asm"
        run("assemble", "--bam", i / "shared.bam", "--genome", ro / "genome.fa", "--gtf", i / "annot.gtf",
            "--out", out, check=0)
        log = json.load(open(out / "adapter_log.json"))
        assert log["genome_index"]["action"] == "built_in_work_dir"
        tx = [c for c in run.calls() if c[0] == "transcriptome"][0]
        assert argval(tx, "-g", "--genome") == str(out / ".work" / "genome.fa")
    finally:
        os.chmod(str(ro), 0o755)


def test_no_samtools_warns_and_continues(fake, tmp_path):
    pybin = tmp_path / "pybin"
    pybin.mkdir()
    os.symlink(sys.executable, str(pybin / "python3"))
    run = fake("v301", extra_env={"PATH": "%s:%s" % (tmp_path / "env_v301" / "bin", pybin)})
    os.remove(str(run.tmp / "env_v301" / "bin" / "samtools"))
    os.remove(str(run.inputs / "genome.fa.fai"))
    out = run.tmp / "asm"
    p = assemble(run, out)
    if p.returncode != 0 or "samtools_faidx" in [c["label"] for c in commands(out)]:
        pytest.skip("a system samtools is on PATH")
    log = json.load(open(out / "adapter_log.json"))
    assert log["genome_index"]["action"] == "missing"
    assert any("no samtools" in w for w in log["warnings"])
