# flair-eval

Compare FLAIR runs (different versions, branches or settings) on the same long-read data.
The result is one self-contained HTML report. You can open any single run in it, or tick
several runs to compare them against a baseline.

- **Supported FLAIR versions:** 2.2.0, 3.0.1 and `dev` (the tip of BrooksLabUCSC/flair `dev`,
  resolved to a commit at launch). You can also name other releases, branches, conda envs or
  containers, but they are not tested ([docs/flair_versions.md](docs/flair_versions.md)).
- **Where it runs:** the pipeline is set up for the lab's **phoenix** cluster (UCSC PRISM). There
  you need no `-profile`, `--genome hg38` picks the lab reference, and FLAIR envs come from a
  shared lab cache. Elsewhere, see [Other machines](docs/usage.md#other-machines-and-clusters).
- **Demo:** [a report built from real FLAIR 2.2.0, 3.0.1 and dev runs](https://hdheath.github.io/flair-eval/demo/report.html),
  and a [web page that writes the launch command for you](https://hdheath.github.io/flair-eval/builder/).

## Quick start

You need Nextflow 25.04+ with Java 17+, and Python 3.10+ for the `flair-eval` helper CLI.

```bash
conda create -y -n nextflow -c conda-forge -c bioconda nextflow openjdk=17
conda activate nextflow
pip install "flair-eval @ git+https://github.com/hdheath/flair-eval@main#subdirectory=cli"
flair-eval doctor
```

`flair-eval doctor` checks Java, Nextflow, conda and disk space, and says whether it found phoenix.

Try the bundled test data, a 2 Mb human region with 20,000 PacBio and 20,000 ONT reads, which
FLAIR 3.0.1, 2.2.0 and dev each finish in minutes:

```bash
flair-eval run --test --outdir /private/groups/brookslab/$USER/flair-eval/dd_mm_yy_test
```

Then compare FLAIR versions on your own reads:

```bash
flair-eval run --reads /data/sample.fq.gz --genome hg38 --flair 2.2.0,3.0.1,dev \
  --outdir /private/groups/brookslab/$USER/flair-eval/dd_mm_yy_name
```

On phoenix, `flair-eval run` submits the Nextflow driver as a SLURM job and prints the job id
and where the report will be: `<outdir>/report/report.html`. `--outdir` is always required. Put
it on groups or warm-archive storage, never in your home folder (about 100 GB quota).

## Documentation

| Read | For |
|---|---|
| [docs/usage.md](docs/usage.md) | Setting up and launching: `--flair` tokens, reads or BAM, `samples.csv` and `runs.yaml`, the CLI, the web builder, plain `nextflow run`, profiles, report-only mode |
| [docs/output.md](docs/output.md) | What lands in `--outdir`, what each bundle holds, and how to read the report |
| [docs/flair_versions.md](docs/flair_versions.md) | How each FLAIR version is run and what differs between them |
| [docs/test_data.md](docs/test_data.md) | The bundled test data: sources, licences, how it was built |

**For developers:** [docs/CONTRACTS.md](docs/CONTRACTS.md) defines the interfaces between the
pipeline, the FLAIR adapter, the report and the front ends. [docs/DESIGN.md](docs/DESIGN.md)
describes the architecture. Run the tests with:

```bash
pytest tests/python cli/tests
```

```bash
nf-test test -c tests/nf-test/nf-test.config
```

The previous multi-assembler pipeline is kept in `legacy/` for reference only.
