# flair-eval

Compare FLAIR runs (different versions, branches or settings) on the same long-read data.
Results come out as one self-contained HTML report: open any single run, or tick several
runs to compare them.

> **Status:** v2 is under construction on the `v2` branch. See [docs/DESIGN.md](docs/DESIGN.md)
> and [docs/CONTRACTS.md](docs/CONTRACTS.md). The previous multi-assembler pipeline is in `legacy/`
> for reference, and on the `wip/2026-11-09-snapshot` branch.

## Usage

You need Nextflow ≥ 25.04 and Java 17+. FLAIR itself is fetched per run (see "FLAIR sources").

### One command, no files

```bash
nextflow run hdheath/flair-eval -profile local,conda \
  --reads sample.fq.gz --genome genome.fa --gtf annotation.gtf \
  --flair 3.0.1,dev --outdir results
# -> results/report/report.html
```

`--flair` takes comma-separated tokens; each one becomes a run:

| token | FLAIR used | run id |
|---|---|---|
| `3.0.1`, `release:3.0.1` | that release | `flair-3.0.1` |
| `dev`, `main`, `git:<ref>` | BrooksLabUCSC/flair at that branch, tag or commit | `flair-<ref>` |
| `git:<url>@<ref>` | another repo at that ref | `flair-<ref>` |
| `env:/abs/path` | an existing conda env | `env-<basename>` |
| `container:<uri>` | a Docker/Singularity image | `ctr-<name>` |

Add `--flair_args='--stringent'` (with `=`; the space form is mis-parsed by Nextflow) to pass extra arguments to every run's transcriptome step.
Use `--bam sorted.bam` (with its `.bai`) instead of `--reads` to skip alignment, and
`--sample`, `--platform pacbio|ont`, `--library cDNA|dRNA` to describe the sample.

### Config files

For several samples, or runs that differ in settings, write `samples.csv` and `runs.yaml`:

```csv
sample,reads,bam,genome,gtf,platform,library
A549_pb,A549.fq.gz,,hg38.fa,gencode.gtf,pacbio,cDNA
WTC11_ont,,WTC11.bam,hg38.fa,gencode.gtf,ont,cDNA
```

```yaml
runs:
  - id: flair-3.0.1
    flair: { release: "3.0.1" }
  - id: dev-stringent
    flair: { git: "https://github.com/BrooksLabUCSC/flair", ref: "dev" }
    transcriptome_args: "--stringent"
    label: "dev, stringent"
    tags: { branch: dev, setting: stringent }
```

```bash
nextflow run hdheath/flair-eval -profile local,conda \
  --input samples.csv --runs runs.yaml --outdir results
```

Every run is applied to every sample. Relative paths in `samples.csv` are read next to the CSV
first, then from the launch directory. Columns and keys are listed in
[docs/CONTRACTS.md](docs/CONTRACTS.md) section 1 and checked against
`assets/schema/*.json`; all problems are reported at once, in plain English.
`--input` cannot be mixed with the one-command sample flags, nor `--runs` with `--flair`.

### What runs

1. One FLAIR environment per distinct FLAIR source (cached under `--env_cache`).
2. A shared alignment per sample with a pinned FLAIR (`--align_flair`, default 3.0.1), skipped
   when a BAM is given. A run with `own_align: true` aligns with its own FLAIR instead.
3. Per sample × run: assemble, then quantify and predictProductivity (switch off per run with
   `quantify: false` / `productivity: false`), then structure evaluation against the GTF.
   Quantify also runs for BAM-only samples; a FLAIR version that needs reads records the step
   as skipped instead of failing.
4. One bundle per job, and one report over all bundles.

`--regions regions.bed` restricts assembly and evaluation to a set of regions. FLAIR then gets a
genome + GTF cut to the chromosomes in the BED (same coordinates), which keeps test regions fast.

### Reports over many launches

Bundles accumulate. Add earlier ones to a new report with `--bundles`, or skip compute entirely:

```bash
nextflow run hdheath/flair-eval --report_only \
  --bundles 'old_results/bundles,other/bundles/*/flair-dev' --outdir combined
```

`--baseline <run_id>` picks the run that differences are shown against. By default it is the
first run in runs.yaml or `--flair` order; with `--bundles` it is the first run of the earliest
launch, so adding a run later keeps the old baseline. The report lets you switch it.
`--report_title` sets the heading. Output layout: [docs/output.md](docs/output.md).

### Profiles

| profile | what it does |
|---|---|
| `local` | run on this machine |
| `slurm` | SLURM on PRISM: partition short/medium/long picked from each task's time; sizes from measured use |
| `conda` | conda on: activates each run's FLAIR env and builds a small Python env for evaluation/report |
| `docker`, `singularity` | container engine on (needed for `container:` runs; evaluation/report use a pandas image) |
| `test` | the bundled test data (`tests/data/`) |

Combine an executor with an environment, e.g. `-profile slurm,conda` or `-profile local,singularity`.
Without `conda`/`docker`/`singularity`, evaluation and report use the `python3` on `PATH`, which
needs pandas and numpy.

On PRISM, keep results and caches off home (quota ~100 GB):

```bash
nextflow run hdheath/flair-eval -profile slurm,conda \
  --input samples.csv --runs runs.yaml \
  --outdir /private/groups/brookslab/$USER/results/flair-eval/dd_mm_yy_name \
  --env_cache /private/groups/brookslab/$USER/flair-env-cache \
  -work-dir /private/groups/brookslab/$USER/work/flair-eval
```

### Help, checks and tests

```bash
nextflow run hdheath/flair-eval --help
```

```bash
nextflow run . -stub -profile local --input samples.csv --runs runs.yaml --outdir /tmp/stub
```

```bash
nf-test test -c tests/nf-test/nf-test.config
```

`-stub` runs the whole graph with placeholder outputs, which is a quick way to check configs.
