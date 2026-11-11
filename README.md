# flair-eval

Compare FLAIR runs (different versions, branches or settings) on the same long-read data.
Results come out as one self-contained HTML report: open any single run, or tick several
runs to compare them.

> **Status:** v2, the FLAIR-only rewrite, is on `main` (development continues on `v2`). See
> [docs/DESIGN.md](docs/DESIGN.md) and [docs/CONTRACTS.md](docs/CONTRACTS.md). The previous
> multi-assembler pipeline is in `legacy/` for reference.

**Supported FLAIR versions:** 2.2.0, 3.0.1 and `dev` (the tip of BrooksLabUCSC/flair `dev`,
resolved to a commit at launch). Other releases, branches, envs and containers can be named too,
but are not tested. See [docs/flair_versions.md](docs/flair_versions.md).

## Usage

flair-eval is set up for the lab's **phoenix** cluster (UCSC PRISM): with no `-profile` it runs
every step as a SLURM job, uses conda, and shares one FLAIR env cache with the whole lab
(`/private/groups/brookslab/flair-eval/env_cache`), so 2.2.0 and 3.0.1 are built once for everyone.
You need Nextflow ≥ 25.04 and Java 17+ on `PATH`. FLAIR itself is fetched per run.

`--outdir` is required. Put it on groups or warm-archive storage, never in home (~100 GB quota);
on phoenix Nextflow's work dir goes inside it (`<outdir>/work`).
Launch with `umask 002` so the shared caches stay writable for everyone (`flair-eval run` does this).

### One command, no files

```bash
nextflow run hdheath/flair-eval \
  --reads sample.fq.gz --genome hg38 \
  --flair 2.2.0,3.0.1,dev \
  --outdir /private/groups/brookslab/$USER/flair-eval/dd_mm_yy_name
# -> <outdir>/report/report.html
```

`--genome hg38` is a lab preset: GRCh38 primary assembly + GENCODE v49 basic from
`/private/groups/brookslab/reference_*`, and it fills `--gtf` unless you give one. It also brings
the end references and the TUSCO gene list used by the report's Ends and TUSCO panels (refTSS v4.1,
PolyA_DB v4.1 and SQANTI3's `tusco_human.tsv`, in `/private/groups/brookslab/flair-eval/reference/`). Any other genome:
`--genome genome.fa --gtf annotation.gtf`. The Nextflow driver itself runs for the whole
comparison, so start it inside a SLURM job or `tmux`, not as a bare login-node command; the
`flair-eval run` CLI ([docs/usage_cli.md](docs/usage_cli.md)) submits it for you.

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
  - id: flair-2.2.0
    flair: { release: "2.2.0" }
  - id: flair-3.0.1
    flair: { release: "3.0.1" }
  - id: dev-stringent
    flair: { git: "https://github.com/BrooksLabUCSC/flair", ref: "dev" }
    transcriptome_args: "--stringent"
    label: "dev, stringent"
    tags: { branch: dev, setting: stringent }
```

```bash
nextflow run hdheath/flair-eval --input samples.csv --runs runs.yaml \
  --outdir /private/groups/brookslab/$USER/flair-eval/dd_mm_yy_name
```

A `genome` cell may name a preset (`hg38`) and leave `gtf` empty.

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
   Quantify also runs for BAM-only samples: dev counts from the BAM, and 2.2.0 / 3.0.1 count
   the BAM's primary reads (extracted with `samtools fastq`).
4. One bundle per job, and one report over all bundles.

`--regions regions.bed` restricts assembly and evaluation to a set of regions. FLAIR then gets a
genome + GTF cut to the chromosomes in the BED (same coordinates), which keeps test regions fast.

`--ref_tss sites.bed`, `--ref_pas sites.bed` (BED6; poly(A) sites as [pos-1,pos)) and
`--tusco tusco_human.tsv` turn on the Ends and TUSCO evaluation for genomes without a preset, or
replace the preset's files. Each is optional; a missing one only hides its panel.

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
| none (= `phoenix`) | phoenix: `slurm` + `conda` + the shared lab env cache; work dir in `<outdir>/work` |
| `local` | run on this machine, requests capped at its CPUs and memory |
| `slurm` | SLURM on PRISM: partition short/medium/long picked from each task's time; sizes from measured use |
| `conda` | conda on: activates each run's FLAIR env and builds a small Python env for evaluation/report |
| `docker`, `singularity` | container engine on (needed for `container:` runs; evaluation/report use a pandas image) |
| `test` | the bundled test data (`tests/data/`) |

Off phoenix, combine an executor with an environment, e.g. `-profile local,conda` or
`-profile local,singularity`.
Without `conda`/`docker`/`singularity`, evaluation and report use the `python3` on `PATH`, which
needs pandas and numpy.

`--env_cache` and `-work-dir` override the phoenix defaults if you need a private cache or a
work dir outside the outdir.

### Help, checks and tests

```bash
nextflow run hdheath/flair-eval --help
```

```bash
nextflow run . -stub -profile local --input samples.csv --runs runs.yaml --outdir stub-out
```

```bash
nf-test test -c tests/nf-test/nf-test.config
```

`-stub` runs the whole graph with placeholder outputs, which is a quick way to check configs.
