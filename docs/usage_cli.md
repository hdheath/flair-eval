# Setting up and launching flair-eval

flair-eval compares FLAIR runs (different versions, branches or settings) on the same long-read
data and writes one HTML report. There are three ways to set up a launch. All of them produce the
same `nextflow run` command and the same config files, and all follow the rules in
[CONTRACTS.md](CONTRACTS.md) §1.

| Way in | Best for |
|---|---|
| [One command](#one-command-no-files) | One sample, a few FLAIR versions with the same settings |
| [Config-builder page](#config-builder-page) | Filling in a form, then copying the command or downloading the files |
| [`flair-eval` CLI](#the-flair-eval-cli) | A guided setup in the terminal, checking configs, launching |

## First time? Start here
```bash
conda create -y -n nextflow -c conda-forge -c bioconda nextflow openjdk=17
conda activate nextflow
pip install "flair-eval @ git+https://github.com/hdheath/flair-eval@v2#subdirectory=cli"
flair-eval doctor          # checks Java, Nextflow, conda, containers, disk; suggests a profile
flair-eval run --test      # a small run on the bundled test data
flair-eval init            # then set up your own data
```
Nextflow needs Java 17+; flair-eval needs Nextflow 25.04+. `flair-eval doctor` says what is
missing and how to fix it. To install the CLI from a clone instead: `pip install -e cli/`.

## One command, no files
```bash
nextflow run hdheath/flair-eval -r v2 \
  --reads /data/HG002.fq.gz --genome /ref/GRCh38.fa --gtf /ref/gencode.gtf \
  --flair 3.0.1,dev \
  --outdir /big/storage/flair-eval-results
```
This makes one sample (named after the reads file, or `--sample NAME`) and one run per
`--flair` token. Add `--platform pacbio|ont` and `--library cDNA|dRNA` if you know them.

### `--flair` tokens
| token | FLAIR used | run id |
|---|---|---|
| `3.0.1` or `release:3.0.1` | release 3.0.1 (bioconda) | `flair-3.0.1` |
| `dev`, `main`, `git:<ref>` | that branch, tag or commit of `https://github.com/BrooksLabUCSC/flair` | `flair-<ref>` |
| `git:<url>@<ref>` | `<ref>` of another repository | `flair-<ref>` |
| `env:/path/to/env` | an existing conda env, used as is | `env-<env folder name>` |
| `container:<image>` | a Docker/Singularity image | `ctr-<image name>` |

Details:
- A bare token that starts with a digit is a release. Any other bare word except `dev` and
  `main` is an error; write `git:<branch>` for other branches.
- In `git:<url>@<ref>` the URL is everything before the last `@`. It can be `https://…`,
  `git@host:org/repo.git` or a local path such as `/home/me/flair`.
- Characters other than letters, digits, `.`, `_` and `-` become `-` in the run id, so
  `git:feature/ends` gets the id `flair-feature-ends`.
- The container name drops the registry path, the `:tag` and any `@sha256:` digest:
  `container:docker://quay.io/biocontainers/flair:3.0.1--pyhdfd78af_0` gets the id `ctr-flair`.
  A local image keeps its file name: `container:/images/flair.sif` gets `ctr-flair.sif`.
- Two tokens that give the same run id are an error, and so are two runs in `runs.yaml` with the
  same id. The message names both; write those runs in `runs.yaml` with ids you choose.

### Passing FLAIR arguments
`--flair_args` goes to `flair transcriptome` for every `--flair` run. **Write it with `=`**:
```bash
--flair_args='--stringent --check_splice'
```
Nextflow reads a value that starts with `-` as a new option. `--flair_args --stringent` sets
`flair_args` to `true` and makes a stray `stringent` param. The builder page and the CLI always
write the `=` form for you. (`flair-eval run --flair_args "--stringent"` is fine: the CLI fixes it.)

Use `runs.yaml` when runs need different arguments, labels, tags, or other step settings, and
`samples.csv` when you have more than one sample or extra tracks.

## Config files
`samples.csv` has one row per sample:
```csv
sample,reads,bam,genome,gtf,platform,library
HG002,/data/HG002.fq.gz,,/ref/GRCh38.fa,/ref/gencode.gtf,pacbio,cDNA
WTC11,,/data/WTC11.bam,/ref/GRCh38.fa,/ref/gencode.gtf,ont,dRNA
```
Each row needs `sample`, `genome`, `gtf`, and `reads` or `bam` (a sorted BAM with its `.bai`
next to it). Optional columns `cage`, `pas`, `junctions` and `spikein_gtf` are read now and used
by later report sections. A relative path is looked up next to `samples.csv` first, then in the
folder you launch from. Full paths avoid any doubt.

`runs.yaml` has one entry per FLAIR configuration:
```yaml
runs:
  - id: "flair-3.0.1"
    flair:
      release: "3.0.1"
    transcriptome_args: ""
  - id: "dev-stringent"
    label: "dev, stringent"
    flair:
      git: "https://github.com/BrooksLabUCSC/flair"
      ref: "dev"
    transcriptome_args: "--stringent"
    tags:
      setting: "stringent"
```
Each run gives exactly one of `release`, `git` (with `ref`), `env` or `container`. Optional keys:
`transcriptome_args` (default empty), `quantify` and `productivity` (default `true`),
`own_align` (default `false`: all runs share one alignment per sample), `label` (name in the
report) and `tags` (text key/value pairs used as report filters).
Put version numbers in quotes: YAML reads `release: 3.0` as the number 3, not the text "3.0".

Launch with:
```bash
nextflow run hdheath/flair-eval -r v2 --input samples.csv --runs runs.yaml --outdir /big/storage/out
```

## Config-builder page
`frontend/builder/index.html` is a single static page (no build step, nothing sent anywhere).
Open it from a clone, or from GitHub Pages once that is enabled
(repository Settings → Pages → Source: GitHub Actions). The `Builder page` workflow then
publishes `https://hdheath.github.io/flair-eval/`: a start page with a three-step quick start,
the builder at `/builder/`, and a demo report built from `tests/fixtures/bundles` at
`/demo/report.html`.

- Add samples, then one card per run. Each card has a FLAIR source picker (release, git branch,
  conda env, container), common setting presets, a label and tags.
- The panel on the right checks everything against the pipeline's schemas as you type and shows
  the exact command. It uses the one-command form whenever that says exactly the same thing,
  and otherwise tells you why it needs `samples.csv` and `runs.yaml`.
- Pick where jobs run (this computer or a SLURM cluster) and how FLAIR gets installed (conda,
  Singularity/Apptainer, Docker); the page turns that into `-profile`.
- Copy the command, and download or copy the two files.

## The `flair-eval` CLI
Python 3.10+. Install it straight from GitHub, or from a clone:
```bash
pip install "flair-eval @ git+https://github.com/hdheath/flair-eval@v2#subdirectory=cli"
pip install -e cli/
```

### `flair-eval doctor`
Checks this machine and prints a fix for each problem:
```text
$ flair-eval doctor --outdir /big/storage/out
  ok    Java: version 17
  FAIL  Nextflow: not found
        Nextflow is not on PATH. Install one of:
            conda create -y -n nextflow -c conda-forge -c bioconda nextflow openjdk=17
  ok    conda: found: FLAIR releases and git branches can be built
  WARN  Containers: no singularity, apptainer or docker: container: runs will not work
  ok    Cluster: SLURM found (sbatch): jobs can run on the cluster
  ok    Disk: 2,310 GB free at /big/storage

Suggested profile for this machine: -profile slurm,conda
```
It also checks the Nextflow version (25.04+), git (for git runs), and warns when `--outdir` is in
your home folder.

### `flair-eval init`
Asks a few questions and writes `samples.csv` and `runs.yaml`:
```text
$ flair-eval init
Sample 1
  Reads (FASTQ/FASTA, .gz is fine) or a sorted, indexed BAM: /data/HG002.fq.gz
  Sample name [HG002]:
  Genome FASTA: /ref/GRCh38.fa
  ...
Runs: which FLAIR builds do you want to compare?
  FLAIR [3.0.1]: 3.0.1,dev
  Settings for flair-3.0.1:
    1) FLAIR defaults
    2) Stringent  (--stringent)
    ...
```
The first question offers the bundled test data instead (answer `y` and nothing is written).
Platform and library default to a guess from the file name (`…_ONT_cDNA.fq.gz` gives ont and
cDNA); type `-` to leave one unset.
Without questions: `flair-eval init --reads r.fq.gz --genome g.fa --gtf a.gtf --flair 3.0.1,dev`.
`--dir DIR` writes elsewhere; `--force` overwrites existing files.

### `flair-eval check`
Validates `samples.csv` and `runs.yaml` (in the current folder, or `--samples`/`--runs`) and
explains each problem:
```text
$ flair-eval check
samples.csv: 1 problem(s)
runs.yaml: 1 problem(s)
  error: samples.csv row 2: platform: 'nanopore' is not allowed; use one of pacbio, ont
  error: runs.yaml run 1: flair: git needs a ref: a branch, tag or commit SHA (e.g. ref: dev)
  warning: samples.csv row 1: reads: file not found: /data/HG002.fq.gz
```
Input files that do not exist here are warnings, since configs are often written on one machine
and run on another. `--strict` turns warnings into errors; `--no-files` skips the file checks.
`--bundles DIR …` validates every `bundle.json` under those folders.

### `flair-eval run`
Checks the configs and Java/Nextflow, prints the command, and runs it:
```bash
flair-eval run --outdir /big/storage/out --profile slurm             # uses ./samples.csv + ./runs.yaml
flair-eval run --samples cfg/samples.csv --runs cfg/runs.yaml --outdir out
flair-eval run --reads r.fq.gz --genome g.fa --gtf a.gtf --flair 3.0.1,dev --outdir out
flair-eval run --test                                                 # bundled test data, -profile test
flair-eval run ... --dry-run                                          # print the command only
flair-eval run ... -- -with-tower                                     # after --: passed to nextflow
```
Without `--profile`, `run` and `report` use the profile `flair-eval doctor` suggests for this
machine (SLURM if `sbatch` exists, else local; plus conda, else Singularity/Apptainer, else
Docker) and say so. `--test` adds `test` in front and writes to `./flair-eval-test`.
Options: `-r/--ref` (pipeline revision, default `v2`), `--profile`, `--outdir`, `--work-dir`,
`--resume`, and the pipeline params `--baseline`, `--report_title`, `--align_flair`,
`--regions`, `--bundles`, `--env_cache`. Pipeline params equal to their default are left off.

### `flair-eval report`
Builds a report from bundles you already have, without running FLAIR:
```bash
flair-eval report out_nov/bundles out_dec/bundles --title "Nov vs Dec" --outdir report_out
```
This runs `nextflow run hdheath/flair-eval --report_only true --bundles out_nov/bundles,out_dec/bundles …`.
If two bundles share a sample and run id, the newest one is used.

## On a shared cluster
- Point `--outdir` (and `--work-dir`) at project or scratch storage, not your home folder.
- Use the cluster profile (`--profile slurm`) so each step is its own job.
- Set `--env_cache` to shared storage so each FLAIR version is built once and reused.

## Keeping the front ends in step
The page (`frontend/builder/logic.js`) and the CLI (`cli/src/flair_eval/`) implement the same
rules. `cli/tests/test_builder_parity.py` runs both on the same cases (needs Node) and fails if
they disagree. Both read copies of `assets/schema/*.json`; after a schema change run
`python cli/tools/sync_schemas.py` (CI runs it with `--check`).
