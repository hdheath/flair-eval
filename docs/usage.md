# Setting up and launching flair-eval

flair-eval runs every FLAIR configuration you list (a *run*) on every dataset you list (a
*sample*), then builds one report over all of them. You can set up a launch in three ways. All
three produce the same `nextflow run` command and the same config files.

| Way in | Best for |
|---|---|
| [`flair-eval` CLI](#the-flair-eval-cli) | Guided setup in the terminal, config checks, launching (submits the job for you on phoenix) |
| [Config-builder page](#config-builder-page) | Filling in a web form, then copying the command or downloading the files |
| [Plain `nextflow run`](#running-nextflow-directly) | Scripts, or if you already know Nextflow |

Install the CLI and run the test data first: see the [README quick start](../README.md#quick-start).

## Before you start

flair-eval assumes the lab's **phoenix** cluster (UCSC PRISM; login node `mustard`). There:
- you need no `-profile`. The default profile runs every step as a SLURM job, installs FLAIR with
  conda, and takes FLAIR envs from the shared lab cache
  (`/private/groups/brookslab/flair-eval/env_cache`), so each version is built once for everyone;
- `--genome hg38` picks the lab reference: GRCh38 primary assembly + GENCODE v49 basic. It also
  fills `--gtf` unless you give one, and brings the end references and gene list that the
  report's Ends and TUSCO panels use (refTSS v4.1, PolyA_DB v4.1, SQANTI3 `tusco_human.tsv`, in
  `/private/groups/brookslab/flair-eval/reference/`);
- `flair-eval run` submits the Nextflow driver as a SLURM job, so nothing long runs on the login node.

**`--outdir` is always required.** On phoenix put it on groups or warm-archive storage, never in
your home folder (about 100 GB quota), e.g.
`--outdir /private/groups/brookslab/$USER/flair-eval/dd_mm_yy_name`. Nextflow's work directory
goes inside it (`<outdir>/work`).

Not on phoenix? See [Other machines and clusters](#other-machines-and-clusters).

## One sample, several FLAIR versions

```bash
flair-eval run --reads /data/HG002.fq.gz --genome hg38 \
  --flair 2.2.0,3.0.1,dev \
  --outdir /private/groups/brookslab/$USER/flair-eval/dd_mm_yy_hg002
```

This makes one sample and one run per `--flair` token. The sample is named after the reads or
BAM file unless you give `--sample NAME`. Add `--platform pacbio|ont` and `--library cDNA|dRNA`
if you know them. To start from an aligned BAM, give `--bam /data/HG002.bam` instead of
`--reads` (see [Reads or a BAM](#reads-or-a-bam)). With your own reference, give
`--genome /ref/genome.fa --gtf /ref/annotation.gtf` instead of `hg38`.

### `--flair` tokens
| token | FLAIR used | run id |
|---|---|---|
| `3.0.1` or `release:3.0.1` | that release | `flair-3.0.1` |
| `dev`, `main`, `git:<ref>` | that branch, tag or commit of `https://github.com/BrooksLabUCSC/flair` | `flair-<ref>` |
| `git:<url>@<ref>` | `<ref>` of another repository | `flair-<ref>` |
| `env:/path/to/env` | an existing conda env, used as is | `env-<env folder name>` |
| `container:<image>` | a Docker/Singularity image | `ctr-<image name>` |

- A bare token that starts with a digit is a release. Any other bare word except `dev` and
  `main` is an error: write `git:<branch>` for other branches.
- In `git:<url>@<ref>` the URL is everything before the last `@`. It can be `https://…`,
  `git@host:org/repo.git` or a local path such as `/home/me/flair`.
- Characters other than letters, digits, `.`, `_` and `-` become `-` in the run id, so
  `git:feature/ends` gets the id `flair-feature-ends`.
- The container name drops the registry path, the `:tag` and any `@sha256:` digest:
  `container:docker://quay.io/biocontainers/flair:3.0.1--pyhdfd78af_0` gets the id `ctr-flair`.
  A local image keeps its file name: `container:/images/flair.sif` gets `ctr-flair.sif`.
- Two runs with the same id are an error, whether they come from tokens or from `runs.yaml`. The
  message names both. Write such runs in `runs.yaml` with ids you choose.

### Passing FLAIR arguments
`--flair_args` goes to `flair transcriptome` (2.2.0: `collapse`) for every `--flair` run. With
plain `nextflow run`, **write it with `=`**:
```bash
--flair_args='--stringent --check_splice'
```
Nextflow reads a value that starts with `-` as a new option, so `--flair_args --stringent` would
set `flair_args` to `true` (the pipeline stops and says how to quote it). The CLI accepts either
form, and the CLI and the builder page always write the `=` form for you.

Each version runs at its own FLAIR defaults unless you pass arguments, and those defaults differ a
lot between versions ([docs/flair_versions.md](flair_versions.md#2-assemble-transcriptome--collapse)).
When runs need different arguments, labels or tags, use `runs.yaml`.

### Reads or a BAM
| sample gives | what happens |
|---|---|
| reads only | aligned **once per sample** with `flair align` from FLAIR `--align_flair` (default 3.0.1, FLAIR's own minimap2 settings; `--nvrna` when `library` is dRNA). Every run assembles from that one BAM, so differences between runs come from FLAIR, not from the alignment. |
| BAM only | no alignment; every run assembles from your BAM. Quantify still runs: dev counts from the BAM, and 2.2.0 / 3.0.1 count the BAM's primary alignments, extracted with `samtools fastq -F 0x904` (so unmapped reads are not counted). |
| both | your BAM is used for assembly; the reads file is used where a version wants reads (2.2.0 `collapse`, quantify). |

The BAM must be coordinate-sorted, have its `.bai` (or `.csi`) next to it, and be aligned to the
same genome FASTA you give (`hg38` = GRCh38 primary assembly). Any spliced aligner is fine
(minimap2 `-ax splice`, pbmm2, an earlier `flair align`).

A run with `own_align: true` in `runs.yaml` ("Align its own BAM" on the builder page) aligns the
reads with its own FLAIR instead. It needs a reads file.

## Several samples or different settings: config files

Write `samples.csv` and `runs.yaml` (or let `flair-eval init` or the builder page write them).
Every run is applied to every sample.

`samples.csv` has one row per sample:
```csv
sample,reads,bam,genome,gtf,platform,library
HG002,/data/HG002.fq.gz,,/ref/GRCh38.fa,/ref/gencode.gtf,pacbio,cDNA
WTC11,,/data/WTC11.bam,hg38,,ont,dRNA
```
- Required: `sample`, `genome`, `gtf`, and `reads` or `bam` (a sorted BAM with its index next to it).
- `genome` may be a preset (`hg38`); the `gtf` cell may then be left empty.
- Optional: `platform` (`pacbio`, `ont`) and `library` (`cDNA`, `dRNA`). The columns `cage`, `pas`,
  `junctions` and `spikein_gtf` are accepted and checked, but no report section uses them yet.
- A relative path is looked up next to `samples.csv` first, then in the folder you launch from.
  Full paths avoid any doubt. `flair-eval run` on phoenix starts the job in `--outdir`, so it
  refuses a path that only works from your current folder.

`runs.yaml` has one entry per FLAIR configuration:
```yaml
runs:
  - id: "flair-3.0.1"
    flair:
      release: "3.0.1"
  - id: "dev-stringent"
    label: "dev, stringent"
    flair:
      git: "https://github.com/BrooksLabUCSC/flair"
      ref: "dev"
    transcriptome_args: "--stringent"
    tags:
      setting: "stringent"
```
- Each run gives exactly one of `release`, `git` (with `ref`), `env` or `container`.
- Optional keys: `transcriptome_args` (default empty), `quantify` and `productivity` (default
  `true`), `own_align` (default `false`: all runs share one alignment per sample), `label` (the
  name shown in the report) and `tags` (text key/value pairs, used as report filters).
- Put version numbers in quotes: YAML reads `release: 3.0` as the number 3, not the text "3.0".
- The first run is the report's default baseline.

Launch with:
```bash
flair-eval run --samples samples.csv --runs runs.yaml --outdir /private/groups/brookslab/$USER/flair-eval/dd_mm_yy_name
```
`--input` (samples.csv) cannot be mixed with the one-sample flags, nor `--runs` with `--flair`.
Every problem in the files is reported at once, in plain English.

## What runs

1. One FLAIR environment per distinct FLAIR source, built once and cached (`--env_cache`).
2. A shared alignment per sample, skipped when the sample has a BAM.
3. Per sample × run: assemble, then quantify and predictProductivity (turn them off per run with
   `quantify: false` / `productivity: false`), then structure evaluation against the GTF, end
   precision against refTSS / PolyA_DB, TUSCO, internal priming and isoform variation.
4. One results *bundle* per sample × run, and one report over all bundles.

What each step writes: [docs/output.md](output.md).

## Optional settings

| setting | CLI / builder | what it does |
|---|---|---|
| `--regions regions.bed` | yes | restrict assembly and evaluation to these regions. FLAIR gets the genome and GTF cut to the BED's chromosomes (same coordinates), which keeps test regions fast. |
| `--align_flair 3.0.1` | yes | FLAIR release used for the shared alignment |
| `--baseline <run_id>` | yes | the run that differences are shown against (default: the first run; the report can switch it) |
| `--report_title "…"` | yes | the report heading |
| `--bundles <dirs>` | yes | add bundles from earlier launches to this report ([below](#reports-over-many-launches)) |
| `--env_cache <dir>` | yes | where FLAIR envs are built (phoenix default: the shared lab cache) |
| `--ref_tss`, `--ref_pas`, `--tusco` | **no**: plain `nextflow run` only, hidden from `--help` | end references (BED6; poly(A) sites as `[pos-1,pos)`) and the SQANTI3 TUSCO list, for genomes without a preset or to replace the preset's files. Each is optional; a missing one only hides its report panels. |

## The `flair-eval` CLI
Python 3.10+. Install from GitHub (or `pip install -e cli/` from a clone):
```bash
pip install "flair-eval @ git+https://github.com/hdheath/flair-eval@main#subdirectory=cli"
```

### `flair-eval doctor`
Checks this machine and prints a fix for each problem:
```text
$ flair-eval doctor --outdir /private/groups/brookslab/$USER/flair-eval
  ok    Java: version 17
  FAIL  Nextflow: not found
        Nextflow is not on PATH. Install one of:
            conda create -y -n nextflow -c conda-forge -c bioconda nextflow openjdk=17
  ok    conda: found: FLAIR releases and git branches can be built
  WARN  Containers: no singularity, apptainer or docker: container: runs will not work
  ok    Cluster: phoenix (UCSC PRISM) detected: SLURM jobs, conda, no -profile needed
  ok    FLAIR envs: shared lab cache /private/groups/brookslab/flair-eval/env_cache
  ok    Disk: 61,200 GB free at /private/groups/brookslab
```
Phoenix is recognised by `sbatch` on `PATH` plus `/private/groups/brookslab`. Elsewhere doctor
suggests a profile instead (`slurm` if `sbatch` exists, else `local`; plus conda, else
Singularity/Apptainer, else Docker). It also checks the Nextflow version (25.04+), git (for git
runs), and warns when `--outdir` is in your home folder.

### `flair-eval init`
Asks a few questions and writes `samples.csv` and `runs.yaml`:
```text
$ flair-eval init
Sample 1
  Reads (FASTQ/FASTA, .gz is fine) or a sorted, indexed BAM: /data/HG002.fq.gz
  Sample name [HG002]:
  Genome FASTA (or a lab preset: hg38) [hg38]:
  ...
Runs: which FLAIR builds do you want to compare?
  FLAIR [3.0.1]: 2.2.0,3.0.1,dev
  Settings for flair-3.0.1:
    1) FLAIR defaults
    2) Stringent  (--stringent)
    ...
```
- The first question offers the bundled test data instead (answer `y` and nothing is written).
- Platform and library default to a guess from the file name (`…_ONT_cDNA.fq.gz` gives ont and
  cDNA); type `-` to leave one unset.
- On phoenix the genome defaults to `hg38` and the output folder to a dated folder on groups
  storage; elsewhere you type both. The output folder is not written to a file: init prints the
  `flair-eval run … --outdir …` command to use next.
- Without questions: `flair-eval init --reads r.fq.gz --genome hg38 --flair 2.2.0,3.0.1,dev`
  (`--outdir` is optional and only used in the printed command). `--dir DIR` writes elsewhere;
  `--force` overwrites existing files.

### `flair-eval check`
Validates `samples.csv` and `runs.yaml` (in the current folder, or `--samples`/`--runs`) against
the pipeline's schemas and explains each problem:
```text
$ flair-eval check
samples.csv: 1 problem(s)
runs.yaml: 1 problem(s)
  error: samples.csv row 2: platform: 'nanopore' is not allowed; use one of pacbio, ont
  error: runs.yaml run 1: flair: git needs a ref: a branch, tag or commit SHA (e.g. ref: dev)
  warning: samples.csv row 1: reads: file not found: /data/HG002.fq.gz
```
Input files that do not exist here are warnings, since configs are often written on one machine
and run on another. A genome preset is checked as its lab files (a warning off phoenix).
`--strict` turns warnings into errors, and `--no-files` skips the file checks.
`--bundles DIR …` validates every `bundle.json` under those folders.

### `flair-eval run`
Checks the configs and Java/Nextflow, then launches. `--outdir` is required.
```bash
OUT=/private/groups/brookslab/$USER/flair-eval/dd_mm_yy_hg002
flair-eval run --outdir $OUT                                      # uses ./samples.csv + ./runs.yaml
flair-eval run --samples cfg/samples.csv --runs cfg/runs.yaml --outdir $OUT
flair-eval run --reads r.fq.gz --genome hg38 --flair 2.2.0,3.0.1,dev --outdir $OUT
flair-eval run --test --outdir $OUT                               # bundled test data
flair-eval run ... --dry-run                                      # print only (on phoenix: the job script)
flair-eval run ... --foreground                                   # on phoenix: run in this shell
flair-eval run ... --resume                                       # continue an earlier launch
flair-eval run ... -- -stub                                       # after --: passed to nextflow
```
**On phoenix** `run` writes `<outdir>/pipeline_info/launch.sbatch` and submits it:
```text
Submitted the Nextflow driver as SLURM job 812345 (partition long, 2 CPUs, 8G, up to 7-00:00:00).
  script:  /private/groups/brookslab/<you>/flair-eval/dd_mm_yy_hg002/pipeline_info/launch.sbatch
  log:     /private/groups/brookslab/<you>/flair-eval/dd_mm_yy_hg002/pipeline_info/nextflow-812345.log
  status:  squeue -j 812345
  cancel:  scancel 812345
The report will be at /private/groups/brookslab/<you>/flair-eval/dd_mm_yy_hg002/report/report.html when the job finishes.
```
- The driver job only schedules the pipeline; every step is its own SLURM job.
- The job sets `umask 002`, so other lab members can write lock files in the shared caches.
- It starts in `--outdir`, so `.nextflow/` (the `-resume` cache) and `.nextflow.log` live with
  the results. Every path in the script is made absolute.
- It inherits your shell's environment: activate the Nextflow conda env before launching.
- To resume, run the same `flair-eval run … --resume` again.

Off phoenix, `run` runs Nextflow in the current shell, so start it from a batch job or `tmux`.

Without `--profile`, `run` and `report` pick the profile: none on phoenix, otherwise the one
`flair-eval doctor` suggests, and they say which. `--test` adds `test` in front (`test,phoenix`
on phoenix, since naming any profile replaces the default one).

Other options:
- `-r/--ref`: pipeline revision (default `main`; `-r v2` runs the development branch).
- `--pipeline`: another copy of the pipeline, a GitHub fork (`--pipeline me/flair-eval -r my-branch`)
  or a local checkout (`--pipeline ~/flair-eval`; `-r` is not used then).
- `--profile`, `--work-dir`, `--skip-checks` (do not check for Java and Nextflow).
- The pipeline settings in [Optional settings](#optional-settings). Settings equal to their
  default are left off the command.
- `-- -stub` runs the whole graph with placeholder outputs, a quick way to test configs end to end.

### `flair-eval report`
Builds a report from bundles you already have, without running FLAIR:
```bash
flair-eval report out_nov/bundles out_dec/bundles --title "Nov vs Dec" \
  --outdir /private/groups/brookslab/$USER/flair-eval/dd_mm_yy_nov_vs_dec
```
On phoenix it is submitted like `run`, as a smaller job (partition `medium`, 2 CPUs, 8 GB,
12 hours). `--foreground` and `--dry-run` work as for `run`. See
[Reports over many launches](#reports-over-many-launches).

## Config-builder page
<https://hdheath.github.io/flair-eval/builder/> is a single static page: nothing you type is sent
anywhere. You can also open `frontend/builder/index.html` from a clone.

- Add samples, then one card per run. Each card has a FLAIR source picker (release, git branch,
  conda env, container), common setting presets, a label and tags.
- Each sample's genome is `hg38 (lab preset)` (the default; the GTF is then optional) or a custom
  FASTA with its GTF.
- The panel on the right checks everything against the pipeline's schemas as you type and shows
  the exact command. It uses the one-command form when that says exactly the same thing, and
  otherwise says why it needs `samples.csv` and `runs.yaml`.
- The output folder is required. Type your user name: the page quotes paths, so `$USER` would
  not be expanded.
- Where it runs: with `phoenix (default)` and no custom profile, the page writes a
  `flair-eval run …` command, which checks the inputs and submits the Nextflow driver as a SLURM
  job. `This computer (local)` writes `nextflow run … -profile local,<conda|singularity|docker>`,
  and a custom profile writes `nextflow run … -profile <yours>` (on phoenix, run that inside an
  `sbatch` job).
- Copy the command and download the files it needs. Save them in the folder you launch from.

## Running Nextflow directly
```bash
nextflow run hdheath/flair-eval -r main \
  --reads /data/HG002.fq.gz --genome hg38 --flair 2.2.0,3.0.1,dev \
  --outdir /private/groups/brookslab/$USER/flair-eval/dd_mm_yy_hg002
```
```bash
nextflow run hdheath/flair-eval -r main --input samples.csv --runs runs.yaml \
  --outdir /private/groups/brookslab/$USER/flair-eval/dd_mm_yy_name
```
The Nextflow driver runs for the whole comparison. On phoenix start it inside an `sbatch` job
(not on `mustard`), and launch with `umask 002` so the shared caches stay writable for everyone.
Give `-r main` so Nextflow does not reuse an older revision it pulled before.

`nextflow run hdheath/flair-eval --help` lists the parameters; `--help <param>` explains one.

### Profiles
| profile | what it does |
|---|---|
| none (= `phoenix`) | phoenix: `slurm` + `conda` + the shared lab env cache; work dir in `<outdir>/work` |
| `local` | run on this machine, with requests capped at its CPUs and memory |
| `slurm` | SLURM with PRISM's partitions (short/medium/long, picked from each task's time) |
| `conda` | builds a small Python env for evaluation and report, and activates each run's FLAIR env |
| `docker`, `singularity` | container engine on (needed for `container:` runs; evaluation and report use a pandas image) |
| `test` | the bundled test data ([docs/test_data.md](test_data.md)); combine with an executor, e.g. `test,phoenix` or `test,local,conda` |

Off phoenix, combine an executor with an environment, e.g. `-profile local,conda` or
`-profile local,singularity`. Without `conda`, `docker` or `singularity`, evaluation and the
report use the `python3` on `PATH`, which needs pandas and numpy.

`-stub` runs the whole graph with placeholder outputs, which checks configs in seconds:
```bash
nextflow run hdheath/flair-eval -r main -stub -profile local --input samples.csv --runs runs.yaml --outdir stub-out
```

## Reports over many launches
Bundles are self-contained, so a report can combine launches. Add earlier bundles to a new launch
with `--bundles`, or skip compute entirely:
```bash
nextflow run hdheath/flair-eval -r main --report_only \
  --bundles 'old_results/bundles,other/bundles/*/flair-dev' --outdir combined
```
(`flair-eval report` does the same.) If two bundles share a sample and run id, the newest one is
used. With `--bundles` and no `--baseline`, the baseline is the first run of the earliest launch,
so adding a run later keeps the old baseline.

## Other machines and clusters
- Give a profile (`--profile slurm,conda` or `-profile local,conda`): with no profile the pipeline
  uses the phoenix one. `conf/slurm.config` uses phoenix's partition names (short, medium, long);
  adapt it if yours differ.
- Point `--outdir` (and `--work-dir`) at project or scratch storage, not your home folder.
- Set `--env_cache` to shared storage so each FLAIR version is built once and reused.
- The `hg38` preset is phoenix paths: give `--genome` and `--gtf` files instead, and
  `--ref_tss`/`--ref_pas`/`--tusco` if you want the Ends and TUSCO panels.
- `flair-eval run` only submits the driver itself on phoenix; elsewhere start it from a batch job
  or `tmux`.
