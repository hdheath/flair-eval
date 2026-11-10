# flair-eval v2: shared contracts

These are the interfaces that the parallel work streams build against. Change one only by
editing this file on the `v2` branch, and say so in the change. Every stream must keep to it.
Design rationale is in [DESIGN.md](DESIGN.md).

## 0. Ownership (who writes which files)
| Stream | Owns |
|---|---|
| **core** | `main.nf`, `nextflow.config`, `conf/`, `nextflow_schema.json`, `modules/**` except `modules/flair_env/`, `subworkflows/**`, `tests/nf-test/`, `docs/output.md`, `README.md` "Usage" section, `.github/workflows/ci.yml` |
| **report** | `bin/structure_metrics.py`, `bin/build_report.py`, `bin/flair_eval_lib/` (shared Python), `assets/report/`, `tests/python/test_structure_metrics.py`, `tests/python/test_report*.py`, `tests/fixtures/bundles/` |
| **adapter** | `bin/flair_adapter.py`, `bin/flair_env.py`, `modules/flair_env/`, `tests/python/test_flair_adapter*.py`, `tests/python/test_flair_env*.py`, `docs/flair_versions.md` |
| **testdata** | `tests/data/**`, `bin/fetch_test_data.sh`, `docs/test_data.md` |
| **frontend** | `frontend/builder/` (static config-builder page), `cli/` (pip package `flair-eval`), `docs/usage_cli.md`, `.github/workflows/pages.yml` |
| **integration (parent session, on `v2`)** | `docs/CONTRACTS.md`, `assets/schema/*.json`, `.gitignore`, `docs/DESIGN.md`, merges |

The placeholders are gone: the core, adapter, report and frontend branches were merged into `v2` on
2026-11-09. Before starting new work, merge `v2` into your stream branch. **§7 (Amendments) overrides
any earlier section it touches.**

`legacy/` holds the old multi-assembler pipeline. Read it, port from it, and never edit it.
Reusable pieces: `legacy/bin/evaluation/flair_structural.py` (SQANTI-style classifier) and
`legacy/bin/evaluation/ted_end_precision.py` (later end-accuracy section).

## 1. Inputs
### samples.csv  (schema: `assets/schema/samples.schema.json`)
| column | req | notes |
|---|---|---|
| `sample` | yes | `^[A-Za-z0-9][A-Za-z0-9._-]*$`, unique |
| `reads` | one of reads/bam | FASTQ/FASTA(.gz) |
| `bam` | one of reads/bam | sorted + indexed (`.bai` next to it) |
| `genome` | yes | FASTA |
| `gtf` | yes | reference annotation |
| `platform` | no | `pacbio` \| `ont` |
| `library` | no | `cDNA` \| `dRNA` |
| `cage`, `pas`, `junctions`, `spikein_gtf` | no | reserved for later sections; parsed and passed through, unused in v1 |

### runs.yaml  (schema: `assets/schema/runs.schema.json`)
```yaml
runs:
  - id: flair-3.0.1                # required, ^[A-Za-z0-9][A-Za-z0-9._-]*$, unique
    flair:                         # required: EXACTLY ONE of release | git | env | container
      release: "3.0.1"
      # git: "https://github.com/BrooksLabUCSC/flair"   ref: "dev"  (branch, tag or SHA)
      # env: "/abs/path/to/conda/env"
      # container: "docker://..." | "/abs/path/image.sif"
    transcriptome_args: ""         # optional, appended verbatim to flair transcriptome/collapse
    quantify: true                 # optional, default true
    productivity: true             # optional, default true
    own_align: false               # optional, default false (true = this run aligns its own BAM)
    label: "FLAIR 3.0.1"           # optional display name (default = id)
    tags: { version: "3.0.1" }     # optional str→str, used as report filters
```

### One-command mode (no files)
`--reads`/`--bam` + `--genome` + `--gtf` [+ `--sample`, `--platform`, `--library`] makes one
in-memory sample (sample name defaults to the reads/BAM basename without its extensions).
`--flair` takes a comma-separated list of tokens. Each token becomes one run:

| token | run.flair | run.id |
|---|---|---|
| `3.0.1` or `release:3.0.1` | `{release: "3.0.1"}` | `flair-3.0.1` |
| `dev`, `main`, `git:<ref>` | `{git: "https://github.com/BrooksLabUCSC/flair", ref: <ref>}` | `flair-<ref>` |
| `git:<url>@<ref>` | `{git: <url>, ref: <ref>}` | `flair-<ref>` |
| `env:/path` | `{env: "/path"}` | `env-<basename>` |
| `container:<uri>` | `{container: <uri>}` | `ctr-<basename-without-tag>` |

`--flair_args "..."` is applied to every run built from `--flair`. Within one source, `--input`
and the quick sample flags are mutually exclusive, and so are `--runs` and `--flair`.

### Pipeline params (names are fixed; core writes them into `nextflow_schema.json`)
| param | default | meaning |
|---|---|---|
| `input` | null | samples.csv |
| `runs` | null | runs.yaml |
| `reads`, `bam`, `genome`, `gtf`, `sample`, `platform`, `library` | null | one-command sample |
| `flair` | null | one-command runs (see above) |
| `flair_args` | `""` | args for one-command runs |
| `outdir` | `./flair-eval-results` | output root |
| `bundles` | null | comma-separated dirs/globs of earlier bundles to include in the report |
| `baseline` | null | run id to compute report differences against (default = first run) |
| `report_only` | false | skip compute; build the report from `bundles` only |
| `align_flair` | `"3.0.1"` | FLAIR release used for the shared align step |
| `regions` | null | optional BED to restrict assembly/eval (test profile uses it) |
| `report_title` | `"FLAIR comparison"` | report heading |
| `env_cache` | `"${launchDir}/flair-env-cache"` | where FLAIR envs are built/cached (PRISM: set on groups storage) |

## 2. FLAIR environment + adapter CLIs (adapter stream implements; core calls)
**`modules/flair_env/main.nf`** (adapter owns it) exports process `FLAIR_ENV`.
- input: `val(run)`, the run map from runs.yaml.
- output: `tuple val(run_id), path("flair_env.json")`.

`flair_env.json`:
```json
{"run_id":"…","source":"release|git|env|container","requested":"3.0.1 | <url>@<ref> | <path> | <uri>",
 "resolved_commit":"<sha or null>","flair_version":"<output of flair --version>",
 "conda_prefix":"<abs path or null>","container":"<uri or null>",
 "python":"3.x.y","pip_freeze_sha256":"…","created":"ISO-8601"}
```
Core runs FLAIR tasks inside that env: `conda_prefix` via the `conda` directive, or `container`.
envs are cached under `params.env_cache` and keyed by
`release` / `resolved_commit` / path. If the same key is requested twice, it is built once.

**`bin/flair_adapter.py`** runs *inside* the FLAIR env and is the only thing that calls `flair`:
```
flair_adapter.py align        --reads R --genome G [--gtf A] --threads N --out DIR
flair_adapter.py assemble     --bam B --reads R --genome G --gtf A --args "…" --threads N --out DIR [--regions BED]
flair_adapter.py quantify     --isoforms-fa F --reads R --sample S --threads N --out DIR
flair_adapter.py productivity --isoforms-bed B --genome G --gtf A --out DIR
```
- It detects the FLAIR version and picks the CLI: ≥3 uses `transcriptome`, 2.x uses `correct` + `collapse`.
- It never patches FLAIR.
- An unsupported step on a version writes `skipped.json` and exits 0.

Normalised outputs in `DIR`:
| step | files |
|---|---|
| align | `aligned.bam`, `aligned.bam.bai` |
| assemble | `isoforms.bed` (BED12, name = isoform id), `isoforms.gtf`, `isoforms.fa`, `read_map.tsv` (read_id, isoform_id; optional) |
| quantify | `counts.tsv` (isoform_id, count, cpm) |
| productivity | `productivity.tsv` (isoform_id, cds_start, cds_end, class) |
| every step | `adapter_log.json`: exact commands, FLAIR version, exit codes, wall time |

Gene ids: `isoforms.gtf` must carry `gene_id` + `transcript_id`. When FLAIR gives no gene, use
`novel_gene_<n>`.

## 3. Bundle (core writes it, report reads it)
Each job = (sample, run). Layout: `<outdir>/bundles/<sample>/<run_id>/`
```
bundle.json        # schema below (assets/schema/bundle.schema.json)
isoforms.bed       # copied from assemble
isoforms.gtf
tables/isoforms.tsv          # one row per isoform (columns below)
tables/categories.tsv        # category, subcategory, n
counts.tsv / productivity.tsv  (when produced)
provenance/flair_env.json  provenance/adapter_*.json  provenance/trace.tsv (this job's tasks)
```
`structure_metrics.py --isoforms-bed B --isoforms-gtf T --gtf A --sample S --run-json RUN.json --env-json flair_env.json --out BUNDLE_DIR [--regions BED]`
writes `bundle.json` + `tables/*`. RUN.json is the run map as JSON. The core stream calls it.

`tables/isoforms.tsv` columns, in this order:
`isoform_id gene_id chrom start end strand n_exons length category subcategory ref_transcript ref_gene`
- `category` ∈ `FSM ISM NIC NNC genic genic_intron antisense intergenic fusion`
- `subcategory` (FSM only) ∈ `CIM FSM5 FSM3 FSMTE`, otherwise empty

`bundle.json` (v1):
```json
{
  "schema": "flair-eval/bundle/1",
  "sample": "…", "run_id": "…", "label": "…", "tags": {"k": "v"},
  "flair": {"source": "release", "requested": "3.0.1", "version": "3.0.1", "commit": null},
  "pipeline": {"version": "…", "commit": "…", "created": "ISO-8601"},
  "inputs": {"genome": "…", "gtf": "…", "reads": "…", "bam": "…", "regions": null},
  "summary": {"isoforms": 0, "genes": 0, "loci": 0, "mono_exon": 0, "multi_exon": 0,
              "ref_transcripts_total": 0, "ref_transcripts_fsm": 0, "ref_transcripts_cim": 0,
              "ref_genes_detected": 0},
  "categories": {"FSM": 0, "ISM": 0, "NIC": 0, "NNC": 0, "genic": 0, "genic_intron": 0,
                 "antisense": 0, "intergenic": 0, "fusion": 0},
  "fsm_subcategories": {"CIM": 0, "FSM5": 0, "FSM3": 0, "FSMTE": 0},
  "junctions": {"total": 0, "known": 0, "novel_canonical": 0, "novel_noncanonical": 0},
  "histograms": {
    "isoforms_per_gene":   {"edges": [1,2,3,4,5,6,8,11,16,21], "counts": []},
    "transcript_length":   {"edges": [0,250,500,750,1000,1500,2000,3000,4000,6000,10000], "counts": []},
    "exon_length":         {"edges": [0,50,100,150,200,300,500,1000,2000,5000], "counts": []},
    "exons_per_isoform":   {"edges": [1,2,3,4,5,6,8,11,16,21,31], "counts": []}
  },
  "compute": {"cpu_hours": null, "peak_rss_gb": null, "wall_minutes": null},
  "sections": {}
}
```
- Histogram `edges` are fixed (as above) so runs compare bin by bin. The last bin is open-ended
  (≥ last edge), so `len(counts) == len(edges)`.
- Later sections (end accuracy, concordance, quant) add keys under `sections` only. The report
  must ignore section keys it doesn't know, and hide any panel whose section is missing.
- `compute` is filled by core from the Nextflow trace. Null is allowed.

## 4. Report
`build_report.py --bundles DIR_OR_GLOB [...] --out report.html [--title T] [--baseline RUN_ID]`
- It walks every `bundle.json` under the given paths. If two bundles share (sample, run_id), the newest `pipeline.created` wins.
- The output is one HTML file with no network access needed: JS/CSS inlined, data embedded as JSON.
- Views: an overview table (run × metric), a single-run page, and a comparison view (tick runs;
  filter by sample/tag/version; differences against the baseline).
- Default target is under 5 MB for 50 bundles. Per-isoform tables are not embedded in v1.

## 5. Outputs
```
<outdir>/bundles/<sample>/<run_id>/…      (section 3)
<outdir>/report/report.html
<outdir>/pipeline_info/{execution_trace.txt, report.html, timeline.html, params.json}
```

## 6. Conventions
- Python 3.10+, stdlib + pandas/numpy only in `bin/` (structure metrics, report). The adapter is stdlib only,
  because it runs inside arbitrary FLAIR envs.
- Every bin script has `--help`, exits non-zero on bad input, and is covered by pytest.
- Commit to your own stream branch (`v2-<stream>`). Never push or merge unless the user says to.

## 7. Amendments — round 1 (2026-11-09, after the first merge)
These settle the gaps the streams reported. Each item is tagged with the stream(s) that act on it.

**Inputs and run ids**
- A1 [core, frontend] **`--flair` token rules.** `cli/src/flair_eval/tokens.py` is the reference, and core's Groovy parser must match it exactly:
  - A bare token is a release only if it starts with a digit. `dev` and `main` mean git refs. Any other bare word is an error that suggests `git:<ref>`.
  - In `git:<url>@<ref>`, the URL is everything before the *last* `@`.
  - In ids, any character outside `[A-Za-z0-9._-]` becomes `-` (`git:feature/x` gives `flair-feature-x`).
  - Container ids drop the `:tag` and `@digest`. A local `.sif` gives `ctr-<basename>` (e.g. `ctr-flair.sif`).
- A2 [core, frontend] **Duplicate run ids** (from tokens or from runs.yaml) are an error. The message names both sources and tells the user to switch to runs.yaml with explicit ids.
- A3 [core, testdata, frontend] **Relative paths in samples.csv** are resolved against the CSV's directory first, then the launch directory.
- A4 [core, frontend, all docs] **Quote `--flair_args` with `=`:** `--flair_args='--stringent'`. Nextflow mis-parses the space form `--flair_args "--stringent"` (it sets `flair_args=true` plus a stray param). Docs and examples must use the `=` form. Core should also fail when `flair_args` is the Boolean `true`, explaining how to quote it.
- A5 [core] **Default baseline:** the first run in runs.yaml / `--flair` order. Core always passes `--baseline`. `build_report.py` falls back to name order only when no baseline is given (report-only mode over old bundles).

**Adapter CLI (§2)**
- A6 [core, adapter] **`assemble --reads` is optional.** A BAM-only sample skips QUANTIFY unless the version can quantify from the BAM (see A7).
- A7 [core, adapter] **`quantify` gains optional `--bam B --genome G --isoforms-bed I`.** Core always passes the ones it has. The adapter picks what the version needs and writes `skipped.json` if a required one is missing.
- A8 [core, adapter] **`assemble` may also write `isoforms.flair.bed`** (FLAIR's native BED, all columns). Core declares it as an optional output and stages it into QUANTIFY / PRODUCTIVITY. The adapter takes it via an optional `--flair-bed` flag rather than finding it through symlinks.
- A9 [core, adapter] **All adapter commands accept optional `--platform` and `--library`.** The adapter records them in `adapter_log.json` and may map them to version flags, documenting each mapping in `docs/flair_versions.md`.
- A10 [adapter, report] **`counts.tsv`, `productivity.tsv` and `read_map.tsv` have one header line,** with the column names in §2.
- A11 [core, adapter] **`FLAIR_ENV`** gets `{id, flair}` from the first run that uses each source; `run_id` in `flair_env.json` is that env key's id. It has label `flair_env`, runs on the host or base env (never inside a FLAIR env), and needs network, git and mamba/conda. Its stub must write `source` as one of `release|git|env|container`, taken from the run's `flair` key; `stub` is not allowed.

**Bundle (§3)**
- A12 [core, report] **Who writes what in `bundle.json`.** `structure_metrics.py` writes the structural fields. Core's BUNDLE step fills `pipeline`, `inputs` and `compute`. `compute` counts only this job's FLAIR tasks (assemble, quantify, productivity), not the shared align or the evaluation step.
- A13 [core, report] **`structure_metrics.py` gains optional `--genome G --reads R --bam B`.** Core always passes the genome. Without the genome, `junctions.motif_checked = false` and novel junctions are reported only as `junctions.novel`.
- A14 [report] **New optional keys in the bundle schema:** `summary.ref_genes_total`, `junctions.novel`, `junctions.motif_checked`.
- A15 [core] **Provenance:**
  - An adapter `skipped.json` is stored as `provenance/adapter_<step>_skipped.json`.
  - `provenance/trace.tsv` holds the Nextflow trace rows for this job's tasks, with columns `task_id name status exit realtime %cpu peak_rss rss duration`.

**Report sections (round 2, see §8)**
- A16 [report] Sections listed in `sections` follow §8. The report must keep working with bundles that don't have them.

**Repo housekeeping (done on `v2` by integration)**
- `.gitignore` has a `!tests/data/` exception, and the stale tracked legacy `.nf-test/` folder has been removed.
- `-entry REPORT` from DESIGN.md is replaced by `--report_only`.
- `nf-test` config stays at `tests/nf-test/nf-test.config`.

## 8. Round-2 report sections (`sections.<name>` in bundle.json)
Each section is optional and has its own panel. It is computed only when its inputs exist.
- `sections.concordance`: computed by the report from run pairs over a sample's bundles, not stored per bundle.
  `build_report.py` reads `isoforms.bed` beside each `bundle.json` for this. It covers:
  - shared / unique isoforms, matched by intron chain (mono-exon matched by overlap ≥50% on the same strand);
  - Jaccard;
  - genes gained or lost relative to the baseline;
  - the top 20 loci where runs differ most (gene, coordinates, counts per run), shown as a table with copyable coordinates.
- `sections.quant`: added by core from `counts.tsv`. It holds `{"isoforms_quantified", "reads_assigned", "cpm_hist": {edges, counts}}`. The report shows CPM correlation between runs that share a sample by reading `counts.tsv` beside the bundle.
- `sections.productivity`: added by core from `productivity.tsv`. It holds `{"class_counts": {class: n}}`.
- `sections.compute`: the report shows the existing `compute` block as a runtime panel (CPU-h, peak RSS, wall time per run). It needs no new key.

## 9. Amendments — round 2 (2026-11-10, after the round-2 merge)
- A17 [core, report] **Run order.** `bundle.json` gains optional `run_index` (integer ≥0): the run's position in runs.yaml or `--flair` order. Core's BUNDLE writes it. The report orders runs as baseline first, then `run_index`, then id. Bundles without it (old ones) sort after indexed runs, by id.
- A18 [report, core] **`tables/genes.tsv`** (`gene_id gene_name`, one header line) is an optional output of `structure_metrics.py`. Core publishes it in the bundle when present.
- A19 [core, report] **One implementation of sections.** `python3 -m flair_eval_lib.sections <bundle_dir>` is the single place that builds `sections.quant` and `sections.productivity` from the TSVs. Core's BUNDLE calls it instead of keeping its own copy.
- A20 [adapter] **dev `quantify` must use `--genome`, not `-g`.** dev builds 57040d2 and 425b260 reject `-g`. Choose flags from `--help`, as is already done elsewhere. When a version supports `--parallel_mode`, pass `--parallel_mode byregion` to dev `transcriptome` too, never bychrom.
- A21 [core, adapter, report] **"Reads assigned" depends on how a version quantifies:**
  - 2.x, 3.0.1 and TED count from reads;
  - dev counts from the genome BAM.

  The adapter records `quantify_input: reads|bam` in `adapter_log.json`, core copies it to `sections.quant.input`, and the report notes it wherever reads-assigned numbers are compared.
- A22 [core, report] **Baseline across launches.**
  - When `--baseline` is not given and `--bundles` is non-empty, core does **not** pass `--baseline`.
  - `build_report.py` then picks the baseline as the run with the lowest `run_index` among the bundles with the oldest `pipeline.created`, i.e. the first run of the first launch. Ties are broken by id.
  - With no `--bundles`, A5 still applies: the first run of this launch.
  - The report shows which run is the baseline and why ("first run of the earliest launch"), and the user can switch it in the page.
- A23 [core, adapter] **Stage `<genome>.fai`** next to the genome for every FLAIR and evaluation step. dev quantify workers race to build it otherwise (seen as `KeyError: sequence 'chr17' not present`). That is a FLAIR dev bug; document it in `docs/flair_versions.md`.
- A24 [core, report] **Launch id.**
  - `bundle.json` gains optional `pipeline.launch`, which is Nextflow's `workflow.sessionId`. Core's BUNDLE writes it, and it is the same for every bundle from one launch. A `-resume`d launch reuses its session id, so it counts as the same launch.
  - The report groups bundles into launches by `pipeline.launch` when present. It falls back to the A22 time-and-`run_index` heuristic for bundles without it.
  - "Earliest launch" means the launch whose earliest `pipeline.created` is oldest.

## 10. Amendments — phoenix defaults (2026-11-10)
- A25 [core, adapter, frontend] **phoenix is the default platform.** All users run on the lab's phoenix cluster (UCSC PRISM).
  - No `-profile` = `standard` = `conf/phoenix.config`: `conf/slurm.config` + conda + `params.env_cache = /private/groups/brookslab/flair-eval/env_cache` (one cache for the lab) + `conda.cacheDir = /private/groups/brookslab/flair-eval/conda_cache` (Nextflow's own envs, also shared) + `workDir = <outdir>/work`. Launch with umask 002. `phoenix` is an alias. `local`, `docker`, `singularity` and `test` are unchanged.
  - `--outdir` is required everywhere (no default). Without it the run stops before any task.
  - `params.genomes` (hidden) maps preset keys to `{fasta, gtf, note}`; today only `hg38` = GRCh38 primary assembly + GENCODE v49 basic from the lab reference folders. A `genome` value (flag or samples.csv cell) equal to a key becomes the preset FASTA, and an empty `gtf` becomes the preset GTF. Presets are explicit: no genome is still an error.
  - `flair_env.py resolve` sets umask 002 so a shared cache stays group-writable. The cache root is setgid, group brookslab.
  - Front ends assume phoenix: no profile, `--outdir` required, `hg38` offered as the genome, and the CLI submits the Nextflow driver as a SLURM job.
