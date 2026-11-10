# flair-eval v2: shared contracts

These are the interfaces that the parallel work streams build against. Change one only by
editing this file on the `v2` branch, and say so in the change. Every stream must keep to it.
Design rationale is in [DESIGN.md](DESIGN.md).

## 0. Ownership (who writes which files)
| Stream | Owns |
|---|---|
| **core** | `main.nf`, `nextflow.config`, `conf/`, `nextflow_schema.json`, `modules/**` except `modules/flair_env/`, `subworkflows/**`, `tests/nf-test/` |
| **report** | `bin/structure_metrics.py`, `bin/build_report.py`, `bin/flair_eval_lib/` (shared Python), `assets/report/`, `tests/python/test_structure_metrics.py`, `tests/python/test_report*.py`, `tests/fixtures/bundles/` |
| **adapter** | `bin/flair_adapter.py`, `bin/flair_env.py`, `modules/flair_env/`, `tests/python/test_flair_adapter*.py`, `tests/python/test_flair_env*.py`, `docs/flair_versions.md` |
| **testdata** | `tests/data/**`, `bin/fetch_test_data.sh`, `docs/test_data.md` |
| **frontend** | `frontend/builder/` (static config-builder page), `cli/` (pip package `flair-eval`), `docs/usage_cli.md` |
| **shared, edit here only** | `docs/CONTRACTS.md`, `assets/schema/*.json` |

Placeholder stubs with the contract CLIs exist for `bin/flair_env.py`, `bin/flair_adapter.py`,
`bin/structure_metrics.py`, `bin/build_report.py` and `modules/flair_env/main.nf`. The owning stream
replaces them. Core may rely on their stub blocks / `--help` only.

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
