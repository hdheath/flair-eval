# flair-eval v2: design

Status: draft, 2026-11-09. Decisions were agreed in the setup session. The old
multi-assembler pipeline stays on `wip/2026-11-09-snapshot`.

## Goal
A Nextflow pipeline that other people can run to **compare FLAIR runs**, meaning different
versions, branches or settings, on the same data. Results land in **one self-contained HTML
report**. In the report you can open any single run, or tick several runs to compare them
side by side.

## Vocabulary
| Term | Meaning |
|---|---|
| **sample** | One dataset: reads or BAM + genome + GTF, plus optional orthogonal tracks. |
| **run** | One FLAIR configuration: *which FLAIR* (version source) + *which arguments*. |
| **job** | One run applied to one sample. The pipeline executes runs × samples. |
| **bundle** | A job's results folder: `metrics.json` + a few standard tables + `provenance.json`. Bundles are immutable and reusable. |
| **report** | `report.html`, built from any set of bundles, whether from one launch or many over time. |

## Inputs
`samples.csv` has one row per sample:
- Required: `sample`, `genome`, `gtf`, and either `reads` or `bam`.
- Optional: `platform` (pacbio/ont), `library` (cDNA/dRNA).
- Optional tracks for later report sections: `cage`, `pas`, `junctions` (STAR SJ.out.tab),
  `spikein_gtf`.

`runs.yaml` has one entry per run. Each run gives exactly one FLAIR source:
```yaml
runs:
  - id: flair-3.0.1
    flair: { release: "3.0.1" }                     # bioconda/pip; the env is built and cached
    transcriptome_args: ""
  - id: dev-tip
    flair: { git: "https://github.com/BrooksLabUCSC/flair", ref: "dev" }   # pip-installed, cached per commit
    transcriptome_args: "--stringent"
  - id: ted-local
    flair: { env: "/path/to/conda/env" }            # or { container: "docker://..." }
    transcriptome_args: "--ted --ted_chemistry pacbio_cdna"
    label: "TED (local build)"                       # optional display name
    tags: { branch: ted-port, setting: ted }          # optional free-form tags for report filters
```

## Pipeline steps
1. **Resolve FLAIR environments.** There is one env per distinct source. git refs resolve to a commit SHA, so the
   cache key is the exact code. Each env records `flair --version` and its pip freeze.
2. **Shared alignment, once per sample.** A pinned FLAIR + minimap2 runs `flair align` (+ `correct` when a
   run needs it). Every run sees the same BAM. A run can set `own_align: true` to test alignment
   itself. Pre-aligned BAMs skip this step.
3. **Per job: assemble.** FLAIR ≥3 runs `flair transcriptome`; FLAIR 2.x runs `correct` → `collapse`. A small
   **version adapter** maps a run onto the right CLI and normalises outputs to
   `isoforms.bed/.gtf/.fa` + `read_map`. Known 3.0.1 quirks are documented in the adapter.
   Then `quantify` and `predictProductivity` run when the version supports them.
4. **Per job: evaluate → bundle.** Structure vs GTF (v1), described below. Compute cost comes from the Nextflow trace.
5. **Report.** Collects bundles from this launch plus any `--bundles <dir|list>` from earlier ones.
   It writes `report.html` with the data embedded, so it opens offline, can be emailed, or can be hosted on private.prism.

FLAIR source code is never patched. Version differences are handled in the adapter.

## Report v1: isoform structure vs GTF
The classifier is in-house and SQANTI-compatible. It reuses `bin/evaluation/flair_structural.py`
(FSM/ISM/NIC/NNC/genic/intergenic/antisense/fusion; FSM subcategories CIM/FSM5/FSM3/FSMTE; 50-bp end window).
Category names match SQANTI3.

Each run's page shows:
- isoforms, genes and loci detected
- category counts and fractions
- isoforms per gene
- transcript and exon length distributions
- mono- vs multi-exon split
- junction novelty (known / novel canonical / novel non-canonical)
- reference transcripts recovered (FSM/CIM)

Comparison view:
- Tick runs and filter by tag, version or sample.
- The same panels are shown overlaid or side by side, plus a run × metric table with differences
  against a chosen baseline run.

Later sections are only shown when their input exists:
- TSS/TTS end accuracy (CAGE/PAS, `ted_end_precision.py`)
- run-vs-run concordance (shared/unique isoforms, changed genes, example loci)
- quantification / spike-in / CDS / runtime

## Front end for new users
There are four ways in, all mapping onto the same params and schemas (see CONTRACTS.md §1):
1. **One command, no files.**
   `nextflow run hdheath/flair-eval --reads r.fq.gz --genome g.fa --gtf a.gtf --flair 2.2.0,3.0.1,dev --outdir out`.
2. **`nextflow_schema.json`.** Gives `--help`, input validation with plain-English errors, a web
   launch form via `nf-core pipelines launch`, and a Seqera Platform launchpad.
3. **Config-builder page** (`frontend/builder/`, GitHub Pages). It has a form for samples, plus run cards
   with a version picker and presets. It validates against `assets/schema/*.json`, shows the
   exact command, and downloads `samples.csv` / `runs.yaml`.
4. **pip CLI** (`cli/`, package `flair-eval`):
   - `flair-eval init` is an interactive wizard that writes the config files.
   - `flair-eval run` wraps `nextflow run`.
   - `flair-eval report` gives report-only mode.
   - `flair-eval check` validates configs.

## Layout (new, lean)
```
main.nf                      # one entry; -entry REPORT rebuilds the report from bundles only
nextflow.config              # profiles: test, slurm (PRISM), local, docker/singularity/conda
conf/                        # per-profile resources
modules/  flair_env/  align/  assemble/  quantify/  productivity/  evaluate_structure/  report/
bin/      flair_adapter.py  structure_metrics.py (from flair_structural.py)  build_report.py
assets/   report template (HTML + JS, charting lib inlined)
tests/    pytest for bin/, nf-test for modules; data/ = tiny test sets
docs/     DESIGN.md, usage.md, output.md
```

## Test data (`-profile test`)
- **Human region.** About 1–2 Mb of hg38 + its GENCODE slice + about 20k subsampled public reads
  (PacBio and ONT, from public sources such as ENCODE or PacBio downloads).
- **SIRV.** SIRV genome + annotation + spike-in reads (known truth).
- Both are small enough to commit or fetch from a fixed URL. CI runs `-profile test` on both.

## Milestones
1. **Skeleton.** samples.csv + runs.yaml parsing. Env resolution for `env` and `release`. Shared align → `flair transcriptome` →
   bundle with structure metrics. Report v1 for a single run. Runs on the human test region.
2. **Comparison report.** Multi-run toggles, baseline differences, tag filters, `--bundles` accumulation.
3. **Version breadth.** `git` sources cached by SHA. 2.x collapse adapter. quantify + productivity. SIRV test set.
4. **Usability.** Containers, README/usage docs, GitHub Actions CI on `-profile test`.
5. **Extra sections.** End accuracy, concordance, quant/CDS/runtime.

## Open points (decide during M1)
- Exact pinned versions for the shared align step.
- Whether `report.html` embeds the full per-isoform tables or only aggregates. Aggregates
  keep the file small; per-isoform tables enable drill-down.
- Which public reads to use for the test region (needs a source with a clear licence).
