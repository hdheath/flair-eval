# flair-eval: design

## Goal
A Nextflow pipeline that anyone in the lab can run to **compare FLAIR runs** (different
versions, branches or settings) on the same data. Results land in **one self-contained HTML
report** where you can look at any single run, or compare runs against a baseline.

Principles:
- **Isolate FLAIR.** Every run of a sample assembles from the same alignment, and every FLAIR env
  that flair-eval builds gets the same minimap2/samtools/bedtools, so differences come from FLAIR
  itself.
- **Never patch FLAIR.** Version differences are handled in one adapter; FLAIR bugs are
  documented ([flair_versions.md](flair_versions.md)), not fixed.
- **Results are reusable.** Each job writes an immutable, self-contained bundle, and a report can
  be built from any set of bundles, from one launch or many.
- **Phoenix first.** Defaults suit the lab's phoenix cluster; other machines need a profile.

## Vocabulary
| Term | Meaning |
|---|---|
| **sample** | One dataset: reads and/or BAM + genome + GTF. |
| **run** | One FLAIR configuration: which FLAIR (release, git ref, conda env or container) + which arguments. |
| **job** | One run applied to one sample. The pipeline executes runs × samples. |
| **bundle** | A job's results folder: `bundle.json`, a few standard tables, provenance. |
| **report** | `report.html`, built from any set of bundles. |
| **baseline** | The run that every comparison in the report is made against. |

## Pipeline
```
FLAIR_ENV            once per distinct FLAIR source (+ align_flair); cached by release / commit / path
PREP_REFERENCE       with --regions: genome + GTF cut to the regions' chromosomes
ALIGN_SHARED         once per sample without a BAM, with the pinned align_flair
ALIGN_OWN            per job, for runs with own_align: true
ASSEMBLE             per job: flair transcriptome (≥3) or correct + collapse (2.x), via the adapter
QUANTIFY             per job, when run.quantify (reads, or reads extracted from the BAM, or the BAM for dev)
PRODUCTIVITY         per job, when run.productivity (FLAIR's inline call, else predictProductivity)
EVALUATE_STRUCTURE   per job: SQANTI-style classes vs the GTF, end precision, TUSCO, priming, variation
BUNDLE               per job -> <outdir>/bundles/<sample>/<run_id>/
REPORT               all bundles of this launch + --bundles -> <outdir>/report/report.html
```
`--report_only` runs REPORT alone over `--bundles`.

## Code layout
```
main.nf, nextflow.config, nextflow_schema.json   entry, params, profiles (phoenix is the default)
conf/                    base, modules (publishing), phoenix, slurm, test
subworkflows/inputs/     samples.csv / runs.yaml / one-command parsing and validation
subworkflows/flair_eval/ the job graph above
modules/<step>/          one process per step
bin/flair_env.py         builds and caches FLAIR envs (stdlib only)
bin/flair_adapter.py     the only caller of `flair`; maps each version's CLI, normalises outputs (stdlib only)
bin/structure_metrics.py classification + bundle.json (SQANTI3-compatible category names)
bin/build_report.py      the report: bundles + cross-run data, inlined into assets/report/template.html
bin/flair_eval_lib/      shared Python: classify, metrics, sections, ends, tusco, variation, concordance, crossrun
assets/schema/           JSON schemas for samples.csv, runs.yaml, bundle.json
assets/report/           report template, CSS and the viewer (report.js)
assets/demo/bundles/     real bundles for the GitHub Pages demo report
cli/                     pip package `flair-eval`: doctor, init, check, run, report
frontend/builder/        static config-builder page + Pages start page
tests/                   pytest (tests/python, cli/tests), nf-test, test data, report fixtures
```

The interfaces between these parts are in [CONTRACTS.md](CONTRACTS.md).

## Ways in
All map onto the same params and schemas:
1. **One command:** `--reads/--bam`, `--genome`, `--flair 2.2.0,3.0.1,dev`.
2. **Config files:** `samples.csv` + `runs.yaml`.
3. **`flair-eval` CLI:** writes and checks the files, builds the command, and on phoenix submits
   the Nextflow driver as a SLURM job.
4. **Config-builder page:** the same rules in the browser (`logic.js`, tested for parity with the CLI).
5. **`nextflow_schema.json`:** `--help` and validation (it also enables `nf-core pipelines launch`).

## Not done yet
- Report sections that use a sample's own `cage`, `pas`, `junctions` or `spikein_gtf` (the
  columns are parsed but unused).
- A Loci tab and a read 5′-end heatmap (need read alignments in the bundle).
- `--ref_tss`, `--ref_pas` and `--tusco` in the CLI and the builder.
- A SIRV spike-in test set (on a separate branch).
