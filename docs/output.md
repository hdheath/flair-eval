# flair-eval v2: outputs

Everything lands under `--outdir` (default `./flair-eval-results`; on PRISM put it on groups or
warm-archive storage). The layout is fixed by [CONTRACTS.md](CONTRACTS.md) sections 3 and 5.

```
<outdir>/
├── bundles/<sample>/<run_id>/     one folder per job (sample × run)
├── report/report.html             the comparison report
└── pipeline_info/
    ├── execution_trace.txt        Nextflow trace of every task
    ├── report.html                Nextflow resource report
    ├── timeline.html              Nextflow timeline
    └── params.json                the parameters of this launch
```

## Bundles

A bundle is everything the report needs about one FLAIR run on one sample. Bundles are
self-contained (plain copies, no symlinks into `work/`), so you can move, archive or share them
and rebuild a report from any set of them later.

```
bundles/<sample>/<run_id>/
├── bundle.json                    summary numbers (assets/schema/bundle.schema.json)
├── isoforms.bed                   FLAIR isoforms, BED12, name = isoform id
├── isoforms.gtf                   FLAIR isoforms with gene_id + transcript_id
├── counts.tsv                     isoform_id, count, cpm         (when quantify ran)
├── productivity.tsv               isoform_id, cds_start, cds_end, class (when productivity ran)
├── tables/
│   ├── isoforms.tsv               one row per isoform with its structural category
│   ├── categories.tsv             category, subcategory, n
│   └── genes.tsv                  gene_id, gene_name of reference genes (when written)
└── provenance/
    ├── flair_env.json             which FLAIR ran: source, version, commit, env/container
    ├── adapter_<step>.json        exact FLAIR commands, version, exit codes, wall time
    ├── adapter_<step>_skipped.json  present when that FLAIR version cannot do the step
    └── trace.tsv                  resource use of this job's tasks
```

`<step>` is `assemble`, `quantify`, `productivity`, `eval`, or `align` for runs with
`own_align: true`. A FLAIR version that cannot run a step leaves `adapter_<step>_skipped.json`
instead of that step's table.

### What fills `bundle.json`

`bin/structure_metrics.py` writes the structure numbers (`summary`, `categories`,
`fsm_subcategories`, `junctions`, `histograms`, `flair`). The BUNDLE step then adds what only the
pipeline knows:

| key | filled with |
|---|---|
| `pipeline.version` / `pipeline.commit` | the pipeline manifest version and git commit (null for a local checkout) |
| `pipeline.launch` | Nextflow's session id: the same for every bundle of one launch, and kept by `-resume` (A24) |
| `pipeline.created` | when the bundle was written (UTC). The newest bundle wins when the report sees the same (sample, run_id) twice |
| `inputs` | absolute paths of the original `genome`, `gtf`, `reads`, `bam`, `regions`, plus `platform`, `library` and `reference_cut_to_regions` (below); other keys written by structure_metrics are kept |
| `label`, `tags` | from the run, when structure_metrics left them empty |
| `compute` | from this job's Nextflow task traces (below) |

### `compute` and `provenance/trace.tsv`

`compute` is the cost of the FLAIR run itself: the job's `assemble`, `quantify` and
`productivity` tasks (CONTRACTS section 7, A12). Alignment (shared or `own_align`) and the
evaluation step are not counted, so runs on the same sample compare fairly.

| field | meaning |
|---|---|
| `cpu_hours` | Σ realtime × %CPU over the counted tasks |
| `peak_rss_gb` | largest peak RSS of any counted task (GiB) |
| `wall_minutes` | Σ realtime of the counted tasks |

Values come from each task's `.command.trace`, which Nextflow writes because `nextflow.config`
turns tracing on. A value is `null` when a trace could not be read.

`trace.tsv` has one row per task of the job (all steps, including `align` for own_align runs and
the evaluation), with the columns of CONTRACTS section 7, A15, in raw units:

| column | meaning |
|---|---|
| `task_id` | the task's work-directory hash, as Nextflow prints it (`ab/cdef12`) |
| `name` | process and tag, e.g. `FLAIR_EVAL:ASSEMBLE (sample/run)` |
| `status`, `exit` | `COMPLETED` / `FAILED` and the exit code |
| `realtime`, `duration` | ms; duration runs from task submission to completion |
| `%cpu` | percent of one CPU |
| `peak_rss`, `rss` | bytes |

### `sections`

BUNDLE adds two report sections (CONTRACTS sections 8 and 9) when their files exist, by calling
the one shared implementation, `python3 -m flair_eval_lib.sections <bundle_dir>` (A19):

| key | from | content |
|---|---|---|
| `sections.quant` | `counts.tsv` | `isoforms_quantified`, `reads_assigned`, `cpm_hist`, and `input` (`reads` or `bam`: how this FLAIR version quantified, copied from the adapter log, A21) |
| `sections.productivity` | `productivity.tsv` | `class_counts`: isoforms per productivity class |

Both count only isoforms listed in `tables/isoforms.tsv`, so a `--regions` run counts the same
isoforms in every panel. The CPM bin edges are defined in `bin/flair_eval_lib/sections.py`.

### Other `bundle.json` keys set by BUNDLE

- `run_index` (A17): the run's position in runs.yaml or `--flair` order, starting at 0. The
  report lists runs as baseline first, then by `run_index`.

### `--regions` and the reference FLAIR sees

With `--regions`, a PREP_REFERENCE step copies the genome and GTF cut down to the chromosomes
named in the BED (whole chromosomes, so all coordinates stay the same). Assembly, quantify,
productivity and the structure evaluation use that copy; the shared alignment still uses the
full genome. This keeps a small test region fast: on the first real run some FLAIR builds
otherwise re-read a 3 GB genome over shared storage for every region, or spent 20+ minutes
building whole-genome annotation indexes. `inputs.reference_cut_to_regions` records it.

## Report

`report/report.html` is built by `bin/build_report.py` from this launch's bundles plus any
`--bundles` folders. It is a single file with all data and code inlined, so it opens offline and
can be emailed or hosted. `--report_title` sets its heading and `--baseline` the run that
differences are computed against.

To rebuild the report later from bundles only, without running FLAIR:

```bash
nextflow run hdheath/flair-eval --report_only \
  --bundles 'results_a/bundles,results_b/bundles/*/flair-dev' --outdir results_combined
```

## What is not published

FLAIR's raw outputs (alignments, read maps, isoform FASTA, intermediate files) stay in the
Nextflow `work/` directory. The shared alignment is not copied either. Clean `work/` when you no
longer need `-resume`.

Re-running into the same `--outdir` overwrites `pipeline_info/` and the report, and overwrites
the files of any bundle that is produced again. Bundles of other runs are left alone. A file that
an earlier launch wrote but this one does not (for example `counts.tsv` after setting
`quantify: false`) is not deleted, so use a fresh `--outdir` when you change a run's steps.
