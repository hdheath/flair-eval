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
│   └── categories.tsv             category, subcategory, n
└── provenance/
    ├── flair_env.json             which FLAIR ran: source, version, commit, env/container
    ├── adapter_<step>.json        exact FLAIR commands, version, exit codes, wall time
    ├── adapter_<step>_skipped.json  present when that FLAIR version cannot do the step
    └── trace.tsv                  resource use of this job's tasks
```

`<step>` is `assemble`, `quantify`, `productivity`, `eval`, or `align` for runs with
`own_align: true`.

### What fills `bundle.json`

`bin/structure_metrics.py` writes the structure numbers (`summary`, `categories`,
`fsm_subcategories`, `junctions`, `histograms`, `flair`). The BUNDLE step then adds what only the
pipeline knows:

| key | filled with |
|---|---|
| `pipeline.version` / `pipeline.commit` | the pipeline manifest version and git commit (null for a local checkout) |
| `pipeline.created` | when the bundle was written (UTC). The newest bundle wins when the report sees the same (sample, run_id) twice |
| `inputs` | absolute paths of `genome`, `gtf`, `reads`, `bam`, `regions`, plus `platform` and `library` |
| `label`, `tags` | from the run, when structure_metrics left them empty |
| `compute` | from this job's Nextflow task traces (below) |

### `compute` and `provenance/trace.tsv`

`compute` is the cost of the FLAIR run itself: the job's `align` (own_align runs only),
`assemble`, `quantify` and `productivity` tasks. The shared per-sample alignment and the
evaluation step are not counted, so runs on the same sample compare fairly.

| field | meaning |
|---|---|
| `cpu_hours` | Σ realtime × %CPU over the counted tasks |
| `peak_rss_gb` | largest peak RSS of any counted task (GiB) |
| `wall_minutes` | Σ realtime of the counted tasks |

Values come from each task's `.command.trace`, which Nextflow writes because `nextflow.config`
turns tracing on. A value is `null` when a trace could not be read.

`trace.tsv` has one row per task of the job, with columns `step, counted, realtime_ms, pct_cpu,
peak_rss_bytes, peak_vmem_bytes, rchar, wchar`.

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
