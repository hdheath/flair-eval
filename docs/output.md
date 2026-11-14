# Outputs and the report

Everything lands under `--outdir` (required; on phoenix put it on groups or warm-archive storage).
The layout is fixed by [CONTRACTS.md](CONTRACTS.md) sections 3 and 5.

```
<outdir>/
├── bundles/<sample>/<run_id>/     one folder per job (sample × run)
├── report/report.html             the comparison report
├── work/                          Nextflow work dir (phoenix profile; FLAIR's raw outputs)
└── pipeline_info/
    ├── execution_trace.txt        Nextflow trace of every task
    ├── report.html                Nextflow resource report
    ├── timeline.html              Nextflow timeline
    ├── params.json                the parameters of this launch
    └── launch.sbatch, nextflow-<job>.log   the driver job (flair-eval run on phoenix)
```

## Reading the report

`report/report.html` is one file with all data and code inlined: it opens offline, and can be
emailed or hosted. It is built from this launch's bundles plus any `--bundles` folders
([report-only mode](usage.md#reports-over-many-launches)). See the
[demo report](https://hdheath.github.io/flair-eval/demo/report.html).

- **One sample at a time.** The Sample picker drives every tab; nothing is pooled across samples.
- **Every comparison is against the baseline run.** It defaults to the first run (`--baseline`
  changes it), and the control bar switches it. Each run keeps its colour when you switch.
- **Overlay or side by side.** With two or more runs, charts overlay; the Layout toggle draws one
  panel per run. Every chart has a Table toggle with its numbers.
- **Missing inputs hide panels.** A panel whose data is missing says what it needs (for
  example "needs refTSS sites (--genome hg38, or --ref_tss)"). When only some runs have the data,
  the panel is drawn for those and names the others.

| Tab | Answers |
|---|---|
| What differs | Where do the runs differ most? Headline metrics side by side, each run's biggest differences from the baseline (at most one per metric family), and every metric in one table ("All metrics") |
| Structure | How do the isoforms compare with the reference annotation? Structural categories (FSM, ISM, NIC, NNC, …, each defined in one line), isoform counts, lengths, exons per isoform |
| Recovery | Which reference genes and transcripts are found? TUSCO benchmark, genes detected, gained and lost against the baseline, FSM isoforms with exact ends |
| Diversity | How many isoforms per gene, how many are coding or NMD, and do genes vary by ends, splicing or both? |
| Agreement | How similar are the runs? Shared and unique isoforms, expression correlation, 5′ and 3′ end agreement, TUSCO genes vs the baseline |
| Ends | How close are isoform 5′ and 3′ ends to refTSS and PolyA_DB sites (precision only) and to the reference transcript's ends, and which 3′ ends sit on genomic A-runs (internal-priming risk)? |
| Quantification | How many reads were assigned, and how many isoforms have reads? |
| Compute | Wall time, CPU time and peak memory of each run's FLAIR steps |

**Read counts are not always like for like across versions:** dev counts alignments in the
genome BAM, while 2.2.0 and 3.0.1 count reads from the reads file. The Quantification tab says so
when the runs being compared counted differently.

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
│   ├── genes.tsv                  gene_id, gene_name of reference genes (when written)
│   └── ends.tsv                   per isoform: 5'/3' ends, distances to end references, poly(A) motif, downstream A
└── provenance/
    ├── flair_env.json             which FLAIR ran: source, version, commit, env/container
    ├── adapter_<step>.json        exact FLAIR commands, version, exit codes, wall time
    ├── adapter_<step>_skipped.json  present when that FLAIR version cannot do the step
    └── trace.tsv                  resource use of this job's tasks
```

`<step>` is `assemble`, `quantify`, `productivity`, or `align` for runs with `own_align: true`. A FLAIR version that cannot run a step leaves `adapter_<step>_skipped.json`
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
| `inputs` | absolute paths of the original `genome`, `gtf`, `reads`, `bam`, `regions`, `ref_tss`, `ref_pas`, `tusco` (null when not used), plus `platform`, `library` and `reference_cut_to_regions` (below); other keys written by structure_metrics are kept |
| `label`, `tags` | from the run, when structure_metrics left them empty |
| `compute` | from this job's Nextflow task traces (below) |

### `compute` and `provenance/trace.tsv`

`compute` is the cost of the FLAIR run itself: the job's `assemble`, `quantify` and
`productivity` tasks (CONTRACTS A12). Alignment (shared or `own_align`) and the
evaluation step are not counted, so runs on the same sample compare fairly.

| field | meaning |
|---|---|
| `cpu_hours` | Σ realtime × %CPU over the counted tasks |
| `peak_rss_gb` | largest peak RSS of any counted task (GiB) |
| `wall_minutes` | Σ realtime of the counted tasks |

Values come from each task's `.command.trace`, which Nextflow writes because `nextflow.config`
turns tracing on. A value is `null` when a trace could not be read.

`trace.tsv` has one row per task of the job (all steps, including `align` for own_align runs and
the evaluation), with the columns of CONTRACTS A15, in raw units:

| column | meaning |
|---|---|
| `task_id` | the task's work-directory hash, as Nextflow prints it (`ab/cdef12`) |
| `name` | process and tag, e.g. `FLAIR_EVAL:ASSEMBLE (sample/run)` |
| `status`, `exit` | `COMPLETED` / `FAILED` and the exit code |
| `realtime`, `duration` | ms; duration runs from task submission to completion |
| `%cpu` | percent of one CPU |
| `peak_rss`, `rss` | bytes |

### `sections`

BUNDLE adds two report sections (CONTRACTS section 8) when their files exist, by calling
the one shared implementation, `python3 -m flair_eval_lib.sections <bundle_dir>` (A19):

| key | from | content |
|---|---|---|
| `sections.quant` | `counts.tsv` | `isoforms_quantified`, `reads_assigned`, `cpm_hist`, and `input` (`reads` or `bam`: how this FLAIR version quantified, copied from the adapter log, A21) |
| `sections.productivity` | `productivity.tsv` | `class_counts`: isoforms per productivity class |

Both count only isoforms listed in `tables/isoforms.tsv`, so a `--regions` run counts the same
isoforms in every panel. The CPM bin edges are defined in `bin/flair_eval_lib/sections.py`.

### End, TUSCO, priming and variation sections (CONTRACTS section 10)

EVALUATE_STRUCTURE passes the sample's end references and TUSCO list (`--ref_tss`, `--ref_pas`,
`--tusco`, or the genome preset's, A26). These sections cover main chromosomes only (chr1-22, X, Y).

| key | needs | content |
|---|---|---|
| `sections.ends` | `ref_tss` and/or `ref_pas` | per reference (`refs`: file, md5, sites), the isoforms within each window of `windows` bp of a site, signed distance histograms, FSM distances to their reference transcript, the same numbers by category, and the poly(A) motif count (A29). The half without its file is null. |
| `sections.tusco` | `tusco` | SQANTI3 TUSCO benchmark: TP/PTP/FP/FN genes, sensitivity, precision, `per_gene` labels, plus SQANTI3's transcript-level counts in `transcripts` (A30). Rules: `bin/flair_eval_lib/tusco.py`. |
| `sections.priming` | the genome | isoforms whose 20 nt downstream of the 3' end are >= 60% A (internal-priming risk), and the A-fraction histogram (A31) |
| `sections.variation` | nothing | genes with >= 2 isoforms that differ by ends only, splicing only, both, or neither (A32). Rules: `bin/flair_eval_lib/variation.py`. |

`tables/ends.tsv` columns (A28): `isoform_id category subcategory strand tss tes d_reftss d_pas
d_ref_tss d_ref_tes polya_motif down_a_frac`. `tss`/`tes` are 0-based end bases. Distances are
genomic and signed in transcript orientation: negative means the isoform end is upstream (5′) of
the site. A site farther than 10 kb leaves the cell empty; `d_ref_*` are for FSM and ISM only.

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

## What is not published

FLAIR's raw outputs (alignments, read maps, isoform FASTA, intermediate files) stay in the
Nextflow `work/` directory. The shared alignment is not copied either. Clean `work/` when you no
longer need `-resume`.

Re-running into the same `--outdir` overwrites `pipeline_info/` and the report, and overwrites
the files of any bundle that is produced again. Bundles of other runs are left alone. A file that
an earlier launch wrote but this one does not (for example `counts.tsv` after setting
`quantify: false`) is not deleted, so use a fresh `--outdir` when you change a run's steps.
