# Integration notes (core stream, round 2)

Issues found while wiring the merged `v2` streams together and running the first real
end-to-end job. Each item names the owning stream and file. Core did not edit those files.

## Open items for other streams

### Round 3 status
- Fixed upstream in `v2` (a40c81f): the adapter accepts `--flair-bed`, `--platform` and `--library`
  (core now passes all of them; `tests/nf-test/adapter_flags.nf.test` checks QUANTIFY and
  PRODUCTIVITY against the merged adapter's real `--help`), and the FLAIR_ENV stub writes a real
  `source`. The report draws the section 8 panels.
- Fixed in round 3: the adapter writes `quantify_input`; the final run shows `sections.quant.input`
  = `bam` for dev and `reads` for the others, and the report's not-like-for-like note.
- Still open, adapter: `nextflow lint` warns on `modules/flair_env/main.nf:33`
  (`containsKey(it)`, implicit closure parameter, deprecated in strict syntax). Suggested:
  `kinds.find { k -> run.flair instanceof Map && run.flair.containsKey(k) }`.

## FLAIR behaviour seen in the first real run (not fixable in flair-eval; FLAIR is never patched)

### `flair-ted-52826` (TED fork): intron-support database build is ~quadratic in annotated introns
Run: e2e WTC11 PacBio, chr17 TP53 slice, `env: ~/miniforge3/envs/flair-ted-52826`. The command
the adapter runs (from `assemble/logs/03_flair_transcriptome.log`):

```
flair transcriptome -b regions.bam -g genome.fa -f annotation.gtf -t 12 -o flair \
  --sample_name wtc11-pacbio --parallel_mode byregion --generate_map
```

With the full GENCODE v49 basic GTF it sat for 20+ min at one core, 100 % CPU, no I/O, nothing in
`flair.intermediate/` (the same picture as the "never forked" genome-wide runs in the user's
notes). `--log-level INFO` shows it stuck at `building intron support database`. A py-spy profile
(120 s, run as a child so ptrace is allowed) puts ~90 % of samples in:

```
_find_intron (flair/intron_support.py:50-52) <- _add_support (intron_support.py:66)
  -> _find_point_strand (intron_support.py:44) -> overlap (flair/interval_index.py:46-49)
     -> _freeze (interval_index.py:37-38) + ruranges overlaps
```

i.e. every annotated intron triggers a re-freeze of the interval index plus a full overlap
query, so the build grows roughly with the square of the number of introns. The released 3.0.1
and the current `dev` build the same database in seconds. Not a deadlock, just very slow.
Core's mitigation: with `--regions`, PREP_REFERENCE hands FLAIR a reference cut to the regions
(see docs/output.md). Genome-wide TED runs will still pay this cost; it belongs in a FLAIR issue.

### adapter: no `--parallel_mode byregion` for the dev build
In the same run the adapter passed `--parallel_mode byregion` to the TED env but not to
`git: dev` (`flair transcriptome ... --sample_name wtc11-pacbio --generate_map --temp_dir ...`).
CLAUDE.md asks for byregion on every FLAIR. The dev assembly still finished, so this only
affects consistency between runs; worth checking `profile.transcriptome.parallel_flag` detection
in `bin/flair_adapter.py` for dev.

### report: `sections.quant` / `sections.productivity` are not drawn yet
BUNDLE writes both sections (CONTRACTS section 8) and every bundle of the real run has them,
e.g. `flair-3.0.1`: `quant.isoforms_quantified = 253`, `quant.reads_assigned = 26422`,
`productivity.class_counts = {NGO: 28, PRO: 168, PTC: 58}`. The merged `bin/build_report.py`
ignores them (it should, per A16, until the panels exist), so the report has no quant or
productivity panel yet. Concordance and these panels are the report stream's round-2 work.

## First real end-to-end run (2026-11-09)

| | |
|---|---|
| data | LRGASP WTC11 PacBio cDNA (`human_minimap_OUTPUT.sorted.bam`), sliced to chr17:7,400,000-9,400,000: 27,945 alignments, 27,799 reads |
| reference | GRCh38 primary assembly, GENCODE v49 basic |
| runs | `release 2.2.0`, `release 3.0.1`, `git dev` (@389748e), `env flair-ted-52826` (@20c0f47) |
| command | `nextflow run <repo> -profile slurm,conda --input samples.csv --runs runs.yaml --regions regions.bed` |
| location | `/private/groups/brookslab/hdheath/claude-tmp/09_11_26_flair_eval_e2e/` (`inputs/prep.sh`, `run_e2e.sbatch`, `results/`) |
| tasks | 27: 4 FLAIR_ENV, 1 PREP_REFERENCE, 4 × (ASSEMBLE, QUANTIFY, PRODUCTIVITY, EVALUATE_STRUCTURE, BUNDLE), 1 REPORT; no ALIGN (BAM given) |
| env builds | release 2.2.0 2.2 min, release 3.0.1 2.7 min, git dev 3.3 min, env 1 s |
| assemble | 2.2.0 2.0 min, 3.0.1 1.6 min, dev 1.6 min, TED 11.3 min (see above) |
| wall time | launch 1: 21 min, cancelled (dev and TED slow on the full reference); launch 2 (`-resume`, PREP_REFERENCE): 15.3 min; launch 3 (`-resume`, BUNDLE fix): 1.8 min |
| result | `results/report/report.html`: 4 runs, every bundle filled; Overview, Single run and Compare render with no console errors (headless Chromium) |

## Final real end-to-end run (round 4, 2026-11-09, merged `v2` dc8f54d + core)

Same data, reference and region as above; `-profile slurm,conda --regions`; fresh launch (no
`-resume`; FLAIR envs reused from `env_cache`).

### Core fix found on the way: stage the genome `.fai` with the genome
The first round-4 launch failed in `QUANTIFY (wtc11_pacbio/dev-tip)`:

```
flair quantify --manifest quantify/.work/manifest.tsv --genome genome.fa \
  --isoform_bed assemble/isoforms.flair.bed --threads 8 --output quantify/.work/quant ...
  File "flair/flair_quantify.py", line 390, in write_batch_reference
    seq = ''.join(genome.fetch(bed.chrom, block.start, block.end) ...)
KeyError: "sequence 'chr17' not present"
```

`genome.fa` did hold chr17. Core staged the genome without its `.fai`, so dev quantify's pool
workers built the index concurrently and some read it half-written. Re-running the identical
command in the same work dir once the `.fai` existed exited 0. Core now stages `<genome>.fai`
next to the genome for ASSEMBLE, QUANTIFY, PRODUCTIVITY and EVALUATE_STRUCTURE (the original
genome's `.fai` when it exists, else PREP_REFERENCE's). The underlying race is in FLAIR dev
quantify (it should build the index once before forking); worth a FLAIR issue.

### Runs and results (`results_final/`)

| run_index | run id | label | FLAIR | isoforms | FSM | quant input | assemble |
|---|---|---|---|---|---|---|---|
| 0 | flair-2.2.0 | FLAIR 2.2.0 | release 2.2.0 | 280 | 129 | reads | 48 s |
| 1 | flair-3.0.1 | FLAIR 3.0.1 | release 3.0.1 | 253 | 138 | reads | 29 s |
| 2 | dev-tip | dev (tip) | git dev @389748e (the tip; a separate @389748e run would duplicate it) | 143 | 97 | bam | 71 s |
| 3 | dev-425b260 | dev @425b260 | git 425b260, `quantify: false` (its quantify crashes in FLAIR) | 128 | 89 | (none) | 93 s |
| 4 | ted | TED | env flair-ted-52826 @20c0f47 | 134 | 91 | reads | 6.0 min |

- 31 tasks: 5 FLAIR_ENV (all cache hits, 1-7 s), 1 PREP_REFERENCE (14 s), 5 ASSEMBLE, 4 QUANTIFY,
  5 PRODUCTIVITY, 5 EVALUATE_STRUCTURE, 5 BUNDLE, 1 REPORT. No ALIGN (BAM given).
- Wall time of the launch: 11.0 min (SLURM job 46962781), dominated by TED's assemble.
- Report: `results_final/report/report.html` (127 KB). Checked in headless Chromium: every tab
  (Overview, Structure, Agreement, Expression, Runtime, Run detail, Guide) renders with no console
  errors and no NaN/undefined; runs listed in runs.yaml order; Expression shows the
  "Not like for like ... dev (tip) counts from the genome BAM" note; quant, productivity
  ("Coding potential"), agreement and runtime panels are populated.

### Second-run reuse (`results_add/`)
Launch with only `release: 3.0.0` (`inputs/runs_add.yaml`) and
`--bundles results_final/bundles`, into a separate outdir:
- 8 tasks, all for the new run (FLAIR_ENV cold build 2.6 min, PREP_REFERENCE, ASSEMBLE 29 s,
  QUANTIFY, PRODUCTIVITY, EVALUATE_STRUCTURE, BUNDLE) plus REPORT; none of the five earlier runs
  was recomputed. Wall time 6.8 min (SLURM job 46963330).
- `results_add/report/report.html` (131 KB) shows all 6 runs, no console errors.
- FLAIR 3.0.0 and 3.0.1 give identical isoform structures on this region (separate envs,
  different pip freezes; `isoforms.bed` columns 1-12 identical).
- Design note for integration: per A5 core passes this launch's first run as `--baseline`, so in
  an "add a run" launch the new run (3.0.0) becomes the baseline and is listed first. A user adding
  to earlier results probably wants the old baseline; they can pass `--baseline flair-2.2.0`. A
  possible amendment: when `--bundles` is given and no `--baseline`, use the baseline recorded in
  the earlier bundles, or let the report choose by `run_index`.

Paths: `/private/groups/brookslab/hdheath/claude-tmp/09_11_26_flair_eval_e2e/`
(`inputs/runs_final.yaml`, `inputs/runs_add.yaml`, `run_final.sbatch`, `results_final/`,
`results_add/`, `driver_final_<jobid>.log`).

## Round 5: A22 (baseline across launches), core side

- Core passes `--baseline` only when the user gave one, or when there are no `--bundles`
  (then this launch's first run, A5). With `--bundles` and no `--baseline` it passes none.
- Done in round 6, see below.

## Round 6: A24 launch id, and the add-a-run report re-checked

- BUNDLE writes `pipeline.launch` = `workflow.sessionId` (A24). nf-test checks that every bundle
  of one launch carries the same UUID.
- `results_add/` relaunched with `-resume` (SLURM job 46969356, 4.8 min). The FLAIR steps of the
  3.0.0 job (ASSEMBLE, QUANTIFY, PRODUCTIVITY, EVALUATE_STRUCTURE) re-ran as well as BUNDLE and
  REPORT, not only REPORT: the merged adapter's `.fai` guard changed `bin/flair_adapter.py`, and
  Nextflow hashes the bin scripts a task calls. FLAIR_ENV and PREP_REFERENCE came from cache.
- The new `flair-3.0.0` bundle has a launch id; the five `results_final/` bundles have none, so
  the report uses its fallback rule. Its embedded data shows
  `baseline_rule = {rule: earliest_launch, launches: 2}`, baseline `flair-2.2.0`.
- Headless Chromium: header reads "baseline FLAIR 2.2.0, first run of the earliest of 2
  launches"; runs listed FLAIR 2.2.0, FLAIR 3.0.0, FLAIR 3.0.1, dev (tip), dev @425b260, TED
  (baseline first, then run_index, then id); no console errors.
- Report: `/private/groups/brookslab/hdheath/claude-tmp/09_11_26_flair_eval_e2e/results_add/report/report.html`
  (the round-4 version is kept beside it as `report.round4.html`).
