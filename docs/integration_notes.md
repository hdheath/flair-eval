# Integration notes (core stream, round 2)

Issues found while wiring the merged `v2` streams together and running the first real
end-to-end job. Each item names the owning stream and file. Core did not edit those files.

## Open items for other streams

### adapter: `bin/flair_adapter.py` does not accept `--platform`, `--library` or `--flair-bed` yet
CONTRACTS section 7, A8 and A9. Core passes each of these flags only when
`flair_adapter.py <step> --help` lists it (`modules/adapter_opts.nf`), so nothing breaks today,
but the values are silently not used.

```
$ bin/flair_adapter.py assemble --bam x.bam --genome g.fa --gtf a.gtf --out o --platform pacbio
flair_adapter.py: error: unrecognized arguments: --platform pacbio
```

Once the adapter adds them, core picks them up with no change.

### adapter: `modules/flair_env/main.nf` stub writes `"source":"stub"`
CONTRACTS section 7, A11 says the stub must write `source` as the run's `flair` key
(`release|git|env|container`). The current stub block prints `"source":"stub"`, which is not
allowed by the bundle schema. It only matters in `-stub` runs (core's EVALUATE_STRUCTURE stub
writes its own `flair` block), so nothing fails, but the stub `flair_env.json` copied into
`provenance/` is not contract-shaped. Suggested stub line:

```groovy
def src = (run.flair.keySet() as List).find { k -> k in ['release', 'git', 'env', 'container'] }
printf '{"run_id":"%s","source":"%s",...}' '${run.id}' '${src}'
```

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
