# flair-eval: interface contracts

The interfaces between the pipeline (Nextflow), the FLAIR adapter (`bin/flair_adapter.py`,
`bin/flair_env.py`), the evaluation and report code (`bin/structure_metrics.py`,
`bin/build_report.py`, `bin/flair_eval_lib/`) and the front ends (`cli/`, `frontend/builder/`).
Machine-readable parts are in `assets/schema/*.json`. Architecture: [DESIGN.md](DESIGN.md).
User guide: [usage.md](usage.md); outputs: [output.md](output.md).

Change an interface here and in the schemas together, in the same commit. Rules carry tags
`A1`–`A40`, `A42`–`A44` from the rounds in which they were agreed; code comments cite them, and
[section 7](#7-rule-index) lists where each one now lives.

`legacy/` holds the old multi-assembler pipeline. Read it and port from it; never edit it.

## 1. Inputs

### samples.csv (schema: `assets/schema/samples.schema.json`)
| column | req | notes |
|---|---|---|
| `sample` | yes | `^[A-Za-z0-9][A-Za-z0-9._-]*$`, unique |
| `reads` | one of reads/bam | FASTQ/FASTA(.gz) |
| `bam` | one of reads/bam | coordinate-sorted, index (`.bai`/`.csi`) next to it |
| `genome` | yes | FASTA, or a genome preset key (A25) |
| `gtf` | yes, unless `genome` is a preset | reference annotation |
| `platform` | no | `pacbio` \| `ont` |
| `library` | no | `cDNA` \| `dRNA` |
| `junctions` | no | splice junctions for this sample: STAR `SJ.out.tab` or BED (BED6+, read count in the score column). EVALUATE_STRUCTURE scores every run's junctions against it (section 13, A42) |
| `cage`, `pas`, `spikein_gtf` | no | reserved: parsed, checked and passed through; no step uses them yet. Truth sets for end accuracy come from `--truth` (section 12), not from `cage`/`pas` |

Relative paths are resolved against the CSV's directory first, then the launch directory (A3).

### runs.yaml (schema: `assets/schema/runs.schema.json`)
```yaml
runs:
  - id: flair-3.0.1                # required, ^[A-Za-z0-9][A-Za-z0-9._-]*$, unique (A2)
    flair:                         # required: EXACTLY ONE of release | git | env | container
      release: "3.0.1"
      # git: "https://github.com/BrooksLabUCSC/flair"   ref: "dev"  (branch, tag or SHA)
      # env: "/abs/path/to/conda/env"
      # container: "docker://..." | "/abs/path/image.sif"
    transcriptome_args: ""         # optional, appended verbatim to flair transcriptome / 2.x collapse
    quantify: true                 # optional, default true
    productivity: true             # optional, default true
    own_align: false               # optional, default false (true = this run aligns its own BAM)
    label: "FLAIR 3.0.1"           # optional display name (default = id)
    tags: { version: "3.0.1" }     # optional str→str, used as report filters
```

### One-command mode (no files)
`--reads`/`--bam` + `--genome` [+ `--gtf`, `--sample`, `--platform`, `--library`] makes one
in-memory sample; the sample name defaults to the reads/BAM basename without its extensions.
`--flair` takes comma-separated tokens; each becomes one run:

| token | run.flair | run.id |
|---|---|---|
| `3.0.1` or `release:3.0.1` | `{release: "3.0.1"}` | `flair-3.0.1` |
| `dev`, `main`, `git:<ref>` | `{git: "https://github.com/BrooksLabUCSC/flair", ref: <ref>}` | `flair-<ref>` |
| `git:<url>@<ref>` | `{git: <url>, ref: <ref>}` | `flair-<ref>` |
| `env:/path` | `{env: "/path"}` | `env-<basename>` |
| `container:<uri>` | `{container: <uri>}` | `ctr-<basename-without-tag>` |

Token rules (A1). `cli/src/flair_eval/tokens.py` is the reference; the pipeline's Groovy parser
and `frontend/builder/logic.js` must match it exactly (`cli/tests/test_builder_parity.py` and
`tests/nf-test/functions.nf.test` check this):
- A bare token is a release only if it starts with a digit. `dev` and `main` are git refs. Any
  other bare word is an error that suggests `git:<ref>`.
- In `git:<url>@<ref>`, the URL is everything before the *last* `@`.
- In ids, any character outside `[A-Za-z0-9._-]` becomes `-` (`git:feature/x` → `flair-feature-x`).
- Container ids drop the registry path, `:tag` and `@digest`; a local `.sif` keeps its file name
  (`ctr-flair.sif`).
- Duplicate run ids, from tokens or runs.yaml, are an error that names both sources and suggests
  runs.yaml with explicit ids (A2).

`--flair_args` applies to every run built from `--flair`. Docs and front ends write it with `=`
(`--flair_args='--stringent'`), because Nextflow reads `--flair_args "--stringent"` as
`flair_args=true` plus a stray param; the pipeline stops with a quoting hint when
`flair_args` is the Boolean `true` (A4). `--input` and the one-command sample flags are mutually
exclusive, and so are `--runs` and `--flair`.

### Pipeline params (`nextflow.config`, `nextflow_schema.json`)
| param | default | meaning |
|---|---|---|
| `input`, `runs` | null | samples.csv, runs.yaml |
| `reads`, `bam`, `genome`, `gtf`, `sample`, `platform`, `library` | null | one-command sample |
| `flair` | null | one-command runs (tokens above) |
| `flair_args` | `""` | args for one-command runs |
| `outdir` | none: **required** (A25) | output root; the run stops before any task without it |
| `bundles` | null | comma-separated dirs/globs of earlier bundles to add to the report |
| `baseline` | null | run id for report differences (section 4) |
| `report_only` | false | skip compute; build the report from `bundles` only (needs `bundles`) |
| `report_title` | `"FLAIR comparison"` | report heading |
| `align_flair` | `"3.0.1"` | FLAIR release used for the shared align step |
| `regions` | null | BED restricting assembly and evaluation |
| `ref_tss`, `ref_pas`, `tusco` | null; filled from a genome preset (A26) | end references and TUSCO list (section 10). Hidden in `nextflow_schema.json` until the CLI and builder expose them (A35) |
| `truth` | null | truth sheet CSV (section 12, A36). Hidden like `ref_tss` (A35) |
| `env_cache` | `"${launchDir}/flair-env-cache"`; phoenix: the shared lab cache (A25) | where FLAIR envs are built and cached |
| `genomes` | the `hg38` preset (hidden) | genome presets (section 9) |

## 2. FLAIR environment and adapter

### `FLAIR_ENV` (`modules/flair_env/main.nf`)
- Input: `{id, flair}` from the first run that uses each distinct source. `run_id` in
  `flair_env.json` is that env key's id (A11).
- Output: `tuple val(run_id), path("flair_env.json")`.
- Label `flair_env`; runs on the host or base env, never inside a FLAIR env; needs network, git
  and conda/mamba. Its stub writes `source` as one of `release|git|env|container`, from the run's
  `flair` key (A11).
- Envs are cached under `params.env_cache`, keyed by `release` / `resolved_commit` / path; a key
  requested twice is built once (file lock).

```json
{"run_id":"…","source":"release|git|env|container","requested":"3.0.1 | <url>@<ref> | <path> | <uri>",
 "resolved_commit":"<sha or null>","flair_version":"<output of flair --version>",
 "conda_prefix":"<abs path or null>","container":"<uri or null>",
 "python":"3.x.y","pip_freeze_sha256":"…","created":"ISO-8601"}
```
FLAIR tasks run inside that env: `conda_prefix/bin` first on `PATH` (and the `conda` directive
with `-profile conda`), or the `container`.

### `bin/flair_adapter.py`
Runs *inside* the FLAIR env (stdlib only) and is the only thing that calls `flair`. It never
patches FLAIR.
```
flair_adapter.py align        --reads R --genome G [--gtf A] --threads N --out DIR
flair_adapter.py assemble     --bam B [--reads R] --genome G --gtf A --args "…" --threads N --out DIR [--regions BED] [--sample S]
flair_adapter.py quantify     --isoforms-fa F [--reads R] --sample S [--bam B --genome G --isoforms-bed I --flair-bed FB] --threads N --out DIR
flair_adapter.py productivity --isoforms-bed B --genome G --gtf A [--flair-bed FB] --out DIR
flair_adapter.py detect       [--json]
```
The four step commands also take `--flair-bin` and `--keep-work`, and all of them take `--platform` and
`--library` (A9).
- The adapter picks the CLI by probing `--help` (not by version number): 2.x uses `correct` +
  `collapse`, ≥3 uses `transcriptome`. Option names come from the same help text (A20: e.g.
  dev `quantify` gets `--genome`, never `-g`; `--parallel_mode byregion` wherever offered, never
  bychrom).
- `assemble --reads` is optional (A6). For a version that needs reads (2.x `collapse`) and has
  none, the adapter extracts the BAM's primary alignments with `samtools fastq -F 0x904`.
- `quantify` gets everything the pipeline has (A7) and uses what the version needs. A version
  that counts from reads but has no `--reads` gets the BAM's primary alignments as reads (A6b;
  `adapter_log.json` records `quantify_reads: "extracted_from_bam"`).
- `--platform` and `--library` are recorded in `adapter_log.json` and mapped to version flags
  only where FLAIR's own help ties a flag to the library type (A9); each mapping is in
  [flair_versions.md](flair_versions.md#5---platform-and---library-a9).
- `assemble` and BAM-based `quantify` run `samtools faidx` first when the genome has no `.fai`
  (second guard for A23).
- An unsupported step writes `skipped.json` and exits 0.

Normalised outputs in `DIR`:
| step | files |
|---|---|
| align | `aligned.bam`, `aligned.bam.bai` |
| assemble | `isoforms.bed` (BED12, name = isoform id), `isoforms.gtf`, `isoforms.fa`, `read_map.tsv` (read_id, isoform_id; optional), `isoforms.flair.bed` (FLAIR's own BED, all columns; optional, A8) |
| quantify | `counts.tsv` (isoform_id, count, cpm) |
| productivity | `productivity.tsv` (isoform_id, cds_start, cds_end, class) |
| every step | `adapter_log.json`: exact commands, FLAIR version, exit codes, wall time; `quantify_input: reads\|bam` for quantify (A21) |

- `isoforms.gtf` carries `gene_id` + `transcript_id`; with no gene from FLAIR, `novel_gene_<n>`.
- `counts.tsv`, `productivity.tsv` and `read_map.tsv` have exactly one header line with the column
  names above (A10).
- The pipeline stages `isoforms.flair.bed` into QUANTIFY and PRODUCTIVITY and passes it as
  `--flair-bed` (A8).
- The pipeline stages `<genome>.fai` next to the genome for every FLAIR and evaluation step:
  dev quantify workers otherwise race to build it (A23).

## 3. Bundle

One bundle per job (sample × run): `<outdir>/bundles/<sample>/<run_id>/`.
```
bundle.json                    schema: assets/schema/bundle.schema.json
isoforms.bed, isoforms.gtf     copied from assemble
counts.tsv, productivity.tsv   when produced
tables/isoforms.tsv            one row per isoform (columns below)
tables/categories.tsv          category, subcategory, n
tables/genes.tsv               gene_id, gene_name; optional (A18)
tables/ends.tsv                section 10 (A28)
tables/end_accuracy.tsv        section 12 (A39)
provenance/flair_env.json
provenance/adapter_<step>.json, adapter_<step>_skipped.json (A15)
provenance/trace.tsv           this job's Nextflow trace rows (A15)
provenance/end_accuracy/<call>/{precision_recall_summary.tsv, judge.log}   section 12 (A37)
```

### `structure_metrics.py` (EVALUATE_STRUCTURE)
```
structure_metrics.py --isoforms-bed B --isoforms-gtf T --gtf A --sample S --run-json RUN.json \
  --env-json flair_env.json --out BUNDLE_DIR [--regions BED] [--genome G] [--reads R] [--bam B] \
  [--ref-tss BED] [--ref-pas BED] [--tusco TSV] [--junction-file SJ]
```
Writes `bundle.json` and `tables/*`. RUN.json is the run map as JSON. The pipeline always passes
the genome (A13); without it `junctions.motif_checked = false` and novel junctions are reported
only as `junctions.novel`. `--ref-tss/--ref-pas/--tusco` add the section 10 output (A27).
`sections.structure_qc` is always written; `--junction-file` adds its `junction_support` key (section 13).

`tables/isoforms.tsv` columns, in order:
`isoform_id gene_id chrom start end strand n_exons length category subcategory ref_transcript ref_gene`
- `category` ∈ `FSM ISM NIC NNC genic genic_intron antisense intergenic fusion` (SQANTI3 names)
- `subcategory` (FSM only) ∈ `CIM FSM5 FSM3 FSMTE`, otherwise empty

### `bundle.json`
```json
{
  "schema": "flair-eval/bundle/1",
  "sample": "…", "run_id": "…", "label": "…", "tags": {"k": "v"}, "run_index": 0,
  "flair": {"source": "release", "requested": "3.0.1", "version": "3.0.1", "commit": null},
  "pipeline": {"version": "…", "commit": "…", "created": "ISO-8601", "launch": "<session id>"},
  "inputs": {"genome": "…", "gtf": "…", "reads": "…", "bam": "…", "regions": null,
             "ref_tss": null, "ref_pas": null, "tusco": null, "platform": null, "library": null,
             "truth": null, "truth_sets": [],
             "reference_cut_to_regions": false},
  "summary": {"isoforms": 0, "genes": 0, "loci": 0, "mono_exon": 0, "multi_exon": 0,
              "ref_transcripts_total": 0, "ref_transcripts_fsm": 0, "ref_transcripts_cim": 0,
              "ref_genes_detected": 0, "ref_genes_total": 0},
  "categories": {"FSM": 0, "ISM": 0, "NIC": 0, "NNC": 0, "genic": 0, "genic_intron": 0,
                 "antisense": 0, "intergenic": 0, "fusion": 0},
  "fsm_subcategories": {"CIM": 0, "FSM5": 0, "FSM3": 0, "FSMTE": 0},
  "junctions": {"total": 0, "known": 0, "novel_canonical": 0, "novel_noncanonical": 0,
                "novel": 0, "motif_checked": true},
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
- Histogram `edges` are fixed so runs compare bin by bin. The last bin is open-ended (≥ last
  edge), so `len(counts) == len(edges)`.
- **Who writes what (A12):** `structure_metrics.py` writes the structural fields and section 10.
  The BUNDLE step fills `pipeline`, `inputs`, `compute`, `run_index` and `sections.quant` /
  `sections.productivity` (section 8), and `label`/`tags` when left empty.
- `compute` counts only this job's FLAIR tasks (assemble, quantify, productivity), not the
  shared align or the evaluation step (A12). Null is allowed.
- Optional keys (A14): `summary.ref_genes_total`, `junctions.novel`, `junctions.motif_checked`.
- `run_index` (A17): the run's position in runs.yaml or `--flair` order, from 0.
- `pipeline.launch` (A24): Nextflow's `workflow.sessionId`, the same for every bundle of one
  launch; `-resume` keeps it.
- New report data goes under `sections` only. The report ignores section keys it does not know
  and hides a panel whose section is missing (A16).
- `provenance/trace.tsv` columns (A15): `task_id name status exit realtime %cpu peak_rss rss duration`.

## 4. Report

```
build_report.py --bundles DIR_OR_GLOB [...] --out report.html [--title T] [--baseline RUN_ID]
```
- Reads every `bundle.json` under the given paths. If two bundles share (sample, run_id), the
  newest `pipeline.created` wins.
- One HTML file that needs no network: JS/CSS inlined, data embedded as JSON. Per-isoform tables
  are not embedded; target under 5 MB for 50 bundles.
- **Baseline.** Without `--bundles`, the pipeline passes this launch's first run (A5). With
  `--bundles` and no `--baseline`, it passes none (A22), and `build_report.py` picks the run with
  the lowest `run_index` in the earliest launch, ties broken by id. Launches are grouped by
  `pipeline.launch` (A24), else by time and `run_index`; "earliest" = the launch whose earliest
  `pipeline.created` is oldest. When no bundle has a `run_index` it falls back to name order. The embedded
  `baseline_rule` (`given | earliest_launch | name_order`) says which rule applied, and the page
  can switch the baseline.
- **Run order (A17):** baseline first, then `run_index`, then id; bundles without `run_index`
  sort after indexed runs, by id.
- Cross-run data is computed per sample by `build_report.py` (sections 8 and 11).
- The page itself: section 11.

## 5. Outputs
```
<outdir>/bundles/<sample>/<run_id>/…      (section 3)
<outdir>/report/report.html
<outdir>/pipeline_info/{execution_trace.txt, report.html, timeline.html, params.json}
```
FLAIR's raw outputs stay in the Nextflow work dir.

## 6. Conventions
- Python 3.10+. `bin/` uses the stdlib plus pandas/numpy only. The one exception is the vendored
  end judge under `modules/end_accuracy/resources/usr/bin/` (section 12), which imports scipy; the
  END_ACCURACY process has its own environment for it. The adapter and `flair_env.py`
  are stdlib only, because they run inside arbitrary FLAIR envs.
- Every bin script has `--help`, exits non-zero on bad input, and is covered by pytest.
- FLAIR is never patched: version differences live in the adapter, bugs are documented in
  [flair_versions.md](flair_versions.md).
- `-entry REPORT` from early designs is `--report_only`. The nf-test config is
  `tests/nf-test/nf-test.config`.

## 7. Rule index
Where each tagged rule lives now.

| rules | section |
|---|---|
| A1–A4 (tokens, duplicate ids, relative paths, `--flair_args`) | 1 |
| A5 (default baseline) | 4 |
| A6, A6b, A7–A11 (adapter CLI, BAM-only, FLAIR_ENV) | 2 |
| A12–A15, A17, A18, A24 (bundle contents, provenance, run order, launch id) | 3 (A17, A24 also 4) |
| A16, A19, A21 (sections, one implementation, reads vs BAM counts) | 8 |
| A20, A23 (adapter flags, genome `.fai`) | 2 |
| A22 (baseline across launches) | 4 |
| A25 (phoenix defaults) | 9 |
| A26–A32, A35 (end references, TUSCO, priming, variation; as built) | 10 |
| A33, A34 (cross-run data, report rules) | 11 |
| A36–A40 (truth sheet, end judge, end accuracy section, table, report) | 12 |
| A42–A44 (structure QC section, junction-file support, report) | 13 |

## 8. Sections from counts and productivity; run concordance
Each section is optional, computed only when its inputs exist, and the report keeps working
without it (A16).
- One implementation (A19): `python3 -m flair_eval_lib.sections <bundle_dir>` builds
  `sections.quant` and `sections.productivity` from the bundle's TSVs; BUNDLE calls it. Both count
  only isoforms listed in `tables/isoforms.tsv`.
- `sections.quant`: `{"isoforms_quantified", "reads_assigned", "cpm_hist": {edges, counts}, "input"}`.
- `sections.productivity`: `{"class_counts": {class: n}}`.
- **Reads assigned depends on how a version quantifies (A21):** 2.x, 3.0.1 and TED count from
  reads; dev counts alignments in the genome BAM. The adapter records `quantify_input`, BUNDLE
  copies it to `sections.quant.input`, and the report notes it wherever reads-assigned numbers
  are compared.
- **Concordance** is computed by `build_report.py` from run pairs over one sample's bundles
  (`by_sample.<s>.concordance`, `.quant`; `bin/flair_eval_lib/concordance.py`), not stored per
  bundle. It reads `isoforms.bed`, `tables/isoforms.tsv`, `tables/genes.tsv` and `counts.tsv`
  beside each `bundle.json`: shared/unique isoforms matched by intron chain (mono-exon by ≥50%
  same-strand overlap), Jaccard, genes gained/lost against the baseline, CPM correlation.
- The runtime panel shows the existing `compute` block; no new key.

## 9. Phoenix defaults (A25)
All lab users run on phoenix (UCSC PRISM), so it is the default platform.
- No `-profile` = `standard` = `conf/phoenix.config` (`phoenix` is an alias): `conf/slurm.config`
  + conda + `params.env_cache = /private/groups/brookslab/flair-eval/env_cache` (one FLAIR env
  cache for the lab) + `conda.cacheDir = /private/groups/brookslab/flair-eval/conda_cache`
  (Nextflow's own envs) + `workDir = <outdir>/work`. Launch with umask 002.
- `--outdir` is required everywhere (no default).
- `params.genomes` (hidden) maps preset keys to `{fasta, gtf, ref_tss, ref_pas, tusco, note}`;
  today only `hg38` = GRCh38 primary assembly + GENCODE v49 basic from the lab reference folders.
  A `genome` value (flag or samples.csv cell) equal to a key becomes the preset FASTA, and an
  empty `gtf` becomes the preset GTF. No genome is still an error. The same presets are in
  `GENOME_PRESETS` in `cli/src/flair_eval/phoenix.py` and `frontend/builder/logic.js`; tests check
  all three match.
- `flair_env.py resolve` sets umask 002 so the shared cache stays group-writable; the cache root
  is setgid, group brookslab.
- Front ends assume phoenix: no profile, `--outdir` required, `hg38` offered as the genome, and
  the CLI submits the Nextflow driver as a SLURM job.

## 10. End references, TUSCO, priming and variation
### Inputs (A26)
- `params.ref_tss` (BED6), `params.ref_pas` (BED6, one row per cleavage site, `[pos-1,pos)`) and
  `params.tusco` (SQANTI3 `tusco_human.tsv` format), all null by default.
- The `hg38` preset brings `/private/groups/brookslab/flair-eval/reference/{end_refs/refTSS_v4.1.hg38.bed,
  end_refs/PolyA_DB_v4.1.main.hg38.bed, tusco/tusco_human.tsv}`; each folder has a `SOURCE.tsv`.
  When the genome is the preset, a null param is filled from it; an explicit param wins. A
  sample with a custom genome gets no preset refs.
- EVALUATE_STRUCTURE stages whichever exist and passes them to `structure_metrics.py`. A missing
  preset file warns and turns off its section; a missing explicit file stops the run before any
  task (A35).
- Code: `bin/flair_eval_lib/ends.py`, `tusco.py`, `variation.py` (pandas/numpy). Main chromosomes
  only (chr1–22, X, Y, with or without `chr`) (A27).
- The genome is always passed, so `tables/ends.tsv`, `sections.priming` and `sections.variation`
  are always written (A35).

### `tables/ends.tsv` (A28, A35)
One row per main-chromosome isoform, one header line:
`isoform_id category subcategory strand tss tes d_reftss d_pas d_ref_tss d_ref_tes polya_motif down_a_frac`
- `tss`/`tes`: the isoform's 5′/3′ end bases (0-based, strand-aware).
- Distances are genomic (not spliced), signed in transcript orientation as isoform minus site:
  negative = the isoform end is upstream (5′). This equals SQANTI3 `diff_to_TTS` and has the
  opposite sign to SQANTI3 `diff_to_TSS`.
  - `d_reftss` / `d_pas`: to the nearest same-strand refTSS / PolyA_DB site (0 inside a site
    interval); empty when no site lies within 10 kb.
  - `d_ref_tss` / `d_ref_tes`: to the matched reference transcript's ends; FSM/ISM only; uncapped.
- `polya_motif`: the first hexamer, in priority order, from {AATAAA, ATTAAA, AGTAAA, TATAAA,
  CATAAA, GATAAA, AATATA, AATACA, AATAGA, AAAAAG, ACTAAA, AAGAAA} found in the 50 nt upstream of
  the 3′ end, or empty.
- `down_a_frac`: fraction of A in the 20 nt of genome just downstream of the 3′ end
  (strand-aware). Internal-priming risk = ≥ 0.6.
- An unstranded isoform gets `tss`/`tes` and empty cells elsewhere, and is left out of section counts.

### `sections.ends` (A29, A35; needs `ref_tss` and/or `ref_pas`)
```json
{"refs": {"tss": {"file": "refTSS_v4.1.hg38.bed", "md5": "…", "sites": 0}, "pas": {…}},
 "windows": [0, 10, 25, 50, 100, 200, 300, 500],
 "tss": {"n": 0, "no_site": 0, "within": [n per window], "dist_hist": {"edges": E, "counts": []}},
 "tes": {"n": 0, "no_site": 0, "within": [...], "dist_hist": {"edges": E, "counts": []}},
 "fsm_tss_hist": {"edges": E, "counts": []}, "fsm_tes_hist": {"edges": E, "counts": []},
 "by_category": {"FSM": {"n": 0, "tss_within_50": 0, "tes_within_50": 0}, "…": {}},
 "polya_motif": {"n_with_motif": 0, "n": 0}}
```
- The half without its file is null. `by_category` lists all 9 categories, null where a
  reference is missing. `polya_motif` is null without a genome.
- E = `[-1000, -200, -50, -10, 0, 10, 50, 200, 1000]`; bins are `[edge_i, edge_i+1)` with open
  bins below the first and above the last edge, so `len(counts) == len(E) + 1`; exactly 0 falls
  in `[0, 10)`.
- `within[i]` = isoforms with |d| ≤ windows[i]; precision = within / n. `n` includes `no_site`
  (stranded isoforms with no site within 10 kb). With database references this is precision
  only: the report never shows recall from them.
- `fsm_*_hist` bin `d_ref_tss` / `d_ref_tes` of FSM isoforms.
- `refs.*.sites` counts sites after removing duplicate position+strand rows.

### `sections.tusco` (A30, A35; needs `tusco`)
`{"genes_total": 46, "genes_in_regions": n, "TP": n, "PTP": n, "FP": n, "FN": n, "sensitivity": x,
"precision": x, "per_gene": {gene_id: "TP|PTP|FP|FN"}}` plus `precision_redundant`,
`positive_detection_rate`, `transcripts {n, TP, PTP, FP}`, `id_type`, `file`, `md5`, `gene_names`.
- Definitions follow SQANTI3 `src/utilities/report_qc/TUSCO_report.R` (commit 589dfd4);
  `tusco.py` documents the rules and cites the lines.
- `TP/PTP/FP/FN` are gene counts summing to `genes_in_regions`. FP = a TUSCO gene with
  transcripts but none FSM/ISM. `sensitivity` = TP genes / genes in scope; `precision` = TP
  transcripts / TUSCO-gene transcripts.
- With `--regions`, only TUSCO genes inside the regions count.

### `sections.priming` (A31, A35)
`{"n": 0, "risk": 0, "down_a_hist": {"edges": [0,0.1,…,1.0], "counts": []}}` from `down_a_frac`;
11 edges, `len(counts) == len(edges)`, the last bin is exactly 1.0.

### `sections.variation` (A32, A35)
For genes with ≥ 2 isoforms (the isoforms' own gene ids): `{"genes_multi", "ends_only",
"splice_only", "both", "neither"}`, with `genes_multi` = the sum of the other four.
- ends only: every isoform shares one intron chain (or all are mono-exon) and ends differ by > 50 bp;
- splice only: more than one chain, and isoforms sharing a chain have ends within 50 bp;
- neither: one chain, ends within 50 bp; both: the rest. Rules: `variation.py`.

## 11. The report viewer
Tabs: What differs, Structure, Recovery, Diversity, Agreement, Ends, Quantification (view id
`reads`), Compute
(`assets/report/report.js`; developer notes in `assets/report/README.md`). Out of scope for now:
a Loci tab (read-to-isoform locus plots) and a read 5′-end heatmap, which need read alignments
in the bundle; signal metaplots around ends. Recall against a sample's own truth sets is section 12.

### Cross-run data (A33, A35)
Computed by `build_report.py` per sample (`bin/flair_eval_lib/crossrun.py`), not stored per
bundle, and precomputed for every ordered run pair or baseline so switching the baseline needs no
rebuild:
- `by_sample.<s>.end_agreement`: Jaccard of TSS/TES end clusters (same strand, single-linkage
  within 50 bp) between run pairs, from `isoforms.bed`.
- `by_sample.<s>.tusco_incremental`: genes TP/PTP in a run but FN in the baseline, from
  `sections.tusco.per_gene`.
- `by_sample.<s>.lost_deciles`: the baseline's detected genes in ten CPM deciles (from
  `counts.tsv`), with each other run's missing genes per decile.

### Rules (A34)
- One sample at a time: a Sample picker drives every tab; nothing is pooled across samples.
- Colour = run: a fixed `--series-N` slot in run order, so colours never move when the baseline
  changes. Shared and "other" are neutral grey.
- The baseline is drawn first, with a larger outlined dot, is switchable in the control bar, and
  every comparison is "vs baseline".
- With ≥ 2 runs, plots overlay by default, with a side-by-side toggle.
- "All metrics" and "Biggest differences" come from a metric registry: each section registers
  its rows, with a lower-is-worse flag. "Biggest differences" takes at most one row per family.
- The Quantification tab says once that dev counts from the genome BAM (A21), as a callout only
  when the compared runs counted reads differently.
- A panel whose section is missing is replaced by a one-line reason in the reader's terms (the
  param and preset that fill it, never section keys). Old bundles without the newer sections
  still render.

## 12. End accuracy against truth sets
Each FLAIR run's isoform 5′ and 3′ ends are scored against the sample's own orthogonal truth
sets (CAGE, PRO-cap, RAMPAGE for the 5′ end; dRNA for the 3′ end), next to the end databases of
section 10. The ends half of section 10 (`sections.ends`, nearest-site distances from
`ends.py`) stays as it is. This section adds a second, independent score from the judge below.

### The truth sheet (A36; schema `assets/schema/truth.schema.json`)
`--truth` is a CSV with one row per truth set:

| column | req | notes |
|---|---|---|
| `sample` | yes | a sample of this launch. Rows for other samples are skipped with a warning, so one sheet can serve many launches |
| `end` | yes | `tss` (5′ end) \| `tes` (3′ end) |
| `assay` | yes | the set's name in the report: `^[A-Za-z0-9][A-Za-z0-9._-]*$`, unique per sample. `refTSS` and `PolyA_DB` are reserved for the databases |
| `path` | yes | BED6 or wider, stranded sites; columns 1–3 are the interval matched, column 6 the strand. Relative paths resolve as in A3 |
| `provisional` | no | `true` \| `false` (empty = false). True marks the set as provisional in the report |

- A sample may have any number of sets per end. Every problem in the sheet is listed before any
  task runs, together with the samples.csv and runs.yaml problems.
- Truth data stay outside the repository; the sheet only points at them.
- The lab's sheet for six cell lines (A549, K562, HCT116, GM12878, HEK293T, WTC11) is
  `/private/groups/brookslab/flair-eval/reference/ortho_truth/truth_6lines.csv` (with `SOURCE.tsv`).
  Its sample column holds the cell line names; every row is provisional.

### END_ACCURACY and the judge (A37)
- One task per job (sample × run) whose sample has at least one end set: its truth sets and/or
  `ref_tss` / `ref_pas` (section 10). Without any, the task does not run and the section is absent.
- **The judge.** `bin/end_accuracy.py` scores nothing itself. It runs flair4_ted's end judge,
  `ted_end_precision.py`, which is copied unchanged with its helpers into
  `modules/end_accuracy/resources/usr/bin/` (provenance, md5s and documented quirks:
  `modules/end_accuracy/README.md`). It then reads each call's one-row
  `precision_recall_summary.tsv`.
- **Rule.** An isoform end is a hit when it lies on the same strand within 50 bp of any part of a
  site's interval. Only chr1–22, X and Y are scored. Single-end matching is many-to-one.
  - precision = hits / isoform ends
  - recall = sites with a hit / all sites in scope
  - both ends = isoforms whose two ends both hit / all isoforms, counting isoforms of one
    junction chain on the same (5′ site, 3′ site) pair once
  - With `--regions`, only isoforms and sites inside the regions are scored.
- **Calls.** Each kind is planned separately (`truth`, `reference`): every 5′ set with every
  3′ set when the kind has both ends, otherwise each set alone. A set's single-end numbers are
  the same in every call it is in.
- The judge's own outputs and log are copied to `provenance/end_accuracy/<call>/`, where
  `<call>` is the set ids joined by `__`.

### `sections.end_accuracy` (A38)
```json
{"scorer": {"script": "ted_end_precision.py", "source": "…", "md5": "…", "window": 50, "regions": null},
 "provisional": true,
 "sets": [{"id": "CAGE", "kind": "truth", "end": "tss", "assay": "CAGE", "file": "…", "path": "…",
           "md5": "…", "provisional": true}],
 "ends": [{"set": "CAGE", "kind": "truth", "end": "tss", "assay": "CAGE", "isoforms": 0,
           "isoform_ends": 0, "sites": 0, "sites_hit": 0, "precision": 0.0, "naive_precision": 0.0,
           "dedup_tp": 0, "redundant_calls": 0, "recall": 0.0, "f1": 0.0}],
 "pairs": [{"kind": "truth", "tss_set": "CAGE", "tes_set": "dRNA", "isoforms": 0, "unique_pairs": 0,
            "both_hit": 0, "both_hit_raw": 0, "precision": 0.0, "naive_precision": 0.0}]}
```
- `sets`: references first (ids `refTSS`, `PolyA_DB`), then truth sets in sheet order. `path` is
  the original path and `md5` is the md5 of the staged file.
- `ends`: one row per set, ordered by kind (reference, truth), then end (tss, tes), then set order.
  Fractions are 0–1 and null when the denominator is 0.
- **Reference rows are precision only.** `recall`, `f1` and `sites_hit` are null, because the
  databases pool many tissues (A29). Truth rows have all three.
- `provisional` is true when any truth set is.
- `dedup_tp` and `redundant_calls` come from the judge's junction-chain grouping. They are kept
  as diagnostics. They are not the numerator of `precision`.
- `regions` lists the judge's `--region` values (`chr:start-end`), or null.

### `tables/end_accuracy.tsv` (A39)
One header line, then one row per `ends` entry, then one row per `pairs` entry:
```
kind end set assay partner isoforms isoform_ends sites sites_hit precision naive_precision
dedup_tp redundant_calls recall f1 unique_pairs both_hit both_hit_raw pair_precision pair_naive_precision
```
- Pair rows have `end` = `pair`, `set` = the 5′ set and `partner` = the 3′ set.
- Empty cells are null. Floats are written with 6 significant digits.

### Report (A40)
On the Ends tab, when any shown run has `sections.end_accuracy`:
- The **Precision at 50 bp** table takes its numbers from this section. It has one column per
  end set:
  - databases: precision
  - truth sets: precision and recall, with F1 and sites hit
- Each truth set gets a panel with precision, recall and F1.
- A **Both ends on a site** panel shows the pair precision for each 5′ × 3′ pair.
- A **Provisional truth** callout names the provisional sets.
- A note says that the distance charts (section 10) measure each end against the nearest site
  position, so their 50 bp values can differ a little from the table.
- Without this section the table falls back to `sections.ends`. A hidden **Ends against truth
  sets** card then names `--truth`.
- Metric registry rows are not added for truth sets, because their names differ from sample to
  sample.

## 13. Structure QC (SQANTI3-style checks)
`sections.structure_qc`, written by `structure_metrics.py` for every run
(`bin/flair_eval_lib/structure_qc.py`). Each key follows a figure or column of SQANTI3's QC
report. The rules are taken from SQANTI3's definitions and written anew: no SQANTI3 code is
copied (SQANTI3 is GPL-3).

### Subcategories (A42)
`subcategories` = `{"ISM": {...}, "NIC": {...}, "NNC": {...}}`, counts of SQANTI3's names:
- ISM: `3prime_fragment` (lacks the 5′ exons; the 5′-truncation signal), `5prime_fragment`,
  `internal_fragment`, `intron_retention` (an exon spans more than one reference exon),
  `mono-exon`. Multi-exon ISMs compare their first and last junction with the matched
  reference transcript, on the transcript strand.
- NIC: `combination_of_known_junctions`, `combination_of_known_splicesites`,
  `intron_retention` (an exon contains a reference intron), `mono-exon_by_intron_retention`.
- NNC: `at_least_one_novel_splicesite`, `intron_retention`.

### The other keys (A42)
- `redundancy`: FSM and FSM+ISM isoforms per matched reference transcript.
  `fsm_isoforms`, `fsm_refs`, `fsm_refs_multi` (references with more than one FSM), the same for
  FSM+ISM, and `fsm_hist` / `fsm_ism_hist` over edges `[1,2,3,4,5,6,8,11]`.
- `by_category`: `length` and `exons` histograms per category, with the edges of
  `histograms.transcript_length` and `histograms.exons_per_isoform`.
- `ref_recovery`: the same two histograms for reference transcripts in scope (`all`) and for
  those with an FSM (`fsm`).
- With the genome (the pipeline always passes it):
  - `canonical`: `multi_exon`, `all_canonical` (every junction GT-AG, GC-AG or AT-AC), and the
    same per category.
  - `rts`: distinct stranded junctions with an RT-switching repeat. The two windows are the 12 nt
    ending one base into the intron at the donor and the 12 nt ending one base into the
    downstream exon at the acceptor, read on the transcript strand. A junction is RT-switching
    when the windows share an 8-nt stretch with at most one mismatch (SQANTI3's
    `rt_switching.py` defaults). `by_type` splits known/novel × canonical/non-canonical;
    `isoforms_with_rts` counts multi-exon isoforms with at least one. `params` records the values.
  - `priming`: the internal-priming risk of section 10 (`down_a_frac` ≥ 0.6) by category
    (`by_category`) and by `mono`/`multi` exon (`by_exons`), each `{"n", "risk"}`. When
    productivity ran, `sections.add_sections` adds `by_coding` (`coding` = PRO/PTC/NST,
    `noncoding` = NGO).

### Support in the sample's junctions file (A43)
The sample's `junctions` column (section 1), staged as `refs/sj/*` and passed as
`--junction-file`. It is the same file for every run of the sample, so runs compare against the
same evidence. A run's own runs.yaml junctions are not used here.
- STAR `SJ.out.tab` (9 columns; 1-based first and last intron base; unique reads, column 7) or
  BED (intron `[start, end)`, reads in the score column). Strand is ignored.
- `junction_support` = `{"file", "md5", "file_junctions" (distinct junctions in the file),
  "junctions", "supported", "known": {"n", "supported"}, "novel": {"n", "supported"},
  "multi_exon", "all_supported", "by_category"}`. Supported = at least one read.
- It is independent evidence only when the file comes from short reads (STAR on Illumina).
  A file made from the same long reads (intronProspector on the sample's BAM), or one FLAIR was
  given to correct junctions, supports the runs partly by construction. The report says so and
  names the file.

### Report (A44)
- Structure tab:
  - an ISM/NIC/NNC subcategory chart;
  - a junction QC panel: all-canonical, RT-switching and junction-file support.
- Recovery tab: reference recovery by transcript length.
- Diversity tab: FSM redundancy.
- Ends tab: priming risk by group.
- Metric-registry rows, each needing its key: ISM 3′ fragments (share of ISM), isoforms with all
  junctions canonical, RT-switching junctions, junctions in the junctions file, multi-exon
  isoforms fully supported, references with more than one FSM.
