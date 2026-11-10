# FLAIR versions: CLI and output differences

flair-eval never patches FLAIR. Every version difference is handled in `bin/flair_adapter.py`
(the CLI) or `bin/flair_env.py` (the install). This page lists what differs, what the adapter does
about it, and the FLAIR bugs we document instead of fixing.

Versions checked on 2026-11-09:

| id here | what | `flair --version` |
|---|---|---|
| **2.2.0** | PyPI `flair-brookslab==2.2.0` | `FLAIR 2.2.0` |
| **3.0.1** | PyPI `flair-brookslab==3.0.1` | `FLAIR 3.0.1` |
| **dev** | BrooksLabUCSC `dev` at d033600 (2026-10-26) and at 389748e (2026-11-09) | `FLAIR 3.0.0+dev` |
| **TED** | TED fork of dev (base 51fa32f, 2026-05-22), HEAD 20c0f47 + uncommitted edits | `FLAIR 3.0.0+dev` |

dev and TED report the same version string, so the adapter **detects the CLI by probing**
`flair --help` and each subcommand's `--help` rather than by version number.
`flair_adapter.py detect --json` prints what it found.

## 1. Subcommands and the CLI the adapter picks

| | 2.2.0 | 3.0.1 | TED | dev |
|---|---|---|---|---|
| `flair --help` format | `module: align, correct, collapse, quantify` | argparse `{align,correct,transcriptome,...}` | argparse `{align,transcriptome,...}` | `subcommands:` block |
| `align` | yes | yes | yes | **no** → `align` writes `skipped.json` |
| `correct` / `collapse` | yes | yes (not used) | no | no |
| `transcriptome` | no | yes | yes | yes |
| **assemble runs** | BAM → `bam2Bed12` → `correct` → `collapse` | `transcriptome` | `transcriptome` | `transcriptome` |
| `predictProductivity` script | yes | yes | yes, other CLI | **no** (CDS is inline) |
| `quantify` input | reads FASTA/Q manifest + `-i isoforms.fa` | same | same + **required** `--isoform_bed` | **genome BAM** manifest + `-g genome` + `--isoform_bed` |

## 2. `assemble` (transcriptome / collapse)

| | 2.2.0 collapse | 3.0.1 transcriptome | TED transcriptome | dev transcriptome |
|---|---|---|---|---|
| alignment input | reads BED12 from `-q` (adapter makes it from the shared BAM with FLAIR's own `flair.bam2Bed12`, which keeps primary flag 0/16 like 2.x align) + raw reads `-r` | `-b BAM` | `-b BAM` | `-b BAM` |
| `--sample_name` | — | — (silently ignored if given) | **required** (adapter passes `--sample flaireval`) | **required** (same) |
| parallel flag | — | `--parallelmode` | `--parallel_mode` | `--parallel_mode` |
| adapter adds | `--generate_map --temp_dir` | `--parallelmode byregion` | `--parallel_mode byregion --generate_map` | `--parallel_mode byregion --generate_map --temp_dir` |
| unknown flags | error | **ignored silently** (`parse_known_args`, logged only at `--log-level INFO`) | error | error |
| temp files | `--temp_dir` | uuid dir in cwd | cwd | `--temp_dir` (default `$TMPDIR`) |
| isoform BED | 12 col | 12 col | **30 col** (31 with `--ted`) | **25 col** (`flair_bed.as`) |
| isoform names | `<transcript>_<gene>` (`<read>_chr:pos` when novel) | `<transcript>_<gene>` | `FLT00000001…` | `FLT00000001…` |
| gene ids | in the name / GTF `gene_id` (ENSG… or `chr:pos`) | same | BED col 13 `FLG…` | BED col 13 `FLG…` |
| GTF `transcript_id` | name **without** `_<gene>` | same | = BED name | = BED name |
| read map | `<o>.isoform.read.map.txt` (only with `--generate_map`) | `<o>.isoform.read.map.txt` (always) | `.txt` with `--generate_map` | `<o>.isoform.read.map.tsv`, header `name reads` |
| CDS / productivity | no | only with `--predictCDS` (rewrites the name field) | inline: thick = CDS, col 18 | inline: thick = CDS, col 18, `<o>.aaseq.tsv` |
| also writes | — | `.isoform.counts.txt` | `.isoform.counts.txt`, `.bb` | `.isoform.counts.tsv` (header), `.bb` |

Default settings differ a lot between versions. The pipeline runs each version at its own defaults
unless a run sets `transcriptome_args`, so these differences are part of what a comparison measures:

| default | 2.2.0 collapse | 3.0.1 | TED | dev (d033600) |
|---|---|---|---|---|
| min reads per isoform | `-s 3` | `-s 3` | `--sjc_support` / `--se_support` | `--sjc_support 1`, `--single_exon_support 3` |
| relative support | none (`-s <1` = fraction) | none | `--frac_support 0.05` | `--frac_support 0.05` |
| stringent / check_splice | off | off | on (`--no_stringent`, `--no_check_splice`) | on (same) |
| `--max_ends` | 2 | 2 | 1 | 1 |
| `--filter` | default | default | (TED-specific) | nosubset |
| MAPQ | `--quality 0` | fixed 0 | `--quality` | `--quality 1` |
| `--ss_window` | 15 (correct) | 15 | 15 | 10 |

### What the adapter normalises (CONTRACTS §2)
- `isoforms.bed`: the first 12 columns of FLAIR's BED, one row per name (a repeated name keeps
  its first row and adds a warning), sorted by chrom/start/end/name, score clamped to 0–1000, no
  trailing commas in the block lists. thickStart/thickEnd are kept (dev and TED put the CDS there).
- `isoforms.gtf`: rebuilt from `isoforms.bed` so it always matches it: `transcript` + `exon` lines,
  `transcript_id` = the BED name (FLAIR 2.2.0/3.0.1 drop the `_<gene>` suffix in their own GTF;
  the normalised one does not). `gene_id` comes from, in order: BED column 13 (dev/TED), FLAIR's
  GTF (matched as `<transcript_id>_<gene_id>` = BED name), the part of the name after its last `_`,
  else `novel_gene_<n>`. `adapter_log.json` counts each source (`outputs.gene_id_sources`).
- `isoforms.fa`: FLAIR's FASTA, headers cut at the first space, one record per BED isoform, in BED order.
- `read_map.tsv`: long format with header `read_id	isoform_id`, from whichever read map the version wrote.
- `isoforms.flair.bed`: FLAIR's own BED, verbatim (all columns). dev quantify and the inline
  productivity read it. `native/` keeps FLAIR's GTF and counts file for provenance.
- `--regions BED`: the BAM is cut with `samtools view -L` (reads overlapping the regions) before
  FLAIR sees it; for 2.2.0 the raw reads given to `collapse` are cut to the reads in that BAM.

## 3. `quantify`

| | 2.2.0 | 3.0.1 | TED | dev |
|---|---|---|---|---|
| manifest line | `sid cond batch reads.fq` | same | same | `sid cond batch genome.bam` (`--manifest`) |
| isoforms | `-i isoforms.fa` | `-i isoforms.fa` | `-i isoforms.fa` + `--isoform_bed` (required) | `--isoform_bed` FLAIR 25-col BED (BED12 rejected) + `-g genome` |
| `_` in sample/condition/batch | not checked | **rejected** | not checked | — |
| stringent / check_splice | flags, off | flags, off | **hardwired on** | — (always strict) |
| zero-read isoforms | omitted | **omitted** | **omitted** | omitted |
| counts file | `<o>.counts.tsv`: `ids  s_c_b` (floats) | same (ints) | same | `<o>.counts.tsv` (d033600) / `<o>.isoforms.counts.tsv` (389748e): `gene_id isoform_id sid` |
| `--tpm` | CPM | **CPM, not TPM** (no lengths passed) | CPM | — |

The adapter writes a fixed manifest `s1 c1 b1 <reads|bam>` (the real sample name goes into
`adapter_log.json`), adds every isoform missing from FLAIR's table with count 0, and computes
`cpm = count / Σcount × 10^6` itself. `counts.tsv` has the header `isoform_id count cpm`, one row
per `isoforms.fa` record in FASTA order. FLAIR's own table is kept as `flair.counts.tsv`.

**dev needs more than the contract passes.** The contract `quantify` call gives the isoform FASTA
and reads; dev quantify counts from the genome BAM. The adapter accepts optional
`--bam B --genome G --isoforms-bed I`; it finds FLAIR's own BED as `isoforms.flair.bed` next to
`--isoforms-bed` or `--isoforms-fa` (also behind Nextflow's staging symlinks). Without `--bam` and
`--genome`, dev quantify writes `skipped.json`. TED quantify needs FLAIR's BED too and finds it the
same way.

## 4. `productivity`

| | 2.2.0 | 3.0.1 | TED | dev |
|---|---|---|---|---|
| source used | `predictProductivity` | `predictProductivity` | inline (BED col 18) | inline (BED col 18) |
| standalone CLI | `-i -g -f -o [--append_column] (--firstTIS \| --longestORF)` | same | `-i -g -f -o [--is_flair_bed]` (no `--append_column`) | none |

The adapter prefers the version's inline call (column 18 of `isoforms.flair.bed`, PRO/PTC/NST/NGO,
CDS in thickStart/thickEnd). That is what `flair transcriptome` itself reported. Otherwise it runs
`predictProductivity --append_column --longestORF`, the same ORF mode FLAIR's own pipeline uses
(`transcriptome --predictCDS`, 2.x `collapse`). The input BED is copied to a plain name first,
because TED's script splits paths on the first `.bed`. Output `productivity.tsv`:
`isoform_id cds_start cds_end class`, 0-based genomic thick coordinates, `NA` when FLAIR gave no
call. `adapter_log.json` records `outputs.method` (`inline` or `predictProductivity`).

## 5. `align`

| | 2.2.0 | 3.0.1 | TED | dev |
|---|---|---|---|---|
| available | yes | yes | yes | **no** |
| outputs | `<o>.bam`, `.bam.bai`, `.bed` | same | same | — |
| note | | `-secondary=no` is passed in a form minimap2 does not take, so minimap2 runs with `-s 0` (flair_align.py:156) | | the shared align step must use a version with `align` (`align_flair`, default 3.0.1) |

The adapter runs FLAIR's defaults (`flair align -g G -r R -t N`); `--gtf` is accepted for
interface symmetry and not passed on.

## 6. Known FLAIR bugs and traps (documented, not patched)
- **dev runtime tools outside pip**: dev `transcriptome` calls `bedToBigBed` (conda
  `ucsc-bedtobigbed`) and fails without it; dev's `misc/flair_conda_env.yaml` also lists
  `intron-prospector` and a `pysam` range. pip cannot install these, so `flair_env.py` adds the
  runtime conda packages of the commit's own `misc/flair_conda_env.yaml` (found on 389748e).
- **dev output names move**: between d033600 and 389748e, quantify's `<o>.counts.tsv` became
  `<o>.isoforms.counts.tsv` (+ `<o>.isoforms.sample_info.tsv`). The adapter looks for the known
  names, then for any `<o>.*counts*` table, and logs a warning when it had to guess.
- **3.0.1 transcriptome** ignores unknown flags (e.g. dev-only `--sample_name`, `--generate_map`).
  The adapter compares `--args` with each version's `--help` and records unknown flags as warnings.
- **3.0.1 combine** (not used in v1): with filtering on, `flair_combine.py:318` reads `theseisos`
  before assignment and crashes when the first locus has all chains rejected.
- **3.0.1 predictProductivity**: GENCODE `NMD_exception` never matches combine-style names; an
  open ORF that is a multiple of 3 drops its final codon (NST).
- **`--tpm`** in 3.0.1 and dev-pr684 writes CPM (counts_to_tpm.py is called without lengths).
- **TED**: `FlairBed.toRow` crashes when `frac_support` is None; standalone `predictProductivity`
  takes start codons from any gene (differs from the inline column 18) and splits paths on `.bed`;
  a genome-wide run in one job can deadlock even with the OMP exports. Run it per chromosome or
  with `--regions`, and always with `byregion`.
- **TED env provenance**: it is an editable install with uncommitted edits, so `resolved_commit`
  (20c0f47) is its HEAD, not the exact code; `flair_env.py` prints a note when this happens.
- Every adapter call exports `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=XGBOOST_NTHREAD=1`
  (TED byregion deadlocks without them).

## 7. Installing a version (`flair_env.py`)
- **release** → a conda env (`conda-forge`, `bioconda`) with Python chosen from the release's
  `requires-python` (3.12 when allowed), the uniform tool pins `minimap2=2.30 samtools=1.22.1
  bedtools=2.31.1`, and the runtime conda extras of the tag's `misc/flair_conda_env.yaml` (none
  for 2.2.0 and 3.0.1), then `pip install flair-brookslab==<v>`. PyPI is the route because it is
  FLAIR's primary release channel and it lets every version share the same aligner and tools, so a
  comparison isolates FLAIR. Bioconda's `flair=<v>` (it has 1.4–3.0.1, including 2.1.1, which
  PyPI lacks) is used only when PyPI has no such version; it brings its own tool versions unless
  pinned. `--route pip|bioconda` forces one.
- **git** → `git ls-remote` resolves the ref (annotated tags dereferenced; a short SHA needs a
  clone), the commit is fetched shallow by SHA into `<cache>/src/<sha>`, Python is chosen from its
  `pyproject.toml`, the runtime conda extras come from its `misc/flair_conda_env.yaml` (the R
  stack, longshot/whatshap and doc/build tools are left out), then `pip install <src>` (not
  editable) into `<cache>/envs/git-<sha12>`.
- **env** / **container** → used as is after `flair --version` works (`apptainer exec` for
  containers). For an env, `resolved_commit` comes from pip's `direct_url.json` (VCS or local
  checkout) when there is one.
- Builds happen under `<cache>/locks/<key>.lock` (`fcntl.lockf`), so concurrent tasks build a key
  once; an env without its `.flair-eval-env.json` marker is a broken build and is rebuilt. Build
  logs are in `<cache>/logs/`. Conda's package cache goes to `<cache>/conda-pkgs` unless
  `CONDA_PKGS_DIRS` is set (keeps it off `$HOME`). The key does not include the recipe: to
  rebuild a cached env (e.g. after changing tool pins), delete `<cache>/envs/<key>`.
- Overrides: `FLAIR_EVAL_CONDA`, `FLAIR_EVAL_TOOL_PINS`, `FLAIR_EVAL_PYPI`,
  `FLAIR_EVAL_CONTAINER_RUNTIME` (or the matching `--conda`, `--tool-pins`, `--route`,
  `--container-runtime` flags).

## 8. Tests
- `tests/python/test_flair_adapter.py` and `test_flair_env.py`: fake `flair`, `predictProductivity`,
  `samtools`, conda, pip, git and apptainer for each CLI generation; no FLAIR needed.
- `tests/python/test_flair_adapter_integration.py`: opt-in (`FLAIR_EVAL_IT_ENVS`, see its
  docstring) end-to-end run of every step with real envs, checking contract-shaped outputs.

Integration results, 2026-11-09, on FLAIR's own test inputs (`test/input/basic.reads.fa`, 1,013
reads; `basic.annotation.gtf`; hg38 chr12+chr17+chr20), shared BAM from 3.0.1 `align`, 4 threads
on a SLURM short node:

| env | `align` | isoforms (all / `--regions`) | genes | quantify (contract call) | quantify (`--bam`) assigned reads | productivity |
|---|---|---|---|---|---|---|
| 2.2.0 (existing env) | ok | 14 / 13 | 3 | ok | 763 | predictProductivity: 9 PRO, 5 PTC |
| 3.0.1 (existing env, Python 3.13) | ok | 15 / 11 | 3 | ok | 974 | predictProductivity: 7 PRO, 5 PTC, 3 NGO |
| 3.0.1 (`release: 3.0.1`, built from scratch, Python 3.12) | ok | 15 / 11 | 3 | ok | 974 | predictProductivity |
| TED 20c0f47+ (existing env) | ok | 7 / 4 | 3 | ok | 774 | inline: 5 PRO, 2 NGO |
| dev d033600 (existing env) | skipped | 7 / 4 | 3 | skipped (needs `--bam`) | 751 | inline: 5 PRO, 2 NGO |
| dev 389748e (`git: …flair, ref: dev`, built from scratch) | skipped | 7 / 4 | 3 | skipped (needs `--bam`) | 751 | inline |

The freshly built `release: 3.0.1` env gives byte-identical `isoforms.bed/.gtf`, `counts.tsv` and
`productivity.tsv` to the existing 3.0.1 env, and dev 389748e matches dev d033600 on this data.
Both builds took about 2 min each on a short node. A second resolve of the same release, started
5 s later, waited 121 s on the lock and then reused the env.
