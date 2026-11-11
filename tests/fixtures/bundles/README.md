# Report fixtures: real FLAIR bundles

Ten bundles from existing FLAIR runs, restricted to **chr17:7,400,000-9,400,000** (hg38; TP53
locus, 68 protein-coding genes). Input paths in `bundle.json` are reduced to file names.

Each bundle folder holds what the report reads:

| file | made by |
|---|---|
| `bundle.json`, `tables/isoforms.tsv`, `tables/categories.tsv`, `tables/genes.tsv` | `bin/structure_metrics.py` |
| `isoforms.bed` | the run's FLAIR isoform BED (first 12 columns), region isoforms only |
| `counts.tsv` | FLAIR quantify on the region's reads (see below) |
| `productivity.tsv` | `bin/flair_adapter.py productivity` (FLAIR predictProductivity) |
| `sections.quant`, `sections.productivity` | `python3 -m flair_eval_lib.sections <bundle dir>` |
| `compute` | `/usr/bin/time -v` of the quantify + productivity steps on the slice |
| `run_index` | the run order above (A17) |
| `sections.quant.input` | `reads` for 2.2.0 / 3.0.1, `bam` for the dev runs (A21) |
| `tables/ends.tsv`, `sections.ends`, `.priming`, `.variation` | round 3, see below |
| `sections.tusco` (`wtc11_pacbio` only) | round 3, **synthetic**, see below |

## Samples and runs

| sample | data |
|---|---|
| `wtc11_pacbio` | LRGASP WTC11 PacBio cDNA (ENCODE), one minimap2 alignment shared by all runs |
| `wtc11_ont` | LRGASP WTC11 ONT cDNA replicate 3 (ENCODE), shared alignment |

| run_id | FLAIR | settings |
|---|---|---|
| `flair-2.2.0` | release 2.2.0 | `correct` + `collapse` |
| `flair-3.0.1` | release 3.0.1 | `--junction_support 3 --stringent --check_splice --max_ends 1` |
| `dev-57040d2` | dev branch @ 57040d2 (2026-08-19) | `--junction_support 2 --max_ends 1 --quality 0` |
| `dev-425b260` | dev branch @ 425b260 (2026-09-03) | `--junction_support 2 --max_ends 1 --quality 0` |
| `dev-425b260-normends` | dev branch @ 425b260 | `--junction_support 2 --max_ends 1 --normalize_ends` |

Reference: GENCODE v49 basic (chr17 slice). Genome: GRCh38 primary assembly (junction motifs).
The assemblies were genome-wide runs; quantify and productivity were run on the slice only
(reads from the shared BAMs over chr17:7.4-9.4 Mb: 27,799 PacBio, 110,973 ONT).

Caveats, so the numbers are not over-read:
- **Quantify differs by version** (A21, recorded in `sections.quant.input`; the report notes it).
  FLAIR 2.2.0 and 3.0.1 quantify from reads (via the adapter). The dev branch quantifies from
  the genome BAM. The adapter's dev path passed `-g`, which this
  dev `flair quantify` rejects (it only accepts `--genome`), so the dev counts were made with the
  same command and `--genome` spelled out. FLAIR itself was not changed.
- CPM is relative to the slice's reads, so values are much higher than genome-wide CPM.
- `compute` covers quantify + productivity on the slice, not the genome-wide assembly.

## Round-3 sections (CONTRACTS A28-A32)

Added by hand for the report tests, before `structure_metrics.py` wrote them, so they follow the
contract shapes but were not made by the pipeline. Regenerate them with
`structure_metrics.py --ref-tss --ref-pas --tusco` once that lands.

- `tables/ends.tsv`, `sections.ends`, `sections.priming` and `sections.variation` were computed
  from each bundle's `isoforms.bed` and `tables/isoforms.tsv` by a one-off script that reads
  A28-A32 independently of `bin/flair_eval_lib/ends.py`: refTSS v4.1 and PolyA_DB v4.1 (hg38),
  GENCODE v49 basic for the FSM/ISM reference ends, and the GRCh38 genome for the polyA motif and
  the downstream-A fraction. `refs` records the files' names, md5 and genome-wide site counts.
  `priming.down_a_hist` has 11 bins (the last is exactly 1.0, as in section 3).
- `sections.tusco` is **made up**: the fixture region holds no TUSCO gene. `sensitivity` = TP / genes as in `tusco.py`; `precision` is a made-up number. `wtc11_pacbio` runs
  carry random calls over the 46 real TUSCO gene ids (seeded by run id, `genes_in_regions` = 46)
  so the TUSCO panels have something to draw. `wtc11_ont` has no `tusco` section, which tests the
  hidden-panel path.

## Cross-check of the classifier against SQANTI3

On the same isoforms (chr17 slice, 8 of the 10 runs, 1,518 isoforms), category agreement is
99.0% and multi-exon FSM subcategory agreement is 100%. The 15 differences are 10 sub-chain
isoforms with a retained intron (ISM here, NIC in SQANTI3) and 5 mono-exon isoforms that SQANTI3
calls antisense although they overlap a same-strand exon. Mono-exon FSMs get CIM/FSM5/FSM3/FSMTE
here; SQANTI3 labels them only "mono-exon".

## Rebuilding

Per bundle:
```
structure_metrics.py --isoforms-bed <run>.isoforms.bed --isoforms-gtf <run>.isoforms.gtf \
  --gtf gencode.v49.basic.annotation.gtf --sample <sample> --run-json run.json --env-json env.json \
  --regions regions.bed --genome GRCh38.primary_assembly.genome.fa --out <sample>/<run_id>
flair_adapter.py quantify --isoforms-fa <slice>.fa --isoforms-bed <slice>.bed --reads <slice>.fq.gz \
  --bam <slice>.bam --genome GRCh38.primary_assembly.genome.fa --sample <sample> --out quantify
flair_adapter.py productivity --isoforms-bed <slice>.bed --genome GRCh38.primary_assembly.genome.fa \
  --gtf gencode.v49.basic.chr17.gtf --out productivity
python3 -m flair_eval_lib.sections <sample>/<run_id>
```
with `regions.bed` = `chr17  7400000  9400000`.
