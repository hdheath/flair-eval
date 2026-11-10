# Report fixtures: real FLAIR bundles

Ten bundles made by `bin/structure_metrics.py` from existing FLAIR runs, restricted to
**chr17:7,400,000-9,400,000** (hg38; TP53 locus, 68 protein-coding genes). Only `bundle.json`
and `tables/` are kept. Input paths in `bundle.json` are reduced to file names.

| sample | data |
|---|---|
| `wtc11_pacbio` | LRGASP WTC11 PacBio cDNA (ENCODE), minimap2 alignment shared by all runs |
| `wtc11_ont` | LRGASP WTC11 ONT cDNA replicate 3 (ENCODE), shared alignment |

| run_id | FLAIR | settings |
|---|---|---|
| `flair-2.2.0` | release 2.2.0 | `correct` + `collapse` |
| `flair-3.0.1` | release 3.0.1 | `--junction_support 3 --stringent --check_splice --max_ends 1` |
| `dev-57040d2` | dev branch @ 57040d2 (2026-08-19) | `--junction_support 2 --max_ends 1 --quality 0` |
| `dev-425b260` | dev branch @ 425b260 (2026-09-03) | `--junction_support 2 --max_ends 1 --quality 0` |
| `dev-425b260-normends` | dev branch @ 425b260 | `--junction_support 2 --max_ends 1 --normalize_ends` |

Reference: GENCODE v49 basic (chr17 slice); genome GRCh38 primary assembly (for junction motifs).
All runs used the same short-read junction file. Run labels/tags are in each bundle.

Cross-check against SQANTI3 run on the same isoforms (chr17 slice, 8 of the 10 runs, 1,518
isoforms): category agreement 99.0%, multi-exon FSM subcategory agreement 100%. The 15
differences are 10 sub-chain isoforms with a retained intron (ISM here, NIC in SQANTI3) and
5 mono-exon isoforms that SQANTI3 calls antisense while they overlap a same-strand exon.
Mono-exon FSMs get CIM/FSM5/FSM3/FSMTE here; SQANTI3 labels them only "mono-exon".

To rebuild, run per bundle:
```
structure_metrics.py --isoforms-bed <run>.isoforms.bed --isoforms-gtf <run>.isoforms.gtf \
  --gtf gencode.v49.basic.annotation.gtf --sample <sample> --run-json run.json --env-json env.json \
  --regions regions.bed --genome GRCh38.primary_assembly.genome.fa --out <sample>/<run_id>
```
with `regions.bed` = `chr17  7400000  9400000`.
