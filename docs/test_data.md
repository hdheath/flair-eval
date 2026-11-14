# Test data (`-profile test`)

One small dataset is committed under `tests/data/`. It lets anyone run the pipeline in minutes:

```bash
flair-eval run --test --outdir /private/groups/brookslab/$USER/flair-eval/dd_mm_yy_test
```

or `nextflow run hdheath/flair-eval -r main -profile test,phoenix --outdir …` (off phoenix:
`-profile test,local,conda`). It compares FLAIR 3.0.1 (the baseline), 2.2.0 and dev on two samples.

| set | genome | annotation | reads |
|---|---|---|---|
| `human_region` | hg38 `chr19:1-2,000,000` | GENCODE v49 basic, genes wholly inside the window | 20,000 PacBio + 20,000 ONT (WTC11 cDNA) |

The reads come from the LRGASP WTC11 libraries in ENCODE. Spike-in (SIRV) test data belongs to
a separate spike-in pipeline branch and is not part of this test set.

## Files

All paths are relative to `tests/data/`. In `test_samples.csv`, paths are relative to the CSV's
own folder (CONTRACTS A3).

| file | bytes | content |
|---|---:|---|
| `test_samples.csv` | | 2 samples: `human_region_pacbio` and `human_region_ont` |
| `test_runs.yaml` | | `flair-3.0.1` (release; first, so it is the baseline), `flair-2.2.0` (release) and `flair-dev` (git `dev`) |
| `test_regions.bed` | | `chr19 0 2000000`, the whole contig. Optional: the genome already is the region. |
| `md5sums.txt` | | md5 of every file in `human_region/`; check with `cd tests/data && md5sum -c md5sums.txt` |
| `build_test_data.sh` | | rebuilds `human_region/` and `md5sums.txt` from the sources (deterministic: same md5s every run) |
| `human_region/human_region.genome.fa` (+ `.fai`) | 2,033,341 | hg38 chr19 bases 1–2,000,000 as one contig named `chr19` (length 2,000,000). Coordinates are identical to hg38. |
| `human_region/human_region.gtf` | 6,948,133 | GENCODE v49 basic: all lines of the 165 genes that lie wholly inside the window (758 transcripts, 71 protein-coding genes, e.g. PTBP1, TCF3, BSG, CIRBP, DAZAP1, MBD3, GPX4, STK11, ABCA7) |
| `human_region/human_region.regions.bed` | 16 | `chr19 0 2000000` |
| `human_region/human_region.pacbio.fa.gz` | 1,752,650 | 20,000 PacBio reads (of 52,626 in the window) |
| `human_region/human_region.ont.fa.gz` | 2,879,511 | 20,000 ONT reads (of 34,492 in the window) |
| `human_region/human_region.refTSS.bed` | 24,278 | refTSS v4.1 human hg38 TSS clusters wholly inside the window (620 sites). refTSS is CC BY 4.0: Abugessaisa et al., *J Mol Biol* 2019 and later releases, https://reftss.riken.jp |
| `human_region/human_region.polyA_DB.bed` | 31,190 | PolyA_DB v4.1 main polyA sites inside the window (648 sites), one row per cleavage site `[pos-1, pos)`. Source: PolyA_DB v4.1 `HumanPas` release (`hg38.PAS.main.tsv`), redistributed as a small subset with citation |
| `human_region/human_region.tusco.tsv` | 3,158 | SQANTI3 `tusco_human.tsv` (TUSCO single-isoform genes, 46 genes), copied whole from SQANTI3 commit 589dfd4 (GPL-3.0). No TUSCO gene lies in the window, so the test profile's TUSCO panel is empty by design. |
| `build_test_refs.sh` | | cuts the three reference files above from the lab copies (the `hg38` preset files in `/private/groups/brookslab/flair-eval/reference/`); `build_test_data.sh` calls it |

Total: 13.6 MB on disk. Git stores it compressed, mostly because the GTF and genome compress well.

**Reads are gzipped FASTA, with no quality strings.** That keeps the set small enough to commit,
with no hosting and no network access needed in CI. ONT FASTQ would be about 21 MB. Neither
minimap2 nor FLAIR needs the qualities. PacBio read names keep the source ENCODE file accession
as a prefix (`ENCFF245IPAm54284U_…/ccs`). ONT read names are the original UUIDs.

## How the reads were chosen

- Reads were taken from existing alignments of the full ENCODE files: ONT with `minimap2 -ax splice`
  (the `flair align` settings), PacBio with `minimap2 -ax splice`, both with minimap2 2.30 against
  GRCh38 primary.
- Only primary alignments were used (`-F 0x904`): reads whose primary alignment lies wholly inside
  chr19:1–2,000,000. Because reads were chosen by where their *primary* alignment falls, reads
  from paralogs elsewhere in the genome are left out, and every read realigns to the sliced genome.
- The 20,000 reads per set are the ones with the smallest `md5("flair-eval-testdata-v1:" + read_name)`.
  This is uniform and reproducible, and keeps natural expression levels.
- The reads are written back as the full original sequences (`samtools fasta` on the primary
  record, so reverse-strand reads are reverse-complemented back). They are in the BAM's
  coordinate order.

## Sources and licences

| piece | source | terms |
|---|---|---|
| PacBio reads | ENCODE experiment [ENCSR507JOF](https://www.encodeproject.org/experiments/ENCSR507JOF/) (LRGASP WTC11, PacBio Sequel II, Mortazavi lab, UCI). FASTQ files ENCFF563QZR (rep 1), ENCFF370NFS (rep 2) and ENCFF245IPA (rep 3); all three reps are pooled. | ENCODE data use policy |
| ONT reads | ENCODE experiment [ENCSR539ZXJ](https://www.encodeproject.org/experiments/ENCSR539ZXJ/) (LRGASP Nanopore PCR-cDNA WTC11, Brooks lab, UCSC), file ENCFF263YFG (rep 1) | ENCODE data use policy |
| hg38 | GRCh38 primary assembly (GENCODE `GRCh38.primary_assembly.genome.fa`) | Genome Reference Consortium reference, publicly available |
| GENCODE | GENCODE v49 basic annotation (`gencode.v49.basic.annotation.gtf`) | gencodegenes.org: "All the GENCODE project data is open access" (EMBL-EBI terms of use) |

**ENCODE terms** (https://www.encodeproject.org/help/citing-encode/, read 2026-11-09): "External data
users may freely download, analyze and publish results based on any ENCODE data without
restrictions." ENCODE asks users to cite the ENCODE Consortium, acknowledge the production labs,
and give the ENCSR/ENCFF accessions used. Those accessions are listed above. The policy does not
use the word "redistribute".

Please cite: the ENCODE Consortium; LRGASP (Pardo-Palacios et al., *Nature Methods* 2024); GENCODE
(Frankish et al.).

## Runtime check (FLAIR 3.0.1)

Each sample was run with `flair align -t 4` then
`flair transcriptome -f <gtf> -t 4 --parallelmode byregion`. Env `flair-3.0.1` (FLAIR 3.0.1, minimap2 2.30).
SLURM `short` partition, 4 CPUs, 8 GB. Run on 2026-11-10.

| sample | align | transcriptome | total | reads aligned | isoforms | peak RSS (align) |
|---|---:|---:|---:|---:|---:|---:|
| `human_region_pacbio` | 105 s | 38 s | 2.4 min | 19,998 | 221 | 2.8 GB |
| `human_region_ont` | 66 s | 22 s | 1.5 min | 20,000 | 120 | 1.9 GB |

Both samples finish in under 5 minutes. 4 CPUs and 4 GB are enough for each FLAIR step on this data.

## Rebuilding

```bash
conda activate flair-3.0.1          # any env with samtools >= 1.12, python3, gzip
bash tests/data/build_test_data.sh  # writes tests/data/human_region/ and md5sums.txt
```

The source paths at the top of the script point at the lab's local copies on PRISM. Each one can
be overridden with an environment variable of the same name, e.g. after downloading the ENCODE
files and aligning them yourself. `N_READS` (default 20000) and `SEED` change the subsample.
Changing either changes the md5s.
