# flair-eval

Compare FLAIR runs (different versions, branches or settings) on the same long-read data.
Results come out as one self-contained HTML report: open any single run, or tick several
runs to compare them.

> **Status:** v2 is under construction on the `v2` branch. See [docs/DESIGN.md](docs/DESIGN.md)
> and [docs/CONTRACTS.md](docs/CONTRACTS.md). The previous multi-assembler pipeline is in `legacy/`
> for reference, and on the `wip/2026-11-09-snapshot` branch.

Planned quick start:
```bash
nextflow run hdheath/flair-eval \
  --reads sample.fq.gz --genome genome.fa --gtf annotation.gtf \
  --flair 3.0.1,dev --outdir results
# -> results/report/report.html
```
