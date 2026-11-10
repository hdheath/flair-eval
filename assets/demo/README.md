# Demo bundles

The Pages demo report (`/demo/`) is built from these three bundles: FLAIR 2.2.0, FLAIR 3.0.1 and
the `dev` tip (389748e, 2026-11-09). All are on LRGASP WTC11 PacBio cDNA (ENCODE),
hg38 chr17:7,400,000-9,400,000 (TP53 locus), GENCODE v49 basic, with default FLAIR settings.

They are the pipeline's own outputs, with input paths cut to file names and `provenance/` left out.
The report's tests use `tests/fixtures/bundles/` instead.
