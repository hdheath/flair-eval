# flair-eval (CLI)

Set up, check and launch [flair-eval](https://github.com/hdheath/flair-eval), the Nextflow
pipeline that compares FLAIR runs (versions, branches, settings) on the same data.

```bash
pip install "flair-eval @ git+https://github.com/hdheath/flair-eval@v2#subdirectory=cli"
OUT=/private/groups/brookslab/$USER/flair-eval/dd_mm_yy_name   # --outdir is required
flair-eval doctor            # check Java, Nextflow, conda, disk; detects phoenix
flair-eval run --test --outdir $OUT-test   # try the bundled test data
flair-eval init              # questions -> samples.csv + runs.yaml
flair-eval check             # validate them against the pipeline's schemas
flair-eval run --outdir $OUT
flair-eval run --reads r.fq.gz --genome hg38 --flair 2.2.0,3.0.1,dev --outdir $OUT
flair-eval report out1/bundles out2/bundles --outdir $OUT-report
```

On the UCSC PRISM phoenix cluster no `-profile` is needed, `--genome hg38` is the lab reference,
and `run`/`report` submit the Nextflow driver as a SLURM job (`--foreground` to run it here).

Full guide: [docs/usage_cli.md](../docs/usage_cli.md).

Development: `pip install -e "cli[test]" && pytest cli/tests`. The parity tests need Node.js.
After changing `assets/schema/*.json`, run `python cli/tools/sync_schemas.py`.
