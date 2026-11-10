# flair-eval (CLI)

Set up, check and launch [flair-eval](https://github.com/hdheath/flair-eval), the Nextflow
pipeline that compares FLAIR runs (versions, branches, settings) on the same data.

```bash
pip install -e cli/          # from a clone of the repository
flair-eval init              # questions -> samples.csv + runs.yaml
flair-eval check             # validate them against the pipeline's schemas
flair-eval run --outdir /big/storage/out
flair-eval report out1/bundles out2/bundles --outdir report_out
```

Full guide: [docs/usage_cli.md](../docs/usage_cli.md).

Development: `pip install -e "cli[test]" && pytest cli/tests`. The parity tests need Node.js.
After changing `assets/schema/*.json`, run `python cli/tools/sync_schemas.py`.
