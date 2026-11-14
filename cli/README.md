# flair-eval (CLI)

Set up, check and launch [flair-eval](https://github.com/hdheath/flair-eval), the Nextflow
pipeline that compares FLAIR runs (versions, branches, settings) on the same data.

```bash
pip install "flair-eval @ git+https://github.com/hdheath/flair-eval@main#subdirectory=cli"
flair-eval doctor            # check Java, Nextflow, conda, disk; detects phoenix
flair-eval init              # questions -> samples.csv + runs.yaml
flair-eval check             # validate them against the pipeline's schemas
flair-eval run --outdir /private/groups/brookslab/$USER/flair-eval/dd_mm_yy_name
flair-eval report out1/bundles out2/bundles --outdir /private/groups/brookslab/$USER/flair-eval/dd_mm_yy_report
```

Guide: [docs/usage.md](https://github.com/hdheath/flair-eval/blob/main/docs/usage.md#the-flair-eval-cli).

## Development

```bash
pip install -e "cli[test]" && pytest cli/tests
```

The builder parity tests need Node.js. The page (`frontend/builder/logic.js`) and this package
implement the same rules; see [frontend/builder/README.md](../frontend/builder/README.md).
After changing `assets/schema/*.json`, run `python cli/tools/sync_schemas.py`.
