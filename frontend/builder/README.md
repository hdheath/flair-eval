# Config-builder page

A static page for setting up a flair-eval launch: samples, FLAIR run cards, and the exact
`nextflow run` command, with `samples.csv` / `runs.yaml` to download. No build step.
How to use it: [docs/usage.md](../../docs/usage.md#config-builder-page).

| file | what |
|---|---|
| `index.html` | the builder: markup and styles |
| `welcome.html` | GitHub Pages start page (served as `/`; links `builder/` and `demo/report.html`) |
| `app.js` | page wiring (form, panel, copy/download, browser-only autosave) |
| `logic.js` | all rules: `--flair` tokens, validation, file writers, command builder. Mirrors `cli/src/flair_eval/` |
| `schemas.js` | generated from `assets/schema/*.json` by `cli/tools/sync_schemas.py`. Do not edit |

## Publishing

`.github/workflows/pages.yml` runs the CLI and page tests, then publishes
<https://hdheath.github.io/flair-eval/>: `welcome.html` at `/`, this page at `/builder/`, and a
demo report built from `assets/demo/bundles` at `/demo/report.html`. It runs on pushes to `main`
that touch the front ends, the schemas, the report code or the demo bundles (the full path list
is in the workflow).

## Keeping the page and the CLI in step

`logic.js` and `cli/src/flair_eval/` implement the same rules.
- Genome presets: `GENOME_PRESETS` in `cli/src/flair_eval/phoenix.py`, in `logic.js`, and
  `params.genomes` in `nextflow.config`. Tests check that all three match.
- `cli/tests/test_builder_parity.py` runs both on the same cases (needs Node) and fails if they
  disagree.
- Both read copies of `assets/schema/*.json`. After a schema change run
  `python cli/tools/sync_schemas.py` (CI runs it with `--check`).
