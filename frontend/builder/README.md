# Config-builder page

A static page for setting up a flair-eval launch: samples, FLAIR run cards, and the exact
`nextflow run` command, with `samples.csv` / `runs.yaml` to download. No build step.

| file | what |
|---|---|
| `index.html` | the builder: markup and styles |
| `welcome.html` | GitHub Pages start page (served as `/`; links `builder/` and `demo/report.html`) |
| `app.js` | page wiring (form, panel, copy/download, browser-only autosave) |
| `logic.js` | all rules: `--flair` tokens, validation, file writers, command builder. Mirrors `cli/src/flair_eval/` |
| `schemas.js` | generated from `assets/schema/*.json` by `cli/tools/sync_schemas.py`. Do not edit |

Open `index.html` directly in a browser, or let `.github/workflows/pages.yml` publish it to
GitHub Pages as `/builder/`, next to `welcome.html` at `/` and a demo report at `/demo/`. `cli/tests/test_builder_parity.py` checks `logic.js` against the Python package.
Usage: [docs/usage_cli.md](../../docs/usage_cli.md#config-builder-page).
