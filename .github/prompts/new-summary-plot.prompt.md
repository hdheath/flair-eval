---
description: "Create a new summary visualization process: standalone Python plot script + Nextflow process + subworkflow wiring"
---

# New Summary Plot

Create a new cross-method summary plot for the flair-eval pipeline.

## What does the plot show?
${input:description}

## What data does it need?
${input:data_source}

## Steps

1. Create standalone Python script at `bin/evaluation/${input:script_name}.py`
   - CLI: `--bed label:path [...] --output dir/`
   - Import parsing from `signal_utils`
   - Use pub_style conventions or inline matplotlib with Okabe-Ito palette
   - Self-contained plotting (no import from plots.py)

2. Add Nextflow process `${input:process_name}` to `modules/visualization/summary/main.nf`
   - Uses same input channel as `IsoformsPerGeneHist` (label:path BED pairs)
   - publishDir to `summary/${input:output_folder}/`

3. Wire in `subworkflows/summary_and_viz.nf`
   - Add `include` statement
   - Call process with appropriate channel

4. Verify syntax: `python -c "import py_compile; py_compile.compile('bin/evaluation/${input:script_name}.py', doraise=True)"`
