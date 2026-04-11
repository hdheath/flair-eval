---
description: "Run the flair-eval Nextflow pipeline with resume, check for failures, and re-run failed processes"
---

# Run Pipeline

Run or resume the flair-eval pipeline.

## Configuration
${input:params_file}

## Commands

```bash
# Activate environment
conda activate nextflow_env

# Run with resume (reuses cached results)
cd /private/groups/brookslab/hdheath/projects/flair-eval
nextflow run main.nf -params-file params/${input:params_file} -resume

# Check for failures
grep -r 'ERROR\|FAILED' .nextflow.log | tail -20

# List cached work directories
nextflow log last -f name,status,hash | grep -v COMPLETED
```

## Cache Invalidation

If a Python script was modified and you need to force re-run:
1. Find the process work dir: `nextflow log last -f name,hash | grep ProcessName`
2. Delete the cache: `rm -rf work/<hash_prefix>*`
3. Re-run with `-resume`
