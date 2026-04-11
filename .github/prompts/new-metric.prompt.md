---
description: "Add a new evaluation metric to the pipeline: compute in Python, wire through ted_core.py, add to evaluation TSV output"
---

# New Evaluation Metric

Add a new metric to the flair-eval evaluation pipeline.

## What does the metric measure?
${input:description}

## Steps

1. **Implement computation** in appropriate module under `bin/evaluation/`
   - Import shared parsing from `signal_utils.py`
   - Return a metrics dict (key → numeric value)
   - Add type annotations and docstring

2. **Integrate in `ted_core.py`**
   - Import the new function
   - Call within `timed_section()` block
   - `metrics.update(new_metrics)` to include in TSV output

3. **Export** from `bin/evaluation/__init__.py` (import + `__all__`)

4. **Add plot** (optional) — create plot function in `plots.py` or standalone script

5. **Verify**
   ```bash
   python -c "import py_compile; py_compile.compile('bin/evaluation/module.py', doraise=True)"
   conda run -n flair-refactor-ted python -c "import sys; sys.path.insert(0,'bin'); from evaluation.module import func; print('OK')"
   ```
