---
description: "Use when creating or editing Nextflow process definitions. Covers publishDir patterns, errorStrategy, input/output tuples, and Python script calling conventions."
applyTo: "modules/**/main.nf"
---

# Nextflow Process Conventions

- Process names: CamelCase
- Visualization/summary processes: `errorStrategy 'ignore'`, `optional: true` on plot outputs
- Always include `tag` for log identification
- Python scripts: `python ${projectDir}/bin/evaluation/script.py`
- BED file args use `label:path` format — build in Groovy loop:
  ```groovy
  def bed_args = []
  for (int i = 0; i < bed_labels.size(); i++) {
      bed_args << "${bed_labels[i]}:${bed_files[i]}"
  }
  ```
- End script blocks with `|| true` for non-critical processes
