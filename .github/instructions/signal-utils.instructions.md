---
description: "Use when modifying signal_utils.py — the canonical shared parser for BED12, GTF, signal tracks, and junction chain grouping. All evaluation scripts must import from here."
applyTo: "bin/evaluation/signal_utils.py"
---

# signal_utils.py — Shared Parsing API

This is the **single source of truth** for BED12/GTF parsing across the entire pipeline. Never duplicate this logic elsewhere.

## When Adding New Functions

1. Add full type annotations and docstring
2. Update the module docstring "Provides:" list at the top
3. Export from `bin/evaluation/__init__.py` (both the import block and `__all__` list)
4. Verify: `python -c "import py_compile; py_compile.compile('bin/evaluation/signal_utils.py', doraise=True)"`

## Junction Format

Junctions are **2-tuples** `(donor, acceptor)` — genomic coordinates only. Chrom/strand context comes from the parent isoform dict. Code that needs globally unique junctions should prefix with `(chrom, strand)`.

## Dict Schema

Every `parse_bed12()` dict must include: chrom, start, end, name, score, strand, tss, tts, junctions, n_exons, spliced_len.
