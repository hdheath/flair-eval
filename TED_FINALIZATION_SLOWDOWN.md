# TED `--ted` finalization slowdown — FIXED 2026-06-02

**Status:** Fix APPLIED (working tree, branch `ted-port-52826`, uncommitted). Author/owner: hdheath.
**Severity:** High — silently dropped TED outputs on full-genome runs.

**History note:** The `diagnose_zero_denominator_genes` / `guard_zero_denominators` machinery
(`ted_diag.py` + the `--ted` block) was an **uncommitted** working-tree patch on top of `fad8b34`;
it exists in no committed TED version. The last committed version (`fad8b34`) had *no* guard and would
have hit a `ZeroDivisionError` in stock `_iso_passes_support_filter` on this data (TED snapping creates
the zero-denominator condition). So the prior version didn't "work" — it would have crashed. The patch
fixed the crash but did it with a per-partition full-file re-read (the diagnostic), which is the
runtime regression this change removes while KEEPING the cheap crash guard.

---

## Symptom

Full-genome `ASSEMBLY` run `full_4lib_0531` (default workflow, all caching). Run summary reported
`Succeeded: 8 / Failed: 3`, but **3 of the 4 TED transcriptomes were missing** from the results. Only
`WTC11_cDNA` produced a TED isoform set; `A549_cDNA`, `A549_ONT_cDNA`, `A549_dRNA` did not.

The 3 missing TED jobs were `flair transcriptome --ted` SLURM jobs that were **killed at the 24 h wall
limit (exit 140)**. Because `FlairTranscriptome` has `errorStrategy = 'ignore'`
(`nextflow.config:98`), the failures were swallowed and the rows silently dropped — hence the
misleading "Succeeded" summary.

| Dataset | chemistry | input BAM | mapped reads | TED result |
|---|---|---|---|---|
| A549_dRNA | ont_dRNA | 23.3 GB | (≫20M) | killed @24h — **HANG** (Mode B) |
| A549_ONT_cDNA | ont_cDNA | 7.0 GB | 23.7 M | killed @24h — slow (Mode A) |
| A549_cDNA | pacbio_cDNA | 2.4 GB | 11.1 M | killed @24h — slow (Mode A) |
| WTC11_cDNA | pacbio_cDNA | 3.1 GB | 8.8 M | ✅ ~22h |

Note WTC11 is **not** the largest input — it has the *fewest* reads, yet it's the one that finished.
So this is **not** a raw-data-scale problem. It is TED-path specific: plain (non-`--ted`) FLAIR runs
complete in far less time because they never enter the code path below.

---

## Root cause

In `_run_region`, after `calc_final_iso_support` already reads `isoform.ends.tsv` once, the `--ted`
block runs a **per-partition diagnostic that re-reads the entire `isoform.ends.tsv` again**:

`flair_transcriptome.py:1570-1588`
```python
if getattr(args, 'ted', False):
    try:
        from flair.ted_diag import (diagnose_zero_denominator_genes,
                                    guard_zero_denominators)
        diagnose_zero_denominator_genes(            # <-- pure logging; re-reads the whole ends file
            partition.output_path('isoform.ends.tsv'),
            final_transcript_objs, iso_to_counts, gene_to_tot)
        guard_zero_denominators(final_transcript_objs, gene_to_tot)   # <-- the real, cheap fix
    except Exception as _diag_err:
        logging.warning('[ted_diag] diagnostic/guard skipped: %s', _diag_err)
```

`diagnose_zero_denominator_genes` (`ted_diag.py:20-122`) only does expensive work when there is at
least one zero-denominator gene (`if not offenders: return`, line 52). It was written on the
assumption that **zero-denominator genes are rare**. They are not: TED end-snapping reassigns reads so
their terminal junction indices no longer match their isoform, which produces a zero full-length
denominator for many genes. Observed offending partitions per job: **A549_ONT 247, WTC11 222,
A549_cDNA 91, A549_dRNA 295**. So the "rare" branch fires on **nearly every partition**, and each time
it re-opens and re-scans the full `isoform.ends.tsv` from disk (`ted_diag.py:60-86`) just to print a
report capped at 25 genes / 12 reads.

Net effect with `--ted`: `isoform.ends.tsv` is read **twice per partition** instead of once. Noisier
chemistries (ONT/dRNA: ~3,615 / ~1,527 TED peaks on chr22 vs WTC11 pacbio's ~180) generate more
fragmented reads → more zero-denominator genes → more partitions hitting the redundant second read.
This is the ~2× finalization I/O that pushes the slower datasets past 24 h.

**This is a diagnostic-only function.** Its own docstring: *"Pure logging; nothing here changes
pipeline output."* The call-site comment even flags it as temporary: *"This is the
investigation-safe default; revisit once the diag log is read."* The diag log has now been read.

### Secondary: `ERROR: not full-length` spam (same hot path)

`calc_final_iso_support` prints `ERROR: not full-length` once **per non-full-length read**
(`flair_transcriptome.py:1444`). Verified pure side-effect: that `else:` branch only prints — it
increments no tally and does not `continue`/`break`; lines 1445-1446 run regardless. It produced
0.5M–2.5M lines per job (18–57 MB logs). Not the primary timeout driver (block-buffered to NFS), but
it's noise on the same hot loop and should be silenced in the same change.

---

## Fix (behavior-neutral)

### 1. Disable the per-partition diagnostic; keep the guard

`flair_transcriptome.py:1570-1588` → replace the `--ted` block with:
```python
if getattr(args, 'ted', False):
    try:
        from flair.ted_diag import guard_zero_denominators
        # diagnose_zero_denominator_genes() is diagnostic-only and re-reads the
        # ENTIRE isoform.ends.tsv per partition. Zero-denominator genes are NOT
        # rare under TED snapping (91-295 partitions/job), so this fired on nearly
        # every partition and ~doubled finalization I/O -> 24h SLURM kills.
        # Gate it behind an explicit debug flag; the guard below is the real fix.
        if getattr(args, 'ted_debug_zero_denom', False):
            from flair.ted_diag import diagnose_zero_denominator_genes
            diagnose_zero_denominator_genes(
                partition.output_path('isoform.ends.tsv'),
                final_transcript_objs, iso_to_counts, gene_to_tot)
        guard_zero_denominators(final_transcript_objs, gene_to_tot)
    except Exception as _diag_err:
        logging.warning('[ted_diag] guard skipped: %s', _diag_err)
```
(If no debug flag is wanted, just delete the `diagnose_*` import and call entirely.)

**Why safe:** `guard_zero_denominators` is unchanged and still prevents the `ZeroDivisionError`. The
removed function is pure logging and alters no output.

### 2. Silence the per-read print

`flair_transcriptome.py:1444` → delete the `print('ERROR: not full-length')` (or replace the whole
loop's spam with a single per-partition counter logged once). The `else:` branch then becomes a
no-op `pass`, which is correct — these reads are intentionally not counted as full-length.

**Why safe:** the branch had no effect beyond printing; tallies are unchanged.

---

## Expected impact

- Eliminates the redundant second full scan of `isoform.ends.tsv` on every offending partition
  (~2× → 1× finalization read). Expected to bring **A549_cDNA and A549_ONT_cDNA under 24 h**.
- Removes 0.5M–2.5M log lines / 18–57 MB per job.
- **Outputs byte-for-byte identical** (both changes are logging-only).

## Not fixed by this change

- **A549_dRNA HANG (Mode B):** its log froze ~11.75 h in and emitted nothing for ~12 h until the wall
  kill — a stuck/oversized partition, not steady slowness. Removing the diagnostic reduces its load
  but the hang needs separate diagnosis (instrument `_run_region` to log per-partition start/end +
  `isoform.ends.tsv` size, rerun with a short `-t`). More wall-time alone would just hang for 48 h.

## Cheapest confirmation before committing 24h×12-core compute

Rerun **only A549_dRNA and A549_cDNA** with `-t 3:00:00` and `PYTHONUNBUFFERED=1`, plus one
instrumentation line in `_run_region` (~`flair_transcriptome.py:1564`) logging
`partition X: ends.tsv=N lines, start/end timestamp`. Costs a few core-hours and decisively shows:
(a) A549_dRNA's exact stall partition (Mode B), (b) A549_cDNA steadily clearing partitions (Mode A).

---

## Pipeline-side follow-ups (reporter will handle; no TED source change)

- Change `FlairTranscriptome` `errorStrategy` from `'ignore'` to `'retry'`/`'terminate'`
  (`nextflow.config:98`) so TED failures surface instead of silently dropping rows.
- Add `PYTHONUNBUFFERED=1` to the process env so future hangs show the true last point of progress.
