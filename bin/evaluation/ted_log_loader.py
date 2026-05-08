"""
ted_log_loader.py — single source of truth for loading flair_ted.ted_log.tsv.

Every diagnostic that consumes flair_ted.ted_log.tsv (ted_rejection_analysis,
ted_confusion_matrix, ted_log_analysis, ...) must go through these helpers.
They handle the schema differences between per-SJC and --ted_global mode:

  - per-SJC: junc_id has real (start, end) bounds, cluster_id is 1:1 with
    a chain-scoped cluster, stage in {ted_cluster, tss_snap, fallback},
    tss_pos and tts_pos always >= 0.
  - --ted_global: emits additional `global_peak` rows that are
    PARTITION-LEVEL per-axis summaries — junc_id like
    `chr*:None-None:strand:0j`, with one of (tss_pos, tts_pos) = -1
    (sentinel for "TSS-only peak" or "TTS-only peak"). These rows track
    which CAGE/dRNA-style peaks passed partition-level scoring; they
    are NOT per-isoform candidates and must be filtered before entering
    recall/precision/distance/violin plots, otherwise they leak into
    outputs as phantom FPs/FNs and NaN-filled distributions.

Reference implementation: analysis/ted_pr_sweep/ted_pr_sweep.py:178-190
(filters the same set of stages and rejects sentinel positions). The
logic is duplicated here rather than imported to keep this module
self-contained — it has zero relative imports so standalone scripts in
this directory (which Python invokes with `python /path/to/script.py`)
can sibling-import it without dragging in the rest of the `evaluation`
package's heavyweight dependencies.
"""

from __future__ import annotations


# Per-isoform candidate stages — anything else (currently just `global_peak`)
# is partition-level summary metadata, not a candidate isoform.
PER_ISOFORM_STAGES = ("ted_cluster", "global_assign",
                      "fallback", "tss_snap")


def is_global_mode(path) -> bool:
    """Return True if the ted_log.tsv was produced by a --ted_global run.

    Detection: scan the leading `# hdbscan_params: ...` comment line that
    ted.py prepends to the TSV; if it contains `ted_global=True` (case
    sensitive — that's how flair_transcriptome.py writes it) the run was
    in global mode.
    """
    try:
        with open(path) as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                if "ted_global=True" in line:
                    return True
    except (OSError, FileNotFoundError):
        pass
    return False


def load_ted_log(path, *, drop_partition_summary: bool = True):
    """Schema-aware loader for flair_ted.ted_log.tsv.

    drop_partition_summary=True (default) strips global-mode `global_peak`
    rows and any sentinel-coordinate rows. Pass False only when you
    explicitly want the partition-level peak-promotion outcomes (e.g.,
    auditing which CAGE peaks the validator rejected at the partition
    level — distinct from per-isoform reject reasons).
    """
    import pandas as pd  # lazy: callers may not need pandas at import time

    df = pd.read_csv(path, sep="\t", comment="#")

    if drop_partition_summary:
        if "stage" in df.columns:
            df = df[df["stage"].isin(PER_ISOFORM_STAGES)].copy()
        # Strip sentinel-coordinate rows even if `stage` column was absent
        # (older logs predate the column). Both axes must be set on a
        # per-isoform row by construction.
        if "tss_pos" in df.columns:
            df = df[df["tss_pos"].notna() & (df["tss_pos"] >= 0)].copy()
        if "tts_pos" in df.columns:
            df = df[df["tts_pos"].notna() & (df["tts_pos"] >= 0)].copy()

    # Cast numeric columns once so downstream callers don't have to
    # repeat the boilerplate.
    for c in ("n_reads", "jc_n_reads_total",
              "tss_pos", "tts_pos",
              "tss_spread_iqr", "tts_spread_iqr",
              "TED_confidence",
              "TED_tss_reality", "TED_tts_reality",
              "TED_tss_depth", "TED_tts_depth", "TED_depth",
              "TED_tss_model", "TED_tts_model",
              "TED_tss_annot", "TED_tts_annot",
              "TED_tss_annot_dist", "TED_tts_annot_dist",
              "threshold_tss", "threshold_tts"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df
