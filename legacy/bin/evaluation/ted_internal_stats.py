"""
TED internal algorithm statistics aggregator.

Reads FLAIR's per-locus TED decision log (TSV produced by --ted_log) and
aggregates the locus_summary rows into run-level metrics.  This moves the
summary computation OUT of FLAIR and INTO the eval pipeline so all
metrics live in one place.

Returns a flat dict with ``ted_internal_`` prefixed keys, ready to merge
into the main evaluation metrics dict.
"""

import csv
from pathlib import Path
from typing import Dict, Optional


def _safe_float(val, default=0.0):
    """Convert a string to float, returning *default* on failure."""
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def _safe_int(val, default=0):
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return default


def _safe_bool(val):
    return str(val).strip().lower() in ('true', '1', 'yes')


def _pct(num, denom):
    return f'{100 * num / denom:.2f}' if denom > 0 else '0.00'


def aggregate_ted_log(ted_log_path: str) -> Dict[str, str]:
    """Read a FLAIR per-locus TED decision log and aggregate run-level stats.

    Parameters
    ----------
    ted_log_path : str or Path
        Path to the per-locus TSV (produced by ``flair transcriptome --ted --ted_log``).

    Returns
    -------
    dict
        Flat dict with ``ted_internal_`` prefixed keys.  All values are strings
        (ready for TSV output).  Returns empty dict if the file is missing or
        has no locus_summary rows.
    """
    path = Path(ted_log_path)
    if not path.exists() or path.stat().st_size == 0:
        return {}

    # ── Accumulators (mirrors _ted_run_stats in flair_transcriptome.py) ──
    n_loci_processed = 0
    n_loci_with_output = 0
    n_loci_all_discarded = 0
    n_loci_fallback = 0
    n_loci_single_cluster_tss = 0
    n_loci_single_cluster_tts = 0
    total_tss_clusters_found = 0
    total_tts_clusters_found = 0
    total_tss_clusters_real = 0
    total_tts_clusters_real = 0
    total_reads_processed = 0
    total_reads_real_tss = 0
    total_reads_artifact_tss = 0
    total_reads_real_tts = 0
    total_reads_artifact_tts = 0
    total_noise_tss = 0
    total_noise_tts = 0
    total_reads_kept = 0
    total_reads_reassigned = 0
    total_reads_discarded = 0
    total_ref_extensions = 0
    total_ref_ext_5prime = 0
    total_ref_ext_3prime = 0
    total_isoforms_output = 0

    # ── Also extract config from the first locus_summary row ──
    config = {}

    with open(path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row.get('step') != 'locus_summary':
                continue

            n_loci_processed += 1
            # Use n_reads_total (individual reads) when available,
            # falling back to n_reads (input items) for older logs.
            n_reads = _safe_int(row.get('n_reads_total',
                                        row.get('n_reads', 0)))
            total_reads_processed += n_reads

            n_iso_out = _safe_int(row.get('n_isoforms_output', 0))
            if n_iso_out > 0:
                n_loci_with_output += 1
            else:
                n_loci_all_discarded += 1

            # Detect fallback: use explicit flag (newer logs) or n_reads < 4
            if _safe_bool(row.get('_is_fallback', 'False')) or n_reads < 4:
                n_loci_fallback += 1

            if _safe_bool(row.get('single_cluster_tss', 'False')):
                n_loci_single_cluster_tss += 1
            if _safe_bool(row.get('single_cluster_tts', 'False')):
                n_loci_single_cluster_tts += 1

            total_tss_clusters_found += _safe_int(row.get('n_tss_clusters_found', 0))
            total_tts_clusters_found += _safe_int(row.get('n_tts_clusters_found', 0))
            total_tss_clusters_real += _safe_int(row.get('n_tss_clusters_real', 0))
            total_tts_clusters_real += _safe_int(row.get('n_tts_clusters_real', 0))

            total_reads_real_tss += _safe_int(row.get('n_reads_real_tss', 0))
            total_reads_artifact_tss += _safe_int(row.get('n_reads_artifact_tss', 0))
            total_reads_real_tts += _safe_int(row.get('n_reads_real_tts', 0))
            total_reads_artifact_tts += _safe_int(row.get('n_reads_artifact_tts', 0))
            total_noise_tss += _safe_int(row.get('n_noise_tss', 0))
            total_noise_tts += _safe_int(row.get('n_noise_tts', 0))

            total_reads_kept += _safe_int(row.get('n_reads_kept', 0))
            total_reads_reassigned += _safe_int(row.get('n_reads_reassigned', 0))
            total_reads_discarded += _safe_int(row.get('n_reads_discarded', 0))

            total_ref_extensions += _safe_int(row.get('n_ref_extensions', 0))
            total_ref_ext_5prime += _safe_int(row.get('n_ref_ext_5prime', 0))
            total_ref_ext_3prime += _safe_int(row.get('n_ref_ext_3prime', 0))

            total_isoforms_output += n_iso_out

            # Capture config from first row
            if not config:
                config = {
                    'cluster_mode': row.get('cluster_mode', ''),
                    'normalization': row.get('normalization', ''),
                    'noise_handling': row.get('noise_handling', ''),
                    'allow_single_cluster': row.get('allow_single_cluster', ''),
                    'cluster_selection': row.get('cluster_selection', ''),
                    'epsilon': row.get('epsilon', ''),
                    'threshold': row.get('threshold', ''),
                }

    if n_loci_processed == 0:
        return {}

    # ── Build output dict ────────────────────────────────────────────────
    n_loci = n_loci_processed
    total_fate = total_reads_kept + total_reads_reassigned + total_reads_discarded
    tss_total = total_reads_real_tss + total_reads_artifact_tss + total_noise_tss
    tts_total = total_reads_real_tts + total_reads_artifact_tts + total_noise_tts

    result = {}

    # Config
    for k, v in config.items():
        result[f'ted_internal_config_{k}'] = v

    # Locus statistics
    result['ted_internal_n_loci_processed'] = str(n_loci)
    result['ted_internal_n_loci_with_output'] = str(n_loci_with_output)
    result['ted_internal_n_loci_all_discarded'] = str(n_loci_all_discarded)
    result['ted_internal_n_loci_fallback'] = str(n_loci_fallback)
    result['ted_internal_pct_loci_with_output'] = _pct(n_loci_with_output, n_loci)
    result['ted_internal_pct_loci_all_discarded'] = _pct(n_loci_all_discarded, n_loci)
    result['ted_internal_pct_loci_fallback'] = _pct(n_loci_fallback, n_loci)

    # Cluster topology
    result['ted_internal_n_loci_single_cluster_tss'] = str(n_loci_single_cluster_tss)
    result['ted_internal_n_loci_single_cluster_tts'] = str(n_loci_single_cluster_tts)
    result['ted_internal_pct_loci_single_cluster_tss'] = _pct(n_loci_single_cluster_tss, n_loci)
    result['ted_internal_pct_loci_single_cluster_tts'] = _pct(n_loci_single_cluster_tts, n_loci)
    result['ted_internal_total_tss_clusters_found'] = str(total_tss_clusters_found)
    result['ted_internal_total_tts_clusters_found'] = str(total_tts_clusters_found)
    result['ted_internal_total_tss_clusters_real'] = str(total_tss_clusters_real)
    result['ted_internal_total_tts_clusters_real'] = str(total_tts_clusters_real)
    result['ted_internal_avg_tss_clusters_per_locus'] = f'{total_tss_clusters_found / n_loci:.2f}'
    result['ted_internal_avg_tts_clusters_per_locus'] = f'{total_tts_clusters_found / n_loci:.2f}'
    result['ted_internal_pct_tss_clusters_validated'] = _pct(total_tss_clusters_real, total_tss_clusters_found)
    result['ted_internal_pct_tts_clusters_validated'] = _pct(total_tts_clusters_real, total_tts_clusters_found)

    # Read classification by end type
    result['ted_internal_total_reads_processed'] = str(total_reads_processed)
    result['ted_internal_total_reads_real_tss'] = str(total_reads_real_tss)
    result['ted_internal_total_reads_artifact_tss'] = str(total_reads_artifact_tss)
    result['ted_internal_total_noise_tss'] = str(total_noise_tss)
    result['ted_internal_pct_reads_real_tss'] = _pct(total_reads_real_tss, tss_total)
    result['ted_internal_pct_reads_artifact_tss'] = _pct(total_reads_artifact_tss, tss_total)
    result['ted_internal_pct_noise_tss'] = _pct(total_noise_tss, tss_total)
    result['ted_internal_total_reads_real_tts'] = str(total_reads_real_tts)
    result['ted_internal_total_reads_artifact_tts'] = str(total_reads_artifact_tts)
    result['ted_internal_total_noise_tts'] = str(total_noise_tts)
    result['ted_internal_pct_reads_real_tts'] = _pct(total_reads_real_tts, tts_total)
    result['ted_internal_pct_reads_artifact_tts'] = _pct(total_reads_artifact_tts, tts_total)
    result['ted_internal_pct_noise_tts'] = _pct(total_noise_tts, tts_total)

    # Read fate
    result['ted_internal_total_reads_kept'] = str(total_reads_kept)
    result['ted_internal_total_reads_reassigned'] = str(total_reads_reassigned)
    result['ted_internal_total_reads_discarded'] = str(total_reads_discarded)
    result['ted_internal_pct_reads_kept'] = _pct(total_reads_kept, total_fate)
    result['ted_internal_pct_reads_reassigned'] = _pct(total_reads_reassigned, total_fate)
    result['ted_internal_pct_reads_discarded'] = _pct(total_reads_discarded, total_fate)

    # Reference extension
    result['ted_internal_total_ref_extensions'] = str(total_ref_extensions)
    result['ted_internal_total_ref_ext_5prime'] = str(total_ref_ext_5prime)
    result['ted_internal_total_ref_ext_3prime'] = str(total_ref_ext_3prime)
    result['ted_internal_pct_loci_with_ref_extension'] = _pct(total_ref_extensions, n_loci)

    # Output
    result['ted_internal_total_isoforms_output'] = str(total_isoforms_output)
    result['ted_internal_avg_isoforms_per_locus'] = f'{total_isoforms_output / n_loci:.2f}'

    return result
