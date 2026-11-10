"""
BED file I/O and bedtools wrapper functions.

Provides utilities for reading, writing, and manipulating BED format files,
as well as wrappers for common bedtools operations.
"""

import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from .utils import run, timed_section, get_logger

logger = get_logger()


def read_bed6(path: Path) -> List[dict]:
    """Read peaks as BED6 if possible; tolerate BED3 by adding Strand='.'

    Captures Name (col 4) and Score (col 5) when available.
    Score is stored as float to support TPM values.
    """
    rows = []
    try:
        with open(path) as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if not row or row[0].startswith('#'):
                    continue
                if len(row) >= 6:
                    try:
                        score = float(row[4])
                    except (ValueError, IndexError):
                        score = 0.0
                    rows.append({
                        'Chrom': row[0],
                        'Start': int(row[1]),
                        'End': int(row[2]),
                        'Name': row[3],
                        'Score': score,
                        'Strand': row[5],
                    })
                elif len(row) >= 3:
                    rows.append({
                        'Chrom': row[0],
                        'Start': int(row[1]),
                        'End': int(row[2]),
                        'Name': '',
                        'Score': 0.0,
                        'Strand': '.',
                    })
        return rows
    except Exception as e:
        logger.warning(f"Error reading {path}: {e}")
        return []


def prepare_bed6_sorted(path: Path, tmps: list) -> Path:
    """Trim to BED6 and sort."""
    trimmed = path.with_suffix(path.suffix + ".trimmed.tmp")
    sorted_p = path.with_suffix(path.suffix + ".sorted.tmp")

    res = run(["cut", "-f", "1-6", str(path)])
    trimmed.write_text(res.stdout)

    res = run(["bedtools", "sort", "-i", str(trimmed)])
    sorted_p.write_text(res.stdout)

    try:
        trimmed.unlink()
    except Exception:
        pass
    tmps.append(sorted_p)
    return sorted_p


def extract_tss_tts_bed(iso_bed: Path, tmps: list) -> Tuple[Path, Path]:
    """
    Extract single-position TSS and TTS coordinates from isoforms BED.

    For + strand: TSS = start, TTS = end
    For - strand: TSS = end, TTS = start

    Returns (tss_bed, tts_bed) as single-position BED6 files.
    """
    tss_rows = []
    tts_rows = []

    with open(iso_bed) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            if len(row) < 6:
                continue

            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            name = row[3]
            score = row[4]
            strand = row[5]

            # Calculate single-position TSS/TTS
            if strand == '+':
                tss_pos = start
                tts_pos = end
            else:  # strand == '-'
                tss_pos = end
                tts_pos = start

            # Create 1bp windows for TSS
            tss_rows.append([chrom, str(tss_pos), str(tss_pos + 1), name, score, strand])
            # Create 1bp windows for TTS
            tts_rows.append([chrom, str(tts_pos), str(tts_pos + 1), name, score, strand])

    # Write TSS BED
    tss_bed = iso_bed.with_suffix(".tss.tmp")
    with open(tss_bed, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(tss_rows)
    tmps.append(tss_bed)

    # Write TTS BED
    tts_bed = iso_bed.with_suffix(".tts.tmp")
    with open(tts_bed, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(tts_rows)
    tmps.append(tts_bed)

    logger.debug(f"Extracted {len(tss_rows)} TSS and TTS positions from {iso_bed}")
    return tss_bed, tts_bed


def run_closest(a: Path, b: Path) -> List[List[str]]:
    """Run bedtools closest."""
    with timed_section(f"bedtools_closest"):
        res = run(["bedtools", "closest", "-a", str(a), "-b", str(b), "-s", "-d"])
        rows = []
        for line in res.stdout.strip().split('\n'):
            if line and not line.startswith('#'):
                rows.append(line.split('\t'))
        return rows


def extract_distance_and_peak(rows: List[List[str]], label: str, max_dist: int) -> Dict[str, Tuple[int, str]]:
    """Extract distances to nearest peak. Returns {isoform_id: (dist, peak_id)}"""
    result = {}
    for row in rows:
        if not row or len(row) < 13:
            continue
        tx_id = row[3]
        try:
            dist = int(row[-1])
        except (ValueError, IndexError):
            dist = max_dist + 1

        try:
            peak_id = f"{row[6]}_{row[7]}_{row[8]}"
        except IndexError:
            peak_id = "."

        # Keep first occurrence
        if tx_id not in result:
            result[tx_id] = (dist, peak_id)

    return result


def extract_signed_distances(rows: List[List[str]]) -> List[int]:
    """
    Extract signed distances from bedtools closest output for ALL transcripts.

    Signed distance convention:
    - Negative = transcript end is upstream (5' direction) of experimental peak
    - Positive = transcript end is downstream (3' direction) of experimental peak

    For + strand: signed_dist = transcript_pos - peak_pos
    For - strand: signed_dist = peak_pos - transcript_pos (flipped due to orientation)

    Args:
        rows: Output from bedtools closest (columns: transcript BED6, peak BED6, distance)

    Returns:
        List of signed distances for all transcripts (no filtering applied)
    """
    distances = []
    for row in rows:
        if not row or len(row) < 13:
            continue

        # Parse transcript position (columns 0-5 are transcript BED6)
        try:
            tx_start = int(row[1])
            tx_end = int(row[2])
            tx_strand = row[5]
        except (ValueError, IndexError):
            continue

        # Parse peak position (columns 6-11 are peak BED6)
        try:
            peak_start = int(row[7])
            peak_end = int(row[8])
        except (ValueError, IndexError):
            continue

        # Check if no peak was found (bedtools closest reports -1 or "." for no match)
        try:
            abs_dist = int(row[-1])
            if abs_dist == -1:
                continue
        except (ValueError, IndexError):
            continue

        # Calculate midpoints for more accurate distance
        tx_pos = (tx_start + tx_end) // 2
        peak_pos = (peak_start + peak_end) // 2

        # Calculate signed distance based on strand
        if tx_strand == '+':
            signed_dist = tx_pos - peak_pos
        else:  # '-' strand
            signed_dist = peak_pos - tx_pos

        distances.append(signed_dist)

    return distances


def vectorized_overlap_counts(
    bed_rows: List[dict],
    peaks_rows: List[dict],
    window: int,
    signal=None,
) -> Tuple[int, int]:
    """Return (#peak summits matched one-to-one within window, total peaks).

    Thin back-compat wrapper around :func:`one_to_one_match_counts` that
    returns just the recall (numerator, denominator) pair.  Prefer
    one_to_one_match_counts when you also need the precision numerator /
    call count (W1 — so precision and recall share the same TP).
    """
    tp, _n_calls, n_peaks = one_to_one_match_counts(
        bed_rows, peaks_rows, window, signal=signal)
    return tp, n_peaks


def one_to_one_match_counts(
    bed_rows: List[dict],
    peaks_rows: List[dict],
    window: int,
    signal=None,
) -> Tuple[int, int, int]:
    """Return (TP, n_calls, n_peaks) under summit + one-to-one matching.

    Contract-compliant accounting (the R3 fix), shared by precision AND
    recall in ted_core.tss_tts_metrics (the W1 unification):

      * each peak interval collapses to ONE summit — the dominant-signal
        position when a `signal` track is supplied, else the peak centre,
      * called ends match summits ONE-TO-ONE (greedy nearest), so a fat
        cluster of isoform ends near one peak credits exactly one TP, not
        one per call.

    Because matching is one-to-one, the number of matched calls equals the
    number of matched summits — so a single TP count feeds both:
        precision = TP / n_calls      recall = TP / n_peaks

    Strand-aware; strand-agnostic ('.') peaks match both strands.  The
    single-position TSS/TTS endpoint convention from extract_tss_tts_bed
    means bed_row Start == the called-end coordinate (End == Start+1).

    `signal`, if given, maps (chrom, strand) -> sorted [(pos, value)] so the
    summit is the dominant CTSS / density-peak; otherwise the interval
    centre is used (graceful fallback).
    """
    n_calls = len(bed_rows)
    if not peaks_rows:
        return 0, n_calls, 0

    from .peak_matching import match_positions_one_to_one, _interval_summit

    # Collapse each peak to a single summit; bucket by (chrom, strand).
    #
    # A strand-agnostic ('.') peak is placed in EXACTLY ONE bucket — its own
    # '.' bucket — so it contributes a single summit (one recall target).  The
    # old code placed it in BOTH the '+' and '-' buckets but counted it once in
    # the denominator, so a '+' call and a '-' call near one '.' peak each
    # claimed it (TP=2 for a single peak) and recall could exceed 1.0.  Here
    # each peak contributes exactly one summit, so n_summits == #peak intervals.
    summits_by_cs = defaultdict(list)  # (chrom, strand) -> [summit_pos, ...]
    for row in peaks_rows:
        strand = row['Strand']
        sig = signal.get((row['Chrom'], strand)) if signal else None
        # '.' peaks have no per-strand signal track; fall back to '+' then '-'.
        if sig is None and signal and strand == '.':
            sig = (signal.get((row['Chrom'], '+'))
                   or signal.get((row['Chrom'], '-')))
        summit = _interval_summit(int(row['Start']), int(row['End']), sig)
        summits_by_cs[(row['Chrom'], strand)].append(summit)
    for key in summits_by_cs:
        summits_by_cs[key].sort()
    # One summit per peak interval — the recall denominator.
    n_summits = sum(len(v) for v in summits_by_cs.values())

    # Group called ends (single-position endpoints) by (chrom, strand).
    ends_by_cs = defaultdict(list)
    for row in bed_rows:
        ends_by_cs[(row['Chrom'], row['Strand'])].append(int(row['Start']))

    matched_keys = set()  # (chrom, strand, summit) — distinct matched summits
    for (chrom, strand), positions in ends_by_cs.items():
        # A call may match a summit on its own strand OR a strand-agnostic
        # ('.') summit on the same chrom.  Merge both candidate summit sets
        # (tagging each with its source bucket) and run ONE one-to-one match
        # so the call is consumed once and each summit is credited once —
        # preserving precision = TP/n_calls <= 1 and recall = TP/n_summits <= 1.
        candidates: List[Tuple[int, str]] = []
        for bucket_strand in (strand, '.'):
            if bucket_strand == '.' and strand == '.':
                continue  # '.' calls already covered by their own bucket below
            for s in summits_by_cs.get((chrom, bucket_strand), ()):
                candidates.append((s, bucket_strand))
        # '.' calls match only '.' summits (no strand to disambiguate).
        if strand == '.':
            candidates = [(s, '.') for s in summits_by_cs.get((chrom, '.'), ())]
        if not candidates:
            continue
        candidates.sort()
        cand_pos = [p for p, _ in candidates]
        assign = match_positions_one_to_one(positions, cand_pos, window)
        for j in assign:
            if j >= 0:
                summit_pos, bucket_strand = candidates[j]
                matched_keys.add((chrom, bucket_strand, summit_pos))

    tp = len(matched_keys)
    # Clamp TP to the summit count as a final guard so recall (TP/n_summits)
    # can never exceed 1.0 regardless of bucketing.
    if tp > n_summits:
        tp = n_summits
    return tp, n_calls, n_summits


def near_match_counts(bed_rows, peaks_rows, window, signal=None):
    """MANY-TO-ONE 'near a peak' counts (NO one-to-one consumption).

    Returns (ends_near, n_calls, peaks_recovered, n_peaks):
      ends_near       = # called ends within `window` of ANY peak summit  (precision numerator)
      peaks_recovered = # peak summits with ANY called end within `window` (recall numerator)
    A promoter producing several isoforms therefore credits EVERY nearby end
    (hdheath 2026-06-15: alt-TSS within one promoter are expected, not duplicates).
    Summit construction mirrors one_to_one_match_counts (peak centre / dominant
    signal, strand-aware) so only the consumption rule differs.
    """
    from bisect import bisect_left
    from collections import defaultdict
    from .peak_matching import _interval_summit
    n_calls = len(bed_rows)
    if not peaks_rows:
        return 0, n_calls, 0, 0
    summits_by_cs = defaultdict(list)
    for row in peaks_rows:
        strand = row['Strand']
        sig = signal.get((row['Chrom'], strand)) if signal else None
        if sig is None and signal and strand == '.':
            sig = (signal.get((row['Chrom'], '+')) or signal.get((row['Chrom'], '-')))
        summits_by_cs[(row['Chrom'], strand)].append(
            _interval_summit(int(row['Start']), int(row['End']), sig))
    for k in summits_by_cs:
        summits_by_cs[k].sort()
    calls_by_cs = defaultdict(list)
    for row in bed_rows:
        calls_by_cs[(row['Chrom'], row['Strand'])].append(int(row['Start']))
    for k in calls_by_cs:
        calls_by_cs[k].sort()
    def _near(pos, arr):
        if not arr:
            return False
        i = bisect_left(arr, pos)
        for j in (i - 1, i):
            if 0 <= j < len(arr) and abs(pos - arr[j]) <= window:
                return True
        return False
    ends_near = 0
    for (c, st), poss in calls_by_cs.items():
        summ = summits_by_cs.get((c, st), [])
        ends_near += sum(1 for p in poss if _near(p, summ))
    n_peaks = sum(len(v) for v in summits_by_cs.values())
    peaks_recovered = 0
    for (c, st), summ in summits_by_cs.items():
        cl = calls_by_cs.get((c, st), [])
        peaks_recovered += sum(1 for p in summ if _near(p, cl))
    return ends_near, n_calls, peaks_recovered, n_peaks
