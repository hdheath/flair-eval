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
    """Read peaks as BED6 if possible; tolerate BED3 by adding Strand='.'"""
    rows = []
    try:
        with open(path) as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if not row or row[0].startswith('#'):
                    continue
                if len(row) >= 6:
                    rows.append({
                        'Chrom': row[0],
                        'Start': int(row[1]),
                        'End': int(row[2]),
                        'Strand': row[5]
                    })
                elif len(row) >= 3:
                    rows.append({
                        'Chrom': row[0],
                        'Start': int(row[1]),
                        'End': int(row[2]),
                        'Strand': '.'
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


def vectorized_overlap_counts(bed_rows: List[dict], peaks_rows: List[dict], window: int) -> Tuple[int, int]:
    """Return (#peaks matched by ≥1 isoform within window, total peaks)."""
    if not peaks_rows:
        return 0, 0

    consumed = set()

    # Group bed rows by chrom/strand
    bed_by_cs = defaultdict(list)
    for row in bed_rows:
        key = (row['Chrom'], row['Strand'])
        bed_by_cs[key].append(row)

    # Group peaks by chrom/strand
    peaks_by_cs = defaultdict(list)
    for idx, row in enumerate(peaks_rows):
        key = (row['Chrom'], row['Strand'])
        peaks_by_cs[key].append((idx, row))
        if row['Strand'] == '.':
            # Strand-agnostic peaks go to both strands
            for strand in ['+', '-']:
                peaks_by_cs[(row['Chrom'], strand)].append((idx, row))

    # Check overlaps
    for (c, s), bed_chunk in bed_by_cs.items():
        if (c, s) not in peaks_by_cs:
            continue

        peaks_chunk = peaks_by_cs[(c, s)]
        # Sort peaks by start
        peaks_chunk_sorted = sorted(peaks_chunk, key=lambda x: x[1]['Start'])

        for bed_row in bed_chunk:
            lo = bed_row['Start'] - window
            hi = bed_row['End'] + window

            for idx, peak in peaks_chunk_sorted:
                if peak['Start'] > hi:
                    break
                if peak['End'] >= lo:
                    consumed.add(idx)

    return len(consumed), len(peaks_rows)
