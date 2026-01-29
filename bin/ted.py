#!/usr/bin/env python3
"""
TED (Transcriptome End Distribution) - Standard Library Version
Evaluates TSS/TTS accuracy and read assignment metrics for FLAIR outputs.
Uses only Python standard library + bedtools/samtools.
"""

import argparse
import csv
import io
import logging
import math
import os
import statistics
import subprocess
import sys
import tempfile
import time
import numpy as np
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# Optional matplotlib import for plotting
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend for server environments
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s'
)
logger = logging.getLogger(__name__)

# Performance timing storage
_TIMING_DATA = defaultdict(list)

@contextmanager
def timed_section(section_name: str):
    """Context manager to time code sections and store results."""
    start = time.time()
    try:
        yield
    finally:
        elapsed = time.time() - start
        _TIMING_DATA[section_name].append(elapsed)
        logger.debug(f"[TIMING] {section_name}: {elapsed:.2f}s")


# ────────────────────────── Utility Functions ──────────────────────────

def _which(exe: str) -> bool:
    """Check if executable exists on PATH."""
    from shutil import which
    return which(exe) is not None


def _run(cmd: list) -> subprocess.CompletedProcess:
    """Run command and return result."""
    logger.debug(f"RUN: {' '.join(cmd)}")
    try:
        return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed ({' '.join(cmd)}): {e.stderr.strip()}")
        raise


def _read_bed6(path: Path) -> List[dict]:
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


def _prepare_bed6_sorted(path: Path, tmps: list) -> Path:
    """Trim to BED6 and sort."""
    trimmed = path.with_suffix(path.suffix + ".trimmed.tmp")
    sorted_p = path.with_suffix(path.suffix + ".sorted.tmp")
    
    res = _run(["cut", "-f", "1-6", str(path)])
    trimmed.write_text(res.stdout)
    
    res = _run(["bedtools", "sort", "-i", str(trimmed)])
    sorted_p.write_text(res.stdout)
    
    try:
        trimmed.unlink()
    except Exception:
        pass
    tmps.append(sorted_p)
    return sorted_p


def _extract_tss_tts_bed(iso_bed: Path, tmps: list) -> Tuple[Path, Path]:
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


def _run_closest(a: Path, b: Path) -> List[List[str]]:
    """Run bedtools closest."""
    with timed_section(f"bedtools_closest"):
        res = _run(["bedtools", "closest", "-a", str(a), "-b", str(b), "-s", "-d"])
        rows = []
        for line in res.stdout.strip().split('\n'):
            if line and not line.startswith('#'):
                rows.append(line.split('\t'))
        return rows


def _extract_distance_and_peak(rows: List[List[str]], label: str, max_dist: int) -> Dict[str, Tuple[int, str]]:
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


def _extract_signed_distances(rows: List[List[str]]) -> List[int]:
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


def _plot_distance_histogram(
    distances: List[int],
    output_path: Path,
    title: str,
    bin_size: int = 50,
    min_dist: int = -1000,
    max_dist: int = 1000,
) -> bool:
    """
    Create a histogram of signed distances to experimental peaks.

    Args:
        distances: List of signed distances
        output_path: Path to save the plot
        title: Plot title
        bin_size: Size of each bin in bp (default: 50)
        min_dist: Minimum distance for x-axis (default: -1000)
        max_dist: Maximum distance for x-axis (default: 1000)

    Returns:
        True if plot was created successfully, False otherwise
    """
    with timed_section("plot_distance_histogram"):
        if not HAS_MATPLOTLIB:
            logger.warning("matplotlib not available; cannot create distance histogram")
            return False

    if not distances:
        logger.warning(f"No distances to plot for {title}")
        return False

    # Clamp distances into overflow bins at the edges
    # Distances beyond min_dist/max_dist are collected into the outermost bins
    clamped = [max(min_dist, min(max_dist, d)) for d in distances]
    n_clamped_low = sum(1 for d in distances if d < min_dist)
    n_clamped_high = sum(1 for d in distances if d > max_dist)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # Calculate bins
    bins = list(range(min_dist, max_dist + bin_size, bin_size))

    # Create histogram with clamped distances so outliers appear in edge bins
    ax.hist(clamped, bins=bins, edgecolor='black', alpha=0.7, color='steelblue')

    # Add vertical line at 0
    ax.axvline(x=0, color='red', linestyle='--', linewidth=1.5, label='Perfect alignment')

    # Labels and title
    ax.set_xlabel('Distance to Nearest Peak (bp)', fontsize=12)
    ax.set_ylabel('Number of Transcripts', fontsize=12)
    ax.set_title(title, fontsize=14)

    # Add statistics annotation
    mean_dist = statistics.mean(distances)
    median_dist = statistics.median(distances)
    std_dist = statistics.stdev(distances) if len(distances) > 1 else 0

    stats_text = f'n = {len(distances)}\nMean = {mean_dist:.1f} bp\nMedian = {median_dist:.1f} bp\nStd = {std_dist:.1f} bp'
    if n_clamped_low or n_clamped_high:
        stats_text += f'\n<{min_dist}: {n_clamped_low}  >{max_dist}: {n_clamped_high}'
    ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Set x-axis limits
    ax.set_xlim(min_dist, max_dist)

    # Add legend
    ax.legend(loc='upper left')

    # Tight layout
    plt.tight_layout()

    # Save figure
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved distance histogram to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to save histogram to {output_path}: {e}")
        plt.close(fig)
        return False


def _vectorized_overlap_counts(bed_rows: List[dict], peaks_rows: List[dict], window: int) -> Tuple[int, int]:
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


def _safe_f1(p: Optional[float], r: Optional[float]) -> Optional[float]:
    """Calculate F1 score."""
    if p is None or r is None:
        return None
    s = p + r
    return (2 * p * r / s) if s > 0 else 0.0


# ────────────────────────── Read Assignment Functions ──────────────────────────

def _parse_map_file_comprehensive(map_path: Path) -> Tuple[Set[str], dict]:
    """Parse map file once to extract both unique read IDs and per-isoform stats.
    
    Returns:
        Tuple of (unique_read_ids, reads_per_isoform_stats)
    """
    if not map_path.exists() or map_path.stat().st_size == 0:
        empty_stats = {
            "reads_per_isoform_mean": None,
            "reads_per_isoform_median": None,
            "reads_per_isoform_min": None,
            "reads_per_isoform_max": None,
        }
        return set(), empty_stats
    
    uniq = set()
    read_counts = []
    
    with open(map_path) as f:
        for line in f:
            line = line.strip()
            if not line or "\t" not in line:
                continue
            _, rhs = line.split("\t", 1)
            reads = [r.strip() for r in rhs.split(",") if r.strip()]
            read_counts.append(len(reads))
            uniq.update(reads)
    
    if not read_counts:
        stats = {
            "reads_per_isoform_mean": None,
            "reads_per_isoform_median": None,
            "reads_per_isoform_min": None,
            "reads_per_isoform_max": None,
        }
    else:
        stats = {
            "reads_per_isoform_mean": statistics.mean(read_counts),
            "reads_per_isoform_median": statistics.median(read_counts),
            "reads_per_isoform_min": min(read_counts),
            "reads_per_isoform_max": max(read_counts),
        }
    
    return uniq, stats


def _get_assigned_read_ids(map_path: Path) -> Set[str]:
    """Extract set of unique read IDs from isoform.read.map.txt.
    
    Note: For better performance, use _parse_map_file_comprehensive if you also
    need reads_per_isoform statistics.
    """
    if not map_path.exists() or map_path.stat().st_size == 0:
        return set()
    
    uniq = set()
    with open(map_path) as f:
        for line in f:
            line = line.strip()
            if not line or "\t" not in line:
                continue
            _, rhs = line.split("\t", 1)
            for rid in rhs.split(","):
                rid = rid.strip()
                if rid:
                    uniq.add(rid)
    return uniq


def _count_reads_per_isoform(map_path: Path) -> dict:
    """
    Calculate reads per isoform statistics from isoform.read.map.txt.
    
    Returns dict with: mean, median, min, max reads per isoform.
    """
    if not map_path.exists() or map_path.stat().st_size == 0:
        return {
            "reads_per_isoform_mean": None,
            "reads_per_isoform_median": None,
            "reads_per_isoform_min": None,
            "reads_per_isoform_max": None,
        }
    
    read_counts = []
    with open(map_path) as f:
        for line in f:
            line = line.strip()
            if not line or "\t" not in line:
                continue
            _, rhs = line.split("\t", 1)
            reads = [r.strip() for r in rhs.split(",") if r.strip()]
            read_counts.append(len(reads))
    
    if not read_counts:
        return {
            "reads_per_isoform_mean": None,
            "reads_per_isoform_median": None,
            "reads_per_isoform_min": None,
            "reads_per_isoform_max": None,
        }
    
    return {
        "reads_per_isoform_mean": statistics.mean(read_counts),
        "reads_per_isoform_median": statistics.median(read_counts),
        "reads_per_isoform_min": min(read_counts),
        "reads_per_isoform_max": max(read_counts),
    }


def _count_assigned_reads_by_type(bam_path: Path, assigned_read_ids: Set[str]) -> dict:
    """Count primary and supplementary alignments for assigned reads using samtools.

    Uses streaming to avoid loading entire BAM output into memory.
    """
    if not bam_path or not bam_path.exists() or not assigned_read_ids:
        return {
            "assigned_primary": 0,
            "assigned_supplementary": 0,
            "assigned_total": 0,
        }

    if not _which("samtools"):
        logger.warning("samtools not available; cannot count assigned reads by type")
        return {
            "assigned_primary": 0,
            "assigned_supplementary": 0,
            "assigned_total": 0,
        }

    primary_count = 0
    supplementary_count = 0

    try:
        # Stream BAM output line by line to avoid memory issues with large files
        logger.debug(f"RUN (streaming): samtools view -F 4 {bam_path}")
        proc = subprocess.Popen(
            ["samtools", "view", "-F", "4", str(bam_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1  # Line buffered
        )

        # Process output line by line (streaming)
        for line in proc.stdout:
            line = line.rstrip('\n')
            if not line:
                continue

            # Only parse the first two fields (read name and flag)
            # This is faster than splitting the entire line
            tab1 = line.find('\t')
            if tab1 == -1:
                continue
            tab2 = line.find('\t', tab1 + 1)
            if tab2 == -1:
                tab2 = len(line)

            read_name = line[:tab1]
            if read_name not in assigned_read_ids:
                continue

            try:
                flag = int(line[tab1+1:tab2])
            except ValueError:
                continue

            is_supplementary = bool(flag & 0x800)
            is_secondary = bool(flag & 0x100)

            if is_supplementary:
                supplementary_count += 1
            elif not is_secondary:
                primary_count += 1

        # Wait for process to finish and check return code
        proc.wait()
        if proc.returncode != 0:
            stderr = proc.stderr.read()
            logger.warning(f"samtools view returned non-zero exit code: {stderr}")

    except Exception as e:
        logger.warning(f"Error counting assigned reads by type from BAM: {e}")
        return {
            "assigned_primary": 0,
            "assigned_supplementary": 0,
            "assigned_total": 0,
        }

    return {
        "assigned_primary": primary_count,
        "assigned_supplementary": supplementary_count,
        "assigned_total": primary_count + supplementary_count,
    }


def _count_total_alignments_bam(bam: Path) -> dict:
    """Count total primary and supplementary alignments in BAM."""
    if not _which("samtools"):
        logger.warning("samtools not available; cannot count total alignments")
        return {
            "total_primary": 0,
            "total_supplementary": 0,
            "total_alignments": 0,
        }
    
    if not bam.exists() or bam.stat().st_size == 0:
        return {
            "total_primary": 0,
            "total_supplementary": 0,
            "total_alignments": 0,
        }
    
    # Count primary alignments (exclude unmapped=0x4, secondary=0x100, supplementary=0x800)
    try:
        res = _run(["samtools", "view", "-c", "-F", "2308", str(bam)])
        primary = int(res.stdout.strip() or "0")
    except Exception:
        primary = 0
    
    # Count supplementary alignments (include supplementary=0x800, exclude unmapped=0x4, secondary=0x100)
    try:
        res = _run(["samtools", "view", "-c", "-f", "2048", "-F", "260", str(bam)])
        supplementary = int(res.stdout.strip() or "0")
    except Exception:
        supplementary = 0
    
    return {
        "total_primary": primary,
        "total_supplementary": supplementary,
        "total_alignments": primary + supplementary,
    }


def count_lines(path: Path) -> int:
    """Count lines in a file."""
    if not path.exists():
        return 0
    try:
        with open(path) as f:
            return sum(1 for _ in f)
    except Exception:
        return 0


# ────────────────────────── Read-End Entropy Functions ──────────────────────────

def _parse_isoform_ends(iso_bed: Path) -> Dict[str, dict]:
    """Parse isoform BED to get model TSS/TTS positions per isoform.

    Returns {isoform_id: {'chrom': str, 'start': int, 'end': int, 'strand': str,
                          'tss': int, 'tts': int}}
    """
    with timed_section("parse_isoform_ends"):
        isoforms = {}
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
                strand = row[5]
                if strand == '+':
                    tss, tts = start, end
                else:
                    tss, tts = end, start
                isoforms[name] = {
                    'chrom': chrom, 'start': start, 'end': end,
                    'strand': strand, 'tss': tss, 'tts': tts,
                }
        return isoforms


def _parse_read_map(map_path: Path) -> Dict[str, List[str]]:
    """Parse isoform.read.map.txt to get {isoform_id: [read_id, ...]}."""
    with timed_section("parse_read_map"):
        iso_to_reads = {}
        if not map_path.exists() or map_path.stat().st_size == 0:
            return iso_to_reads
        with open(map_path) as f:
            for line in f:
                line = line.strip()
                if not line or '\t' not in line:
                    continue
                iso_id, rhs = line.split('\t', 1)
                reads = [r.strip() for r in rhs.split(',') if r.strip()]
                if reads:
                    iso_to_reads[iso_id] = reads
        return iso_to_reads


def _parse_reads_bed_ends(reads_bed: Path) -> Dict[str, dict]:
    """Parse reads BED12 to get 5'/3' end positions per read.

    Returns {read_name: {'chrom': str, 'start': int, 'end': int, 'strand': str,
                         'tss': int, 'tts': int}}
    For reads with duplicate names (supplementary alignments), keeps only the first
    (primary) occurrence.
    """
    reads = {}
    with open(reads_bed) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            if len(row) < 6:
                continue
            name = row[3]
            if name in reads:
                continue  # keep first (primary) alignment
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            strand = row[5]
            if strand == '+':
                tss, tts = start, end
            else:
                tss, tts = end, start
            reads[name] = {
                'chrom': chrom, 'start': start, 'end': end,
                'strand': strand, 'tss': tss, 'tts': tts,
            }
    return reads


def _shannon_entropy(values: List[int], bin_size: int = 10) -> float:
    """Compute Shannon entropy (bits) of a list of integer values, binned."""
    if not values:
        return 0.0
    # Bin the values
    binned = [v // bin_size for v in values]
    counts = defaultdict(int)
    for b in binned:
        counts[b] += 1
    n = len(binned)
    entropy = 0.0
    for c in counts.values():
        p = c / n
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy


def _compute_read_end_entropy(
    iso_bed: Path,
    map_path: Path,
    reads_bed: Path,
    bin_size: int = 10,
) -> dict:
    """
    Compute per-isoform entropy of read-end positions relative to the isoform model.

    For each isoform, collects the 5' and 3' end positions of all assigned reads,
    computes the offset from the isoform model's TSS/TTS, and measures Shannon entropy
    of these offsets (binned). High entropy = reads are spread out; low = tightly clustered.

    Returns dict with:
        - tss_entropy_per_isoform: list of (isoform_id, entropy, n_reads)
        - tts_entropy_per_isoform: list of (isoform_id, entropy, n_reads)
        - all_tss_offsets: list of all read TSS offsets from model TSS
        - all_tts_offsets: list of all read TTS offsets from model TTS
        - summary metrics (mean/median entropy)
    """
    isoforms = _parse_isoform_ends(iso_bed)
    iso_to_reads = _parse_read_map(map_path)
    read_ends = _parse_reads_bed_ends(reads_bed)

    tss_entropy_list = []  # (iso_id, entropy, n_reads)
    tts_entropy_list = []
    all_tss_offsets = []
    all_tts_offsets = []

    for iso_id, read_ids in iso_to_reads.items():
        if iso_id not in isoforms:
            continue
        iso = isoforms[iso_id]
        model_tss = iso['tss']
        model_tts = iso['tts']
        strand = iso['strand']

        tss_offsets = []
        tts_offsets = []

        for rid in read_ids:
            if rid not in read_ends:
                continue
            r = read_ends[rid]
            # Signed offset: positive = downstream of model end
            if strand == '+':
                tss_offsets.append(r['tss'] - model_tss)
                tts_offsets.append(r['tts'] - model_tts)
            else:
                # For minus strand, flip sign so positive = downstream (towards 5')
                tss_offsets.append(model_tss - r['tss'])
                tts_offsets.append(model_tts - r['tts'])

        if tss_offsets:
            tss_ent = _shannon_entropy(tss_offsets, bin_size=bin_size)
            tss_entropy_list.append((iso_id, tss_ent, len(tss_offsets)))
            all_tss_offsets.extend(tss_offsets)

        if tts_offsets:
            tts_ent = _shannon_entropy(tts_offsets, bin_size=bin_size)
            tts_entropy_list.append((iso_id, tts_ent, len(tts_offsets)))
            all_tts_offsets.extend(tts_offsets)

    # Summary metrics
    tss_entropies = [e for _, e, _ in tss_entropy_list]
    tts_entropies = [e for _, e, _ in tts_entropy_list]

    result = {
        "tss_entropy_per_isoform": tss_entropy_list,
        "tts_entropy_per_isoform": tts_entropy_list,
        "all_tss_offsets": all_tss_offsets,
        "all_tts_offsets": all_tts_offsets,
        "tss_entropy_mean": statistics.mean(tss_entropies) if tss_entropies else None,
        "tss_entropy_median": statistics.median(tss_entropies) if tss_entropies else None,
        "tts_entropy_mean": statistics.mean(tts_entropies) if tts_entropies else None,
        "tts_entropy_median": statistics.median(tts_entropies) if tts_entropies else None,
    }

    logger.info(f"Read-end entropy: {len(tss_entropy_list)} isoforms with TSS data, "
                f"{len(tts_entropy_list)} with TTS data")

    return result


def _plot_read_end_entropy(
    entropy_data: dict,
    plot_output_dir: Path,
    plot_prefix: str,
) -> None:
    """Generate read-end entropy plots:
    1. Aggregate histogram of read-to-model offsets (TSS and TTS)
    2. Distribution of per-isoform entropies (TSS and TTS)
    """
    if not HAS_MATPLOTLIB:
        return

    plot_output_dir.mkdir(parents=True, exist_ok=True)

    # Plot 1a: Aggregate TSS offset histogram
    tss_offsets = entropy_data["all_tss_offsets"]
    if tss_offsets:
        _plot_distance_histogram(
            distances=tss_offsets,
            output_path=plot_output_dir / f"{plot_prefix}_read_tss_offset_histogram.png",
            title=f"Read 5' End Offset from Isoform Model TSS\n{plot_prefix}",
        )

    # Plot 1b: Aggregate TTS offset histogram
    tts_offsets = entropy_data["all_tts_offsets"]
    if tts_offsets:
        _plot_distance_histogram(
            distances=tts_offsets,
            output_path=plot_output_dir / f"{plot_prefix}_read_tts_offset_histogram.png",
            title=f"Read 3' End Offset from Isoform Model TTS\n{plot_prefix}",
        )

    # Plot 2a: Per-isoform TSS entropy distribution
    tss_entropies = [e for _, e, _ in entropy_data["tss_entropy_per_isoform"]]
    if tss_entropies:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(tss_entropies, bins=50, edgecolor='black', alpha=0.7, color='coral')
        ax.set_xlabel("Shannon Entropy (bits, 10bp bins)", fontsize=12)
        ax.set_ylabel("Number of Isoforms", fontsize=12)
        ax.set_title(f"Per-Isoform TSS Read-End Entropy (5' End)\n{plot_prefix}", fontsize=14)
        mean_e = statistics.mean(tss_entropies)
        median_e = statistics.median(tss_entropies)
        stats_text = f'n = {len(tss_entropies)}\nMean = {mean_e:.2f} bits\nMedian = {median_e:.2f} bits'
        ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.tight_layout()
        plt.savefig(plot_output_dir / f"{plot_prefix}_tss_entropy_distribution.png",
                    dpi=150, bbox_inches='tight')
        plt.close(fig)

    # Plot 2b: Per-isoform TTS entropy distribution
    tts_entropies = [e for _, e, _ in entropy_data["tts_entropy_per_isoform"]]
    if tts_entropies:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(tts_entropies, bins=50, edgecolor='black', alpha=0.7, color='mediumpurple')
        ax.set_xlabel("Shannon Entropy (bits, 10bp bins)", fontsize=12)
        ax.set_ylabel("Number of Isoforms", fontsize=12)
        ax.set_title(f"Per-Isoform TTS Read-End Entropy (3' End)\n{plot_prefix}", fontsize=14)
        mean_e = statistics.mean(tts_entropies)
        median_e = statistics.median(tts_entropies)
        stats_text = f'n = {len(tts_entropies)}\nMean = {mean_e:.2f} bits\nMedian = {median_e:.2f} bits'
        ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.tight_layout()
        plt.savefig(plot_output_dir / f"{plot_prefix}_tts_entropy_distribution.png",
                    dpi=150, bbox_inches='tight')
        plt.close(fig)


# ────────────────────────── Peak Recoverability Functions ──────────────────────────

def _find_recoverable_peaks(
    peaks_path: Path,
    read_ends: List[dict],
    window: int,
) -> Dict[str, int]:
    """Determine which peaks have at least one read end within the window.

    Args:
        peaks_path: BED file of experimental peaks (CAGE or QuantSeq)
        read_ends: List of dicts with 'Chrom', 'Start', 'End', 'Strand' for each read end
        window: Distance window in bp

    Returns:
        Dict mapping peak IDs (chrom_start_end) to the number of supporting reads
    """
    peaks_rows = _read_bed6(peaks_path)
    if not peaks_rows or not read_ends:
        return {}

    # Group read ends by chrom/strand
    reads_by_cs = defaultdict(list)
    for r in read_ends:
        key = (r['Chrom'], r['Strand'])
        reads_by_cs[key].append(r)

    recoverable: Dict[str, int] = {}
    for peak in peaks_rows:
        peak_id = f"{peak['Chrom']}_{peak['Start']}_{peak['End']}"
        lo = peak['Start'] - window
        hi = peak['End'] + window

        count = 0
        strands_to_check = [peak['Strand']]
        if peak['Strand'] == '.':
            strands_to_check = ['+', '-']

        for strand in strands_to_check:
            key = (peak['Chrom'], strand)
            for r in reads_by_cs.get(key, []):
                if r['End'] >= lo and r['Start'] <= hi:
                    count += 1

        if count > 0:
            recoverable[peak_id] = count

    logger.debug(f"Recoverable peaks: {len(recoverable)}/{len(peaks_rows)}")
    return recoverable


def _classify_isoform_recoverability(
    closest_rows: List[List[str]],
    recoverable_peaks: Dict[str, int],
) -> List[bool]:
    """For each row from bedtools closest, classify whether the nearest peak is recoverable.

    Returns list of booleans parallel to the signed distances extracted from the same rows.
    """
    classifications = []
    for row in closest_rows:
        if not row or len(row) < 13:
            continue

        # Check if no peak was found
        try:
            abs_dist = int(row[-1])
            if abs_dist == -1:
                continue
        except (ValueError, IndexError):
            continue

        # Build peak ID from columns 6-8 (peak BED6)
        try:
            peak_id = f"{row[6]}_{row[7]}_{row[8]}"
        except IndexError:
            classifications.append(False)
            continue

        classifications.append(peak_id in recoverable_peaks)

    return classifications


def _extract_read_end_positions(reads_bed: Path, end_type: str) -> List[dict]:
    """Extract single-position BED-like dicts for read TSS or TTS positions.

    Args:
        reads_bed: Reads BED12 file
        end_type: 'tss' for 5' ends, 'tts' for 3' ends

    Returns:
        List of dicts with 'Chrom', 'Start', 'End', 'Strand'
    """
    positions = []
    with open(reads_bed) as f:
        reader = csv.reader(f, delimiter='\t')
        seen = set()
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            if len(row) < 6:
                continue
            name = row[3]
            if name in seen:
                continue
            seen.add(name)

            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            strand = row[5]

            if end_type == 'tss':
                pos = start if strand == '+' else end
            else:  # tts
                pos = end if strand == '+' else start

            positions.append({
                'Chrom': chrom,
                'Start': pos,
                'End': pos + 1,
                'Strand': strand,
            })
    return positions


def _write_recoverable_peaks_bed(
    peaks_path: Path,
    recoverable_ids: Dict[str, int],
    output_path: Path,
) -> int:
    """Write recoverable peak coordinates to a BED file.

    Returns number of peaks written.
    """
    peaks_rows = _read_bed6(peaks_path)
    count = 0
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        for peak in peaks_rows:
            peak_id = f"{peak['Chrom']}_{peak['Start']}_{peak['End']}"
            if peak_id in recoverable_ids:
                writer.writerow([peak['Chrom'], peak['Start'], peak['End'],
                                 peak_id, '.', peak['Strand']])
                count += 1
    logger.info(f"Wrote {count} recoverable peaks to {output_path}")
    return count


def _plot_distance_histogram_colored(
    distances: List[int],
    recoverable_mask: List[bool],
    output_path: Path,
    title: str,
    bin_size: int = 50,
    min_dist: int = -1000,
    max_dist: int = 1000,
) -> bool:
    """Create a histogram with bars colored by peak recoverability.

    Recoverable = at least one long read end is within the window of the peak.
    """
    if not HAS_MATPLOTLIB:
        return False
    if not distances:
        return False
    if len(distances) != len(recoverable_mask):
        logger.warning(f"Distance/mask length mismatch: {len(distances)} vs {len(recoverable_mask)}, "
                       f"falling back to uncolored plot")
        return _plot_distance_histogram(distances, output_path, title, bin_size, min_dist, max_dist)

    # Split distances by recoverability
    dists_recoverable = [d for d, r in zip(distances, recoverable_mask) if r]
    dists_unrecoverable = [d for d, r in zip(distances, recoverable_mask) if not r]

    # Clamp
    clamp = lambda d: max(min_dist, min(max_dist, d))
    clamped_rec = [clamp(d) for d in dists_recoverable]
    clamped_unrec = [clamp(d) for d in dists_unrecoverable]

    n_clamped_low = sum(1 for d in distances if d < min_dist)
    n_clamped_high = sum(1 for d in distances if d > max_dist)

    fig, ax = plt.subplots(figsize=(10, 6))
    bins = list(range(min_dist, max_dist + bin_size, bin_size))

    # Stacked histogram: recoverable on bottom, unrecoverable on top
    ax.hist([clamped_rec, clamped_unrec], bins=bins, stacked=True,
            edgecolor='black', alpha=0.7,
            color=['steelblue', 'lightcoral'],
            label=[f'Read-supported peak ({len(dists_recoverable)})',
                   f'No read support ({len(dists_unrecoverable)})'])

    ax.axvline(x=0, color='red', linestyle='--', linewidth=1.5, label='Perfect alignment')
    ax.set_xlabel('Distance to Nearest Peak (bp)', fontsize=12)
    ax.set_ylabel('Number of Transcripts', fontsize=12)
    ax.set_title(title, fontsize=14)

    # Stats on full distribution
    mean_dist = statistics.mean(distances)
    median_dist = statistics.median(distances)
    std_dist = statistics.stdev(distances) if len(distances) > 1 else 0
    stats_text = (f'n = {len(distances)}\nMean = {mean_dist:.1f} bp\n'
                  f'Median = {median_dist:.1f} bp\nStd = {std_dist:.1f} bp')
    if n_clamped_low or n_clamped_high:
        stats_text += f'\n<{min_dist}: {n_clamped_low}  >{max_dist}: {n_clamped_high}'
    ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.set_xlim(min_dist, max_dist)
    ax.legend(loc='upper left')
    plt.tight_layout()

    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return True
    except Exception as e:
        logger.error(f"Failed to save colored histogram: {e}")
        plt.close(fig)
        return False


def _find_captured_peaks(
    closest_rows: List[List[str]],
    window: int,
) -> Set[str]:
    """Find peak IDs that have at least one isoform end within the window.

    These are peaks already "captured" by the isoform model.
    """
    captured = set()
    for row in closest_rows:
        if not row or len(row) < 13:
            continue
        try:
            abs_dist = int(row[-1])
            if abs_dist == -1:
                continue
        except (ValueError, IndexError):
            continue
        if abs_dist <= window:
            try:
                peak_id = f"{row[6]}_{row[7]}_{row[8]}"
                captured.add(peak_id)
            except IndexError:
                continue
    return captured


def _plot_read_support_distribution(
    recoverable_counts: Dict[str, int],
    output_path: Path,
    title: str,
) -> bool:
    """Plot read support counts for recoverable peaks.

    Displays:
      - 0–100: unit-width bins (one bar per integer)
      - 100+: a single overflow bar whose height equals the number of peaks with support >100
    """
    if not HAS_MATPLOTLIB:
        return False
    if not recoverable_counts:
        return False

    counts = list(recoverable_counts.values())

    fig, ax = plt.subplots(figsize=(12, 6))  # wider helps readability

    # Build frequency table for 0..100 and overflow (100+)
    max_bin = 100
    freq = [0] * (max_bin + 1)  # index i = count i
    overflow = 0

    for c in counts:
        if c < 0:
            continue
        if c <= max_bin:
            freq[c] += 1
        else:
            overflow += 1

    # X positions: 0..100, plus overflow at 101 labeled "100+"
    x_main = list(range(0, max_bin + 1))
    y_main = freq

    overflow_x = max_bin + 1  # 101
    ax.bar(x_main, y_main, width=1.0, edgecolor='black', alpha=0.7)
    ax.bar([overflow_x], [overflow], width=1.0, edgecolor='black', alpha=0.7)

    ax.set_xlabel('Number of Supporting Long Reads', fontsize=12)
    ax.set_ylabel('Number of Peaks', fontsize=12)
    ax.set_title(title, fontsize=14)

    # Ticks: keep sparse for readability, and label overflow as 100+
    tick_positions = [0, 1, 2, 3, 4, 5, 10, 20, 30, 50, 75, 100, overflow_x]
    tick_labels = [str(t) for t in tick_positions[:-1]] + ['100+']
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)

    # Limits so the overflow bar is fully visible
    ax.set_xlim(-0.5, overflow_x + 1.5)

    # Summary stats
    mean_c = statistics.mean(counts)
    median_c = statistics.median(counts)
    std_c = statistics.stdev(counts) if len(counts) > 1 else 0
    stats_text = (
        f'n = {len(counts)} peaks\n'
        f'Mean = {mean_c:.1f} reads\n'
        f'Median = {median_c:.1f} reads\n'
        f'Std = {std_c:.1f}\n'
        f'>100: {overflow}'
    )
    ax.text(
        0.98, 0.98, stats_text,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment='top',
        horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
    )

    # Optional: if one bin dominates and hides structure, uncomment:
    # ax.set_yscale('log')

    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return True
    except Exception as e:
        logger.error(f"Failed to save read support histogram: {e}")
        plt.close(fig)
        return False


# ────────────────────────── Missed Peak Analysis Functions ──────────────────────────

def _get_reads_to_isoforms(iso_to_reads: Dict[str, List[str]]) -> Dict[str, str]:
    """Invert isoform->reads map to get read->isoform mapping."""
    read_to_iso = {}
    for iso_id, read_ids in iso_to_reads.items():
        for rid in read_ids:
            if rid not in read_to_iso:
                read_to_iso[rid] = iso_id
    return read_to_iso


def _classify_missed_peak_reads(
    peak_id: str,
    peak_info: dict,
    reads_supporting_peak: List[str],
    iso_to_reads: Dict[str, List[str]],
    isoforms: Dict[str, dict],
    read_ends: Dict[str, dict],
    end_type: str,  # 'tss' (CAGE/5') or 'tts' (QuantSeq/3')
    window: int = 50,
    nearby_threshold: int = 200,
) -> dict:
    """
    For reads that support a missed peak, classify where they actually got assigned.

    Categories:
    - unassigned: read wasn't assigned to any isoform
    - assigned_nearby: read assigned to isoform with the relevant end (TSS or TTS)
                      within nearby_threshold bp of the peak
    - assigned_distant: read assigned to isoform with relevant end far from peak
    - assigned_wrong_strand: read assigned to isoform on different strand

    Args:
        peak_id: Peak identifier (chrom_start_end)
        peak_info: Dict with 'Chrom', 'Start', 'End', 'Strand'
        reads_supporting_peak: List of read IDs that have ends near this peak
        iso_to_reads: Dict mapping isoform_id -> list of read_ids
        isoforms: Dict mapping isoform_id -> {'chrom','start','end','strand','tss','tts'}
        read_ends: Dict mapping read_id -> {'chrom','start','end','strand','tss','tts'} (currently unused here)
        end_type: 'tss' to classify relative to isoform TSS (CAGE) or 'tts' for isoform TTS (QuantSeq)
        window: Window used for peak detection (not used directly in this classifier)
        nearby_threshold: Distance threshold for "nearby" classification

    Returns:
        Dict with counts for each category and example read IDs
    """
    if end_type not in ("tss", "tts"):
        raise ValueError(f"end_type must be 'tss' or 'tts', got: {end_type}")

    read_to_iso = _get_reads_to_isoforms(iso_to_reads)

    peak_pos = (peak_info["Start"] + peak_info["End"]) // 2
    peak_strand = peak_info["Strand"]

    categories = {
        "unassigned": [],
        "assigned_nearby": [],
        "assigned_distant": [],
        "assigned_wrong_strand": [],
    }

    for read_id in reads_supporting_peak:
        iso_id = read_to_iso.get(read_id)
        if not iso_id:
            categories["unassigned"].append(read_id)
            continue

        iso = isoforms.get(iso_id)
        if not iso:
            categories["unassigned"].append(read_id)
            continue

        # Strand mismatch (ignore if peak strand is '.')
        if peak_strand != "." and iso.get("strand") != peak_strand:
            categories["assigned_wrong_strand"].append(read_id)
            continue

        # Compare to the relevant isoform endpoint
        iso_end_pos = iso.get(end_type)
        if iso_end_pos is None:
            categories["unassigned"].append(read_id)
            continue

        dist = abs(int(iso_end_pos) - int(peak_pos))
        if dist <= nearby_threshold:
            categories["assigned_nearby"].append(read_id)
        else:
            categories["assigned_distant"].append(read_id)

    return {
        "peak_id": peak_id,
        "end_type": end_type,
        "total_supporting_reads": len(reads_supporting_peak),
        "unassigned_count": len(categories["unassigned"]),
        "assigned_nearby_count": len(categories["assigned_nearby"]),
        "assigned_distant_count": len(categories["assigned_distant"]),
        "assigned_wrong_strand_count": len(categories["assigned_wrong_strand"]),
        "unassigned_reads": categories["unassigned"][:5],
        "assigned_nearby_reads": categories["assigned_nearby"][:5],
        "assigned_distant_reads": categories["assigned_distant"][:5],
        "assigned_wrong_strand_reads": categories["assigned_wrong_strand"][:5],
    }



def _characterize_truncation_pattern(
    read_starts: List[int],
    peak_pos: int,
    strand: str,
    sharp_threshold: int = 20,
    trailing_ratio_threshold: float = 0.3,
) -> Tuple[str, dict]:
    """
    Characterize the shape of read start distribution relative to peak.

    For 5' ends, examines whether reads cluster at the peak or trail off,
    which can indicate RT truncation patterns.

    Patterns:
    - "sharp": reads cluster tightly at peak (good signal, low truncation)
    - "trailing": reads trail off in 5' direction (truncation gradient)
    - "bimodal": two distinct clusters (possible alternative TSS)
    - "sparse": too few reads to characterize

    Args:
        read_starts: List of read start positions (5' ends)
        peak_pos: Position of the experimental peak
        strand: Strand of the peak ('+' or '-')
        sharp_threshold: Max distance from peak for "sharp" classification
        trailing_ratio_threshold: Ratio of upstream to downstream reads for "trailing"

    Returns:
        Tuple of (pattern_name, details_dict)
    """
    if len(read_starts) < 5:
        return "sparse", {"n_reads": len(read_starts)}

    # Calculate signed distances (positive = downstream of peak in transcription direction)
    if strand == '+':
        distances = [pos - peak_pos for pos in read_starts]
    else:
        distances = [peak_pos - pos for pos in read_starts]

    # Upstream (5' of peak, negative values) vs downstream (3' of peak, positive values)
    upstream = [d for d in distances if d < 0]
    downstream = [d for d in distances if d >= 0]
    at_peak = [d for d in distances if abs(d) <= sharp_threshold]

    details = {
        "n_reads": len(read_starts),
        "n_upstream": len(upstream),
        "n_downstream": len(downstream),
        "n_at_peak": len(at_peak),
        "mean_distance": statistics.mean(distances) if distances else 0,
        "median_distance": statistics.median(distances) if distances else 0,
        "std_distance": statistics.stdev(distances) if len(distances) > 1 else 0,
    }

    # Check for sharp clustering at peak
    at_peak_ratio = len(at_peak) / len(read_starts)
    if at_peak_ratio >= 0.7:
        return "sharp", details

    # Check for trailing pattern (many reads upstream/5' of peak = truncation)
    if len(upstream) > 0 and len(downstream) > 0:
        upstream_ratio = len(upstream) / len(read_starts)
        if upstream_ratio >= trailing_ratio_threshold:
            # Further check: is there a gradient (more reads closer to peak)?
            if upstream:
                upstream_dists = sorted([abs(d) for d in upstream])
                # Check if closer reads are more numerous
                close_upstream = sum(1 for d in upstream if abs(d) <= 100)
                far_upstream = sum(1 for d in upstream if abs(d) > 100)
                if close_upstream > far_upstream:
                    details["trailing_gradient"] = True
            return "trailing", details

    # Check for bimodal distribution
    if len(distances) >= 10:
        sorted_dists = sorted(distances)
        # Simple bimodal check: look for gap in middle
        mid_start = len(sorted_dists) // 3
        mid_end = 2 * len(sorted_dists) // 3
        mid_range = sorted_dists[mid_end] - sorted_dists[mid_start]
        full_range = sorted_dists[-1] - sorted_dists[0]
        if full_range > 0:
            # If middle third spans less than 20% of range, might be bimodal
            if mid_range / full_range < 0.2 and full_range > 100:
                return "bimodal", details

    # Default: dispersed pattern
    return "dispersed", details


def _extract_sequences_batch(
    genome_path: Path,
    regions: List[Tuple[str, int, str, int, int]],
) -> Dict[int, str]:
    """
    Extract multiple sequences in a single samtools call (50-100x faster than individual calls).

    Args:
        genome_path: Path to genome FASTA file
        regions: List of (chrom, pos, strand, upstream, downstream) tuples

    Returns:
        Dict mapping index -> sequence string (strand-aware)
    """
    if not regions or not _which("samtools"):
        return {}

    # Create temp file with all regions
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.regions') as f:
        regions_file = f.name
        for idx, (chrom, pos, strand, upstream, downstream) in enumerate(regions):
            # Calculate genomic coordinates
            if strand == '+':
                start = max(0, pos - upstream)
                end = pos + downstream
            else:
                start = max(0, pos - downstream)
                end = pos + upstream
            
            # Write region (1-based for samtools)
            f.write(f"{chrom}:{start + 1}-{end}\n")

    try:
        # Single samtools call for all regions
        res = _run(["samtools", "faidx", str(genome_path), "-r", regions_file])
        
        # Parse sequences from FASTA output
        sequences = {}
        current_seq = []
        current_idx = -1
        
        for line in res.stdout.strip().split('\n'):
            if line.startswith('>'):
                # Save previous sequence
                if current_idx >= 0 and current_seq:
                    seq = ''.join(current_seq).upper()
                    _, _, strand, _, _ = regions[current_idx]
                    
                    # Reverse complement for minus strand
                    if strand == '-':
                        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
                        seq = ''.join(complement.get(b, 'N') for b in reversed(seq))
                    
                    sequences[current_idx] = seq
                
                # Start new sequence
                current_idx += 1
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_idx >= 0 and current_seq:
            seq = ''.join(current_seq).upper()
            _, _, strand, _, _ = regions[current_idx]
            
            if strand == '-':
                complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
                seq = ''.join(complement.get(b, 'N') for b in reversed(seq))
            
            sequences[current_idx] = seq
        
        return sequences
        
    except Exception as e:
        logger.warning(f"Batch sequence extraction failed: {e}")
        return {}
    finally:
        # Clean up temp file
        try:
            os.unlink(regions_file)
        except Exception:
            pass


def _extract_sequence_context(
    genome_path: Path,
    chrom: str,
    pos: int,
    strand: str,
    upstream: int = 50,
    downstream: int = 50,
) -> Optional[str]:
    """
    Extract sequence context around a position using samtools faidx.
    
    Note: For multiple sequences, use _extract_sequences_batch() for 50-100x speedup.

    Args:
        genome_path: Path to genome FASTA file
        chrom: Chromosome name
        pos: Position (0-based)
        strand: Strand ('+' or '-')
        upstream: bp upstream (5' direction on strand)
        downstream: bp downstream (3' direction on strand)

    Returns:
        Sequence string (strand-aware), or None if extraction fails
    """
    # Use batch extraction for single sequence (simpler than duplicating logic)
    result = _extract_sequences_batch(genome_path, [(chrom, pos, strand, upstream, downstream)])
    return result.get(0)


def _compute_position_frequency_matrix(
    sequences: List[str],
) -> Optional[List[Dict[str, float]]]:
    """
    Compute position frequency matrix from aligned sequences.

    Args:
        sequences: List of equal-length sequences

    Returns:
        List of dicts, one per position, with nucleotide frequencies
    """
    if not sequences:
        return None

    seq_len = len(sequences[0])
    if not all(len(s) == seq_len for s in sequences):
        logger.warning("Sequences have unequal lengths; cannot compute PFM")
        return None

    pfm = []
    for i in range(seq_len):
        counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
        for seq in sequences:
            base = seq[i] if i < len(seq) else 'N'
            counts[base] = counts.get(base, 0) + 1

        total = sum(counts.values())
        freqs = {base: count / total for base, count in counts.items()}
        pfm.append(freqs)

    return pfm


def _compute_information_content(pfm: List[Dict[str, float]]) -> List[float]:
    """Compute information content (bits) per position."""
    ic = []
    for pos_freqs in pfm:
        # Shannon entropy
        entropy = 0
        for base in ['A', 'T', 'G', 'C']:
            p = pos_freqs.get(base, 0)
            if p > 0:
                entropy -= p * math.log2(p)
        # Information content = max entropy (2 bits) - observed entropy
        ic.append(2.0 - entropy)
    return ic


def _analyze_motifs_at_ends(
    read_ends: Dict[str, dict],
    genome_path: Path,
    end_type: str,
    upstream: int = 50,
    downstream: int = 50,
    max_sequences: int = 5000,
) -> dict:
    """
    Analyze sequence motifs around read transcript ends using batch extraction.

    Args:
        read_ends: Dict mapping read_id -> {'chrom', 'start', 'end', 'strand', 'tss', 'tts'}
        genome_path: Path to genome FASTA file
        end_type: 'tss' for 5' ends, 'tts' for 3' ends
        upstream: bp upstream of end to extract
        downstream: bp downstream of end to extract
        max_sequences: Maximum sequences to analyze (for performance)

    Returns:
        Dict with PFM, information content, and summary statistics
    """
    if not genome_path or not genome_path.exists():
        logger.warning("Genome file not available for motif analysis")
        return {}

    # Sample reads if too many
    read_list = list(read_ends.items())
    if len(read_list) > max_sequences:
        import random
        random.seed(42)
        read_list = random.sample(read_list, max_sequences)

    # Build regions list for batch extraction
    regions = []
    valid_reads = []
    for read_id, read_info in read_list:
        # Note: _parse_reads_bed_ends returns lowercase keys ('chrom', 'strand', 'tss', 'tts')
        chrom = read_info['chrom']
        strand = read_info['strand']
        pos = read_info.get(end_type)

        if pos is not None:
            regions.append((chrom, pos, strand, upstream, downstream))
            valid_reads.append((read_id, chrom, pos, strand))

    if not regions:
        logger.warning(f"No valid positions for {end_type} motif analysis")
        return {}

    # Batch extract all sequences (50-100x faster than loop)
    logger.debug(f"Extracting {len(regions)} sequences for {end_type} motif analysis (batch mode)")
    sequence_dict = _extract_sequences_batch(genome_path, regions)

    # Filter to expected length
    sequences = []
    positions_used = []
    expected_len = upstream + downstream
    for idx, (read_id, chrom, pos, strand) in enumerate(valid_reads):
        seq = sequence_dict.get(idx)
        if seq and len(seq) == expected_len:
            sequences.append(seq)
            positions_used.append((chrom, pos, strand))

    if not sequences:
        logger.warning(f"No sequences extracted for {end_type} motif analysis")
        return {}

    # Compute position frequency matrix
    pfm = _compute_position_frequency_matrix(sequences)
    if not pfm:
        return {}

    ic = _compute_information_content(pfm)

    # Find positions with high information content
    high_ic_positions = [(i - upstream, ic[i]) for i in range(len(ic)) if ic[i] > 0.5]

    result = {
        "n_sequences": len(sequences),
        "sequence_length": len(sequences[0]) if sequences else 0,
        "pfm": pfm,
        "information_content": ic,
        "mean_ic": statistics.mean(ic) if ic else 0,
        "high_ic_positions": high_ic_positions,
    }

    return result


def _plot_sequence_logo(
    pfm: List[Dict[str, float]],
    output_path: Path,
    title: str,
    upstream: int = 50,
) -> bool:
    """
    Traditional sequence logo:
      - total stack height at each position = information content (bits)
      - each letter height = p(base) * information_content
      - letters are drawn as scaled glyphs (TextPath) so height truly matches bits
    """
    if not HAS_MATPLOTLIB or not pfm:
        return False

    try:
        import matplotlib.pyplot as plt
        from matplotlib.textpath import TextPath
        from matplotlib.patches import PathPatch
        from matplotlib.font_manager import FontProperties
        from matplotlib.transforms import Affine2D
    except Exception:
        return False

    # Information content per position (bits), computed from A/T/G/C frequencies
    ic = _compute_information_content(pfm)
    if not ic:
        return False

    # X positions centered around 0 at the end site
    positions = list(range(-upstream, len(pfm) - upstream))

    # Classic DNA colors (feel free to change)
    colors = {'A': '#00CC00', 'C': '#0000FF', 'G': '#FFB300', 'T': '#FF0000'}
    bases = ['A', 'C', 'G', 'T']

    fig, ax = plt.subplots(figsize=(16, 4))

    # Use a bold monospace-ish font for a "logo" look
    fp = FontProperties(family="DejaVu Sans Mono", weight="bold")
    # TextPath is defined in "points"; we normalize each glyph to height=1 in data units
    glyph_cache = {b: TextPath((0, 0), b, size=1, prop=fp) for b in bases}

    # Helper: add one letter as a true scaled glyph in data coordinates
    def _add_letter(letter: str, x: float, y: float, height_bits: float) -> None:
        if height_bits <= 0:
            return

        tp = glyph_cache[letter]
        bb = tp.get_extents()
        glyph_h = bb.height if bb.height > 0 else 1.0
        glyph_w = bb.width if bb.width > 0 else 1.0

        # Scale so glyph height == height_bits in data units.
        sy = height_bits / glyph_h

        # Make letters fit reasonably within each column; width ~ 0.9 data units.
        target_w = 0.9
        sx = (target_w / glyph_w)

        # Center the glyph at x (column), and place its bottom at y
        # We shift by -bb.x0/-bb.y0 to start at origin, then scale, then translate.
        trans = (
            Affine2D()
            .translate(-bb.x0, -bb.y0)      # move glyph bbox to start at (0,0)
            .scale(sx, sy)                  # scale to target width + exact bits height
            .translate(x - target_w / 2.0, y)  # center in column; bottom at y
        )

        patch = PathPatch(
            tp,
            transform=trans + ax.transData,
            facecolor=colors.get(letter, "black"),
            edgecolor="none",
        )
        ax.add_patch(patch)

    # Draw stacks
    for i, (pos_freqs, x) in enumerate(zip(pfm, positions)):
        # Total bits at this position
        R = ic[i]
        if R <= 0:
            continue

        # Compute letter heights in bits: p(base) * R
        heights = []
        for b in bases:
            p = float(pos_freqs.get(b, 0.0) or 0.0)
            if p < 0:
                p = 0.0
            heights.append((b, p * R))

        # Draw from smallest to largest so largest ends up on top
        heights.sort(key=lambda t: t[1])

        y0 = 0.0
        for b, h in heights:
            # Skip extremely tiny letters to reduce clutter
            if h < 1e-3:
                continue
            _add_letter(b, x, y0, h)
            y0 += h

    # Axis / annotations
    ax.axvline(x=0, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='End position')
    ax.set_xlabel('Position relative to transcript end (bp)', fontsize=12)
    ax.set_ylabel('Information Content (bits)', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_xlim(-upstream - 0.5, len(pfm) - upstream + 0.5)
    ax.set_ylim(0, 2.05)

    # Tidy look
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, linestyle='--', alpha=0.25)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=colors[b], label=b) for b in bases]
    legend_elements.append(plt.Line2D([0], [0], color='red', linestyle='--', label='End position'))
    ax.legend(handles=legend_elements, loc='upper right', ncol=5, fontsize=9, frameon=False)

    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved sequence logo to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to save sequence logo: {e}")
        plt.close(fig)
        return False



def _analyze_missed_peaks_comprehensive(
    missed_peaks: Dict[str, int],
    peaks_path: Path,
    read_end_positions: List[dict],
    iso_to_reads: Dict[str, List[str]],
    isoforms: Dict[str, dict],
    read_ends: Dict[str, dict],
    window: int,
    end_type: str,
    genome_path: Optional[Path] = None,
) -> dict:
    """
    Comprehensive analysis of missed peaks including read classification and truncation patterns.

    Args:
        missed_peaks: Dict of peak_id -> read_count for missed recoverable peaks
        peaks_path: Path to peaks BED file
        read_end_positions: List of read end position dicts (Chrom/Start/End/Strand)
        iso_to_reads: Isoform to reads mapping
        isoforms: Isoform info dict
        read_ends: Full read ends dict with TSS/TTS positions (from _parse_reads_bed_ends)
        window: Distance window used to define "recoverable"
        end_type: 'tss' for CAGE missed peaks, 'tts' for QuantSeq missed peaks
        genome_path: (optional) genome FASTA path (currently unused here)

    Returns:
        Dict with analysis results
    """
    with timed_section("analyze_missed_peaks_comprehensive"):
        if end_type not in ("tss", "tts"):
            raise ValueError(f"end_type must be 'tss' or 'tts', got: {end_type}")

    peaks_rows = _read_bed6(peaks_path)
    peak_id_to_info = {f"{p['Chrom']}_{p['Start']}_{p['End']}": p for p in peaks_rows}

    # Index read ends by chrom/strand for quick lookup
    reads_by_chrom_strand = defaultdict(list)
    for r in read_end_positions:
        reads_by_chrom_strand[(r["Chrom"], r["Strand"])].append(r)

    # Map read_id -> endpoint position for the requested end_type
    # NOTE: _parse_reads_bed_ends returns lowercase keys for chrom/strand, but keeps 'tss'/'tts' numeric fields
    read_id_to_end = {}
    # Create reverse index: (chrom, strand, pos) -> read_id for O(1) lookup
    pos_to_read_id = {}
    for rid, rinfo in read_ends.items():
        pos = rinfo.get(end_type)
        if pos is None:
            continue
        read_id_to_end[rid] = {
            "Chrom": rinfo["chrom"],
            "pos": pos,
            "Strand": rinfo["strand"],
        }
        # Index by (chrom, strand, pos) for fast reverse lookup
        key = (rinfo["chrom"], rinfo["strand"], pos)
        pos_to_read_id[key] = rid

    classification_summary = defaultdict(int)
    truncation_patterns = defaultdict(int)
    all_classifications = []

    for peak_id, read_count in missed_peaks.items():
        peak_info = peak_id_to_info.get(peak_id)
        if not peak_info:
            continue

        peak_pos = (peak_info["Start"] + peak_info["End"]) // 2

        # Find read IDs supporting this peak (within window around peak_pos)
        supporting_reads = []
        strands = [peak_info["Strand"]] if peak_info["Strand"] != "." else ["+", "-"]

        for strand in strands:
            key = (peak_info["Chrom"], strand)
            for r in reads_by_chrom_strand.get(key, []):
                # r is a 1bp interval [Start, End)
                if abs(r["Start"] - peak_pos) <= window:
                    # Use reverse index for O(1) lookup instead of O(N) iteration
                    lookup_key = (r["Chrom"], r["Strand"], r["Start"])
                    rid = pos_to_read_id.get(lookup_key)
                    if rid:
                        supporting_reads.append(rid)

        # Deduplicate + cap for performance
        supporting_reads = list(set(supporting_reads))[:100]

        # Calculate average read length for supporting reads
        read_lengths = []
        for rid in supporting_reads:
            if rid in read_ends:
                rinfo = read_ends[rid]
                read_length = abs(rinfo.get('end', 0) - rinfo.get('start', 0))
                if read_length > 0:
                    read_lengths.append(read_length)
        avg_read_length = statistics.mean(read_lengths) if read_lengths else 0

        # Classify reads (THIS is the key: pass end_type through)
        classification = _classify_missed_peak_reads(
            peak_id=peak_id,
            peak_info=peak_info,
            reads_supporting_peak=supporting_reads,
            iso_to_reads=iso_to_reads,
            isoforms=isoforms,
            read_ends=read_ends,
            end_type=end_type,
            window=window,
        )

        classification_summary["unassigned"] += classification.get("unassigned_count", 0)
        classification_summary["assigned_nearby"] += classification.get("assigned_nearby_count", 0)
        classification_summary["assigned_distant"] += classification.get("assigned_distant_count", 0)
        classification_summary["assigned_wrong_strand"] += classification.get("assigned_wrong_strand_count", 0)

        # Add average read length to classification
        classification['avg_read_length'] = avg_read_length

        # Truncation pattern analysis only makes sense for 5' ends
        if end_type == "tss" and supporting_reads:
            read_positions = [read_id_to_end[rid]["pos"] for rid in supporting_reads if rid in read_id_to_end]
            if read_positions:
                pattern, details = _characterize_truncation_pattern(
                    read_starts=read_positions,
                    peak_pos=peak_pos,
                    strand=peak_info["Strand"] if peak_info["Strand"] != "." else "+",
                )
                truncation_patterns[pattern] += 1
                classification["truncation_pattern"] = pattern
                classification["truncation_details"] = details

        all_classifications.append(classification)

    return {
        "n_missed_peaks": len(missed_peaks),
        "end_type": end_type,
        "classification_summary": dict(classification_summary),
        "truncation_patterns": dict(truncation_patterns) if end_type == "tss" else {},
        "peak_classifications": all_classifications[:100],
    }



def _plot_read_classification_summary(
    classification_summary: Dict[str, int],
    output_path: Path,
    title: str,
) -> bool:
    """Plot pie chart of read classification categories for missed peaks."""
    if not HAS_MATPLOTLIB:
        return False

    if not classification_summary or sum(classification_summary.values()) == 0:
        return False

    fig, ax = plt.subplots(figsize=(8, 8))

    labels = []
    sizes = []
    colors = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99']
    category_labels = {
        'unassigned': 'Unassigned',
        'assigned_nearby': 'Assigned Nearby',
        'assigned_distant': 'Assigned Distant',
        'assigned_wrong_strand': 'Wrong Strand',
    }

    for cat in ['unassigned', 'assigned_nearby', 'assigned_distant', 'assigned_wrong_strand']:
        count = classification_summary.get(cat, 0)
        if count > 0:
            labels.append(f"{category_labels[cat]}\n({count})")
            sizes.append(count)

    if not sizes:
        plt.close(fig)
        return False

    ax.pie(sizes, labels=labels, colors=colors[:len(sizes)], autopct='%1.1f%%',
           startangle=90, textprops={'fontsize': 10})
    ax.set_title(title, fontsize=14)

    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return True
    except Exception as e:
        logger.error(f"Failed to save classification pie chart: {e}")
        plt.close(fig)
        return False


def _plot_truncation_patterns(
    truncation_patterns: Dict[str, int],
    output_path: Path,
    title: str,
) -> bool:
    """Plot bar chart of truncation pattern categories."""
    if not HAS_MATPLOTLIB:
        return False

    if not truncation_patterns or sum(truncation_patterns.values()) == 0:
        return False

    fig, ax = plt.subplots(figsize=(10, 6))

    pattern_order = ['sharp', 'trailing', 'bimodal', 'dispersed', 'sparse']
    colors = ['#2ecc71', '#e74c3c', '#9b59b6', '#3498db', '#95a5a6']

    patterns = []
    counts = []
    bar_colors = []

    for i, pattern in enumerate(pattern_order):
        if pattern in truncation_patterns:
            patterns.append(pattern.capitalize())
            counts.append(truncation_patterns[pattern])
            bar_colors.append(colors[i])

    if not patterns:
        plt.close(fig)
        return False

    bars = ax.bar(patterns, counts, color=bar_colors, edgecolor='black', alpha=0.8)

    # Add count labels on bars
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(count), ha='center', va='bottom', fontsize=10)

    ax.set_xlabel('Truncation Pattern', fontsize=12)
    ax.set_ylabel('Number of Peaks', fontsize=12)
    ax.set_title(title, fontsize=14)

    # Add descriptions
    descriptions = {
        'Sharp': 'Reads cluster at peak',
        'Trailing': 'RT truncation gradient',
        'Bimodal': 'Two clusters (alt TSS?)',
        'Dispersed': 'Spread out reads',
        'Sparse': '<5 reads',
    }
    desc_text = '\n'.join([f"{k}: {v}" for k, v in descriptions.items() if k in patterns])
    ax.text(0.98, 0.98, desc_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return True
    except Exception as e:
        logger.error(f"Failed to save truncation patterns plot: {e}")
        plt.close(fig)
        return False


def _write_annotated_peaks_bed(
    peaks_path: Path,
    missed_peaks: Dict[str, int],
    analysis_results: dict,
    output_path: Path,
    end_type: str,
) -> int:
    """
    Write annotated CSV file for missed recoverable peaks with classification info.

    CSV format: chrom, start, end, igv_id, read_support, strand, read_classification, truncation_pattern
    """
    peaks_rows = _read_bed6(peaks_path)

    # Build lookup from analysis results
    classification_lookup = {}
    for cls in analysis_results.get('peak_classifications', []):
        peak_id = cls.get('peak_id')
        if peak_id:
            # Determine dominant classification
            counts = [
                ('unassigned', cls.get('unassigned_count', 0)),
                ('nearby', cls.get('assigned_nearby_count', 0)),
                ('distant', cls.get('assigned_distant_count', 0)),
                ('wrong_strand', cls.get('assigned_wrong_strand_count', 0)),
            ]
            dominant = max(counts, key=lambda x: x[1])[0]
            pattern = cls.get('truncation_pattern', 'NA')
            avg_read_length = cls.get('avg_read_length', 0)
            classification_lookup[peak_id] = (dominant, pattern, avg_read_length)

    count = 0
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        # Header
        writer.writerow(['#chrom', 'start', 'end', 'igv_id', 'read_support',
                        'strand', 'read_classification', 'truncation_pattern', 'avg_read_length'])

        for peak in peaks_rows:
            peak_id = f"{peak['Chrom']}_{peak['Start']}_{peak['End']}"
            if peak_id in missed_peaks:
                read_count = missed_peaks[peak_id]
                classification, pattern, avg_read_length = classification_lookup.get(peak_id, ('unknown', 'NA', 0))
                # Create IGV coordinate format
                igv_id = f"{peak['Chrom']}:{peak['Start']}-{peak['End']}"
                writer.writerow([
                    peak['Chrom'], peak['Start'], peak['End'],
                    igv_id, read_count, peak['Strand'],
                    classification, pattern if end_type == 'tss' else 'NA',
                    int(avg_read_length)
                ])
                count += 1

    logger.info(f"Wrote {count} annotated missed peaks to {output_path}")
    return count


def _analyze_all_recoverable_peaks_truncation(
    recoverable_peaks: Dict[str, int],
    captured_peaks: set,
    peaks_path: Path,
    read_end_positions: List[dict],
    read_ends: Dict[str, dict],
    window: int,
) -> Dict[str, dict]:
    """
    Analyze truncation patterns for ALL recoverable CAGE peaks (not just missed).

    Args:
        recoverable_peaks: Dict of peak_id -> read_count for recoverable peaks
        captured_peaks: Set of peak IDs that were captured by isoform ends
        peaks_path: Path to peaks BED file
        read_end_positions: List of read end position dicts
        read_ends: Full read ends dict with TSS/TTS positions
        window: Distance window

    Returns:
        Dict mapping peak_id -> {pattern, details, read_count, is_captured}
    """
    with timed_section("analyze_all_recoverable_peaks_truncation"):
        peaks_rows = _read_bed6(peaks_path)
        peak_id_to_info = {
            f"{p['Chrom']}_{p['Start']}_{p['End']}": p for p in peaks_rows
        }

    # Index read ends by position for quick lookup
    reads_by_chrom_strand = defaultdict(list)
    for r in read_end_positions:
        key = (r['Chrom'], r['Strand'])
        reads_by_chrom_strand[key].append(r)

    # Build read ID to end position mapping
    # Note: _parse_reads_bed_ends returns lowercase keys ('chrom', 'strand', 'tss', 'tts')
    read_id_to_end = {}
    # Create reverse index: (chrom, strand, pos) -> read_id for O(1) lookup
    pos_to_read_id = {}
    for rid, rinfo in read_ends.items():
        pos = rinfo.get('tss')
        if pos is not None:
            read_id_to_end[rid] = {
                'Chrom': rinfo['chrom'],
                'pos': pos,
                'Strand': rinfo['strand'],
            }
            # Index by (chrom, strand, pos) for fast reverse lookup
            key = (rinfo['chrom'], rinfo['strand'], pos)
            pos_to_read_id[key] = rid

    peak_patterns = {}

    for peak_id, read_count in recoverable_peaks.items():
        if peak_id not in peak_id_to_info:
            continue

        peak_info = peak_id_to_info[peak_id]
        peak_pos = (peak_info['Start'] + peak_info['End']) // 2

        # Find reads supporting this peak
        supporting_reads = []
        strands = [peak_info['Strand']] if peak_info['Strand'] != '.' else ['+', '-']
        for strand in strands:
            key = (peak_info['Chrom'], strand)
            for r in reads_by_chrom_strand.get(key, []):
                if abs(r['Start'] - peak_pos) <= window or abs(r['End'] - peak_pos) <= window:
                    # Use reverse index for O(1) lookup instead of O(N) nested loop
                    lookup_key = (r['Chrom'], r['Strand'], r['Start'])
                    rid = pos_to_read_id.get(lookup_key)
                    if rid:
                        supporting_reads.append(rid)

        supporting_reads = list(set(supporting_reads))[:100]  # Limit for performance

        # Get read positions for truncation pattern analysis
        read_positions = []
        for rid in supporting_reads:
            if rid in read_id_to_end:
                read_positions.append(read_id_to_end[rid]['pos'])

        pattern = 'sparse'
        details = {}
        if read_positions:
            pattern, details = _characterize_truncation_pattern(
                read_starts=read_positions,
                peak_pos=peak_pos,
                strand=peak_info['Strand'] if peak_info['Strand'] != '.' else '+',
            )

        peak_patterns[peak_id] = {
            'pattern': pattern,
            'details': details,
            'read_count': read_count,
            'is_captured': peak_id in captured_peaks,
            'chrom': peak_info['Chrom'],
            'start': peak_info['Start'],
            'end': peak_info['End'],
            'strand': peak_info['Strand'],
        }

    return peak_patterns


def _write_troubled_regions_tsv(
    missed_peaks: Dict[str, int],
    peak_patterns: Dict[str, dict],
    analysis_results: dict,
    output_path: Path,
) -> int:
    """
    Write TSV file for missed recoverable peaks with IGV-compatible coordinates.

    This TSV is peak-level (not per-read). It uses the classification results
    from _analyze_missed_peaks_comprehensive to explain why each peak was missed.

    Args:
        missed_peaks: Dict of peak_id -> read_count for missed recoverable peaks
        peak_patterns: Optional dict with chrom/start/end/strand/read_count/pattern for each peak.
                       (For CAGE you populate this via _analyze_all_recoverable_peaks_truncation.
                        For QuantSeq you can pass {} and coordinates will be derived from peak_id.)
        analysis_results: Output of _analyze_missed_peaks_comprehensive
        output_path: Path to output TSV file

    Returns:
        Number of rows written
    """
    # Build lookup from analysis results
    classification_lookup = {}
    for cls in analysis_results.get("peak_classifications", []):
        pid = cls.get("peak_id")
        if not pid:
            continue
        counts = {
            "unassigned": cls.get("unassigned_count", 0),
            "assigned_nearby": cls.get("assigned_nearby_count", 0),
            "assigned_distant": cls.get("assigned_distant_count", 0),
            "assigned_wrong_strand": cls.get("assigned_wrong_strand_count", 0),
        }
        dominant = max(counts.items(), key=lambda x: x[1])[0] if sum(counts.values()) > 0 else "unknown"
        classification_lookup[pid] = {"dominant": dominant, "counts": counts}

    def _parse_peak_id(peak_id: str) -> Tuple[str, int, int]:
        # expects chrom_start_end
        parts = peak_id.split("_")
        if len(parts) < 3:
            return ("", 0, 0)
        chrom = "_".join(parts[:-2])  # tolerate chrom names that include underscores
        try:
            start = int(parts[-2])
            end = int(parts[-1])
        except ValueError:
            start, end = 0, 0
        return chrom, start, end

    count = 0
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "igv_coords",
            "chrom", "start", "end", "strand",
            "read_support",
            "unassigned_reads", "assigned_nearby", "assigned_distant", "wrong_strand",
            "dominant_class",
            "truncation_pattern",
            "peak_id",
        ])

        for peak_id in missed_peaks:
            cls_info = classification_lookup.get(peak_id, {})
            counts = cls_info.get("counts", {})

            # Prefer coords from peak_patterns (CAGE case), else derive from peak_id (QuantSeq ok)
            pattern_info = peak_patterns.get(peak_id, {})
            if pattern_info.get("chrom"):
                chrom = pattern_info.get("chrom", "")
                start = int(pattern_info.get("start", 0))
                end = int(pattern_info.get("end", 0))
                strand = pattern_info.get("strand", ".")
                read_support = pattern_info.get("read_count", missed_peaks.get(peak_id, 0))
                trunc_pattern = pattern_info.get("pattern", "NA")
            else:
                chrom, start, end = _parse_peak_id(peak_id)
                strand = "."
                read_support = missed_peaks.get(peak_id, 0)
                trunc_pattern = "NA"

            igv_coords = f"{chrom}:{start}-{end}" if chrom else ""

            writer.writerow([
                igv_coords,
                chrom, start, end, strand,
                read_support,
                counts.get("unassigned", 0),
                counts.get("assigned_nearby", 0),
                counts.get("assigned_distant", 0),
                counts.get("assigned_wrong_strand", 0),
                cls_info.get("dominant", "unknown"),
                trunc_pattern,
                peak_id,
            ])
            count += 1

    logger.info(f"Wrote {count} troubled regions to {output_path}")
    return count



def _plot_all_truncation_patterns(
    peak_patterns: Dict[str, dict],
    output_path: Path,
    title: str,
) -> bool:
    """Plot bar chart of truncation patterns for ALL recoverable peaks, split by captured/missed."""
    if not HAS_MATPLOTLIB:
        return False

    if not peak_patterns:
        return False

    # Count patterns separately for captured vs missed
    captured_patterns = defaultdict(int)
    missed_patterns = defaultdict(int)

    for peak_id, info in peak_patterns.items():
        pattern = info.get('pattern', 'unknown')
        if info.get('is_captured'):
            captured_patterns[pattern] += 1
        else:
            missed_patterns[pattern] += 1

    pattern_order = ['sharp', 'trailing', 'bimodal', 'dispersed', 'sparse']
    colors_captured = '#2ecc71'  # Green
    colors_missed = '#e74c3c'    # Red

    fig, ax = plt.subplots(figsize=(12, 6))

    x = np.arange(len(pattern_order))
    width = 0.35

    captured_counts = [captured_patterns.get(p, 0) for p in pattern_order]
    missed_counts = [missed_patterns.get(p, 0) for p in pattern_order]

    bars1 = ax.bar(x - width/2, captured_counts, width, label='Captured', color=colors_captured, alpha=0.8)
    bars2 = ax.bar(x + width/2, missed_counts, width, label='Missed', color=colors_missed, alpha=0.8)

    # Add count labels
    for bar, count in zip(bars1, captured_counts):
        if count > 0:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    str(count), ha='center', va='bottom', fontsize=9)
    for bar, count in zip(bars2, missed_counts):
        if count > 0:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    str(count), ha='center', va='bottom', fontsize=9)

    ax.set_xlabel('Truncation Pattern', fontsize=12)
    ax.set_ylabel('Number of Peaks', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels([p.capitalize() for p in pattern_order])
    ax.legend()

    # Add descriptions
    descriptions = {
        'Sharp': 'Reads cluster at peak',
        'Trailing': 'RT truncation gradient',
        'Bimodal': 'Two clusters (alt TSS?)',
        'Dispersed': 'Spread out reads',
        'Sparse': '<5 reads',
    }
    desc_text = '\n'.join([f"{k}: {v}" for k, v in descriptions.items()])
    ax.text(0.98, 0.98, desc_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return True
    except Exception as e:
        logger.error(f"Failed to save truncation patterns plot: {e}")
        plt.close(fig)
        return False


# ────────────────────────── Core Metrics Calculation ──────────────────────────

def _tss_tts_metrics(
    iso_bed: Path,
    peaks: Dict[str, Optional[Path]],
    window: int,
    reads_bed: Optional[Path] = None,
    plot_output_dir: Optional[Path] = None,
    plot_prefix: str = "",
    test_regions_dir: Optional[Path] = None,
    genome_path: Optional[Path] = None,
    map_file: Optional[Path] = None,
    n_isoforms: Optional[int] = None,
) -> Tuple[dict, Dict[str, List[int]]]:
    """
    Compute precision/recall/F1 for TSS/TTS and extract signed distances for plotting.

    Args:
        iso_bed: Isoforms BED file
        peaks: Dict with keys 'prime5', 'prime3', 'ref_prime5', 'ref_prime3'
        window: Distance window for matching (default 50bp)
        reads_bed: Reads BED12 file (for recoverability and motif analysis)
        plot_output_dir: Directory to save distance histogram plots (optional)
        plot_prefix: Prefix for plot filenames (e.g., sample name)
        test_regions_dir: Directory to save recoverable peak BED files (optional)
        genome_path: Path to genome FASTA for motif analysis (optional)
        map_file: Path to isoform.read.map.txt for read classification (optional)
        n_isoforms: Pre-computed isoform count (optional, to avoid re-reading iso_bed)

    Returns:
        Tuple of (metrics dict, signed_distances dict)
        - metrics: precision/recall/F1 for each peak type
        - signed_distances: dict with keys 'prime5', 'prime3' containing distance lists
    """
    metrics: Dict[str, Optional[float]] = {
        "5prime_precision": None, "5prime_recall": None, "5prime_f1": None,
        "3prime_precision": None, "3prime_recall": None, "3prime_f1": None,
        "ref5prime_precision": None, "ref5prime_recall": None, "ref5prime_f1": None,
        "ref3prime_precision": None, "ref3prime_recall": None, "ref3prime_f1": None,
    }
    signed_distances: Dict[str, List[int]] = {
        "prime5": [],
        "prime3": [],
        "ref_prime5": [],
        "ref_prime3": [],
    }

    have_any = any(peaks.values())
    if not have_any:
        logger.warning("No peak files provided; TSS/TTS metrics will be None")
        return metrics, signed_distances

    if not _which("bedtools"):
        raise RuntimeError("bedtools is required for TSS/TTS metrics")

    # Count transcripts (use cached count if provided)
    if n_isoforms is not None:
        n_tx = n_isoforms
    else:
        n_tx = 0
        with open(iso_bed) as f:
            for line in f:
                if line and not line.startswith('#'):
                    n_tx += 1

    if n_tx == 0:
        logger.warning(f"Isoforms BED '{iso_bed}' has 0 transcripts")
        return metrics, signed_distances

    # Pre-compute read end positions for recoverability analysis
    read_tss_positions = []
    read_tts_positions = []
    if reads_bed and reads_bed.exists() and reads_bed.stat().st_size > 0:
        read_tss_positions = _extract_read_end_positions(reads_bed, 'tss')
        read_tts_positions = _extract_read_end_positions(reads_bed, 'tts')
        logger.debug(f"Read end positions: {len(read_tss_positions)} TSS, {len(read_tts_positions)} TTS")

    tmp_local = []
    try:
        # Extract single-position TSS and TTS BED files
        tss_bed, tts_bed = _extract_tss_tts_bed(iso_bed, tmp_local)
        tss_sorted = _prepare_bed6_sorted(tss_bed, tmp_local)
        tts_sorted = _prepare_bed6_sorted(tts_bed, tmp_local)

        def _side(key: str, label: str, endpoint_bed: Path, read_end_positions: List[dict]
                  ) -> Tuple[Optional[float], Optional[float], Optional[float], List[int], List[bool], Dict[str, int], List[List[str]]]:
            pth = peaks.get(key)
            if pth is None or not pth.exists() or pth.stat().st_size == 0:
                logger.warning(f"Peaks file for '{key}' missing or empty")
                return None, None, None, [], [], {}, []

            peaks_sorted = _prepare_bed6_sorted(pth, tmp_local)
            closest_rows = _run_closest(endpoint_bed, peaks_sorted)
            dist_map = _extract_distance_and_peak(closest_rows, label, window)

            # Extract signed distances for histogram plotting
            signed_dists = _extract_signed_distances(closest_rows)

            # Compute peak recoverability if read positions available
            recoverable_ids: Dict[str, int] = {}
            recov_mask = []
            if read_end_positions:
                recoverable_ids = _find_recoverable_peaks(pth, read_end_positions, window)
                recov_mask = _classify_isoform_recoverability(closest_rows, recoverable_ids)

            # Count isoforms within window
            m = sum(1 for dist, _ in dist_map.values() if dist <= window)
            precision = (m / n_tx) if n_tx else None

            # For recall
            peaks_rows = _read_bed6(pth)
            if not peaks_rows:
                logger.warning(f"Peaks file '{pth}' parsed to 0 intervals")

            # Read endpoint BED as rows
            endpoint_rows = _read_bed6(endpoint_bed)
            c, t = _vectorized_overlap_counts(endpoint_rows, peaks_rows, window)
            recall = (c / t) if t else None
            f1 = _safe_f1(precision, recall)

            return precision, recall, f1, signed_dists, recov_mask, recoverable_ids, closest_rows

        # Use TSS bed for 5' peaks, TTS bed for 3' peaks
        # CAGE/QuantSeq get read end positions for recoverability; ref peaks do not
        p5, r5, f5, dists_5prime, mask_5prime, recov_5prime, closest_5prime = _side("prime5", "5", tss_sorted, read_tss_positions)
        p3, r3, f3, dists_3prime, mask_3prime, recov_3prime, closest_3prime = _side("prime3", "3", tts_sorted, read_tts_positions)
        rp5, rr5, rf5, dists_ref5prime, _, _, _ = _side("ref_prime5", "ref5", tss_sorted, [])
        rp3, rr3, rf3, dists_ref3prime, _, _, _ = _side("ref_prime3", "ref3", tts_sorted, [])

        signed_distances["prime5"] = dists_5prime
        signed_distances["prime3"] = dists_3prime
        signed_distances["ref_prime5"] = dists_ref5prime
        signed_distances["ref_prime3"] = dists_ref3prime

        # Add recoverability counts to metrics
        if mask_5prime:
            metrics["5prime_recoverable_peaks"] = len(recov_5prime)
            n_total_5 = len(_read_bed6(peaks["prime5"])) if peaks.get("prime5") else 0
            metrics["5prime_total_peaks"] = n_total_5
            metrics["5prime_peak_recovery_rate"] = len(recov_5prime) / n_total_5 if n_total_5 else None
        if mask_3prime:
            metrics["3prime_recoverable_peaks"] = len(recov_3prime)
            n_total_3 = len(_read_bed6(peaks["prime3"])) if peaks.get("prime3") else 0
            metrics["3prime_total_peaks"] = n_total_3
            metrics["3prime_peak_recovery_rate"] = len(recov_3prime) / n_total_3 if n_total_3 else None

        metrics.update({
            "5prime_precision": p5, "5prime_recall": r5, "5prime_f1": f5,
            "3prime_precision": p3, "3prime_recall": r3, "3prime_f1": f3,
            "ref5prime_precision": rp5, "ref5prime_recall": rr5, "ref5prime_f1": rf5,
            "ref3prime_precision": rp3, "ref3prime_recall": rr3, "ref3prime_f1": rf3,
        })

        # Write recoverable peaks BED files
        if test_regions_dir:
            test_regions_dir.mkdir(parents=True, exist_ok=True)
            if recov_5prime and peaks.get("prime5"):
                _write_recoverable_peaks_bed(
                    peaks["prime5"], recov_5prime,
                    test_regions_dir / f"{plot_prefix}_recoverable_cage_peaks.bed")
            if recov_3prime and peaks.get("prime3"):
                _write_recoverable_peaks_bed(
                    peaks["prime3"], recov_3prime,
                    test_regions_dir / f"{plot_prefix}_recoverable_quantseq_peaks.bed")

        # Generate plots if output directory is specified
        if plot_output_dir and HAS_MATPLOTLIB:
            plot_output_dir.mkdir(parents=True, exist_ok=True)

            # Plot CAGE (5' end) distance histogram — colored by recoverability
            if dists_5prime:
                cage_plot_path = plot_output_dir / f"{plot_prefix}_cage_distance_histogram.png"
                if mask_5prime and len(mask_5prime) == len(dists_5prime):
                    _plot_distance_histogram_colored(
                        distances=dists_5prime,
                        recoverable_mask=mask_5prime,
                        output_path=cage_plot_path,
                        title=f"Distance to Nearest CAGE Peak (5' End)\n{plot_prefix}",
                    )
                else:
                    _plot_distance_histogram(
                        distances=dists_5prime,
                        output_path=cage_plot_path,
                        title=f"Distance to Nearest CAGE Peak (5' End)\n{plot_prefix}",
                    )

            # Plot QuantSeq (3' end) distance histogram — colored by recoverability
            if dists_3prime:
                quantseq_plot_path = plot_output_dir / f"{plot_prefix}_quantseq_distance_histogram.png"
                if mask_3prime and len(mask_3prime) == len(dists_3prime):
                    _plot_distance_histogram_colored(
                        distances=dists_3prime,
                        recoverable_mask=mask_3prime,
                        output_path=quantseq_plot_path,
                        title=f"Distance to Nearest QuantSeq Peak (3' End)\n{plot_prefix}",
                    )
                else:
                    _plot_distance_histogram(
                        distances=dists_3prime,
                        output_path=quantseq_plot_path,
                        title=f"Distance to Nearest QuantSeq Peak (3' End)\n{plot_prefix}",
                    )

            # Plot Reference TSS (5' end) distance histogram
            if dists_ref5prime:
                ref_tss_plot_path = plot_output_dir / f"{plot_prefix}_ref_tss_distance_histogram.png"
                _plot_distance_histogram(
                    distances=dists_ref5prime,
                    output_path=ref_tss_plot_path,
                    title=f"Distance to Nearest Reference TSS (5' End)\n{plot_prefix}",
                )

            # Plot Reference TTS (3' end) distance histogram
            if dists_ref3prime:
                ref_tts_plot_path = plot_output_dir / f"{plot_prefix}_ref_tts_distance_histogram.png"
                _plot_distance_histogram(
                    distances=dists_ref3prime,
                    output_path=ref_tts_plot_path,
                    title=f"Distance to Nearest Reference TTS (3' End)\n{plot_prefix}",
                )

            # Plot read support distribution for missed recoverable peaks
            # (recoverable by reads but not captured by any isoform end within window)
            missed_5 = {}
            missed_3 = {}
            captured_5 = set()
            captured_3 = set()
            if recov_5prime and closest_5prime:
                captured_5 = _find_captured_peaks(closest_5prime, window)
                missed_5 = {pid: cnt for pid, cnt in recov_5prime.items() if pid not in captured_5}
                if missed_5:
                    _plot_read_support_distribution(
                        recoverable_counts=missed_5,
                        output_path=plot_output_dir / f"{plot_prefix}_cage_peak_read_support.png",
                        title=f"Read Support for Missed Recoverable CAGE Peaks\n{plot_prefix}",
                    )
            if recov_3prime and closest_3prime:
                captured_3 = _find_captured_peaks(closest_3prime, window)
                missed_3 = {pid: cnt for pid, cnt in recov_3prime.items() if pid not in captured_3}
                if missed_3:
                    _plot_read_support_distribution(
                        recoverable_counts=missed_3,
                        output_path=plot_output_dir / f"{plot_prefix}_quantseq_peak_read_support.png",
                        title=f"Read Support for Missed Recoverable QuantSeq Peaks\n{plot_prefix}",
                    )

            # Comprehensive missed peak analysis with read classification and truncation patterns
            if map_file and map_file.exists() and reads_bed and reads_bed.exists():
                iso_to_reads = _parse_read_map(map_file)
                isoforms = _parse_isoform_ends(iso_bed)
                read_ends_full = _parse_reads_bed_ends(reads_bed)

                # Analyze ALL recoverable CAGE peaks (5' ends) for truncation patterns
                peak_patterns_5 = {}
                if recov_5prime and peaks.get("prime5"):
                    peak_patterns_5 = _analyze_all_recoverable_peaks_truncation(
                        recoverable_peaks=recov_5prime,
                        captured_peaks=captured_5,
                        peaks_path=peaks["prime5"],
                        read_end_positions=read_tss_positions,
                        read_ends=read_ends_full,
                        window=window,
                    )

                    # Count truncation patterns for all recoverable peaks
                    all_trunc_5 = defaultdict(int)
                    for pid, info in peak_patterns_5.items():
                        all_trunc_5[info.get('pattern', 'unknown')] += 1

                    if all_trunc_5:
                        metrics["5prime_all_pattern_sharp"] = all_trunc_5.get('sharp', 0)
                        metrics["5prime_all_pattern_trailing"] = all_trunc_5.get('trailing', 0)
                        metrics["5prime_all_pattern_bimodal"] = all_trunc_5.get('bimodal', 0)
                        metrics["5prime_all_pattern_dispersed"] = all_trunc_5.get('dispersed', 0)
                        metrics["5prime_all_pattern_sparse"] = all_trunc_5.get('sparse', 0)

                    # Plot truncation patterns for ALL recoverable peaks (captured vs missed)
                    _plot_all_truncation_patterns(
                        peak_patterns=peak_patterns_5,
                        output_path=plot_output_dir / f"{plot_prefix}_cage_all_truncation_patterns.png",
                        title=f"Truncation Patterns at ALL Recoverable CAGE Peaks\n{plot_prefix}",
                    )

                # Analyze missed CAGE peaks (5' ends) - classification and detailed analysis
                analysis_5 = {}
                if missed_5 and peaks.get("prime5"):
                    analysis_5 = _analyze_missed_peaks_comprehensive(
                        missed_peaks=missed_5,
                        peaks_path=peaks["prime5"],
                        read_end_positions=read_tss_positions,
                        iso_to_reads=iso_to_reads,
                        isoforms=isoforms,
                        read_ends=read_ends_full,
                        window=window,
                        end_type='tss',
                        genome_path=genome_path,
                    )

                    # Add summary metrics
                    if analysis_5.get('classification_summary'):
                        metrics["5prime_missed_unassigned"] = analysis_5['classification_summary'].get('unassigned', 0)
                        metrics["5prime_missed_nearby"] = analysis_5['classification_summary'].get('assigned_nearby', 0)
                        metrics["5prime_missed_distant"] = analysis_5['classification_summary'].get('assigned_distant', 0)

                    if analysis_5.get('truncation_patterns'):
                        metrics["5prime_missed_pattern_sharp"] = analysis_5['truncation_patterns'].get('sharp', 0)
                        metrics["5prime_missed_pattern_trailing"] = analysis_5['truncation_patterns'].get('trailing', 0)
                        metrics["5prime_missed_pattern_bimodal"] = analysis_5['truncation_patterns'].get('bimodal', 0)

                    # Plot classification summary
                    if analysis_5.get('classification_summary'):
                        _plot_read_classification_summary(
                            classification_summary=analysis_5['classification_summary'],
                            output_path=plot_output_dir / f"{plot_prefix}_cage_read_classification.png",
                            title=f"Read Assignment for Missed CAGE Peaks\n{plot_prefix}",
                        )

                    # Plot truncation patterns (missed only - legacy plot)
                    if analysis_5.get('truncation_patterns'):
                        _plot_truncation_patterns(
                            truncation_patterns=analysis_5['truncation_patterns'],
                            output_path=plot_output_dir / f"{plot_prefix}_cage_missed_truncation_patterns.png",
                            title=f"Truncation Patterns at Missed CAGE Peaks\n{plot_prefix}",
                        )

                    # Write annotated BED with classifications
                    if test_regions_dir:
                        _write_annotated_peaks_bed(
                            peaks_path=peaks["prime5"],
                            missed_peaks=missed_5,
                            analysis_results=analysis_5,
                            output_path=test_regions_dir / f"{plot_prefix}_missed_cage_peaks_annotated.csv",
                            end_type='tss',
                        )

                        # Write troubled regions TSV with IGV coordinates
                        _write_troubled_regions_tsv(
                            missed_peaks=missed_5,
                            peak_patterns=peak_patterns_5,
                            analysis_results=analysis_5,
                            output_path=test_regions_dir / f"{plot_prefix}_missed_cage_peaks.tsv",
                        )

                # Analyze missed QuantSeq peaks (3' ends)
                if missed_3 and peaks.get("prime3"):
                    analysis_3 = _analyze_missed_peaks_comprehensive(
                        missed_peaks=missed_3,
                        peaks_path=peaks["prime3"],
                        read_end_positions=read_tts_positions,
                        iso_to_reads=iso_to_reads,
                        isoforms=isoforms,
                        read_ends=read_ends_full,
                        window=window,
                        end_type='tts',
                        genome_path=genome_path,
                    )

                    # Add summary metrics
                    if analysis_3.get('classification_summary'):
                        metrics["3prime_missed_unassigned"] = analysis_3['classification_summary'].get('unassigned', 0)
                        metrics["3prime_missed_nearby"] = analysis_3['classification_summary'].get('assigned_nearby', 0)
                        metrics["3prime_missed_distant"] = analysis_3['classification_summary'].get('assigned_distant', 0)

                    # Plot classification summary
                    if analysis_3.get('classification_summary'):
                        _plot_read_classification_summary(
                            classification_summary=analysis_3['classification_summary'],
                            output_path=plot_output_dir / f"{plot_prefix}_quantseq_read_classification.png",
                            title=f"Read Assignment for Missed QuantSeq Peaks\n{plot_prefix}",
                        )

                    # Write annotated BED
                    if test_regions_dir:
                        _write_annotated_peaks_bed(
                            peaks_path=peaks["prime3"],
                            missed_peaks=missed_3,
                            analysis_results=analysis_3,
                            output_path=test_regions_dir / f"{plot_prefix}_missed_quantseq_peaks_annotated.csv",
                            end_type='tts',
                        )

                        # Write troubled regions TSV with IGV coordinates (no truncation patterns for QuantSeq)
                        _write_troubled_regions_tsv(
                            missed_peaks=missed_3,
                            peak_patterns={},  # No truncation analysis for QuantSeq
                            analysis_results=analysis_3,
                            output_path=test_regions_dir / f"{plot_prefix}_missed_quantseq_peaks.tsv",
                        )

                # Motif analysis around read ends (if genome available)
                if genome_path and genome_path.exists():
                    # TSS motif analysis
                    tss_motifs = _analyze_motifs_at_ends(
                        read_ends=read_ends_full,
                        genome_path=genome_path,
                        end_type='tss',
                    )
                    if tss_motifs.get('pfm'):
                        _plot_sequence_logo(
                            pfm=tss_motifs['pfm'],
                            output_path=plot_output_dir / f"{plot_prefix}_tss_motif_logo.png",
                            title=f"Sequence Context at Read 5' Ends\n{plot_prefix}",
                        )
                        metrics["tss_motif_mean_ic"] = tss_motifs.get('mean_ic')

                    # TTS motif analysis
                    tts_motifs = _analyze_motifs_at_ends(
                        read_ends=read_ends_full,
                        genome_path=genome_path,
                        end_type='tts',
                    )
                    if tts_motifs.get('pfm'):
                        _plot_sequence_logo(
                            pfm=tts_motifs['pfm'],
                            output_path=plot_output_dir / f"{plot_prefix}_tts_motif_logo.png",
                            title=f"Sequence Context at Read 3' Ends\n{plot_prefix}",
                        )
                        metrics["tts_motif_mean_ic"] = tts_motifs.get('mean_ic')

        return metrics, signed_distances
    finally:
        for f in tmp_local:
            try:
                f.unlink()
            except Exception:
                pass


def calculate_ted_metrics(
    iso_bed: Path,
    map_file: Path,
    bam_file: Optional[Path] = None,
    corrected_bed: Optional[Path] = None,
    reads_bed: Optional[Path] = None,
    prime5_peaks: Optional[Path] = None,
    prime3_peaks: Optional[Path] = None,
    ref_prime5_peaks: Optional[Path] = None,
    ref_prime3_peaks: Optional[Path] = None,
    window: int = 50,
    stage: str = "collapse",
    plot_output_dir: Optional[Path] = None,
    plot_prefix: str = "",
    test_regions_dir: Optional[Path] = None,
    genome_path: Optional[Path] = None,
) -> dict:
    """
    Calculate all TED metrics for a single sample.

    Args:
        iso_bed: Isoforms BED file (*.isoforms.bed)
        map_file: Read map file (*.isoform.read.map.txt)
        bam_file: BAM file (for alignment metrics)
        corrected_bed: Corrected BED (for collapse denominator)
        reads_bed: Reads BED12 file (for read-end entropy and recoverability analysis)
        prime5_peaks: Experimental 5' peaks (CAGE)
        prime3_peaks: Experimental 3' peaks (QuantSeq/PolyA)
        ref_prime5_peaks: Reference 5' peaks
        ref_prime3_peaks: Reference 3' peaks
        window: Distance window for TSS/TTS matching (default 50)
        stage: "collapse" or "transcriptome"
        plot_output_dir: Directory to save distance histogram plots (optional)
        plot_prefix: Prefix for plot filenames (e.g., sample name)
        test_regions_dir: Directory to save recoverable peak BED files (optional)
        genome_path: Path to genome FASTA for motif analysis (optional)

    Returns:
        Dict with all metrics
    """
    # Basic counts - read isoform names and extract gene IDs in single pass
    isoform_names = []
    genes = set()
    with open(iso_bed) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            if len(row) >= 4:
                name = row[3]
                isoform_names.append(name)
                # Extract gene IDs in same loop
                parts = str(name).split("_")
                for part in parts:
                    # Only count gene IDs (ENSG*, ENSMUSG*), not transcript IDs (ENST*, ENSMUST*)
                    if part.startswith(("ENSG", "ENSMUSG")):
                        genes.add(part.split(".")[0])
    
    n_iso = len(isoform_names)
    n_genes = len(genes)
    
    # Read assignment - parse map file once for both unique IDs and per-isoform stats
    assigned_read_ids, reads_per_iso_stats = _parse_map_file_comprehensive(map_file)
    assigned_unique_ids = len(assigned_read_ids)
    
    # Count by type if BAM available
    if bam_file and bam_file.exists():
        assigned_by_type = _count_assigned_reads_by_type(bam_file, assigned_read_ids)
        total_alignments = _count_total_alignments_bam(bam_file)
    else:
        assigned_by_type = {
            "assigned_primary": 0,
            "assigned_supplementary": 0,
            "assigned_total": 0,
        }
        total_alignments = {
            "total_primary": 0,
            "total_supplementary": 0,
            "total_alignments": 0,
        }
    
    # Denominator for assignment_rate
    if stage == "collapse" and corrected_bed and corrected_bed.exists():
        input_molecules = count_lines(corrected_bed)
    else:
        input_molecules = total_alignments.get("total_primary", 0) or None
    
    # TSS/TTS metrics
    peaks = {
        "prime5": prime5_peaks,
        "prime3": prime3_peaks,
        "ref_prime5": ref_prime5_peaks,
        "ref_prime3": ref_prime3_peaks,
    }
    tss_tts, signed_distances = _tss_tts_metrics(
        iso_bed, peaks, window,
        reads_bed=reads_bed,
        plot_output_dir=plot_output_dir,
        plot_prefix=plot_prefix,
        test_regions_dir=test_regions_dir,
        genome_path=genome_path,
        map_file=map_file,
        n_isoforms=n_iso,  # Pass cached count to avoid re-reading
    )

    # Read-end entropy analysis (skip if no reads or no assigned reads)
    entropy_metrics = {}
    if reads_bed and reads_bed.exists() and reads_bed.stat().st_size > 0 and assigned_unique_ids > 0:
        entropy_data = _compute_read_end_entropy(iso_bed, map_file, reads_bed)
        entropy_metrics = {
            "tss_entropy_mean": entropy_data["tss_entropy_mean"],
            "tss_entropy_median": entropy_data["tss_entropy_median"],
            "tts_entropy_mean": entropy_data["tts_entropy_mean"],
            "tts_entropy_median": entropy_data["tts_entropy_median"],
        }
        if plot_output_dir:
            _plot_read_end_entropy(entropy_data, plot_output_dir, plot_prefix)

    # Build result
    result = {
        "stage": stage,
        "isoforms_observed": n_iso,
        "genes_observed": n_genes,
        "assigned_unique_read_ids": assigned_unique_ids,
        "assigned_primary_alignments": assigned_by_type.get("assigned_primary", 0),
        "assigned_supplementary_alignments": assigned_by_type.get("assigned_supplementary", 0),
        "assigned_total_alignments": assigned_by_type.get("assigned_total", 0),
        "reads_per_isoform_mean": reads_per_iso_stats.get("reads_per_isoform_mean"),
        "reads_per_isoform_median": reads_per_iso_stats.get("reads_per_isoform_median"),
        "reads_per_isoform_min": reads_per_iso_stats.get("reads_per_isoform_min"),
        "reads_per_isoform_max": reads_per_iso_stats.get("reads_per_isoform_max"),
        "input_primary_alignments": total_alignments.get("total_primary", 0) or None,
        "input_supplementary_alignments": total_alignments.get("total_supplementary", 0) or None,
        "input_total_alignments": total_alignments.get("total_alignments", 0) or None,
        "assignment_rate": (assigned_unique_ids / input_molecules) if (input_molecules and input_molecules > 0) else None,
        "primary_alignment_utilization": (assigned_by_type.get("assigned_primary", 0) / total_alignments.get("total_primary", 1)) if total_alignments.get("total_primary", 0) > 0 else None,
        "total_alignment_utilization": (assigned_by_type.get("assigned_total", 0) / total_alignments.get("total_alignments", 1)) if total_alignments.get("total_alignments", 0) > 0 else None,
        **tss_tts,
        **entropy_metrics,
    }

    return result


def main():
    parser = argparse.ArgumentParser(
        description="TED: Transcriptome End Distribution metrics for FLAIR outputs"
    )
    parser.add_argument("--isoforms-bed", required=True, type=Path, help="Isoforms BED file (*.isoforms.bed)")
    parser.add_argument("--map-file", required=True, type=Path, help="Read map file (*.isoform.read.map.txt)")
    parser.add_argument("--bam", type=Path, help="BAM file for alignment metrics")
    parser.add_argument("--corrected-bed", type=Path, help="Corrected BED (for collapse stage)")
    parser.add_argument("--reads-bed", type=Path, help="Reads BED12 file (for read-end entropy analysis)")
    parser.add_argument("--prime5-peaks", type=Path, help="Experimental 5' peaks (CAGE)")
    parser.add_argument("--prime3-peaks", type=Path, help="Experimental 3' peaks (QuantSeq/PolyA)")
    parser.add_argument("--ref-prime5-peaks", type=Path, help="Reference 5' peaks")
    parser.add_argument("--ref-prime3-peaks", type=Path, help="Reference 3' peaks")
    parser.add_argument("--window", type=int, default=50, help="Distance window for TSS/TTS matching (default: 50)")
    parser.add_argument("--stage", choices=["collapse", "transcriptome"], default="collapse", help="Pipeline stage")
    parser.add_argument("--output", type=Path, required=True, help="Output TSV file")
    parser.add_argument("--plot-output-dir", type=Path, help="Directory to save distance histogram plots (CAGE and QuantSeq)")
    parser.add_argument("--test-regions-dir", type=Path, help="Directory to save recoverable peak BED files")
    parser.add_argument("--genome", type=Path, help="Genome FASTA file (for motif analysis)")
    parser.add_argument("--timing-output", type=Path, help="Output file for performance timing report")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    # Metadata arguments for result tracking
    parser.add_argument("--test-name", type=str, help="Test set name")
    parser.add_argument("--dataset-name", type=str, help="Dataset name")
    parser.add_argument("--align-mode", type=str, help="Alignment mode")
    parser.add_argument("--partition-mode", type=str, help="Partition mode")
    parser.add_argument("--pipeline-mode", type=str, help="Pipeline mode (process label)")
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Build plot prefix from metadata
    plot_prefix_parts = []
    if args.dataset_name:
        plot_prefix_parts.append(args.dataset_name)
    if args.align_mode:
        plot_prefix_parts.append(args.align_mode)
    if args.partition_mode:
        plot_prefix_parts.append(args.partition_mode)
    if args.pipeline_mode:
        plot_prefix_parts.append(args.pipeline_mode)
    if args.stage:
        plot_prefix_parts.append(args.stage)
    plot_prefix = "_".join(plot_prefix_parts) if plot_prefix_parts else "ted"

    # Calculate metrics
    logger.info(f"Calculating TED metrics for {args.isoforms_bed}")
    metrics = calculate_ted_metrics(
        iso_bed=args.isoforms_bed,
        map_file=args.map_file,
        bam_file=args.bam,
        corrected_bed=args.corrected_bed,
        reads_bed=args.reads_bed,
        prime5_peaks=args.prime5_peaks,
        prime3_peaks=args.prime3_peaks,
        ref_prime5_peaks=args.ref_prime5_peaks,
        ref_prime3_peaks=args.ref_prime3_peaks,
        window=args.window,
        stage=args.stage,
        plot_output_dir=args.plot_output_dir,
        plot_prefix=plot_prefix,
        test_regions_dir=args.test_regions_dir,
        genome_path=args.genome,
    )
    
    # Add metadata to metrics if provided
    if args.test_name:
        metrics = {'test_name': args.test_name, **metrics}
    if args.dataset_name:
        metrics = {'dataset': args.dataset_name, **metrics}
    if args.align_mode:
        metrics = {'align_mode': args.align_mode, **metrics}
    if args.partition_mode:
        metrics = {'partition_mode': args.partition_mode, **metrics}
    
    # Extract transcriptome_mode from pipeline_mode if provided
    # pipeline_mode format is typically "transcriptome_<mode>" or "collapse_<mode>"
    if args.pipeline_mode:
        if args.pipeline_mode.startswith('transcriptome_'):
            transcriptome_mode = args.pipeline_mode.replace('transcriptome_', '', 1)
        elif args.pipeline_mode.startswith('collapse_'):
            transcriptome_mode = args.pipeline_mode.replace('collapse_', '', 1)
        else:
            transcriptome_mode = args.pipeline_mode
        metrics = {'transcriptome_mode': transcriptome_mode, **metrics}
    
    # Write output as TSV
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=metrics.keys(), delimiter='\t')
        writer.writeheader()
        writer.writerow(metrics)
    
    logger.info(f"Wrote metrics to {args.output}")
    
    # Write timing report if requested
    if args.timing_output:
        with open(args.timing_output, 'w') as tf:
            tf.write("TED Performance Timing Report\n")
            tf.write("=" * 80 + "\n\n")
            
            # Sort by total time (sum of all calls)
            sorted_sections = sorted(
                _TIMING_DATA.items(),
                key=lambda x: sum(x[1]),
                reverse=True
            )
            
            total_time = sum(sum(times) for times in _TIMING_DATA.values())
            tf.write(f"Total instrumented time: {total_time:.2f}s\n\n")
            
            tf.write(f"{'Section':<50} {'Calls':>8} {'Total (s)':>12} {'Avg (s)':>12} {'% of Total':>12}\n")
            tf.write("-" * 100 + "\n")
            
            for section, times in sorted_sections:
                n_calls = len(times)
                total = sum(times)
                avg = total / n_calls
                pct = (total / total_time * 100) if total_time > 0 else 0
                tf.write(f"{section:<50} {n_calls:>8} {total:>12.2f} {avg:>12.2f} {pct:>11.1f}%\n")
        
        logger.info(f"Wrote timing report to {args.timing_output}")
    
    # Print summary
    logger.info(f"Isoforms: {metrics['isoforms_observed']}")
    logger.info(f"Genes: {metrics['genes_observed']}")
    if metrics['assignment_rate']:
        logger.info(f"Assignment rate: {metrics['assignment_rate']:.2%}")
    else:
        logger.info("Assignment rate: N/A")


if __name__ == "__main__":
    main()
