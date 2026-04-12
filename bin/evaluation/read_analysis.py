"""
Read assignment analysis functions.

Provides utilities for parsing isoform-to-read mappings, counting alignments,
and extracting read position information from BAM and BED files.
"""

import csv
import statistics
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from .utils import run, which, timed_section, get_logger

logger = get_logger()


def parse_map_file_comprehensive(map_path: Path) -> Tuple[Set[str], dict]:
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


def get_assigned_read_ids(map_path: Path) -> Set[str]:
    """Extract set of unique read IDs from isoform.read.map.txt.

    Note: For better performance, use parse_map_file_comprehensive if you also
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


def count_reads_per_isoform(map_path: Path) -> dict:
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


def count_assigned_reads_by_type(bam_path: Path, assigned_read_ids: Set[str]) -> dict:
    """Count primary and supplementary alignments for assigned reads using samtools.

    Uses streaming to avoid loading entire BAM output into memory.
    """
    if not bam_path or not bam_path.exists() or not assigned_read_ids:
        return {
            "assigned_primary": 0,
            "assigned_supplementary": 0,
            "assigned_total": 0,
        }

    if not which("samtools"):
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


def count_total_alignments_bam(bam: Path) -> dict:
    """Count total primary and supplementary alignments in BAM."""
    if not which("samtools"):
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
        res = run(["samtools", "view", "-c", "-F", "2308", str(bam)])
        primary = int(res.stdout.strip() or "0")
    except Exception:
        primary = 0

    # Count supplementary alignments (include supplementary=0x800, exclude unmapped=0x4, secondary=0x100)
    try:
        res = run(["samtools", "view", "-c", "-f", "2048", "-F", "260", str(bam)])
        supplementary = int(res.stdout.strip() or "0")
    except Exception:
        supplementary = 0

    return {
        "total_primary": primary,
        "total_supplementary": supplementary,
        "total_alignments": primary + supplementary,
    }


def parse_isoform_ends(iso_bed: Path) -> Dict[str, dict]:
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


def parse_read_map(map_path: Path) -> Dict[str, List[str]]:
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


def parse_reads_bed_ends(reads_bed: Path) -> Dict[str, dict]:
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
