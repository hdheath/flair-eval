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
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s'
)
logger = logging.getLogger(__name__)


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

def _get_assigned_read_ids(map_path: Path) -> Set[str]:
    """Extract set of unique read IDs from isoform.read.map.txt."""
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
    """Count primary and supplementary alignments for assigned reads using samtools."""
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
        # Get read names and flags
        res = _run(["samtools", "view", "-F", "4", str(bam_path)])
        
        for line in res.stdout.strip().split('\n'):
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 2:
                continue
            
            read_name = fields[0]
            if read_name not in assigned_read_ids:
                continue
            
            flag = int(fields[1])
            is_supplementary = bool(flag & 0x800)
            is_secondary = bool(flag & 0x100)
            
            if is_supplementary:
                supplementary_count += 1
            elif not is_secondary:
                primary_count += 1
                
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


# ────────────────────────── Core Metrics Calculation ──────────────────────────

def _tss_tts_metrics(iso_bed: Path, peaks: Dict[str, Optional[Path]], window: int) -> dict:
    """
    Compute precision/recall/F1 for TSS/TTS.
    
    Args:
        iso_bed: Isoforms BED file
        peaks: Dict with keys 'prime5', 'prime3', 'ref_prime5', 'ref_prime3'
        window: Distance window for matching (default 50bp)
    
    Returns:
        Dict with precision/recall/F1 for each peak type
    """
    metrics: Dict[str, Optional[float]] = {
        "5prime_precision": None, "5prime_recall": None, "5prime_f1": None,
        "3prime_precision": None, "3prime_recall": None, "3prime_f1": None,
        "ref5prime_precision": None, "ref5prime_recall": None, "ref5prime_f1": None,
        "ref3prime_precision": None, "ref3prime_recall": None, "ref3prime_f1": None,
    }
    
    have_any = any(peaks.values())
    if not have_any:
        logger.warning("No peak files provided; TSS/TTS metrics will be None")
        return metrics

    if not _which("bedtools"):
        raise RuntimeError("bedtools is required for TSS/TTS metrics")

    # Count transcripts
    n_tx = 0
    with open(iso_bed) as f:
        for line in f:
            if line and not line.startswith('#'):
                n_tx += 1
    
    if n_tx == 0:
        logger.warning(f"Isoforms BED '{iso_bed}' has 0 transcripts")
        return metrics

    tmp_local = []
    try:
        # Extract single-position TSS and TTS BED files
        tss_bed, tts_bed = _extract_tss_tts_bed(iso_bed, tmp_local)
        tss_sorted = _prepare_bed6_sorted(tss_bed, tmp_local)
        tts_sorted = _prepare_bed6_sorted(tts_bed, tmp_local)

        def _side(key: str, label: str, endpoint_bed: Path) -> Tuple[Optional[float], Optional[float], Optional[float]]:
            pth = peaks.get(key)
            if pth is None or not pth.exists() or pth.stat().st_size == 0:
                logger.warning(f"Peaks file for '{key}' missing or empty")
                return None, None, None

            peaks_sorted = _prepare_bed6_sorted(pth, tmp_local)
            closest_rows = _run_closest(endpoint_bed, peaks_sorted)
            dist_map = _extract_distance_and_peak(closest_rows, label, window)
            
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

            return precision, recall, f1

        # Use TSS bed for 5' peaks, TTS bed for 3' peaks
        p5, r5, f5 = _side("prime5", "5", tss_sorted)
        p3, r3, f3 = _side("prime3", "3", tts_sorted)
        rp5, rr5, rf5 = _side("ref_prime5", "ref5", tss_sorted)
        rp3, rr3, rf3 = _side("ref_prime3", "ref3", tts_sorted)

        metrics.update({
            "5prime_precision": p5, "5prime_recall": r5, "5prime_f1": f5,
            "3prime_precision": p3, "3prime_recall": r3, "3prime_f1": f3,
            "ref5prime_precision": rp5, "ref5prime_recall": rr5, "ref5prime_f1": rf5,
            "ref3prime_precision": rp3, "ref3prime_recall": rr3, "ref3prime_f1": rf3,
        })
        return metrics
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
    prime5_peaks: Optional[Path] = None,
    prime3_peaks: Optional[Path] = None,
    ref_prime5_peaks: Optional[Path] = None,
    ref_prime3_peaks: Optional[Path] = None,
    window: int = 50,
    stage: str = "collapse",
) -> dict:
    """
    Calculate all TED metrics for a single sample.
    
    Args:
        iso_bed: Isoforms BED file (*.isoforms.bed)
        map_file: Read map file (*.isoform.read.map.txt)
        bam_file: BAM file (for alignment metrics)
        corrected_bed: Corrected BED (for collapse denominator)
        prime5_peaks: Experimental 5' peaks (CAGE)
        prime3_peaks: Experimental 3' peaks (QuantSeq/PolyA)
        ref_prime5_peaks: Reference 5' peaks
        ref_prime3_peaks: Reference 3' peaks
        window: Distance window for TSS/TTS matching (default 50)
        stage: "collapse" or "transcriptome"
    
    Returns:
        Dict with all metrics
    """
    # Basic counts - read isoform names
    isoform_names = []
    with open(iso_bed) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            if len(row) >= 4:
                isoform_names.append(row[3])
    
    n_iso = len(isoform_names)
    
    # Extract gene IDs (only ENSG* and ENSMUSG* patterns, not transcript IDs)
    genes = set()
    for name in isoform_names:
        parts = str(name).split("_")
        for part in parts:
            # Only count gene IDs (ENSG*, ENSMUSG*), not transcript IDs (ENST*, ENSMUST*)
            if part.startswith(("ENSG", "ENSMUSG")):
                genes.add(part.split(".")[0])
    n_genes = len(genes)
    
    # Read assignment
    assigned_read_ids = _get_assigned_read_ids(map_file)
    assigned_unique_ids = len(assigned_read_ids)
    
    # Reads per isoform statistics
    reads_per_iso_stats = _count_reads_per_isoform(map_file)
    
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
    tss_tts = _tss_tts_metrics(iso_bed, peaks, window)
    
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
    parser.add_argument("--prime5-peaks", type=Path, help="Experimental 5' peaks (CAGE)")
    parser.add_argument("--prime3-peaks", type=Path, help="Experimental 3' peaks (QuantSeq/PolyA)")
    parser.add_argument("--ref-prime5-peaks", type=Path, help="Reference 5' peaks")
    parser.add_argument("--ref-prime3-peaks", type=Path, help="Reference 3' peaks")
    parser.add_argument("--window", type=int, default=50, help="Distance window for TSS/TTS matching (default: 50)")
    parser.add_argument("--stage", choices=["collapse", "transcriptome"], default="collapse", help="Pipeline stage")
    parser.add_argument("--output", type=Path, required=True, help="Output TSV file")
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
    
    # Calculate metrics
    logger.info(f"Calculating TED metrics for {args.isoforms_bed}")
    metrics = calculate_ted_metrics(
        iso_bed=args.isoforms_bed,
        map_file=args.map_file,
        bam_file=args.bam,
        corrected_bed=args.corrected_bed,
        prime5_peaks=args.prime5_peaks,
        prime3_peaks=args.prime3_peaks,
        ref_prime5_peaks=args.ref_prime5_peaks,
        ref_prime3_peaks=args.ref_prime3_peaks,
        window=args.window,
        stage=args.stage,
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
    if args.pipeline_mode:
        metrics = {'pipeline_mode': args.pipeline_mode, **metrics}
    # Always include stage
    metrics = {**metrics, 'stage': args.stage}
    
    # Write output as TSV
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=metrics.keys(), delimiter='\t')
        writer.writeheader()
        writer.writerow(metrics)
    
    logger.info(f"Wrote metrics to {args.output}")
    
    # Print summary
    logger.info(f"Isoforms: {metrics['isoforms_observed']}")
    logger.info(f"Genes: {metrics['genes_observed']}")
    if metrics['assignment_rate']:
        logger.info(f"Assignment rate: {metrics['assignment_rate']:.2%}")
    else:
        logger.info("Assignment rate: N/A")


if __name__ == "__main__":
    main()
