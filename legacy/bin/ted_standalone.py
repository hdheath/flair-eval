#!/usr/bin/env python3
"""
TED (Transcriptome End Distribution) - Standalone Version
Evaluates TSS/TTS accuracy and read assignment metrics for FLAIR outputs.
"""

import argparse
import io
import logging
import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd


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


def _read_bed6(path: Path) -> pd.DataFrame:
    """Read peaks as BED6 if possible; tolerate BED3 by adding Strand='.'"""
    try:
        return pd.read_csv(
            path, sep="\t", header=None, usecols=[0, 1, 2, 5],
            names=["Chrom", "Start", "End", "Strand"],
            dtype={"Chrom": str, "Start": np.int32, "End": np.int32, "Strand": str},
            comment="#",
        )
    except Exception:
        logger.warning(f"{path} does not look like BED6; treating as BED3 (strand='.')")
        df = pd.read_csv(
            path, sep="\t", header=None, usecols=[0, 1, 2],
            names=["Chrom", "Start", "End"],
            dtype={"Chrom": str, "Start": np.int32, "End": np.int32},
            comment="#",
        )
        df["Strand"] = "."
        return df


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
    df = pd.read_csv(
        iso_bed, sep="\t", header=None, comment="#",
        usecols=[0, 1, 2, 3, 4, 5],
        names=["Chrom", "Start", "End", "Name", "Score", "Strand"],
        dtype={"Chrom": str, "Start": np.int32, "End": np.int32, "Name": str, "Score": str, "Strand": str},
    )
    
    # Calculate single-position TSS/TTS
    df["TSS"] = np.where(df["Strand"] == "+", df["Start"], df["End"])
    df["TTS"] = np.where(df["Strand"] == "+", df["End"], df["Start"])
    
    # Create 1bp windows for TSS
    tss_bed = iso_bed.with_suffix(".tss.tmp")
    df_tss = df.assign(Start=df["TSS"], End=df["TSS"] + 1)[["Chrom", "Start", "End", "Name", "Score", "Strand"]]
    df_tss.to_csv(tss_bed, sep="\t", header=False, index=False)
    tmps.append(tss_bed)
    
    # Create 1bp windows for TTS
    tts_bed = iso_bed.with_suffix(".tts.tmp")
    df_tts = df.assign(Start=df["TTS"], End=df["TTS"] + 1)[["Chrom", "Start", "End", "Name", "Score", "Strand"]]
    df_tts.to_csv(tts_bed, sep="\t", header=False, index=False)
    tmps.append(tts_bed)
    
    logger.debug(f"Extracted {len(df)} TSS and TTS positions from {iso_bed}")
    return tss_bed, tts_bed


def _run_closest(a: Path, b: Path) -> pd.DataFrame:
    """Run bedtools closest."""
    res = _run(["bedtools", "closest", "-a", str(a), "-b", str(b), "-s", "-d"])
    return pd.read_csv(io.StringIO(res.stdout), sep="\t", header=None, comment="#")


def _extract_distance_and_peak(df: pd.DataFrame, label: str, max_dist: int) -> pd.DataFrame:
    """Extract distances to nearest peak."""
    orig_cols = df.shape[1]
    tx_series = df.iloc[:, 3].astype(str)
    dist_series = pd.to_numeric(df.iloc[:, orig_cols - 1], errors="coerce").fillna(max_dist + 1).astype(int)
    peak_series = (df.iloc[:, 6].astype(str) + "_" + df.iloc[:, 7].astype(str) + "_" + df.iloc[:, 8].astype(str))
    out = pd.DataFrame({"isoform_id": tx_series, f"{label}_dist": dist_series, f"{label}_peak": peak_series})
    return out.groupby("isoform_id", as_index=False).first()


def _vectorized_overlap_counts(bed_df: pd.DataFrame, peaks_df: pd.DataFrame, window: int) -> Tuple[int, int]:
    """Return (#peaks matched by ≥1 isoform within window, total peaks)."""
    consumed = np.zeros(len(peaks_df), dtype=bool)
    for (c, s), chunk in bed_df.groupby(["Chrom", "Strand"], sort=False):
        mask = (peaks_df["Chrom"] == c) & (
            (peaks_df["Strand"] == s) | (peaks_df["Strand"] == ".")
        )
        idxs = np.flatnonzero(mask)
        if idxs.size == 0:
            continue
        sub = peaks_df.iloc[idxs]
        starts = sub["Start"].to_numpy()
        ends = sub["End"].to_numpy()
        order = np.argsort(starts)
        starts = starts[order]
        ends = ends[order]
        global_i = idxs[order]
        for _, row in chunk.iterrows():
            lo, hi = int(row["Start"]) - window, int(row["End"]) + window
            r = np.searchsorted(starts, hi, side="right")
            if r == 0:
                continue
            hits = np.nonzero(ends[:r] >= lo)[0]
            if hits.size == 0:
                continue
            matched = global_i[hits]
            new_hits = matched[~consumed[matched]]
            consumed[new_hits] = True
    return int(consumed.sum()), int(len(peaks_df))


def _safe_f1(p: Optional[float], r: Optional[float]) -> Optional[float]:
    """Calculate F1 score."""
    if p is None or r is None:
        return None
    s = p + r
    return (2 * p * r / s) if s > 0 else 0.0


# ────────────────────────── Read Assignment Functions ──────────────────────────

def _get_assigned_read_ids(map_path: Path) -> set:
    """Extract set of unique read IDs from isoform.read.map.txt."""
    if not map_path.exists() or map_path.stat().st_size == 0:
        return set()
    uniq = set()
    with open(map_path) as fh:
        for line in fh:
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
    with open(map_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or "\t" not in line:
                continue
            _, rhs = line.split("\t", 1)
            # Count unique reads for this isoform
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
        "reads_per_isoform_mean": float(np.mean(read_counts)),
        "reads_per_isoform_median": float(np.median(read_counts)),
        "reads_per_isoform_min": int(np.min(read_counts)),
        "reads_per_isoform_max": int(np.max(read_counts)),
    }


def _count_assigned_reads_by_type(bam_path: Path, assigned_read_ids: set) -> dict:
    """Count primary and supplementary alignments for assigned reads."""
    if not bam_path or not bam_path.exists() or not assigned_read_ids:
        return {
            "assigned_primary": 0,
            "assigned_supplementary": 0,
            "assigned_total": 0,
        }
    
    try:
        import pysam
    except ImportError:
        logger.warning("pysam not available; cannot count assigned reads by type")
        return {
            "assigned_primary": 0,
            "assigned_supplementary": 0,
            "assigned_total": 0,
        }
    
    primary_count = 0
    supplementary_count = 0
    
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bf:
            for read in bf:
                if read.is_unmapped:
                    continue
                if read.query_name not in assigned_read_ids:
                    continue
                
                if read.is_supplementary:
                    supplementary_count += 1
                elif not read.is_secondary:
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

    full_bed_df = pd.read_csv(
        iso_bed, sep="\t", header=None, comment="#",
        usecols=[0, 1, 2, 3, 5],
        names=["Chrom", "Start", "End", "Name", "Strand"],
        dtype={"Chrom": str, "Start": np.int32, "End": np.int32, "Name": str, "Strand": str},
    )
    n_tx = len(full_bed_df)
    if n_tx == 0:
        logger.warning(f"Isoforms BED '{iso_bed}' has 0 transcripts")
        return metrics
        
    bed_df = full_bed_df[["Chrom", "Start", "End", "Strand"]]

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
            dfc = _run_closest(endpoint_bed, peaks_sorted)
            ddf = _extract_distance_and_peak(dfc, label, window)
            m = int((ddf[f"{label}_dist"] <= window).sum())
            precision = (m / n_tx) if n_tx else None

            peaks_df = _read_bed6(pth)
            if peaks_df.empty:
                logger.warning(f"Peaks file '{pth}' parsed to 0 intervals")
            
            # For recall, use single-position endpoint BED
            endpoint_df = pd.read_csv(
                endpoint_bed, sep="\t", header=None, comment="#",
                usecols=[0, 1, 2, 5],
                names=["Chrom", "Start", "End", "Strand"],
                dtype={"Chrom": str, "Start": np.int32, "End": np.int32, "Strand": str},
            )
            c, t = _vectorized_overlap_counts(endpoint_df, peaks_df, window)
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
    # Basic counts
    df_iso = pd.read_csv(iso_bed, sep="\t", header=None, comment="#", usecols=[3], names=["Name"], dtype=str)
    n_iso = len(df_iso)
    
    # Extract gene IDs (simple version - looks for ENS* patterns)
    genes = set()
    for name in df_iso["Name"]:
        parts = str(name).split("_")
        for part in parts:
            if part.startswith(("ENS", "ENSG", "ENST", "ENSMUSG", "ENSMUST")):
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
    
    # Denominator
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
        "input_molecules": int(input_molecules) if (input_molecules is not None) else None,
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
    
    # Write output
    df = pd.DataFrame([metrics])
    df.to_csv(args.output, sep="\t", index=False)
    logger.info(f"Wrote metrics to {args.output}")
    
    # Print summary
    logger.info(f"Isoforms: {metrics['isoforms_observed']}")
    logger.info(f"Genes: {metrics['genes_observed']}")
    logger.info(f"Assignment rate: {metrics['assignment_rate']:.2%}" if metrics['assignment_rate'] else "Assignment rate: N/A")


if __name__ == "__main__":
    main()
