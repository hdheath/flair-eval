"""
Core TED (Transcript End Distance) metrics calculation.

Provides the main functions for calculating TSS/TTS accuracy metrics
and the top-level calculate_ted_metrics function.
"""

import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .utils import timed_section, which, safe_f1, count_lines, get_logger
from .bed_utils import (
    read_bed6,
    prepare_bed6_sorted,
    extract_tss_tts_bed,
    run_closest,
    extract_distance_and_peak,
    extract_signed_distances,
    vectorized_overlap_counts,
)
from .read_analysis import (
    parse_map_file_comprehensive,
    count_assigned_reads_by_type,
    count_total_alignments_bam,
    parse_isoform_ends,
    parse_read_map,
    parse_reads_bed_ends,
)
from .entropy import compute_read_end_entropy
from .peak_analysis import (
    find_recoverable_peaks,
    classify_isoform_recoverability,
    extract_read_end_positions,
    write_recoverable_peaks_bed,
    find_captured_peaks,
    analyze_missed_peaks_comprehensive,
    analyze_all_recoverable_peaks_truncation,
    write_annotated_peaks_bed,
    write_troubled_regions_tsv,
)
from .motif import analyze_motifs_at_ends
from .plots import (
    HAS_MATPLOTLIB,
    plot_distance_histogram,
    plot_distance_histogram_colored,
    plot_read_end_entropy,
    plot_read_support_distribution,
    plot_read_classification_summary,
    plot_truncation_patterns,
    plot_all_truncation_patterns,
    plot_sequence_logo,
)

logger = get_logger()


def tss_tts_metrics(
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
    max_count_cage: Optional[int] = None,
    max_count_quantseq: Optional[int] = None,
    max_count_ref_tss: Optional[int] = None,
    max_count_ref_tts: Optional[int] = None,
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

    if not which("bedtools"):
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
        read_tss_positions = extract_read_end_positions(reads_bed, 'tss')
        read_tts_positions = extract_read_end_positions(reads_bed, 'tts')
        logger.debug(f"Read end positions: {len(read_tss_positions)} TSS, {len(read_tts_positions)} TTS")

    tmp_local = []
    try:
        # Extract single-position TSS and TTS BED files
        with timed_section("extract_tss_tts_bed"):
            tss_bed, tts_bed = extract_tss_tts_bed(iso_bed, tmp_local)
            tss_sorted = prepare_bed6_sorted(tss_bed, tmp_local)
            tts_sorted = prepare_bed6_sorted(tts_bed, tmp_local)

        def _side(key: str, label: str, endpoint_bed: Path, read_end_positions: List[dict]
                  ) -> Tuple[Optional[float], Optional[float], Optional[float], List[int], List[bool], Dict[str, int], List[List[str]]]:
            pth = peaks.get(key)
            if pth is None or not pth.exists() or pth.stat().st_size == 0:
                logger.warning(f"Peaks file for '{key}' missing or empty")
                return None, None, None, [], [], {}, []

            peaks_sorted = prepare_bed6_sorted(pth, tmp_local)
            closest_rows = run_closest(endpoint_bed, peaks_sorted)
            dist_map = extract_distance_and_peak(closest_rows, label, window)

            # Extract signed distances for histogram plotting
            signed_dists = extract_signed_distances(closest_rows)

            # Compute peak recoverability if read positions available
            recoverable_ids: Dict[str, int] = {}
            recov_mask = []
            if read_end_positions:
                with timed_section(f"find_recoverable_peaks_{label}"):
                    recoverable_ids = find_recoverable_peaks(pth, read_end_positions, window)
                with timed_section(f"classify_isoform_recoverability_{label}"):
                    recov_mask = classify_isoform_recoverability(closest_rows, recoverable_ids)

            # Count isoforms within window
            m = sum(1 for dist, _ in dist_map.values() if dist <= window)
            precision = (m / n_tx) if n_tx else None

            # For recall
            peaks_rows = read_bed6(pth)
            if not peaks_rows:
                logger.warning(f"Peaks file '{pth}' parsed to 0 intervals")

            # Read endpoint BED as rows
            endpoint_rows = read_bed6(endpoint_bed)
            c, t = vectorized_overlap_counts(endpoint_rows, peaks_rows, window)
            recall = (c / t) if t else None
            f1 = safe_f1(precision, recall)

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
            n_total_5 = len(read_bed6(peaks["prime5"])) if peaks.get("prime5") else 0
            metrics["5prime_total_peaks"] = n_total_5
            metrics["5prime_peak_recovery_rate"] = len(recov_5prime) / n_total_5 if n_total_5 else None
        if mask_3prime:
            metrics["3prime_recoverable_peaks"] = len(recov_3prime)
            n_total_3 = len(read_bed6(peaks["prime3"])) if peaks.get("prime3") else 0
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
                write_recoverable_peaks_bed(
                    peaks["prime5"], recov_5prime,
                    test_regions_dir / f"{plot_prefix}_recoverable_cage_peaks.bed")
            if recov_3prime and peaks.get("prime3"):
                write_recoverable_peaks_bed(
                    peaks["prime3"], recov_3prime,
                    test_regions_dir / f"{plot_prefix}_recoverable_quantseq_peaks.bed")

        # Generate plots if output directory is specified
        if plot_output_dir and HAS_MATPLOTLIB:
            plot_output_dir.mkdir(parents=True, exist_ok=True)

            # Calculate maximum bin counts for each histogram type to ensure consistent y-axes
            # Use provided max_count values if available (for cross-run consistency),
            # otherwise auto-calculate from current data
            from evaluation.plots import calculate_histogram_max_count
            
            # For CAGE (experimental 5')
            if max_count_cage is not None:
                final_max_cage = max_count_cage
            else:
                final_max_cage = calculate_histogram_max_count(dists_5prime) if dists_5prime else None
            
            # For QuantSeq (experimental 3')
            if max_count_quantseq is not None:
                final_max_quantseq = max_count_quantseq
            else:
                final_max_quantseq = calculate_histogram_max_count(dists_3prime) if dists_3prime else None
            
            # For Reference TSS (5')
            if max_count_ref_tss is not None:
                final_max_ref_tss = max_count_ref_tss
            else:
                final_max_ref_tss = calculate_histogram_max_count(dists_ref5prime) if dists_ref5prime else None
            
            # For Reference TTS (3')
            if max_count_ref_tts is not None:
                final_max_ref_tts = max_count_ref_tts
            else:
                final_max_ref_tts = calculate_histogram_max_count(dists_ref3prime) if dists_ref3prime else None

            # Plot CAGE (5' end) distance histogram -- colored by recoverability
            if dists_5prime:
                cage_plot_path = plot_output_dir / f"{plot_prefix}_cage_distance_histogram.png"
                if mask_5prime and len(mask_5prime) == len(dists_5prime):
                    plot_distance_histogram_colored(
                        distances=dists_5prime,
                        recoverable_mask=mask_5prime,
                        output_path=cage_plot_path,
                        title=f"Distance to Nearest CAGE Peak (5' End)\n{plot_prefix}",
                        max_count=final_max_cage,
                    )
                else:
                    plot_distance_histogram(
                        distances=dists_5prime,
                        output_path=cage_plot_path,
                        title=f"Distance to Nearest CAGE Peak (5' End)\n{plot_prefix}",
                        max_count=final_max_cage,
                    )

            # Plot QuantSeq (3' end) distance histogram -- colored by recoverability
            if dists_3prime:
                quantseq_plot_path = plot_output_dir / f"{plot_prefix}_quantseq_distance_histogram.png"
                if mask_3prime and len(mask_3prime) == len(dists_3prime):
                    plot_distance_histogram_colored(
                        distances=dists_3prime,
                        recoverable_mask=mask_3prime,
                        output_path=quantseq_plot_path,
                        title=f"Distance to Nearest QuantSeq Peak (3' End)\n{plot_prefix}",
                        max_count=final_max_quantseq,
                    )
                else:
                    plot_distance_histogram(
                        distances=dists_3prime,
                        output_path=quantseq_plot_path,
                        title=f"Distance to Nearest QuantSeq Peak (3' End)\n{plot_prefix}",
                        max_count=final_max_quantseq,
                    )

            # Plot Reference TSS (5' end) distance histogram
            if dists_ref5prime:
                ref_tss_plot_path = plot_output_dir / f"{plot_prefix}_ref_tss_distance_histogram.png"
                plot_distance_histogram(
                    distances=dists_ref5prime,
                    output_path=ref_tss_plot_path,
                    title=f"Distance to Nearest Reference TSS (5' End)\n{plot_prefix}",
                    max_count=final_max_ref_tss,
                )

            # Plot Reference TTS (3' end) distance histogram
            if dists_ref3prime:
                ref_tts_plot_path = plot_output_dir / f"{plot_prefix}_ref_tts_distance_histogram.png"
                plot_distance_histogram(
                    distances=dists_ref3prime,
                    output_path=ref_tts_plot_path,
                    title=f"Distance to Nearest Reference TTS (3' End)\n{plot_prefix}",
                    max_count=final_max_ref_tts,
                )

            # Plot read support distribution for missed recoverable peaks
            # (recoverable by reads but not captured by any isoform end within window)
            missed_5 = {}
            missed_3 = {}
            captured_5 = set()
            captured_3 = set()
            if recov_5prime and closest_5prime:
                captured_5 = find_captured_peaks(closest_5prime, window)
                missed_5 = {pid: cnt for pid, cnt in recov_5prime.items() if pid not in captured_5}
                if missed_5:
                    plot_read_support_distribution(
                        recoverable_counts=missed_5,
                        output_path=plot_output_dir / f"{plot_prefix}_cage_peak_read_support.png",
                        title=f"Read Support for Missed Recoverable CAGE Peaks\n{plot_prefix}",
                    )
            if recov_3prime and closest_3prime:
                captured_3 = find_captured_peaks(closest_3prime, window)
                missed_3 = {pid: cnt for pid, cnt in recov_3prime.items() if pid not in captured_3}
                if missed_3:
                    plot_read_support_distribution(
                        recoverable_counts=missed_3,
                        output_path=plot_output_dir / f"{plot_prefix}_quantseq_peak_read_support.png",
                        title=f"Read Support for Missed Recoverable QuantSeq Peaks\n{plot_prefix}",
                    )

            # Comprehensive missed peak analysis with read classification and truncation patterns
            if map_file and map_file.exists() and reads_bed and reads_bed.exists():
                iso_to_reads = parse_read_map(map_file)
                isoforms = parse_isoform_ends(iso_bed)
                read_ends_full = parse_reads_bed_ends(reads_bed)

                # Analyze ALL recoverable CAGE peaks (5' ends) for truncation patterns
                peak_patterns_5 = {}
                if recov_5prime and peaks.get("prime5"):
                    peak_patterns_5 = analyze_all_recoverable_peaks_truncation(
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
                    plot_all_truncation_patterns(
                        peak_patterns=peak_patterns_5,
                        output_path=plot_output_dir / f"{plot_prefix}_cage_all_truncation_patterns.png",
                        title=f"Truncation Patterns at ALL Recoverable CAGE Peaks\n{plot_prefix}",
                    )

                # Analyze missed CAGE peaks (5' ends) - classification and detailed analysis
                analysis_5 = {}
                if missed_5 and peaks.get("prime5"):
                    analysis_5 = analyze_missed_peaks_comprehensive(
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
                        plot_read_classification_summary(
                            classification_summary=analysis_5['classification_summary'],
                            output_path=plot_output_dir / f"{plot_prefix}_cage_read_classification.png",
                            title=f"Read Assignment for Missed CAGE Peaks\n{plot_prefix}",
                        )

                    # Plot truncation patterns (missed only - legacy plot)
                    if analysis_5.get('truncation_patterns'):
                        plot_truncation_patterns(
                            truncation_patterns=analysis_5['truncation_patterns'],
                            output_path=plot_output_dir / f"{plot_prefix}_cage_missed_truncation_patterns.png",
                            title=f"Truncation Patterns at Missed CAGE Peaks\n{plot_prefix}",
                        )

                    # Write annotated BED with classifications
                    if test_regions_dir:
                        write_annotated_peaks_bed(
                            peaks_path=peaks["prime5"],
                            missed_peaks=missed_5,
                            analysis_results=analysis_5,
                            output_path=test_regions_dir / f"{plot_prefix}_missed_cage_peaks_annotated.csv",
                            end_type='tss',
                        )

                        # Write troubled regions TSV with IGV coordinates
                        write_troubled_regions_tsv(
                            missed_peaks=missed_5,
                            peak_patterns=peak_patterns_5,
                            analysis_results=analysis_5,
                            output_path=test_regions_dir / f"{plot_prefix}_missed_cage_peaks.tsv",
                        )

                # Analyze missed QuantSeq peaks (3' ends)
                if missed_3 and peaks.get("prime3"):
                    analysis_3 = analyze_missed_peaks_comprehensive(
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
                        plot_read_classification_summary(
                            classification_summary=analysis_3['classification_summary'],
                            output_path=plot_output_dir / f"{plot_prefix}_quantseq_read_classification.png",
                            title=f"Read Assignment for Missed QuantSeq Peaks\n{plot_prefix}",
                        )

                    # Write annotated BED
                    if test_regions_dir:
                        write_annotated_peaks_bed(
                            peaks_path=peaks["prime3"],
                            missed_peaks=missed_3,
                            analysis_results=analysis_3,
                            output_path=test_regions_dir / f"{plot_prefix}_missed_quantseq_peaks_annotated.csv",
                            end_type='tts',
                        )

                        # Write troubled regions TSV with IGV coordinates (no truncation patterns for QuantSeq)
                        write_troubled_regions_tsv(
                            missed_peaks=missed_3,
                            peak_patterns={},  # No truncation analysis for QuantSeq
                            analysis_results=analysis_3,
                            output_path=test_regions_dir / f"{plot_prefix}_missed_quantseq_peaks.tsv",
                        )

                # Motif analysis around read ends (if genome available)
                if genome_path and genome_path.exists():
                    # TSS motif analysis
                    tss_motifs = analyze_motifs_at_ends(
                        read_ends=read_ends_full,
                        genome_path=genome_path,
                        end_type='tss',
                    )
                    if tss_motifs.get('pfm'):
                        plot_sequence_logo(
                            pfm=tss_motifs['pfm'],
                            output_path=plot_output_dir / f"{plot_prefix}_tss_motif_logo.png",
                            title=f"Sequence Context at Read 5' Ends\n{plot_prefix}",
                        )
                        metrics["tss_motif_mean_ic"] = tss_motifs.get('mean_ic')

                    # TTS motif analysis
                    tts_motifs = analyze_motifs_at_ends(
                        read_ends=read_ends_full,
                        genome_path=genome_path,
                        end_type='tts',
                    )
                    if tts_motifs.get('pfm'):
                        plot_sequence_logo(
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
    max_count_cage: Optional[int] = None,
    max_count_quantseq: Optional[int] = None,
    max_count_ref_tss: Optional[int] = None,
    max_count_ref_tts: Optional[int] = None,
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
    with timed_section("parse_map_file_comprehensive"):
        assigned_read_ids, reads_per_iso_stats = parse_map_file_comprehensive(map_file)
        assigned_unique_ids = len(assigned_read_ids)

    # Count by type if BAM available
    if bam_file and bam_file.exists():
        with timed_section("count_alignments_by_type"):
            assigned_by_type = count_assigned_reads_by_type(bam_file, assigned_read_ids)
            total_alignments = count_total_alignments_bam(bam_file)
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
    with timed_section("tss_tts_metrics"):
        tss_tts, signed_distances = tss_tts_metrics(
            iso_bed, peaks, window,
            reads_bed=reads_bed,
            plot_output_dir=plot_output_dir,
            plot_prefix=plot_prefix,
            test_regions_dir=test_regions_dir,
            genome_path=genome_path,
            map_file=map_file,
            n_isoforms=n_iso,  # Pass cached count to avoid re-reading
            max_count_cage=max_count_cage,
            max_count_quantseq=max_count_quantseq,
            max_count_ref_tss=max_count_ref_tss,
            max_count_ref_tts=max_count_ref_tts,
        )

    # Read-end entropy analysis (skip if no reads or no assigned reads)
    entropy_metrics = {}
    if reads_bed and reads_bed.exists() and reads_bed.stat().st_size > 0 and assigned_unique_ids > 0:
        with timed_section("compute_read_end_entropy"):
            entropy_data = compute_read_end_entropy(iso_bed, map_file, reads_bed)
            entropy_metrics = {
                "tss_entropy_mean": entropy_data["tss_entropy_mean"],
                "tss_entropy_median": entropy_data["tss_entropy_median"],
                "tts_entropy_mean": entropy_data["tts_entropy_mean"],
                "tts_entropy_median": entropy_data["tts_entropy_median"],
            }
            if plot_output_dir:
                plot_read_end_entropy(entropy_data, plot_output_dir, plot_prefix)

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
