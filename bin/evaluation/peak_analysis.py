"""
Peak analysis functions.

Provides utilities for analyzing peak recoverability, classifying missed peaks,
and generating annotated output files for troubleshooting.
"""

import csv
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from .utils import timed_section, get_logger
from .bed_utils import read_bed6
from .truncation import characterize_truncation_pattern

logger = get_logger()


def find_recoverable_peaks(
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
    peaks_rows = read_bed6(peaks_path)
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


def classify_isoform_recoverability(
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


def extract_read_end_positions(reads_bed: Path, end_type: str) -> List[dict]:
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


def write_recoverable_peaks_bed(
    peaks_path: Path,
    recoverable_ids: Dict[str, int],
    output_path: Path,
) -> int:
    """Write recoverable peak coordinates to a BED file.

    Returns number of peaks written.
    """
    peaks_rows = read_bed6(peaks_path)
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


def find_captured_peaks(
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


def get_reads_to_isoforms(iso_to_reads: Dict[str, List[str]]) -> Dict[str, str]:
    """Invert isoform->reads map to get read->isoform mapping."""
    read_to_iso = {}
    for iso_id, read_ids in iso_to_reads.items():
        for rid in read_ids:
            if rid not in read_to_iso:
                read_to_iso[rid] = iso_id
    return read_to_iso


def classify_missed_peak_reads(
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

    read_to_iso = get_reads_to_isoforms(iso_to_reads)

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


def analyze_missed_peaks_comprehensive(
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
        read_ends: Full read ends dict with TSS/TTS positions (from parse_reads_bed_ends)
        window: Distance window used to define "recoverable"
        end_type: 'tss' for CAGE missed peaks, 'tts' for QuantSeq missed peaks
        genome_path: (optional) genome FASTA path (currently unused here)

    Returns:
        Dict with analysis results
    """
    with timed_section("analyze_missed_peaks_comprehensive"):
        if end_type not in ("tss", "tts"):
            raise ValueError(f"end_type must be 'tss' or 'tts', got: {end_type}")

    peaks_rows = read_bed6(peaks_path)
    peak_id_to_info = {f"{p['Chrom']}_{p['Start']}_{p['End']}": p for p in peaks_rows}

    # Index read ends by chrom/strand for quick lookup
    reads_by_chrom_strand = defaultdict(list)
    for r in read_end_positions:
        reads_by_chrom_strand[(r["Chrom"], r["Strand"])].append(r)

    # Map read_id -> endpoint position for the requested end_type
    # NOTE: parse_reads_bed_ends returns lowercase keys for chrom/strand, but keeps 'tss'/'tts' numeric fields
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

    with timed_section(f"analyze_missed_peaks_loop_{end_type}"):
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
            classification = classify_missed_peak_reads(
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
                    pattern, details = characterize_truncation_pattern(
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
        "peak_classifications": all_classifications,
    }


def analyze_all_recoverable_peaks_truncation(
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
        peaks_rows = read_bed6(peaks_path)
        peak_id_to_info = {
            f"{p['Chrom']}_{p['Start']}_{p['End']}": p for p in peaks_rows
        }

    # Index read ends by position for quick lookup
    reads_by_chrom_strand = defaultdict(list)
    for r in read_end_positions:
        key = (r['Chrom'], r['Strand'])
        reads_by_chrom_strand[key].append(r)

    # Build read ID to end position mapping
    # Note: parse_reads_bed_ends returns lowercase keys ('chrom', 'strand', 'tss', 'tts')
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
            pattern, details = characterize_truncation_pattern(
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


def write_annotated_peaks_bed(
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
    peaks_rows = read_bed6(peaks_path)

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


def write_troubled_regions_tsv(
    missed_peaks: Dict[str, int],
    peak_patterns: Dict[str, dict],
    analysis_results: dict,
    output_path: Path,
) -> int:
    """
    Write TSV file for missed recoverable peaks with IGV-compatible coordinates.

    This TSV is peak-level (not per-read). It uses the classification results
    from analyze_missed_peaks_comprehensive to explain why each peak was missed.

    Args:
        missed_peaks: Dict of peak_id -> read_count for missed recoverable peaks
        peak_patterns: Optional dict with chrom/start/end/strand/read_count/pattern for each peak.
                       (For CAGE you populate this via analyze_all_recoverable_peaks_truncation.
                        For QuantSeq you can pass {} and coordinates will be derived from peak_id.)
        analysis_results: Output of analyze_missed_peaks_comprehensive
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
