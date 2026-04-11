#!/usr/bin/env python3
"""
Plot transcript models with optional strict region clipping and orthogonal signal.

Panels:
1) Called transcript models
2) Main long-read alignments (isoform-color matched)
3) Optional CAGE signal panel (if bedGraph provided)
4) Optional QuantSeq signal panel (if bedGraph provided)
"""

from __future__ import annotations

import argparse
import colorsys
import hashlib
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / "evaluation"))

import matplotlib.colors as mcolors
import matplotlib.patches as mplpatches
import matplotlib.pyplot as plt
import pysam

from pub_style import apply_rc  # Nature Portfolio rcParams
apply_rc()

# Shared vertical geometry defaults (in axis y-units).
READ_BLOCK_HEIGHT = 0.70
# Isoform block scaling relative to read block height.
MODEL_BLOCK_MIN_READ_RATIO = 1.4
MODEL_BLOCK_TARGET_READ_RATIO = 1.5
MODEL_BLOCK_MAX_READ_RATIO = 2.0
MODEL_ROW_GAP_FRACTION = 0.10
MODEL_PANEL_MIN_SPAN = 3.0
UNASSIGNED_GROUP_ID = "__unassigned__"


def assign_rows(intervals: Sequence[Tuple[int, int]]) -> List[int]:
    """Greedy non-overlapping row assignment."""
    indexed = sorted(enumerate(intervals), key=lambda x: x[1][0])
    row_assignment = [0] * len(intervals)
    row_ends: List[int] = []

    for orig_idx, (start, end) in indexed:
        placed = False
        for row_idx, last_end in enumerate(row_ends):
            if start >= last_end:
                row_assignment[orig_idx] = row_idx
                row_ends[row_idx] = end
                placed = True
                break
        if not placed:
            row_assignment[orig_idx] = len(row_ends)
            row_ends.append(end)
    return row_assignment


def parse_region(region_text: str) -> Tuple[str, int, int]:
    """Parse 'chr:start-end' into (chrom, start, end)."""
    if ":" not in region_text or "-" not in region_text:
        raise ValueError(f"Invalid region format '{region_text}'. Expected chr:start-end")
    chrom, coords = region_text.split(":", 1)
    start_s, end_s = coords.split("-", 1)
    start = int(start_s.replace(",", ""))
    end = int(end_s.replace(",", ""))
    if end <= start:
        raise ValueError(f"Invalid region bounds '{region_text}'. end must be greater than start")
    return chrom, start, end


def load_isoform_read_map(readmap_path: Path) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    iso_to_reads: Dict[str, List[str]] = {}
    read_to_iso: Dict[str, str] = {}
    with readmap_path.open() as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            parts = text.split("\t", 1)
            if len(parts) != 2:
                continue
            iso_id, reads_raw = parts
            reads = [r for r in reads_raw.split(",") if r]
            iso_to_reads[iso_id] = reads
            for read_id in reads:
                read_to_iso[read_id] = iso_id
    return iso_to_reads, read_to_iso


def load_isoforms_bed(isoforms_path: Path) -> Dict[str, dict]:
    isoforms: Dict[str, dict] = {}
    with isoforms_path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            iso_id = fields[3]
            strand = fields[5]
            block_sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
            block_starts = [int(x) for x in fields[11].rstrip(",").split(",") if x]
            blocks = [
                (start + block_starts[i], start + block_starts[i] + block_sizes[i])
                for i in range(min(len(block_sizes), len(block_starts)))
            ]
            isoforms[iso_id] = {
                "chrom": chrom,
                "strand": strand,
                "start": start,
                "end": end,
                "blocks": blocks,
            }
    return isoforms


def filter_isoforms_to_region(
    isoforms: Dict[str, dict],
    chrom: str,
    region_start: int,
    region_end: int,
) -> Dict[str, dict]:
    filtered: Dict[str, dict] = {}
    for iso_id, info in isoforms.items():
        if info["chrom"] != chrom:
            continue
        if info["end"] <= region_start or info["start"] >= region_end:
            continue
        filtered[iso_id] = info
    return filtered


def load_peaks(
    bed_path: Optional[Path],
    chrom: str,
    region_start: int,
    region_end: int,
) -> List[Tuple[int, int]]:
    if bed_path is None or not bed_path.exists():
        return []
    peaks: List[Tuple[int, int]] = []
    with bed_path.open() as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            fields = text.split("\t")
            if len(fields) < 3:
                continue
            if fields[0] != chrom:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            if end <= region_start or start >= region_end:
                continue
            peaks.append((max(region_start, start), min(region_end, end)))
    return peaks


def pick_region(
    isoforms: Dict[str, dict],
    bam_path: Path,
    strict_region: Optional[Tuple[str, int, int]] = None,
) -> Tuple[str, int, int]:
    if strict_region is not None:
        return strict_region

    if isoforms:
        chroms = sorted({v["chrom"] for v in isoforms.values()})
        chrom = chroms[0]
        start = min(v["start"] for v in isoforms.values())
        end = max(v["end"] for v in isoforms.values())
        return chrom, start, end

    # Fallback when no isoforms are present.
    with pysam.AlignmentFile(str(bam_path), "rb") as sam:
        first = None
        for aln in sam:
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            first = aln
            break
    if first is None:
        raise RuntimeError("Could not infer plotting region: no isoforms and no mapped reads found.")
    chrom = first.reference_name
    start = max(0, first.reference_start - 1000)
    end = first.reference_end + 1000
    return chrom, start, end


def load_primary_alignments(
    bam_path: Path,
    chrom: str,
    region_start: int,
    region_end: int,
    max_reads: int,
) -> Dict[str, pysam.AlignedSegment]:
    alignments: Dict[str, pysam.AlignedSegment] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as sam:
        for aln in sam.fetch(chrom, region_start, region_end):
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            # keep one alignment per read name
            if aln.query_name not in alignments:
                alignments[aln.query_name] = aln
            # max_reads is kept for backward compatibility; balancing is applied later.
            _ = max_reads
    return alignments


def _stable_read_rank(read_id: str) -> int:
    digest = hashlib.md5(read_id.encode("utf-8")).hexdigest()
    return int(digest[:12], 16)


def _count_reads_by_group(
    read_ids: Iterable[str],
    read_to_iso: Dict[str, str],
) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for read_id in read_ids:
        group_id = read_to_iso.get(read_id, UNASSIGNED_GROUP_ID)
        counts[group_id] = counts.get(group_id, 0) + 1
    return counts


def _junction_chain_from_alignment(aln: pysam.AlignedSegment) -> Tuple[Tuple[int, int], ...]:
    """
    Extract ordered intron chain from CIGAR (N operations) in reference coordinates.
    """
    if aln.reference_start is None:
        return tuple()
    pos = int(aln.reference_start)
    chain: List[Tuple[int, int]] = []
    for op, length in (aln.cigartuples or []):
        if op in (0, 2, 3, 7, 8):  # M, D, N, =, X consume reference
            if op == 3:
                chain.append((pos, pos + int(length)))
            pos += int(length)
        elif op in (1, 4, 5, 6):  # I, S, H, P do not consume reference
            continue
    return tuple(chain)


def _cap_redundant_group_reads(
    read_ids: Sequence[str],
    alignments: Dict[str, pysam.AlignedSegment],
    max_redundant_pattern_reads: int,
    pattern_end_tolerance: int,
) -> Tuple[List[str], int]:
    """
    Limit near-duplicate read patterns within one group.

    Pattern definition:
      - identical splice-junction chain
      - start/end coordinates within +/- pattern_end_tolerance bp
    """
    max_per_pattern = int(max_redundant_pattern_reads)
    if max_per_pattern <= 0:
        return list(read_ids), 0

    tol = max(0, int(pattern_end_tolerance))
    bin_size = max(1, tol + 1)

    # Per-chain cluster state
    # clusters_by_chain[chain] = [[rep_start, rep_end, kept_count], ...]
    clusters_by_chain: Dict[Tuple[Tuple[int, int], ...], List[List[int]]] = {}
    # grid index to speed up tolerance lookups:
    # index_by_chain[chain][(start_bin, end_bin)] -> [cluster_idx, ...]
    index_by_chain: Dict[Tuple[Tuple[int, int], ...], Dict[Tuple[int, int], List[int]]] = {}

    kept: List[str] = []
    dropped = 0

    for read_id in read_ids:
        aln = alignments.get(read_id)
        if aln is None:
            continue
        if aln.reference_start is None or aln.reference_end is None:
            kept.append(read_id)
            continue

        start = int(aln.reference_start)
        end = int(aln.reference_end)
        chain = _junction_chain_from_alignment(aln)

        clusters = clusters_by_chain.setdefault(chain, [])
        grid = index_by_chain.setdefault(chain, {})

        sb = start // bin_size
        eb = end // bin_size

        chosen_idx: Optional[int] = None
        for ds in (-1, 0, 1):
            for de in (-1, 0, 1):
                for idx in grid.get((sb + ds, eb + de), []):
                    rep_start, rep_end, _ = clusters[idx]
                    if abs(start - rep_start) <= tol and abs(end - rep_end) <= tol:
                        chosen_idx = idx
                        break
                if chosen_idx is not None:
                    break
            if chosen_idx is not None:
                break

        if chosen_idx is None:
            chosen_idx = len(clusters)
            clusters.append([start, end, 0])
            grid.setdefault((sb, eb), []).append(chosen_idx)

        if clusters[chosen_idx][2] < max_per_pattern:
            clusters[chosen_idx][2] += 1
            kept.append(read_id)
        else:
            dropped += 1

    return kept, dropped


def select_balanced_main_alignments(
    alignments: Dict[str, pysam.AlignedSegment],
    read_to_iso: Dict[str, str],
    max_reads: int,
    max_redundant_pattern_reads: int,
    pattern_end_tolerance: int,
) -> Tuple[Dict[str, pysam.AlignedSegment], Dict[str, int], Dict[str, int], int, int]:
    """
    Apply a balanced cap over read groups (isoforms + unassigned).

    Strategy:
      1) Group reads by assigned isoform id, with explicit unassigned group.
      2) Deterministically order reads within each group.
      3) Round-robin pick across groups until max_reads is reached.
    """
    full_counts = _count_reads_by_group(alignments.keys(), read_to_iso)
    if not alignments:
        return {}, full_counts, {}, 0, 0

    grouped: Dict[str, List[str]] = {}
    for read_id, aln in alignments.items():
        group_id = read_to_iso.get(read_id, UNASSIGNED_GROUP_ID)
        grouped.setdefault(group_id, []).append(read_id)

    redundant_dropped = 0
    for group_id in grouped:
        ordered = sorted(
            grouped[group_id],
            key=lambda rid: (
                _stable_read_rank(rid),
                alignments[rid].reference_start if alignments[rid].reference_start is not None else -1,
            )
        )
        kept, dropped = _cap_redundant_group_reads(
            ordered,
            alignments=alignments,
            max_redundant_pattern_reads=max_redundant_pattern_reads,
            pattern_end_tolerance=pattern_end_tolerance,
        )
        grouped[group_id] = kept
        redundant_dropped += dropped

    post_redundancy_ids = [rid for group_reads in grouped.values() for rid in group_reads]
    post_redundancy_total = len(post_redundancy_ids)

    if max_reads <= 0 or post_redundancy_total <= max_reads:
        selected_alignments = {read_id: alignments[read_id] for read_id in post_redundancy_ids}
        selected_counts = _count_reads_by_group(post_redundancy_ids, read_to_iso)
        return selected_alignments, full_counts, selected_counts, redundant_dropped, post_redundancy_total

    group_order = sorted(
        grouped.keys(),
        key=lambda gid: (gid == UNASSIGNED_GROUP_ID, gid),
    )
    next_idx = {gid: 0 for gid in group_order}
    active_groups = [gid for gid in group_order if grouped[gid]]
    selected_ids: List[str] = []

    while active_groups and len(selected_ids) < max_reads:
        next_active: List[str] = []
        for gid in active_groups:
            idx = next_idx[gid]
            if idx >= len(grouped[gid]):
                continue
            selected_ids.append(grouped[gid][idx])
            next_idx[gid] = idx + 1
            if next_idx[gid] < len(grouped[gid]):
                next_active.append(gid)
            if len(selected_ids) >= max_reads:
                break
        active_groups = next_active

    selected_alignments = {read_id: alignments[read_id] for read_id in selected_ids}
    selected_counts = _count_reads_by_group(selected_ids, read_to_iso)
    return selected_alignments, full_counts, selected_counts, redundant_dropped, post_redundancy_total


def draw_alignment(
    ax,
    aln: pysam.AlignedSegment,
    y: float,
    color: str,
    alpha: float = 0.92,
    height: float = 0.70,
):
    blocks = aln.get_blocks()
    if not blocks:
        return

    ax.add_patch(
        mplpatches.Rectangle(
            (blocks[0][0], y - 0.06),
            blocks[-1][1] - blocks[0][0],
            0.12,
            facecolor=color,
            edgecolor="none",
            alpha=alpha * 0.9,
        )
    )

    for start, end in blocks:
        ax.add_patch(
            mplpatches.Rectangle(
                (start, y - height / 2.0),
                max(1, end - start),
                height,
                facecolor=color,
                edgecolor="none",
                alpha=alpha,
            )
        )


def draw_isoform_model(
    ax,
    blocks: Sequence[Tuple[int, int]],
    y: float,
    color: str,
    strand: str,
    region_start: int,
    region_end: int,
    exon_height: float = READ_BLOCK_HEIGHT * MODEL_BLOCK_MIN_READ_RATIO,
    backbone_height: float = 0.10,
):
    if not blocks:
        return
    tx_start = blocks[0][0]
    tx_end = blocks[-1][1]

    ax.add_patch(
        mplpatches.Rectangle(
            (tx_start, y - backbone_height / 2.0),
            tx_end - tx_start,
            backbone_height,
            facecolor=color,
            edgecolor="none",
            alpha=0.85,
        )
    )
    for start, end in blocks:
        ax.add_patch(
            mplpatches.Rectangle(
                (start, y - exon_height / 2.0),
                max(1, end - start),
                exon_height,
                facecolor=color,
                edgecolor="none",
                alpha=0.95,
            )
        )

    # Boundary flags: model extends beyond strict plotting window.
    clipped_left = tx_start < region_start
    clipped_right = tx_end > region_end
    tri_w = max(10, int((region_end - region_start) * 0.0025))
    tri_half_h = max(0.10, min(0.22, exon_height * 0.18))
    if clipped_left:
        tri_left = [
            (region_start, y),
            (region_start + tri_w, y + tri_half_h),
            (region_start + tri_w, y - tri_half_h),
        ]
        ax.add_patch(mplpatches.Polygon(tri_left, closed=True, facecolor=color, edgecolor="none", alpha=0.9))
    if clipped_right:
        tri_right = [
            (region_end, y),
            (region_end - tri_w, y + tri_half_h),
            (region_end - tri_w, y - tri_half_h),
        ]
        ax.add_patch(mplpatches.Polygon(tri_right, closed=True, facecolor=color, edgecolor="none", alpha=0.9))
    vis_start = max(region_start, tx_start)
    vis_end = min(region_end, tx_end)

    if not clipped_left and not clipped_right:
        draw_endpoint_direction_arrows(
            ax,
            tx_start=vis_start,
            tx_end=vis_end,
            y=y,
            strand=strand,
            color="#FFFFFF",
        )
    else:
        draw_direction_arrows(
            ax,
            tx_start=vis_start,
            tx_end=vis_end,
            y=y,
            strand=strand,
            color=color,
        )


def draw_endpoint_direction_arrows(
    ax,
    tx_start: int,
    tx_end: int,
    y: float,
    strand: str,
    color: str = "#FFFFFF",
):
    """Minimal endpoint strand arrows for fully captured models."""
    width = tx_end - tx_start
    if width < 100:
        return

    pad = max(16.0, min(120.0, width * 0.08))
    arrow_len = max(18.0, min(80.0, width * 0.10))
    anchors = [tx_start + pad, tx_end - pad]

    for anchor in anchors:
        if strand == "-":
            x0, x1 = anchor + arrow_len / 2.0, anchor - arrow_len / 2.0
        else:
            x0, x1 = anchor - arrow_len / 2.0, anchor + arrow_len / 2.0
        ax.annotate(
            "",
            xy=(x1, y),
            xytext=(x0, y),
            arrowprops=dict(arrowstyle="-|>", color=color, lw=0.9, alpha=0.95, mutation_scale=7),
            zorder=6,
        )


def draw_direction_arrows(
    ax,
    tx_start: int,
    tx_end: int,
    y: float,
    strand: str,
    color: str,
):
    """Draw small strand-direction arrows on transcript models."""
    width = tx_end - tx_start
    if width < 80:
        return

    n_arrows = max(1, min(4, width // 3000 + 1))
    spacing = width / float(n_arrows + 1)
    arrow_len = max(24.0, min(120.0, width * 0.10))

    for idx in range(1, n_arrows + 1):
        center = tx_start + idx * spacing
        if strand == "-":
            x0, x1 = center + arrow_len / 2.0, center - arrow_len / 2.0
        else:
            x0, x1 = center - arrow_len / 2.0, center + arrow_len / 2.0
        ax.annotate(
            "",
            xy=(x1, y),
            xytext=(x0, y),
            arrowprops=dict(arrowstyle="-|>", color=color, lw=0.9, alpha=0.92, mutation_scale=8),
            zorder=5,
        )


def add_vertical_peak_highlights(
    axes: Iterable,
    peaks: Sequence[Tuple[int, int]],
    color: str,
    alpha_span: float = 0.10,
    alpha_line: float = 0.22,
):
    for start, end in peaks:
        width = end - start
        for ax in axes:
            if width <= 1:
                ax.axvline(start, color=color, linewidth=1.0, alpha=alpha_line, zorder=0)
            else:
                ax.axvspan(start, end, color=color, alpha=alpha_span, zorder=0)


def _accumulate_signal_interval(
    bin_sums: List[float],
    bin_nonzero_bases: List[int],
    region_start: int,
    region_end: int,
    n_bins: int,
    interval_start: int,
    interval_end: int,
    value: float,
) -> None:
    ov_start = max(region_start, interval_start)
    ov_end = min(region_end, interval_end)
    if ov_end <= ov_start:
        return
    span = region_end - region_start
    if span <= 0 or n_bins <= 0:
        return
    bin_w = span / float(n_bins)
    left_bin = max(0, int((ov_start - region_start) / bin_w))
    right_bin = min(n_bins - 1, int((ov_end - 1 - region_start) / bin_w))
    for i in range(left_bin, right_bin + 1):
        b_start = region_start + int(i * bin_w)
        b_end = region_start + int((i + 1) * bin_w) if i < n_bins - 1 else region_end
        if b_end <= b_start:
            continue
        seg_start = max(ov_start, b_start)
        seg_end = min(ov_end, b_end)
        if seg_end <= seg_start:
            continue
        seg_len = seg_end - seg_start
        bin_sums[i] += value * seg_len
        if value > 0:
            bin_nonzero_bases[i] += seg_len


def compute_binned_signal_from_bedgraphs(
    bedgraph_paths: Sequence[Path],
    chrom: str,
    region_start: int,
    region_end: int,
    n_bins: int,
) -> Tuple[List[float], List[float], int]:
    if region_end <= region_start:
        return [], [], 0
    valid_paths = [p for p in bedgraph_paths if p and p.exists()]
    if not valid_paths:
        return [], [], 0

    n_bins = max(1, n_bins)
    span = region_end - region_start
    bin_sums = [0.0] * n_bins
    bin_nonzero = [0] * n_bins

    for bg_path in valid_paths:
        with bg_path.open() as handle:
            for line in handle:
                text = line.strip()
                if not text or text.startswith("#") or text.startswith("track"):
                    continue
                fields = text.split("\t")
                if len(fields) < 4 or fields[0] != chrom:
                    continue
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                    value = float(fields[3])
                except ValueError:
                    continue
                if end <= region_start:
                    continue
                if start >= region_end:
                    # bedGraph is generally sorted; safe early break by chrom.
                    break
                _accumulate_signal_interval(
                    bin_sums=bin_sums,
                    bin_nonzero_bases=bin_nonzero,
                    region_start=region_start,
                    region_end=region_end,
                    n_bins=n_bins,
                    interval_start=start,
                    interval_end=end,
                    value=value,
                )

    # Convert weighted sums back to mean signal per bin.
    ys: List[float] = []
    xs: List[float] = []
    nonzero_bins = 0
    for i in range(n_bins):
        b_start = region_start + int(i * (span / float(n_bins)))
        b_end = region_start + int((i + 1) * (span / float(n_bins))) if i < n_bins - 1 else region_end
        b_len = max(1, b_end - b_start)
        mean_val = bin_sums[i] / float(b_len)
        ys.append(mean_val)
        xs.append((b_start + b_end) / 2.0)
        if mean_val > 0:
            nonzero_bins += 1

    return xs, ys, nonzero_bins


def build_isoform_colors(iso_order: Sequence[str]) -> Dict[str, str]:
    palette = [
        "#F94144",
        "#F3722C",
        "#F8961E",
        "#F9C74F",
        "#90BE6D",
        "#43AA8B",
        "#4D96FF",
        "#5E60CE",
        "#FF66C4",
        "#00C2A8",
        "#B8DE6F",
        "#FF7B9C",
    ]
    color_map: Dict[str, str] = {}
    for idx, iso_id in enumerate(iso_order):
        if idx < len(palette):
            color_map[iso_id] = palette[idx]
            continue
        hue = (idx * 0.61803398875) % 1.0
        rgb = colorsys.hsv_to_rgb(hue, 0.72, 0.95)
        color_map[iso_id] = mcolors.rgb2hex(rgb)
    return color_map


def main():
    parser = argparse.ArgumentParser(description="Plot isoforms and orthogonal signal tracks")
    parser.add_argument("--bam", required=True, help="Main long-read BAM")
    parser.add_argument("--readmap", required=True, help="Isoform read map file")
    parser.add_argument("--isoforms", required=True, help="Isoforms BED12 file")
    parser.add_argument("--output", required=True, help="Output file prefix")
    parser.add_argument("--region", default=None, help="Strict plotting region, format chr:start-end")

    # Peaks (support legacy arg names and explicit *-peaks names)
    parser.add_argument("--cage", dest="cage_peaks", default=None, help="CAGE peak BED (legacy arg name)")
    parser.add_argument("--quantseq", dest="quantseq_peaks", default=None, help="QuantSeq peak BED (legacy arg name)")
    parser.add_argument("--cage-peaks", dest="cage_peaks", default=None, help="CAGE peak BED")
    parser.add_argument("--quantseq-peaks", dest="quantseq_peaks", default=None, help="QuantSeq peak BED")

    parser.add_argument("--cage-signal-plus", default=None, help="Optional CAGE plus-strand bedGraph signal")
    parser.add_argument("--cage-signal-minus", default=None, help="Optional CAGE minus-strand bedGraph signal")
    parser.add_argument("--quantseq-signal-plus", default=None, help="Optional QuantSeq plus-strand bedGraph signal")
    parser.add_argument("--quantseq-signal-minus", default=None, help="Optional QuantSeq minus-strand bedGraph signal")
    parser.add_argument(
        "--max-main-reads",
        type=int,
        default=3500,
        help="Max main long-read alignments to render after balanced isoform/unassigned sampling (<=0 means all)",
    )
    parser.add_argument(
        "--max-redundant-pattern-reads",
        type=int,
        default=50,
        help="Within each isoform/unassigned group, cap near-identical chain+end patterns to this many reads (<=0 disables)",
    )
    parser.add_argument(
        "--pattern-end-tolerance",
        type=int,
        default=1,
        help="Start/end tolerance (bp) for near-identical pattern grouping",
    )
    parser.add_argument("--orth-bins", type=int, default=300, help="Number of bins for orthogonal pileups")
    args = parser.parse_args()

    bam_path = Path(args.bam)
    readmap_path = Path(args.readmap)
    isoforms_path = Path(args.isoforms)
    out_prefix = Path(args.output)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    print("[plot] Loading isoform assignments and models")
    iso_to_reads, read_to_iso = load_isoform_read_map(readmap_path)
    isoforms = load_isoforms_bed(isoforms_path)
    if not isoforms:
        print("[plot] Warning: no isoforms loaded; continuing with read-only plotting")

    strict_region = parse_region(args.region) if args.region else None
    chrom, region_start, region_end = pick_region(isoforms, bam_path, strict_region=strict_region)
    if strict_region is not None:
        isoforms = filter_isoforms_to_region(isoforms, chrom, region_start, region_end)
    print(f"[plot] Region: {chrom}:{region_start}-{region_end}")

    # Limit mapping dictionaries to isoforms currently in scope
    iso_ids_in_scope = set(isoforms.keys())
    if iso_ids_in_scope:
        iso_to_reads = {k: v for k, v in iso_to_reads.items() if k in iso_ids_in_scope}
        read_to_iso = {read: iso for read, iso in read_to_iso.items() if iso in iso_ids_in_scope}

    print("[plot] Loading main BAM alignments in region")
    all_main_alignments = load_primary_alignments(
        bam_path,
        chrom,
        region_start,
        region_end,
        max_reads=0,
    )
    main_alignments, full_group_counts, kept_group_counts, redundant_dropped, post_redundancy_total = select_balanced_main_alignments(
        alignments=all_main_alignments,
        read_to_iso=read_to_iso,
        max_reads=args.max_main_reads,
        max_redundant_pattern_reads=args.max_redundant_pattern_reads,
        pattern_end_tolerance=args.pattern_end_tolerance,
    )
    n_iso_groups_total = sum(1 for gid in full_group_counts if gid != UNASSIGNED_GROUP_ID)
    n_iso_groups_kept = sum(1 for gid, n in kept_group_counts.items() if gid != UNASSIGNED_GROUP_ID and n > 0)
    unassigned_total = full_group_counts.get(UNASSIGNED_GROUP_ID, 0)
    unassigned_kept = kept_group_counts.get(UNASSIGNED_GROUP_ID, 0)
    print(
        "[plot] Main long-read groups: "
        f"isoforms={n_iso_groups_kept}/{n_iso_groups_total}, "
        f"unassigned={unassigned_kept}/{unassigned_total} reads kept"
    )
    if redundant_dropped > 0:
        print(
            "[plot] Diversity cap applied: "
            f"{post_redundancy_total}/{len(all_main_alignments)} reads kept after "
            f"per-pattern cap (max={args.max_redundant_pattern_reads}, tolerance={args.pattern_end_tolerance} bp)"
        )
    if args.max_main_reads > 0 and post_redundancy_total > args.max_main_reads:
        print(
            "[plot] Balanced cap applied: "
            f"{len(main_alignments)}/{post_redundancy_total} reads "
            "(round-robin across isoform groups + unassigned)"
        )

    filtered_iso_to_reads: Dict[str, List[str]] = {}
    for iso_id, reads in iso_to_reads.items():
        keep = [r for r in reads if r in main_alignments]
        if keep:
            filtered_iso_to_reads[iso_id] = keep

    unassigned_reads = [r for r in sorted(main_alignments.keys()) if r not in read_to_iso]

    iso_order = sorted(
        [iso_id for iso_id in isoforms.keys() if iso_id in filtered_iso_to_reads or iso_id in iso_to_reads],
        key=lambda x: isoforms[x]["start"] if x in isoforms else 10**18,
    )
    if not iso_order:
        iso_order = sorted(isoforms.keys(), key=lambda x: isoforms[x]["start"])
    iso_color_map = build_isoform_colors(iso_order)

    cage_signal_paths = [Path(p) for p in (args.cage_signal_plus, args.cage_signal_minus) if p]
    quantseq_signal_paths = [Path(p) for p in (args.quantseq_signal_plus, args.quantseq_signal_minus) if p]

    cage_peaks = load_peaks(Path(args.cage_peaks), chrom, region_start, region_end) if args.cage_peaks else []
    quantseq_peaks = load_peaks(Path(args.quantseq_peaks), chrom, region_start, region_end) if args.quantseq_peaks else []

    x_cage, y_cage, n_cage = compute_binned_signal_from_bedgraphs(
        bedgraph_paths=cage_signal_paths,
        chrom=chrom,
        region_start=region_start,
        region_end=region_end,
        n_bins=args.orth_bins,
    )
    has_cage_signal = bool(x_cage) and max(y_cage) > 0
    if cage_signal_paths and not has_cage_signal:
        print("[plot] CAGE bedGraph provided but no non-zero signal in region; skipping CAGE panel")

    x_quant, y_quant, n_quant = compute_binned_signal_from_bedgraphs(
        bedgraph_paths=quantseq_signal_paths,
        chrom=chrom,
        region_start=region_start,
        region_end=region_end,
        n_bins=args.orth_bins,
    )
    has_quantseq_signal = bool(x_quant) and max(y_quant) > 0
    if quantseq_signal_paths and not has_quantseq_signal:
        print("[plot] QuantSeq bedGraph provided but no non-zero signal in region; skipping QuantSeq panel")
    if not has_cage_signal and not has_quantseq_signal:
        print("[plot] No orthogonal signal panels for this region; plotting models + long-read alignments only")

    n_models = max(1, len(iso_order))
    n_main_reads = max(1, len(main_alignments))
    orth_panel_count = int(has_cage_signal) + int(has_quantseq_signal)
    fig_h = max(6.5, min(22.0, 2.6 + 0.20 * n_models + 0.05 * n_main_reads + 1.2 * orth_panel_count))

    panel_layout: List[Tuple[str, float]] = [("models", 1.2)]
    if has_cage_signal:
        panel_layout.append(("cage", 1.0))
    panel_layout.append(("reads", 2.9))
    if has_quantseq_signal:
        panel_layout.append(("quant", 1.0))

    fig = plt.figure(figsize=(7.2, fig_h))
    gs = fig.add_gridspec(len(panel_layout), 1, height_ratios=[h for _, h in panel_layout], hspace=0.10)
    axes_by_name: Dict[str, plt.Axes] = {}
    previous_ax = None
    for idx, (panel_name, _) in enumerate(panel_layout):
        ax = fig.add_subplot(gs[idx, 0], sharex=previous_ax) if previous_ax is not None else fig.add_subplot(gs[idx, 0])
        axes_by_name[panel_name] = ax
        previous_ax = ax

    ax_models = axes_by_name["models"]
    ax_reads = axes_by_name["reads"]
    ax_cage = axes_by_name.get("cage")
    ax_quant = axes_by_name.get("quant")
    all_axes = list(axes_by_name.values())
    add_vertical_peak_highlights(all_axes, cage_peaks, color="#d62728")
    add_vertical_peak_highlights(all_axes, quantseq_peaks, color="#1f77b4")

    # Panel 1: transcript models
    read_block_height = READ_BLOCK_HEIGHT
    model_min_h = MODEL_BLOCK_MIN_READ_RATIO * read_block_height
    model_target_h = MODEL_BLOCK_TARGET_READ_RATIO * read_block_height
    model_max_h = MODEL_BLOCK_MAX_READ_RATIO * read_block_height
    model_block_height = min(max(model_target_h, model_min_h), model_max_h)
    model_backbone_height = max(0.06, min(model_block_height * 0.10, 0.14))
    model_row_step = model_block_height * (1.0 + MODEL_ROW_GAP_FRACTION)
    first_model_y = model_block_height / 2.0 + 0.25
    y = first_model_y
    n_models_drawn = 0
    for iso_id in iso_order:
        info = isoforms.get(iso_id)
        if not info:
            continue
        draw_isoform_model(
            ax_models,
            info["blocks"],
            y,
            iso_color_map.get(iso_id, "#555555"),
            strand=info.get("strand", "+"),
            region_start=region_start,
            region_end=region_end,
            exon_height=model_block_height,
            backbone_height=model_backbone_height,
        )
        y += model_row_step
        n_models_drawn += 1
    if n_models_drawn == 0:
        ax_models.text(0.5, 0.5, "No isoform models in region", transform=ax_models.transAxes,
                       ha="center", va="center", fontsize=7)
        model_upper = first_model_y + 0.8
    else:
        last_center = first_model_y + (n_models_drawn - 1) * model_row_step
        model_upper = last_center + model_block_height / 2.0 + 0.35
    model_top = max(model_upper, -0.8 + MODEL_PANEL_MIN_SPAN)
    ax_models.set_ylim(-0.8, model_top)
    ax_models.set_yticks([])
    ax_models.set_title("Transcript Models", fontsize=8, fontweight="normal", loc="left")

    # Optional CAGE signal panel
    if ax_cage is not None:
        ax_cage.fill_between(x_cage, y_cage, color="#d62728", alpha=0.55, linewidth=0)
        ax_cage.plot(x_cage, y_cage, color="#d62728", linewidth=0.9)
        ax_cage.set_ylim(0, max(y_cage) * 1.12)
        ax_cage.set_ylabel("signal", fontsize=8)
        ax_cage.set_title(f"CAGE (signal, nonzero bins={n_cage})", fontsize=7, fontweight="normal", loc="left")

    # Panel 3: long-read alignments (isoform color matched)
    read_y = 0.0
    if unassigned_reads:
        for read_id in unassigned_reads:
            aln = main_alignments.get(read_id)
            if aln is None:
                continue
            draw_alignment(ax_reads, aln, read_y, color="#D0D0D0", alpha=0.88, height=read_block_height)
            read_y += 1.0
        read_y += 0.55

    for iso_id in iso_order:
        read_ids = filtered_iso_to_reads.get(iso_id, [])
        if not read_ids:
            continue
        col = iso_color_map.get(iso_id, "#777777")
        for read_id in read_ids:
            aln = main_alignments.get(read_id)
            if aln is None:
                continue
            draw_alignment(ax_reads, aln, read_y, color=col, alpha=0.92, height=read_block_height)
            read_y += 1.0
        read_y += 0.55
    if read_y == 0:
        ax_reads.text(0.5, 0.5, "No long-read alignments in region", transform=ax_reads.transAxes,
                      ha="center", va="center", fontsize=7)
        read_y = 1.0
    ax_reads.set_ylim(-0.8, read_y + 0.6)
    ax_reads.set_yticks([])
    ax_reads.set_title("Long Read Alignments (Isoform-colored)", fontsize=8, fontweight="normal", loc="left")

    # Optional QuantSeq signal panel
    if ax_quant is not None:
        ax_quant.fill_between(x_quant, y_quant, color="#1f77b4", alpha=0.55, linewidth=0)
        ax_quant.plot(x_quant, y_quant, color="#1f77b4", linewidth=0.9)
        ax_quant.set_ylim(0, max(y_quant) * 1.12)
        ax_quant.set_ylabel("signal", fontsize=8)
        ax_quant.set_title(f"Quantseq (signal, nonzero bins={n_quant})", fontsize=7, fontweight="normal", loc="left")

    for ax in all_axes:
        ax.set_xlim(region_start, region_end)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
    for ax in all_axes[:-1]:
        ax.tick_params(axis="x", labelbottom=False)
    all_axes[-1].set_xlabel(f"{chrom}:{region_start:,}-{region_end:,}")

    out_png = str(out_prefix) + ".png"
    print(f"[plot] Saving {out_png}")
    fig.savefig(out_png, dpi=2200, bbox_inches="tight")
    plt.close(fig)

    print("[plot] Done")


if __name__ == "__main__":
    main()
