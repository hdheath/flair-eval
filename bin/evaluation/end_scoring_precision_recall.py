#!/usr/bin/env python3
"""
end_scoring_precision_recall.py — Precision/recall stratified by SQANTI
category, with F1 boxplots, transcripts-per-gene, and end redundancy.

Three figure outputs:
  1. precision_recall_f1_tpg.png    — 4-panel (like aim-3): P/R scatter TSS,
                                     P/R scatter TTS, F1 boxplot, TPG boxplot.
  2. precision_recall_by_sqanti.png — P/R scatter faceted by SQANTI category
                                     (FSM, ISM, NIC, NNC, SEM, SEN).
  3. end_redundancy.png             — Stacked bar of peak-deduplicated vs
                                     redundant end-site calls per junction chain.

Input:
  --evaluation-tsv: merged evaluation TSV (from synthesize_evaluations.py)
                    Must have 5prime_precision, 5prime_recall, 3prime_precision,
                    3prime_recall, 5prime_f1, 3prime_f1, isoforms_per_gene_mean,
                    FSM, ISM, NIC, NNC, SEM, SEN, transcriptome_mode, dataset.
  --isoforms-bed:   FLAIR isoforms BED12 (for per-transcript classification)
  --gtf:            Reference annotation GTF
  --peaks-5prime:   Experimental 5' peaks BED6 (e.g. CAGE; optional)
  --peaks-3prime:   Experimental 3' peaks BED6 (e.g. dRNA, dRNA; optional)
  --outdir:         Output directory


All figures follow the pub_style conventions:
  - Okabe-Ito colorblind-safe palette
  - DPI 300, no top/right spines
  - legend_outside, savefig with tight bbox
"""

import argparse
import csv
import logging
import sys
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np

# Import project style + structural classification
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pub_style import (
    PALETTE,
    ASSEMBLER_COLORS,
    LIBRARY_SHAPES,
    MODE_COLORS,
    apply_rc,
    style_ax,
    legend_outside,
    savefig,
)
from flair_structural import (
    parse_gtf_transcripts,
    build_reference_structures,
    classify_transcripts_per_isoform,
)

apply_rc()

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s")
log = logging.getLogger(__name__)


# ── Colour scheme ───────────────────────────────────────────────────────────

CATEGORY_COLORS = {
    "FSM": "#2ecc71",   # green
    "ISM": "#3498db",   # blue
    "NIC": "#f39c12",   # yellow/orange
    "NNC": "#e74c3c",   # red
    "SEM": "#9b59b6",   # purple
    "SEN": "#95a5a6",   # grey
}
CATEGORY_ORDER = ["FSM", "ISM", "NIC", "NNC", "SEM", "SEN"]


# ── Utilities ───────────────────────────────────────────────────────────────

def _safe_f1(precision, recall):
    if precision is None or recall is None:
        return None
    if precision + recall == 0:
        return 0.0
    return 2 * precision * recall / (precision + recall)


def _nearest_distance_abs(query: int, sorted_positions: List[int]) -> int:
    """Absolute distance to nearest position in a sorted list."""
    if not sorted_positions:
        return 999999
    idx = bisect_left(sorted_positions, query)
    best = 999999
    for i in (idx - 1, idx):
        if 0 <= i < len(sorted_positions):
            d = abs(query - sorted_positions[i])
            if d < best:
                best = d
    return best


def parse_peaks_bed(path: str) -> Dict[Tuple[str, str], List[int]]:
    """Parse BED6 peaks into {(chrom, strand): sorted positions}."""
    peaks: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            chrom, start, end, strand = cols[0], int(cols[1]), int(cols[2]), cols[5]
            mid = (start + end) // 2
            peaks[(chrom, strand)].append(mid)
    for k in peaks:
        peaks[k].sort()
    return peaks


def parse_gtf_ends(gtf_path: str) -> Dict[str, Dict[str, List[int]]]:
    """Extract annotated TSS/TTS positions from GTF, keyed by chrom."""
    ends: Dict[str, Dict[str, set]] = defaultdict(lambda: {"tss": set(), "tts": set()})
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            chrom = cols[0]
            start = int(cols[3]) - 1
            end = int(cols[4])
            strand = cols[6]
            if strand == "+":
                ends[chrom]["tss"].add(start)
                ends[chrom]["tts"].add(end)
            else:
                ends[chrom]["tss"].add(end)
                ends[chrom]["tts"].add(start)
    return {c: {"tss": sorted(ends[c]["tss"]), "tts": sorted(ends[c]["tts"])}
            for c in ends}


# ── Per-Category Precision/Recall ───────────────────────────────────────────

def compute_per_category_precision_recall(
    classified_isoforms: List[dict],
    annotated_ends: Dict[str, Dict[str, List[int]]],
    window: int = 50,
) -> Dict[str, dict]:
    """
    For each SQANTI category, compute 5' and 3' precision/recall.

    Precision = fraction of isoform ends within `window` bp of an annotation.
    Recall = fraction of annotated ends matched by at least one isoform.

    Returns:
        {category: {
            "n": count,
            "5prime_precision": float, "5prime_recall": float, "5prime_f1": float,
            "3prime_precision": float, "3prime_recall": float, "3prime_f1": float,
        }}
    """
    # Group isoforms by category
    by_cat: Dict[str, List[dict]] = defaultdict(list)
    for iso in classified_isoforms:
        by_cat[iso["category"]].append(iso)

    results = {}
    for cat in CATEGORY_ORDER:
        isos = by_cat.get(cat, [])
        n = len(isos)
        if n == 0:
            results[cat] = {
                "n": 0,
                "5prime_precision": None, "5prime_recall": None, "5prime_f1": None,
                "3prime_precision": None, "3prime_recall": None, "3prime_f1": None,
            }
            continue

        for end_type, end_label in [("tss", "5prime"), ("tts", "3prime")]:
            matched_isos = 0
            matched_annots = set()
            for iso in isos:
                chrom = iso["chrom"]
                strand = iso["strand"]
                if strand == "+":
                    pos = iso["start"] if end_type == "tss" else iso["end"]
                else:
                    pos = iso["end"] if end_type == "tss" else iso["start"]

                chrom_ends = annotated_ends.get(chrom, {}).get(end_type, [])
                dist = _nearest_distance_abs(pos, chrom_ends)
                if dist <= window:
                    matched_isos += 1
                    # Find which annotation was matched
                    idx = bisect_left(chrom_ends, pos)
                    for i in (idx - 1, idx):
                        if 0 <= i < len(chrom_ends) and abs(pos - chrom_ends[i]) <= window:
                            matched_annots.add((chrom, chrom_ends[i]))

            precision = matched_isos / n if n > 0 else None

            # For recall: how many unique annotated ends were hit
            all_annots = set()
            for chrom in annotated_ends:
                for pos in annotated_ends[chrom].get(end_type, []):
                    all_annots.add((chrom, pos))
            recall = len(matched_annots) / len(all_annots) if all_annots else None

            f1 = _safe_f1(precision, recall)

            if cat not in results:
                results[cat] = {"n": n}
            results[cat][f"{end_label}_precision"] = precision
            results[cat][f"{end_label}_recall"] = recall
            results[cat][f"{end_label}_f1"] = f1

    return results


# ── End Redundancy Metric ──────────────────────────────────────────────────

def compute_end_redundancy(
    classified_isoforms: List[dict],
    isoforms_bed_path: str,
    peaks_5prime: Optional[Dict[Tuple[str, str], List[int]]] = None,
    peaks_3prime: Optional[Dict[Tuple[str, str], List[int]]] = None,
    annotated_ends: Optional[Dict[str, Dict[str, List[int]]]] = None,
    window: int = 50,
) -> dict:
    """
    Compute end redundancy: multiple isoform ends calling the same peak/annotation.

    For each junction chain group (isoforms sharing all splice junctions):
      - Count unique TSS/TTS positions
      - Count unique peaks (or annotation sites) those positions map to
      - Redundant = unique_positions - unique_peaks

    Returns:
        {
            "groups_total": int,
            "groups_with_redundancy": int,
            "redundancy_rate": float,
            "tss_unique_positions": int,
            "tss_unique_peaks": int,
            "tss_redundant": int,
            "tts_unique_positions": int,
            "tts_unique_peaks": int,
            "tts_redundant": int,
            "per_group": [  # for plotting
                {"jc_key": str, "n_isoforms": int,
                 "tss_unique": int, "tss_peaks": int, "tss_redundant": int,
                 "tts_unique": int, "tts_peaks": int, "tts_redundant": int,
                 "category_counts": {cat: int}},
            ]
        }
    """
    # Parse BED12 to get junction chains
    isoforms_by_jc: Dict[Tuple[str, str, tuple], List[dict]] = defaultdict(list)
    iso_lookup = {iso["name"]: iso for iso in classified_isoforms}

    with open(isoforms_bed_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            chrom, strand = parts[0], parts[5]
            name = parts[3]
            start, end = int(parts[1]), int(parts[2])
            esizes = [int(x) for x in parts[10].rstrip(",").split(",")]
            estarts = [int(x) for x in parts[11].rstrip(",").split(",")]
            exons = [(start + estarts[i], start + estarts[i] + esizes[i])
                     for i in range(len(esizes))]
            introns = tuple((exons[x][1], exons[x + 1][0])
                            for x in range(len(exons) - 1))
            if not introns:
                continue  # Skip single-exon for junction chain analysis

            if strand == "+":
                tss, tts = start, end
            else:
                tss, tts = end, start

            isoforms_by_jc[(chrom, strand, introns)].append({
                "name": name,
                "tss": tss,
                "tts": tts,
                "category": iso_lookup.get(name, {}).get("category", "unknown"),
            })

    def _find_peak(pos, peaks_dict, chrom, strand, window):
        """Map a position to its nearest peak."""
        positions = peaks_dict.get((chrom, strand), [])
        if not positions:
            return None
        idx = bisect_left(positions, pos)
        best = None
        for i in (idx - 1, idx):
            if 0 <= i < len(positions) and abs(pos - positions[i]) <= window:
                if best is None or abs(pos - positions[i]) < abs(pos - best):
                    best = positions[i]
        return best

    def _find_annot(pos, annot_ends, chrom, end_type, window):
        """Map a position to nearest annotation site."""
        positions = annot_ends.get(chrom, {}).get(end_type, [])
        if not positions:
            return None
        idx = bisect_left(positions, pos)
        best = None
        for i in (idx - 1, idx):
            if 0 <= i < len(positions) and abs(pos - positions[i]) <= window:
                if best is None or abs(pos - positions[i]) < abs(pos - best):
                    best = positions[i]
        return best

    total_groups = 0
    groups_with_redundancy = 0
    total_tss_unique = 0
    total_tss_peaks = 0
    total_tts_unique = 0
    total_tts_peaks = 0
    per_group = []

    for (chrom, strand, jc), isos in isoforms_by_jc.items():
        if len(isos) < 2:
            continue
        total_groups += 1

        # Category distribution within this group
        cat_counts = defaultdict(int)
        for iso in isos:
            cat_counts[iso["category"]] += 1

        # Unique TSS/TTS positions
        tss_positions = sorted(set(iso["tss"] for iso in isos))
        tts_positions = sorted(set(iso["tts"] for iso in isos))
        n_tss_unique = len(tss_positions)
        n_tts_unique = len(tts_positions)

        # Map to peaks (or annotations)
        tss_peaks = set()
        for pos in tss_positions:
            if peaks_5prime:
                pk = _find_peak(pos, peaks_5prime, chrom, strand, window)
            elif annotated_ends:
                pk = _find_annot(pos, annotated_ends, chrom, "tss", window)
            else:
                pk = None
            if pk is not None:
                tss_peaks.add(pk)

        tts_peaks = set()
        for pos in tts_positions:
            if peaks_3prime:
                pk = _find_peak(pos, peaks_3prime, chrom, strand, window)
            elif annotated_ends:
                pk = _find_annot(pos, annotated_ends, chrom, "tts", window)
            else:
                pk = None
            if pk is not None:
                tts_peaks.add(pk)

        n_tss_peaks = len(tss_peaks)
        n_tts_peaks = len(tts_peaks)
        tss_redundant = max(0, n_tss_unique - max(1, n_tss_peaks))
        tts_redundant = max(0, n_tts_unique - max(1, n_tts_peaks))

        if tss_redundant > 0 or tts_redundant > 0:
            groups_with_redundancy += 1

        total_tss_unique += n_tss_unique
        total_tss_peaks += n_tss_peaks
        total_tts_unique += n_tts_unique
        total_tts_peaks += n_tts_peaks

        per_group.append({
            "jc_key": f"{chrom}:{strand}:{hash(jc) & 0xFFFF:04x}",
            "n_isoforms": len(isos),
            "tss_unique": n_tss_unique,
            "tss_peaks": n_tss_peaks,
            "tss_redundant": tss_redundant,
            "tts_unique": n_tts_unique,
            "tts_peaks": n_tts_peaks,
            "tts_redundant": tts_redundant,
            "category_counts": dict(cat_counts),
        })

    return {
        "groups_total": total_groups,
        "groups_with_redundancy": groups_with_redundancy,
        "redundancy_rate": groups_with_redundancy / max(1, total_groups),
        "tss_unique_positions": total_tss_unique,
        "tss_unique_peaks": total_tss_peaks,
        "tss_redundant": max(0, total_tss_unique - total_tss_peaks),
        "tts_unique_positions": total_tts_unique,
        "tts_unique_peaks": total_tts_peaks,
        "tts_redundant": max(0, total_tts_unique - total_tts_peaks),
        "per_group": per_group,
    }


# ═══════════════════════════════════════════════════════════════════════════
# PLOTS
# ═══════════════════════════════════════════════════════════════════════════

# ── Plot 1: Precision/Recall + F1 + TPG (aim-3 style) ──────────────────────

def plot_precision_recall_f1_tpg(
    eval_rows: List[dict],
    outdir: Path,
    mode_key: str = "transcriptome_mode",
    prefix: str = "",
):
    """
    4-panel figure (like aim-3_precision-recall_v.1.py):
      TL: 5' Precision vs Recall scatter
      TR: 3' Precision vs Recall scatter
      BL: Mean F1 boxplot by mode
      BR: Isoforms-per-gene boxplot by mode
    """
    # Filter rows with valid data
    valid = [r for r in eval_rows
             if r.get("5prime_precision") and r.get("5prime_recall")]
    if not valid:
        log.warning("No valid precision/recall data for plot_precision_recall_f1_tpg")
        return

    fig, axes = plt.subplots(2, 2, figsize=(10, 9), constrained_layout=True)

    # Discover modes
    modes = sorted(set(r.get(mode_key, "unknown") for r in valid))
    mode_colors = {}
    for i, m in enumerate(modes):
        mode_colors[m] = MODE_COLORS.get(m, PALETTE[i % len(PALETTE)])

    # ── TL: 5' P/R scatter ──
    ax = axes[0, 0]
    for r in valid:
        m = r.get(mode_key, "unknown")
        p5 = float(r["5prime_precision"]) * 100
        r5 = float(r["5prime_recall"]) * 100
        ax.scatter(r5, p5, s=30, c=mode_colors[m], edgecolors="black",
                   linewidth=0.4, alpha=0.7, zorder=2)
    style_ax(ax, xlabel="5′ Recall (%)", ylabel="5′ Precision (%)",
             title="TSS Precision vs Recall")
    ax.set_xlim(0, 105)
    ax.set_ylim(0, 105)
    ax.plot([0, 100], [0, 100], "--", color="#cccccc", linewidth=0.8, zorder=0)

    # ── TR: 3' P/R scatter ──
    ax = axes[0, 1]
    for r in valid:
        m = r.get(mode_key, "unknown")
        p3 = float(r["3prime_precision"]) * 100
        r3 = float(r["3prime_recall"]) * 100
        ax.scatter(r3, p3, s=30, c=mode_colors[m], edgecolors="black",
                   linewidth=0.4, alpha=0.7, zorder=2)
    style_ax(ax, xlabel="3′ Recall (%)", ylabel="3′ Precision (%)",
             title="TTS Precision vs Recall")
    ax.set_xlim(0, 105)
    ax.set_ylim(0, 105)
    ax.plot([0, 100], [0, 100], "--", color="#cccccc", linewidth=0.8, zorder=0)

    # ── BL: Mean F1 boxplot ──
    ax = axes[1, 0]
    f1_by_mode = defaultdict(list)
    for r in valid:
        m = r.get(mode_key, "unknown")
        f1_5 = float(r.get("5prime_f1", 0) or 0)
        f1_3 = float(r.get("3prime_f1", 0) or 0)
        mean_f1 = (f1_5 + f1_3) / 2 if (f1_5 and f1_3) else max(f1_5, f1_3)
        f1_by_mode[m].append(mean_f1)

    positions = list(range(len(modes)))
    bp_data = [f1_by_mode.get(m, [0]) for m in modes]
    bp = ax.boxplot(bp_data, positions=positions, widths=0.6, patch_artist=True,
                    showfliers=True, flierprops=dict(markersize=3))
    for patch, m in zip(bp["boxes"], modes):
        patch.set_facecolor(mode_colors[m])
        patch.set_alpha(0.6)
    ax.set_xticks(positions)
    ax.set_xticklabels([m.replace("-", "\n").replace("_", "\n") for m in modes],
                       fontsize=7)
    style_ax(ax, ylabel="Mean F1 Score", title="F1 by Pipeline Mode")

    # ── BR: Isoforms/gene boxplot ──
    ax = axes[1, 1]
    tpg_by_mode = defaultdict(list)
    for r in valid:
        m = r.get(mode_key, "unknown")
        tpg = float(r.get("isoforms_per_gene_mean", 0) or 0)
        if tpg > 0:
            tpg_by_mode[m].append(tpg)

    bp_data = [tpg_by_mode.get(m, [0]) for m in modes]
    bp = ax.boxplot(bp_data, positions=positions, widths=0.6, patch_artist=True,
                    showfliers=True, flierprops=dict(markersize=3))
    for patch, m in zip(bp["boxes"], modes):
        patch.set_facecolor(mode_colors[m])
        patch.set_alpha(0.6)
    ax.set_xticks(positions)
    ax.set_xticklabels([m.replace("-", "\n").replace("_", "\n") for m in modes],
                       fontsize=7)
    style_ax(ax, ylabel="Isoforms per Gene (mean)", title="Transcripts per Gene")

    # Panel labels
    for idx, (label, ax_) in enumerate(zip("ABCD", axes.flat)):
        ax_.text(-0.08, 1.05, label, transform=ax_.transAxes,
                 fontsize=8, fontweight="normal", va="top")

    # Legend
    handles = [mpatches.Patch(color=mode_colors[m], label=m.replace("_", " "),
                              alpha=0.7)
               for m in modes]
    legend_outside(axes[0, 1], handles=handles)

    savefig(fig, outdir / f"{prefix}precision_recall_f1_tpg.png")
    log.info(f"  → {prefix}precision_recall_f1_tpg.png")


# ── Plot 2: Precision/Recall by SQANTI Category ────────────────────────────

def plot_precision_recall_by_sqanti(
    per_category: Dict[str, dict],
    outdir: Path,
    library_type: Optional[str] = None,
    prefix: str = "",
):
    """
    2-panel (TSS, TTS): grouped bar chart of precision & recall per category.
    """
    cats_with_data = [c for c in CATEGORY_ORDER if per_category.get(c, {}).get("n", 0) > 0]
    if not cats_with_data:
        log.warning("No SQANTI categories with data for P/R plot")
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    if library_type:
        fig.suptitle(f"Library: {library_type}", fontsize=8, fontstyle="italic", y=1.02)

    for panel_idx, (end_label, end_title) in enumerate([
        ("5prime", "TSS (5′)"), ("3prime", "TTS (3′)")
    ]):
        ax = axes[panel_idx]
        x = np.arange(len(cats_with_data))
        width = 0.35

        precisions = []
        recalls = []
        counts = []
        for cat in cats_with_data:
            d = per_category[cat]
            p = d.get(f"{end_label}_precision")
            r = d.get(f"{end_label}_recall")
            precisions.append(100 * p if p is not None else 0)
            recalls.append(100 * r if r is not None else 0)
            counts.append(d["n"])

        bars_p = ax.bar(x - width / 2, precisions, width, label="Precision",
                        color=[CATEGORY_COLORS[c] for c in cats_with_data],
                        edgecolor="black", linewidth=0.5, alpha=0.85)
        bars_r = ax.bar(x + width / 2, recalls, width, label="Recall",
                        color=[CATEGORY_COLORS[c] for c in cats_with_data],
                        edgecolor="black", linewidth=0.5, alpha=0.45,
                        hatch="///")

        # Count annotations
        for i, (bar_p, cnt) in enumerate(zip(bars_p, counts)):
            ax.text(bar_p.get_x() + width, -4, f"n={cnt}",
                    ha="center", va="top", fontsize=7, fontstyle="italic")

        ax.set_xticks(x)
        ax.set_xticklabels(cats_with_data, fontsize=7)
        style_ax(ax,
                 xlabel="SQANTI Category",
                 ylabel=f"{end_title} %" if panel_idx == 0 else "",
                 title=f"{end_title} Precision & Recall by Category")
        ax.set_ylim(-8, 110)

        if panel_idx == 0:
            # Custom legend: filled = precision, hatched = recall
            prec_patch = mpatches.Patch(facecolor=PALETTE[7], edgecolor="black",
                                        linewidth=0.5, alpha=0.85, label="Precision")
            recall_patch = mpatches.Patch(facecolor=PALETTE[7], edgecolor="black",
                                          linewidth=0.5, alpha=0.45, hatch="///",
                                          label="Recall")
            ax.legend(handles=[prec_patch, recall_patch], fontsize=8,
                      frameon=False, loc="upper right")

    # Panel labels
    for idx, (label, ax_) in enumerate(zip("AB", axes)):
        ax_.text(-0.06, 1.05, label, transform=ax_.transAxes,
                 fontsize=8, fontweight="normal", va="top")

    savefig(fig, outdir / f"{prefix}precision_recall_by_sqanti.png")
    log.info(f"  → {prefix}precision_recall_by_sqanti.png")


# ── Plot 3: End Redundancy ─────────────────────────────────────────────────

def plot_end_redundancy(redundancy: dict, outdir: Path,
                       library_type: Optional[str] = None,
                       prefix: str = ""):
    """
    3-panel plot:
      A: Stacked bar — distinct peaks vs redundant calls (TSS & TTS)
      B: Histogram — redundant positions per junction-chain group
      C: Scatter — group size vs redundancy, colored by dominant category
    """
    per_group = redundancy["per_group"]
    if not per_group:
        log.warning("No junction-chain groups for end redundancy plot")
        return

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    if library_type:
        fig.suptitle(f"Library: {library_type}", fontsize=8, fontstyle="italic", y=1.02)

    # ── Panel A: Summary stacked bar ──
    ax = axes[0]
    labels = ["TSS (5′)", "TTS (3′)"]
    peaks = [redundancy["tss_unique_peaks"], redundancy["tts_unique_peaks"]]
    redundant = [redundancy["tss_redundant"], redundancy["tts_redundant"]]

    x = np.arange(len(labels))
    ax.bar(x, peaks, width=0.5, color=PALETTE[2], edgecolor="none",
           label="Distinct peaks", alpha=0.85)
    ax.bar(x, redundant, width=0.5, bottom=peaks, color=PALETTE[5],
           edgecolor="none", label="Redundant calls", alpha=0.85)

    for i, (p, r) in enumerate(zip(peaks, redundant)):
        total = p + r
        if total > 0:
            pct_red = 100 * r / total
            ax.text(i, total + 0.5, f"{pct_red:.0f}% redundant",
                    ha="center", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    style_ax(ax, ylabel="End-site positions",
             title=f"End Redundancy\n({redundancy['groups_total']} junction-chain groups)")
    ax.legend(fontsize=8, frameon=False)

    # ── Panel B: Redundancy histogram ──
    ax = axes[1]
    tss_red = [g["tss_redundant"] for g in per_group]
    tts_red = [g["tts_redundant"] for g in per_group]
    max_red = max(max(tss_red, default=0), max(tts_red, default=0))
    bins = np.arange(-0.5, max_red + 1.5, 1)

    ax.hist(tss_red, bins=bins, alpha=0.55, color=PALETTE[4],
            edgecolor="none", label="TSS", density=False)
    ax.hist(tts_red, bins=bins, alpha=0.55, color=PALETTE[0],
            edgecolor="none", label="TTS", density=False)
    style_ax(ax, xlabel="Redundant positions per group",
             ylabel="Junction-chain groups",
             title="Redundancy Distribution")
    ax.legend(fontsize=8, frameon=False)

    # ── Panel C: Group size vs redundancy scatter ──
    ax = axes[2]
    for g in per_group:
        total_redundant = g["tss_redundant"] + g["tts_redundant"]
        # Dominant category
        if g["category_counts"]:
            dom_cat = max(g["category_counts"], key=g["category_counts"].get)
            color = CATEGORY_COLORS.get(dom_cat, PALETTE[7])
        else:
            color = PALETTE[7]
        ax.scatter(g["n_isoforms"], total_redundant, s=15, c=color,
                   edgecolors="black", linewidth=0.3, alpha=0.6, zorder=2)

    style_ax(ax, xlabel="Isoforms per junction chain",
             ylabel="Total redundant end positions",
             title="Group Size vs Redundancy")

    # Category legend
    handles = [mpatches.Patch(color=CATEGORY_COLORS[c], label=c)
               for c in CATEGORY_ORDER
               if any(c in g.get("category_counts", {}) for g in per_group)]
    if handles:
        ax.legend(handles=handles, fontsize=7, frameon=False, loc="upper left")

    # Panel labels
    for idx, (label, ax_) in enumerate(zip("ABC", axes)):
        ax_.text(-0.06, 1.05, label, transform=ax_.transAxes,
                 fontsize=8, fontweight="normal", va="top")

    savefig(fig, outdir / f"{prefix}end_redundancy.png")
    log.info(f"  → {prefix}end_redundancy.png")


# ── Write TSV outputs ──────────────────────────────────────────────────────

def write_per_category_tsv(per_category: Dict[str, dict], outpath: Path,
                          library_type: Optional[str] = None):
    base_fields = ["category", "n",
                   "5prime_precision", "5prime_recall", "5prime_f1",
                   "3prime_precision", "3prime_recall", "3prime_f1"]
    fieldnames = (["library_type"] + base_fields) if library_type else base_fields
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for cat in CATEGORY_ORDER:
            d = per_category.get(cat, {})
            row = {"category": cat, "n": d.get("n", 0)}
            if library_type:
                row["library_type"] = library_type
            for k in base_fields[2:]:
                v = d.get(k)
                row[k] = f"{v:.4f}" if v is not None else ""
            w.writerow(row)
    log.info(f"  → {outpath.name}")


def write_redundancy_tsv(redundancy: dict, outpath: Path,
                        library_type: Optional[str] = None):
    base_fields = ["groups_total", "groups_with_redundancy", "redundancy_rate",
                   "tss_unique_positions", "tss_unique_peaks", "tss_redundant",
                   "tts_unique_positions", "tts_unique_peaks", "tts_redundant"]
    fieldnames = (["library_type"] + base_fields) if library_type else base_fields
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        row = {k: redundancy[k] for k in base_fields}
        row["redundancy_rate"] = f"{row['redundancy_rate']:.4f}"
        if library_type:
            row["library_type"] = library_type
        w.writerow(row)

    # Also write per-group detail
    detail_path = outpath.with_suffix(".detail.tsv")
    detail_fields = ["jc_key", "n_isoforms",
                     "tss_unique", "tss_peaks", "tss_redundant",
                     "tts_unique", "tts_peaks", "tts_redundant",
                     "dominant_category"]
    with open(detail_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=detail_fields, delimiter="\t")
        w.writeheader()
        for g in redundancy["per_group"]:
            dom = max(g["category_counts"], key=g["category_counts"].get) \
                if g["category_counts"] else ""
            w.writerow({
                **{k: g[k] for k in detail_fields[:8]},
                "dominant_category": dom,
            })
    log.info(f"  → {outpath.name}, {detail_path.name}")


# ── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Precision/recall by SQANTI category + end redundancy metric.",
    )
    parser.add_argument("--evaluation-tsv", default=None,
                        help="Merged evaluation TSV (for P/R/F1/TPG plot)")
    parser.add_argument("--isoforms-bed", required=True,
                        help="FLAIR isoforms BED12")
    parser.add_argument("--gtf", required=True,
                        help="Reference annotation GTF")
    parser.add_argument("--peaks-5prime", default=None,
                        help="Experimental 5' peaks BED6 (for redundancy; falls back to annotation)")
    parser.add_argument("--peaks-3prime", default=None,
                        help="Experimental 3' peaks BED6 (for redundancy; falls back to annotation)")
    parser.add_argument("--window", type=int, default=50,
                        help="Window for precision/recall and peak matching (bp)")
    parser.add_argument("--library-type", default=None,
                        help="Library type label (e.g. ont_cDNA, pacbio_isoseq). "
                             "Added as column in output TSVs and subtitle in plots.")
    parser.add_argument("--prefix", default="",
                        help="Prefix prepended to output filenames (e.g. 'A549_default')")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for plots and TSVs")

    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    pfx = f"{args.prefix}_" if args.prefix else ""

    # ── Build reference structures for SQANTI classification ──
    log.info("Parsing reference GTF...")
    transcripttexons = parse_gtf_transcripts(args.gtf)
    refjuncs, refjuncchains, refseends = build_reference_structures(transcripttexons)

    log.info("Classifying isoforms...")
    classified = classify_transcripts_per_isoform(
        args.isoforms_bed, refjuncs, refjuncchains, refseends
    )
    log.info(f"  {len(classified)} isoforms classified")

    cat_counts = defaultdict(int)
    for iso in classified:
        cat_counts[iso["category"]] += 1
    log.info(f"  Categories: {dict(cat_counts)}")

    # ── Parse annotation ends for P/R ──
    annotated_ends = parse_gtf_ends(args.gtf)

    # ── Per-category precision/recall ──
    log.info("Computing per-category precision/recall...")
    per_category = compute_per_category_precision_recall(
        classified, annotated_ends, window=args.window
    )
    lib_tag = f".{args.library_type}" if args.library_type else ""
    write_per_category_tsv(per_category,
                          outdir / f"{pfx}precision_recall_by_category{lib_tag}.tsv",
                          library_type=args.library_type)

    # ── End redundancy ──
    log.info("Computing end redundancy...")
    peaks_5 = parse_peaks_bed(args.peaks_5prime) if args.peaks_5prime else None
    peaks_3 = parse_peaks_bed(args.peaks_3prime) if args.peaks_3prime else None

    redundancy = compute_end_redundancy(
        classified_isoforms=classified,
        isoforms_bed_path=args.isoforms_bed,
        peaks_5prime=peaks_5,
        peaks_3prime=peaks_3,
        annotated_ends=annotated_ends,
        window=args.window,
    )
    write_redundancy_tsv(redundancy, outdir / f"{pfx}end_redundancy{lib_tag}.tsv",
                        library_type=args.library_type)

    # ── Plots ──
    log.info("Generating plots...")

    # Plot 1: P/R + F1 + TPG (needs evaluation TSV)
    if args.evaluation_tsv:
        eval_rows = []
        with open(args.evaluation_tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                eval_rows.append(row)
        plot_precision_recall_f1_tpg(eval_rows, outdir, prefix=pfx)

    # Plot 2: P/R by SQANTI category
    plot_precision_recall_by_sqanti(per_category, outdir,
                                    library_type=args.library_type, prefix=pfx)

    # Plot 3: End redundancy
    plot_end_redundancy(redundancy, outdir, library_type=args.library_type, prefix=pfx)

    log.info("Done.")


if __name__ == "__main__":
    main()
