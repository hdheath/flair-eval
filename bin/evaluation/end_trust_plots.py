#!/usr/bin/env python3
"""
end_trust_plots.py — Publication-quality plots for end-scoring evaluation.

Reads the per-boundary TSV produced by evaluate_end_trust.py and produces
six figure panels answering:

  1. When should we trust 5' vs 3' ends?
  2. How does library type change end trust?
  3. When do sequence features rescue boundaries?
  4. How does alpha blending affect displacement?
  5. What is the annotation vs confidence relationship?
  6. Weight decomposition: which signal dimension matters most?

All figures follow the pub_style conventions:
  - Okabe-Ito colorblind-safe palette
  - DPI 300, no top/right spines
  - legend_outside, savefig with tight bbox
"""

import argparse
import csv
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Import the project's shared pub_style
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pub_style import (
    PALETTE,
    LIBRARY_SHAPES,
    MODE_COLORS,
    apply_rc,
    style_ax,
    legend_outside,
    savefig,
)

apply_rc()

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-8s  %(message)s")
log = logging.getLogger(__name__)


# ── Colour scheme for end-scoring plots ─────────────────────────────────────

PROFILE_COLORS = {
    "default":       PALETTE[7],   # grey
    "ont_cDNA":      PALETTE[0],   # orange
    "ont_dRNA":      PALETTE[5],   # vermillion
    "pacbio_isoseq": PALETTE[4],   # blue
    "pacbio_masseq": PALETTE[2],   # green
}

END_COLORS = {
    "tss": PALETTE[4],   # blue
    "tts": PALETTE[0],   # orange
}

TRUST_COLORS = {
    "trusted":   PALETTE[2],   # green
    "uncertain": PALETTE[3],   # yellow
    "untrusted": PALETTE[5],   # vermillion
}

ALPHA_CMAP = plt.cm.viridis


# ── Data loading ────────────────────────────────────────────────────────────

def load_boundary_data(tsv_path: str) -> List[dict]:
    """Load per-boundary TSV into list of dicts."""
    rows = []
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            row["pos"] = int(row["pos"])
            row["read_depth"] = int(row["read_depth"])
            row["annot_displacement"] = int(row["annot_displacement"])
            row["alpha"] = float(row["alpha"])
            row["seq_score"] = float(row["seq_score"])
            row["depth_score"] = float(row["depth_score"])
            row["tech_penalty"] = float(row["tech_penalty"])
            row["confidence"] = float(row["confidence"])
            row["tss_trust"] = float(row["tss_trust"])
            row["tts_trust"] = float(row["tts_trust"])
            row["seq_weight"] = float(row["seq_weight"])
            row["depth_weight"] = float(row["depth_weight"])
            row["tech_weight"] = float(row["tech_weight"])
            rows.append(row)
    return rows


def load_summary_data(tsv_path: str) -> List[dict]:
    """Load aggregate summary TSV."""
    rows = []
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            for k in row:
                try:
                    row[k] = float(row[k])
                except (ValueError, TypeError):
                    pass
            rows.append(row)
    return rows


# ── Plot 1: End Trust by Library Type ───────────────────────────────────────

def plot_end_trust_by_library(rows: List[dict], outdir: Path):
    """
    Grouped bar chart: % boundaries within 50bp of annotation,
    grouped by library profile, separated by TSS vs TTS.

    Answers: "Which library type gives us the most trustworthy 5' vs 3' ends?"
    """
    # Filter to alpha=0 (raw scoring, no blending)
    data = [r for r in rows if r["alpha"] == 0.0]

    # Group by (profile, end_type)
    groups = defaultdict(list)
    for r in data:
        groups[(r["profile"], r["end_type"])].append(abs(r["annot_displacement"]))

    profiles = [p for p in PROFILE_COLORS if any(k[0] == p for k in groups)]
    if not profiles:
        log.warning("No data for plot_end_trust_by_library")
        return

    x = np.arange(len(profiles))
    width = 0.35

    fig, ax = plt.subplots(figsize=(max(6, len(profiles) * 1.8), 4.5))

    tss_pcts = []
    tts_pcts = []
    for p in profiles:
        tss_disps = groups.get((p, "tss"), [])
        tts_disps = groups.get((p, "tts"), [])
        tss_pcts.append(100 * sum(1 for d in tss_disps if d <= 50) / max(1, len(tss_disps)))
        tts_pcts.append(100 * sum(1 for d in tts_disps if d <= 50) / max(1, len(tts_disps)))

    bars_tss = ax.bar(x - width / 2, tss_pcts, width, label="TSS (5′)",
                      color=END_COLORS["tss"], edgecolor="none", alpha=0.85)
    bars_tts = ax.bar(x + width / 2, tts_pcts, width, label="TTS (3′)",
                      color=END_COLORS["tts"], edgecolor="none", alpha=0.85)

    # Value labels
    for bar in list(bars_tss) + list(bars_tts):
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 1, f"{h:.0f}%",
                ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels([p.replace("_", "\n") for p in profiles], fontsize=7)
    style_ax(ax, ylabel="Boundaries within 50 bp\nof annotation (%)",
             title="End trust by library type")
    ax.set_ylim(0, 105)
    legend_outside(ax)

    savefig(fig, outdir / "end_trust_by_library.png")
    log.info("  → end_trust_by_library.png")


# ── Plot 2: Displacement Distributions by End Type ──────────────────────────

def plot_displacement_distributions(rows: List[dict], outdir: Path):
    """
    Overlaid density histograms of annotation displacement for TSS vs TTS
    across all profiles.  One panel per library profile (2×3 grid).

    Answers: "How accurate are raw boundaries, and which end is worse?"
    """
    data = [r for r in rows if r["alpha"] == 0.0]

    profiles = sorted(set(r["profile"] for r in data))
    n = len(profiles)
    ncols = min(3, n)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.5 * nrows),
                             squeeze=False, constrained_layout=True)

    for idx, profile in enumerate(profiles):
        ax = axes[idx // ncols][idx % ncols]
        for end_type in ["tss", "tts"]:
            disps = [r["annot_displacement"] for r in data
                     if r["profile"] == profile and r["end_type"] == end_type]
            if not disps:
                continue
            # Clip for display
            clipped = [max(-500, min(500, d)) for d in disps]
            ax.hist(clipped, bins=50, alpha=0.55, density=True,
                    color=END_COLORS[end_type], edgecolor="none",
                    label=f"{'TSS (5′)' if end_type == 'tss' else 'TTS (3′)'}")

        ax.axvline(0, color="#333333", linestyle="--", linewidth=0.8, alpha=0.6)
        style_ax(ax, xlabel="Displacement from annotation (bp)",
                 ylabel="Density" if idx % ncols == 0 else "",
                 title=profile.replace("_", " "))
        if idx == 0:
            ax.legend(fontsize=8, frameon=False)

    # Hide empty panels
    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    savefig(fig, outdir / "displacement_distributions.png")
    log.info("  → displacement_distributions.png")


# ── Plot 3: Confidence vs Displacement Scatter ──────────────────────────────

def plot_confidence_vs_displacement(rows: List[dict], outdir: Path):
    """
    Scatter: x = confidence score, y = |displacement from annotation|,
    colored by library profile, panels for TSS vs TTS.

    Answers: "Does higher confidence actually mean better boundary placement?"
    """
    # Use the scorer's native confidence (alpha doesn't matter for this)
    data = [r for r in rows if r["alpha"] == 0.0]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)

    for panel_idx, end_type in enumerate(["tss", "tts"]):
        ax = axes[panel_idx]
        for profile in PROFILE_COLORS:
            subset = [r for r in data
                      if r["end_type"] == end_type and r["profile"] == profile]
            if not subset:
                continue
            confs = [r["confidence"] for r in subset]
            abs_disps = [min(500, abs(r["annot_displacement"])) for r in subset]
            ax.scatter(confs, abs_disps, s=12, alpha=0.35, rasterized=True,
                       color=PROFILE_COLORS[profile], edgecolors="none",
                       label=profile.replace("_", " "))

        style_ax(ax,
                 xlabel="Confidence score",
                 ylabel="|Displacement| (bp)" if panel_idx == 0 else "",
                 title=f"{'TSS (5′)' if end_type == 'tss' else 'TTS (3′)'}")
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-10, 520)

    legend_outside(axes[1])
    savefig(fig, outdir / "confidence_vs_displacement.png")
    log.info("  → confidence_vs_displacement.png")


# ── Plot 4: Alpha Sweep — Displacement Improvement ─────────────────────────

def plot_alpha_sweep(summary_rows: List[dict], outdir: Path):
    """
    Line plot: mean |displacement| vs alpha for each (profile, end_type).

    Answers: "How much does increasing alpha reduce boundary error?"
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True,
                             sharey=True)

    for panel_idx, end_type in enumerate(["tss", "tts"]):
        ax = axes[panel_idx]
        for profile in PROFILE_COLORS:
            subset = [r for r in summary_rows
                      if r["end_type"] == end_type and r["profile"] == profile]
            if not subset:
                continue
            subset.sort(key=lambda r: r["alpha"])
            alphas = [r["alpha"] for r in subset]
            mean_abs = [r["mean_abs_displacement"] for r in subset]
            ax.plot(alphas, mean_abs, "o-", color=PROFILE_COLORS[profile],
                    markersize=5, linewidth=1.5,
                    label=profile.replace("_", " "))

        style_ax(ax,
                 xlabel="Scoring alpha (0 = heuristic, 1 = sequence-aware)",
                 ylabel="Mean |displacement| (bp)" if panel_idx == 0 else "",
                 title=f"{'TSS (5′)' if end_type == 'tss' else 'TTS (3′)'}")
        ax.set_xlim(-0.05, 1.05)

    legend_outside(axes[1])
    savefig(fig, outdir / "alpha_sweep_displacement.png")
    log.info("  → alpha_sweep_displacement.png")


# ── Plot 5: Trust Category Distribution ─────────────────────────────────────

def plot_trust_categories(rows: List[dict], outdir: Path):
    """
    Stacked horizontal bar: proportion of boundaries in each trust category
    (trusted / uncertain / untrusted), one bar per (profile, end_type).

    Answers: "What fraction of boundaries should we actually trust?"
    """
    data = [r for r in rows if r["alpha"] == 0.0]

    # Group by (profile, end_type)
    groups = defaultdict(lambda: defaultdict(int))
    totals = defaultdict(int)
    for r in data:
        key = f"{r['profile']}  {'5′' if r['end_type'] == 'tss' else '3′'}"
        groups[key][r["trust_category"]] += 1
        totals[key] += 1

    labels = sorted(groups.keys())
    if not labels:
        return

    fig, ax = plt.subplots(figsize=(8, max(3.5, len(labels) * 0.45)),
                           constrained_layout=True)

    y = np.arange(len(labels))
    left = np.zeros(len(labels))
    cats = ["trusted", "uncertain", "untrusted"]

    for cat in cats:
        widths = [100 * groups[l].get(cat, 0) / max(1, totals[l]) for l in labels]
        ax.barh(y, widths, left=left, height=0.65,
                color=TRUST_COLORS[cat], edgecolor="white", linewidth=0.5,
                label=cat.capitalize())
        left += np.array(widths)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    style_ax(ax, xlabel="Boundaries (%)", title="End trust distribution")
    ax.set_xlim(0, 100)
    legend_outside(ax)

    savefig(fig, outdir / "trust_category_distribution.png")
    log.info("  → trust_category_distribution.png")


# ── Plot 6: Weight Decomposition — Which Signal Matters Most ────────────────

def plot_weight_decomposition(rows: List[dict], outdir: Path):
    """
    Grouped bar: for each (profile, end_type), show the mean contribution
    of seq_score, depth_score, and tech_penalty to the final confidence.

    Answers: "How much does sequence feature vs read depth vs tech bias
    contribute to the confidence score under each library profile?"
    """
    data = [r for r in rows if r["alpha"] == 0.0]

    groups = defaultdict(lambda: {"seq": [], "depth": [], "tech": []})
    for r in data:
        key = (r["profile"], r["end_type"])
        groups[key]["seq"].append(r["seq_score"] * r["seq_weight"])
        groups[key]["depth"].append(r["depth_score"] * r["depth_weight"])
        groups[key]["tech"].append(r["tech_penalty"] * r["tech_weight"])

    profiles_in_data = sorted(set(k[0] for k in groups))
    if not profiles_in_data:
        return

    fig, axes = plt.subplots(1, 2, figsize=(max(8, len(profiles_in_data) * 2.5), 4.5),
                             constrained_layout=True, sharey=True)

    for panel_idx, end_type in enumerate(["tss", "tts"]):
        ax = axes[panel_idx]
        profiles = [p for p in profiles_in_data if (p, end_type) in groups]
        x = np.arange(len(profiles))
        width = 0.25

        means_seq = [np.mean(groups[(p, end_type)]["seq"]) for p in profiles]
        means_depth = [np.mean(groups[(p, end_type)]["depth"]) for p in profiles]
        means_tech = [np.mean(groups[(p, end_type)]["tech"]) for p in profiles]

        ax.bar(x - width, means_seq, width, label="Sequence feature",
               color=PALETTE[2], edgecolor="none", alpha=0.85)
        ax.bar(x, means_depth, width, label="Read depth",
               color=PALETTE[4], edgecolor="none", alpha=0.85)
        ax.bar(x + width, means_tech, width, label="Tech penalty",
               color=PALETTE[5], edgecolor="none", alpha=0.85)

        ax.set_xticks(x)
        ax.set_xticklabels([p.replace("_", "\n") for p in profiles], fontsize=8)
        style_ax(ax,
                 ylabel="Mean weighted contribution" if panel_idx == 0 else "",
                 title=f"{'TSS (5′)' if end_type == 'tss' else 'TTS (3′)'}")

    axes[0].legend(fontsize=8, frameon=False, loc="upper right")
    savefig(fig, outdir / "weight_decomposition.png")
    log.info("  → weight_decomposition.png")


# ── Plot 7: Rescue Effectiveness ────────────────────────────────────────────

def plot_rescue_effectiveness(rows: List[dict], outdir: Path):
    """
    For boundaries that were rescued (rescue_reason != ""),
    compare displacement to annotation vs non-rescued boundaries.

    Answers: "When sequence features rescue a low-depth boundary,
    is it actually closer to the true annotation site?"
    """
    data = [r for r in rows if r["alpha"] == 0.0]

    rescued = [abs(r["annot_displacement"]) for r in data if r["rescue_reason"]]
    not_rescued = [abs(r["annot_displacement"]) for r in data if not r["rescue_reason"]]

    if not rescued or not not_rescued:
        log.warning("Not enough rescued boundaries for rescue_effectiveness plot")
        return

    fig, ax = plt.subplots(figsize=(6, 4.5))

    bins = np.linspace(0, 500, 50)
    ax.hist(not_rescued, bins=bins, alpha=0.55, density=True,
            color=PALETTE[7], edgecolor="none", label="Not rescued")
    ax.hist(rescued, bins=bins, alpha=0.65, density=True,
            color=PALETTE[2], edgecolor="none", label="Rescued by\nsequence features")

    # Add median lines
    med_r = np.median(rescued)
    med_nr = np.median(not_rescued)
    ax.axvline(med_r, color=PALETTE[2], linestyle="--", linewidth=1.2, alpha=0.7)
    ax.axvline(med_nr, color=PALETTE[7], linestyle="--", linewidth=1.2, alpha=0.7)

    # Stats annotation
    stats_text = (
        f"Rescued:  n={len(rescued):,}, med={med_r:.0f} bp\n"
        f"Not rescued: n={len(not_rescued):,}, med={med_nr:.0f} bp"
    )
    ax.text(0.97, 0.97, stats_text, transform=ax.transAxes,
            va="top", ha="right", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.5))

    style_ax(ax, xlabel="|Displacement from annotation| (bp)",
             ylabel="Density",
             title="Rescue effectiveness")
    ax.legend(fontsize=8, frameon=False)

    savefig(fig, outdir / "rescue_effectiveness.png")
    log.info("  → rescue_effectiveness.png")


# ── Plot 8: Combined Dashboard ──────────────────────────────────────────────

def plot_end_trust_dashboard(rows: List[dict], summary_rows: List[dict],
                             outdir: Path):
    """
    2×2 publication-quality dashboard combining the four most important views:
      (A) End trust by library type (bar)
      (B) Confidence vs displacement (scatter)
      (C) Alpha sweep (line)
      (D) Trust category distribution (stacked bar)
    """
    data_a0 = [r for r in rows if r["alpha"] == 0.0]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    panel_labels = ["A", "B", "C", "D"]

    # ── Panel A: End trust by library ──
    ax = axes[0, 0]
    profiles = [p for p in PROFILE_COLORS
                if any(r["profile"] == p for r in data_a0)]
    x = np.arange(len(profiles))
    width = 0.35
    groups = defaultdict(list)
    for r in data_a0:
        groups[(r["profile"], r["end_type"])].append(abs(r["annot_displacement"]))

    tss_pcts = [100 * sum(1 for d in groups.get((p, "tss"), []) if d <= 50)
                / max(1, len(groups.get((p, "tss"), [1]))) for p in profiles]
    tts_pcts = [100 * sum(1 for d in groups.get((p, "tts"), []) if d <= 50)
                / max(1, len(groups.get((p, "tts"), [1]))) for p in profiles]

    ax.bar(x - width / 2, tss_pcts, width, color=END_COLORS["tss"],
           edgecolor="none", alpha=0.85, label="TSS (5′)")
    ax.bar(x + width / 2, tts_pcts, width, color=END_COLORS["tts"],
           edgecolor="none", alpha=0.85, label="TTS (3′)")
    ax.set_xticks(x)
    ax.set_xticklabels([p.replace("_", "\n") for p in profiles], fontsize=7)
    style_ax(ax, ylabel="Within 50 bp (%)", title="End trust by library")
    ax.set_ylim(0, 105)
    ax.legend(fontsize=7, frameon=False, loc="lower right")

    # ── Panel B: Confidence vs displacement ──
    ax = axes[0, 1]
    for profile in PROFILE_COLORS:
        subset = [r for r in data_a0 if r["profile"] == profile]
        if not subset:
            continue
        confs = [r["confidence"] for r in subset]
        abs_disps = [min(500, abs(r["annot_displacement"])) for r in subset]
        ax.scatter(confs, abs_disps, s=8, alpha=0.25, rasterized=True,
                   color=PROFILE_COLORS[profile], edgecolors="none",
                   label=profile.replace("_", " "))
    style_ax(ax, xlabel="Confidence", ylabel="|Displacement| (bp)",
             title="Confidence vs displacement")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-10, 520)

    # ── Panel C: Alpha sweep ──
    ax = axes[1, 0]
    for profile in PROFILE_COLORS:
        for end_type, ls in [("tss", "-"), ("tts", "--")]:
            subset = [r for r in summary_rows
                      if r["end_type"] == end_type and r["profile"] == profile]
            if not subset:
                continue
            subset.sort(key=lambda r: r["alpha"])
            alphas_v = [r["alpha"] for r in subset]
            mean_abs = [r["mean_abs_displacement"] for r in subset]
            lbl = f"{profile.replace('_', ' ')} {'5′' if end_type == 'tss' else '3′'}"
            ax.plot(alphas_v, mean_abs, ls, color=PROFILE_COLORS[profile],
                    markersize=3, linewidth=1.2, alpha=0.7)

    # Simplified legend: just profiles (solid) + end type key
    handles_prof = [mpatches.Patch(color=PROFILE_COLORS[p], label=p.replace("_", " "))
                    for p in PROFILE_COLORS if any(r["profile"] == p for r in summary_rows)]
    handles_end = [plt.Line2D([0], [0], color="k", linestyle="-", label="TSS (5′)"),
                   plt.Line2D([0], [0], color="k", linestyle="--", label="TTS (3′)")]
    ax.legend(handles=handles_prof + handles_end, fontsize=6, frameon=False,
              loc="upper right", ncol=2)
    style_ax(ax, xlabel="Alpha", ylabel="Mean |displacement| (bp)",
             title="Alpha sweep")
    ax.set_xlim(-0.05, 1.05)

    # ── Panel D: Trust categories ──
    ax = axes[1, 1]
    cat_groups = defaultdict(lambda: defaultdict(int))
    cat_totals = defaultdict(int)
    for r in data_a0:
        key = f"{r['profile']}  {'5′' if r['end_type'] == 'tss' else '3′'}"
        cat_groups[key][r["trust_category"]] += 1
        cat_totals[key] += 1
    labels = sorted(cat_groups.keys())
    y = np.arange(len(labels))
    left = np.zeros(len(labels))
    for cat in ["trusted", "uncertain", "untrusted"]:
        widths = [100 * cat_groups[l].get(cat, 0) / max(1, cat_totals[l]) for l in labels]
        ax.barh(y, widths, left=left, height=0.65,
                color=TRUST_COLORS[cat], edgecolor="white", linewidth=0.5,
                label=cat.capitalize())
        left += np.array(widths)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6)
    ax.invert_yaxis()
    style_ax(ax, xlabel="Boundaries (%)", title="Trust distribution")
    ax.set_xlim(0, 100)
    ax.legend(fontsize=7, frameon=False, loc="lower right")

    # Panel labels
    for idx, (label, ax_) in enumerate(zip(panel_labels,
                                            [axes[0, 0], axes[0, 1],
                                             axes[1, 0], axes[1, 1]])):
        ax_.text(-0.08, 1.05, label, transform=ax_.transAxes,
                 fontsize=8, fontweight="normal", va="top")

    savefig(fig, outdir / "end_trust_dashboard.png")
    log.info("  → end_trust_dashboard.png")


# ── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate publication-quality end-trust evaluation plots.",
    )
    parser.add_argument("--boundaries", required=True,
                        help="Per-boundary TSV from evaluate_end_trust.py")
    parser.add_argument("--summary", required=True,
                        help="Summary TSV from evaluate_end_trust.py --summary")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for PNG files")

    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info(f"Loading boundary data: {args.boundaries}")
    rows = load_boundary_data(args.boundaries)
    log.info(f"  {len(rows)} rows")

    log.info(f"Loading summary data: {args.summary}")
    summary_rows = load_summary_data(args.summary)
    log.info(f"  {len(summary_rows)} rows")

    log.info("Generating plots...")
    plot_end_trust_by_library(rows, outdir)
    plot_displacement_distributions(rows, outdir)
    plot_confidence_vs_displacement(rows, outdir)
    plot_alpha_sweep(summary_rows, outdir)
    plot_trust_categories(rows, outdir)
    plot_weight_decomposition(rows, outdir)
    plot_rescue_effectiveness(rows, outdir)
    plot_end_trust_dashboard(rows, summary_rows, outdir)

    log.info("Done — all plots saved.")


if __name__ == "__main__":
    main()
