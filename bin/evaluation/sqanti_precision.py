#!/usr/bin/env python3
"""
sqanti_precision.py — End precision broken down by SQANTI structural category.

For each mode and each SQANTI category (FSM, ISM, NIC, NNC), computes:
  - 5′ end precision: JC-deduplicated TSS precision vs annotated GTF ends
  - 3′ end precision: JC-deduplicated TTS precision vs annotated GTF ends

Uses the same compute_jc_deduplicated_precision_recall() as TedEndPrecision
so metrics are consistent across the pipeline.  Multiple isoforms in the same
junction chain that map to the same annotation end count as ONE true positive.

This answers: "do NIC isoforms from TED modes land on better TSS/TTS positions
than NIC isoforms from FLAIR?" — separating end quality from splice-chain quality.

Produces:
  category_end_precision_bar.png
      Grouped bar chart. X = SQANTI category, bars = modes, Y = end precision.
      Two panels: TSS (top) and TTS (bottom). Count (n) shown in each bar.

  category_end_precision_heatmap.png
      Heatmap: modes (rows) × categories (cols), colour = end precision.
      Two panels side by side for TSS and TTS. Missing = grey.

Usage:
    python sqanti_precision.py \\
        --bed label1:bed1.bed label2:bed2.bed ... \\
        --gtf annotation.gtf \\
        --output outdir/ \\
        [--window 50] [--verbose]
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import numpy as np
import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler
    from flair_structural import (
        classify_transcripts_per_isoform,
        parse_reference,
    )
    from ted_end_precision import (
        parse_gtf_ends,
        parse_gtf_transcripts,
        compute_jc_deduplicated_precision_recall,
    )
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler
    from evaluation.flair_structural import (
        classify_transcripts_per_isoform,
        parse_reference,
    )
    from evaluation.ted_end_precision import (
        parse_gtf_ends,
        parse_gtf_transcripts,
        compute_jc_deduplicated_precision_recall,
    )

apply_rc()

# ── Constants ─────────────────────────────────────────────────────────────────

CATEGORY_ORDER  = ["FSM", "ISM", "NIC", "NNC"]
CATEGORY_COLORS = {
    "FSM": "#009E73",
    "ISM": "#56B4E9",
    "NIC": "#E69F00",
    "NNC": "#D55E00",
}
END_WINDOW = 50   # bp tolerance for end match


# ── Core computation ──────────────────────────────────────────────────────────

def compute_category_precision(
    isoforms: List[dict],
    annotated_ends: Dict[str, Dict[str, List[int]]],
    annot_transcripts: List[dict],
    window: int = END_WINDOW,
) -> Dict[str, Dict]:
    """Per-category 5′ and 3′ JC-deduplicated end precision for one mode.

    Uses compute_jc_deduplicated_precision_recall (GTF path, no orthogonal peaks)
    on the subset of isoforms in each SQANTI category.  This matches the
    authoritative metric used by TedEndPrecision in the Nextflow pipeline.

    Returns:
        {category: {"n": int, "tss_prec": float|None, "tts_prec": float|None}}
    """
    by_cat: Dict[str, List[dict]] = defaultdict(list)
    for iso in isoforms:
        cat = iso.get("category")
        if cat in CATEGORY_ORDER:
            by_cat[cat].append(iso)

    result = {}
    for cat in CATEGORY_ORDER:
        isos = by_cat.get(cat, [])
        n = len(isos)
        if n == 0:
            result[cat] = {"n": 0, "tss_prec": None, "tts_prec": None}
            continue

        metrics = compute_jc_deduplicated_precision_recall(
            isos, annotated_ends, annot_transcripts,
            window=window,
            peaks_5prime=None,   # GTF-only path — matches ted_end_precision GTF mode
            peaks_3prime=None,
        )
        result[cat] = {
            "n":        n,
            "tss_prec": metrics.get("5prime_dedup_precision"),
            "tts_prec": metrics.get("3prime_dedup_precision"),
        }
    return result


# ── Plot 1: Grouped bar chart ─────────────────────────────────────────────────

def plot_bar(
    data: Dict[str, Dict[str, Dict]],
    styler: ModeStyler,
    output_path: Path,
) -> None:
    """Two-panel grouped bar: X = category, grouped bars = modes, Y = end precision."""
    modes = list(data.keys())
    cats  = [c for c in CATEGORY_ORDER if any(
        data[m].get(c, {}).get("n", 0) > 0 for m in modes
    )]
    if not cats or not modes:
        return

    n_cats  = len(cats)
    n_modes = len(modes)
    width   = 0.8 / n_modes
    x       = np.arange(n_cats)

    fig, axes = plt.subplots(2, 1, figsize=(W2 * 0.75, W1 * 1.5), sharex=True)

    for ax_idx, (ax, end_key, end_label) in enumerate([
        (axes[0], "tss_prec", "5′ TSS  end precision"),
        (axes[1], "tts_prec", "3′ TTS  end precision"),
    ]):
        handles = []
        for i, mode in enumerate(modes):
            vals, ns = [], []
            for cat in cats:
                d = data[mode].get(cat, {"n": 0, "tss_prec": None, "tts_prec": None})
                v = d.get(end_key)
                vals.append(v if v is not None else 0.0)
                ns.append(d.get("n", 0))

            offset = (i - n_modes / 2 + 0.5) * width
            bars = ax.bar(
                x + offset, vals, width * 0.90,
                color=styler.color(mode), edgecolor="none", alpha=0.88,
            )
            handles.append(styler.legend_handle(mode, label=_short(mode), markersize=4))

            # n= count inside taller bars
            for bar, n in zip(bars, ns):
                if n > 0 and bar.get_height() > 0.12:
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        bar.get_height() / 2,
                        f"n={n}", ha="center", va="center",
                        fontsize=3.5, color="white", rotation=90,
                    )

        ax.set_ylim(0, 1.08)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:.0%}")
        )
        ax.axhline(1.0, color="#cccccc", lw=0.5, ls="--")
        style_ax(ax, ylabel=end_label)

        if ax_idx == 0:
            ax.legend(handles=handles, fontsize=5, frameon=False,
                      loc="lower right", ncol=max(1, n_modes // 3))

    # Category colour band along x-axis on the bottom panel
    for j, cat in enumerate(cats):
        axes[1].add_patch(plt.Rectangle(
            (j - 0.5, -0.13), 1, 0.06,
            color=CATEGORY_COLORS[cat], clip_on=False, transform=axes[1].get_xaxis_transform(),
        ))

    axes[1].set_xticks(x)
    axes[1].set_xticklabels(cats, fontsize=7)
    axes[1].set_xlabel("SQANTI structural category", fontsize=7, labelpad=10)

    fig.tight_layout(pad=0.4, h_pad=0.8)
    savefig(fig, output_path)


# ── Plot 2: Heatmap ───────────────────────────────────────────────────────────

def plot_heatmap(
    data: Dict[str, Dict[str, Dict]],
    output_path: Path,
) -> None:
    """Modes × categories heatmap, colour = end precision, two panels TSS/TTS."""
    modes = list(data.keys())
    cats  = [c for c in CATEGORY_ORDER if any(
        data[m].get(c, {}).get("n", 0) > 0 for m in modes
    )]
    if not cats or not modes:
        return

    fig, axes = plt.subplots(1, 2, figsize=(W2, W1 * 1.15))

    cmap = plt.get_cmap("RdYlGn").copy()
    cmap.set_bad(color="#e0e0e0")

    for ax, end_key, end_label in [
        (axes[0], "tss_prec", "5′ TSS end precision"),
        (axes[1], "tts_prec", "3′ TTS end precision"),
    ]:
        mat  = np.full((len(modes), len(cats)), np.nan)
        nmat = np.zeros((len(modes), len(cats)), dtype=int)

        for i, mode in enumerate(modes):
            for j, cat in enumerate(cats):
                d = data[mode].get(cat, {})
                v = d.get(end_key)
                if v is not None and d.get("n", 0) > 0:
                    mat[i, j]  = v
                    nmat[i, j] = d["n"]

        im = ax.imshow(mat, aspect="auto", vmin=0, vmax=1,
                       cmap=cmap, interpolation="nearest")

        for i in range(len(modes)):
            for j in range(len(cats)):
                if not np.isnan(mat[i, j]):
                    txt_col = "white" if (mat[i, j] < 0.35 or mat[i, j] > 0.82) else "#222222"
                    ax.text(j, i, f"{mat[i,j]:.0%}\n(n={nmat[i,j]})",
                            ha="center", va="center",
                            fontsize=4.5, color=txt_col, linespacing=1.3)

        # Category colour band across the top
        for j, cat in enumerate(cats):
            ax.add_patch(plt.Rectangle(
                (j - 0.5, -1.0), 1, 0.6,
                color=CATEGORY_COLORS[cat], clip_on=False,
            ))

        ax.set_xticks(range(len(cats)))
        ax.set_xticklabels(cats, fontsize=7, labelpad=8)
        ax.set_yticks(range(len(modes)))
        ax.set_yticklabels([_short(m) for m in modes], fontsize=6)
        ax.set_title(end_label, fontsize=7, pad=4)
        ax.tick_params(length=0)

        cb = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
        cb.set_label("End precision", fontsize=6)
        cb.ax.tick_params(labelsize=5)
        cb.ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:.0%}")
        )

    fig.tight_layout(pad=0.5)
    savefig(fig, output_path)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _short(mode: str) -> str:
    return (mode
            .replace("TED-", "")
            .replace("FLAIR-", "FL-")
            .replace("isoquant_", "IQ-"))


def _parse_label_path(entries: List[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for e in entries:
        if ":" not in e:
            continue
        label, path = e.split(":", 1)
        out[label] = path
    return out


def _load_categories_tsv(path: str) -> List[dict]:
    """Load a pre-classified isoform categories TSV written by flair_eval.py.

    Expected format (tab-separated, header row):
        isoform_name  category

    Returns list of dicts compatible with classify_transcripts_per_isoform output:
        [{"name": str, "category": str}, ...]
    """
    isoforms = []
    with open(path) as fh:
        header = fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            isoforms.append({"name": parts[0], "category": parts[1]})
    return isoforms


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bed", nargs="+", required=True,
                        help="label:path pairs for BED12 isoform files")
    parser.add_argument("--categories-tsv", nargs="+", default=[],
                        help="label:path pairs for pre-classified isoform category TSVs "
                             "(produced by flair_eval.py --categories-output). "
                             "When provided, skips GTF re-parsing and re-classification.")
    parser.add_argument("--gtf", required=True,
                        help="Reference GTF for annotated end positions")
    parser.add_argument("--window", type=int, default=END_WINDOW,
                        help=f"bp tolerance for end match (default {END_WINDOW})")
    parser.add_argument("--output", required=True)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    # Always need annotated ends + transcripts for compute_jc_deduplicated_precision_recall
    if args.verbose:
        print("  Parsing reference GTF for end positions...", file=sys.stderr)
    annotated_ends = parse_gtf_ends(args.gtf)
    annot_transcripts = parse_gtf_transcripts(args.gtf)
    if args.verbose:
        n_tss = sum(len(v["tss"]) for v in annotated_ends.values())
        n_tts = sum(len(v["tts"]) for v in annotated_ends.values())
        print(f"  Reference: {n_tss} annotated TSS, {n_tts} annotated TTS, "
              f"{len(annot_transcripts)} transcripts", file=sys.stderr)

    # Pre-classified categories TSVs (from flair_eval.py) — skip GTF parse + classification
    cat_paths = _parse_label_path(args.categories_tsv) if args.categories_tsv else {}

    # Only parse reference structures when at least one label lacks pre-classified categories
    bed_paths = _parse_label_path(args.bed)
    needs_classification = [lbl for lbl in bed_paths if lbl not in cat_paths]
    if needs_classification:
        if args.verbose:
            print("  Parsing reference GTF for classification "
                  f"(needed for: {', '.join(needs_classification)})...", file=sys.stderr)
        refjuncs, refjuncchains, refseends = parse_reference(args.gtf)
    else:
        refjuncs = refjuncchains = refseends = None

    data: Dict[str, Dict[str, Dict]] = {}

    for label, path in bed_paths.items():
        if not Path(path).exists():
            print(f"WARNING: {path} not found — skipping {label}", file=sys.stderr)
            continue

        if label in cat_paths:
            # Fast path: load pre-classified categories from TSV
            if args.verbose:
                print(f"  Loading pre-classified categories for {label}...", file=sys.stderr)
            try:
                isoforms = _load_categories_tsv(cat_paths[label])
            except Exception as e:
                print(f"WARNING: failed to load categories TSV for {label}: {e} — "
                      "falling back to re-classification", file=sys.stderr)
                isoforms = None
        else:
            isoforms = None

        if isoforms is None:
            # Slow path: classify from BED + GTF reference structures
            if args.verbose:
                print(f"  Classifying {label} from BED...", file=sys.stderr)
            try:
                isoforms = classify_transcripts_per_isoform(
                    path, refjuncs, refjuncchains, refseends
                )
            except Exception as e:
                print(f"WARNING: classification failed for {label}: {e}", file=sys.stderr)
                continue

        if not isoforms:
            continue
        data[label] = compute_category_precision(
            isoforms, annotated_ends, annot_transcripts, args.window
        )
        if args.verbose:
            for cat, d in data[label].items():
                if d["n"] > 0:
                    tss_str = f"{d['tss_prec']:.1%}" if d['tss_prec'] is not None else "N/A"
                    tts_str = f"{d['tts_prec']:.1%}" if d['tts_prec'] is not None else "N/A"
                    print(f"    {cat:4s}  n={d['n']:4d}  TSS={tss_str}  TTS={tts_str}",
                          file=sys.stderr)

    if not data:
        print("No data — exiting", file=sys.stderr)
        sys.exit(1)

    styler = ModeStyler(list(data.keys()))
    plot_bar(data, styler, out / "category_end_precision_bar.png")
    plot_heatmap(data, out / "category_end_precision_heatmap.png")
    print(f"Saved category end precision plots to {args.output}")


if __name__ == "__main__":
    main()
