#!/usr/bin/env python3
"""
sjc_alt_end_analysis.py — Within-SJC joint end redundancy analysis.

For each assembler, groups multi-exon isoforms by junction chain (SJC) and
asks: within each SJC group, does every isoform have a unique combination of
5' and 3' orthogonal peak support?

The core question: are any two isoforms in the same SJC group hitting the
*same peak on both ends* — meaning one is fully redundant with the other?
If they differ on even one end, they may be biologically justified.

Each isoform in a multi-isoform SJC group is assigned a
(5prime_peak_iv, 3prime_peak_iv) pair.  Two isoforms are fully redundant if
their pair is identical.  If only one end matches, they are partially redundant.

Classification per isoform (relative to all others in its SJC group):

  Fully unique        : (5p, 3p) pair not shared with any other group member
  5p-only unique      : 5p peak distinct, 3p peak shared with another member
  3p-only unique      : 3p peak distinct, 5p peak shared with another member
  Fully redundant     : both ends shared with the same other group member
  Fully unsupported   : no peak on either end

Three panels:

Panel A — Per-assembler stacked bar of isoform classifications
  Shows what fraction of isoforms in multi-isoform SJC groups fall into
  each category.

Panel B — Group-level redundancy: how many groups are "fully justified"?
  For each SJC group, count how many isoforms are fully unique.
  Histogram of (n_fully_unique / n_isoforms) ratio — ratio = 1.0 means
  every isoform in the group has its own unique (5p, 3p) peak pair.

Panel C — Pairwise end-distance ECDF for fully-redundant vs fully-unique pairs
  For each pair of isoforms in the same SJC group, compute the minimum of
  (5p_distance, 3p_distance) and the sum.  Split by whether the pair is
  fully redundant, partially redundant, or fully distinct.

Inputs:
  --bed       label:path  (BED12 isoform files, one per assembler)
  --cage-peaks   BED6 CAGE peaks
  --qs-peaks     BED6 dRNA peaks
  --cage-plus/--cage-minus/--qs-plus/--qs-minus  bedGraph signal tracks
  --output       output directory
"""

from __future__ import annotations

import argparse
import logging
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, PALETTE
    from signal_utils import (
        parse_isoforms, group_by_junction_chain,
        load_signal_tracks,
    )
    from ted_end_precision import parse_peaks_bed, _nearest_annot
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, PALETTE
    from evaluation.signal_utils import (
        parse_isoforms, group_by_junction_chain,
        load_signal_tracks,
    )
    from evaluation.ted_end_precision import parse_peaks_bed, _nearest_annot

apply_rc()
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s")
log = logging.getLogger(__name__)

WINDOW = 50  # peak-match window in bp


# ── Per-assembler SJC analysis ───────────────────────────────────────────────

def _analyse_assembler(
    label: str,
    isoforms: List[dict],
    cage_peaks: Dict[Tuple[str, str], List[Tuple[int, int]]],
    qs_peaks:   Dict[Tuple[str, str], List[Tuple[int, int]]],
) -> dict:
    """Analyse joint (5p, 3p) peak-pair redundancy within SJC groups.

    For each multi-isoform SJC group, assigns every isoform a
    (5prime_peak_iv, 3prime_peak_iv) tuple and classifies it as:

      fully_unique   : peak-pair not shared with any other group member
      unique_5p_only : 5p peak is unique in the group; 3p is shared
      unique_3p_only : 3p peak is unique in the group; 5p is shared
      fully_redundant: both ends shared with at least one other member
      fully_unsupported: no peak on either end

    Also records per-group (n_isoforms, n_fully_unique) for Panel B and
    pairwise distances split by redundancy class for Panel C.
    """
    groups = group_by_junction_chain(isoforms)

    n_fully_unique    = 0
    n_unique_5p_only  = 0
    n_unique_3p_only  = 0
    n_fully_redundant = 0
    n_fully_unsup     = 0

    # Panel B: (n_isoforms, n_fully_unique) per group
    group_stats: List[Tuple[int, int]] = []

    # Panel C: pairwise end distances split by redundancy class
    # Each entry: (5p_dist, 3p_dist)
    pairs_fully_redundant: List[Tuple[float, float]] = []
    pairs_partial:         List[Tuple[float, float]] = []
    pairs_fully_distinct:  List[Tuple[float, float]] = []

    for jc_key, members in groups.items():
        if len(members) < 2:
            continue
        chrom, strand, _ = jc_key
        cage_list = cage_peaks.get((chrom, strand), [])
        qs_list   = qs_peaks.get((chrom, strand), [])

        # Assign each member its (5p_peak, 3p_peak) pair
        enriched = []
        for iso in members:
            p5 = _nearest_annot(iso["tss"], cage_list, WINDOW)
            p3 = _nearest_annot(iso["tts"], qs_list,  WINDOW)
            enriched.append({
                "name": iso["name"],
                "tss":  iso["tss"],
                "tts":  iso["tts"],
                "p5":   p5,
                "p3":   p3,
            })

        # Count how many times each peak-pair appears in this group
        pair_counts: Dict[Tuple, int] = defaultdict(int)
        for e in enriched:
            pair_counts[(e["p5"], e["p3"])] += 1

        # Count occurrences of each individual end peak across the group
        p5_counts: Dict[Optional[Tuple[int,int]], int] = defaultdict(int)
        p3_counts: Dict[Optional[Tuple[int,int]], int] = defaultdict(int)
        for e in enriched:
            p5_counts[e["p5"]] += 1
            p3_counts[e["p3"]] += 1

        # Classify each isoform
        n_fu_this_group = 0
        for e in enriched:
            p5, p3 = e["p5"], e["p3"]
            no_5p = p5 is None
            no_3p = p3 is None

            if no_5p and no_3p:
                n_fully_unsup += 1
                continue

            # Is the full pair unique in this group?
            pair_unique = pair_counts[(p5, p3)] == 1
            # Is each individual end unique in this group?
            p5_unique = (not no_5p) and p5_counts[p5] == 1
            p3_unique = (not no_3p) and p3_counts[p3] == 1

            if pair_unique:
                # Both ends together are unique — fully justified
                n_fully_unique += 1
                n_fu_this_group += 1
            elif p5_unique and not p3_unique:
                n_unique_5p_only += 1
            elif p3_unique and not p5_unique:
                n_unique_3p_only += 1
            else:
                # Neither end is unique in this group → fully redundant
                n_fully_redundant += 1

        group_stats.append((len(enriched), n_fu_this_group))

        # Panel C: pairwise distances
        for a, b in combinations(enriched, 2):
            d5 = abs(a["tss"] - b["tss"])
            d3 = abs(a["tts"] - b["tts"])
            pa5, pa3 = a["p5"], a["p3"]
            pb5, pb3 = b["p5"], b["p3"]
            same_5p = (pa5 is not None and pa5 == pb5)
            same_3p = (pa3 is not None and pa3 == pb3)
            if same_5p and same_3p:
                pairs_fully_redundant.append((d5, d3))
            elif same_5p or same_3p:
                pairs_partial.append((d5, d3))
            else:
                pairs_fully_distinct.append((d5, d3))

    return {
        "n_fully_unique":    n_fully_unique,
        "n_unique_5p_only":  n_unique_5p_only,
        "n_unique_3p_only":  n_unique_3p_only,
        "n_fully_redundant": n_fully_redundant,
        "n_fully_unsup":     n_fully_unsup,
        "group_stats":               group_stats,
        "pairs_fully_redundant":     pairs_fully_redundant,
        "pairs_partial":             pairs_partial,
        "pairs_fully_distinct":      pairs_fully_distinct,
    }


# ── Panel A — Per-assembler isoform classification (stacked bar) ─────────────

def plot_panel_a(results: dict, outdir: Path) -> None:
    """Horizontal stacked bar: classification of isoforms in multi-isoform
    SJC groups by joint (5p, 3p) peak-pair uniqueness.
    """
    labels = [l for l in results if results[l] is not None]
    labels = [l for l in labels if (
        results[l]["n_fully_unique"] + results[l]["n_unique_5p_only"] +
        results[l]["n_unique_3p_only"] + results[l]["n_fully_redundant"] +
        results[l]["n_fully_unsup"]
    ) > 0]
    if not labels:
        log.warning("Panel A: no multi-isoform SJC groups, skipping")
        return

    fig, ax = plt.subplots(figsize=(W2, max(W1, 0.28 * len(labels) + 0.5)))
    y = np.arange(len(labels))
    h = 0.65

    c_fu   = PALETTE[2]   # green  — fully unique (5p and 3p both distinct)
    c_5p   = PALETTE[0]   # blue   — 5p unique, 3p shared
    c_3p   = PALETTE[1]   # orange — 3p unique, 5p shared
    c_red  = "#D55E00"    # vermillion — fully redundant (both ends shared)
    c_unsup = "#999999"   # grey   — no peak on either end

    totals = [
        results[l]["n_fully_unique"] + results[l]["n_unique_5p_only"] +
        results[l]["n_unique_3p_only"] + results[l]["n_fully_redundant"] +
        results[l]["n_fully_unsup"]
        for l in labels
    ]

    def frac(key, l):
        t = totals[labels.index(l)]
        return results[l][key] / t if t > 0 else 0.0

    specs = [
        ("n_fully_unique",    c_fu,    "Fully unique (distinct 5p AND 3p peak)"),
        ("n_unique_5p_only",  c_5p,    "5p unique, 3p shared"),
        ("n_unique_3p_only",  c_3p,    "3p unique, 5p shared"),
        ("n_fully_redundant", c_red,   "Fully redundant (same 5p AND 3p peak)"),
        ("n_fully_unsup",     c_unsup, "Fully unsupported (no peak on either end)"),
    ]

    left = np.zeros(len(labels))
    for key, col, lbl in specs:
        vals = [frac(key, l) for l in labels]
        ax.barh(y, vals, h, left=left, label=lbl, color=col)
        left += np.array(vals)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:.0%}")
    )
    ax.legend(fontsize=5.5, loc="lower right", frameon=False,
              handlelength=1.2, ncol=1)
    for i, (l, t) in enumerate(zip(labels, totals)):
        ax.text(1.02, i, f"n={t:,}", va="center", fontsize=5,
                transform=ax.get_yaxis_transform())

    style_ax(ax, xlabel="Fraction of isoforms in multi-isoform SJC groups",
             title="Joint (5p, 3p) peak-pair uniqueness within SJC groups")
    ax.axvline(0.0, color="#888888", linewidth=0.4)
    fig.tight_layout()
    savefig(fig, outdir / "sjc_alt_end_panel_a")


# ── Panel B — Group-level: fraction of isoforms that are fully unique ────────

def plot_panel_b(results: dict, outdir: Path) -> None:
    """Histogram of (n_fully_unique / n_isoforms) per SJC group.

    Ratio = 1.0: every isoform in the group has its own unique (5p, 3p) pair.
    Ratio = 0.0: no isoform in the group is uniquely justified.
    """
    labels = [l for l in results if results[l] is not None]
    labels = [l for l in labels if results[l]["group_stats"]]
    if not labels:
        log.warning("Panel B: no group data, skipping")
        return

    n = len(labels)
    fig, axes = plt.subplots(1, n, figsize=(max(W1, 1.8 * n), W1 * 1.3),
                             sharey=True, sharex=True, squeeze=False)
    axes = axes[0]

    bin_edges = np.linspace(0, 1, 21)

    strata = [(2, 2), (3, 3), (4, 4), (5, 9999)]
    strata_labels = ["2 isoforms", "3 isoforms", "4 isoforms", "5+ isoforms"]
    strata_colors = [PALETTE[0], PALETTE[1], PALETTE[2], PALETTE[3]]

    for i, label in enumerate(labels):
        ax = axes[i]
        ratios_by_stratum = [[] for _ in strata]
        for n_iso, n_fu in results[label]["group_stats"]:
            ratio = n_fu / n_iso if n_iso > 0 else 0.0
            for si, (lo, hi) in enumerate(strata):
                if lo <= n_iso <= hi:
                    ratios_by_stratum[si].append(ratio)
                    break

        bottom = np.zeros(len(bin_edges) - 1)
        for ratios, slbl, col in zip(ratios_by_stratum, strata_labels, strata_colors):
            if not ratios:
                continue
            counts, _ = np.histogram(ratios, bins=bin_edges)
            ax.bar(bin_edges[:-1], counts, width=np.diff(bin_edges),
                   bottom=bottom, align="edge",
                   label=f"{slbl} (n={len(ratios):,})",
                   color=col, alpha=0.8, linewidth=0)
            bottom += counts

        ax.axvline(1.0, color="#333333", linewidth=0.8, linestyle="--", alpha=0.7)
        ax.set_title(label, fontsize=6)
        ax.set_xlabel("Fully-unique isoforms / group size", fontsize=6)
        if i == 0:
            ax.set_ylabel("SJC groups", fontsize=6)
        ax.legend(fontsize=4.5, frameon=False, loc="upper left")
        style_ax(ax)

    fig.suptitle(
        "Fraction of isoforms per SJC group with a unique (5p, 3p) peak pair\n"
        "Ratio = 1: all isoforms justified  |  Ratio < 1: at least one redundant",
        fontsize=6, y=1.03,
    )
    fig.tight_layout()
    savefig(fig, outdir / "sjc_alt_end_panel_b")


# ── Panel C — Pairwise end distances by redundancy class ─────────────────────

def plot_panel_c(results: dict, outdir: Path) -> None:
    """ECDF of pairwise (5p_dist, 3p_dist) within SJC groups, split by
    whether the pair is fully redundant, partially sharing, or fully distinct.

    Uses the sum of 5p and 3p distances as the x-axis (total end divergence).
    """
    labels = [l for l in results if results[l] is not None]
    labels = [l for l in labels if (
        results[l]["pairs_fully_redundant"] or
        results[l]["pairs_partial"] or
        results[l]["pairs_fully_distinct"]
    )]
    if not labels:
        log.warning("Panel C: no pairwise distance data, skipping")
        return

    n = len(labels)
    fig, axes = plt.subplots(1, n, figsize=(max(W1, 1.6 * n), W1 * 1.1),
                             sharey=True, sharex=True, squeeze=False)
    axes = axes[0]

    def _total_dist(pairs):
        return [d5 + d3 for d5, d3 in pairs]

    x_max = 0
    for l in labels:
        r = results[l]
        all_d = (
            _total_dist(r["pairs_fully_redundant"]) +
            _total_dist(r["pairs_partial"]) +
            _total_dist(r["pairs_fully_distinct"])
        )
        if all_d:
            x_max = max(x_max, np.percentile(all_d, 99))
    x_max = max(x_max, 300)

    def _ecdf(data):
        xs = np.sort(data)
        ys = np.arange(1, len(xs) + 1) / len(xs)
        return xs, ys

    grp_spec = [
        ("pairs_fully_redundant", "Fully redundant (same 5p+3p peak)", "#D55E00", "--"),
        ("pairs_partial",         "Partially sharing (one end shared)", "#E69F00", ":"),
        ("pairs_fully_distinct",  "Fully distinct (different 5p+3p)",   PALETTE[2], "-"),
    ]

    for i, label in enumerate(labels):
        ax = axes[i]
        r  = results[label]
        for key, lbl, col, ls in grp_spec:
            d = _total_dist(r[key])
            if not d:
                continue
            xs, ys = _ecdf(d)
            ax.plot(xs, ys, color=col, linestyle=ls, linewidth=1.0,
                    label=f"{lbl} (n={len(d):,})")

        ax.set_xlim(0, x_max)
        ax.set_ylim(0, 1.05)
        ax.set_title(label, fontsize=6)
        ax.set_xlabel("Sum of 5p + 3p end distance (bp)", fontsize=6)
        ax.legend(fontsize=5, frameon=False, loc="lower right")
        style_ax(ax, ylabel=("Cumulative fraction" if i == 0 else None))

    fig.suptitle(
        "Pairwise end divergence within SJC groups\n"
        "x = 5p_dist + 3p_dist between each pair of isoforms in the same group",
        fontsize=6, y=1.02,
    )
    fig.tight_layout()
    savefig(fig, outdir / "sjc_alt_end_panel_c")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bed", nargs="+", required=True,
                        help="label:path pairs for BED12 isoform files")
    parser.add_argument("--cage-peaks", required=True, help="CAGE peaks BED6")
    parser.add_argument("--qs-peaks",   required=True, help="dRNA peaks BED6")
    parser.add_argument("--cage-plus",  required=True, help="CAGE bedGraph (+ strand)")
    parser.add_argument("--cage-minus", required=True, help="CAGE bedGraph (- strand)")
    parser.add_argument("--qs-plus",    required=True, help="dRNA bedGraph (+ strand)")
    parser.add_argument("--qs-minus",   required=True, help="dRNA bedGraph (- strand)")
    parser.add_argument("--output",     required=True, help="Output directory")
    parser.add_argument("--verbose",    action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Parse inputs ─────────────────────────────────────────────────────
    beds_by_method: dict = {}
    for entry in args.bed:
        label, path = entry.split(":", 1)
        isos = parse_isoforms(path)
        beds_by_method[label] = isos
        log.info("Loaded %d isoforms for %s", len(isos), label)

    cage_peaks = parse_peaks_bed(args.cage_peaks)
    qs_peaks   = parse_peaks_bed(args.qs_peaks)
    log.info("Loaded %d / %d CAGE / dRNA peak groups",
             len(cage_peaks), len(qs_peaks))

    cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
        args.cage_plus, args.cage_minus, args.qs_plus, args.qs_minus,
    )

    # ── Run analysis per assembler ────────────────────────────────────────
    results: dict = {}
    for label, isoforms in beds_by_method.items():
        r = _analyse_assembler(label, isoforms, cage_peaks, qs_peaks)
        total = (
            r["n_fully_unique"] + r["n_unique_5p_only"] +
            r["n_unique_3p_only"] + r["n_fully_redundant"] + r["n_fully_unsup"]
        )
        if total == 0:
            log.info("%s: no multi-isoform SJC groups, skipping", label)
            results[label] = None
        else:
            n_groups = len(r["group_stats"])
            n_fully_justified = sum(1 for ni, nfu in r["group_stats"] if nfu == ni)
            log.info(
                "%s: %d isoforms in %d multi-iso groups — "
                "fully_unique=%d  5p_only=%d  3p_only=%d  "
                "fully_redundant=%d  unsupported=%d  "
                "groups_fully_justified=%d/%d (%.0f%%)",
                label, total, n_groups,
                r["n_fully_unique"], r["n_unique_5p_only"], r["n_unique_3p_only"],
                r["n_fully_redundant"], r["n_fully_unsup"],
                n_fully_justified, n_groups,
                100 * n_fully_justified / n_groups if n_groups else 0,
            )
            results[label] = r

    plot_panel_a(results, outdir)
    plot_panel_b(results, outdir)
    plot_panel_c(results, outdir)

    log.info("Done — output in %s", outdir)


if __name__ == "__main__":
    main()
