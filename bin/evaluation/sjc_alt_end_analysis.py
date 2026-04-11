#!/usr/bin/env python3
"""
sjc_alt_end_analysis.py — Analyse alternative ends within splice-junction-chain groups.

For each assembler, groups multi-exon isoforms by junction chain (SJC) and
examines groups with >1 isoform.  Two panels per end type (5′ and 3′):

Panel A — Alt-end TP/FP breakdown (stacked bar)
  For alt-end isoforms: fraction that are TP-unique (hit a peak not matched by
  another member), TP-redundant (hit a peak already matched), FP (no peak hit).

Panel B — Boundary signal at representative vs alternative ends (violin)
  Within each multi-end SJC group the representative end (highest read count)
  vs alternative ends, coloured by TP/FP status.

Inputs mirror the cumulative-signal-plot pattern:
  --bed       label:path  (BED12 isoform files, one per assembler)
  --read-map  label:path  (read-map files, optional per assembler)
  --cage-peaks   BED6 CAGE peaks
  --qs-peaks     BED6 QuantSeq peaks
  --cage-plus/--cage-minus/--qs-plus/--qs-minus  bedGraph signal tracks
  --output       output directory
"""

from __future__ import annotations

import argparse
import logging
import sys
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, PALETTE
    from signal_utils import (
        parse_isoforms, group_by_junction_chain, isoform_signal,
        load_signal_tracks,
    )
    from cumulative_signal_plot import load_read_map, _lookup_read_count
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, PALETTE
    from evaluation.signal_utils import (
        parse_isoforms, group_by_junction_chain, isoform_signal,
        load_signal_tracks,
    )
    from evaluation.cumulative_signal_plot import load_read_map, _lookup_read_count

apply_rc()
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s")
log = logging.getLogger(__name__)

WINDOW = 50  # match distance window for peak-based TP


# ── Peak helpers (same logic as ted_end_precision.py) ───────────────────────

def _parse_peaks_bed(path: str) -> Dict[Tuple[str, str], List[int]]:
    """Parse BED6 peaks → {(chrom, strand): sorted midpoints}."""
    peaks: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            chrom, start, end, strand = cols[0], int(cols[1]), int(cols[2]), cols[5]
            peaks[(chrom, strand)].append((start + end) // 2)
    for k in peaks:
        peaks[k].sort()
    return peaks


def _nearest_peak(pos: int, sorted_peaks: List[int], window: int = WINDOW
                  ) -> Optional[int]:
    """Return nearest peak midpoint within *window*, or None."""
    if not sorted_peaks:
        return None
    idx = bisect_left(sorted_peaks, pos)
    best, best_d = None, window + 1
    for i in (idx - 1, idx):
        if 0 <= i < len(sorted_peaks):
            d = abs(pos - sorted_peaks[i])
            if d < best_d:
                best_d, best = d, sorted_peaks[i]
    return best if best_d <= window else None


# ── Per-assembler SJC analysis ──────────────────────────────────────────────

def _analyse_assembler(
    label: str,
    isoforms: List[dict],
    read_counts: dict[str, int],
    peaks: Dict[Tuple[str, str], List[int]],
    cage_p, cage_m, qs_p, qs_m,
    end_type: str,  # "tss" or "tts"
) -> dict:
    """Analyse alt-end SJC groups for one assembler and one end type.

    Returns dict with:
      tp_unique, tp_redundant, fp  — counts for Panel A
      rep_tp_signals, rep_fp_signals — representative-end signals for Panel B
      alt_tp_signals, alt_fp_signals — alternative-end signals for Panel B
    """
    groups = group_by_junction_chain(isoforms)

    tp_unique = 0
    tp_redundant = 0
    fp = 0
    rep_tp_signals: list[float] = []
    rep_fp_signals: list[float] = []
    alt_tp_signals: list[float] = []
    alt_fp_signals: list[float] = []

    for jc_key, members in groups.items():
        if len(members) < 2:
            continue
        chrom, strand, _ = jc_key
        peak_list = peaks.get((chrom, strand), [])

        # Compute per-member: end position, signal, read count, peak match
        enriched = []
        for iso in members:
            pos = iso["tss"] if end_type == "tss" else iso["tts"]
            sig_tss, sig_tts = isoform_signal(iso, cage_p, cage_m, qs_p, qs_m)
            sig = sig_tss if end_type == "tss" else sig_tts
            rc = _lookup_read_count(iso["name"], read_counts)
            peak_hit = _nearest_peak(pos, peak_list)
            enriched.append({
                "name": iso["name"],
                "pos": pos,
                "signal": sig,
                "read_count": rc if rc is not None else 0,
                "peak_hit": peak_hit,
            })

        # Representative = highest read count (tie-break: first)
        enriched.sort(key=lambda x: -x["read_count"])
        representative = enriched[0]
        alternatives = enriched[1:]

        # Track peaks already matched by any member to detect redundancy
        # Process all members in read-count order; first match to a peak is "unique"
        seen_peaks: Set[Optional[int]] = set()

        # Classify representative
        rp = representative["peak_hit"]
        if rp is not None:
            seen_peaks.add(rp)
            # representative is always "unique" for its peak
        # We only count alternatives for the alt-end breakdown

        # Classify each alternative
        for alt in alternatives:
            ap = alt["peak_hit"]
            if ap is None:
                fp += 1
                alt_fp_signals.append(alt["signal"])
            elif ap in seen_peaks:
                tp_redundant += 1
                alt_tp_signals.append(alt["signal"])
            else:
                tp_unique += 1
                seen_peaks.add(ap)
                alt_tp_signals.append(alt["signal"])

        # Panel B: representative signal
        if representative["peak_hit"] is not None:
            rep_tp_signals.append(representative["signal"])
        else:
            rep_fp_signals.append(representative["signal"])

    return {
        "tp_unique": tp_unique,
        "tp_redundant": tp_redundant,
        "fp": fp,
        "rep_tp_signals": rep_tp_signals,
        "rep_fp_signals": rep_fp_signals,
        "alt_tp_signals": alt_tp_signals,
        "alt_fp_signals": alt_fp_signals,
    }


# ── Panel A — Stacked bar ──────────────────────────────────────────────────

def plot_panel_a(
    results: dict[str, dict],
    end_label: str,
    outdir: Path,
) -> None:
    """Stacked bar: alt-end TP-unique / TP-redundant / FP per assembler."""
    labels = [l for l in results if results[l] is not None]
    # Filter labels with at least one alt-end isoform
    labels = [l for l in labels
              if (results[l]["tp_unique"] + results[l]["tp_redundant"] + results[l]["fp"]) > 0]
    if not labels:
        log.warning("Panel A (%s): no assemblers with alt-end groups, skipping", end_label)
        return

    tp_u = [results[l]["tp_unique"] for l in labels]
    tp_r = [results[l]["tp_redundant"] for l in labels]
    fp_  = [results[l]["fp"] for l in labels]

    # Convert to fractions
    totals = [u + r + f for u, r, f in zip(tp_u, tp_r, fp_)]
    frac_u = [u / t if t > 0 else 0 for u, t in zip(tp_u, totals)]
    frac_r = [r / t if t > 0 else 0 for r, t in zip(tp_r, totals)]
    frac_f = [f / t if t > 0 else 0 for f, t in zip(fp_, totals)]

    fig, ax = plt.subplots(figsize=(W2, W1))
    y = np.arange(len(labels))
    h = 0.6

    ax.barh(y, frac_u, h, label="TP-unique", color=PALETTE[2])
    ax.barh(y, frac_r, h, left=frac_u, label="TP-redundant", color=PALETTE[0])
    ax.barh(y, frac_f, h, left=[u + r for u, r in zip(frac_u, frac_r)],
            label="FP (no peak)", color=PALETTE[3])

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlim(0, 1)
    ax.legend(fontsize=6, loc="lower right", frameon=False)
    # Annotate absolute counts on each bar
    for i, l in enumerate(labels):
        ax.text(1.02, i, f"n={totals[i]}", va="center", fontsize=5,
                transform=ax.get_yaxis_transform())
    style_ax(ax, xlabel="Fraction of alternative ends",
             title=f"Alt-end classification within SJC groups ({end_label})")
    fig.tight_layout()
    savefig(fig, outdir / f"sjc_alt_end_panel_a_{end_label}")


# ── Panel B — Violin / strip plot ───────────────────────────────────────────

def plot_panel_b(
    results: dict[str, dict],
    end_label: str,
    outdir: Path,
) -> None:
    """Violin: boundary signal for representative vs alternative ends."""
    labels = [l for l in results if results[l] is not None]
    labels = [l for l in labels
              if (len(results[l]["rep_tp_signals"]) + len(results[l]["rep_fp_signals"])
                  + len(results[l]["alt_tp_signals"]) + len(results[l]["alt_fp_signals"])) > 0]
    if not labels:
        log.warning("Panel B (%s): no data, skipping", end_label)
        return

    n = len(labels)
    fig, axes = plt.subplots(1, n, figsize=(W2, W1), sharey=True, squeeze=False)
    axes = axes[0]

    for i, label in enumerate(labels):
        ax = axes[i]
        r = results[label]
        # Combine representative and alternative signals with category labels
        categories = []
        signals = []
        colors = []

        for sig, cat, col in [
            (r["rep_tp_signals"],  "Rep\n(TP)",  PALETTE[2]),
            (r["rep_fp_signals"],  "Rep\n(FP)",  PALETTE[3]),
            (r["alt_tp_signals"],  "Alt\n(TP)",  PALETTE[0]),
            (r["alt_fp_signals"],  "Alt\n(FP)",  PALETTE[3]),
        ]:
            if sig:
                categories.append(cat)
                signals.append(sig)
                colors.append(col)

        if not signals:
            ax.set_visible(False)
            continue

        # Box plot with strip overlay
        positions = list(range(len(signals)))
        bp = ax.boxplot(signals, positions=positions, widths=0.5,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color="black", linewidth=1))
        for patch, col in zip(bp["boxes"], colors):
            patch.set_facecolor(col)
            patch.set_alpha(0.4)

        # Strip (jittered points)
        for j, (sigs, col) in enumerate(zip(signals, colors)):
            jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(sigs))
            ax.scatter(j + jitter, sigs, s=4, alpha=0.5, color=col,
                       edgecolors="none", zorder=3)

        ax.set_xticks(positions)
        ax.set_xticklabels(categories, fontsize=5)
        ax.set_title(label, fontsize=6)

        if i == 0:
            style_ax(ax, ylabel=f"{end_label} boundary signal")
        else:
            style_ax(ax)

    fig.suptitle(f"Boundary signal: representative vs alternative ends ({end_label})",
                 fontsize=7, y=1.02)
    fig.tight_layout()
    savefig(fig, outdir / f"sjc_alt_end_panel_b_{end_label}")


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bed", nargs="+", required=True,
                        help="label:path pairs for BED12 isoform files")
    parser.add_argument("--read-map", nargs="+", default=[],
                        help="label:path pairs for isoform read-map files")
    parser.add_argument("--cage-peaks", required=True, help="CAGE peaks BED6")
    parser.add_argument("--qs-peaks", required=True, help="QuantSeq peaks BED6")
    parser.add_argument("--cage-plus", required=True, help="CAGE bedGraph (+ strand)")
    parser.add_argument("--cage-minus", required=True, help="CAGE bedGraph (- strand)")
    parser.add_argument("--qs-plus", required=True, help="QuantSeq bedGraph (+ strand)")
    parser.add_argument("--qs-minus", required=True, help="QuantSeq bedGraph (- strand)")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Parse inputs ────────────────────────────────────────────────────
    beds_by_method: dict[str, list[dict]] = {}
    for entry in args.bed:
        label, path = entry.split(":", 1)
        beds_by_method[label] = parse_isoforms(path)
        log.info("Loaded %d isoforms for %s", len(beds_by_method[label]), label)

    read_maps_by_method: dict[str, dict[str, int]] = {}
    for entry in args.read_map:
        label, path = entry.split(":", 1)
        rm = load_read_map(path)
        read_maps_by_method[label] = rm
        log.info("Loaded read-map for %s: %d entries%s",
                 label, len(rm), " (self-ref, skipped)" if not rm else "")

    cage_peaks = _parse_peaks_bed(args.cage_peaks)
    qs_peaks = _parse_peaks_bed(args.qs_peaks)
    log.info("Loaded %d / %d CAGE / QuantSeq peak groups",
             len(cage_peaks), len(qs_peaks))

    cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
        args.cage_plus, args.cage_minus, args.qs_plus, args.qs_minus,
    )

    # ── Run analysis per end type ───────────────────────────────────────
    for end_type, end_label, peaks in [
        ("tss", "5prime", cage_peaks),
        ("tts", "3prime", qs_peaks),
    ]:
        results: dict[str, dict | None] = {}
        for label, isoforms in beds_by_method.items():
            rc = read_maps_by_method.get(label, {})
            r = _analyse_assembler(
                label, isoforms, rc, peaks,
                cage_p, cage_m, qs_p, qs_m,
                end_type,
            )
            total = r["tp_unique"] + r["tp_redundant"] + r["fp"]
            if total == 0:
                log.info("%s %s: no multi-end SJC groups, skipping", label, end_label)
                results[label] = None
            else:
                log.info("%s %s: %d alt ends (TP-unique=%d, TP-redundant=%d, FP=%d)",
                         label, end_label, total,
                         r["tp_unique"], r["tp_redundant"], r["fp"])
                results[label] = r

        plot_panel_a(results, end_label, outdir)
        plot_panel_b(results, end_label, outdir)

    log.info("Done — output in %s", outdir)


if __name__ == "__main__":
    main()
