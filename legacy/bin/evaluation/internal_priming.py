#!/usr/bin/env python3
"""
internal_priming.py — Internal priming analysis at isoform TTS positions.

For each isoform, computes the A-content in a window downstream of the TTS
using the reference genome sequence as a proxy for internal priming likelihood.
Also computes distance to the nearest annotated polyA site from a GTF.

Produces three plots:
  A1. internal_priming_summary.png
      Fraction of isoforms with high A-content (>= threshold) at TTS,
      stratified by read-support bins. Per-method bars, all methods on same axes.

  A2. internal_priming_cross_tool.png
      Same A-content metric across ALL tools (FLAIR, FLAIR-permissive,
      IsoQuant, TED modes) — single grouped bar chart for direct comparison.

  A3. proximal_apa_distance_support.png
      Scatter: distance from TTS to nearest annotated polyA site (y-axis)
      vs read support (x-axis), coloured by within-50bp (likely real) vs
      distant (likely internal priming). Faceted by method.

Usage:
    python internal_priming.py \\
        --bed label1:bed1.bed label2:bed2.bed ... \\
        --read-map label1:map1.txt label2:map2.txt ... \\
        --genome genome.fa \\
        --gtf annotation.gtf \\
        --output output_dir/ \\
        [--window 30] \\
        [--a-threshold 0.60] \\
        [--apa-window 50]
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
    from signal_utils import parse_isoforms, tss_tts, load_read_map, lookup_read_count
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
    from evaluation.signal_utils import parse_isoforms, tss_tts, load_read_map, lookup_read_count

apply_rc()

# ── Constants ────────────────────────────────────────────────────────────────

DEFAULT_WINDOW    = 30    # bp downstream of TTS to scan for A-content
DEFAULT_A_THRESH  = 0.60  # A-fraction threshold for "likely internal priming"
DEFAULT_APA_WIN   = 50    # bp window for "proximal to polyA site"
SUPPORT_BINS      = [1, 2, 5, 10, 20, float("inf")]
SUPPORT_LABELS    = ["1", "2–4", "5–9", "10–19", "20+"]

# ── Motif-concordance-by-confidence-bin (Anvar An-F3 analog) ──────────────────
# Canonical poly(A)-signal hexamers (AATAAA-class). AATAAA + ATTAAA cover the
# large majority of mammalian PAS; the remainder are the well-characterised
# single-variant set (Beaudoing 2000 / Gruber-Zavolan 2019).
PAS_HEXAMERS = [
    "AATAAA", "ATTAAA", "AGTAAA", "TATAAA", "CATAAA", "GATAAA",
    "AATATA", "AATACA", "AATAGA", "AAAAAG", "ACTAAA",
]
# 5' TSS initiator/core-promoter signatures: Inr (consensus BBCABW, we use the
# diagnostic CA dinucleotide at the +1 site) and the TATA box upstream.
TSS_HEXAMERS_INR = ["TCAGTT", "TCAGTC", "CCAGTT", "TCATTC", "CCATTT"]  # Inr-like CA core
TATA_BOX = ["TATAAA", "TATAAT", "TATATA", "TATAWA"]
# TED-reality confidence bins (left-closed). reality in [0,1].
REALITY_BINS   = [0.0, 0.25, 0.5, 0.75, 1.0001]
REALITY_LABELS = ["[0,0.25)", "[0.25,0.5)", "[0.5,0.75)", "[0.75,1]"]
DEFAULT_PAS_HEX_WINDOW = 50   # bp upstream of PAS to scan for the hexamer
DEFAULT_TSS_MOTIF_WINDOW = 50  # bp around TSS for Inr/TATA


# ── Genome FASTA ─────────────────────────────────────────────────────────────

def load_fasta(path: str | Path) -> Dict[str, str]:
    """Load a FASTA file into {chrom: sequence}."""
    seqs: Dict[str, str] = {}
    chrom = None
    buf: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if chrom is not None:
                    seqs[chrom] = "".join(buf).upper()
                chrom = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if chrom is not None:
        seqs[chrom] = "".join(buf).upper()
    return seqs


def a_content_downstream(chrom: str, tts: int, strand: str,
                          genome: Dict[str, str], window: int) -> Optional[float]:
    """Return fraction of A bases in a window downstream of the TTS.

    Downstream means: 3' of the TTS in the direction of transcription.
      + strand: sequence [tts, tts+window)
      - strand: sequence [tts-window, tts), then reverse-complement → A becomes T
    """
    seq = genome.get(chrom)
    if seq is None:
        return None
    if strand == "+":
        region = seq[max(0, tts):min(len(seq), tts + window)]
        return region.count("A") / len(region) if region else None
    else:
        region = seq[max(0, tts - window):min(len(seq), tts)]
        # downstream on minus strand → look for T (complement of A)
        return region.count("T") / len(region) if region else None


# ── GTF polyA site extraction ─────────────────────────────────────────────────

def extract_polya_sites(gtf_path: str | Path) -> Dict[str, List[int]]:
    """Extract annotated TTS positions per chrom from a GTF as polyA site proxies."""
    sites: Dict[str, List[int]] = defaultdict(list)
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip().split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            chrom  = cols[0]
            start  = int(cols[3]) - 1
            end    = int(cols[4])
            strand = cols[6]
            _, tts = tss_tts(start, end, strand)
            sites[chrom].append(tts)
    # sort each list
    for ch in sites:
        sites[ch] = sorted(set(sites[ch]))
    return sites


def nearest_apa_distance(chrom: str, tts: int,
                          polya_sites: Dict[str, List[int]]) -> Optional[int]:
    """Return distance (bp) from tts to nearest annotated polyA site."""
    import bisect
    positions = polya_sites.get(chrom)
    if not positions:
        return None
    idx = bisect.bisect_left(positions, tts)
    candidates = []
    if idx < len(positions):
        candidates.append(abs(positions[idx] - tts))
    if idx > 0:
        candidates.append(abs(positions[idx - 1] - tts))
    return min(candidates) if candidates else None


# ── Support bin assignment ───────────────────────────────────────────────────

def support_bin(n: int) -> int:
    for i, upper in enumerate(SUPPORT_BINS[1:]):
        if n < upper:
            return i
    return len(SUPPORT_BINS) - 2


# ── Per-isoform data building ─────────────────────────────────────────────────

def build_isoform_data(
    beds_by_method: Dict[str, List[dict]],
    read_maps:       Dict[str, Dict[str, int]],
    genome:          Dict[str, str],
    polya_sites:     Dict[str, List[int]],
    window:          int,
    a_thresh:        float,
) -> Dict[str, List[dict]]:
    """Return {method: [{reads, a_content, apa_dist, likely_ip}]}."""
    data: Dict[str, List[dict]] = {}
    for method, isoforms in beds_by_method.items():
        rc = read_maps.get(method, {})
        rows = []
        for iso in isoforms:
            reads  = lookup_read_count(iso["name"], rc) if rc else 0
            a_frac = a_content_downstream(
                iso["chrom"], iso["tts"], iso["strand"], genome, window)
            apa_d  = nearest_apa_distance(iso["chrom"], iso["tts"], polya_sites)
            if a_frac is None:
                continue
            rows.append({
                "reads":     reads,
                "a_content": a_frac,
                "apa_dist":  apa_d,
                "likely_ip": a_frac >= a_thresh,
                "bin":       support_bin(reads),
            })
        data[method] = rows
    return data


# ── Plot A1: internal priming fraction by read-support bin ───────────────────

def plot_internal_priming_summary(
    data: Dict[str, List[dict]],
    output_dir: Path,
    a_thresh: float,
):
    methods = list(data.keys())
    if not methods:
        return
    styler = ModeStyler(methods)
    n_bins = len(SUPPORT_LABELS)
    x      = np.arange(n_bins)
    width  = 0.8 / max(len(methods), 1)

    fig, ax = plt.subplots(figsize=(W2, W1 * 0.75))
    handles = []
    for i, m in enumerate(methods):
        rows = data[m]
        fracs = []
        for b in range(n_bins):
            bin_rows = [r for r in rows if r["bin"] == b]
            if bin_rows:
                fracs.append(sum(r["likely_ip"] for r in bin_rows) / len(bin_rows))
            else:
                fracs.append(0.0)
        offset = (i - len(methods) / 2 + 0.5) * width
        bars = ax.bar(x + offset, fracs, width * 0.9,
                      color=styler.color(m), edgecolor="none", alpha=0.85)
        handles.append(styler.legend_handle(m, label=m, markersize=5))

    ax.set_xticks(x)
    ax.set_xticklabels(SUPPORT_LABELS, fontsize=7)
    ax.set_ylim(0, 1)
    style_ax(ax,
             xlabel="Read support (reads per isoform)",
             ylabel=f"Fraction with TTS A-content ≥ {a_thresh:.0%}")
    ax.axhline(a_thresh, color="0.5", lw=0.5, ls="--", zorder=0)
    legend_outside(fig, handles=handles, loc="upper right",
                   bbox_to_anchor=(1.0, 1.0), ncol=1, fontsize=6)
    fig.tight_layout(pad=0.3)
    savefig(fig, output_dir / "internal_priming_summary.png")


# ── Plot A2: cross-tool single grouped bar ───────────────────────────────────

def plot_cross_tool_comparison(
    data: Dict[str, List[dict]],
    output_dir: Path,
    a_thresh: float,
):
    methods = list(data.keys())
    if not methods:
        return
    styler  = ModeStyler(methods)
    fracs   = [
        (sum(r["likely_ip"] for r in data[m]) / len(data[m]))
        if data[m] else 0.0
        for m in methods
    ]
    x = np.arange(len(methods))

    fig, ax = plt.subplots(figsize=(W2, W1 * 0.65))
    colors  = [styler.color(m) for m in methods]
    ax.bar(x, fracs, 0.7, color=colors, edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=5, rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.axhline(a_thresh, color="0.5", lw=0.5, ls="--", zorder=0)
    style_ax(ax, ylabel=f"Fraction with TTS A-content ≥ {a_thresh:.0%}")
    fig.tight_layout(pad=0.3)
    savefig(fig, output_dir / "internal_priming_cross_tool.png")


# ── Plot A3: APA distance vs read support ────────────────────────────────────

def plot_proximal_apa(
    data: Dict[str, List[dict]],
    output_dir: Path,
    apa_window: int,
):
    methods = [m for m in data if any(r["apa_dist"] is not None for r in data[m])]
    if not methods:
        return

    styler = ModeStyler(methods)
    n = len(methods)
    ncols = min(n, 4)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(W2, W1 * 0.75 * nrows / max(nrows, 1)),
        squeeze=False,
    )

    for idx, m in enumerate(methods):
        ax   = axes[idx // ncols][idx % ncols]
        rows = [r for r in data[m] if r["apa_dist"] is not None and r["reads"] > 0]
        if not rows:
            ax.set_visible(False)
            continue

        xs       = np.array([r["reads"]   for r in rows], dtype=float)
        ys       = np.array([r["apa_dist"] for r in rows], dtype=float)
        proximal = ys <= apa_window

        ax.scatter(xs[~proximal], ys[~proximal],
                   s=1.5, alpha=0.25, color="#CC79A7",
                   edgecolors="none", rasterized=True, label="Distant")
        ax.scatter(xs[proximal], ys[proximal],
                   s=1.5, alpha=0.4, color="#009E73",
                   edgecolors="none", rasterized=True, label="Proximal")
        ax.axhline(apa_window, color="0.5", lw=0.5, ls="--")
        ax.set_xscale("log")
        style_ax(ax)
        ax.text(0.04, 0.97, m, transform=ax.transAxes,
                ha="left", va="top", fontsize=5, fontweight="bold")
        pct = proximal.sum() / len(proximal) * 100 if len(proximal) else 0
        ax.text(0.04, 0.88, f"{pct:.0f}% proximal",
                transform=ax.transAxes, ha="left", va="top", fontsize=5,
                color="#009E73")

        if idx >= n - ncols:
            ax.set_xlabel("Read support", fontsize=6)
        if idx % ncols == 0:
            ax.set_ylabel("Distance to nearest polyA site (bp)", fontsize=6)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    # Shared legend
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#009E73",
               markersize=5, label=f"≤ {apa_window} bp (proximal)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#CC79A7",
               markersize=5, label=f"> {apa_window} bp (likely IP)"),
    ]
    fig.legend(handles=handles, loc="lower right", fontsize=6, frameon=False,
               bbox_to_anchor=(1.0, 0.0))
    fig.tight_layout(pad=0.3)
    savefig(fig, output_dir / "proximal_apa_distance_support.png")


# ── Motif-concordance-by-confidence-bin (Rank 1, Anvar An-F3 analog) ──────────

def _seq_upstream_of_pas(chrom: str, tts: int, strand: str,
                          genome: Dict[str, str], window: int) -> Optional[str]:
    """Return the sequence in *window* bp ending at the PAS, in txn orientation.

    The canonical AATAAA hexamer sits ~10–30 bp 5' (upstream) of the cleavage
    site, so we scan upstream of the TTS in the direction of transcription.
      + strand: genomic [tts-window, tts)
      - strand: genomic (tts, tts+window], reverse-complemented
    """
    seq = genome.get(chrom)
    if seq is None:
        return None
    if strand == "+":
        region = seq[max(0, tts - window):min(len(seq), tts)]
        return region or None
    else:
        region = seq[max(0, tts):min(len(seq), tts + window)]
        if not region:
            return None
        comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        return "".join(comp.get(b, "N") for b in reversed(region))


def _seq_around_tss(chrom: str, tss: int, strand: str,
                     genome: Dict[str, str], window: int) -> Optional[str]:
    """Return sequence ±window around the TSS in txn orientation."""
    seq = genome.get(chrom)
    if seq is None:
        return None
    if strand == "+":
        region = seq[max(0, tss - window):min(len(seq), tss + window)]
        return region or None
    else:
        region = seq[max(0, tss - window):min(len(seq), tss + window)]
        if not region:
            return None
        comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        return "".join(comp.get(b, "N") for b in reversed(region))


def _matches_any(seq: str, motifs: List[str]) -> bool:
    """True if *seq* contains any motif. 'W' = A/T IUPAC handled."""
    if not seq:
        return False
    for m in motifs:
        if "W" in m:
            import re as _re
            pat = m.replace("W", "[AT]")
            if _re.search(pat, seq):
                return True
        elif m in seq:
            return True
    return False


def has_pas_hexamer(chrom: str, tts: int, strand: str,
                    genome: Dict[str, str], window: int) -> Optional[bool]:
    seq = _seq_upstream_of_pas(chrom, tts, strand, genome, window)
    if seq is None:
        return None
    return _matches_any(seq, PAS_HEXAMERS)


def has_tss_motif(chrom: str, tss: int, strand: str,
                  genome: Dict[str, str], window: int) -> Optional[bool]:
    seq = _seq_around_tss(chrom, tss, strand, genome, window)
    if seq is None:
        return None
    return _matches_any(seq, TSS_HEXAMERS_INR + TATA_BOX)


def reality_bin_index(reality: float) -> int:
    for i in range(len(REALITY_BINS) - 1):
        if REALITY_BINS[i] <= reality < REALITY_BINS[i + 1]:
            return i
    return len(REALITY_BINS) - 2


def build_motif_concordance(
    isoforms: List[dict],
    reality_index: dict,
    genome: Dict[str, str],
    pas_window: int,
    tss_window: int,
    a_window: int,
    a_thresh: float,
    end_tol: int,
) -> List[dict]:
    """Per-isoform-end concordance rows, tagged with TED reality bin.

    For each isoform end matched to a TED-reality score (via the ted_log),
    record: the reality bin, whether the 3' end carries a canonical PAS hexamer,
    whether it has the A-rich internal-priming signature, and whether the 5'
    end carries a TSS (Inr/TATA) motif. Ends with no reality match are dropped
    (they have no confidence to bin by).
    """
    try:
        from alt_end_enumeration import reality_at
    except ImportError:
        from evaluation.alt_end_enumeration import reality_at

    rows: List[dict] = []
    for iso in isoforms:
        chrom = iso["chrom"]
        # 3' end
        rty_t = reality_at("tts", chrom, int(iso["tts"]), reality_index, end_tol)
        if rty_t is not None:
            hexp = has_pas_hexamer(chrom, int(iso["tts"]), iso["strand"],
                                   genome, pas_window)
            a_frac = a_content_downstream(chrom, int(iso["tts"]), iso["strand"],
                                          genome, a_window)
            rows.append({
                "axis": "tts", "reality": rty_t,
                "bin": reality_bin_index(rty_t),
                "has_motif": hexp,
                "a_rich": (a_frac is not None and a_frac >= a_thresh),
            })
        # 5' end
        rty_s = reality_at("tss", chrom, int(iso["tss"]), reality_index, end_tol)
        if rty_s is not None:
            mot = has_tss_motif(chrom, int(iso["tss"]), iso["strand"],
                                genome, tss_window)
            rows.append({
                "axis": "tss", "reality": rty_s,
                "bin": reality_bin_index(rty_s),
                "has_motif": mot, "a_rich": None,
            })
    return rows


def write_concordance_tsv(rows_by_method: Dict[str, List[dict]], path: Path):
    """Per-(method, axis, reality bin) motif-concordance summary."""
    import csv as _csv
    cols = ["method", "axis", "reality_bin", "n_ends",
            "frac_with_motif", "frac_a_rich"]
    with open(path, "w", newline="") as f:
        w = _csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for method, rows in rows_by_method.items():
            for axis in ("tss", "tts"):
                for b, label in enumerate(REALITY_LABELS):
                    sub = [r for r in rows if r["axis"] == axis and r["bin"] == b]
                    n = len(sub)
                    motif_vals = [r["has_motif"] for r in sub
                                  if r["has_motif"] is not None]
                    a_vals = [r["a_rich"] for r in sub if r["a_rich"] is not None]
                    frac_motif = (sum(motif_vals) / len(motif_vals)
                                  if motif_vals else "")
                    frac_a = (sum(a_vals) / len(a_vals)) if a_vals else ""
                    w.writerow([
                        method, axis, label, n,
                        round(frac_motif, 5) if motif_vals else "",
                        round(frac_a, 5) if a_vals else "",
                    ])


def plot_motif_concordance(rows_by_method: Dict[str, List[dict]],
                           output_dir: Path, title_prefix: str = ""):
    methods = [m for m in rows_by_method if rows_by_method[m]]
    if not methods:
        return
    styler = ModeStyler(methods)
    x = np.arange(len(REALITY_LABELS))

    fig, axes = plt.subplots(1, 2, figsize=(W2, W1 * 0.7))

    # PAS hexamer fraction (and A-rich) vs reality bin
    ax = axes[0]
    for m in methods:
        rows = [r for r in rows_by_method[m] if r["axis"] == "tts"]
        motif_fracs, a_fracs = [], []
        for b in range(len(REALITY_LABELS)):
            sub = [r for r in rows if r["bin"] == b]
            mv = [r["has_motif"] for r in sub if r["has_motif"] is not None]
            av = [r["a_rich"] for r in sub if r["a_rich"] is not None]
            motif_fracs.append(sum(mv) / len(mv) if mv else np.nan)
            a_fracs.append(sum(av) / len(av) if av else np.nan)
        ax.plot(x, motif_fracs, "-o", color=styler.color(m), ms=3, lw=1.2,
                label=f"{m} AATAAA")
        ax.plot(x, a_fracs, "--s", color=styler.color(m), ms=2.5, lw=0.9,
                alpha=0.6, label=f"{m} A-rich")
    ax.set_xticks(x); ax.set_xticklabels(REALITY_LABELS, fontsize=5, rotation=30)
    ax.set_ylim(0, 1)
    style_ax(ax, xlabel="TED 3' reality bin",
             ylabel="fraction of 3' ends")
    ax.legend(fontsize=4.5, frameon=False, loc="best")

    # TSS motif fraction vs reality bin
    ax = axes[1]
    for m in methods:
        rows = [r for r in rows_by_method[m] if r["axis"] == "tss"]
        fr = []
        for b in range(len(REALITY_LABELS)):
            sub = [r for r in rows if r["bin"] == b]
            mv = [r["has_motif"] for r in sub if r["has_motif"] is not None]
            fr.append(sum(mv) / len(mv) if mv else np.nan)
        ax.plot(x, fr, "-o", color=styler.color(m), ms=3, lw=1.2, label=m)
    ax.set_xticks(x); ax.set_xticklabels(REALITY_LABELS, fontsize=5, rotation=30)
    ax.set_ylim(0, 1)
    style_ax(ax, xlabel="TED 5' reality bin",
             ylabel="fraction of 5' ends with Inr/TATA")
    ax.legend(fontsize=5, frameon=False, loc="best")

    if title_prefix:
        fig.suptitle(f"{title_prefix}motif concordance by TED confidence",
                     fontsize=8)
    fig.tight_layout(pad=0.4)
    savefig(fig, output_dir / "motif_concordance_by_confidence.png")


def run_motif_concordance(args) -> int:
    """CLI entry for --by-confidence mode. Degrades gracefully on empty input."""
    try:
        from alt_end_enumeration import load_end_reality
    except ImportError:
        from evaluation.alt_end_enumeration import load_end_reality

    bed_pairs, log_pairs = {}, {}
    for e in args.bed:
        if ":" in e:
            k, v = e.split(":", 1); bed_pairs[k] = v
    for e in args.ted_log:
        if ":" in e:
            k, v = e.split(":", 1); log_pairs[k] = v

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = output_dir / "motif_concordance_by_confidence.tsv"

    genome = load_fasta(args.genome) if Path(args.genome).exists() else {}

    rows_by_method: Dict[str, List[dict]] = {}
    for method, bed_path in bed_pairs.items():
        if not Path(bed_path).exists():
            rows_by_method[method] = []
            continue
        isos = parse_isoforms(bed_path)
        lp = log_pairs.get(method)
        reality_index = (load_end_reality(Path(lp))
                         if lp and lp.lower() != "none" and Path(lp).exists()
                         else {})
        if not reality_index:
            print(f"WARNING: no TED reality for {method}; "
                  "motif-concordance-by-confidence needs a ted_log — skipping "
                  "this method.", file=sys.stderr)
            rows_by_method[method] = []
            continue
        rows_by_method[method] = build_motif_concordance(
            isos, reality_index, genome,
            args.pas_hex_window, args.tss_motif_window,
            args.window, args.a_threshold, args.apa_window)

    # Always write the TSV (header-only if empty) — never exit non-zero.
    write_concordance_tsv(rows_by_method, tsv_path)
    print(f"Wrote {tsv_path}")
    if any(rows_by_method.values()):
        plot_motif_concordance(rows_by_method, output_dir, args.title_prefix)
        print(f"Wrote {output_dir / 'motif_concordance_by_confidence.png'}")
    else:
        print("No reality-matched ends — skipping plot "
              "(empty input handled gracefully).", file=sys.stderr)
    return 0


# ── CLI ──────────────────────────────────────────────────────────────────────

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bed",       nargs="+", required=True,
                        help="label:path pairs for BED12 isoform files")
    parser.add_argument("--read-map",  nargs="+", default=[],
                        help="label:path pairs for isoform read-map files")
    parser.add_argument("--genome",    required=True,
                        help="Reference genome FASTA")
    parser.add_argument("--gtf",       default=None,
                        help="GTF annotation (for polyA site proxies; "
                             "required unless --by-confidence)")
    parser.add_argument("--output",    required=True,
                        help="Output directory")
    parser.add_argument("--window",    type=int, default=DEFAULT_WINDOW,
                        help=f"bp downstream of TTS for A-content scan (default: {DEFAULT_WINDOW})")
    parser.add_argument("--a-threshold", type=float, default=DEFAULT_A_THRESH,
                        help=f"A-fraction threshold for internal priming call (default: {DEFAULT_A_THRESH})")
    parser.add_argument("--apa-window", type=int, default=DEFAULT_APA_WIN,
                        help=f"bp for proximal polyA classification (default: {DEFAULT_APA_WIN})")
    # ── Rank-1 motif-concordance-by-confidence-bin mode ──
    parser.add_argument("--by-confidence", action="store_true",
                        help="Run motif-concordance-by-TED-confidence-bin "
                             "(Anvar An-F3 analog) instead of the IP plots. "
                             "Requires --ted-log.")
    parser.add_argument("--ted-log", nargs="+", default=[],
                        help="label:path pairs for flair_ted.ted_log TSVs "
                             "(for --by-confidence)")
    parser.add_argument("--pas-hex-window", type=int, default=DEFAULT_PAS_HEX_WINDOW,
                        help="bp upstream of PAS to scan for AATAAA-class hexamer")
    parser.add_argument("--tss-motif-window", type=int, default=DEFAULT_TSS_MOTIF_WINDOW,
                        help="bp around TSS to scan for Inr/TATA")
    parser.add_argument("--title-prefix", default="")
    parser.add_argument("--verbose",   action="store_true")
    args = parser.parse_args(argv)

    if args.by_confidence:
        return run_motif_concordance(args)

    if not args.gtf:
        parser.error("--gtf is required (unless --by-confidence)")

    # --- Load BED files ---
    beds_by_method: Dict[str, List[dict]] = {}
    for entry in args.bed:
        if ":" not in entry:
            print(f"WARNING: skipping malformed --bed entry '{entry}'", file=sys.stderr)
            continue
        label, path = entry.split(":", 1)
        if not Path(path).exists():
            print(f"WARNING: file not found: {path}", file=sys.stderr)
            continue
        isos = parse_isoforms(path)
        if isos:
            beds_by_method[label] = isos
            if args.verbose:
                print(f"  {label}: {len(isos)} isoforms", file=sys.stderr)

    if not beds_by_method:
        print("No isoform data loaded — skipping", file=sys.stderr)
        sys.exit(1)

    # --- Load read maps ---
    read_maps: Dict[str, Dict[str, int]] = {}
    for entry in args.read_map:
        if ":" not in entry:
            continue
        label, path = entry.split(":", 1)
        if not Path(path).exists():
            continue
        rc = load_read_map(path)
        if rc:
            read_maps[label] = rc

    # --- Load genome ---
    if args.verbose:
        print("  Loading genome FASTA...", file=sys.stderr)
    genome = load_fasta(args.genome)

    # --- Load polyA sites from GTF ---
    if args.verbose:
        print("  Extracting polyA sites from GTF...", file=sys.stderr)
    polya_sites = extract_polya_sites(args.gtf)

    # --- Build per-isoform data ---
    if args.verbose:
        print("  Computing A-content and APA distances...", file=sys.stderr)
    data = build_isoform_data(
        beds_by_method, read_maps, genome, polya_sites,
        args.window, args.a_threshold,
    )

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_internal_priming_summary(data, output_dir, args.a_threshold)
    plot_cross_tool_comparison(data, output_dir, args.a_threshold)
    plot_proximal_apa(data, output_dir, args.apa_window)

    print(f"Saved internal priming plots to {args.output}")


if __name__ == "__main__":
    sys.exit(main() or 0)
