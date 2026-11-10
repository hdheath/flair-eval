#!/usr/bin/env python3
"""
sqanti_category_precise_bar.py — Fig C: per-SQANTI-category bar, single tool.

For ONE tool, classify every isoform into a SQANTI structural category
(FSM, ISM, NIC, NNC, SEM, SEN) via
flair_structural.classify_transcripts_per_isoform, then mark each isoform as
"paired-precise" using the SAME interval-edge predicate as the headline paired
metric (ted_end_precision._nearest_annot against the CAGE/dRNA peak intervals):

    paired-precise  iff  tss_match is not None  AND  tts_match is not None

i.e. both the 5' end and the 3' end land within --window of (or inside) a truth
peak interval.  This is exactly the per-isoform (tss_match, tts_match) used by
compute_jc_deduplicated_precision_recall's pair logic, so the segment heights
reconcile with the paired-dedup precision headline.

Output: one bar per category, height = #isoforms in that category, the FILLED
(darker) segment = #paired-precise isoforms; the unfilled remainder is the
not-both-hit isoforms.  Annotated with n and the paired-precise fraction.

Intended to be run once per tool (e.g. TED-hybrid-me1 and FLAIR-default get
SEPARATE plots).

Read-support is NOT filtered here (every emitted isoform is shown) unless
--min-support > 1 is given with --counts.

Usage:
    python sqanti_category_precise_bar.py \\
        --bed isoforms.bed --counts counts.txt \\
        --cage CAGE.bed --drna dRNA.bed --gtf REF.gtf \\
        --label TED-hybrid-me1 --outpng out/figC_ted.png \\
        [--region chr22] [--window 50] [--min-support 1]
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

sys.path.insert(0, str(Path(__file__).resolve().parent))

try:
    from signal_utils import tss_tts, lookup_read_count, load_support_table
except ImportError:
    from evaluation.signal_utils import (
        tss_tts, lookup_read_count, load_support_table)
from flair_structural import (
    classify_transcripts_per_isoform,
    build_reference_structures,
    parse_gtf_transcripts as parse_gtf_transcripts_nested,
)
from ted_end_precision import (
    parse_peaks_bed,
    parse_region_values,
    _nearest_annot,
)
try:
    from peak_matching import is_canonical_chrom
except ImportError:
    from evaluation.peak_matching import is_canonical_chrom

CATEGORY_ORDER = ["FSM", "ISM", "NIC", "NNC", "SEM", "SEN"]
CATEGORY_COLORS = {
    "FSM": "#009E73", "ISM": "#56B4E9", "NIC": "#E69F00",
    "NNC": "#D55E00", "SEM": "#CC79A7", "SEN": "#999999",
}


def _render_bar(label, cat_total, cat_precise, out):
    """Draw the per-category bar (total light + paired-precise solid).

    Shared by the compute path and --replot so the redraw is identical.
    cat_total / cat_precise are {category: int}.
    """
    cats = [c for c in CATEGORY_ORDER if cat_total.get(c, 0) > 0]
    x = range(len(cats))
    totals = [cat_total[c] for c in cats]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for i, c in enumerate(cats):
        col = CATEGORY_COLORS[c]
        ax.bar(i, cat_total[c], width=0.72, color=col, alpha=0.30,
               edgecolor=col, linewidth=0.6)
        ax.bar(i, cat_precise.get(c, 0), width=0.72, color=col, alpha=0.95,
               edgecolor="none")
    for i, c in enumerate(cats):
        nt = cat_total[c]; npp = cat_precise.get(c, 0)
        frac = npp / nt if nt else 0
        ax.text(i, nt, f"n={nt}\n{frac:.0%}", ha="center", va="bottom",
                fontsize=7, linespacing=1.2)

    ax.set_xticks(list(x))
    ax.set_xticklabels(cats, fontsize=9)
    ax.set_ylabel("# isoforms (solid = paired-precise)", fontsize=10)
    ax.set_xlabel("SQANTI structural category", fontsize=10)
    ax.set_title(f"{label}: paired-precise isoforms by SQANTI category", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ymax = max(totals) if totals else 1
    ax.set_ylim(0, ymax * 1.18)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def _read_bar_tsv(tsv):
    """Read <outpng>.tsv -> (label, cat_total, cat_precise) for --replot."""
    cat_total, cat_precise = {}, {}
    label = ""
    with open(tsv) as fh:
        header = fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 4:
                continue
            label = c[0]
            cat = c[1]
            try:
                cat_total[cat] = int(c[2])
                cat_precise[cat] = int(c[3])
            except ValueError:
                continue
    return label, cat_total, cat_precise


def _parse_bed_geometry(bed_path, regions, canonical_only):
    """BED12 -> {name: {chrom,strand,tss,tts,start,end}}, region/canonical filtered."""
    geom = {}
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("track"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            chrom, start, end = p[0], int(p[1]), int(p[2])
            name, strand = p[3], p[5]
            if canonical_only and not is_canonical_chrom(chrom):
                continue
            if regions:
                keep = False
                for rc, rs, re_ in regions:
                    if chrom != rc:
                        continue
                    if rs is None or re_ is None or not (end < rs or start > re_):
                        keep = True
                        break
                if not keep:
                    continue
            tss, tts = tss_tts(start, end, strand)
            geom[name] = {"chrom": chrom, "strand": strand,
                          "tss": tss, "tts": tts, "start": start, "end": end}
    return geom


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", help="BED12 isoforms for ONE tool")
    ap.add_argument("--counts", default=None, help="per-transcript read-count table")
    ap.add_argument("--cage", help="CAGE peaks BED6 (5')")
    ap.add_argument("--drna", help="dRNA peaks BED6 (3')")
    ap.add_argument("--gtf", help="reference GTF for SQANTI classification")
    ap.add_argument("--label", default="")
    ap.add_argument("--outpng", required=True)
    ap.add_argument("--region", nargs="+", default=None)
    ap.add_argument("--window", type=int, default=50)
    ap.add_argument("--min-support", type=int, default=1)
    ap.add_argument("--all-contigs", action="store_true")
    ap.add_argument("--replot", action="store_true",
                    help="re-draw the bar from the existing <outpng>.tsv checkpoint "
                         "with ZERO recompute (no classification/peak parsing).")
    a = ap.parse_args()

    out = Path(a.outpng)
    out.parent.mkdir(parents=True, exist_ok=True)

    # ── REPLOT: redraw the bar from the existing <outpng>.tsv ───────────────
    if a.replot:
        tsv = out.with_suffix(".tsv")
        if not tsv.exists():
            sys.exit(f"--replot: checkpoint not found: {tsv}")
        print(f"[replot] reading checkpoint {tsv} (no recompute)", flush=True)
        label, cat_total, cat_precise = _read_bar_tsv(tsv)
        _render_bar(label or a.label, cat_total, cat_precise, out)
        print("SQANTI_CATEGORY_PRECISE_DONE")
        return

    if not (a.bed and a.cage and a.drna and a.gtf and a.label):
        sys.exit("compute mode requires --bed --cage --drna --gtf --label "
                 "(or pass --replot to redraw from <outpng>.tsv)")
    canonical_only = not a.all_contigs
    regions = parse_region_values(a.region)

    print(f"[{a.label}] parsing reference GTF for classification ...", flush=True)
    transcripttoexons = parse_gtf_transcripts_nested(a.gtf)
    refjuncs, refjuncchains, refseends = build_reference_structures(transcripttoexons)

    print(f"[{a.label}] classifying isoforms ...", flush=True)
    classified = classify_transcripts_per_isoform(
        a.bed, refjuncs, refjuncchains, refseends)
    name_to_cat = {c["name"]: c["category"] for c in classified}

    geom = _parse_bed_geometry(a.bed, regions, canonical_only)

    counts = load_support_table(a.counts)

    print(f"[{a.label}] loading truth peaks ...", flush=True)
    peaks_5 = parse_peaks_bed(a.cage, regions, canonical_only=canonical_only)
    peaks_3 = parse_peaks_bed(a.drna, regions, canonical_only=canonical_only)

    # per-category counts: total + paired-precise
    cat_total = defaultdict(int)
    cat_precise = defaultdict(int)
    n_kept = 0
    for name, g in geom.items():
        cat = name_to_cat.get(name)
        if cat not in CATEGORY_ORDER:
            continue
        # support filter (optional)
        if a.min_support > 1:
            c = lookup_read_count(name, counts) if counts else 0
            if c < a.min_support:
                continue
        n_kept += 1
        cat_total[cat] += 1
        key = (g["chrom"], g["strand"])
        p5 = peaks_5.get(key, [])
        p3 = peaks_3.get(key, [])
        tss_match = _nearest_annot(g["tss"], p5, a.window) if p5 else None
        tts_match = _nearest_annot(g["tts"], p3, a.window) if p3 else None
        if tss_match is not None and tts_match is not None:
            cat_precise[cat] += 1

    print(f"[{a.label}] {n_kept} isoforms; per-category:")
    for cat in CATEGORY_ORDER:
        if cat_total[cat]:
            frac = cat_precise[cat] / cat_total[cat]
            print(f"    {cat:4s}  n={cat_total[cat]:6d}  precise={cat_precise[cat]:6d}  "
                  f"({frac:.1%})")

    # write tidy TSV alongside the PNG
    tsv = out.with_suffix(".tsv")
    with open(tsv, "w") as fh:
        fh.write("label\tcategory\tn_isoforms\tn_paired_precise\tprecise_fraction\n")
        for cat in CATEGORY_ORDER:
            nt = cat_total[cat]
            npp = cat_precise[cat]
            fr = (npp / nt) if nt else ""
            fh.write(f"{a.label}\t{cat}\t{nt}\t{npp}\t"
                     f"{fr if fr == '' else f'{fr:.6f}'}\n")
    print("wrote", tsv)

    # plot: one bar per category; total (light) + paired-precise (solid)
    _render_bar(a.label, cat_total, cat_precise, out)
    print("SQANTI_CATEGORY_PRECISE_DONE")


if __name__ == "__main__":
    main()
