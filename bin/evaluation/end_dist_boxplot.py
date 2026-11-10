#!/usr/bin/env python3
"""
end_dist_boxplot.py — Fig D: per-isoform end-to-nearest-ORTHOGONAL-PEAK distance.

For every isoform of every tool (strand-aware, distance to the peak INTERVAL —
0 if the end is inside a peak, else distance to the nearest peak edge):
  * TSS distance = isoform 5' TSS -> nearest CAGE peak
  * TTS distance = isoform 3' TTS -> nearest dRNA peak

Isoforms with read-support <= 1 are dropped (counts table; FLAIR/TED have inline
support via FlairBed col16 / BED col5).  One box per tool, two panels (5' TSS,
3' TTS), all tools overlaid.

This shows how tightly each tool's ends snap to the orthogonal end-support data
(a low, tight box = ends on CAGE/dRNA peaks; a long tail = ends with no nearby
orthogonal support).

Usage:
    python end_dist_boxplot.py \\
        --beds L1:bed1 L2:gtf2 ... \\
        --labels L1 L2 ...                 # (optional; derived from --beds tokens)
        --counts L1:counts1 L2:counts2 ... \\
        --cage CAGE.bed --drna DRNA.bed --outdir out/ \\
        [--region chr22] [--min-support 1]

(--beds entries may be plain paths aligned with --labels, OR LABEL:path tokens.)
"""
from __future__ import annotations

import argparse
import sys
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))

try:
    from signal_utils import (
        parse_isoforms, tss_tts, lookup_read_count, load_support_table,
        load_isoforms_any, is_gtf_path)
except ImportError:
    from evaluation.signal_utils import (
        parse_isoforms, tss_tts, lookup_read_count, load_support_table,
        load_isoforms_any, is_gtf_path)
from ted_end_precision import parse_gtf_ends, parse_region_values, parse_peaks_bed
try:
    from peak_matching import is_canonical_chrom
except ImportError:
    from evaluation.peak_matching import is_canonical_chrom

PALETTE = ["#1b1b6f", "#185fa5", "#1d9e75", "#63992f", "#b19a00", "#d4537e", "#993c1d"]


def _dist_to_interval(pos, lo, hi):
    """Distance from pos to half-open peak interval [lo, hi): 0 inside, else to
    the nearest edge (mirrors peak_matching._match_distance)."""
    if pos < lo:
        return lo - pos
    if pos >= hi:
        return pos - (hi - 1)
    return 0


def _build_peak_index(peaks):
    """{(chrom,strand):[(lo,hi)...]} -> {(chrom,strand):(starts, intervals, maxw)}."""
    idx = {}
    for k, ivs in peaks.items():
        ivs = sorted(ivs)
        starts = [lo for lo, _ in ivs]
        maxw = max((hi - lo for lo, hi in ivs), default=0)
        idx[k] = (starts, ivs, maxw)
    return idx


def _nearest_peak_dist(pos, key, idx):
    """Min distance from pos to the nearest peak interval on (chrom, strand).

    Returns None when that (chrom, strand) has no peaks (so the end contributes
    no point — e.g. an isoform on a strand with no CAGE/dRNA peaks).  Scans
    outward from the bisect point, bounded by the widest interval so a broad
    upstream peak that still contains/abuts `pos` is not missed.
    """
    entry = idx.get(key)
    if not entry:
        return None
    starts, ivs, maxw = entry
    i = bisect_left(starts, pos)
    best = None
    j = i                                   # rightward: intervals with lo >= pos
    while j < len(ivs):
        lo, hi = ivs[j]
        if best is not None and lo - pos > best:
            break
        d = _dist_to_interval(pos, lo, hi)
        best = d if best is None else min(best, d)
        if best == 0:
            return 0
        j += 1
    j = i - 1                               # leftward: intervals with lo < pos
    while j >= 0:
        lo, hi = ivs[j]
        d = _dist_to_interval(pos, lo, hi)
        best = d if best is None else min(best, d)
        if best == 0:
            return 0
        if best is not None and (pos - lo) - maxw > best:
            break
        j -= 1
    return best


def _support(iso, counts):
    name = iso.get("name", "")
    if counts:
        c = lookup_read_count(name, counts)
        if c:
            return float(c)
    rs = iso.get("read_support")
    if rs is not None:
        try:
            return float(rs)
        except (TypeError, ValueError):
            pass
    try:
        return float(iso.get("score", 0) or 0)
    except (TypeError, ValueError):
        return 0.0


def _in_regions(chrom, start, end, regions):
    if not regions:
        return True
    for rc, rs, re_ in regions:
        if chrom != rc:
            continue
        if rs is None or re_ is None or not (end < rs or start > re_):
            return True
    return False


def collect_distances(bed, counts, cage_idx, drna_idx, regions, min_support,
                      canonical_only):
    """Return (tss_dists, tts_dists, has_support) for one tool.

    Distances are to the ORTHOGONAL truth peaks (strand-aware): each isoform's
    5' TSS to the nearest CAGE peak interval, each 3' TTS to the nearest dRNA
    peak interval (0 if the end falls inside a peak; else distance to the peak
    edge — the same interval predicate as the precision metric).

    Auto-detects BED12 vs assembler GTF via the SHARED loader
    (signal_utils.load_isoforms_any).  has_support is False when the input is a
    GTF tool with no usable per-isoform read count (the caller decides whether
    the >min-support filter can be applied).
    """
    isos = load_isoforms_any(bed, regions=regions,
                             canonical_only=canonical_only, counts=counts)
    has_support = any(iso.get("has_support") for iso in isos) or not isos
    apply_filter = has_support and min_support >= 1
    tss_d, tts_d = [], []
    for iso in isos:
        if apply_filter:
            rs = iso.get("read_support")
            if rs is None or float(rs) <= min_support:
                continue
        key = (iso["chrom"], iso["strand"])
        dt = _nearest_peak_dist(iso["tss"], key, cage_idx)
        du = _nearest_peak_dist(iso["tts"], key, drna_idx)
        if dt is not None:
            tss_d.append(dt)
        if du is not None:
            tts_d.append(du)
    return tss_d, tts_d, has_support


def _parse_beds_labels(beds, labels):
    """Support both LABEL:path tokens and parallel --beds/--labels lists."""
    if all(":" in b and not Path(b).exists() for b in beds):
        # all LABEL:path
        out = []
        for b in beds:
            lab, path = b.split(":", 1)
            out.append((lab, path))
        return out
    assert labels and len(labels) == len(beds), \
        "provide --labels aligned with --beds (or use LABEL:path tokens)"
    return list(zip(labels, beds))


def _parse_counts_args(entries):
    out = {}
    for e in (entries or []):
        if ":" not in e:
            continue
        lab, path = e.split(":", 1)
        out[lab] = path
    return out


def _box_panel(ax, data_by_tool, labels, title):
    """Boxplot of distances per tool (log y; clip at 1 bp floor).

    Draws the box from the RAW per-isoform distance lists in data_by_tool, so
    the whiskers/quartiles come straight from the values (this is what makes the
    raw-distance checkpoint redraw a faithful box; a median/q25/q75 summary alone
    could not).
    """
    vals = []
    for lab in labels:
        arr = np.asarray(data_by_tool.get(lab, []), dtype=float)
        arr = np.maximum(arr, 1.0)  # log floor
        vals.append(arr if arr.size else np.array([1.0]))
    bp = ax.boxplot(vals, positions=range(len(labels)), widths=0.6,
                    showfliers=False, patch_artist=True,
                    medianprops=dict(color="black", linewidth=1.0))
    for i, box in enumerate(bp["boxes"]):
        box.set_facecolor(PALETTE[i % len(PALETTE)])
        box.set_alpha(0.7)
        box.set_edgecolor("black")
        box.set_linewidth(0.6)
    for i, lab in enumerate(labels):
        n = len(data_by_tool.get(lab, []))
        ax.text(i, ax.get_ylim()[1], f"n={n}", ha="center", va="top",
                fontsize=6, rotation=90)
    ax.set_yscale("log")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("dist to nearest peak (bp)", fontsize=9)
    ax.set_title(title, fontsize=11)
    ax.grid(axis="y", alpha=0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ── raw-distance checkpoint <-> figure (compute + --replot share this) ───────

_RAW_NAME = "end_dist_boxplot.raw.tsv"      # THE checkpoint (raw per-end dists)
_SUMMARY_NAME = "end_dist_boxplot.tsv"       # convenience summary (also written)
_PNG_NAME = "end_dist_boxplot.png"


def _write_raw_checkpoint(path, labels, min_support, tss_by_tool, tts_by_tool):
    """Write the RAW per-isoform-end distances (long form): label, end, distance.

    One row per isoform-end.  A header comment records the full ordered label
    list (incl. tools with zero rows) + min_support so --replot rebuilds the
    SAME panel ordering and title with no recompute.
    """
    tab = "\t"
    with open(path, "w") as fh:
        fh.write("# labels\t" + tab.join(labels) + "\n")
        fh.write(f"# min_support\t{min_support}\n")
        fh.write("label\tend\tdistance\n")
        for lab in labels:
            for d in tss_by_tool.get(lab, []):
                fh.write(f"{lab}\t5prime\t{float(d):.6g}\n")
            for d in tts_by_tool.get(lab, []):
                fh.write(f"{lab}\t3prime\t{float(d):.6g}\n")


def _read_raw_checkpoint(path):
    """Read the raw checkpoint -> (labels, min_support, tss_by_tool, tts_by_tool)."""
    labels = []
    min_support = 1
    tss_by_tool = defaultdict(list)
    tts_by_tool = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                parts = line[1:].strip().split("\t")
                if parts and parts[0] == "labels":
                    labels = parts[1:]
                elif parts and parts[0] == "min_support":
                    try:
                        min_support = int(parts[1])
                    except (IndexError, ValueError):
                        pass
                continue
            if line.startswith("label\t"):   # column header
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            lab, end, dist = c[0], c[1], c[2]
            try:
                d = float(dist)
            except ValueError:
                continue
            (tss_by_tool if end == "5prime" else tts_by_tool)[lab].append(d)
    if not labels:
        # fall back to first-seen order across both ends
        seen = []
        for lab in list(tss_by_tool) + list(tts_by_tool):
            if lab not in seen:
                seen.append(lab)
        labels = seen
    return labels, min_support, dict(tss_by_tool), dict(tts_by_tool)


def _write_summary(path, labels, tss_by_tool, tts_by_tool):
    """Convenience summary TSV (median/q25/q75/mean); NOT the checkpoint."""
    with open(path, "w") as fh:
        fh.write("label\tend\tn\tmedian\tq25\tq75\tmean\n")
        for lab in labels:
            for end, arr in (("5prime", tss_by_tool.get(lab, [])),
                             ("3prime", tts_by_tool.get(lab, []))):
                if arr:
                    a_ = np.asarray(arr, dtype=float)
                    fh.write(f"{lab}\t{end}\t{len(arr)}\t{float(np.median(a_))}\t"
                             f"{float(np.percentile(a_, 25))}\t"
                             f"{float(np.percentile(a_, 75))}\t"
                             f"{float(np.mean(a_))}\n")
                else:
                    fh.write(f"{lab}\t{end}\t0\t\t\t\t\n")


def _render(labels, min_support, tss_by_tool, tts_by_tool, png):
    """Draw the two-panel boxplot from raw distances (compute + replot share this)."""
    fig, ax = plt.subplots(1, 2, figsize=(max(8, 1.6 * len(labels)), 5))
    _box_panel(ax[0], tss_by_tool, labels, "5' TSS dist to nearest CAGE peak")
    _box_panel(ax[1], tts_by_tool, labels, "3' TTS dist to nearest dRNA peak")
    fig.suptitle("Per-isoform end distance to orthogonal peak (CAGE 5' / dRNA 3'; "
                 f"reads / isoform > {min_support})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", png)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--beds", nargs="+",
                    help="isoform BED12/GTF paths (or LABEL:path tokens)")
    ap.add_argument("--labels", nargs="+", default=None)
    ap.add_argument("--counts", nargs="+", default=[],
                    help="LABEL:path per-tool read-count tables")
    ap.add_argument("--cage", help="5' truth peaks BED6 (CAGE) — TSS distance target")
    ap.add_argument("--drna", help="3' truth peaks BED6 (dRNA) — TTS distance target")
    ap.add_argument("--gtf", help="(unused; kept for back-compat)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--region", nargs="+", default=None)
    ap.add_argument("--min-support", type=int, default=1,
                    help="drop isoforms with read-support <= this (default 1)")
    ap.add_argument("--all-contigs", action="store_true")
    ap.add_argument("--replot", action="store_true",
                    help="re-draw the PNG from the existing --outdir raw checkpoint "
                         "(end_dist_boxplot.raw.tsv) with ZERO recompute.")
    a = ap.parse_args()
    out = Path(a.outdir)

    # ── REPLOT: redraw the box from the raw checkpoint ──────────────────────
    if a.replot:
        raw = out / _RAW_NAME
        if not raw.exists():
            sys.exit(f"--replot: raw checkpoint not found: {raw}")
        print(f"[replot] reading raw checkpoint {raw} (no recompute)", flush=True)
        labels, min_support, tss_by_tool, tts_by_tool = _read_raw_checkpoint(raw)
        out.mkdir(parents=True, exist_ok=True)
        _render(labels, min_support, tss_by_tool, tts_by_tool, out / _PNG_NAME)
        print("END_DIST_BOXPLOT_DONE")
        return

    # ── COMPUTE path ────────────────────────────────────────────────────────
    if not (a.beds and a.cage and a.drna):
        sys.exit("compute mode requires --beds --cage --drna "
                 "(or pass --replot to redraw from the raw checkpoint)")
    canonical_only = not a.all_contigs
    regions = parse_region_values(a.region)
    counts_paths = _parse_counts_args(a.counts)
    out.mkdir(parents=True, exist_ok=True)

    pairs = _parse_beds_labels(a.beds, a.labels)

    print("[load] orthogonal truth peaks (CAGE 5' / dRNA 3') ...", flush=True)
    cage_idx = _build_peak_index(parse_peaks_bed(a.cage, regions, canonical_only))
    drna_idx = _build_peak_index(parse_peaks_bed(a.drna, regions, canonical_only))

    tss_by_tool, tts_by_tool = {}, {}
    labels = []
    omitted = []   # (label, reason) — tools dropped from Fig D entirely
    for lab, bed in pairs:
        counts = load_support_table(counts_paths.get(lab))
        td, ud, has_support = collect_distances(bed, counts, cage_idx, drna_idx,
                                                regions, a.min_support, canonical_only)
        if (not has_support) and a.min_support > 1:
            # NO SILENT DROP: a GTF tool with no --counts cannot honour a real
            # >min-support threshold (min_support>1), so it is OMITTED from Fig D.
            reason = ("GTF input with no --counts; cannot apply "
                      f">{a.min_support} read-support threshold")
            omitted.append((lab, reason))
            print(f"  OMIT {lab} from Fig D: {reason}", file=sys.stderr, flush=True)
            continue
        labels.append(lab)
        tss_by_tool[lab] = td
        tts_by_tool[lab] = ud
        sup = "yes" if has_support else "NO"
        print(f"  {lab}: TSS n={len(td)} (median {np.median(td) if td else 'NA'}), "
              f"TTS n={len(ud)} (median {np.median(ud) if ud else 'NA'}), "
              f"support={sup}", flush=True)

    if omitted:
        print("[omissions] Fig D (end-distance >min-support) excludes: "
              + ", ".join(f"{l} ({r})" for l, r in omitted),
              file=sys.stderr, flush=True)

    # RAW checkpoint (the source of truth for the box + --replot)
    raw = out / _RAW_NAME
    _write_raw_checkpoint(raw, labels, a.min_support, tss_by_tool, tts_by_tool)
    print("wrote", raw)
    # convenience summary (NOT the checkpoint)
    _write_summary(out / _SUMMARY_NAME, labels, tss_by_tool, tts_by_tool)
    print("wrote", out / _SUMMARY_NAME)

    _render(labels, a.min_support, tss_by_tool, tts_by_tool, out / _PNG_NAME)
    print("END_DIST_BOXPLOT_DONE")


if __name__ == "__main__":
    main()
