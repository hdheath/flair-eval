#!/usr/bin/env python3
"""Near-miss localization slice for the TSS/TTS sequence classifier.

FIX 3 / ML-review T1 (Q2), 2026-05-28. Eval-side diagnostic — does NOT touch
ted.py or the model. Measures the classifier's score FALLOFF as a function of
distance from a true end: for each true CAGE (TSS) / dRNA (TTS) summit, score
``model.predict_proba`` at the true summit and at decoy positions shifted by
+/- {5,10,25,50,100} bp ("near-misses"). A model that genuinely LOCALIZES the
end produces a sharp peak at offset 0 and a steep falloff; a model that merely
recognizes the REGION (the §3 negative-set ceiling the ML review names) gives a
FLAT curve where decoys score about as high as the true site.

Design:
  * Consumes a model .pkl + a peak BED + a genome FASTA — needs NO chr22 run to
    exist, but is built to consume the chr22 CAGE/dRNA peak BEDs + the calibrated
    .pkl + the hg38 genome directly afterwards.
  * Reuses ``flair.ted.SequenceFeatureExtractor`` so the feature/window/strand
    convention is IDENTICAL to how TED scores a candidate end (train<->TED
    parity, including the BED-edge closing-edge -1, see prepare_data.py / the
    centroid fix in ted.py). The summit is taken from the shared
    ``peak_matching`` summit definition (peak centre, or dominant signal when a
    bedGraph is supplied) so "true position" means the same thing as everywhere
    else in the eval harness.

Outputs:
  * <out_prefix>.near_miss.tsv      : one row per (peak, offset) with the score.
  * <out_prefix>.near_miss_summary.tsv : per-offset aggregate (mean/median score,
    mean delta-from-true, frac decoys >= true) + a single localization_auc-style
    falloff number.

This module is import-safe: the heavy deps (the extractor, the model) load only
inside ``main`` / ``score_near_miss`` so the unit test can exercise the pure
falloff math on synthetic scores without a genome or a pickle.
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# Default near-miss offsets (bp). 0 = the true summit. Symmetric +/- shifts.
DEFAULT_OFFSETS: Tuple[int, ...] = (5, 10, 25, 50, 100)


def signed_offsets(offsets: Sequence[int]) -> List[int]:
    """[5,10] -> [-10,-5,0,5,10] (sorted, deduped, always including 0)."""
    s = {0}
    for o in offsets:
        s.add(int(o))
        s.add(-int(o))
    return sorted(s)


# ── Pure falloff math (unit-testable without a genome / model) ──────────────

def summarize_falloff(
    scores_by_offset: Dict[int, List[float]],
) -> Dict[int, Dict[str, float]]:
    """Aggregate per-offset score lists into falloff statistics.

    For each offset returns: n, mean, median, and (vs the true offset 0)
    mean_delta_from_true (score@offset - score@0, averaged per matched peak)
    and frac_ge_true (fraction of decoys scoring >= the true site). A model
    that LOCALIZES has mean_delta_from_true << 0 and frac_ge_true ~ 0 at large
    offsets; a REGION detector has both ~flat (delta ~ 0, frac_ge_true ~ 0.5+).

    `scores_by_offset` must be peak-aligned: index i is the same peak across
    every offset (so delta-from-true is a paired comparison).
    """
    if 0 not in scores_by_offset:
        raise ValueError("scores_by_offset must contain the true offset 0")
    true_scores = scores_by_offset[0]
    out: Dict[int, Dict[str, float]] = {}
    for offset, vals in scores_by_offset.items():
        finite = [v for v in vals if v == v]  # drop NaN
        n = len(finite)
        mean = statistics.fmean(finite) if finite else float("nan")
        median = statistics.median(finite) if finite else float("nan")
        # Paired delta vs the true site (same peak index).
        deltas, ge = [], 0
        paired = 0
        for i, v in enumerate(vals):
            if i < len(true_scores):
                t = true_scores[i]
                if v == v and t == t:
                    deltas.append(v - t)
                    paired += 1
                    if v >= t:
                        ge += 1
        out[offset] = {
            "n": float(n),
            "mean": mean,
            "median": median,
            "mean_delta_from_true": (statistics.fmean(deltas) if deltas
                                     else float("nan")),
            "frac_ge_true": (ge / paired if paired else float("nan")),
        }
    return out


def localization_score(summary: Dict[int, Dict[str, float]]) -> float:
    """One scalar in [0,1]: how well the model localizes the end.

    Defined as the mean over the largest +/- offset of (1 - frac_ge_true),
    i.e. how rarely a far decoy beats the true site. 1.0 = perfect
    localization (no decoy ever beats true); ~0.5 = flat / region detector.
    """
    far = max((o for o in summary if o > 0), default=0)
    if far == 0:
        return float("nan")
    vals = []
    for o in (far, -far):
        f = summary.get(o, {}).get("frac_ge_true")
        if f is not None and f == f:
            vals.append(1.0 - f)
    return statistics.fmean(vals) if vals else float("nan")


# ── Isoform-end-driven near-miss falloff (no model needed) ──────────────────
# The complementary slice the ML/PI review's A5 asks for: instead of probing the
# MODEL at synthetic decoys, take the ACTUAL isoform ends a TED run emitted, bin
# them by SIGNED distance to the nearest truth summit (one-to-one matched), and
# report mean TED reality + one-to-one precision per distance bin. This is the
# end-side near-miss falloff: a good end-caller concentrates its emitted ends at
# distance 0 with high reality and precision; a region-y caller smears ends out
# to +/- tens of bp with reality/precision barely decaying with distance.

# Default signed-distance bin EDGES (bp). Symmetric; the open tails (< -100,
# >= 100) capture far ends. A point exactly on a summit lands in [0,10).
DEFAULT_DIST_BIN_EDGES: Tuple[int, ...] = (
    -100, -50, -25, -10, -5, 0, 5, 10, 25, 50, 100)


def _bin_label(edges: Sequence[int], d: int) -> str:
    """Human label for the bin a signed distance d falls into."""
    edges = list(edges)
    if d < edges[0]:
        return f"<{edges[0]}"
    for lo, hi in zip(edges[:-1], edges[1:]):
        if lo <= d < hi:
            return f"[{lo},{hi})"
    return f">={edges[-1]}"


def signed_distance_to_nearest_summit(
    positions: Sequence[int],
    summits: Sequence[int],
    strand: str,
) -> List[Optional[int]]:
    """Signed distance from each called end to its NEAREST summit (same bucket).

    Sign is in transcription direction: a positive distance means the called end
    is DOWNSTREAM of the truth summit (call - summit on +; summit - call on -),
    matching the decoy-offset convention in `score_near_miss`. `summits` must be
    sorted. Returns None for a bucket with no summits.
    """
    if not summits:
        return [None] * len(positions)
    out: List[Optional[int]] = []
    for pos in positions:
        lo = bisect_left(summits, pos)
        cands = []
        if lo > 0:
            cands.append(summits[lo - 1])
        if lo < len(summits):
            cands.append(summits[lo])
        nearest = min(cands, key=lambda s: abs(pos - s))
        signed = (pos - nearest) if strand == "+" else (nearest - pos)
        out.append(int(signed))
    return out


def bin_ends_by_signed_distance(
    ends: Sequence[Tuple[str, str, int]],
    summit_index,
    window: int,
    *,
    reality: Optional[Sequence[float]] = None,
    bin_edges: Sequence[int] = DEFAULT_DIST_BIN_EDGES,
) -> Dict[str, Dict[str, float]]:
    """Bin called isoform ends by SIGNED distance to nearest truth summit.

    Parameters
    ----------
    ends : flat list of (chrom, strand, pos) called ends (mixed buckets ok).
    summit_index : {(chrom, strand): sorted[int]} truth summits (peak_matching).
    window : the one-to-one match window (bp) used for the per-bin PRECISION
             (an end is a TP iff it greedily one-to-one matches a summit within
             this window — the SAME matcher the rest of the harness uses).
    reality : optional per-end TED reality/confidence (peak-aligned to `ends`).
    bin_edges : signed-distance bin edges (bp).

    Returns {bin_label: {n, frac_ends, mean_reality, precision, mean_abs_dist}}.
    `precision` per bin = (# ends in that bin that are one-to-one TPs) / (# ends
    in that bin) — the near-miss falloff: precision should drop as |distance|
    grows. `mean_reality` is NaN when no reality vector is given.
    """
    # One-to-one TP mask (shared matcher) so per-bin precision is honest.
    # Lazy import keeps this module import-safe for the pure unit test; works
    # both as a CLI (bin/evaluation on path) and as a package (evaluation.*).
    try:
        from peak_matching import match_ends  # CLI / direct
    except ImportError:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        try:
            from peak_matching import match_ends
        except ImportError:
            from evaluation.peak_matching import match_ends  # package
    matched = match_ends(ends, summit_index, window)  # list of summit-id or None
    # Nearest-summit signed distance per end (independent of the match window).
    buckets: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    for k, (chrom, strand, _pos) in enumerate(ends):
        buckets[(chrom, strand)].append(k)
    signed_dist: List[Optional[int]] = [None] * len(ends)
    for (chrom, strand), idxs in buckets.items():
        summits = summit_index.get((chrom, strand), [])
        positions = [ends[k][2] for k in idxs]
        dists = signed_distance_to_nearest_summit(positions, summits, strand)
        for local_i, k in enumerate(idxs):
            signed_dist[k] = dists[local_i]

    agg: Dict[str, dict] = {}
    total = 0
    for k, d in enumerate(signed_dist):
        if d is None:
            continue  # no summit in this bucket; not a near-miss observation
        total += 1
        label = _bin_label(bin_edges, d)
        a = agg.setdefault(label, {"n": 0, "tp": 0, "abs": 0,
                                   "real": [], "edge_lo": None})
        a["n"] += 1
        a["abs"] += abs(d)
        if matched[k] is not None:
            a["tp"] += 1
        if reality is not None and k < len(reality):
            rv = reality[k]
            if rv == rv:  # not NaN
                a["real"].append(float(rv))

    out: Dict[str, Dict[str, float]] = {}
    for label, a in agg.items():
        n = a["n"]
        out[label] = {
            "n": float(n),
            "frac_ends": (n / total) if total else float("nan"),
            "precision": (a["tp"] / n) if n else float("nan"),
            "mean_abs_dist": (a["abs"] / n) if n else float("nan"),
            "mean_reality": (statistics.fmean(a["real"]) if a["real"]
                             else float("nan")),
        }
    return out


# ── Scoring against a real model + genome (needs the heavy deps) ────────────

def score_near_miss(
    model,
    extractor,
    summits: Sequence[Tuple[str, int, str]],
    offsets: Sequence[int],
    boundary_type: str,
) -> Tuple[List[dict], Dict[int, List[float]]]:
    """Score `model` at each summit and its near-miss decoys.

    `summits` is a list of (chrom, summit_pos, strand). `extractor` is a
    flair.ted.SequenceFeatureExtractor for the matching boundary_type.

    Returns (rows, scores_by_offset). `scores_by_offset[o][i]` is the score for
    peak i at signed offset o (NaN if the window fell off the contig), so the
    lists stay peak-aligned for the paired delta-from-true stats.

    The decoy genomic position is summit + offset on the + strand and
    summit - offset on the - strand, so a positive offset always means
    "downstream of the true end" in transcription direction — matching the
    extractor's strand-aware window and the shifted-negative semantics the
    realistic-negative retrain (B1/F1) will use.
    """
    import numpy as np

    signed = signed_offsets(offsets)
    scores_by_offset: Dict[int, List[float]] = {o: [] for o in signed}
    rows: List[dict] = []
    for chrom, summit, strand in summits:
        for o in signed:
            decoy = summit + o if strand == "+" else summit - o
            try:
                X = extractor.extract_batch([(chrom, decoy, strand)])
                if np.any(np.isnan(X)):
                    score = float("nan")
                else:
                    score = float(model.predict_proba(X)[0, 1])
            except Exception:
                score = float("nan")
            scores_by_offset[o].append(score)
            rows.append({
                "boundary_type": boundary_type,
                "chrom": chrom,
                "summit": summit,
                "strand": strand,
                "offset": o,
                "decoy_pos": decoy,
                "score": score,
            })
    return rows, scores_by_offset


def _summits_from_index(index) -> List[Tuple[str, int, str]]:
    out: List[Tuple[str, int, str]] = []
    for (chrom, strand), positions in index.items():
        for pos in positions:
            out.append((chrom, int(pos), strand))
    return out


def _isoform_ends_to_calls(iso_bed: Path, boundary_type: str,
                           chroms: Optional[Sequence[str]] = None):
    """Read isoform BED6/BED12 -> [(chrom, strand, end_pos)] for the chosen end.

    TSS = 5' end (start on +, end on -); TTS = 3' end (end on +, start on -).
    Returns (ends, names) so an optional ted_log reality column can be aligned
    by isoform name.
    """
    ends: List[Tuple[str, str, int]] = []
    names: List[str] = []
    chrom_set = set(chroms) if chroms else None
    with open(iso_bed) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith(("#", "track", "browser")):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            chrom, start, end, name, _score, strand = f[0], int(f[1]), int(f[2]), f[3], f[4], f[5]
            if chrom_set is not None and chrom not in chrom_set:
                continue
            if boundary_type == "tss":
                pos = start if strand == "+" else end
            else:  # tts
                pos = end if strand == "+" else start
            ends.append((chrom, strand, int(pos)))
            names.append(name)
    return ends, names


def _reality_by_name_from_ted_log(ted_log: Path, boundary_type: str) -> Dict[str, float]:
    """Map isoform name -> TED reality for the chosen end from a ted_log.tsv.

    Uses the global_assign 'pass' rows (one per emitted cluster). Reads the
    TED_tss_reality / TED_tts_reality column; falls back to TED_confidence.
    Best-effort and header-driven so it survives column reordering.
    """
    out: Dict[str, float] = {}
    real_col = f"TED_{boundary_type}_reality"
    with open(ted_log) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        name_i = idx.get("read_name", idx.get("cluster_id", idx.get("junc_id")))
        col = idx.get(real_col, idx.get("TED_confidence"))
        if name_i is None or col is None:
            return out
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) <= max(name_i, col):
                continue
            try:
                out[f[name_i]] = float(f[col])
            except (ValueError, IndexError):
                continue
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    # Two modes: 'model' (default, score model at synthetic decoys) and
    # 'ends' (bin ACTUAL isoform ends by signed distance to nearest summit and
    # report mean reality + one-to-one precision per bin).
    top = argparse.ArgumentParser(
        description="Near-miss slices for TED end calling.")
    top.add_argument("--mode", choices=["model", "ends"], default="model",
                     help="model: classifier score falloff at synthetic decoys "
                          "(needs a .pkl + genome). ends: bin real isoform ends "
                          "by signed distance to nearest truth summit, reporting "
                          "mean TED reality + one-to-one precision per bin.")
    known, rest = top.parse_known_args(argv)
    if known.mode == "ends":
        return _main_ends(rest)
    return _main_model(rest)


def _main_ends(argv: Sequence[str]) -> int:
    p = argparse.ArgumentParser(
        description="Near-miss falloff over ACTUAL isoform ends: bin by signed "
                    "distance to nearest truth summit; report mean reality + "
                    "one-to-one precision per bin.")
    p.add_argument("--isoforms", required=True, type=Path,
                   help="Isoform BED6/BED12 (the TED run's emitted ends).")
    p.add_argument("--peaks", required=True, type=Path,
                   help="Truth peak BED6 (CAGE for TSS, dRNA/QuantSeq for TTS).")
    p.add_argument("--boundary-type", required=True, choices=["tss", "tts"])
    p.add_argument("--out-prefix", required=True, type=Path)
    p.add_argument("--window", type=int, default=50,
                   help="One-to-one match window (bp) for per-bin precision. "
                        "Default 50 (matches the harness end-PR convention).")
    p.add_argument("--ted-log", type=Path, default=None,
                   help="Optional ted_log.tsv to attach TED reality per end "
                        "(by isoform name) so each bin reports mean reality.")
    p.add_argument("--chroms", nargs="+", default=None)
    p.add_argument("--signal-plus", type=Path, default=None)
    p.add_argument("--signal-minus", type=Path, default=None)
    args = p.parse_args(argv)

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from peak_matching import (load_summit_index_from_bed,
                               load_signal_bedgraphs, count_summits)

    signal = None
    if args.signal_plus or args.signal_minus:
        signal = load_signal_bedgraphs(
            str(args.signal_plus) if args.signal_plus else None,
            str(args.signal_minus) if args.signal_minus else None,
            chroms=args.chroms)
    summit_index = load_summit_index_from_bed(
        str(args.peaks), chroms=args.chroms, signal=signal)

    ends, names = _isoform_ends_to_calls(
        args.isoforms, args.boundary_type, chroms=args.chroms)
    if not ends:
        p.error(f"No isoform ends loaded from {args.isoforms}.")

    reality = None
    if args.ted_log:
        rmap = _reality_by_name_from_ted_log(args.ted_log, args.boundary_type)
        reality = [rmap.get(nm, float("nan")) for nm in names]

    bins = bin_ends_by_signed_distance(
        ends, summit_index, args.window, reality=reality)

    out_tsv = Path(str(args.out_prefix) + ".near_miss_ends.tsv")
    # Order rows by bin lower edge for a readable falloff curve.
    def _sort_key(lbl: str) -> int:
        s = lbl.strip("<>=[)")
        try:
            return int(s.split(",")[0])
        except ValueError:
            return 10 ** 9
    with open(out_tsv, "w") as fh:
        fh.write("bin\tn\tfrac_ends\tprecision\tmean_reality\tmean_abs_dist\n")
        for lbl in sorted(bins, key=_sort_key):
            b = bins[lbl]
            fh.write(f"{lbl}\t{int(b['n'])}\t{b['frac_ends']:.6f}\t"
                     f"{b['precision']:.6f}\t{b['mean_reality']:.6f}\t"
                     f"{b['mean_abs_dist']:.6f}\n")

    out_json = Path(str(args.out_prefix) + ".near_miss_ends.json")
    with open(out_json, "w") as fh:
        json.dump({
            "boundary_type": args.boundary_type,
            "n_ends": len(ends),
            "n_summits": count_summits(summit_index),
            "window": args.window,
            "per_bin": bins,
        }, fh, indent=2)

    print(f"[near_miss_slice:ends] {args.boundary_type}: {len(ends)} ends vs "
          f"{count_summits(summit_index)} summits -> {out_tsv}")
    return 0


def _main_model(argv: Sequence[str]) -> int:
    p = argparse.ArgumentParser(
        description="Near-miss localization slice for the TSS/TTS classifier "
                    "(score falloff vs distance from the true end).")
    p.add_argument("--model", required=True, type=Path,
                   help="Classifier .pkl (calibrated or raw); same predict_proba "
                        "interface TED uses.")
    p.add_argument("--genome", required=True, type=Path,
                   help="Genome FASTA (indexed; matches the TED run's genome).")
    p.add_argument("--peaks", required=True, type=Path,
                   help="True-end peak BED6 (CAGE for TSS, dRNA/QuantSeq for TTS).")
    p.add_argument("--boundary-type", required=True, choices=["tss", "tts"],
                   help="Which end / which model window convention to use.")
    p.add_argument("--out-prefix", required=True, type=Path)
    p.add_argument("--offsets", type=int, nargs="+", default=list(DEFAULT_OFFSETS),
                   help="Unsigned near-miss offsets in bp (symmetric +/- applied). "
                        "Default 5 10 25 50 100.")
    p.add_argument("--chroms", nargs="+", default=None,
                   help="Restrict to these chromosomes (e.g. chr22).")
    p.add_argument("--signal-plus", type=Path, default=None,
                   help="Optional signal bedGraph (+) so the summit is the "
                        "dominant-signal position, not the peak centre.")
    p.add_argument("--signal-minus", type=Path, default=None,
                   help="Optional signal bedGraph (-).")
    args = p.parse_args(argv)

    # Heavy imports kept local so the module imports without these deps.
    import pickle
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from peak_matching import (load_summit_index_from_bed,
                               load_signal_bedgraphs, count_summits)
    try:
        from flair.ted import SequenceFeatureExtractor
    except ImportError as e:
        p.error(f"Could not import flair.ted.SequenceFeatureExtractor "
                f"(run under the flair-ted env): {e}")

    signal = None
    if args.signal_plus or args.signal_minus:
        signal = load_signal_bedgraphs(
            str(args.signal_plus) if args.signal_plus else None,
            str(args.signal_minus) if args.signal_minus else None,
            chroms=args.chroms)

    index = load_summit_index_from_bed(
        str(args.peaks), chroms=args.chroms, signal=signal)
    summits = _summits_from_index(index)
    if not summits:
        p.error(f"No summits loaded from {args.peaks} (chroms={args.chroms}).")

    with open(args.model, "rb") as fh:
        model = pickle.load(fh)
    extractor = SequenceFeatureExtractor(str(args.genome),
                                         boundary_type=args.boundary_type)

    rows, scores_by_offset = score_near_miss(
        model, extractor, summits, args.offsets, args.boundary_type)
    summary = summarize_falloff(scores_by_offset)
    loc = localization_score(summary)

    out_tsv = Path(str(args.out_prefix) + ".near_miss.tsv")
    with open(out_tsv, "w") as fh:
        fh.write("boundary_type\tchrom\tsummit\tstrand\toffset\tdecoy_pos\tscore\n")
        for r in rows:
            fh.write(f"{r['boundary_type']}\t{r['chrom']}\t{r['summit']}\t"
                     f"{r['strand']}\t{r['offset']}\t{r['decoy_pos']}\t"
                     f"{r['score']:.6f}\n")

    out_sum = Path(str(args.out_prefix) + ".near_miss_summary.tsv")
    with open(out_sum, "w") as fh:
        fh.write("offset\tn\tmean_score\tmedian_score\t"
                 "mean_delta_from_true\tfrac_ge_true\n")
        for o in sorted(summary):
            s = summary[o]
            fh.write(f"{o}\t{int(s['n'])}\t{s['mean']:.6f}\t{s['median']:.6f}\t"
                     f"{s['mean_delta_from_true']:.6f}\t{s['frac_ge_true']:.6f}\n")

    out_json = Path(str(args.out_prefix) + ".near_miss_summary.json")
    with open(out_json, "w") as fh:
        json.dump({
            "boundary_type": args.boundary_type,
            "n_summits": count_summits(index),
            "offsets": signed_offsets(args.offsets),
            "localization_score": loc,
            "per_offset": {str(o): summary[o] for o in summary},
        }, fh, indent=2)

    print(f"[near_miss_slice] {args.boundary_type}: {len(summits)} summits, "
          f"localization_score={loc:.4f}")
    print(f"  rows  -> {out_tsv}")
    print(f"  summ  -> {out_sum} / {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
