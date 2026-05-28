#!/usr/bin/env python3
"""
peak_matching.py — Shared, contract-compliant truth-peak matching.

This module is the *single source of truth* for what counts as a true
positive when a called transcript end is compared against an orthogonal
truth peak (CAGE for 5'/TSS, direct-RNA / 3'-seq for 3'/TTS).  It is
imported by both the evaluator (bin/evaluation/*) and the logistic-
regression weight fit (analysis/ted_logreg/ted_logreg.py) so that "TP"
means the same thing everywhere.

It implements the EVAL CONTRACT (2026-05-28 direction memo §5, item 2,
"Truth preprocessing — the R3 fix") and the lit-review dedup-to-summit
guidance:

  (a) Collapse each CAGE/dRNA peak *cluster* (one BED interval) to ONE
      summit.  When a signal track (bedGraph) is supplied the summit is
      the dominant-signal position inside the interval; otherwise it is
      the interval centre.  Strand-aware throughout.

  (b) Match each called end to the NEAREST summit within ±W, with each
      summit matched AT MOST ONCE.  Assignment is greedy nearest-first
      across all (call, summit) pairs in a (chrom, strand) bucket, so a
      fat cluster of called ends near one summit yields exactly ONE true
      positive — not one per call.  This removes the many-to-one
      interval-membership recall inflator (R3).

  (c) Window sweep: ±10 / ±25 / ±50 / ±100, with ±50 as the primary
      window (LRGASP/SQANTI3) and ±25 reported as a secondary 3' column.

The OLD behaviour (interval-membership, many-to-one) lived in:
  - analysis/ted_logreg/ted_logreg.py::_hit_peak
  - analysis/ted_pr_sweep/ted_pr_sweep.py::_hit_peak / _matched_peak_id
  - bin/evaluation/ted_end_precision.py::_nearest_annot   (interval id)
  - bin/evaluation/ted_global_pileup_pr.py::_nearest_peak (interval id)
Each of those credited every call that fell inside (or within ±W of) a
peak interval as a TP, and de-duplicated recall only by interval identity
(so precision/TP counts were many-to-one).  This module replaces that
logic with summit dedup + one-to-one accounting.

Pure-stdlib (no numpy/pandas required) so it can be imported anywhere.
"""
from __future__ import annotations

from bisect import bisect_left
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# Contract windows.  ±50 is the primary (LRGASP/SQANTI3); ±25 is the
# secondary tighter 3' column; ±10/±100 bound the sweep.
SWEEP_WINDOWS: Tuple[int, ...] = (10, 25, 50, 100)
PRIMARY_WINDOW: int = 50
SECONDARY_3PRIME_WINDOW: int = 25

# A (chrom, strand) -> sorted list of summit positions.
SummitIndex = Dict[Tuple[str, str], List[int]]


# ── Summit collapse ──────────────────────────────────────────────────────────

def _interval_summit(start: int, end: int,
                     signal: Optional[Sequence[Tuple[int, float]]] = None) -> int:
    """Collapse one peak interval [start, end) to a single summit position.

    Summit definition (eval contract / litreview dedup-to-summit):
      * If a signal track is supplied (list of (position, value) covering
        the interval), the summit is the position of maximum signal inside
        [start, end).  Ties break to the position closest to the interval
        centre (deterministic, strand-agnostic — initiation/cleavage
        wobble is symmetric about the dominant site).
      * Otherwise the summit is the interval centre (start + end) // 2.

    NOTE (flag for ted-pi / hdheath): no signal bedGraphs are wired into
    the eval harness fixtures yet, so production currently falls back to
    the peak centre.  When the per-sample signal bedGraphs from
    05_06_26.csv are available, pass them in via `signal` so the summit is
    the dominant CTSS / density-peak rather than the geometric centre.
    """
    center = (start + end) // 2
    if not signal:
        return center
    best_pos = center
    best_val = None
    best_center_dist = None
    # `signal` is sorted by position when produced by load_signal_bedgraphs,
    # so bisect to the window start and stop at the window end.  Detect
    # sortedness cheaply (one pass over positions) so an ad-hoc unsorted
    # list (e.g. in unit tests) still scans correctly.
    positions = [p for p, _ in signal]
    sorted_track = all(positions[i] <= positions[i + 1]
                       for i in range(len(positions) - 1))
    if sorted_track:
        lo = bisect_left(positions, start)
        scan = signal[lo:]
    else:
        scan = signal
    for pos, val in scan:
        if pos < start:
            continue
        if pos >= end:
            if sorted_track:
                break
            continue
        cdist = abs(pos - center)
        if (best_val is None or val > best_val or
                (val == best_val and cdist < best_center_dist)):
            best_val = val
            best_pos = pos
            best_center_dist = cdist
    return best_pos if best_val is not None else center


def build_summit_index(
    peaks: Iterable[Tuple[str, int, int, str]],
    signal: Optional[Dict[Tuple[str, str], Sequence[Tuple[int, float]]]] = None,
) -> SummitIndex:
    """Build a (chrom, strand) -> sorted [summit_pos, ...] index.

    `peaks` is an iterable of (chrom, start, end, strand) tuples (BED6
    half-open).  Each interval collapses to exactly one summit, so the
    number of summits equals the number of input peak intervals — this is
    the recall denominator and is what one-to-one matching consumes.

    `signal`, if given, maps (chrom, strand) -> sequence of (pos, value)
    so the summit can be the dominant-signal position; otherwise the
    interval centre is used.
    """
    index: SummitIndex = defaultdict(list)
    for chrom, start, end, strand in peaks:
        sig = None
        if signal is not None:
            sig = signal.get((chrom, strand))
        summit = _interval_summit(int(start), int(end), sig)
        index[(chrom, strand)].append(summit)
    for key in index:
        index[key].sort()
    return index


def load_summit_index_from_bed(
    path: str,
    chroms: Optional[Iterable[str]] = None,
    signal: Optional[Dict[Tuple[str, str], Sequence[Tuple[int, float]]]] = None,
) -> SummitIndex:
    """Read a BED6 peak file and return a summit index (one summit/peak).

    `chroms` optionally restricts to a set/iterable of chromosome names.
    Lines beginning with '#' and rows with < 6 columns are skipped.
    """
    chroms_set = set(chroms) if chroms is not None else None
    rows: List[Tuple[str, int, int, str]] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            chrom = cols[0]
            if chroms_set is not None and chrom not in chroms_set:
                continue
            rows.append((chrom, int(cols[1]), int(cols[2]), cols[5]))
    return build_summit_index(rows, signal=signal)


def count_summits(index: SummitIndex) -> int:
    """Total number of truth summits (the recall denominator)."""
    return sum(len(v) for v in index.values())


# ── Signal bedGraph loading ──────────────────────────────────────────────────

# A (chrom, strand) -> list of (position, value) signal samples.  Each
# bedGraph interval [start, end) contributes one sample at its centre with
# the interval's value; _interval_summit picks the dominant-signal position
# inside a peak from these samples.
SignalTrack = Dict[Tuple[str, str], List[Tuple[int, float]]]


def load_signal_bedgraphs(
    plus_path: Optional[str] = None,
    minus_path: Optional[str] = None,
    chroms: Optional[Iterable[str]] = None,
) -> Optional[SignalTrack]:
    """Load per-strand bedGraph(s) into a (chrom, strand) -> [(pos, val)] track.

    Each bedGraph line is "chrom start end value" (4 columns, half-open
    interval).  Every interval contributes a sample at its centre carrying
    its value; build_summit_index/`_interval_summit` then resolve each peak
    to the dominant-signal position covered by these samples.

    `plus_path` provides the '+' strand track, `minus_path` the '-' strand.
    Either may be None (that strand falls back to peak centre).  Returns None
    if neither path is usable, so callers can pass the result straight to
    `signal=` and transparently fall back to the geometric centre.

    `chroms` optionally restricts to a set/iterable of chromosome names
    (cheap pre-filter for whole-genome bedGraphs on tiny/chr22 fixtures).
    """
    chroms_set = set(chroms) if chroms is not None else None
    track: SignalTrack = defaultdict(list)

    def _ingest(path: Optional[str], strand: str) -> None:
        if not path:
            return
        try:
            fh = open(path)
        except OSError:
            return
        with fh:
            for line in fh:
                if not line or line[0] in "#tb":  # skip comments/track/browser
                    if not line[0].isdigit() and not line.startswith("chr"):
                        continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 4:
                    continue
                chrom = cols[0]
                if chroms_set is not None and chrom not in chroms_set:
                    continue
                try:
                    start = int(cols[1])
                    end = int(cols[2])
                    val = float(cols[3])
                except ValueError:
                    continue
                pos = (start + end) // 2
                track[(chrom, strand)].append((pos, val))

    _ingest(plus_path, "+")
    _ingest(minus_path, "-")

    if not track:
        return None
    # Sort each track by position so _interval_summit's scan is well-ordered
    # (it scans all samples; sorting is for determinism / readability).
    for key in track:
        track[key].sort()
    return dict(track)


# ── One-to-one greedy nearest matching ───────────────────────────────────────

def match_positions_one_to_one(
    positions: Sequence[int],
    summits: Sequence[int],
    window: int,
) -> List[int]:
    """Greedy nearest one-to-one assignment of called ends to summits.

    Parameters
    ----------
    positions : called end positions for ONE (chrom, strand) bucket.
    summits   : sorted summit positions for the SAME bucket.
    window    : max |distance| (bp) for a valid match.

    Returns
    -------
    assignment : list, len == len(positions).  assignment[i] is the index
        into `summits` that call i is matched to, or -1 if unmatched.
        Each summit index appears AT MOST ONCE across the returned list.

    Algorithm
    ---------
    Enumerate every (call_i, summit_j) pair with |pos_i - summit_j| <= W,
    sort by (distance, call_i, summit_j) for determinism, then walk the
    sorted pairs assigning greedily: a pair is taken only if neither its
    call nor its summit is already used.  This is the standard greedy
    nearest-match used by CAGEr/CAGEfightR-style peak benchmarking and
    guarantees one summit is never credited to more than one call.

    Greedy-nearest is not globally optimal (that would need Hungarian),
    but it is the field convention, is order-independent, and — crucially
    for the R3 fix — it makes a fat cluster of N calls near a single
    summit yield exactly ONE true positive.
    """
    n = len(positions)
    assignment = [-1] * n
    if n == 0 or len(summits) == 0 or window < 0:
        return assignment

    sorted_summits = summits  # caller guarantees sorted; do not re-sort here
    pairs: List[Tuple[int, int, int]] = []  # (distance, call_idx, summit_idx)
    for i, pos in enumerate(positions):
        lo = bisect_left(sorted_summits, pos - window)
        j = lo
        while j < len(sorted_summits):
            s = sorted_summits[j]
            if s > pos + window:
                break
            pairs.append((abs(pos - s), i, j))
            j += 1

    if not pairs:
        return assignment

    pairs.sort()
    call_used = [False] * n
    summit_used = [False] * len(sorted_summits)
    for dist, i, j in pairs:
        if call_used[i] or summit_used[j]:
            continue
        assignment[i] = j
        call_used[i] = True
        summit_used[j] = True
    return assignment


def match_ends(
    ends: Sequence[Tuple[str, str, int]],
    summit_index: SummitIndex,
    window: int,
) -> List[int]:
    """One-to-one match a flat list of (chrom, strand, pos) called ends.

    Returns a list, len == len(ends), where entry k is a stable summit id
    string "chrom:strand:summit_pos" if call k matched a summit, else None.
    Each summit id appears AT MOST ONCE in the returned list.

    This is the convenience entry point for callers that have ends from
    many (chrom, strand) buckets mixed together (e.g. ted_logreg, the
    pileup PR evaluator).  It buckets internally, runs the greedy one-to-
    one matcher per bucket, and flattens the result back to the input
    order.
    """
    # Bucket call indices by (chrom, strand), preserving original order.
    buckets: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    for k, (chrom, strand, _pos) in enumerate(ends):
        buckets[(chrom, strand)].append(k)

    matched: List[Optional[str]] = [None] * len(ends)
    for key, idxs in buckets.items():
        summits = summit_index.get(key, [])
        if not summits:
            continue
        positions = [ends[k][2] for k in idxs]
        assign = match_positions_one_to_one(positions, summits, window)
        for local_i, k in enumerate(idxs):
            j = assign[local_i]
            if j >= 0:
                matched[k] = f"{key[0]}:{key[1]}:{summits[j]}"
    return matched


# ── Precision / recall on a matched call set ─────────────────────────────────

def precision_recall(
    ends: Sequence[Tuple[str, str, int]],
    summit_index: SummitIndex,
    window: int,
) -> Dict[str, float]:
    """Compute one-to-one precision / recall / F1 for a call set.

    precision = (# calls matched to a summit) / (# calls)
    recall    = (# distinct summits matched)  / (total summits)

    Because matching is one-to-one, #calls-matched == #summits-matched,
    so precision and recall share the same numerator (TP) — exactly the
    contract's "each summit matched at most once" accounting.  A cluster
    of N calls on one summit contributes TP = 1 (not N), removing the
    recall inflation.
    """
    matched = match_ends(ends, summit_index, window)
    tp = sum(1 for m in matched if m is not None)
    n_calls = len(ends)
    n_summits = count_summits(summit_index)
    precision = tp / n_calls if n_calls else float("nan")
    recall = tp / n_summits if n_summits else float("nan")
    f1 = (2 * precision * recall / (precision + recall)
          if precision == precision and recall == recall
          and (precision + recall) > 0 else float("nan"))
    return {
        "tp": tp,
        "n_calls": n_calls,
        "n_summits": n_summits,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "window": window,
    }


def window_sweep_precision_recall(
    ends: Sequence[Tuple[str, str, int]],
    summit_index: SummitIndex,
    windows: Iterable[int] = SWEEP_WINDOWS,
) -> Dict[int, Dict[str, float]]:
    """Run precision_recall at each window in the sweep (default ±10/25/50/100)."""
    return {w: precision_recall(ends, summit_index, w) for w in windows}
