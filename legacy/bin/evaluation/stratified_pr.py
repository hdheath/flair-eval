#!/usr/bin/env python3
"""
stratified_pr.py — precision/recall stratified by:
  (1) reads assigned per isoform  (cumulative min-reads threshold)
  (2) ground-truth peak score / TPM  (recall, quartile bins)
  (3) distance from annotation     (precision binned by CALL dist; recall binned by PEAK dist)

Cross-mode curves, 5' (vs CAGE) and 3' (vs dRNA). One line per transcriptome mode.

CONVENTION (CORRECTED METRIC — matches the pipeline's single-end headline in
ted_end_precision.py, 5prime_precision/3prime_precision):

  * MANY-TO-ONE membership: a CALL (isoform end) is "matched" if it is within
    --window bp of ANY truth-peak INTERVAL (peak_matching._match_distance <= W).
    Duplication is allowed — a promoter/cleavage region that legitimately yields
    several isoforms credits every nearby end. This is the SAME rule
    ted_end_precision._many_to_one_match_counts uses, so the single-end numbers
    here match combined_evaluation.tsv / the ortho TSV headline and cannot
    diverge.
  * a PEAK is "recovered" if it is within --window bp of ANY call (same OR predicate).
  precision = matched calls / calls ;  recall = recovered peaks / peaks  (within a stratum).
  * CANONICAL chroms only (chr1-22, chrX, chrY) — peak_matching.is_canonical_chrom.

The predicate change vs the legacy version: summit±W -> within-W-of-INTERVAL
(via peak_matching.Peak.lo/.hi + _match_distance), plus the canonical filter.
The MATCHING SCHEME stays MANY-TO-ONE (we did NOT switch to one-to-one).

Inputs (module-friendly; all modes share the same cage/drna/gtf):
  --beds B1 B2 ...        per-mode isoforms BED12 (FLAIR; col5 = read support)
                         or assembler GTF (Bambu/IsoQuant/StringTie2/LRAA/IsoSeq)
  --labels L1 L2 ...      mode labels (same length as --beds)
  --counts L:path ...     OPTIONAL per-tool read-count tables (non-FLAIR assemblers
                         need this; FLAIR/TED have inline support).  Each token is
                         LABEL:path.  Falls back to FlairBed col16 / BED col5.
  --cage  CAGE.bed        5' truth peaks (BED6; col5 = score/TPM)
  --drna  DRNA.bed        3' truth peaks (BED6)
  --gtf   REF.gtf         reference annotation (for distance-from-annotation)
  --outdir DIR
  --region chrN[:s-e]     optional chrom/region restrict
  --window 50  --novel-dist 100

Outputs (PNG + tidy TSV):
  stratified_pr.tsv                            tidy long-form table
  pr_vs_support (precision + recall vs reads)  stratified_precision_by_reads.png /
                                               stratified_recall_by_reads.png
  pr_tradeoff_support (5'/3')                  pr_tradeoff_support.png
  pr_vs_annotdist (5'/3')                      stratified_precision_by_annotdist.png /
                                               stratified_recall_by_annotdist.png
  (plus recall-by-gtscore, retained)

PAIRED-END (TSS+TTS) STRATIFIED METRICS (added 2026-06-22)
==========================================================
In addition to the single-end (5'/3') series above, this script now emits a
PAIRED metric per tool, stratified by the SAME read-support KS sweep.  The
paired metric is the headline "both ends right" view: an isoform end-pair is a
true positive only when its TSS matches a CAGE peak AND its TTS matches a dRNA
peak.

  (1) PAIRED PRECISION(k) — JC-DEDUPLICATED.  Take isoforms with
      read_support >= k, group by JUNCTION CHAIN (single-exon = overlap-grouped
      on the same strand), deduplicate to unique (TSS-peak, TTS-peak) pairs,
      and report
          paired_precision = (unique pairs where BOTH ends hit) / (total isoforms).
      This is computed by CALLING ted_end_precision.compute_jc_deduplicated_
      precision_recall on the support-filtered isoform subset and reading its
      `paired_dedup_precision`.  Because it IS that function (not a re-implement),
      the k=1 value reconciles EXACTLY with the canonical paired headline in
      precision_recall_summary.tsv on the same input (verified on chr22 for
      TED-hybrid-me1 + FLAIR-default).

  (2) PAIRED RECALL(k) — vs a tool-independent VALIDATED TRUTH-PAIR set.  Build
      ONCE from the reference GTF: every transcript -> (TSS, TTS) end-pair, KEPT
      only when its TSS is within-50 of a CAGE peak interval AND its TTS is
      within-50 of a dRNA peak interval (orthogonally-validated reference
      end-pairs), then deduplicated by (TSS, TTS) within a small tolerance W.
      For each tool / threshold k,
          paired_recall(k) = (validated truth pairs RECOVERED) / (total validated pairs),
      where a truth pair is RECOVERED when SOME isoform with read_support >= k has
      its TSS within-50 of the truth TSS AND its TTS within-50 of the truth TTS.
      The denominator is the SAME for every tool (tool-independent truth), so the
      recall lines are directly comparable.  Canonical chroms only.

  (3) OUTPUT.  Paired rows are appended to stratified_pr.tsv with
      stratifier=support, end='paired', metric in {paired_precision,
      paired_recall} (the TSV stays tidy).  Two single-panel PNGs are added:
      paired_precision_by_reads.png + paired_recall_by_reads.png (x = reads /
      isoform, log KS; one line per tool).  An OPTIONAL combined paired
      PR-tradeoff (paired_pr_tradeoff.png) is also drawn.  Tools with no usable
      read support (StringTie2 = NO_COUNTS) are omitted from the paired support
      figures, logged the same way as the single-end ones.

The SINGLE-END code path and its metric (within-50-of-interval, canonical,
many-to-one) are UNCHANGED.
"""
import argparse, csv, sys
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
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
try:
    from peak_matching import (
        Peak, _match_distance, is_canonical_chrom, CANONICAL_CHROMS_SET)
except ImportError:
    from evaluation.peak_matching import (
        Peak, _match_distance, is_canonical_chrom, CANONICAL_CHROMS_SET)
# Paired-end machinery: reuse the CANONICAL JC-dedup metric + the SAME
# truth-peak interval matcher (_nearest_annot) the paired headline uses, so the
# paired_precision computed here is byte-identical to precision_recall_summary's
# paired_dedup_precision on the same input.
try:
    from ted_end_precision import (
        compute_jc_deduplicated_precision_recall,
        parse_peaks_bed as _ted_parse_peaks_bed,
        parse_gtf_transcripts as _ted_parse_gtf_transcripts,
        parse_gtf_ends as _ted_parse_gtf_ends,
        _nearest_annot)
except ImportError:
    from evaluation.ted_end_precision import (
        compute_jc_deduplicated_precision_recall,
        parse_peaks_bed as _ted_parse_peaks_bed,
        parse_gtf_transcripts as _ted_parse_gtf_transcripts,
        parse_gtf_ends as _ted_parse_gtf_ends,
        _nearest_annot)

KS = [1, 2, 3, 5, 8, 12, 20, 35, 50]
DIST_BINS = [(0, 50), (50, 100), (100, 500), (500, 2000), (2000, 10**12)]
DIST_LABELS = ["0-50", "50-100", "100-500", "500-2k", ">2k"]
PALETTE = ["#1b1b6f", "#185fa5", "#1d9e75", "#63992f", "#b19a00", "#d4537e", "#993c1d"]


# ── Region handling ───────────────────────────────────────────────────────────

def parse_region(region):
    """Parse 'chr' or 'chr:start-end' -> (chrom, start|None, end|None)."""
    if region is None:
        return None
    if ":" not in region:
        return (region, None, None)
    chrom, coords = region.split(":", 1)
    s, e = coords.split("-", 1)
    s, e = int(s), int(e)
    if s > e:
        s, e = e, s
    return (chrom, s, e)


def _in_region(chrom, pos, reg):
    """chrom (+ optional pos) within region tuple."""
    if reg is None:
        return True
    rc, rs, re_ = reg
    if chrom != rc:
        return False
    if rs is None or re_ is None or pos is None:
        return True
    return rs <= pos <= re_


# ── GTF annotated ends (distance-from-annotation stratifier) ──────────────────

def parse_gtf_ends(gtf_path, reg=None, canonical_only=True):
    """{(chrom,strand): sorted [tss]}, {(chrom,strand): sorted [tts]} from transcripts."""
    tss = defaultdict(set); tts = defaultdict(set)
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] not in ("transcript", "mRNA"):
                continue
            chrom = c[0]
            if canonical_only and not is_canonical_chrom(chrom):
                continue
            start, end, strand = int(c[3]) - 1, int(c[4]), c[6]
            if reg is not None and reg[0] != chrom:
                continue
            if strand == "+":
                tss[(chrom, "+")].add(start); tts[(chrom, "+")].add(end)
            else:
                tss[(chrom, "-")].add(end); tts[(chrom, "-")].add(start)
    return ({k: sorted(v) for k, v in tss.items()},
            {k: sorted(v) for k, v in tts.items()})


# ── Truth peaks as INTERVAL summit index (Peak with .lo/.hi) ──────────────────

def load_peaks(bed, reg=None, canonical_only=True):
    """{(chrom,strand): (sorted [Peak], scores_aligned)} from BED6.

    Each peak interval [start, end) is carried as a Peak (summit=center) with
    .lo/.hi so the matcher credits a call within W of ANY part of the interval.
    """
    tmp = defaultdict(list)
    with open(bed) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            chrom = c[0]
            if canonical_only and not is_canonical_chrom(chrom):
                continue
            start, end, strand = int(c[1]), int(c[2]), c[5]
            center = (start + end) // 2
            if reg is not None and not _in_region(chrom, center, reg):
                continue
            try:
                score = float(c[4])
            except ValueError:
                score = 0.0
            tmp[(chrom, strand)].append((Peak(center, start, end), score))
    out = {}
    for k, v in tmp.items():
        v.sort(key=lambda t: int(t[0]))
        out[k] = ([p for p, _ in v], [s for _, s in v])
    return out


def _peak_half_extents(peaks):
    """(max summit-to-lo, max hi-to-summit) over a peak bucket — for scan bounds."""
    ml = mr = 0
    for p in peaks:
        s = int(p)
        ml = max(ml, s - p.lo)
        mr = max(mr, p.hi - s)
    return ml, mr


def _near_interval(pos, sorted_peaks, ml, mr, window):
    """MANY-TO-ONE: True iff `pos` is within `window` of ANY peak interval."""
    if not sorted_peaks:
        return False
    lo_bound = pos - (window + mr)
    hi_bound = pos + (window + ml)
    j = bisect_left(sorted_peaks, lo_bound)
    while j < len(sorted_peaks):
        p = sorted_peaks[j]
        if int(p) > hi_bound:
            break
        if _match_distance(pos, p) <= window:
            return True
        j += 1
    return False


def _dist(pos, sorted_positions):
    """Nearest scalar distance (annotation distance; plain int positions)."""
    if not sorted_positions:
        return 10**12
    i = bisect_left(sorted_positions, pos)
    best = 10**12
    for j in (i - 1, i):
        if 0 <= j < len(sorted_positions):
            best = min(best, abs(pos - sorted_positions[j]))
    return best


# ── Calls (BED12 / GTF) + per-tool read support ───────────────────────────────

def _resolve_support(iso, counts):
    """Read support for an isoform: counts file -> FlairBed col16 -> BED col5.

    `counts` may be a {name: count} map (from load_read_map / FLAIR/IsoQuant/
    Bambu counts) or None.  parse_isoforms already exposes the inline FlairBed
    read_support (col16) and plain BED score (col5).
    """
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


def load_counts_map(path):
    """Load a per-tool count source -> {transcript_id: count}.

    Auto-detects two file shapes (the manifest's isoform_counts column points to
    a numeric table; isoform_read_map points to a read-list):
      * READ-MAP   (isoform_id<TAB>read1,read2,...): count = #reads (load_read_map).
        Detected when the second column contains a comma.
      * NUMERIC COUNTS (id<TAB>count[<TAB>count]): count = MAX numeric column to
        the right (covers FLAIR/TED `id count count` duplicate columns without
        double-counting, IsoQuant `#feature_id count`, Bambu `TXNAME GENEID
        sample` where the non-numeric GENEID is skipped).
    See signal_utils.load_support_table.
    """
    return load_support_table(path)


def load_calls(bed, reg=None, counts=None, canonical_only=True):
    """Per mode: (calls, has_support).

    calls = list of dicts {chrom,strand,tss,tts,reads}.  has_support is False
    when the input is a GTF tool with NO usable per-isoform read count (so the
    caller can OMIT it from read-support-stratified outputs but keep it in the
    distance-stratified ones).

    Auto-detects BED12 vs assembler GTF via the SHARED loader
    (signal_utils.load_isoforms_any) so a GTF tool produces the same per-isoform
    geometry (chrom/strand/start/end/tss/tts) as a BED tool.  The BED path is
    byte-for-byte the same isoforms parse_isoforms produced before; only the GTF
    branch is new.
    """
    regions = [reg] if reg is not None else None
    isos = load_isoforms_any(bed, regions=regions,
                             canonical_only=canonical_only, counts=counts)
    out = []
    n_with = 0
    for iso in isos:
        rs = iso.get("read_support")
        reads = 0.0 if rs is None else float(rs)
        if iso.get("has_support"):
            n_with += 1
        out.append({"chrom": iso["chrom"], "strand": iso["strand"],
                    "tss": iso["tss"], "tts": iso["tts"], "reads": reads})
    # A tool "has support" only if at least one isoform got a usable count
    # (BED tools always do via inline col16/col5; a GTF tool with no --counts
    # never does).  n=0 with isoforms present => no support source.
    has_support = (n_with > 0) or (not out)
    return out, has_support


def load_calls_full(bed, reg=None, counts=None, canonical_only=True):
    """Per mode: (isoforms, has_support) keeping the FULL JC geometry.

    Unlike :func:`load_calls` (which keeps only chrom/strand/tss/tts/reads for
    the single-end series), this preserves the junction chain, exon count and
    genomic span so the PAIRED metric can JC-group + dedup exactly the way
    ted_end_precision.compute_jc_deduplicated_precision_recall does.

    Each returned isoform dict has the keys that function consumes:
        chrom, strand, tss, tts, junctions, n_exons, start, end,
        read_support (float; 0.0 when the inline support was missing)
    plus has_support (False for a GTF tool with no --counts, so the caller can
    OMIT it from the paired support figures — same policy as the single-end
    path).
    """
    regions = [reg] if reg is not None else None
    isos = load_isoforms_any(bed, regions=regions,
                             canonical_only=canonical_only, counts=counts)
    out = []
    n_with = 0
    for iso in isos:
        rs = iso.get("read_support")
        reads = 0.0 if rs is None else float(rs)
        if iso.get("has_support"):
            n_with += 1
        out.append({
            "chrom": iso["chrom"], "strand": iso["strand"],
            "tss": iso["tss"], "tts": iso["tts"],
            "junctions": tuple(iso.get("junctions", ())),
            "n_exons": int(iso.get("n_exons", 1)),
            "start": iso["start"], "end": iso["end"],
            "read_support": reads,
        })
    has_support = (n_with > 0) or (not out)
    return out, has_support


def call_pos_index(calls, end_type):
    """{(chrom,strand): sorted [positions]} of called ends (for peak recovery)."""
    idx = defaultdict(list)
    for c in calls:
        idx[(c["chrom"], c["strand"])].append(c[end_type])
    for k in idx:
        idx[k].sort()
    return idx


def _peak_recovered_by_calls(peak, sorted_calls, window):
    """MANY-TO-ONE recall: True iff ANY call is within `window` of this peak interval.

    Scan bounds mirror ted_end_precision._many_to_one_match_counts (supported
    zone [lo-W, (hi-1)+W]) so the recall numerator matches the headline.
    """
    if not sorted_calls:
        return False
    s = int(peak)
    scan_lo = min(peak.lo, s) - window
    scan_hi = max(peak.hi, s + 1) + window
    i = bisect_left(sorted_calls, scan_lo)
    while i < len(sorted_calls):
        c = sorted_calls[i]
        if c >= scan_hi:
            break
        if _match_distance(c, peak) <= window:
            return True
        i += 1
    return False


def compute(calls, peaks, annot, end_type, window):
    """Return all stratified series for one mode + end_type.
    peaks: {(chrom,strand): ([Peak], scores)}; annot: {(chrom,strand): sorted ann pos}."""
    callidx = call_pos_index(calls, end_type)
    # pre-extents per peak bucket
    extents = {k: _peak_half_extents(v[0]) for k, v in peaks.items()}

    # per-call: matched? + dist-to-annotation  (MANY-TO-ONE membership)
    call_recs = []
    for c in calls:
        key = (c["chrom"], c["strand"])
        pk = peaks.get(key)
        if pk:
            ml, mr = extents[key]
            matched = _near_interval(c[end_type], pk[0], ml, mr, window)
        else:
            matched = False
        dann = _dist(c[end_type], annot.get(key, []))
        call_recs.append((c["reads"], matched, dann))

    # per-peak: recovered? + score + dist-to-annotation
    peak_recs = []
    for key, (peak_list, score_list) in peaks.items():
        calls_here = callidx.get(key, [])
        ann_here = annot.get(key, [])
        for p, s in zip(peak_list, score_list):
            rec = _peak_recovered_by_calls(p, calls_here, window)
            dann = _dist(int(p), ann_here)
            peak_recs.append((s, rec, dann))

    res = {}
    # (1) support sweep: cumulative min-reads threshold
    Psup, Rsup = [], []
    n_peaks_total = len(peak_recs)
    for k in KS:
        sub = [r for r in call_recs if r[0] >= k]
        Psup.append(sum(1 for r in sub if r[1]) / len(sub) if sub else None)
        # recall: peaks recovered by a call with reads>=k
        sub_idx = defaultdict(list)
        for c in calls:
            if c["reads"] >= k:
                sub_idx[(c["chrom"], c["strand"])].append(c[end_type])
        for kk in sub_idx:
            sub_idx[kk].sort()
        rec_k = 0
        for key, (peak_list, _s) in peaks.items():
            ch = sub_idx.get(key, [])
            rec_k += sum(1 for p in peak_list if _peak_recovered_by_calls(p, ch, window))
        Rsup.append(rec_k / n_peaks_total if n_peaks_total else None)
    res["support"] = {"x": KS, "P": Psup, "R": Rsup}

    # (2) recall by ground-truth peak score (quartile bins)
    scores = sorted(s for s, _, _ in peak_recs)
    if scores:
        qs = [scores[int(q * (len(scores) - 1))] for q in (0.0, 0.25, 0.5, 0.75, 1.0)]
        edges = sorted(set(qs))
        bins = list(zip(edges[:-1], edges[1:])) or [(edges[0], edges[0] + 1)]
        xb, Rb = [], []
        for lo, hi in bins:
            grp = [r for r in peak_recs if lo <= r[0] < hi or (hi == bins[-1][1] and r[0] >= lo)]
            xb.append(f"{lo:.0f}-{hi:.0f}")
            Rb.append(sum(1 for r in grp if r[1]) / len(grp) if grp else None)
        res["gtscore"] = {"x": xb, "R": Rb}
    else:
        res["gtscore"] = {"x": [], "R": []}

    # (3) distance-from-annotation: precision by CALL dist, recall by PEAK dist
    Pd, Rd = [], []
    for lo, hi in DIST_BINS:
        cg = [r for r in call_recs if lo <= r[2] < hi]
        pg = [r for r in peak_recs if lo <= r[2] < hi]
        Pd.append(sum(1 for r in cg if r[1]) / len(cg) if cg else None)
        Rd.append(sum(1 for r in pg if r[1]) / len(pg) if pg else None)
    res["dist"] = {"x": DIST_LABELS, "P": Pd, "R": Rd}
    return res


# ── Paired-end (TSS+TTS) stratified metric ────────────────────────────────────
#
# These reuse the CANONICAL paired machinery from ted_end_precision so the
# numbers reconcile with precision_recall_summary.tsv's paired_dedup_precision.
# Nothing here touches the single-end `compute` path above.

def compute_paired_precision_support(isoforms, peaks_5_intervals,
                                     peaks_3_intervals, annotated_ends,
                                     window):
    """PAIRED PRECISION stratified by read support, per KS threshold.

    For each k in KS: keep isoforms with read_support >= k, then CALL
    ted_end_precision.compute_jc_deduplicated_precision_recall on that subset
    and read its ``paired_dedup_precision``.  This is the JC-deduplicated,
    both-ends-hit precision:

        group by junction chain (single-exon = overlap-grouped on strand),
        dedup to unique (TSS-peak, TTS-peak) pairs,
        paired_precision = (unique pairs with BOTH ends matched) / (total isoforms).

    Because the value comes from THE canonical function (not a re-implementation),
    paired_precision at k=1 == precision_recall_summary's paired_dedup_precision
    on the same (canonical-filtered, region-restricted) input — the consistency
    contract.  The peak matcher is _nearest_annot (within-W-of-interval, the
    SAME predicate the paired headline uses).

    `peaks_5_intervals` / `peaks_3_intervals` are {(chrom,strand): [(start,end)]}
    (ted_end_precision.parse_peaks_bed shape); `annotated_ends` is the GTF-ends
    dict (unused when peaks are provided, but the function signature requires it).

    Returns {"x": KS, "P": [...]}.
    """
    P = []
    for k in KS:
        sub = [iso for iso in isoforms if iso["read_support"] >= k]
        if not sub:
            P.append(None)
            continue
        # compute_jc_deduplicated_precision_recall tags isoforms with a
        # transient "_se" marker and pops it; pass shallow copies so repeated
        # KS calls never see stale state and the caller's dicts stay clean.
        sub_copy = [dict(iso) for iso in sub]
        res = compute_jc_deduplicated_precision_recall(
            sub_copy, annotated_ends, [], window=window,
            peaks_5prime=peaks_5_intervals, peaks_3prime=peaks_3_intervals)
        P.append(res.get("paired_dedup_precision"))
    return {"x": list(KS), "P": P}


def build_validated_truth_pairs(annot_transcripts, peaks_5_intervals,
                                peaks_3_intervals, window, dedup_tol=None):
    """Tool-INDEPENDENT validated reference (TSS, TTS) end-pairs.

    From every reference GTF transcript build its (TSS, TTS) end-pair, then KEEP
    only pairs whose TSS is within `window` of a CAGE peak INTERVAL AND whose TTS
    is within `window` of a dRNA peak INTERVAL (orthogonally-validated reference
    end-pairs — the matcher is _nearest_annot, the same within-W-of-interval
    predicate everywhere else).  Deduplicate the survivors by (chrom, strand,
    TSS-bucket, TTS-bucket) where a bucket = position // (dedup_tol) so two truth
    pairs within `dedup_tol` bp on BOTH ends collapse to one.

    `dedup_tol` defaults to `window` (so pairs within a window on both ends are
    one validated target).  Returns a list of dicts
        {chrom, strand, tss, tts}
    one per UNIQUE validated truth pair.  Canonical filtering is the caller's job
    (parse_gtf_transcripts is fed region/canonical-restricted peaks; the GTF is
    parsed canonical-only at the call site).
    """
    if dedup_tol is None:
        dedup_tol = window
    tol = max(1, int(dedup_tol))
    seen = set()
    pairs = []
    for tx in annot_transcripts:
        chrom, strand = tx["chrom"], tx["strand"]
        if not is_canonical_chrom(chrom):
            continue
        tss, tts = tx["tss"], tx["tts"]
        p5 = peaks_5_intervals.get((chrom, strand), [])
        p3 = peaks_3_intervals.get((chrom, strand), [])
        if _nearest_annot(tss, p5, window) is None:
            continue
        if _nearest_annot(tts, p3, window) is None:
            continue
        key = (chrom, strand, tss // tol, tts // tol)
        if key in seen:
            continue
        seen.add(key)
        pairs.append({"chrom": chrom, "strand": strand, "tss": tss, "tts": tts})
    return pairs


def compute_paired_recall_support(isoforms, truth_pairs, window):
    """PAIRED RECALL stratified by read support, per KS threshold.

    A validated truth pair (chrom, strand, tss, tts) is RECOVERED at threshold k
    when SOME isoform with read_support >= k has BOTH:
        |iso.tss - truth.tss| <= window  AND  |iso.tts - truth.tts| <= window
    (scalar within-W on each end; the truth ends are already validated reference
    positions).  Then
        paired_recall(k) = recovered truth pairs / total validated truth pairs.

    The denominator (len(truth_pairs)) is tool-independent, so recall lines are
    directly comparable across tools.  Returns {"x": KS, "R": [...]}.
    """
    n_total = len(truth_pairs)
    if n_total == 0:
        return {"x": list(KS), "R": [None] * len(KS)}

    # Index isoform TSS positions by (chrom, strand) so each truth pair only
    # probes nearby ends.  Each entry carries (tss, tts, reads) so we can apply
    # the support threshold + the paired (both-ends) predicate in one scan.
    from bisect import bisect_left as _bl
    by_key = defaultdict(list)
    for iso in isoforms:
        by_key[(iso["chrom"], iso["strand"])].append(
            (int(iso["tss"]), int(iso["tts"]), iso["read_support"]))
    sorted_keys = {}
    for key, lst in by_key.items():
        lst.sort(key=lambda t: t[0])           # sort by TSS for bisect
        sorted_keys[key] = [t[0] for t in lst]  # parallel TSS array

    R = []
    for k in KS:
        recovered = 0
        for tp in truth_pairs:
            key = (tp["chrom"], tp["strand"])
            ends = by_key.get(key)
            if not ends:
                continue
            tss_arr = sorted_keys[key]
            lo = _bl(tss_arr, tp["tss"] - window)
            hit = False
            j = lo
            n = len(ends)
            while j < n:
                etss = ends[j][0]
                if etss > tp["tss"] + window:
                    break
                etts, ereads = ends[j][1], ends[j][2]
                if ereads >= k and abs(etts - tp["tts"]) <= window:
                    hit = True
                    break
                j += 1
            if hit:
                recovered += 1
        R.append(recovered / n_total)
    return {"x": list(KS), "R": R}


# ── Plotting ──────────────────────────────────────────────────────────────────

def _plot_grid(data, key, sub, fname, title, xlabel, xlog=False, categorical=False):
    """data: {mode: {'tss': res, 'tts': res}}. sub in ('P','R')."""
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for col, (which, ct) in enumerate([("tss", "5' TSS (vs CAGE)"), ("tts", "3' TTS (vs dRNA)")]):
        a = ax[col]
        xref = None
        for _m, d in data.items():
            xx = d[which].get(key, {}).get("x")
            if xx:
                xref = xx; break
        for i, (mode, d) in enumerate(data.items()):
            series = d[which].get(key, {})
            y = series.get(sub); x = series.get("x")
            if not y or not x:
                continue
            xs = list(range(len(x))) if categorical else x
            a.plot(xs, y, "-o", ms=4, color=PALETTE[i % len(PALETTE)], label=mode)
        a.set_ylim(-0.02, 1.02); a.grid(alpha=0.3)
        if categorical and xref:
            a.set_xticks(range(len(xref))); a.set_xticklabels(xref, fontsize=8, rotation=20)
        elif xlog:
            a.set_xscale("log")
        a.set_title(ct, fontsize=11); a.set_xlabel(xlabel, fontsize=10)
        if col == 0:
            a.set_ylabel("precision" if sub == "P" else "recall", fontsize=11)
            a.legend(fontsize=8, loc="best")
    fig.suptitle(title, fontsize=12); fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(fname, dpi=140); plt.close(fig); print("wrote", fname)


def _plot_pr_tradeoff(data, fname):
    """PR-tradeoff: per tool a (recall, precision) curve traced by support thresholds.

    One point per KS threshold; two panels (5' TSS, 3' TTS).  This is the
    classic precision-recall tradeoff view of the support sweep.
    """
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for col, (which, ct) in enumerate([("tss", "5' TSS (vs CAGE)"), ("tts", "3' TTS (vs dRNA)")]):
        a = ax[col]
        for i, (mode, d) in enumerate(data.items()):
            sup = d[which].get("support", {})
            P = sup.get("P"); R = sup.get("R"); ks = sup.get("x")
            if not P or not R:
                continue
            pts = [(rr, pp) for rr, pp in zip(R, P) if rr is not None and pp is not None]
            if not pts:
                continue
            xr = [p[0] for p in pts]; yp = [p[1] for p in pts]
            a.plot(xr, yp, "-o", ms=4, color=PALETTE[i % len(PALETTE)], label=mode)
        a.set_xlim(-0.02, 1.02); a.set_ylim(-0.02, 1.02); a.grid(alpha=0.3)
        a.set_title(ct, fontsize=11); a.set_xlabel("recall", fontsize=10)
        if col == 0:
            a.set_ylabel("precision", fontsize=11)
            a.legend(fontsize=8, loc="best")
    fig.suptitle("Precision–recall tradeoff across reads / isoform thresholds",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(fname, dpi=140); plt.close(fig); print("wrote", fname)


def _plot_paired_single(data, sub, fname, title, subtitle):
    """SINGLE-panel paired (TSS+TTS) metric vs reads/isoform (log KS).

    data: {mode: {... "paired": {"support": {"x", "P", "R"}}}}; sub in ('P','R').
    One line per tool.  Tools whose paired support series is empty (no read
    support) are simply absent (already omitted upstream).
    """
    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    n_lines = 0
    for i, (mode, d) in enumerate(data.items()):
        series = d.get("paired", {}).get("support", {})
        y = series.get(sub); x = series.get("x")
        if not y or not x:
            continue
        ax.plot(x, y, "-o", ms=4, color=PALETTE[i % len(PALETTE)], label=mode)
        n_lines += 1
    ax.set_ylim(-0.02, 1.02); ax.grid(alpha=0.3); ax.set_xscale("log")
    ax.set_xlabel("reads / isoform", fontsize=11)
    ax.set_ylabel("paired precision" if sub == "P" else "paired recall",
                  fontsize=11)
    if n_lines:
        ax.legend(fontsize=8, loc="best")
    fig.suptitle(title, fontsize=12, y=0.98)
    ax.set_title(subtitle, fontsize=8.5, color="#444")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(fname, dpi=140); plt.close(fig); print("wrote", fname)


def _plot_paired_pr_tradeoff(data, fname):
    """OPTIONAL combined paired PR-tradeoff: per tool a (recall, precision)
    curve traced by the support thresholds (single panel)."""
    fig, ax = plt.subplots(1, 1, figsize=(6.5, 5.5))
    n_lines = 0
    for i, (mode, d) in enumerate(data.items()):
        sup = d.get("paired", {}).get("support", {})
        P = sup.get("P"); R = sup.get("R")
        if not P or not R:
            continue
        pts = [(rr, pp) for rr, pp in zip(R, P) if rr is not None and pp is not None]
        if not pts:
            continue
        xr = [p[0] for p in pts]; yp = [p[1] for p in pts]
        ax.plot(xr, yp, "-o", ms=4, color=PALETTE[i % len(PALETTE)], label=mode)
        n_lines += 1
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02); ax.grid(alpha=0.3)
    ax.set_xlabel("paired recall", fontsize=11)
    ax.set_ylabel("paired precision", fontsize=11)
    if n_lines:
        ax.legend(fontsize=8, loc="best")
    fig.suptitle("Paired (TSS+TTS) precision-recall tradeoff across reads/isoform",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(fname, dpi=140); plt.close(fig); print("wrote", fname)


PAIRED_PREC_SUBTITLE = (
    "JC-deduplicated: per junction chain, unique (TSS-peak,TTS-peak) pairs with "
    "BOTH ends hit / total isoforms (>=k reads). Mirrors paired_dedup_precision.")
PAIRED_RECALL_SUBTITLE = (
    "vs orthogonally-validated reference end-pairs (GTF TSS within-50 of CAGE AND "
    "TTS within-50 of dRNA, deduped); recovered by a >=k-read isoform / total.")


def _parse_counts_args(entries):
    """--counts LABEL:path tokens -> {label: path}."""
    out = {}
    for e in (entries or []):
        if ":" not in e:
            continue
        label, path = e.split(":", 1)
        out[label] = path
    return out


# ── TSV checkpoint <-> data round-trip (compute path + --replot share this) ──

_END_TO_WHICH = {"5prime": "tss", "3prime": "tts"}
_WHICH_TO_END = {v: k for k, v in _END_TO_WHICH.items()}


def _data_to_rows(data):
    """Flatten the plotting `data` dict into checkpoint TSV rows.

    Skips a tool's `support` series when it is empty (no read-support source) so
    the checkpoint records EXACTLY which tools were omitted from Fig A/A2.
    """
    rows = []
    for lab, d in data.items():
        for which in ("tss", "tts"):
            r = d.get(which, {})
            end = _WHICH_TO_END[which]
            sup = r.get("support", {})
            if sup.get("x"):
                for k, P in zip(sup["x"], sup.get("P", [])):
                    rows.append([lab, end, "support", k, "precision", P])
                for k, R in zip(sup["x"], sup.get("R", [])):
                    rows.append([lab, end, "support", k, "recall", R])
            gt = r.get("gtscore", {})
            for x, R in zip(gt.get("x", []), gt.get("R", [])):
                rows.append([lab, end, "gtscore", x, "recall", R])
            di = r.get("dist", {})
            for x, P in zip(di.get("x", []), di.get("P", [])):
                rows.append([lab, end, "dist", x, "precision", P])
            for x, R in zip(di.get("x", []), di.get("R", [])):
                rows.append([lab, end, "dist", x, "recall", R])
        # PAIRED (TSS+TTS) support series — end='paired' keeps the TSV tidy.
        # Skipped entirely when empty (tool omitted from paired support figs).
        psup = d.get("paired", {}).get("support", {})
        if psup.get("x"):
            for k, P in zip(psup["x"], psup.get("P", [])):
                rows.append([lab, "paired", "support", k, "paired_precision", P])
            for k, R in zip(psup["x"], psup.get("R", [])):
                rows.append([lab, "paired", "support", k, "paired_recall", R])
    return rows


def _write_tsv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["mode", "end", "stratifier", "bin", "metric", "value"])
        for r in rows:
            w.writerow([r[0], r[1], r[2], r[3], r[4],
                        "" if r[5] is None else f"{r[5]:.4f}"])


def _rows_to_data(tsv_path):
    """Reconstruct the plotting `data` dict from a checkpoint TSV (replot path).

    Preserves bin order as first-seen and re-derives the support `x` as int KS
    so the regenerated figures match the compute-path figures exactly.
    """
    # nested: data[mode][which][stratifier] -> ordered {bin: {metric: value}}
    raw = defaultdict(lambda: {"tss": defaultdict(lambda: {}),
                               "tts": defaultdict(lambda: {})})
    # paired support series collected separately: paired_raw[mode] ->
    # {"_order": [...], "_d": {bin: {metric: value}}}.
    paired_raw = defaultdict(lambda: {"_order": [], "_d": {}})
    order = []
    with open(tsv_path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for rec in rdr:
            mode = rec["mode"]
            if mode not in order:
                order.append(mode)
            strat = rec["stratifier"]
            b = rec["bin"]
            metric = rec["metric"]
            val = rec["value"]
            v = None if val == "" else float(val)
            if rec["end"] == "paired":
                slot = paired_raw[mode]
                if b not in slot["_d"]:
                    slot["_order"].append(b)
                    slot["_d"][b] = {}
                slot["_d"][b][metric] = v
                continue
            which = _END_TO_WHICH.get(rec["end"])
            if which is None:
                continue
            slot = raw[mode][which].setdefault(strat, {"_order": [], "_d": {}})
            if b not in slot["_d"]:
                slot["_order"].append(b)
                slot["_d"][b] = {}
            slot["_d"][b][metric] = v

    data = {}
    for mode in order:
        d = {}
        for which in ("tss", "tts"):
            res = {}
            buckets = raw[mode][which]
            # support
            sup = buckets.get("support")
            if sup and sup["_order"]:
                xs = [int(b) for b in sup["_order"]]
                res["support"] = {
                    "x": xs,
                    "P": [sup["_d"][b].get("precision") for b in sup["_order"]],
                    "R": [sup["_d"][b].get("recall") for b in sup["_order"]],
                }
            else:
                res["support"] = {"x": [], "P": [], "R": []}
            # gtscore (recall only)
            gt = buckets.get("gtscore")
            if gt and gt["_order"]:
                res["gtscore"] = {
                    "x": list(gt["_order"]),
                    "R": [gt["_d"][b].get("recall") for b in gt["_order"]],
                }
            else:
                res["gtscore"] = {"x": [], "R": []}
            # dist
            di = buckets.get("dist")
            if di and di["_order"]:
                res["dist"] = {
                    "x": list(di["_order"]),
                    "P": [di["_d"][b].get("precision") for b in di["_order"]],
                    "R": [di["_d"][b].get("recall") for b in di["_order"]],
                }
            else:
                res["dist"] = {"x": [], "P": [], "R": []}
            d[which] = res
        # paired support series (int KS x; P=paired_precision, R=paired_recall)
        psup = paired_raw.get(mode)
        if psup and psup["_order"]:
            xs = [int(b) for b in psup["_order"]]
            d["paired"] = {"support": {
                "x": xs,
                "P": [psup["_d"][b].get("paired_precision") for b in psup["_order"]],
                "R": [psup["_d"][b].get("paired_recall") for b in psup["_order"]],
            }}
        else:
            d["paired"] = {"support": {"x": [], "P": [], "R": []}}
        data[mode] = d
    return data


def _render_all(data, out):
    """Draw every PNG from the plotting `data` dict (compute + replot share this)."""
    _plot_grid(data, "support", "P", out / "stratified_precision_by_reads.png",
               "Precision vs reads / isoform", "reads / isoform", xlog=True)
    _plot_grid(data, "support", "R", out / "stratified_recall_by_reads.png",
               "Recall vs reads / isoform", "reads / isoform", xlog=True)
    _plot_pr_tradeoff(data, out / "pr_tradeoff_support.png")
    _plot_grid(data, "gtscore", "R", out / "stratified_recall_by_gtscore.png",
               "Recall vs ground-truth peak score/TPM (quartiles)", "ground-truth peak score bin",
               categorical=True)
    _plot_grid(data, "dist", "P", out / "stratified_precision_by_annotdist.png",
               "Precision vs CALL distance from annotation", "distance from annotation (bp)",
               categorical=True)
    _plot_grid(data, "dist", "R", out / "stratified_recall_by_annotdist.png",
               "Recall vs PEAK distance from annotation", "distance from annotation (bp)",
               categorical=True)
    # PAIRED (TSS+TTS) single-panel figures + optional combined PR-tradeoff.
    _plot_paired_single(data, "P", out / "paired_precision_by_reads.png",
                        "Paired (TSS+TTS) precision vs reads / isoform",
                        PAIRED_PREC_SUBTITLE)
    _plot_paired_single(data, "R", out / "paired_recall_by_reads.png",
                        "Paired (TSS+TTS) recall vs reads / isoform",
                        PAIRED_RECALL_SUBTITLE)
    _plot_paired_pr_tradeoff(data, out / "paired_pr_tradeoff.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--beds", nargs="+",
                    help="per-mode isoforms BED12 (FLAIR) or assembler GTF")
    ap.add_argument("--labels", nargs="+",
                    help="mode labels (aligned with --beds)")
    ap.add_argument("--counts", nargs="+", default=[],
                    help="optional LABEL:path per-tool read-count tables")
    ap.add_argument("--cage")
    ap.add_argument("--drna")
    ap.add_argument("--gtf")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--region", default=None)
    ap.add_argument("--window", type=int, default=50)
    ap.add_argument("--novel-dist", type=int, default=100)
    ap.add_argument("--all-contigs", action="store_true",
                    help="keep ALL contigs (chrM + GL*/KI*); default = canonical only")
    ap.add_argument("--replot", action="store_true",
                    help="re-draw the PNGs from the existing --outdir/stratified_pr.tsv "
                         "checkpoint with ZERO recompute (no isoform/GTF/peak parsing).")
    a = ap.parse_args()
    out = Path(a.outdir)

    # ── REPLOT: read the checkpoint, regenerate identical figures ───────────
    if a.replot:
        tsv = out / "stratified_pr.tsv"
        if not tsv.exists():
            sys.exit(f"--replot: checkpoint not found: {tsv}")
        print(f"[replot] reading checkpoint {tsv} (no recompute)", flush=True)
        data = _rows_to_data(tsv)
        out.mkdir(parents=True, exist_ok=True)
        _render_all(data, out)
        print("STRATIFIED_PR_DONE")
        return

    # ── COMPUTE path ────────────────────────────────────────────────────────
    if not (a.beds and a.labels and a.cage and a.drna and a.gtf):
        sys.exit("compute mode requires --beds --labels --cage --drna --gtf "
                 "(or pass --replot to redraw from the checkpoint)")
    assert len(a.beds) == len(a.labels), "--beds and --labels must align"
    canonical_only = not a.all_contigs
    reg = parse_region(a.region)
    counts_map_paths = _parse_counts_args(a.counts)
    out.mkdir(parents=True, exist_ok=True)

    print("[load] annotation + truth peaks ...", flush=True)
    annTss, annTts = parse_gtf_ends(a.gtf, reg, canonical_only)
    cage = load_peaks(a.cage, reg, canonical_only)
    drna = load_peaks(a.drna, reg, canonical_only)
    n_cage = sum(len(v[0]) for v in cage.values())
    n_drna = sum(len(v[0]) for v in drna.values())
    print(f"  {n_cage} CAGE peaks, {n_drna} dRNA peaks (canonical_only={canonical_only}, "
          f"region={a.region})", flush=True)

    # ── PAIRED prep: truth peaks as ted_end_precision INTERVAL dicts +
    # the tool-independent validated truth-pair set (built ONCE). ──────────
    ted_regions = [(reg[0], reg[1], reg[2])] if reg is not None else None
    print("[load] paired truth: interval peaks + validated truth pairs ...",
          flush=True)
    cage_iv = _ted_parse_peaks_bed(a.cage, ted_regions, canonical_only=canonical_only)
    drna_iv = _ted_parse_peaks_bed(a.drna, ted_regions, canonical_only=canonical_only)
    ted_annot_ends = _ted_parse_gtf_ends(a.gtf, ted_regions)
    annot_tx = _ted_parse_gtf_transcripts(a.gtf, ted_regions)
    truth_pairs = build_validated_truth_pairs(
        annot_tx, cage_iv, drna_iv, a.window)
    print(f"  {len(annot_tx)} GTF transcripts -> {len(truth_pairs)} validated "
          f"(CAGE&dRNA) truth pairs (paired-recall denominator)", flush=True)

    data = {}
    omitted_support = []   # (label, reason) — tools dropped from Fig A/A2
    for bed, lab in zip(a.beds, a.labels):
        counts = load_counts_map(counts_map_paths.get(lab))
        calls, has_support = load_calls(bed, reg, counts, canonical_only)
        rtss = compute(calls, cage, annTss, "tss", a.window)
        rtts = compute(calls, drna, annTts, "tts", a.window)

        # PAIRED: load the FULL JC geometry (junctions kept) and compute the
        # JC-dedup paired precision + the validated-truth-pair paired recall.
        isos_full, _has_full = load_calls_full(bed, reg, counts, canonical_only)
        ppair = compute_paired_precision_support(
            isos_full, cage_iv, drna_iv, ted_annot_ends, a.window)
        rpair = compute_paired_recall_support(isos_full, truth_pairs, a.window)
        paired = {"support": {"x": list(KS), "P": ppair["P"], "R": rpair["R"]}}

        if not has_support:
            # NO SILENT DROP: GTF tool without --counts cannot be thresholded by
            # read support.  Keep its distance/gt-score series (Fig B), but blank
            # the support series so it is omitted from Figs A/A2 + PR-tradeoff and
            # from the PAIRED support figures (same policy).
            is_gtf = is_gtf_path(bed)
            reason = ("GTF input with no --counts (read support unavailable)"
                      if is_gtf else "no per-isoform read support")
            omitted_support.append((lab, reason))
            for r in (rtss, rtts):
                r["support"] = {"x": [], "P": [], "R": []}
            paired = {"support": {"x": [], "P": [], "R": []}}
            print(f"  OMIT {lab} from read-support figures (Fig A/A2/PR-tradeoff "
                  f"+ PAIRED): {reason}; KEPT in distance figures (Fig B)",
                  file=sys.stderr, flush=True)
        data[lab] = {"tss": rtss, "tts": rtts, "paired": paired}
        n_counts = "inline" if counts is None else f"{len(counts)} ids"
        pp1 = paired["support"]["P"][0] if paired["support"]["P"] else None
        pr1 = paired["support"]["R"][0] if paired["support"]["R"] else None
        print(f"  {lab}: {len(calls)} isoforms (counts={n_counts}, "
              f"support={'yes' if has_support else 'NO'}); "
              f"paired@k1 P={pp1} R={pr1}", flush=True)

    if omitted_support:
        print("[omissions] read-support-stratified figures (Fig A/A2/PR-tradeoff) "
              "exclude: " + ", ".join(f"{l} ({r})" for l, r in omitted_support),
              file=sys.stderr, flush=True)

    rows = _data_to_rows(data)
    _write_tsv(out / "stratified_pr.tsv", rows)
    print("wrote", out / "stratified_pr.tsv")

    _render_all(data, out)
    print("STRATIFIED_PR_DONE")


if __name__ == "__main__":
    main()
