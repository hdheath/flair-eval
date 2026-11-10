#!/usr/bin/env python3
"""Paired-end (TSS+TTS) stratified metric for stratified_pr.py (added 2026-06-22).

stratified_pr.py gained, per tool, per read-support threshold k:
  (1) PAIRED PRECISION(k) — JC-deduplicated, computed by CALLING
      ted_end_precision.compute_jc_deduplicated_precision_recall on the
      >=k-read isoform subset and reading paired_dedup_precision.  At k=1 this
      MUST equal that function's paired_dedup_precision on the same input.
  (2) PAIRED RECALL(k) — vs a tool-independent VALIDATED truth-pair set:
      reference GTF transcripts -> (TSS,TTS), kept only when TSS within-50 of a
      CAGE peak AND TTS within-50 of a dRNA peak, deduped; recovered when a
      >=k-read isoform matches BOTH ends within-50.

These tests pin:
  * paired precision @ k=1 == compute_jc_deduplicated_precision_recall's
    paired_dedup_precision (the reconciliation contract), on a synthetic set;
  * paired precision is monotone-non-decreasing as k rises on a clean set;
  * validated truth-pair construction keeps ONLY both-ends-validated pairs and
    dedups within tolerance;
  * paired recall numerator/denominator follow the documented formula, and the
    support threshold gates recovery;
  * a no-support (GTF NO_COUNTS) tool is OMITTED from the paired support series.

conftest.py puts bin/ on sys.path (so `from evaluation.<...>` works).
"""
from __future__ import annotations

from evaluation import stratified_pr as S
from evaluation.ted_end_precision import (
    compute_jc_deduplicated_precision_recall, parse_peaks_bed,
    parse_gtf_ends, parse_gtf_transcripts)

W = 50


def _iso(chrom, strand, tss, tts, junctions, reads, n_exons=2):
    """A full-geometry isoform record (load_calls_full shape)."""
    start = min(tss, tts)
    end = max(tss, tts)
    return {"chrom": chrom, "strand": strand, "tss": tss, "tts": tts,
            "junctions": tuple(junctions), "n_exons": n_exons,
            "start": start, "end": end, "read_support": float(reads)}


# ── peaks: CAGE summit 1000 (+), dRNA summit 2000 (+); a 2nd pair at 5000/6000
def _cage_intervals():
    # ted_end_precision.parse_peaks_bed shape: {(chrom,strand): [(start,end)]}
    return {("chr1", "+"): [(950, 1050), (4950, 5050)]}


def _drna_intervals():
    return {("chr1", "+"): [(1950, 2050), (5950, 6050)]}


# ════════════════════════════════════════════════════════════════════════
# (1) PAIRED PRECISION @ k=1 == canonical paired_dedup_precision
# ════════════════════════════════════════════════════════════════════════

def test_paired_precision_k1_mirrors_canonical():
    cage = _cage_intervals(); drna = _drna_intervals()
    annot_ends = {}   # peaks provided -> GTF ends unused by the paired path
    # 4 isoforms, distinct junction chains:
    #   A: TSS@1000 hits, TTS@2000 hits        -> both hit
    #   B: TSS@1005 hits, TTS@2010 hits        -> both hit (distinct JC)
    #   C: TSS@5000 hits, TTS@9999 (off-peak)  -> only TSS hits -> NOT a pair
    #   D: TSS@8888 (off), TTS@6000 hits       -> only TTS hits -> NOT a pair
    isos = [
        _iso("chr1", "+", 1000, 2000, [(1100, 1200)], 5),
        _iso("chr1", "+", 1005, 2010, [(1100, 1300)], 3),
        _iso("chr1", "+", 5000, 9999, [(5100, 5200)], 2),
        _iso("chr1", "+", 8888, 6000, [(8100, 8200)], 1),
    ]
    # CANONICAL paired_dedup_precision on the SAME isoforms (fresh copies).
    canon = compute_jc_deduplicated_precision_recall(
        [dict(i) for i in isos], annot_ends, [], window=W,
        peaks_5prime=cage, peaks_3prime=drna)["paired_dedup_precision"]

    ppair = S.compute_paired_precision_support(isos, cage, drna, annot_ends, W)
    k1 = ppair["P"][0]   # S.KS[0] == 1
    assert S.KS[0] == 1
    assert k1 == canon
    # 2 of 4 isoforms are both-hit unique pairs -> 0.5
    assert k1 == 0.5


def test_paired_precision_monotone_with_support():
    """On a set where higher-support isoforms are the clean (both-hit) ones,
    paired precision rises with k."""
    cage = _cage_intervals(); drna = _drna_intervals()
    isos = [
        _iso("chr1", "+", 1000, 2000, [(1100, 1200)], 10),   # both hit, 10 reads
        _iso("chr1", "+", 5000, 6000, [(5100, 5200)], 8),    # both hit, 8 reads
        _iso("chr1", "+", 9000, 9500, [(9100, 9200)], 1),    # off-peak, 1 read
        _iso("chr1", "+", 9300, 9800, [(9350, 9400)], 1),    # off-peak, 1 read
    ]
    ppair = S.compute_paired_precision_support(isos, cage, drna, {}, W)
    P = ppair["P"]
    # k=1: 4 isoforms, 2 both-hit -> 0.5 ; k>=2: only the 2 clean ones -> 1.0
    assert P[0] == 0.5
    # KS = [1,2,3,5,8,12,20,35,50]; at k=2 the two 1-read junk drop
    assert P[1] == 1.0
    # non-decreasing
    assert all((b is None or a is None or b >= a - 1e-12)
               for a, b in zip(P, P[1:]))


# ════════════════════════════════════════════════════════════════════════
# (2) VALIDATED TRUTH-PAIR construction
# ════════════════════════════════════════════════════════════════════════

def test_truth_pairs_keep_only_both_validated():
    cage = _cage_intervals(); drna = _drna_intervals()
    # 4 synthetic "transcripts" (parse_gtf_transcripts shape subset):
    annot_tx = [
        {"chrom": "chr1", "strand": "+", "tss": 1000, "tts": 2000},  # both ok
        {"chrom": "chr1", "strand": "+", "tss": 5000, "tts": 6000},  # both ok
        {"chrom": "chr1", "strand": "+", "tss": 1000, "tts": 9999},  # TTS off
        {"chrom": "chr1", "strand": "+", "tss": 8888, "tts": 2000},  # TSS off
        {"chrom": "chrM", "strand": "+", "tss": 1000, "tts": 2000},  # non-canon
    ]
    pairs = S.build_validated_truth_pairs(annot_tx, cage, drna, W)
    # only the 2 both-validated canonical pairs survive
    keys = {(p["chrom"], p["strand"], p["tss"], p["tts"]) for p in pairs}
    assert keys == {("chr1", "+", 1000, 2000), ("chr1", "+", 5000, 6000)}


def test_truth_pairs_dedup_within_tolerance():
    cage = _cage_intervals(); drna = _drna_intervals()
    # two transcripts within tolerance on BOTH ends collapse to one pair
    annot_tx = [
        {"chrom": "chr1", "strand": "+", "tss": 1000, "tts": 2000},
        {"chrom": "chr1", "strand": "+", "tss": 1010, "tts": 2010},  # within W on both
    ]
    pairs = S.build_validated_truth_pairs(annot_tx, cage, drna, W, dedup_tol=W)
    assert len(pairs) == 1


# ════════════════════════════════════════════════════════════════════════
# (3) PAIRED RECALL formula + support gating
# ════════════════════════════════════════════════════════════════════════

def test_paired_recall_formula_and_threshold():
    truth = [
        {"chrom": "chr1", "strand": "+", "tss": 1000, "tts": 2000},
        {"chrom": "chr1", "strand": "+", "tss": 5000, "tts": 6000},
    ]
    # isoform 1 recovers truth pair 1 only with >=5 reads; nothing recovers pair 2
    isos = [
        _iso("chr1", "+", 1010, 1990, [(1100, 1200)], 5),   # both within W of pair1
        _iso("chr1", "+", 5000, 9999, [(5100, 5200)], 9),   # TTS off -> not paired
    ]
    rpair = S.compute_paired_recall_support(isos, truth, W)
    R = rpair["R"]
    assert S.KS[:4] == [1, 2, 3, 5]
    # k=1,2,3: pair1 recovered (1/2 = 0.5); pair2 never
    assert R[0] == 0.5 and R[1] == 0.5 and R[2] == 0.5
    # k=5: still recovered (iso has 5 reads) -> 0.5
    assert R[3] == 0.5
    # k=8: iso1 (5 reads) drops -> 0 recovered
    assert R[4] == 0.0


def test_paired_recall_requires_both_ends():
    truth = [{"chrom": "chr1", "strand": "+", "tss": 1000, "tts": 2000}]
    # an isoform that matches ONLY the TSS must NOT recover the pair
    isos = [_iso("chr1", "+", 1000, 9000, [(1100, 1200)], 10)]
    R = S.compute_paired_recall_support(isos, truth, W)["R"]
    assert all(r == 0.0 for r in R)
    # add a TTS match -> recovered
    isos.append(_iso("chr1", "+", 1000, 2000, [(1100, 1300)], 10))
    R2 = S.compute_paired_recall_support(isos, truth, W)["R"]
    assert R2[0] == 1.0


def test_paired_recall_empty_truth():
    R = S.compute_paired_recall_support([], [], W)["R"]
    assert R == [None] * len(S.KS)


# ════════════════════════════════════════════════════════════════════════
# (4) NO-COUNTS tool OMITTED from paired support series (via _data_to_rows)
# ════════════════════════════════════════════════════════════════════════

def test_no_support_tool_omitted_from_paired_rows():
    data = {
        "FLAIR": {"tss": {"support": {"x": [], "P": [], "R": []},
                          "gtscore": {"x": [], "R": []},
                          "dist": {"x": [], "P": [], "R": []}},
                  "tts": {"support": {"x": [], "P": [], "R": []},
                          "gtscore": {"x": [], "R": []},
                          "dist": {"x": [], "P": [], "R": []}},
                  "paired": {"support": {"x": [1, 2], "P": [0.5, 0.6],
                                         "R": [0.4, 0.3]}}},
        "StringTie2": {"tss": {"support": {"x": [], "P": [], "R": []},
                               "gtscore": {"x": [], "R": []},
                               "dist": {"x": [], "P": [], "R": []}},
                       "tts": {"support": {"x": [], "P": [], "R": []},
                               "gtscore": {"x": [], "R": []},
                               "dist": {"x": [], "P": [], "R": []}},
                       "paired": {"support": {"x": [], "P": [], "R": []}}},
    }
    rows = S._data_to_rows(data)
    paired_modes = {r[0] for r in rows if r[1] == "paired"}
    assert paired_modes == {"FLAIR"}     # StringTie2 (empty) absent
    # FLAIR has paired_precision + paired_recall rows
    metrics = {r[4] for r in rows if r[0] == "FLAIR" and r[1] == "paired"}
    assert metrics == {"paired_precision", "paired_recall"}


# ════════════════════════════════════════════════════════════════════════
# End-to-end on the chr22-style fixture path is covered by the offline
# consistency check; here we pin behaviour via the BED loader on a canonical
# synthetic BED so load_calls_full preserves junctions for the JC-dedup.
# ════════════════════════════════════════════════════════════════════════

def test_load_calls_full_preserves_junctions(tmp_path):
    bed = tmp_path / "t.bed"
    # two-exon + strand isoform on chr1: blocks 200@0, 200@800 over 100..1000
    bed.write_text(
        "chr1\t100\t1000\tt1\t5\t+\t100\t1000\t0\t2\t200,200,\t0,800,\n")
    isos, has = S.load_calls_full(str(bed), None, None, True)
    assert has and len(isos) == 1
    iso = isos[0]
    assert iso["n_exons"] == 2
    assert iso["junctions"] == ((300, 900),)   # donor=100+200, acceptor=100+800
    assert iso["read_support"] == 5.0
    assert iso["tss"] == 100 and iso["tts"] == 1000
