#!/usr/bin/env python3
"""Tests for the THREE assembly-comparison figure fixes (2026-06-21):

  1. GTF isoform input via the SHARED signal_utils.load_isoforms_any (so
     StringTie2/IsoSeq/Bambu/LRAA/IsoQuant GTF tools are no longer silently
     dropped by stratified_pr.py / end_dist_boxplot.py), with read-support for
     GTF tools coming ONLY from --counts, and no-counts tools flagged for
     OMISSION (not silent drop).
  2. A --replot path on all three scripts that reads the script's own checkpoint
     TSV and regenerates identical figures with ZERO recompute.
  3. A RAW per-isoform-end distance checkpoint for end_dist_boxplot.py.

These pin:
  * GTF geometry == equivalent BED geometry (the byte-faithful invariant);
  * the BED code path is unchanged (geometry + inline support + canonical filter);
  * GTF tools with no counts get read_support=None / has_support=False;
  * the stratified_pr TSV <-> data round-trip is loss-free (replot redraws
    the same figure);
  * the end_dist raw checkpoint is long-form (label, end, distance) and
    round-trips through --replot.

conftest.py puts bin/ on sys.path (so `from evaluation.<...>` works).
"""
from __future__ import annotations

from pathlib import Path

import pytest

from evaluation.signal_utils import (
    load_isoforms_any, load_support_table, is_gtf_path)


# ── tiny synthetic fixtures (2 transcripts, + and - strands) ────────────────

_GTF = (
    'chr1\tt\ttranscript\t101\t1000\t.\t+\t.\tgene_id "g1"; transcript_id "iq1";\n'
    'chr1\tt\texon\t101\t300\t.\t+\t.\tgene_id "g1"; transcript_id "iq1";\n'
    'chr1\tt\texon\t801\t1000\t.\t+\t.\tgene_id "g1"; transcript_id "iq1";\n'
    'chr1\tt\ttranscript\t2001\t3000\t.\t-\t.\tgene_id "g2"; transcript_id "iq2";\n'
    'chr1\tt\texon\t2001\t2200\t.\t-\t.\tgene_id "g2"; transcript_id "iq2";\n'
    'chr1\tt\texon\t2801\t3000\t.\t-\t.\tgene_id "g2"; transcript_id "iq2";\n'
)

# BED12 equivalent of the GTF (same genomic spans, exon blocks, strands).
# g1: + 100..1000 blocks 200@0,200@700 ; g2: - 2000..3000 blocks 200@0,200@800
_BED = (
    "chr1\t100\t1000\tiq1\t5\t+\t100\t1000\t0\t2\t200,200,\t0,700,\n"
    "chr1\t2000\t3000\tiq2\t7\t-\t2000\t3000\t0\t2\t200,200,\t0,800,\n"
)


@pytest.fixture
def gtf_path(tmp_path):
    p = tmp_path / "tool.gtf"
    p.write_text(_GTF)
    return str(p)


@pytest.fixture
def bed_path(tmp_path):
    p = tmp_path / "tool.bed"
    p.write_text(_BED)
    return str(p)


@pytest.fixture
def counts_path(tmp_path):
    p = tmp_path / "counts.tsv"
    p.write_text("iq1\t11\niq2\t13\n")
    return str(p)


# ════════════════════════════════════════════════════════════════════════
# FIX 1 — shared GTF loader
# ════════════════════════════════════════════════════════════════════════

def test_is_gtf_path():
    assert is_gtf_path("x.gtf") and is_gtf_path("x.GFF3") and is_gtf_path("x.gff")
    assert not is_gtf_path("x.bed") and not is_gtf_path("x.txt")


def test_gtf_equals_bed_geometry(gtf_path, bed_path, counts_path):
    """A GTF tool and its equivalent BED tool yield IDENTICAL per-isoform
    geometry (chrom/strand/start/end/tss/tts/junctions)."""
    counts = load_support_table(counts_path)
    g = {i["name"]: i for i in load_isoforms_any(gtf_path, counts=counts)}
    b = {i["name"]: i for i in load_isoforms_any(bed_path)}
    assert set(g) == set(b) == {"iq1", "iq2"}
    for name in ("iq1", "iq2"):
        for k in ("chrom", "strand", "start", "end", "tss", "tts", "junctions"):
            assert g[name][k] == b[name][k], f"{name}.{k}: {g[name][k]} != {b[name][k]}"


def test_strand_aware_tss_tts(gtf_path):
    g = {i["name"]: i for i in load_isoforms_any(gtf_path)}
    # + strand: tss=start, tts=end
    assert g["iq1"]["tss"] == g["iq1"]["start"]
    assert g["iq1"]["tts"] == g["iq1"]["end"]
    # - strand: tss=end, tts=start
    assert g["iq2"]["tss"] == g["iq2"]["end"]
    assert g["iq2"]["tts"] == g["iq2"]["start"]


def test_gtf_support_only_from_counts(gtf_path, counts_path):
    """GTF read-support comes ONLY from --counts; no counts -> has_support False."""
    nocounts = load_isoforms_any(gtf_path, counts=None)
    assert all(i["read_support"] is None and i["has_support"] is False
               for i in nocounts)
    counts = load_support_table(counts_path)
    withcounts = {i["name"]: i for i in load_isoforms_any(gtf_path, counts=counts)}
    assert withcounts["iq1"]["read_support"] == 11.0 and withcounts["iq1"]["has_support"]
    assert withcounts["iq2"]["read_support"] == 13.0 and withcounts["iq2"]["has_support"]


def test_bed_support_inline(bed_path):
    """BED tools keep inline col5/col16 support (always has_support)."""
    b = {i["name"]: i for i in load_isoforms_any(bed_path)}
    assert b["iq1"]["read_support"] == 5.0 and b["iq1"]["has_support"]
    assert b["iq2"]["read_support"] == 7.0 and b["iq2"]["has_support"]


def test_bed_path_unchanged_canonical_filter(tmp_path):
    """The BED path drops non-canonical contigs exactly like before."""
    p = tmp_path / "c.bed"
    p.write_text(
        "chr1\t100\t1000\tt1\t5\t+\t100\t1000\t0\t1\t900,\t0,\n"
        "chrM\t10\t100\ttM\t3\t+\t10\t100\t0\t1\t90,\t0,\n"
        "GL000.1\t10\t100\ttG\t3\t+\t10\t100\t0\t1\t90,\t0,\n")
    isos = load_isoforms_any(str(p), canonical_only=True)
    assert {i["name"] for i in isos} == {"t1"}
    # with --all-contigs equivalent (canonical_only=False) everything stays
    isos_all = load_isoforms_any(str(p), canonical_only=False)
    assert {i["name"] for i in isos_all} == {"t1", "tM", "tG"}


def test_region_restriction(gtf_path, bed_path):
    """A region tuple list restricts both GTF and BED identically."""
    reg = [("chr1", 0, 1500)]   # only g1 (100..1000); g2 at 2000..3000 excluded
    g = load_isoforms_any(gtf_path, regions=reg)
    b = load_isoforms_any(bed_path, regions=reg)
    assert {i["name"] for i in g} == {"iq1"}
    assert {i["name"] for i in b} == {"iq1"}


# ════════════════════════════════════════════════════════════════════════
# FIX 1 — stratified_pr omission logic (GTF no-counts dropped from support)
# ════════════════════════════════════════════════════════════════════════

def test_stratified_load_calls_has_support_flag(gtf_path, bed_path, counts_path):
    from evaluation import stratified_pr as S
    # BED -> has_support True
    calls_b, has_b = S.load_calls(bed_path, None, None, True)
    assert has_b and len(calls_b) == 2
    # GTF + counts -> has_support True
    counts = S.load_counts_map(counts_path)
    calls_g, has_g = S.load_calls(gtf_path, None, counts, True)
    assert has_g and len(calls_g) == 2
    # GTF no counts -> has_support False (so caller omits from support figs)
    calls_n, has_n = S.load_calls(gtf_path, None, None, True)
    assert has_n is False and len(calls_n) == 2


# ════════════════════════════════════════════════════════════════════════
# FIX 2 — stratified_pr TSV <-> data round-trip (replot)
# ════════════════════════════════════════════════════════════════════════

def test_stratified_tsv_roundtrip(tmp_path):
    from evaluation import stratified_pr as S
    data = {
        "FLAIR": {
            "tss": {"support": {"x": [1, 2], "P": [0.85, 0.9], "R": [0.5, 0.4]},
                    "gtscore": {"x": ["0-10"], "R": [0.6]},
                    "dist": {"x": ["0-50", "50-100"], "P": [0.8, None],
                             "R": [0.7, 0.3]}},
            "tts": {"support": {"x": [1, 2], "P": [0.7, 0.75], "R": [0.6, 0.5]},
                    "gtscore": {"x": ["0-10"], "R": [0.5]},
                    "dist": {"x": ["0-50"], "P": [0.6], "R": [0.5]}},
            "paired": {"support": {"x": [1, 2], "P": [0.49, 0.63],
                                   "R": [0.62, 0.55]}},
        },
        # GTF no-counts tool: empty support, populated dist
        "StringTie2": {
            "tss": {"support": {"x": [], "P": [], "R": []},
                    "gtscore": {"x": ["0-10"], "R": [0.2]},
                    "dist": {"x": ["0-50"], "P": [0.3], "R": [0.2]}},
            "tts": {"support": {"x": [], "P": [], "R": []},
                    "gtscore": {"x": ["0-10"], "R": [0.1]},
                    "dist": {"x": ["0-50"], "P": [0.2], "R": [0.1]}},
            "paired": {"support": {"x": [], "P": [], "R": []}},
        },
    }
    tsv = tmp_path / "stratified_pr.tsv"
    S._write_tsv(tsv, S._data_to_rows(data))
    back = S._rows_to_data(tsv)

    # FLAIR full series survive (4-dp rounding)
    assert back["FLAIR"]["tss"]["support"]["x"] == [1, 2]
    assert back["FLAIR"]["tss"]["support"]["P"] == [0.85, 0.9]
    assert back["FLAIR"]["tss"]["dist"]["P"] == [0.8, None]   # None preserved
    # FLAIR paired support round-trips (int KS x, paired P + R)
    assert back["FLAIR"]["paired"]["support"]["x"] == [1, 2]
    assert back["FLAIR"]["paired"]["support"]["P"] == [0.49, 0.63]
    assert back["FLAIR"]["paired"]["support"]["R"] == [0.62, 0.55]
    # StringTie2: support empty (omitted), dist present, paired empty (omitted)
    assert back["StringTie2"]["tss"]["support"]["x"] == []
    assert back["StringTie2"]["tss"]["dist"]["x"] == ["0-50"]
    assert back["StringTie2"]["paired"]["support"]["x"] == []


def test_stratified_replot_no_compute(tmp_path):
    """--replot reads ONLY the checkpoint (a bogus --beds path would crash the
    compute path; replot must not touch it)."""
    from evaluation import stratified_pr as S
    out = tmp_path
    data = {"M": {w: {"support": {"x": [1], "P": [0.5], "R": [0.5]},
                      "gtscore": {"x": [], "R": []},
                      "dist": {"x": ["0-50"], "P": [0.5], "R": [0.5]}}
                  for w in ("tss", "tts")}}
    S._write_tsv(out / "stratified_pr.tsv", S._data_to_rows(data))
    import sys
    argv = sys.argv
    sys.argv = ["stratified_pr.py", "--outdir", str(out), "--replot"]
    try:
        S.main()
    finally:
        sys.argv = argv
    assert (out / "stratified_precision_by_reads.png").exists()


# ════════════════════════════════════════════════════════════════════════
# FIX 3 — end_dist raw checkpoint round-trip
# ════════════════════════════════════════════════════════════════════════

def test_end_dist_raw_checkpoint_roundtrip(tmp_path):
    from evaluation import end_dist_boxplot as E
    labels = ["FLAIR", "IsoQuant"]
    tss = {"FLAIR": [0.0, 5.0, 1000.0], "IsoQuant": [0.0, 0.0]}
    tts = {"FLAIR": [0.0, 2000.0], "IsoQuant": [3.0]}
    raw = tmp_path / E._RAW_NAME
    E._write_raw_checkpoint(raw, labels, 1, tss, tts)
    got_labels, got_ms, got_tss, got_tts = E._read_raw_checkpoint(raw)
    assert got_labels == labels
    assert got_ms == 1
    assert got_tss["FLAIR"] == [0.0, 5.0, 1000.0]
    assert got_tts["FLAIR"] == [0.0, 2000.0]
    assert got_tts["IsoQuant"] == [3.0]


def test_end_dist_raw_is_long_form(tmp_path):
    """The raw checkpoint is long-form: one row per isoform-end (label,end,dist)."""
    from evaluation import end_dist_boxplot as E
    raw = tmp_path / E._RAW_NAME
    E._write_raw_checkpoint(raw, ["A"], 1, {"A": [1.0, 2.0]}, {"A": [3.0]})
    body = [ln for ln in raw.read_text().splitlines()
            if ln and not ln.startswith("#") and not ln.startswith("label\t")]
    # 2 TSS rows + 1 TTS row = 3 raw rows
    assert len(body) == 3
    assert body[0].split("\t")[:2] == ["A", "5prime"]
    assert body[-1].split("\t")[:2] == ["A", "3prime"]


def test_end_dist_collect_has_support(gtf_path, bed_path, counts_path):
    from evaluation import end_dist_boxplot as E
    counts = load_support_table(counts_path)
    cage_idx, drna_idx = {}, {}   # has_support depends on read counts, not peaks
    # GTF with counts -> has_support True
    _, _, has_g = E.collect_distances(gtf_path, counts, cage_idx, drna_idx, [], 1, True)
    assert has_g
    # GTF no counts -> has_support False
    _, _, has_n = E.collect_distances(gtf_path, None, cage_idx, drna_idx, [], 1, True)
    assert has_n is False
    # BED -> has_support True
    _, _, has_b = E.collect_distances(bed_path, None, cage_idx, drna_idx, [], 1, True)
    assert has_b


def test_end_dist_nearest_peak_distance():
    """Fig D distance is to the ORTHOGONAL peak INTERVAL: 0 inside, else to edge."""
    from evaluation import end_dist_boxplot as E
    # one (chrom,strand) bucket with two peaks: [100,120) and [500,560)
    idx = E._build_peak_index({("chr1", "+"): [(100, 120), (500, 560)]})
    key = ("chr1", "+")
    assert E._nearest_peak_dist(110, key, idx) == 0      # inside first peak
    assert E._nearest_peak_dist(119, key, idx) == 0      # last base inside [100,120)
    assert E._nearest_peak_dist(90, key, idx) == 10      # 10 bp left of lo=100
    assert E._nearest_peak_dist(130, key, idx) == 11     # 130-(120-1)=11 right of first
    assert E._nearest_peak_dist(540, key, idx) == 0      # inside second peak
    assert E._nearest_peak_dist(300, key, idx) == 181    # 300-(120-1), nearer than 500-300
    assert E._nearest_peak_dist(110, ("chr1", "-"), idx) is None   # no peaks on - strand
    assert E._dist_to_interval(50, 100, 120) == 50
    assert E._dist_to_interval(100, 100, 120) == 0


# ════════════════════════════════════════════════════════════════════════
# FIX 2 — sqanti (Fig C) replot loader
# ════════════════════════════════════════════════════════════════════════

def test_sqanti_bar_tsv_roundtrip(tmp_path):
    from evaluation import sqanti_category_precise_bar as Q
    tsv = tmp_path / "figC_TED.tsv"
    tsv.write_text(
        "label\tcategory\tn_isoforms\tn_paired_precise\tprecise_fraction\n"
        "TED\tFSM\t100\t85\t0.850000\n"
        "TED\tISM\t40\t20\t0.500000\n")
    label, tot, pre = Q._read_bar_tsv(tsv)
    assert label == "TED"
    assert tot["FSM"] == 100 and pre["FSM"] == 85
    assert tot["ISM"] == 40 and pre["ISM"] == 20
