#!/usr/bin/env python3
"""Unit tests for STAGE 1 (two-pass TSS routing) + STAGE 2 (S_chain / S_junction
structural rerankers) in ted.py (flair-ted-52826).

Guarantees:

STAGE 1 — two-pass routing (--ted_two_pass, default OFF):
  * OFF (default): _ted_validate_cluster is byte-identical to the legacy
    single-pass blend (pass_label='', no re-blend), for both TSS and TTS.
  * NEAR  (annot_dist <= novel_dist): keeps the current blend INCLUDING S_a at
    the near threshold; pass_label='near', applied_threshold=near θ.
  * NOVEL (annot_dist >  novel_dist, or no annotation): re-blends reads-only
    WITHOUT S_a using the novel weights {S_d,S_m,S_shape} and the novel θ;
    pass_label='novel', applied_threshold=novel θ.
  * Routing is TSS-only — a TTS end never re-blends regardless of two_pass.
  * The geometry-lock-opening property: a novel TSS that FAILS the near 4-term
    gate can PASS on the reads-only novel blend at the lower novel θ.

STAGE 2 — structural rerankers (NOVEL pass, TSS-only):
  * S_chain (s_chain_score): negated truncation feature derived from the
    dominant same-chain upstream 5'-PEAK (the cT_updist WORKHORSE + peakratio +
    asym; NOT the weak per-read cT_upfrac). High (->1) when no dominant longer
    same-chain molecule peaks upstream, low (->0) when one does. A single stray
    far-upstream read does NOT move it (no peak); a far upstream PEAK does.
  * S_junction (s_junction_score): fj_novel_frac; high when the candidate's
    first junction is NOT annotated; 0 when it is annotated / too few reads.
  * Both enter the NOVEL blend additively via chain_reranker_score; OFF (weight
    0 / not computed) leaves the novel reality unchanged.
  * parse_gtf_introns builds the annotated-intron set S_junction tests against.

Loads ted.py standalone (numpy/scipy/sklearn at call time); skips if the clone
is absent. An end-to-end section (real FLAIR on the tiny fixture) asserts the
two-pass route is taken on a NOVEL candidate and the reranker columns populate.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

_TED_PY = Path(
    "/private/home/hdheath/projects/flair4_ted/flair-ted-52826/src/flair/ted.py")

pytestmark = pytest.mark.skipif(
    not _TED_PY.exists(), reason=f"TED clone not found: {_TED_PY}")


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ted = _load(_TED_PY, "ted_under_test_two_pass")


class _Scorer:
    """Minimal TedScorer stand-in: annotation only, no sequence model."""
    _tss_model = None
    _tts_model = None
    _flair_extractor_tss = None
    _flair_extractor_tts = None

    def __init__(self, tss_annot=(), tts_annot=()):
        self._annotated_ends = {
            'chrT': {'tss': sorted(tss_annot), 'tts': sorted(tts_annot)}}


def _validate(positions, scorer, **kw):
    """_ted_validate_cluster with a saturated-depth (S_d=1) logsat normalizer so
    only the routed weights/terms move the reality."""
    norm = ted.DepthNormalizer([], method='logsat', logsat_constant=1.0)
    base = dict(
        end_type='tss', chrom='chrT', strand='+', scorer=scorer,
        depth_norm=norm,
        w_depth=0.24, w_model=0.19, w_annot=0.17, w_shape=0.40, threshold=0.77,
        two_pass=True, two_pass_novel_dist=500,
        w_depth_novel=0.02, w_model_novel=0.11, w_shape_novel=0.87,
        threshold_novel=0.40)
    base.update(kw)
    return ted._ted_validate_cluster(positions, **base)


# ════════════════════════════════════════════════════════════════════════
# STAGE 1 — two-pass routing
# ════════════════════════════════════════════════════════════════════════
def test_two_pass_off_is_legacy_single_pass():
    """two_pass=False (default) -> pass_label='', no re-blend; reality is the
    legacy 4-term blend; identical whether or not a near anno exists."""
    sc = _Scorer(tss_annot=[1000])
    _, real, scores = _validate([5000] * 10, sc, two_pass=False, shape_score=0.9)
    assert scores['pass_label'] == ''
    # annot_dist 4000 > annot_max_dist (500 default) -> annotation_score 0.
    # legacy 4-term = w_d*1 + w_m*0 + w_a*0 + w_shape*0.9
    assert real == pytest.approx(0.24 * 1 + 0.40 * 0.9)
    # applied_threshold echoes the (unchanged) near threshold.
    assert scores['applied_threshold'] == pytest.approx(0.77)


def test_two_pass_near_keeps_annotation_and_near_theta():
    """NEAR end (at the annotated TSS) keeps S_a and the near θ; pass='near'."""
    sc = _Scorer(tss_annot=[1000])
    is_real, real, scores = _validate([1000] * 10, sc, shape_score=0.9)
    assert scores['pass_label'] == 'near'
    assert scores['applied_threshold'] == pytest.approx(0.77)
    # near = w_d*1 + w_a*exp(0)=1 + w_shape*0.9  (annot_dist 0 -> S_a=1)
    assert real == pytest.approx(0.24 * 1 + 0.17 * 1 + 0.40 * 0.9)
    assert is_real is True


def test_two_pass_novel_drops_annotation_uses_novel_blend():
    """NOVEL end (>novel_dist) re-blends reads-only WITHOUT S_a, novel θ."""
    sc = _Scorer(tss_annot=[1000])
    is_real, real, scores = _validate([5000] * 10, sc, shape_score=0.9)
    assert scores['pass_label'] == 'novel'
    assert scores['applied_threshold'] == pytest.approx(0.40)
    assert scores['annot_dist_bp'] == 4000
    # novel = w_d_novel*1 + w_m_novel*0 + w_shape_novel*0.9 (NO annotation term)
    assert real == pytest.approx(0.02 * 1 + 0.11 * 0 + 0.87 * 0.9)
    assert is_real is True   # 0.803 >= 0.40


def test_two_pass_no_annotation_routes_novel():
    """No annotation on the chrom (annot_dist=-1) -> NOVEL route."""
    sc = _Scorer(tss_annot=[])
    _, _, scores = _validate([5000] * 10, sc, shape_score=0.5)
    assert scores['pass_label'] == 'novel'
    assert scores['annot_dist_bp'] == -1


def test_two_pass_boundary_distance_is_near():
    """annot_dist exactly == novel_dist routes NEAR (<= is near)."""
    sc = _Scorer(tss_annot=[1000])
    _, _, scores = _validate([1500] * 10, sc, two_pass_novel_dist=500,
                             shape_score=0.5)
    assert scores['annot_dist_bp'] == 500
    assert scores['pass_label'] == 'near'


def test_two_pass_tts_never_reblends():
    """Routing is TSS-only: a TTS end keeps the legacy blend even with two_pass
    on (pass_label stays '')."""
    sc = _Scorer(tts_annot=[1000])
    norm = ted.DepthNormalizer([], method='logsat', logsat_constant=1.0)
    _, real, scores = ted._ted_validate_cluster(
        [9000] * 10, end_type='tts', chrom='chrT', strand='+', scorer=sc,
        depth_norm=norm, w_depth=0.67, w_model=0.19, w_annot=0.15,
        threshold=0.84, two_pass=True, two_pass_novel_dist=500,
        w_depth_novel=0.02, w_model_novel=0.11, w_shape_novel=0.87,
        threshold_novel=0.40)
    assert scores['pass_label'] == ''
    # legacy TTS = w_d*1 (model 0, annot 0 at dist 8000) -> 0.67
    assert real == pytest.approx(0.67)


def test_two_pass_opens_geometry_lock():
    """A novel TSS that FAILS the near 4-term gate (θ=0.77) PASSES the reads-only
    novel blend at θ=0.40 — the lever's headline property."""
    sc = _Scorer(tss_annot=[1000])
    # shape 0.45: near blend would be 0.24 + 0.40*0.45 = 0.42 (< 0.77 -> fail);
    # novel blend = 0.02 + 0.87*0.45 = 0.4115 (>= 0.40 -> pass).
    is_real, real, scores = _validate([5000] * 10, sc, shape_score=0.45)
    assert scores['pass_label'] == 'novel'
    assert real == pytest.approx(0.02 + 0.87 * 0.45)
    assert is_real is True
    # the SAME candidate under the near gate fails (sanity on the lock).
    near_blend = 0.24 + 0.40 * 0.45
    assert near_blend < 0.77


def test_two_pass_novel_weights_overridable():
    """Novel weights default to the fitted values but are taken as given."""
    sc = _Scorer(tss_annot=[])
    _, real, _ = _validate([5000] * 10, sc, shape_score=0.8,
                           w_depth_novel=0.0, w_model_novel=0.0,
                           w_shape_novel=1.0)
    assert real == pytest.approx(0.8)   # pure shape


# ════════════════════════════════════════════════════════════════════════
# STAGE 2 — S_chain reranker
# ════════════════════════════════════════════════════════════════════════
def _chain_read(start, end, introns=()):
    """(start, end, introns) tuple as collected in-flow from ReadRec.juncs."""
    return (start, end, tuple(introns))


def test_s_chain_real_no_upstream_support():
    """All same-chain reads start AT the candidate -> no upstream support ->
    cT_upfrac 0 -> S_chain 1 (real 5' boundary). + strand."""
    pos = 10000
    # spliced reads sharing intron (10500, 11000); 5' ends clustered at pos.
    reads = [_chain_read(pos + d, pos + 2500, [(10500, 11000)])
             for d in range(0, 40, 2)]
    s = ted.s_chain_score(reads, pos, '+', T=150)
    assert s == pytest.approx(1.0)


def test_s_chain_truncation_upstream_support():
    """A dominant same-chain PEAK sits ~1 kb UPSTREAM -> high cT_updist with a
    peakratio >= the candidate's own peak -> low S_chain (truncation artifact).
    + strand. This is the WORKHORSE peak signal, not per-read presence."""
    pos = 10000
    introns = [(10500, 11000)]
    reads = []
    # 5 candidate-start reads ...
    reads += [_chain_read(pos + d, pos + 2500, introns) for d in range(0, 10, 2)]
    # ... and 20 less-truncated reads of the SAME chain peaking ~1kb upstream.
    reads += [_chain_read(pos - 1000 + d, pos + 2500, introns)
              for d in range(0, 40, 2)]
    s = ted.s_chain_score(reads, pos, '+', T=150)
    assert s < 0.5   # a dominant upstream PEAK -> truncation-like


def test_s_chain_stray_upstream_read_does_not_move_score():
    """SPEC: a single stray far-upstream read (no PEAK) must NOT lower S_chain.
    A deep candidate 5' peak with one lone read 2 kb upstream stays real-like
    (high), because the workhorse is the pileup PEAK (peakratio ~ 1/cand_peak),
    not per-read upstream presence (which the weak cT_upfrac would have keyed
    on, dragging S_chain down)."""
    pos = 10000
    introns = [(10500, 11000)]
    # 20-read candidate 5' peak AT pos, plus ONE stray read 2 kb upstream.
    reads = [_chain_read(pos + d, pos + 2500, introns) for d in range(0, 40, 2)]
    reads += [_chain_read(pos - 2000, pos + 2500, introns)]
    s = ted.s_chain_score(reads, pos, '+', T=150)
    assert s > 0.75   # stray read forms no dominant peak -> stays real-like
    # contrast: the SAME far displacement as a real CLUSTER (20 reads) tanks it.
    reads_peak = [_chain_read(pos + d, pos + 2500, introns) for d in range(0, 10, 2)]
    reads_peak += [_chain_read(pos - 2000 + d, pos + 2500, introns)
                   for d in range(0, 40, 2)]
    s_peak = ted.s_chain_score(reads_peak, pos, '+', T=150)
    assert s_peak < s   # the cluster (peak) lowers S_chain below the stray case
    assert s_peak < 0.5


def test_s_chain_peak_helpers_directional():
    """_chain_peak_features + _s_chain_from_peak: a FAR dominant upstream peak
    (high updist, peakratio>=1, asym->+1) yields a small S_chain; no upstream
    peak (updist 0) yields S_chain 1; asym NaN (no mass beyond T) -> S_chain 1."""
    # far dominant upstream peak
    assert ted._s_chain_from_peak(0.76, 1.0, 1.0) < 0.5
    # an even farther / bigger cluster is even lower
    assert ted._s_chain_from_peak(0.90, 8.0, 1.0) < ted._s_chain_from_peak(0.76, 1.0, 1.0)
    # no upstream peak -> real
    assert ted._s_chain_from_peak(0.0, 0.0, -1.0) == pytest.approx(1.0)
    # stray-like: tiny peakratio gates the displacement away
    assert ted._s_chain_from_peak(0.96, 0.2, 1.0) > 0.75
    # asym NaN guard
    assert ted._s_chain_from_peak(0.9, 1.0, float('nan')) == pytest.approx(1.0)
    # updist 0 with mass present still real (no peak displacement to penalize)
    assert ted._s_chain_from_peak(0.0, 1.0, 1.0) == pytest.approx(1.0)


def test_s_chain_minus_strand_symmetry():
    """- strand: 5' end = read end; upstream = larger genomic coord. A dominant
    same-chain cluster ending DOWNSTREAM-in-genome (upstream in transcription)
    -> low S_chain."""
    pos = 10000   # 5' end (genomic) on - strand
    introns = [(8000, 9000)]
    reads = []
    reads += [_chain_read(pos - 2500, pos - d, introns) for d in range(0, 10, 2)]
    reads += [_chain_read(pos - 2500, pos + 1000 - d, introns)
              for d in range(0, 40, 2)]
    s = ted.s_chain_score(reads, pos, '-', T=150)
    assert s < 0.5


def test_s_chain_too_few_same_chain_neutral():
    """< MIN_READS same-chain reads -> neutral 0.5 (no evidence)."""
    pos = 10000
    reads = [_chain_read(pos, pos + 2500, [(10500, 11000)])]  # 1 read
    assert ted.s_chain_score(reads, pos, '+', T=150) == 0.5
    assert ted.s_chain_score([], pos, '+') == 0.5


def test_s_chain_different_chain_excluded():
    """Reads with a DIFFERENT downstream chain are not pooled, so their upstream
    starts do not count as truncation (alt-promoter, not a 5' fragment)."""
    pos = 10000
    cand = [_chain_read(pos + d, pos + 2500, [(10500, 11000)])
            for d in range(0, 40, 2)]
    # upstream reads but with a DIFFERENT first intron -> excluded from the pool.
    alt = [_chain_read(pos - 1000 + d, pos + 2500, [(10200, 10400)])
           for d in range(0, 40, 2)]
    s = ted.s_chain_score(cand + alt, pos, '+', T=150)
    assert s == pytest.approx(1.0)   # only the candidate chain counts


# ════════════════════════════════════════════════════════════════════════
# STAGE 2 — S_junction reranker + parse_gtf_introns
# ════════════════════════════════════════════════════════════════════════
def test_s_junction_novel_first_junction():
    """Reads whose first junction is NOT in the annotated set -> fj_novel_frac
    1 (real alt-promoter)."""
    pos = 10000
    annot = {('chrT', '+'): {(11000 // ted._JUNC_TOL, 12000 // ted._JUNC_TOL)}}
    # candidate reads' first junction is (10500, 11000) — NOT annotated.
    reads = [_chain_read(pos, pos + 2500, [(10500, 11000)]) for _ in range(6)]
    s = ted.s_junction_score(reads, pos, '+', 'chrT', annot)
    assert s == pytest.approx(1.0)


def test_s_junction_annotated_first_junction():
    """First junction IS annotated -> fj_novel_frac 0."""
    pos = 10000
    annot = {('chrT', '+'): {(10500 // ted._JUNC_TOL, 11000 // ted._JUNC_TOL)}}
    reads = [_chain_read(pos, pos + 2500, [(10500, 11000)]) for _ in range(6)]
    s = ted.s_junction_score(reads, pos, '+', 'chrT', annot)
    assert s == pytest.approx(0.0)


def test_s_junction_too_few_reads_zero():
    """< _JUNC_MIN_READS starting reads -> 0.0."""
    pos = 10000
    annot = {('chrT', '+'): set()}
    reads = [_chain_read(pos, pos + 2500, [(10500, 11000)])]  # 1 read
    assert ted.s_junction_score(reads, pos, '+', 'chrT', annot) == 0.0


def test_s_junction_no_annotation_set_all_novel():
    """No annotation entry for (chrom, strand) -> every junction reads as novel
    (is_annot False) -> fj_novel_frac 1 for spliced starters."""
    pos = 10000
    reads = [_chain_read(pos, pos + 2500, [(10500, 11000)]) for _ in range(6)]
    s = ted.s_junction_score(reads, pos, '+', 'chrT', {})
    assert s == pytest.approx(1.0)


def test_parse_gtf_introns(tmp_path):
    """parse_gtf_introns yields TOL-binned (donor, acceptor) intron sets per
    (chrom, strand) from consecutive exons of each transcript."""
    gtf = tmp_path / "mini.gtf"
    gtf.write_text(
        'chrT\tx\texon\t101\t200\t.\t+\t.\ttranscript_id "t1";\n'
        'chrT\tx\texon\t301\t400\t.\t+\t.\ttranscript_id "t1";\n'
        'chrT\tx\texon\t501\t600\t.\t+\t.\ttranscript_id "t1";\n')
    introns = ted.parse_gtf_introns(str(gtf))
    s = introns[('chrT', '+')]
    TOL = ted._JUNC_TOL
    # exon1.end=200 -> donor 200 (0-based half-open: GTF end is inclusive 200,
    # so donor coord 200); acceptor = next exon start 301-1 = 300.
    assert (200 // TOL, 300 // TOL) in s
    assert (400 // TOL, 500 // TOL) in s


# ════════════════════════════════════════════════════════════════════════
# STAGE 2 — _ChainReadIndex window
# ════════════════════════════════════════════════════════════════════════
def test_chain_read_index_window():
    """_ChainReadIndex returns reads overlapping the wide ±FETCH_PAD window."""
    pos = 100000
    pad = ted._CHAIN_FETCH_PAD
    chains = [
        (pos - 200, pos + 2000, ((100500, 101000),)),   # overlaps
        (pos + pad + 100, pos + pad + 500, ()),          # starts past hi, no overlap
        (pos - pad - 500, pos - pad - 100, ()),          # ends before lo, no overlap
        (pos - 100, pos + 100, ((100050, 100080),)),     # tiny, overlaps
    ]
    idx = ted._ChainReadIndex(chains, '+')
    got = idx.reads_for(pos)
    starts = sorted(s for (s, _e, _i) in got)
    assert (pos - 200) in starts and (pos - 100) in starts
    assert (pos + pad + 100) not in starts and (pos - pad - 500) not in starts


def test_chain_read_index_caps_deep_locus_deterministically():
    """SCALE GUARD: on an ultra-deep locus reads_for must CAP the returned reads at
    _CHAIN_MAX_READS — so per-peak s_chain/s_junction work is O(cap), not O(local
    depth) (the chr15-class quadratic-in-depth trap) — and be DETERMINISTIC (seeded
    per position). Byte-identical (no subsampling) below the cap."""
    pos = 100000
    introns = ((pos + 500, pos + 600),)
    deep = [(pos - 100 + (k % 50), pos + 2000, introns)
            for k in range(ted._CHAIN_SEEN_CAP + 1000)]   # > both caps, all overlap
    idx = ted._ChainReadIndex(deep, '+')
    got1 = idx.reads_for(pos)
    got2 = idx.reads_for(pos)
    assert len(got1) == ted._CHAIN_MAX_READS              # capped, not len(deep)
    assert got1 == got2                                   # deterministic (same seed)
    # Below the cap -> every overlapper returned unchanged.
    few = [(pos - 100, pos + 2000, introns),
           (pos - 50, pos + 1500, introns),
           (pos + 10, pos + 900, introns)]
    assert len(ted._ChainReadIndex(few, '+').reads_for(pos)) == 3


# ════════════════════════════════════════════════════════════════════════
# STAGE 3 — Pass-2 emit-protection (filter_spliced_iso collapse fix)
# ════════════════════════════════════════════════════════════════════════
# The host's filter_spliced_iso folds a junction-subset isoform into a longer
# co-located one. With --ted_protect_alt_ends a junction-subset whose OUTER
# terminal end is BOTH (a) distinct (>= delta past the superset's end) AND (b)
# TED-approved must NOT be collapsed. These integration tests drive the REAL
# filter_spliced_iso with lightweight stubs to prove the protect map actually
# protects (and that OFF / not-approved / not-distinct still collapse).
def _load_host(name):
    import sys
    clone_src = "/private/home/hdheath/projects/flair4_ted/flair-ted-52826/src"
    if clone_src not in sys.path:
        sys.path.insert(0, clone_src)
    try:
        import flair.flair_transcriptome as ft
    except Exception as e:   # heavy deps (pysam/xgboost) -> skip cleanly
        pytest.skip(f"cannot import flair_transcriptome: {e}")
    return ft


class _Exon(tuple):
    """Minimal Exon(start, end[, name]) — filter_spliced_iso only reads
    .start/.end and asserts isinstance(exons[0], flair Exon). We import the real
    Exon for the isinstance check below."""


class _IsoStub:
    def __init__(self, juncs, exons, score):
        self.juncs = juncs
        self.exons = exons
        self.score = score


class _Annots:
    junc_to_gene = {}


def _alt_tss_scenario(ft, protect_on, approved, distinct):
    """Build a + strand junction-subset QUERY whose genomic-FIRST exon (the TSS)
    is `distinct` bp upstream of a longer SUPERSET, and run filter_spliced_iso.

    superset: exons [(900,1100),(2100,2300),(3000,3500)] juncs J1=(1100,2100),
              J2=(2300,3000)   (3 exons, 2 juncs)
    query   : exons [(qstart,1100),(2100,2300)] juncs [J1]   (2 exons, 1 junc)
              -> J1 is a proper subset of (J1,J2); the query's LAST exon end
              (2300) == superset internal exon end -> internal-overlap subset on
              the 3' side; the query's FIRST exon end (1100) == superset first
              exon end -> terminal-overlap on the 5' (TSS) side.
    distinct True  -> qstart = 900 - 150 (>= delta=100 upstream of 900) so the
                      TSS sticks out; False -> qstart = 950 (inside the superset).
    """
    from flair.isoform_data import Junc, Exon
    J1 = Junc(1100, 2100)
    J2 = Junc(2300, 3000)
    sup_exons = [Exon(900, 1100), Exon(2100, 2300), Exon(3000, 3500)]
    superset = _IsoStub((J1, J2), sup_exons, score=50)
    qstart = (900 - 150) if distinct else 950
    q_exons = [Exon(qstart, 1100), Exon(2100, 2300)]
    query = _IsoStub((J1,), q_exons, score=20)
    all_isoforms = {'SUP': superset, 'Q': query}
    # both juncs of the superset map to its name; J1 also maps to the query.
    junc_to_names = {J1: {'SUP', 'Q'}, J2: {'SUP'}}
    ted_protect_map = None
    if protect_on:
        # + strand: 'first' exon carries the TSS. delta=100 (peak_min_distance_tss).
        ted_protect_map = {'Q': {'first': (100, approved), 'last': (10, approved)}}
    is_not_subset, _ = ft.filter_spliced_iso(
        'nosubset', 3, query.juncs, query.exons, 'Q', query.score,
        _Annots(), junc_to_names, all_isoforms, {}, '+',
        ted_protect_map=ted_protect_map)
    return is_not_subset


def test_protect_off_collapses_alt_tss():
    """Protect OFF (map None): the junction-subset alt-TSS collapses (is_subset)."""
    ft = _load_host("ft_protect_off")
    assert _alt_tss_scenario(ft, protect_on=False, approved=True, distinct=True) is False


def test_protect_on_distinct_approved_survives():
    """Protect ON + distinct + TED-approved: the alt-TSS is NOT collapsed."""
    ft = _load_host("ft_protect_on")
    assert _alt_tss_scenario(ft, protect_on=True, approved=True, distinct=True) is True


def test_protect_on_distinct_not_approved_collapses():
    """Protect ON but the end is NOT TED-approved: still collapses (we only
    rescue ends TED vouches for, not RT-truncation fragments)."""
    ft = _load_host("ft_protect_napprove")
    assert _alt_tss_scenario(ft, protect_on=True, approved=False, distinct=True) is False


def test_protect_on_approved_not_distinct_collapses():
    """Protect ON + approved but the end is INSIDE the superset (not distinct):
    still collapses (a contained fragment is not an alt-TSS)."""
    ft = _load_host("ft_protect_ndist")
    assert _alt_tss_scenario(ft, protect_on=True, approved=True, distinct=False) is False


def test_ted_outer_end_is_distinct_predicate():
    """The directional distinctness predicate: first-exon distinct when the
    query start is >= delta upstream; last-exon when >= delta downstream."""
    ft = _load_host("ft_distinct_pred")
    # first: query_outer <= other_outer - delta
    assert ft._ted_outer_end_is_distinct(800, 1000, 100, True, 'first') is True
    assert ft._ted_outer_end_is_distinct(950, 1000, 100, True, 'first') is False
    # last: query_outer >= other_outer + delta
    assert ft._ted_outer_end_is_distinct(1200, 1000, 100, True, 'last') is True
    assert ft._ted_outer_end_is_distinct(1050, 1000, 100, True, 'last') is False
    # not approved -> never distinct (no protection)
    assert ft._ted_outer_end_is_distinct(800, 1000, 100, False, 'first') is False
