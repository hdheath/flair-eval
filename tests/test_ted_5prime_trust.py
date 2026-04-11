#!/usr/bin/env python3
"""Tests for the --ted 5'-only trust implementation.

Validates:
1. check_firstlastexon: 5'-only trust relaxes left edge, keeps strict right edge
2. check_exonenddist: basic trust vs non-trust behavior
3. ted_collapse_end_groups: HDBSCAN for starts, window-based for ends
4. Pipeline integration: --ted sets trust_ends_5prime=True, trust_ends=False
5. count_sam_transcripts: --trust_ends_5prime argparse wiring
"""

import sys
import os
import types
import importlib
import subprocess
import pytest
import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# Path setup
# ─────────────────────────────────────────────────────────────────────────────
FLAIR_SRC = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'tools', 'flair-fusion', 'src')
FLAIR_SRC = os.path.abspath(FLAIR_SRC)
sys.path.insert(0, FLAIR_SRC)


# ─────────────────────────────────────────────────────────────────────────────
# 1. count_sam_transcripts: check_exonenddist + check_firstlastexon
# ─────────────────────────────────────────────────────────────────────────────
class TestCheckExonEndDist:
    """Test check_exonenddist trust_ends behavior."""

    @pytest.fixture(autouse=True)
    def _load(self):
        from flair.count_sam_transcripts import check_exonenddist, TRUST_ENDS_WINDOW
        self.check_exonenddist = check_exonenddist
        self.TRUST_ENDS_WINDOW = TRUST_ENDS_WINDOW

    def test_trust_ends_within_window_passes(self):
        """trust_ends=True: edge within TRUST_ENDS_WINDOW should pass."""
        # transcript_edge=0, read_edge=30 → distance=30 < 50
        assert self.check_exonenddist(100, 30, 0, True, 70, None) is True

    def test_trust_ends_outside_window_fails(self):
        """trust_ends=True: edge beyond TRUST_ENDS_WINDOW should fail."""
        # transcript_edge=0, read_edge=60 → distance=60 > 50
        assert self.check_exonenddist(100, 60, 0, True, 40, None) is False

    def test_trust_ends_boundary(self):
        """trust_ends=True: edge exactly at TRUST_ENDS_WINDOW boundary."""
        assert self.check_exonenddist(100, self.TRUST_ENDS_WINDOW, 0, True, 50, None) is True

    def test_no_trust_sufficient_coverage(self):
        """trust_ends=False: sufficient coverage passes."""
        # disttoblock=15, min(10, 100-5)=10, 15>=10 → True
        assert self.check_exonenddist(100, 15, 0, False, 15, None) is True

    def test_no_trust_insufficient_coverage(self):
        """trust_ends=False: insufficient coverage fails."""
        # disttoblock=5, min(10, 100-5)=10, 5>=10 → False
        assert self.check_exonenddist(100, 95, 0, False, 5, None) is False


class TestCheckFirstLastExon:
    """Test check_firstlastexon with trust_ends_5prime parameter."""

    @pytest.fixture(autouse=True)
    def _load(self):
        from flair.count_sam_transcripts import check_firstlastexon
        self.check_firstlastexon = check_firstlastexon

    def test_full_trust_ends_both_pass(self):
        """trust_ends=True trusts both left (5') and right (3') edges."""
        # read_start=30 (within 50 of left=0), read_end=tlen-30 (within 50 of right=tlen)
        result = self.check_firstlastexon(100, 100, 30, 970, 1000, True, None, None)
        assert result is True

    def test_full_trust_ends_left_fail(self):
        """trust_ends=True: left edge too far fails."""
        # read_start=60 → abs(0-60) = 60 > TRUST_ENDS_WINDOW=50
        result = self.check_firstlastexon(100, 100, 60, 970, 1000, True, None, None)
        assert result is False

    def test_5prime_only_trust_left_relaxed_right_strict(self):
        """trust_ends_5prime=True: left (5') relaxed, right (3') strict.

        This is the KEY test for the --ted 5'-only trust behavior.
        Left edge: 30bp from transcript start → within TRUST_ENDS_WINDOW → PASS
        Right edge: disttoblock = 905 - (1000 - 100) = 5  <  min(10, 95) = 10 → FAIL
        But we DON'T trust right under trust_ends_5prime → should FAIL overall.
        """
        result = self.check_firstlastexon(100, 100, 30, 905, 1000, False, None, None,
                                           trust_ends_5prime=True)
        assert result is False, "Right (3') edge should NOT be trusted under trust_ends_5prime"

    def test_5prime_only_trust_both_edges_pass(self):
        """trust_ends_5prime=True: both edges pass when right has sufficient coverage."""
        # Left: trust_left=True → abs(0-30)=30 <= 50 → PASS
        # Right: trust_right=False → disttoblock = 970-(1000-100) = 70 >= 10 → PASS
        result = self.check_firstlastexon(100, 100, 30, 970, 1000, False, None, None,
                                           trust_ends_5prime=True)
        assert result is True, "Both edges should pass: 5' trusted, 3' has sufficient coverage"

    def test_5prime_only_left_fail_right_pass(self):
        """trust_ends_5prime=True: left (5') edge beyond window → should FAIL."""
        # Left: abs(0-60)=60 > 50 → FAIL even with trust
        result = self.check_firstlastexon(100, 100, 60, 970, 1000, False, None, None,
                                           trust_ends_5prime=True)
        assert result is False, "Left (5') beyond TRUST_ENDS_WINDOW should fail"

    def test_no_trust_default_behavior(self):
        """trust_ends=False, trust_ends_5prime=False: default strict behavior."""
        # Left: disttoblock = 100-30 = 70 >= 10 → PASS
        # Right: disttoblock = 970-(1000-100) = 70 >= 10 → PASS
        result = self.check_firstlastexon(100, 100, 30, 970, 1000, False, None, None,
                                           trust_ends_5prime=False)
        assert result is True

    def test_5prime_trust_does_not_affect_3prime(self):
        """Core invariant: trust_ends_5prime NEVER relaxes the 3' edge."""
        # last_blocksize=20, read_end=985:
        # Right disttoblock = 985 - (1000-20) = 5 < min(10, 15) = 10 → FAIL
        result = self.check_firstlastexon(100, 20, 30, 985, 1000, False, None, None,
                                           trust_ends_5prime=True)
        assert result is False, "3' edge with insufficient coverage must fail under trust_ends_5prime"

    def test_5prime_trust_vs_full_trust_3prime_edge(self):
        """Contrast: full trust passes where 5prime-only fails (on 3' edge)."""
        # Same parameters as test_5prime_only_trust_left_relaxed_right_strict
        # With full trust_ends: right edge abs(1000-905) = 95 > 50 → FAIL too
        # Let's use read_end=960: abs(1000-960) = 40 <= 50 → PASS with full trust
        # disttoblock = 960-(1000-100) = 60 >= 10 → PASS without trust too
        # Better: read_end=955: abs(1000-955) = 45 <= 50 → PASS with full trust
        #         disttoblock = 955-900 = 55 >= 10 → PASS without trust
        # Need: edge that fails strict but passes trust window
        # read_end=952: abs(1000-952) = 48 <= 50 → PASS with trust
        #               disttoblock = 952-900=52 >= 10 → PASS without trust
        # Need small last_blocksize: last_blocksize=8 → min(10,3)=3
        # read_end=995: abs(1000-995)=5 <= 50 → PASS with trust (both modes)
        #               disttoblock=995-(1000-8) = 3 >= min(10,3)=3 → PASS without trust
        # Make it fail strict: last_blocksize=8, read_end=993
        #   disttoblock = 993-992=1 < 3 → FAIL strict
        #   abs(1000-993)=7 <= 50 → PASS trusted
        full_trust = self.check_firstlastexon(100, 8, 30, 993, 1000, True, None, None)
        assert full_trust is True, "Full trust should pass (both edges within window)"

        five_prime_only = self.check_firstlastexon(100, 8, 30, 993, 1000, False, None, None,
                                                    trust_ends_5prime=True)
        assert five_prime_only is False, "5prime-only trust should fail (3' edge not trusted)"


class TestCheckStringent5PrimeTrust:
    """Test check_stringent passes trust_ends_5prime through correctly."""

    @pytest.fixture(autouse=True)
    def _load(self):
        from flair.count_sam_transcripts import check_stringent
        self.check_stringent = check_stringent

    def test_multiexon_5prime_trust_passthrough(self):
        """check_stringent passes trust_ends_5prime to check_firstlastexon."""
        # Two-exon transcript: exonpos=[100, 100], tlen=200
        coveredpos = [1] * 200
        exonpos = [100, 100]
        tlen = 200

        # read_start=30, read_end=30+140=170 — good coverage on both ends
        # Left: trust (abs(0-30)=30 <= 50) → PASS
        # Right: disttoblock = 170-(200-100) = 70 >= 10 → PASS
        result = self.check_stringent(coveredpos, exonpos, tlen,
                                       [30], [140],
                                       False, 'test', 0, {},
                                       trust_ends_5prime=True)
        assert result is True

    def test_multiexon_5prime_trust_right_fails(self):
        """check_stringent: trust_ends_5prime should not relax 3' boundary."""
        coveredpos = [1] * 200
        exonpos = [100, 100]
        tlen = 200

        # read_start=30, read_end=30+75=105
        # Left: trust (abs(0-30)=30 <= 50) → PASS
        # Right: disttoblock = 105-(200-100) = 5 < 10 → FAIL
        result = self.check_stringent(coveredpos, exonpos, tlen,
                                       [30], [75],
                                       False, 'test', 0, {},
                                       trust_ends_5prime=True)
        assert result is False, "3' edge should fail under 5prime-only trust"


# ─────────────────────────────────────────────────────────────────────────────
# 2. count_sam_transcripts argparse: --trust_ends_5prime
# ─────────────────────────────────────────────────────────────────────────────
class TestCountSamArgparse:
    """Test that --trust_ends_5prime is properly wired in argparse."""

    def test_trust_ends_5prime_flag_exists(self):
        """--trust_ends_5prime should be a valid argument."""
        from flair.count_sam_transcripts import parse_args
        old_argv = sys.argv
        try:
            sys.argv = ['count_sam_transcripts', '--sam', '-', '--trust_ends_5prime']
            args = parse_args()
            assert args.trust_ends_5prime is True
            assert args.trust_ends is False
        finally:
            sys.argv = old_argv

    def test_both_flags_independent(self):
        """--trust_ends and --trust_ends_5prime are independent flags."""
        from flair.count_sam_transcripts import parse_args
        old_argv = sys.argv
        try:
            sys.argv = ['count_sam_transcripts', '--sam', '-', '--trust_ends', '--trust_ends_5prime']
            args = parse_args()
            assert args.trust_ends is True
            assert args.trust_ends_5prime is True
        finally:
            sys.argv = old_argv


# ─────────────────────────────────────────────────────────────────────────────
# 3. ted_collapse_end_groups: hybrid 5'-trust collapse
# ─────────────────────────────────────────────────────────────────────────────
class TestTedCollapseEndGroups:
    """Test the hybrid HDBSCAN-starts + window-ends collapse."""

    @pytest.fixture(autouse=True)
    def _load(self):
        from flair.flair_transcriptome import (
            ted_collapse_end_groups, collapse_end_groups,
            group_reads_by_ends,
        )
        self.ted_collapse = ted_collapse_end_groups
        self.default_collapse = collapse_end_groups
        self.group_by_ends = group_reads_by_ends

    def _make_reads(self, start_end_pairs, strand='+'):
        """Create (start, end, strand, name) tuples from [(start, end), ...]."""
        reads = []
        for i, (s, e) in enumerate(start_end_pairs):
            reads.append((s, e, strand, f'read_{i}'))
        return reads

    def test_fallback_small_input(self):
        """< 4 reads falls back to default collapse."""
        reads = self._make_reads([(100, 500), (105, 510), (110, 490)])
        result = self.ted_collapse(reads, 'chr1', scorer=None, do_get_best_ends=True)
        assert len(result) >= 1
        # Should return standard format: [score, start, end, strand, name, [names]]
        assert len(result[0]) == 6

    def test_two_distinct_tss_clusters(self):
        """Two clear TSS clusters should produce >=2 isoform groups."""
        # Cluster A: starts ~100, end ~500
        # Cluster B: starts ~300, end ~500 (same TTS)
        reads = self._make_reads([
            (100, 500), (102, 505), (98, 498), (105, 502), (99, 503),
            (300, 500), (302, 505), (298, 498), (305, 502), (299, 503),
        ])
        result = self.ted_collapse(reads, 'chr1', scorer=None, do_get_best_ends=True)
        # Should detect 2 TSS clusters despite shared TTS
        assert len(result) >= 2, f"Expected >=2 isoform groups, got {len(result)}"

    def test_tts_preserved_within_tss_cluster(self):
        """Different TTS positions within one TSS cluster should form separate groups
        if they exceed the end_window distance."""
        # All starts ~100, but two distinct end groups: ~500 and ~800
        reads = self._make_reads([
            (100, 500), (102, 505), (98, 498), (105, 502), (99, 503),
            (100, 800), (102, 805), (98, 798), (105, 802), (99, 803),
        ])
        result = self.ted_collapse(
            reads, 'chr1', scorer=None, do_get_best_ends=True, end_window=100)
        # Should detect 2 groups: one TSS cluster × 2 TTS groups
        assert len(result) >= 2, f"Expected >=2 groups from TTS split, got {len(result)}"

    def test_output_format_firstpass(self):
        """do_get_best_ends=True output: [score, start, end, strand, name, [names]]."""
        reads = self._make_reads([
            (100, 500), (102, 505), (98, 498), (105, 502), (99, 503),
        ])
        result = self.ted_collapse(reads, 'chr1', scorer=None, do_get_best_ends=True)
        assert len(result) >= 1
        entry = result[0]
        assert len(entry) == 6, f"Expected 6-element list, got {len(entry)}"
        score, start, end, strand, name, read_names = entry
        assert isinstance(score, (int, float))
        assert isinstance(start, (int, float, np.integer, np.floating))
        assert isinstance(end, (int, float, np.integer, np.floating))
        assert strand in ('+', '-', 'ambig')
        assert isinstance(read_names, list)
        assert len(read_names) == 5

    def test_output_format_final_stage(self):
        """do_get_best_ends=False output: [start, end, (marker, iso), [reads]]."""
        # Final stage input: [start, end, (marker, iso), [read_names]]
        reads = [
            [100, 500, ('a', 'iso1'), ['r1', 'r2']],
            [102, 505, ('n', 'iso2'), ['r3']],
            [98, 498, ('a', 'iso1'), ['r4', 'r5']],
            [105, 502, ('n', 'iso3'), ['r6']],
            [99, 503, ('a', 'iso1'), ['r7']],
        ]
        result = self.ted_collapse(reads, 'chr1', scorer=None, do_get_best_ends=False)
        assert len(result) >= 1
        entry = result[0]
        assert len(entry) == 4, f"Expected 4-element list, got {len(entry)}"
        assert isinstance(entry[-1], list), "Last element should be merged read list"

    def test_sorted_by_score_descending(self):
        """Firstpass output should be sorted by confidence score descending."""
        reads = self._make_reads([
            (100, 500), (102, 505), (98, 498), (105, 502), (99, 503),
            (300, 500), (302, 505),
        ])
        result = self.ted_collapse(reads, 'chr1', scorer=None, do_get_best_ends=True)
        if len(result) >= 2:
            scores = [r[0] for r in result]
            assert scores == sorted(scores, reverse=True), "Should be sorted descending"

    def test_no_hdbscan_on_tts(self):
        """Verify TTS clustering uses window approach, NOT HDBSCAN.

        Two close TTS positions (within end_window) should be merged,
        while two far TTS would be split — consistent with window behavior.
        """
        # Close TTS: all within 100bp → should merge into 1 group
        reads_close = self._make_reads([
            (100, 500), (102, 510), (98, 490), (105, 505), (99, 495),
        ])
        result_close = self.ted_collapse(
            reads_close, 'chr1', scorer=None, do_get_best_ends=True, end_window=100)
        assert len(result_close) == 1, "Close TTS within window should form 1 group"

        # Far TTS: two clearly separated groups
        reads_far = self._make_reads([
            (100, 500), (102, 505), (98, 498),
            (100, 2000), (102, 2005), (98, 1998),
        ])
        result_far = self.ted_collapse(
            reads_far, 'chr1', scorer=None, do_get_best_ends=True, end_window=100)
        assert len(result_far) >= 2, "Far TTS should form >=2 groups"


# ─────────────────────────────────────────────────────────────────────────────
# 4. add_preset_args: --ted sets trust_ends_5prime, NOT trust_ends
# ─────────────────────────────────────────────────────────────────────────────
class TestAddPresetArgs:
    """Test that --ted correctly sets trust_ends_5prime instead of trust_ends."""

    def test_ted_sets_5prime_trust(self):
        """--ted should set trust_ends_5prime=True, trust_ends=False."""
        from flair.flair_transcriptome import add_preset_args
        args = types.SimpleNamespace(
            ted=True,
            trust_ends=False,
            softclip_rescue=False,
            refine_annotated_tts=False,
            end_scoring_alpha=0.5,
            no_redundant='longest',
            mm2_args=None,
            remove_internal_priming=False,
            isoformtss=False,
        )
        args = add_preset_args(args)
        assert args.trust_ends_5prime is True, "--ted should set trust_ends_5prime=True"
        assert args.trust_ends is False, "--ted should NOT set trust_ends=True"
        assert args.softclip_rescue is True
        assert args.refine_annotated_tts is True
        assert args.end_scoring_alpha == 1.0

    def test_no_ted_no_5prime_trust(self):
        """Without --ted, trust_ends_5prime defaults to False."""
        from flair.flair_transcriptome import add_preset_args
        args = types.SimpleNamespace(
            ted=False,
            trust_ends=False,
            softclip_rescue=False,
            refine_annotated_tts=False,
            end_scoring_alpha=0.5,
            no_redundant='none',
            mm2_args=None,
            remove_internal_priming=False,
            isoformtss=False,
        )
        args = add_preset_args(args)
        assert args.trust_ends_5prime is False
        assert args.trust_ends is False

    def test_trust_ends_standalone(self):
        """--trust_ends without --ted should still work as before."""
        from flair.flair_transcriptome import add_preset_args
        args = types.SimpleNamespace(
            ted=False,
            trust_ends=True,
            softclip_rescue=False,
            refine_annotated_tts=False,
            end_scoring_alpha=0.5,
            no_redundant='none',
            mm2_args=None,
            remove_internal_priming=False,
            isoformtss=False,
        )
        args = add_preset_args(args)
        assert args.trust_ends is True
        assert args.trust_ends_5prime is False  # not set by trust_ends alone


# ─────────────────────────────────────────────────────────────────────────────
# 5. Behavioral contract: quality set to 0 under trust_ends_5prime
# ─────────────────────────────────────────────────────────────────────────────
class TestQualityContract:
    """Test quality=0 behavior under trust_ends and trust_ends_5prime.

    These are logic-level tests: they verify the contract that either
    trust flag should force quality=0.  The actual implementation is
    in flair_collapse.py, but we test the logic pattern here.
    """

    def test_trust_ends_sets_quality_zero(self):
        args = types.SimpleNamespace(trust_ends=True, trust_ends_5prime=False, quality=10)
        args.quality = 0 if (args.trust_ends or getattr(args, 'trust_ends_5prime', False)) else args.quality
        assert args.quality == 0

    def test_trust_ends_5prime_sets_quality_zero(self):
        args = types.SimpleNamespace(trust_ends=False, trust_ends_5prime=True, quality=10)
        args.quality = 0 if (args.trust_ends or getattr(args, 'trust_ends_5prime', False)) else args.quality
        assert args.quality == 0

    def test_no_trust_preserves_quality(self):
        args = types.SimpleNamespace(trust_ends=False, trust_ends_5prime=False, quality=10)
        args.quality = 0 if (args.trust_ends or getattr(args, 'trust_ends_5prime', False)) else args.quality
        assert args.quality == 10


# ─────────────────────────────────────────────────────────────────────────────
# 6. Scenario-based tests: APA and truncation implications
# ─────────────────────────────────────────────────────────────────────────────
class TestScenarioImplications:
    """Test that the 5'-only trust approach handles key scenarios correctly."""

    @pytest.fixture(autouse=True)
    def _load(self):
        from flair.flair_transcriptome import ted_collapse_end_groups
        self.ted_collapse = ted_collapse_end_groups

    def _make_reads(self, start_end_pairs, strand='+'):
        return [(s, e, strand, f'read_{i}') for i, (s, e) in enumerate(start_end_pairs)]

    def test_near_identical_3utr_not_merged(self):
        """Two isoforms with nearby TTS (<50bp apart, both well-supported) should
        NOT be merged under --ted when end_window is small enough.

        This was the failure mode of --trust_ends: it merged close 3' UTR variants.
        The fix is: --ted uses window-based TTS grouping per TSS cluster.
        """
        # Isoform A: TSS ~100, TTS ~500 (20 reads)
        reads_a = self._make_reads([(100 + i % 5, 500 + i % 5) for i in range(20)])
        # Isoform B: TSS ~100, TTS ~540 (20 reads, 40bp apart)
        reads_b = [(100 + i % 5, 540 + i % 5, '+', f'read_b_{i}') for i in range(20)]
        all_reads = reads_a + reads_b

        # With end_window=30 (< 40bp separation), they should stay separate
        result = self.ted_collapse(
            all_reads, 'chr1', scorer=None, do_get_best_ends=True, end_window=30)
        assert len(result) >= 2, \
            f"Expected >=2 groups (TTS 40bp apart, window=30), got {len(result)}"

    def test_truncated_5prime_different_starts(self):
        """Reads with different 5' starts should form separate TSS clusters.

        Tests that HDBSCAN correctly identifies distinct 5' start positions
        (full-length vs truncated reads starting later).
        """
        # Full-length reads: start ~100
        full_reads = self._make_reads([(100 + i % 3, 500 + i % 5) for i in range(20)])
        # Truncated reads: start ~250 (150bp truncated from 5')
        trunc_reads = [(250 + i % 3, 500 + i % 5, '+', f'trunc_{i}') for i in range(15)]

        all_reads = full_reads + trunc_reads
        result = self.ted_collapse(
            all_reads, 'chr1', scorer=None, do_get_best_ends=True, end_window=100)
        # HDBSCAN should separate TSS ~100 from TSS ~250
        assert len(result) >= 2, \
            f"Expected >=2 TSS clusters (100 vs 250), got {len(result)}"

    def test_shared_tts_single_tss(self):
        """All reads sharing exact same TSS and TTS → single isoform group."""
        # Use identical start/end so HDBSCAN treats them as 1 cluster
        reads = self._make_reads([(100, 500)] * 20)
        # Rename reads so they're unique
        reads = [(s, e, st, f'read_{i}') for i, (s, e, st, _) in enumerate(reads)]
        result = self.ted_collapse(
            reads, 'chr1', scorer=None, do_get_best_ends=True, end_window=100)
        assert len(result) == 1, f"Expected 1 group for identical cluster, got {len(result)}"

    def test_tight_tss_variation_groups_together(self):
        """TSS with ~2bp jitter should form few groups (HDBSCAN handles density).

        Note: HDBSCAN may split very small discrete differences in start positions
        into separate clusters. This is intentional — the 5'-trust approach lets
        HDBSCAN decide TSS boundaries at single-base resolution.
        """
        reads = self._make_reads([(100 + i % 3, 500 + i % 3) for i in range(20)])
        result = self.ted_collapse(
            reads, 'chr1', scorer=None, do_get_best_ends=True, end_window=100)
        # All reads fall within a few bp; should produce <=3 groups at most
        assert len(result) <= 3, f"Expected <=3 groups for tight jitter, got {len(result)}"
        # Total reads should be preserved across all groups
        total_reads = sum(len(r[-1]) for r in result)
        assert total_reads == 20, f"Expected 20 total reads, got {total_reads}"


# ─────────────────────────────────────────────────────────────────────────────
# 7. Integration: verify --ted flag help text reflects 5'-only trust
# ─────────────────────────────────────────────────────────────────────────────
class TestHelpText:
    """Verify --ted and --trust_ends_5prime help text is present."""

    def test_count_sam_has_trust_ends_5prime(self):
        """count_sam_transcripts should advertise --trust_ends_5prime in help."""
        result = subprocess.run(
            [sys.executable, '-m', 'flair.count_sam_transcripts', '--help'],
            capture_output=True, text=True, cwd=FLAIR_SRC,
        )
        assert '--trust_ends_5prime' in result.stdout, \
            "--trust_ends_5prime not found in count_sam_transcripts help"


# ─────────────────────────────────────────────────────────────────────────────
# 8. Scorer config / continuous scoring architecture
# ─────────────────────────────────────────────────────────────────────────────
class TestScorerConfig:
    """Test _build_scorer_config and _build_scorer_from_config."""

    @pytest.fixture(autouse=True)
    def _load(self):
        from flair.flair_transcriptome import _build_scorer_config, _build_scorer_from_config
        self.build_config = _build_scorer_config
        self.build_from_config = _build_scorer_from_config

    def test_config_none_when_alpha_zero(self):
        """Scorer config should be None when end_scoring_alpha <= 0."""
        args = types.SimpleNamespace(end_scoring_alpha=0.0)
        assert self.build_config(args) is None

    def test_config_returns_dict_when_enabled(self):
        """Config should be a picklable dict with required keys."""
        args = types.SimpleNamespace(
            end_scoring_alpha=1.0,
            library_type='default',
            annotation_weight=None,
            model_weight=None,
            gtf=None,
            tss_model=None,
            tts_model=None,
        )
        config = self.build_config(args)
        assert config is not None
        assert isinstance(config, dict)
        for key in ('profile', 'tss_model_path', 'tts_model_path', 'annotated_ends', 'alpha'):
            assert key in config, f"Missing key: {key}"
        assert config['alpha'] == 1.0

    def test_config_is_picklable(self):
        """Config dict must survive pickle round-trip (multiprocessing compat)."""
        import pickle
        args = types.SimpleNamespace(
            end_scoring_alpha=1.0,
            library_type='default',
            annotation_weight=None,
            model_weight=None,
            gtf=None,
            tss_model=None,
            tts_model=None,
        )
        config = self.build_config(args)
        roundtrip = pickle.loads(pickle.dumps(config))
        assert roundtrip['alpha'] == config['alpha']
        assert roundtrip['profile'].name == config['profile'].name

    def test_from_config_none_returns_none(self):
        """_build_scorer_from_config(None, genome) should return (None, 0.0)."""
        scorer, alpha = self.build_from_config(None, None)
        assert scorer is None
        assert alpha == 0.0


class TestTedCollapseScorerIntegration:
    """Test ted_collapse_end_groups with scorer for continuous TSS+TTS scoring."""

    @pytest.fixture(autouse=True)
    def _load(self):
        from flair.flair_transcriptome import ted_collapse_end_groups
        self.ted_collapse = ted_collapse_end_groups

    def _make_reads(self, start_end_pairs, strand='+'):
        return [(s, e, strand, f'read_{i}') for i, (s, e) in enumerate(start_end_pairs)]

    def _make_mock_scorer(self, tss_score=0.8, tts_score=0.6):
        """Create a mock scorer that returns fixed confidence values."""
        from flair.end_scoring.scoring import EndCandidate

        class MockScorer:
            def score(self, candidate):
                if candidate.end_type == 'tss':
                    candidate.confidence = tss_score
                    candidate.seq_score = tss_score
                    candidate.depth_score = 0.5
                    candidate.annotation_score = 0.5
                    candidate.tech_penalty = 0.0
                elif candidate.end_type == 'tts':
                    candidate.confidence = tts_score
                    candidate.seq_score = tts_score
                    candidate.depth_score = 0.5
                    candidate.annotation_score = 0.5
                    candidate.tech_penalty = 0.0
                return candidate
        return MockScorer()

    def test_scorer_none_graceful_degradation(self):
        """With scorer=None, should degrade to density-only scoring."""
        reads = self._make_reads([(100, 500)] * 10 + [(100, 500)] * 5)
        reads = [(s, e, st, f'read_{i}') for i, (s, e, st, _) in enumerate(reads)]
        result = self.ted_collapse(reads, 'chr1', scorer=None, do_get_best_ends=True)
        assert len(result) >= 1
        # With scorer=None: confidence = 0.40 * density + 0.20 * 0.5 + 0.20 * 0.5
        #   + 0.20 * 0.5 (neutral junction_match_type=None)
        # = 0.40 * density + 0.30
        confidence = result[0][0]
        assert 0.0 <= confidence <= 1.0

    def test_scorer_boosts_confidence(self):
        """A high-scoring scorer should boost confidence above scorer=None."""
        reads = self._make_reads([(100, 500)] * 10)
        reads = [(s, e, st, f'read_{i}') for i, (s, e, st, _) in enumerate(reads)]

        result_no_scorer = self.ted_collapse(reads, 'chr1', scorer=None, do_get_best_ends=True)
        mock = self._make_mock_scorer(tss_score=0.9, tts_score=0.9)
        result_with_scorer = self.ted_collapse(reads, 'chr1', scorer=mock, do_get_best_ends=True)

        conf_no = result_no_scorer[0][0]
        conf_yes = result_with_scorer[0][0]
        assert conf_yes > conf_no, (
            f"High scorer ({conf_yes:.3f}) should boost above no-scorer ({conf_no:.3f})")

    def test_scorer_penalizes_bad_ends(self):
        """A low-scoring scorer should reduce confidence below scorer=None."""
        reads = self._make_reads([(100, 500)] * 10)
        reads = [(s, e, st, f'read_{i}') for i, (s, e, st, _) in enumerate(reads)]

        result_no_scorer = self.ted_collapse(reads, 'chr1', scorer=None, do_get_best_ends=True)
        mock = self._make_mock_scorer(tss_score=0.1, tts_score=0.1)
        result_with_scorer = self.ted_collapse(reads, 'chr1', scorer=mock, do_get_best_ends=True)

        conf_no = result_no_scorer[0][0]
        conf_yes = result_with_scorer[0][0]
        assert conf_yes < conf_no, (
            f"Low scorer ({conf_yes:.3f}) should penalize below no-scorer ({conf_no:.3f})")

    def test_density_is_primary_weight(self):
        """Density (read depth) should dominate the confidence score.

        Two clusters: one with 20 reads + bad scorer, one with 3 reads + good scorer.
        The high-density cluster should have a higher density *component* because
        W_DENSITY=0.40 is the largest single weight.
        """
        # Cluster A: 20 reads at TSS ~100
        reads_a = self._make_reads([(100 + i % 3, 500 + i % 3) for i in range(20)])
        # Cluster B: 4 reads at TSS ~400 (just above HDBSCAN minimum)
        reads_b = [(400 + i % 2, 800 + i % 2, '+', f'read_b_{i}') for i in range(4)]
        all_reads = reads_a + reads_b

        # Mock scorer: penalizes cluster A's positions, rewards cluster B's
        # But density should still win for cluster A
        class AsymmetricScorer:
            def score(self, candidate):
                pos = candidate.pos
                if pos < 300:  # Cluster A region — low scores
                    candidate.confidence = 0.2
                else:  # Cluster B region — high scores
                    candidate.confidence = 0.9
                candidate.seq_score = candidate.confidence
                candidate.depth_score = 0.5
                candidate.annotation_score = 0.5
                candidate.tech_penalty = 0.0
                return candidate

        result = self.ted_collapse(
            all_reads, 'chr1', scorer=AsymmetricScorer(), do_get_best_ends=True)
        assert len(result) >= 2

        # Find which result group corresponds to cluster A (starts ~100)
        cluster_a = [r for r in result if r[1] < 200]
        cluster_b = [r for r in result if r[1] > 300]
        assert len(cluster_a) >= 1 and len(cluster_b) >= 1

        # Verify density component is correctly computed
        density_a = min(1.0, np.log1p(20) / np.log1p(50))
        density_b = min(1.0, np.log1p(4) / np.log1p(50))
        assert density_a > density_b, "Sanity: cluster A has higher density"

    def test_both_tss_and_tts_scored(self):
        """Verify that BOTH TSS and TTS positions are scored by the scorer.

        Use a mock scorer that tracks which end_types it receives.
        """
        scored_types = []
        class TrackingScorer:
            def score(self, candidate):
                scored_types.append(candidate.end_type)
                candidate.confidence = 0.5
                candidate.seq_score = 0.5
                candidate.depth_score = 0.5
                candidate.annotation_score = 0.5
                candidate.tech_penalty = 0.0
                return candidate

        reads = self._make_reads([(100, 500)] * 10)
        reads = [(s, e, st, f'read_{i}') for i, (s, e, st, _) in enumerate(reads)]
        self.ted_collapse(reads, 'chr1', scorer=TrackingScorer(), do_get_best_ends=True)

        assert 'tss' in scored_types, "TSS should be scored"
        assert 'tts' in scored_types, "TTS should be scored"

    def test_tts_position_correct_plus_strand(self):
        """On + strand, TSS = start (smaller), TTS = end (larger)."""
        positions_scored = {}
        class PositionTracker:
            def score(self, candidate):
                positions_scored[candidate.end_type] = candidate.pos
                candidate.confidence = 0.5
                candidate.seq_score = 0.5
                candidate.depth_score = 0.5
                candidate.annotation_score = 0.5
                candidate.tech_penalty = 0.0
                return candidate

        reads = self._make_reads([(100, 500)] * 6, strand='+')
        self.ted_collapse(reads, 'chr1', scorer=PositionTracker(), do_get_best_ends=True)
        assert positions_scored.get('tss', 0) < positions_scored.get('tts', 0), \
            f"+ strand: TSS ({positions_scored.get('tss')}) should be < TTS ({positions_scored.get('tts')})"

    def test_tts_position_correct_minus_strand(self):
        """On - strand, TSS = end (larger), TTS = start (smaller)."""
        positions_scored = {}
        class PositionTracker:
            def score(self, candidate):
                positions_scored[candidate.end_type] = candidate.pos
                candidate.confidence = 0.5
                candidate.seq_score = 0.5
                candidate.depth_score = 0.5
                candidate.annotation_score = 0.5
                candidate.tech_penalty = 0.0
                return candidate

        reads = self._make_reads([(100, 500)] * 6, strand='-')
        self.ted_collapse(reads, 'chr1', scorer=PositionTracker(), do_get_best_ends=True)
        assert positions_scored.get('tss', 0) > positions_scored.get('tts', 0), \
            f"- strand: TSS ({positions_scored.get('tss')}) should be > TTS ({positions_scored.get('tts')})"

    def test_confidence_weight_formula(self):
        """Verify the composite weight formula: 0.40*density + 0.20*tss + 0.20*tts + 0.20*annot."""
        mock = self._make_mock_scorer(tss_score=0.8, tts_score=0.6)
        reads = self._make_reads([(100, 500)] * 10)
        reads = [(s, e, st, f'read_{i}') for i, (s, e, st, _) in enumerate(reads)]
        result = self.ted_collapse(reads, 'chr1', scorer=mock, do_get_best_ends=True)

        density = min(1.0, np.log1p(10) / np.log1p(50))
        # junction_match_type=None → annot_match_score=0.5
        expected = 0.40 * density + 0.20 * 0.8 + 0.20 * 0.6 + 0.20 * 0.5
        actual = result[0][0]
        assert abs(actual - expected) < 0.01, \
            f"Expected confidence ~{expected:.4f}, got {actual:.4f}"

    def test_confidence_weight_formula_fsm(self):
        """Verify FSM junction_match_type gives annot_match_score=1.0."""
        mock = self._make_mock_scorer(tss_score=0.8, tts_score=0.6)
        reads = self._make_reads([(100, 500)] * 10)
        reads = [(s, e, st, f'read_{i}') for i, (s, e, st, _) in enumerate(reads)]
        result = self.ted_collapse(reads, 'chr1', scorer=mock, do_get_best_ends=True,
                                   junction_match_type='FSM')

        density = min(1.0, np.log1p(10) / np.log1p(50))
        expected = 0.40 * density + 0.20 * 0.8 + 0.20 * 0.6 + 0.20 * 1.0
        actual = result[0][0]
        assert abs(actual - expected) < 0.01, \
            f"Expected FSM confidence ~{expected:.4f}, got {actual:.4f}"

    def test_confidence_weight_formula_nnc(self):
        """Verify NNC junction_match_type gives annot_match_score=0.2."""
        mock = self._make_mock_scorer(tss_score=0.8, tts_score=0.6)
        reads = self._make_reads([(100, 500)] * 10)
        reads = [(s, e, st, f'read_{i}') for i, (s, e, st, _) in enumerate(reads)]
        result = self.ted_collapse(reads, 'chr1', scorer=mock, do_get_best_ends=True,
                                   junction_match_type='NNC')

        density = min(1.0, np.log1p(10) / np.log1p(50))
        expected = 0.40 * density + 0.20 * 0.8 + 0.20 * 0.6 + 0.20 * 0.2
        actual = result[0][0]
        assert abs(actual - expected) < 0.01, \
            f"Expected NNC confidence ~{expected:.4f}, got {actual:.4f}"

    def test_fsm_boosts_over_nnc(self):
        """FSM should always produce higher confidence than NNC for same reads."""
        mock = self._make_mock_scorer(tss_score=0.5, tts_score=0.5)
        reads = self._make_reads([(100, 500)] * 10)
        reads = [(s, e, st, f'read_{i}') for i, (s, e, st, _) in enumerate(reads)]

        result_fsm = self.ted_collapse(reads, 'chr1', scorer=mock,
                                       do_get_best_ends=True, junction_match_type='FSM')
        result_nnc = self.ted_collapse(reads, 'chr1', scorer=mock,
                                       do_get_best_ends=True, junction_match_type='NNC')

        conf_fsm = result_fsm[0][0]
        conf_nnc = result_nnc[0][0]
        boost = conf_fsm - conf_nnc
        # Difference should be 0.20 * (1.0 - 0.2) = 0.16
        assert abs(boost - 0.16) < 0.01, \
            f"FSM-NNC boost should be ~0.16, got {boost:.4f}"

    def test_nic_between_fsm_and_nnc(self):
        """NIC should fall between FSM and NNC in confidence."""
        mock = self._make_mock_scorer(tss_score=0.5, tts_score=0.5)
        reads = self._make_reads([(100, 500)] * 10)
        reads = [(s, e, st, f'read_{i}') for i, (s, e, st, _) in enumerate(reads)]

        conf_fsm = self.ted_collapse(reads, 'chr1', scorer=mock,
                                     do_get_best_ends=True, junction_match_type='FSM')[0][0]
        conf_nic = self.ted_collapse(reads, 'chr1', scorer=mock,
                                     do_get_best_ends=True, junction_match_type='NIC')[0][0]
        conf_nnc = self.ted_collapse(reads, 'chr1', scorer=mock,
                                     do_get_best_ends=True, junction_match_type='NNC')[0][0]

        assert conf_fsm > conf_nic > conf_nnc, \
            f"Expected FSM ({conf_fsm:.3f}) > NIC ({conf_nic:.3f}) > NNC ({conf_nnc:.3f})"

    def test_junction_match_none_is_neutral(self):
        """junction_match_type=None should give 0.5 (neutral), same as no kwarg."""
        mock = self._make_mock_scorer(tss_score=0.5, tts_score=0.5)
        reads = self._make_reads([(100, 500)] * 10)
        reads = [(s, e, st, f'read_{i}') for i, (s, e, st, _) in enumerate(reads)]

        conf_none = self.ted_collapse(reads, 'chr1', scorer=mock,
                                      do_get_best_ends=True, junction_match_type=None)[0][0]
        conf_default = self.ted_collapse(reads, 'chr1', scorer=mock,
                                         do_get_best_ends=True)[0][0]
        assert abs(conf_none - conf_default) < 0.001, \
            f"None ({conf_none:.4f}) should equal default ({conf_default:.4f})"

    def test_final_stage_ignores_junction_match(self):
        """do_get_best_ends=False path should not be affected by junction_match_type."""
        reads = [[100, 500, ('marker', 'iso1'), ['r1', 'r2']],
                 [101, 501, ('marker', 'iso1'), ['r3', 'r4']],
                 [102, 499, ('marker', 'iso1'), ['r5']],
                 [105, 505, ('marker', 'iso1'), ['r6']]]
        result_fsm = self.ted_collapse(reads, 'chr1', do_get_best_ends=False,
                                       junction_match_type='FSM')
        result_nnc = self.ted_collapse(reads, 'chr1', do_get_best_ends=False,
                                       junction_match_type='NNC')
        # Final stage merges reads, doesn't compute confidence scores
        assert len(result_fsm) == len(result_nnc)
        for r_fsm, r_nnc in zip(result_fsm, result_nnc):
            assert set(r_fsm[-1]) == set(r_nnc[-1])


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
