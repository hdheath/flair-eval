"""Tests for bin/flair_eval_lib/variation.py (CONTRACTS A32)."""

import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "bin"))

from flair_eval_lib.annotation import Transcript  # noqa: E402
from flair_eval_lib.variation import gene_class, variation_section  # noqa: E402


def t(gene, exons, strand="+", chrom="chr1", name=None):
    return Transcript(name or f"{gene}{exons}", gene, chrom, strand, tuple(exons))


CHAIN1 = [(100, 200), (300, 400)]
CHAIN2 = [(100, 200), (350, 400)]


def test_ends_only():
    assert gene_class([t("g", CHAIN1), t("g", [(30, 200), (300, 400)])]) == "ends_only"
    assert gene_class([t("g", [(0, 500)]), t("g", [(0, 700)])]) == "ends_only"     # mono


def test_splice_only():
    assert gene_class([t("g", CHAIN1), t("g", CHAIN2)]) == "splice_only"
    # same chain twice, ends within 50 bp: still splice only
    assert gene_class([t("g", CHAIN1), t("g", [(140, 200), (300, 440)]),
                       t("g", CHAIN2)]) == "splice_only"
    # mono-exon next to a spliced isoform is a second chain
    assert gene_class([t("g", CHAIN1), t("g", [(100, 400)])]) == "splice_only"


def test_both():
    assert gene_class([t("g", CHAIN1), t("g", [(100, 200), (300, 451)]),
                       t("g", CHAIN2)]) == "both"


def test_neither():
    assert gene_class([t("g", CHAIN1), t("g", [(110, 200), (300, 390)])]) == "neither"


def test_exactly_50_bp_is_not_different():
    assert gene_class([t("g", CHAIN1), t("g", [(50, 200), (300, 400)])]) == "neither"
    assert gene_class([t("g", CHAIN1), t("g", [(49, 200), (300, 400)])]) == "ends_only"


def test_minus_strand_uses_strand_aware_ends():
    a = t("g", CHAIN1, "-")
    b = t("g", [(100, 200), (300, 460)], "-")       # 5' end on - moves by 60
    assert gene_class([a, b]) == "ends_only"


def test_section_counts_and_main_chromosomes():
    isos = [t("a", CHAIN1), t("a", CHAIN2),                         # splice_only
            t("b", CHAIN1), t("b", [(0, 200), (300, 400)]),         # ends_only
            t("c", CHAIN1),                                        # single isoform: skipped
            t("m", CHAIN1, chrom="chrM"), t("m", CHAIN2, chrom="chrM")]
    assert variation_section(isos) == {"genes_multi": 2, "ends_only": 1, "splice_only": 1,
                                       "both": 0, "neither": 0}
