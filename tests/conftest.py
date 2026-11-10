"""Shared pytest fixtures for the flair4_ted test suite.

Centralizes all fixture-data paths so individual test files don't have to
hardcode `BASE = Path("/private/home/...")` — they just request the
fixture they need by name and the path is resolved relative to this file's
location.

Pattern: every fixture is a `Path` object pointing at a file under
`tests/data/`. Tests that need to load the fixture do so themselves via
the appropriate parser; conftest only resolves paths.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest


# ─── Directory anchors ──────────────────────────────────────────────────

REPO_ROOT = Path(__file__).resolve().parent.parent
TESTS_DIR = REPO_ROOT / "tests"
DATA_DIR = TESTS_DIR / "data"
BIN_DIR = REPO_ROOT / "bin"

# Sister repo: /private/home/hdheath/projects/ted/ holds plot_per_locus.py
# and other ted-specific helpers. test_plot_per_locus_loaders imports from
# there. Resolved here once so individual tests don't have to.
TED_REPO_BIN = (REPO_ROOT.parent / "ted" / "bin").resolve()


# Add bin/ to sys.path so `from evaluation.<...> import ...` works in
# test modules without each test having to do it themselves.
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))


# ─── Path fixtures ──────────────────────────────────────────────────────
# Each returns the absolute path to a fixture file. Tests parse them as
# needed (parse_isoforms, parse_peaks_bed, pd.read_csv, etc.).


@pytest.fixture(scope="session")
def data_dir() -> Path:
    """Absolute path to tests/data/."""
    return DATA_DIR


@pytest.fixture(scope="session")
def ted_repo_bin() -> Path:
    """Absolute path to the sister ted/ repo's bin/ — for cross-repo tests
    (e.g. test_plot_per_locus_loaders.py)."""
    return TED_REPO_BIN


# Existing fixtures (already present in tests/data/, owned by older tests)
@pytest.fixture(scope="session")
def tiny_isoforms_bed() -> Path:
    """4 isoforms in 2 SJC groups — used by older tests; reused here for
    sanity-check parses."""
    return DATA_DIR / "tiny_isoforms.bed"


@pytest.fixture(scope="session")
def tiny_cage_peaks() -> Path:
    return DATA_DIR / "tiny_cage_peaks.bed"


@pytest.fixture(scope="session")
def tiny_drna_peaks() -> Path:
    return DATA_DIR / "tiny_drna_peaks.bed"


@pytest.fixture(scope="session")
def tiny_genome_fa() -> Path:
    return DATA_DIR / "tiny_genome.fa"


@pytest.fixture(scope="session")
def tiny_annotation_gtf() -> Path:
    return DATA_DIR / "tiny_annotation.gtf"


# New fixtures added by the test scaffolding plan (see tests/data/).
# These cover global-mode edge cases the older fixtures don't exercise.
@pytest.fixture(scope="session")
def tiny_ted_log_global() -> Path:
    """A 9-row ted_log.tsv mirroring --ted_global mode's multi-stage
    layout. Includes:
      - rows where the same `cluster_id` (c0, c1) appears with DIFFERENT
        (tss_pos, tts_pos) across stages (reproduces the dedup bug fixed
        in plot_per_locus.py:load_ted_log).
      - one row with tss_pos=-1 and another with tts_pos=-1 (reproduces
        the sentinel-handling case for ted_pr_sweep's filter chain).
      - one global_peak partition-summary row (reproduces the stage
        filter ted_pr_sweep applies; load_ted_log itself does NOT filter).
    Used by test_plot_per_locus_loaders.py and
    test_ted_pr_sweep_loaders.py.
    """
    return DATA_DIR / "tiny_ted_log_global.tsv"


@pytest.fixture(scope="session")
def tiny_redundant_isoforms_bed() -> Path:
    """4 isoforms in one SJC group, designed for the
    sjc_alt_end_analysis.py classification regression test:
      - 2 isos sharing (CAGE_X, dRNA_Y) — both ortho peaks present;
        these should classify as `fully_redundant` (post-fix).
      - 2 isos sharing (CAGE_X, None) — TTS misses dRNA window; these
        should classify as `unsupported` (post-fix; pre-fix they were
        mis-classified as `fully_redundant`).
    Reuses the existing tiny_cage_peaks.bed and tiny_drna_peaks.bed
    on chr_test for ortho lookups.
    """
    return DATA_DIR / "tiny_redundant_isoforms.bed"


# Fixtures for dropped_gene_characterization regression test.
# Triple of (GTF with biotypes, FLAIR-baseline BED, TED BED that drops
# the lncRNA + pseudogene genes). Used by
# tests/test_dropped_gene_characterization.py to pin biotype-based
# classification.
@pytest.fixture(scope="session")
def tiny_biotype_annotation_gtf() -> Path:
    """4-gene GTF on chr_test with explicit gene_type tags:
    GENE_PC1 / GENE_PC2 = protein_coding
    GENE_LNC1           = lncRNA
    GENE_PSEUDO1        = transcribed_unprocessed_pseudogene
    Used to verify the analyzer correctly attaches biotype + transcript
    length to each gene.
    """
    return DATA_DIR / "tiny_biotype_annotation.gtf"


@pytest.fixture(scope="session")
def tiny_dropped_flair_bed() -> Path:
    """5 isoforms covering all 4 genes from the biotype GTF — represents
    the FLAIR baseline that has assembled every gene at least once."""
    return DATA_DIR / "tiny_dropped_flair.bed"


@pytest.fixture(scope="session")
def tiny_dropped_ted_bed() -> Path:
    """2 isoforms covering only the protein_coding genes — represents a
    TED mode that has dropped the lncRNA + pseudogene. The expected
    classification: GENE_PC1/GENE_PC2 = kept, GENE_LNC1/GENE_PSEUDO1 =
    dropped."""
    return DATA_DIR / "tiny_dropped_ted.bed"
