#!/usr/bin/env python3
"""Regression tests for evaluation/dropped_gene_characterization.py.

Pins the core analysis steps using a fixed 4-gene fixture:

  Reference GTF (tiny_biotype_annotation.gtf):
    GENEPC1     protein_coding              2 transcripts
    GENEPC2     protein_coding              1 transcript
    GENELNC1    lncRNA                      1 transcript
    GENEPSEUDO1 transcribed_unprocessed_pseudogene  1 transcript

  FLAIR baseline (tiny_dropped_flair.bed): isoforms for all 4 genes.

  TED mode (tiny_dropped_ted.bed): isoforms for only the protein_coding
  genes — the lncRNA and pseudogene are dropped.

Expected classification when running the script with --baseline-label FLAIR:
  GENEPC1     -> kept    (protein_coding)
  GENEPC2     -> kept    (protein_coding)
  GENELNC1    -> dropped (lncRNA)
  GENEPSEUDO1 -> dropped (pseudogene)

Drop rates per biotype:
  protein_coding   = 0/2 = 0.0%
  lncRNA           = 1/1 = 100%
  pseudogene class = 1/1 = 100%
"""
from __future__ import annotations

from pathlib import Path

import pytest

# conftest.py adds bin/ to sys.path
from evaluation.dropped_gene_characterization import (
    parse_gtf_gene_features,
    classify_per_mode,
    build_summary,
)
from evaluation.signal_utils import parse_isoforms


@pytest.fixture(scope="module")
def gtf_features(tiny_biotype_annotation_gtf):
    """Parse the GTF once for the whole test module."""
    return parse_gtf_gene_features(tiny_biotype_annotation_gtf)


@pytest.fixture(scope="module")
def isoforms_per_mode(tiny_dropped_flair_bed, tiny_dropped_ted_bed):
    """Load the FLAIR + TED isoform BEDs into the dict shape the
    analyzer expects."""
    return {
        "FLAIR": parse_isoforms(str(tiny_dropped_flair_bed)),
        "TED":   parse_isoforms(str(tiny_dropped_ted_bed)),
    }


# ── GTF parser sanity checks ──────────────────────────────────────


def test_gtf_has_4_genes_with_biotypes(gtf_features):
    """parse_gtf_gene_features must return all 4 genes from the fixture
    with their gene_type values."""
    gene_to_biotype, gene_to_lens = gtf_features
    assert set(gene_to_biotype) == {
        "GENEPC1", "GENEPC2", "GENELNC1", "GENEPSEUDO1",
    }
    assert gene_to_biotype["GENEPC1"] == "protein_coding"
    assert gene_to_biotype["GENEPC2"] == "protein_coding"
    assert gene_to_biotype["GENELNC1"] == "lncRNA"
    assert gene_to_biotype["GENEPSEUDO1"] == "transcribed_unprocessed_pseudogene"


def test_gtf_transcript_counts_per_gene(gtf_features):
    """Number of transcripts per gene must match the fixture."""
    _, gene_to_lens = gtf_features
    assert len(gene_to_lens["GENEPC1"]) == 2
    assert len(gene_to_lens["GENEPC2"]) == 1
    assert len(gene_to_lens["GENELNC1"]) == 1
    assert len(gene_to_lens["GENEPSEUDO1"]) == 1


def test_gtf_transcript_lengths(gtf_features):
    """Transcript lengths are summed exon sizes per transcript_id.
    Use the actual fixture coordinates — this is the cleanest pin
    against future GTF parser drift."""
    _, gene_to_lens = gtf_features
    # GENEPC1 TX_PC1a: exons 1001-1500 (500bp) + 2501-3500 (1000bp) = 1500
    # GENEPC1 TX_PC1b: exons 1001-1300 (300bp) + 2501-3500 (1000bp) = 1300
    assert sorted(gene_to_lens["GENEPC1"]) == [1300, 1500]
    # GENEPC2 TX_PC2a: 5001-5500 (500) + 6501-7000 (500) = 1000
    assert gene_to_lens["GENEPC2"] == [1000]
    # GENELNC1 TX_LNC1a: 10001-10300 (300) + 10501-10800 (300) = 600
    assert gene_to_lens["GENELNC1"] == [600]
    # GENEPSEUDO1 TX_PSEUDO1a: single exon 15001-15500 = 500
    assert gene_to_lens["GENEPSEUDO1"] == [500]


# ── classify_per_mode regression tests ────────────────────────────


def test_classify_returns_4_rows_for_one_mode(isoforms_per_mode, gtf_features):
    """One mode (TED) and 4 baseline genes -> 4 rows."""
    gene_to_biotype, gene_to_lens = gtf_features
    df = classify_per_mode(
        isoforms_per_mode, "FLAIR", gene_to_biotype, gene_to_lens
    )
    assert len(df) == 4
    assert set(df["mode"].unique()) == {"TED"}


def test_classify_marks_lncRNA_and_pseudogene_dropped(
    isoforms_per_mode, gtf_features
):
    """REGRESSION: TED keeps the protein_coding genes and drops the
    lncRNA + pseudogene. Pin both the status + biotype attachment."""
    gene_to_biotype, gene_to_lens = gtf_features
    df = classify_per_mode(
        isoforms_per_mode, "FLAIR", gene_to_biotype, gene_to_lens
    )
    by_gene = df.set_index("gene_id")
    assert by_gene.loc["GENEPC1", "status"] == "kept"
    assert by_gene.loc["GENEPC1", "biotype"] == "protein_coding"
    assert by_gene.loc["GENEPC2", "status"] == "kept"
    assert by_gene.loc["GENELNC1", "status"] == "dropped"
    assert by_gene.loc["GENELNC1", "biotype"] == "lncRNA"
    assert by_gene.loc["GENEPSEUDO1", "status"] == "dropped"
    assert (
        by_gene.loc["GENEPSEUDO1", "biotype"]
        == "transcribed_unprocessed_pseudogene"
    )


def test_classify_attaches_n_annotated_tx(isoforms_per_mode, gtf_features):
    """The n_annotated_tx column should reflect the GTF's transcript
    count per gene — important because it's used to detect
    'sparsely-annotated genes are dropped more' in real data."""
    gene_to_biotype, gene_to_lens = gtf_features
    df = classify_per_mode(
        isoforms_per_mode, "FLAIR", gene_to_biotype, gene_to_lens
    )
    by_gene = df.set_index("gene_id")
    assert by_gene.loc["GENEPC1", "n_annotated_tx"] == 2
    assert by_gene.loc["GENEPC2", "n_annotated_tx"] == 1
    assert by_gene.loc["GENELNC1", "n_annotated_tx"] == 1


def test_classify_attaches_mean_length(isoforms_per_mode, gtf_features):
    """The annotated_mean_len_bp column averages per-transcript spliced
    lengths. GENEPC1 has two transcripts of 1500/1300 bp -> mean 1400."""
    gene_to_biotype, gene_to_lens = gtf_features
    df = classify_per_mode(
        isoforms_per_mode, "FLAIR", gene_to_biotype, gene_to_lens
    )
    by_gene = df.set_index("gene_id")
    assert by_gene.loc["GENEPC1", "annotated_mean_len_bp"] == 1400
    assert by_gene.loc["GENEPC2", "annotated_mean_len_bp"] == 1000


# ── build_summary regression tests ────────────────────────────────


def test_summary_drop_rates_by_biotype(isoforms_per_mode, gtf_features):
    """REGRESSION: per-biotype drop rates should match expected ratios.
    protein_coding: 0/2 -> 0.0%
    lncRNA:         1/1 -> 100.0%
    pseudogene:     1/1 -> 100.0%
    """
    gene_to_biotype, gene_to_lens = gtf_features
    detail = classify_per_mode(
        isoforms_per_mode, "FLAIR", gene_to_biotype, gene_to_lens
    )
    summary = build_summary(detail)
    assert len(summary) == 1
    row = summary.iloc[0]
    assert row["mode"] == "TED"
    assert row["n_baseline_genes"] == 4
    assert row["n_dropped"] == 2
    assert row["n_kept"] == 2
    assert row["drop_rate_overall_pct"] == 50.0
    assert row["drop_rate_protein_coding_pct"] == 0.0
    assert row["drop_rate_lncRNA_pct"] == 100.0


def test_summary_dropped_kept_n_tx_medians(isoforms_per_mode, gtf_features):
    """REGRESSION: dropped-vs-kept median n_annotated_tx must reflect
    the dropped genes (LNC1=1tx, PSEUDO1=1tx -> median 1) vs kept genes
    (PC1=2tx, PC2=1tx -> median 1.5 floored to 1 by .median()).
    """
    gene_to_biotype, gene_to_lens = gtf_features
    detail = classify_per_mode(
        isoforms_per_mode, "FLAIR", gene_to_biotype, gene_to_lens
    )
    summary = build_summary(detail)
    row = summary.iloc[0]
    # Both groups have median=1 with this tiny fixture; confirm the
    # column populates rather than NaN.
    assert row["dropped_median_n_tx"] == 1
    assert row["kept_median_n_tx"] == 1


# ── 30-column FlairBed (FLG) path regression ──────────────────────
# The legacy fixture above is a 12-col BED whose name suffix happens to
# equal the GTF gene_id, so the FLG↔ENSG join is never exercised. The
# fixture below uses the real production schema: gene_id is an internal
# FLG… id (used for grouping) and the ENSG used for the biotype lookup
# lives in ref_gene_mapping (col 15), versioned (e.g. ENSG001.1). This
# pins the fix for the all-"unknown"/n_annotated_tx=0 collapse on real
# 30-col input.


@pytest.fixture(scope="module")
def flairbed_30col_gtf(tmp_path_factory) -> Path:
    """Tiny GTF keyed by ENSG (versioned) matching the 30-col FlairBed
    fixture's ref_gene_mapping ids. ENSG001 = protein_coding (2 tx),
    ENSG002 = lncRNA (1 tx). Note the ENSG version suffixes differ from
    the BED (.3 vs .1) to exercise version-suffix stripping on both
    sides of the join."""
    gtf = tmp_path_factory.mktemp("flg") / "flg_biotype.gtf"
    gtf.write_text(
        'chr1\tt\tgene\t1001\t5000\t.\t+\t.\tgene_id "ENSG001.3"; gene_type "protein_coding";\n'
        'chr1\tt\texon\t1001\t1200\t.\t+\t.\tgene_id "ENSG001.3"; transcript_id "ENST001.3";\n'
        'chr1\tt\texon\t4701\t5000\t.\t+\t.\tgene_id "ENSG001.3"; transcript_id "ENST001.3";\n'
        'chr1\tt\texon\t1001\t1200\t.\t+\t.\tgene_id "ENSG001.3"; transcript_id "ENST002.3";\n'
        'chr1\tt\texon\t4801\t5000\t.\t+\t.\tgene_id "ENSG001.3"; transcript_id "ENST002.3";\n'
        'chr2\tt\tgene\t8001\t12000\t.\t-\t.\tgene_id "ENSG002.3"; gene_type "lncRNA";\n'
        'chr2\tt\texon\t8001\t8400\t.\t-\t.\tgene_id "ENSG002.3"; transcript_id "ENST004.3";\n'
        'chr2\tt\texon\t11501\t12000\t.\t-\t.\tgene_id "ENSG002.3"; transcript_id "ENST004.3";\n'
    )
    return gtf


@pytest.fixture(scope="module")
def flairbed_30col_ted(tmp_path_factory, data_dir) -> Path:
    """A TED-mode BED derived from the 30-col FlairBed fixture that keeps
    only the FLG0001 (ENSG001, protein_coding) gene and drops FLG0002
    (ENSG002, lncRNA)."""
    src = (data_dir / "tiny_flairbed_30col.bed").read_text().splitlines()
    kept = [ln for ln in src if "\tFLG0001\t" in ln]
    out = tmp_path_factory.mktemp("flg") / "ted_30col.bed"
    out.write_text("\n".join(kept) + "\n")
    return out


def test_flg_path_joins_biotype_via_ref_gene_mapping(
    data_dir, flairbed_30col_gtf, flairbed_30col_ted
):
    """REGRESSION (BUG 1): on a real 30-col FlairBed, grouping is by the
    FLG gene_id but biotype/length must be looked up via the ENSG in
    ref_gene_mapping (version-stripped). Before the fix this collapsed to
    biotype='unknown'/n_annotated_tx=0 for every gene."""
    gene_to_biotype, gene_to_lens = parse_gtf_gene_features(flairbed_30col_gtf)
    isoforms_per_mode = {
        "FLAIR": parse_isoforms(str(data_dir / "tiny_flairbed_30col.bed")),
        "TED": parse_isoforms(str(flairbed_30col_ted)),
    }
    df = classify_per_mode(
        isoforms_per_mode, "FLAIR", gene_to_biotype, gene_to_lens
    )
    by_gene = df.set_index("gene_id")

    # Rows are keyed by the FLG grouping id, not the ENSG.
    assert set(by_gene.index) == {"FLG0001", "FLG0002"}

    # FLG0001 -> ENSG001 -> protein_coding, kept (present in TED).
    assert by_gene.loc["FLG0001", "biotype"] == "protein_coding"
    assert by_gene.loc["FLG0001", "status"] == "kept"
    assert by_gene.loc["FLG0001", "n_annotated_tx"] == 2

    # FLG0002 -> ENSG002 -> lncRNA, dropped (absent from TED).
    assert by_gene.loc["FLG0002", "biotype"] == "lncRNA"
    assert by_gene.loc["FLG0002", "status"] == "dropped"
    assert by_gene.loc["FLG0002", "n_annotated_tx"] == 1

    # The whole point of the fix: nothing collapsed to 'unknown'.
    assert "unknown" not in set(by_gene["biotype"])
