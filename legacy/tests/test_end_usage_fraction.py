"""Tests for bin/evaluation/end_usage_fraction.py (Rank-2 usage-fraction/APA).

Validates: runs + expected TSV columns; read-weighted usage fractions;
APA/ATSS index + PDUI; TED-reality stratification present; EMPTY input
degrades gracefully (header-only TSVs, no crash, exit 0).
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parent.parent / "bin"
if str(BIN) not in sys.path:
    sys.path.insert(0, str(BIN))
if str(BIN / "evaluation") not in sys.path:
    sys.path.insert(0, str(BIN / "evaluation"))

import end_usage_fraction as M  # noqa: E402


def _bed12_line(chrom, start, end, name, strand, blocks):
    n = len(blocks)
    sizes = ",".join(str(s) for _, s in blocks) + ","
    starts = ",".join(str(rs) for rs, _ in blocks) + ","
    return "\t".join(map(str, [
        chrom, start, end, name, 0, strand, start, end, "0,0,0",
        n, sizes, starts])) + "\n"


# ── unit: collapse_sites read-weighting ───────────────────────────────────────

def test_collapse_sites_weighting():
    # two positions within tol → one site, reads summed, rep = higher-support
    sites = M.collapse_sites([(1000, 5), (1020, 20)], 50)
    assert sites == [(1020, 25)]
    # far apart → two sites
    sites = M.collapse_sites([(1000, 5), (2000, 7)], 50)
    assert sorted(sites) == [(1000, 5), (2000, 7)]
    assert M.collapse_sites([], 50) == []


def test_usage_fraction_and_indices(tmp_path):
    """GENEA + strand, shared SJC, 2 PAS at 5000 (30 reads) and 5400 (10 reads).
    Single TSS at 1000. → APA index = 1 - 0.75 = 0.25; PDUI (distal=5400) = 0.25."""
    bed = tmp_path / "iso.bed"
    bed.write_text(
        _bed12_line("chr1", 1000, 5000, "tA1_GENEA", "+", [(0, 100), (3900, 100)]) +
        _bed12_line("chr1", 1000, 5400, "tA2_GENEA", "+", [(0, 100), (3900, 200)])
    )
    rmap = tmp_path / "reads.map.txt"
    rmap.write_text(
        "tA1_GENEA\t" + ",".join(f"r{i}" for i in range(30)) + "\n" +
        "tA2_GENEA\t" + ",".join(f"s{i}" for i in range(10)) + "\n"
    )

    out = tmp_path / "out"
    rc = M.main([
        "--bed", f"ted:{bed}", "--read-map", f"ted:{rmap}",
        "--output-prefix", "t", "--output-dir", str(out),
    ])
    assert rc == 0
    gene_tsv = out / "t_end_usage_fraction.gene.tsv"
    site_tsv = out / "t_end_usage_fraction.site.tsv"
    assert gene_tsv.exists() and site_tsv.exists()

    grows = list(csv.DictReader(gene_tsv.open(), delimiter="\t"))
    assert list(grows[0].keys()) == M.GENE_TSV_COLUMNS
    g = {r["gene_id"]: r for r in grows}["GENEA"]
    assert g["n_pas"] == "2"
    assert g["n_tss"] == "1"
    assert abs(float(g["apa_index"]) - 0.25) < 1e-6
    assert abs(float(g["atss_index"]) - 0.0) < 1e-6
    # distal PAS (5400) carries 10/40 = 0.25 of usage
    assert abs(float(g["pdui"]) - 0.25) < 1e-6

    srows = list(csv.DictReader(site_tsv.open(), delimiter="\t"))
    assert list(srows[0].keys()) == M.SITE_TSV_COLUMNS
    pas_fracs = sorted(float(r["usage_fraction"]) for r in srows
                       if r["axis"] == "pas")
    assert pas_fracs == [0.25, 0.75]
    assert (out / "t_end_usage_fraction.png").exists()


def test_minus_strand_pdui(tmp_path):
    """- strand: distal PAS = smallest genomic coord. TSS at end, TTS at start.

    Both isoforms must GENUINELY share the same junction chain (donor 1300,
    acceptor 4000) so the two distinct PAS (genomic 1000 vs 1100) are
    alternative ends of the SAME spliced structure.  This was the case the
    SJC-grouping fix (FIX 4) corrects: the indices are computed within the
    chain, so distal (1000) usage = 0.5.
    """
    bed = tmp_path / "iso.bed"
    # Shared junction (1300, 4000): exon1 ends at 1300, exon2 = 4000..5000.
    # PAS (genomic start, - strand) at 1000 (distal) and 1100 (proximal).
    bed.write_text(
        _bed12_line("chr1", 1000, 5000, "tA1_GENEA", "-", [(0, 300), (3000, 1000)]) +
        _bed12_line("chr1", 1100, 5000, "tA2_GENEA", "-", [(0, 200), (2900, 1000)])
    )
    out = tmp_path / "out"
    rc = M.main(["--bed", f"ted:{bed}", "--output-prefix", "t",
                 "--output-dir", str(out)])
    assert rc == 0
    grows = list(csv.DictReader(
        (out / "t_end_usage_fraction.gene.tsv").open(), delimiter="\t"))
    g = {r["gene_id"]: r for r in grows}["GENEA"]
    assert g["strand"] == "-"
    assert g["n_pas"] == "2"
    # equal weight (no read map → 1 each) → distal (1000) usage = 0.5
    assert abs(float(g["pdui"]) - 0.5) < 1e-6


def test_sjc_grouping_distinct_from_gene_pooling(tmp_path):
    """FIX 4: ATSS index is computed WITHIN junction chains, then pooled.

    One gene, two distinct junction chains:
      Chain A (junction 2500-4000): 2 isoforms with distinct TSS
                                     (1000, 2000) → within-chain ATSS = 0.5.
      Chain B (junction 3500-4200): 2 isoforms both with TSS 3000
                                     → within-chain ATSS = 0.0.

    SJC-grouped (correct): read-weighted mean = (0.5*2 + 0.0*2)/4 = 0.25.
    Old gene-pooled (wrong): all TSS pooled → sites 1000(1),2000(1),3000(2);
        dominant = 3000 (0.5) → ATSS = 0.5.

    Asserting 0.25 proves SJC grouping is in effect (not the pooled 0.5).
    """
    bed = tmp_path / "iso.bed"
    bed.write_text(
        # Chain A: shared junction (2500, 4000); TSS 1000 and 2000.
        _bed12_line("chr1", 1000, 5000, "tA1_GENEX", "+", [(0, 1500), (3000, 1000)]) +
        _bed12_line("chr1", 2000, 5000, "tA2_GENEX", "+", [(0, 500), (2000, 1000)]) +
        # Chain B: shared junction (3500, 4200); both TSS 3000.
        _bed12_line("chr1", 3000, 5000, "tB1_GENEX", "+", [(0, 500), (1200, 800)]) +
        _bed12_line("chr1", 3000, 5000, "tB2_GENEX", "+", [(0, 500), (1200, 800)])
    )
    out = tmp_path / "out"
    rc = M.main(["--bed", f"ted:{bed}", "--output-prefix", "t",
                 "--output-dir", str(out)])
    assert rc == 0
    grows = list(csv.DictReader(
        (out / "t_end_usage_fraction.gene.tsv").open(), delimiter="\t"))
    g = {r["gene_id"]: r for r in grows}["GENEX"]
    # SJC-grouped ATSS = 0.25 (NOT the gene-pooled 0.5).
    assert abs(float(g["atss_index"]) - 0.25) < 1e-6, (
        f"expected SJC-grouped ATSS 0.25, got {g['atss_index']} "
        "(0.5 would mean gene-level pooling — the bug)")


def test_reality_attached_to_sites(tmp_path):
    bed = tmp_path / "iso.bed"
    bed.write_text(
        _bed12_line("chr1", 1000, 5000, "tA1_GENEA", "+", [(0, 100), (3900, 100)]) +
        _bed12_line("chr1", 1000, 5400, "tA2_GENEA", "+", [(0, 100), (3900, 200)])
    )
    log = tmp_path / "ted_log.tsv"
    log.write_text(
        "# hdbscan_params: ted_global=True\n"
        "junc_id\tcluster_id\tstage\tstatus\tn_reads\tjc_n_reads_total\t"
        "tss_pos\ttts_pos\tTED_tss_reality\tTED_tts_reality\n"
        "chr1:None-None:+:1j\tc0\tglobal_assign\tpass\t30\t40\t1000\t5000\t0.9\t0.88\n"
        "chr1:None-None:+:1j\tc1\tglobal_assign\tpass\t10\t40\t1000\t5400\t0.9\t0.30\n"
    )
    out = tmp_path / "out"
    rc = M.main(["--bed", f"ted:{bed}", "--ted-log", f"ted:{log}",
                 "--output-prefix", "t", "--output-dir", str(out)])
    assert rc == 0
    srows = list(csv.DictReader(
        (out / "t_end_usage_fraction.site.tsv").open(), delimiter="\t"))
    pas = {int(r["site_pos"]): r["ted_reality"] for r in srows if r["axis"] == "pas"}
    # both PAS got a reality value attached (stratification present)
    assert pas[5000] != ""
    assert pas[5400] != ""
    assert abs(float(pas[5400]) - 0.30) < 1e-6


def test_empty_input_graceful(tmp_path):
    bed = tmp_path / "empty.bed"
    bed.write_text("")
    out = tmp_path / "out"
    rc = M.main(["--bed", f"ted:{bed}", "--output-prefix", "e",
                 "--output-dir", str(out)])
    assert rc == 0
    gene_tsv = out / "e_end_usage_fraction.gene.tsv"
    site_tsv = out / "e_end_usage_fraction.site.tsv"
    assert gene_tsv.exists() and site_tsv.exists()
    grows = list(csv.reader(gene_tsv.open(), delimiter="\t"))
    assert grows[0] == M.GENE_TSV_COLUMNS
    assert len(grows) == 1
    assert not (out / "e_end_usage_fraction.png").exists()
