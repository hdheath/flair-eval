"""Tests for bin/evaluation/alt_end_enumeration.py (Rank-3 per-gene #TSS/#PAS).

Validates: runs + expected TSV columns; distinct-site counting with tolerance;
multi-end fraction logic; TED-confidence (reality) stratification is present;
EMPTY input degrades gracefully (header-only TSV, no crash, exit 0).
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

import alt_end_enumeration as M  # noqa: E402


# ── unit: distinct-site counting ──────────────────────────────────────────────

def test_count_distinct_sites_tolerance():
    assert M.count_distinct_sites([], 50) == 0
    assert M.count_distinct_sites([100], 50) == 1
    # within tol → one site
    assert M.count_distinct_sites([100, 130, 145], 50) == 1
    # gap > tol → two sites
    assert M.count_distinct_sites([100, 200], 50) == 2
    # chained single-linkage: 100-140-180 each gap 40 <= 50 → one site
    assert M.count_distinct_sites([100, 140, 180], 50) == 1


# ── helpers to write tiny BED12 ───────────────────────────────────────────────

def _bed12_line(chrom, start, end, name, strand, blocks):
    """blocks = list of (rel_start, size)."""
    n = len(blocks)
    sizes = ",".join(str(s) for _, s in blocks) + ","
    starts = ",".join(str(rs) for rs, _ in blocks) + ","
    return "\t".join(map(str, [
        chrom, start, end, name, 0, strand, start, end, "0,0,0",
        n, sizes, starts])) + "\n"


def _write_bed(path, lines):
    path.write_text("".join(lines))


def test_enumeration_multi_tss_and_pas(tmp_path):
    """GENEA: 2 isoforms, same SJC, 2 distinct TSS → multi_tss.
    GENEB: 2 isoforms, same SJC, 2 distinct PAS → multi_pas.
    GENEC: 1 isoform → not multi-end."""
    bed = tmp_path / "iso.bed"
    lines = [
        # GENEA: + strand, junction chain shared, TSS differ (1000 vs 1300)
        _bed12_line("chr1", 1000, 5000, "tA1_GENEA", "+",
                    [(0, 100), (3900, 100)]),
        _bed12_line("chr1", 1300, 5000, "tA2_GENEA", "+",
                    [(0, 100), (3600, 100)]),
        # GENEB: + strand, same TSS, PAS differ (5000 vs 5400)
        _bed12_line("chr2", 1000, 5000, "tB1_GENEB", "+",
                    [(0, 100), (3900, 100)]),
        _bed12_line("chr2", 1000, 5400, "tB2_GENEB", "+",
                    [(0, 100), (3900, 200)]),
        # GENEC: single isoform
        _bed12_line("chr3", 1000, 5000, "tC1_GENEC", "+",
                    [(0, 100), (3900, 100)]),
    ]
    _write_bed(bed, lines)

    out = tmp_path / "out"
    rc = M.main([
        "--bed", f"ted:{bed}",
        "--output-prefix", "t", "--output-dir", str(out),
    ])
    assert rc == 0
    tsv = out / "t_alt_end_enumeration.tsv"
    assert tsv.exists()
    rows = list(csv.DictReader(tsv.open(), delimiter="\t"))
    assert list(rows[0].keys()) == M.TSV_COLUMNS

    by_gene = {r["gene_id"]: r for r in rows}
    assert by_gene["GENEA"]["n_tss"] == "2"
    assert by_gene["GENEA"]["multi_tss"] == "1"
    assert by_gene["GENEA"]["multi_end"] == "1"
    assert by_gene["GENEB"]["n_pas"] == "2"
    assert by_gene["GENEB"]["multi_pas"] == "1"
    assert by_gene["GENEC"]["multi_end"] == "0"
    # plot produced
    assert (out / "t_alt_end_enumeration.png").exists()


def test_confidence_stratification_present(tmp_path):
    """When a ted_log is supplied, high-conf counts are populated (>= 0),
    and a low-reality end is excluded from the high-conf repertoire."""
    bed = tmp_path / "iso.bed"
    # GENEA + strand: two TSS at 1000 (highconf) and 1300 (lowconf).
    lines = [
        _bed12_line("chr1", 1000, 5000, "tA1_GENEA", "+",
                    [(0, 100), (3900, 100)]),
        _bed12_line("chr1", 1300, 5000, "tA2_GENEA", "+",
                    [(0, 100), (3600, 100)]),
    ]
    _write_bed(bed, lines)

    log = tmp_path / "ted_log.tsv"
    log.write_text(
        "# hdbscan_params: ted_global=True\n"
        "junc_id\tcluster_id\tstage\tstatus\tn_reads\tjc_n_reads_total\t"
        "tss_pos\ttts_pos\tTED_tss_reality\tTED_tts_reality\n"
        # TSS 1000 high reality, TTS 5000 high reality
        "chr1:None-None:+:1j\tc0\tglobal_assign\tpass\t30\t40\t"
        "1000\t5000\t0.95\t0.90\n"
        # TSS 1300 LOW reality, same TTS
        "chr1:None-None:+:1j\tc1\tglobal_assign\tpass\t10\t40\t"
        "1300\t5000\t0.20\t0.90\n"
    )

    out = tmp_path / "out"
    rc = M.main([
        "--bed", f"ted:{bed}", "--ted-log", f"ted:{log}",
        "--output-prefix", "t", "--output-dir", str(out),
        "--reality-bin", "0.5", "--end-tol", "50",
    ])
    assert rc == 0
    rows = list(csv.DictReader(
        (out / "t_alt_end_enumeration.tsv").open(), delimiter="\t"))
    g = {r["gene_id"]: r for r in rows}["GENEA"]
    # full repertoire: 2 distinct TSS
    assert g["n_tss"] == "2"
    # high-confidence repertoire: only the reality>=0.5 TSS (1000) counts → 1
    assert g["n_tss_highconf"] == "1"
    assert int(g["n_pas_highconf"]) >= 1


def test_empty_input_graceful(tmp_path):
    """Empty BED → header-only TSV, no plot, exit 0 (never exit 1)."""
    bed = tmp_path / "empty.bed"
    bed.write_text("")
    out = tmp_path / "out"
    rc = M.main([
        "--bed", f"ted:{bed}",
        "--output-prefix", "e", "--output-dir", str(out),
    ])
    assert rc == 0
    tsv = out / "e_alt_end_enumeration.tsv"
    assert tsv.exists()
    rows = list(csv.reader(tsv.open(), delimiter="\t"))
    assert rows[0] == M.TSV_COLUMNS          # header present
    assert len(rows) == 1                    # header only
    assert not (out / "e_alt_end_enumeration.png").exists()  # no plot


def test_reality_at_nearest_within_tol():
    idx = {("tss", "chr1"): [(1000, 0.9), (1300, 0.2)]}
    # exact match
    assert M.reality_at("tss", "chr1", 1000, idx, 50) == 0.9
    # within tol picks nearest
    assert M.reality_at("tss", "chr1", 1320, idx, 50) == 0.2
    # outside tol → None
    assert M.reality_at("tss", "chr1", 2000, idx, 50) is None
    # missing axis/chrom → None
    assert M.reality_at("tts", "chr1", 1000, idx, 50) is None
