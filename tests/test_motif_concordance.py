"""Tests for Rank-1 motif-concordance-by-confidence-bin
(bin/evaluation/internal_priming.py --by-confidence; EXTENDS internal_priming).

Validates: hexamer/IP detection helpers; the concordance TSV has the expected
columns and bins; high-confidence 3' ends are MORE motif-concordant than the
low-confidence A-rich set (the assay-free "ends are real" signal); confidence
stratification is genuinely present; EMPTY/no-reality input degrades gracefully.
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

import internal_priming as IP  # noqa: E402


# ── unit: motif detection ─────────────────────────────────────────────────────

def test_pas_hexamer_detection():
    # genome with AATAAA at 90-96; PAS (+ strand) at 110 → scan [60,110) finds it
    seq = "C" * 90 + "AATAAA" + "C" * 14   # len 110
    genome = {"chr1": seq}
    assert IP.has_pas_hexamer("chr1", 110, "+", genome, 50) is True
    # PAS where upstream window is all C → no hexamer
    cgenome = {"chr1": "C" * 200}
    assert IP.has_pas_hexamer("chr1", 110, "+", cgenome, 50) is False
    # missing chrom → None
    assert IP.has_pas_hexamer("chrX", 110, "+", genome, 50) is None


def test_pas_hexamer_minus_strand():
    # minus strand: hexamer must appear as TTTATT (revcomp of AATAAA) downstream
    # in genomic coords so that after revcomp it reads AATAAA.
    seq = "C" * 110 + "TTTATT" + "C" * 50
    genome = {"chr1": seq}
    # PAS at genomic 110, - strand → scan genomic (110, 160], revcomp
    assert IP.has_pas_hexamer("chr1", 110, "-", genome, 50) is True


def test_reality_bin_index():
    assert IP.reality_bin_index(0.1) == 0
    assert IP.reality_bin_index(0.3) == 1
    assert IP.reality_bin_index(0.6) == 2
    assert IP.reality_bin_index(0.95) == 3
    assert IP.reality_bin_index(1.0) == 3


def _bed12_line(chrom, start, end, name, strand, blocks):
    n = len(blocks)
    sizes = ",".join(str(s) for _, s in blocks) + ","
    starts = ",".join(str(rs) for rs, _ in blocks) + ","
    return "\t".join(map(str, [
        chrom, start, end, name, 0, strand, start, end, "0,0,0",
        n, sizes, starts])) + "\n"


def _build_genome(length=4000):
    """Genome scaffold: mostly C, with an AATAAA placed just upstream of a
    'real' PAS at 3000, and an A-run (internal-priming) downstream of a
    'fake' PAS at 2000."""
    g = list("C" * length)
    # real PAS at 3000 (+ strand): AATAAA at 2975-2981 (upstream, within 50bp)
    for i, b in enumerate("AATAAA"):
        g[2975 + i] = b
    # fake PAS at 2000 (+ strand): A-run downstream 2000-2030 (internal priming)
    for i in range(2000, 2030):
        g[i] = "A"
    return {"chr1": "".join(g)}


def test_concordance_highconf_more_motif(tmp_path):
    """Two 3' ends on chr1:
      - PAS at 3000, HIGH reality (0.9): has AATAAA upstream → motif-concordant
      - PAS at 2000, LOW reality (0.15): A-run downstream → internal priming
    The TSV must show frac_with_motif higher in the high-reality bin and
    frac_a_rich higher in the low-reality bin."""
    genome_fa = tmp_path / "g.fa"
    g = _build_genome()
    genome_fa.write_text(">chr1\n" + g["chr1"] + "\n")

    bed = tmp_path / "iso.bed"
    bed.write_text(
        # multi-exon so n_exons>=2 (gene grouping not used here, but parse_isoforms
        # works on any BED). End (tts) at 3000 and 2000 on + strand.
        _bed12_line("chr1", 100, 3000, "real_GENEA", "+", [(0, 50), (2850, 50)]) +
        _bed12_line("chr1", 100, 2000, "fake_GENEB", "+", [(0, 50), (1850, 50)])
    )
    log = tmp_path / "ted_log.tsv"
    log.write_text(
        "# hdbscan_params: ted_global=True\n"
        "junc_id\tcluster_id\tstage\tstatus\tn_reads\tjc_n_reads_total\t"
        "tss_pos\ttts_pos\tTED_tss_reality\tTED_tts_reality\n"
        "chr1:None-None:+:1j\tc0\tglobal_assign\tpass\t30\t40\t100\t3000\t0.9\t0.90\n"
        "chr1:None-None:+:1j\tc1\tglobal_assign\tpass\t5\t40\t100\t2000\t0.9\t0.15\n"
    )

    out = tmp_path / "out"
    rc = IP.main([
        "--by-confidence",
        "--bed", f"ted:{bed}", "--ted-log", f"ted:{log}",
        "--genome", str(genome_fa), "--output", str(out),
        "--pas-hex-window", "50", "--window", "30", "--a-threshold", "0.6",
    ])
    assert rc == 0
    tsv = out / "motif_concordance_by_confidence.tsv"
    assert tsv.exists()
    rows = list(csv.DictReader(tsv.open(), delimiter="\t"))
    assert list(rows[0].keys()) == [
        "method", "axis", "reality_bin", "n_ends",
        "frac_with_motif", "frac_a_rich"]

    tts = {r["reality_bin"]: r for r in rows if r["axis"] == "tts"}
    high = tts["[0.75,1]"]
    low = tts["[0,0.25)"]
    # high-conf end: has AATAAA, no A-run
    assert high["n_ends"] == "1"
    assert float(high["frac_with_motif"]) == 1.0
    assert float(high["frac_a_rich"]) == 0.0
    # low-conf end: A-run, no canonical hexamer
    assert low["n_ends"] == "1"
    assert float(low["frac_with_motif"]) == 0.0
    assert float(low["frac_a_rich"]) == 1.0
    # the headline: motif fraction rises with confidence
    assert float(high["frac_with_motif"]) > float(low["frac_with_motif"])
    assert (out / "motif_concordance_by_confidence.png").exists()


def test_empty_input_graceful(tmp_path):
    """Empty BED → header-only TSV, no plot, exit 0."""
    genome_fa = tmp_path / "g.fa"
    genome_fa.write_text(">chr1\n" + "C" * 1000 + "\n")
    bed = tmp_path / "empty.bed"
    bed.write_text("")
    log = tmp_path / "log.tsv"
    log.write_text("# ted_global=True\njunc_id\ttss_pos\ttts_pos\n")
    out = tmp_path / "out"
    rc = IP.main([
        "--by-confidence", "--bed", f"ted:{bed}", "--ted-log", f"ted:{log}",
        "--genome", str(genome_fa), "--output", str(out),
    ])
    assert rc == 0
    tsv = out / "motif_concordance_by_confidence.tsv"
    assert tsv.exists()
    rows = list(csv.reader(tsv.open(), delimiter="\t"))
    assert rows[0] == ["method", "axis", "reality_bin", "n_ends",
                       "frac_with_motif", "frac_a_rich"]
    assert not (out / "motif_concordance_by_confidence.png").exists()


def test_no_ted_log_skips_method_gracefully(tmp_path):
    """No reality available → method skipped, header-only TSV, exit 0 (the
    analysis explicitly requires the per-end confidence)."""
    genome_fa = tmp_path / "g.fa"
    genome_fa.write_text(">chr1\n" + "C" * 4000 + "\n")
    bed = tmp_path / "iso.bed"
    bed.write_text(_bed12_line("chr1", 100, 3000, "x_G", "+",
                               [(0, 50), (2850, 50)]))
    out = tmp_path / "out"
    rc = IP.main([
        "--by-confidence", "--bed", f"ted:{bed}",
        "--genome", str(genome_fa), "--output", str(out),
    ])
    assert rc == 0
    assert (out / "motif_concordance_by_confidence.tsv").exists()
