"""Tests for the judge's per-chain columns (ted_end_precision.py, ted-clustering map tickets 05 and 12).

The judge scores each junction chain's distinct called ends one-to-one against the truth calls,
counts recall over read-supported truth calls (--read-ends, the corrected spliced reads of a
genome-record dump), and scores alternative-end pairs in the 51-200 bp and over-200 bp bins.
Every case is hand-built; coordinates are BED-style (0-based, half-open), and a truth call is a
1-bp reference point."""

import gzip
import os
import sys
import tarfile

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
VENDORED = os.path.join(ROOT, "modules", "end_accuracy", "resources", "usr", "bin")
sys.path.insert(0, VENDORED)

pytest.importorskip("scipy", reason="the judge's signal_utils imports scipy")

import ted_end_precision as judge  # noqa: E402

DUMP_COLUMNS = ["read_name", "chrom", "strand", "start", "end", "aligned_start", "aligned_end",
                "five_prime", "three_prime", "aligned_five_prime", "aligned_three_prime", "introns",
                "junctions_from_annotation"]


def bed12(chrom, strand, exons, score=1, name="iso"):
    """A BED12 line for an isoform with these exons (list of (start, end))."""
    start, end = exons[0][0], exons[-1][1]
    sizes = "".join(f"{e - s}," for s, e in exons)
    offsets = "".join(f"{s - start}," for s, _ in exons)
    return "\t".join(map(str, (chrom, start, end, name, score, strand, start, start, 0, len(exons),
                               sizes, offsets)))


def isoforms(tmp_path, lines):
    path = tmp_path / "isoforms.bed"
    path.write_text("".join(line + "\n" for line in lines))
    return judge.parse_isoforms_bed(str(path))


def points(*calls):
    """A peaks dict of 1-bp truth calls, each (chrom, strand, point)."""
    peaks = {}
    for chrom, strand, p in calls:
        peaks.setdefault((chrom, strand), []).append((p, p + 1))
    for k in peaks:
        peaks[k].sort()
    return peaks


def read(chrom, strand, five, three, introns=()):
    return judge.ReadEnd(chrom, strand, tuple(introns), five, three)


# The two chains used below: both start at exon (1000, 1100) on +.
CHAIN_A = ((1100, 1500),)
CHAIN_B = ((1100, 1700),)


def exons_a(tss, tts=1600):
    return [(tss, 1100), (1500, tts)]


def exons_b(tss, tts=1800):
    return [(tss, 1100), (1700, tts)]


# ------------------------------------------------------------------ headline: one-to-one per chain
def test_two_chains_sharing_a_tss_both_get_credit(tmp_path):
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000)), bed12("chr1", "+", exons_b(1000))])
    m = judge.compute_chain_metrics(iso, points(("chr1", "+", 1000)), None, 50)
    assert m["5prime_chain_ends"] == 2
    assert m["5prime_chain_tp"] == 2
    assert m["5prime_chain_precision"] == 1.0
    assert m["5prime_chain_sites_hit"] == 1
    assert m["5prime_chain_recall_all_calls"] == 1.0


def test_two_calls_on_one_chain_near_one_truth_call_count_once(tmp_path):
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000)), bed12("chr1", "+", exons_a(1030))])
    peaks = points(("chr1", "+", 1010))
    m = judge.compute_chain_metrics(iso, peaks, None, 50)
    assert (m["5prime_chain_ends"], m["5prime_chain_tp"]) == (2, 1)
    assert m["5prime_chain_precision"] == 0.5
    # the old many-to-one headline credits both ends
    old = judge.compute_jc_deduplicated_precision_recall(iso, {}, [], 50, peaks_5prime=peaks)
    assert old["5prime_dedup_precision"] == 1.0


def test_one_chain_end_counted_once_however_many_isoforms_share_it(tmp_path):
    """Distinct called ends: two isoforms with one TSS and different TTSs are one 5' end."""
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000, 1600)),
                              bed12("chr1", "+", exons_a(1000, 1900))])
    m = judge.compute_chain_metrics(iso, points(("chr1", "+", 1000)), points(("chr1", "+", 1900)), 50)
    assert (m["5prime_chain_ends"], m["5prime_chain_tp"]) == (1, 1)
    assert (m["3prime_chain_ends"], m["3prime_chain_tp"]) == (2, 1)
    assert m["3prime_chain_precision"] == 0.5


def test_no_peaks_for_an_axis_leaves_its_columns_empty(tmp_path):
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000))])
    m = judge.compute_chain_metrics(iso, points(("chr1", "+", 1000)), None, 50)
    assert m["3prime_chain_precision"] is None and m["3prime_chain_ends"] is None


# ------------------------------------------------------------------ read-supported recall
def test_recall_counts_only_read_supported_truth_calls(tmp_path):
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000))])
    peaks = points(("chr1", "+", 1000), ("chr1", "+", 5000))
    reads = [read("chr1", "+", 1010, 1600, CHAIN_A),
             read("chr1", "-", 1600, 5000, CHAIN_A)]   # the wrong strand supports nothing
    m = judge.compute_chain_metrics(iso, peaks, None, 50, read_ends=reads)
    assert m["5prime_chain_recall_all_calls"] == 0.5
    assert m["5prime_read_supported_sites"] == 1
    assert m["5prime_chain_sites_hit_read_supported"] == 1
    assert m["5prime_chain_recall_read_supported"] == 1.0


def test_read_support_window_is_the_hit_window(tmp_path):
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000))])
    peaks = points(("chr1", "+", 1000), ("chr1", "+", 3000))
    reads = [read("chr1", "+", 1050, 1600, CHAIN_A), read("chr1", "+", 3051, 3600, CHAIN_A)]
    m = judge.compute_chain_metrics(iso, peaks, None, 50, read_ends=reads)
    assert m["5prime_read_supported_sites"] == 1   # 50 bp supports, 51 bp does not


def test_without_read_ends_read_supported_columns_are_empty(tmp_path):
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000))])
    m = judge.compute_chain_metrics(iso, points(("chr1", "+", 1000)), None, 50)
    assert m["5prime_chain_recall_read_supported"] is None
    assert m["5prime_read_supported_sites"] is None
    assert m["5prime_altpair_51_200_recall"] is None
    assert m["5prime_altpair_51_200_truth_pairs"] is None


# ------------------------------------------------------------------ alternative-end pairs
TANDEM = points(("chr1", "+", 1000), ("chr1", "+", 1120))
TANDEM_READS = [read("chr1", "+", 1000, 1600, CHAIN_A)] * 3 + [read("chr1", "+", 1120, 1600, CHAIN_A)] * 2


def test_tandem_pair_single_linkage_merges_is_missed(tmp_path):
    """Single linkage (200 bp 5' gap) merges TSSs 120 bp apart into one end call."""
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000), score=5)])
    m = judge.compute_chain_metrics(iso, TANDEM, None, 50, read_ends=TANDEM_READS)
    assert m["5prime_altpair_51_200_truth_pairs"] == 1
    assert m["5prime_altpair_51_200_pairs_hit"] == 0
    assert m["5prime_altpair_51_200_recall"] == 0.0
    assert m["5prime_altpair_51_200_alt_ends"] == 0
    assert m["5prime_altpair_51_200_precision"] is None
    assert m["5prime_altpair_over_200_truth_pairs"] == 0
    assert m["5prime_altpair_over_200_recall"] is None


def test_tandem_pair_a_split_recovers(tmp_path):
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000), score=3),
                              bed12("chr1", "+", exons_a(1120), score=2)])
    m = judge.compute_chain_metrics(iso, TANDEM, None, 50, read_ends=TANDEM_READS)
    assert m["5prime_altpair_51_200_pairs_hit"] == 1
    assert m["5prime_altpair_51_200_recall"] == 1.0
    assert (m["5prime_altpair_51_200_alt_ends"], m["5prime_altpair_51_200_alt_ends_hit"]) == (1, 1)
    assert m["5prime_altpair_51_200_precision"] == 1.0


def test_pair_calls_must_be_supported_by_one_chain(tmp_path):
    """Each call of a pair needs a read of the same junction chain ending near it."""
    reads = [read("chr1", "+", 1000, 1600, CHAIN_A), read("chr1", "+", 1120, 1800, CHAIN_B)]
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000)), bed12("chr1", "+", exons_a(1120))])
    m = judge.compute_chain_metrics(iso, TANDEM, None, 50, read_ends=reads)
    assert m["5prime_altpair_51_200_truth_pairs"] == 0


def test_pair_is_hit_only_by_its_own_chain(tmp_path):
    """Two chains each calling one end of the pair do not recover it."""
    reads = TANDEM_READS + [read("chr1", "+", 1120, 1800, CHAIN_B)]
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000)), bed12("chr1", "+", exons_b(1120))])
    m = judge.compute_chain_metrics(iso, TANDEM, None, 50, read_ends=reads)
    assert m["5prime_altpair_51_200_truth_pairs"] == 1
    assert m["5prime_altpair_51_200_pairs_hit"] == 0


def test_one_read_near_both_calls_does_not_make_a_pair(tmp_path):
    """A read ends at one call, its nearest: reads that all end between two calls 70 bp apart
    (each read within 50 bp of both) support only the nearer call."""
    peaks = points(("chr1", "+", 1000), ("chr1", "+", 1070))
    reads = [read("chr1", "+", 1030, 1600, CHAIN_A)] * 3
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000)), bed12("chr1", "+", exons_a(1070))])
    m = judge.compute_chain_metrics(iso, peaks, None, 50, read_ends=reads)
    assert m["5prime_read_supported_sites"] == 2   # a read near a call still supports it for recall
    assert m["5prime_altpair_51_200_truth_pairs"] == 0
    reads.append(read("chr1", "+", 1060, 1600, CHAIN_A))   # a read nearest the other call
    m = judge.compute_chain_metrics(iso, peaks, None, 50, read_ends=reads)
    assert m["5prime_altpair_51_200_truth_pairs"] == 1
    assert m["5prime_altpair_51_200_pairs_hit"] == 1


def test_calls_50_bp_apart_are_not_a_pair(tmp_path):
    peaks = points(("chr1", "+", 1000), ("chr1", "+", 1050))
    reads = [read("chr1", "+", 1000, 1600, CHAIN_A), read("chr1", "+", 1050, 1600, CHAIN_A)]
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000))])
    m = judge.compute_chain_metrics(iso, peaks, None, 50, read_ends=reads)
    assert m["5prime_altpair_51_200_truth_pairs"] == 0
    assert m["5prime_altpair_over_200_truth_pairs"] == 0


def test_minus_strand_over_200_pair_and_tie_to_the_5prime_most_end(tmp_path):
    """On -, the 5' end is the BED end; a read-count tie makes the 5'-most (highest) end main."""
    chain = ((4000, 4500),)
    exons_lo = [(3800, 4000), (4500, 5000)]
    exons_hi = [(3800, 4000), (4500, 5300)]
    iso = isoforms(tmp_path, [bed12("chr1", "-", exons_lo, score=4), bed12("chr1", "-", exons_hi, score=4)])
    peaks = points(("chr1", "-", 5000), ("chr1", "-", 5300))
    reads = [read("chr1", "-", 5000, 3800, chain), read("chr1", "-", 5300, 3800, chain)]
    m = judge.compute_chain_metrics(iso, peaks, None, 50, read_ends=reads)
    assert m["5prime_altpair_over_200_truth_pairs"] == 1
    assert m["5prime_altpair_over_200_pairs_hit"] == 1
    assert m["5prime_altpair_over_200_alt_ends"] == 1   # 5000 is the alternative end, 300 bp away
    assert m["5prime_altpair_over_200_precision"] == 1.0


def test_main_end_has_the_most_reads(tmp_path):
    """The alternative end is the one with fewer reads, here the one that misses."""
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000), score=1),
                              bed12("chr1", "+", exons_a(1300), score=9)])
    m = judge.compute_chain_metrics(iso, points(("chr1", "+", 1300)), None, 50)
    assert m["5prime_altpair_over_200_alt_ends"] == 1
    assert m["5prime_altpair_over_200_precision"] == 0.0


def test_main_end_tie_on_plus_goes_to_the_lowest_coordinate(tmp_path):
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000), score=2),
                              bed12("chr1", "+", exons_a(1300), score=2)])
    m = judge.compute_chain_metrics(iso, points(("chr1", "+", 1300)), None, 50)
    assert m["5prime_altpair_over_200_alt_ends_hit"] == 1   # 1300 is alternative, and hits


def test_main_end_reads_sum_over_isoforms_sharing_it(tmp_path):
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000, 1600), score=2),
                              bed12("chr1", "+", exons_a(1000, 1900), score=2),
                              bed12("chr1", "+", exons_a(1300), score=3)])
    m = judge.compute_chain_metrics(iso, points(("chr1", "+", 1300)), None, 50)
    assert m["5prime_altpair_over_200_alt_ends_hit"] == 1   # 1000 has 4 reads, so 1300 is alternative


def test_alternative_ends_within_50_bp_of_the_main_end_are_in_no_bin(tmp_path):
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000), score=3),
                              bed12("chr1", "+", exons_a(1040), score=1)])
    m = judge.compute_chain_metrics(iso, points(("chr1", "+", 1000)), None, 50)
    assert m["5prime_altpair_51_200_alt_ends"] == 0
    assert m["5prime_altpair_over_200_alt_ends"] == 0


def test_three_prime_pairs(tmp_path):
    iso = isoforms(tmp_path, [bed12("chr1", "+", exons_a(1000, 2000), score=3),
                              bed12("chr1", "+", exons_a(1000, 2100), score=1)])
    peaks3 = points(("chr1", "+", 2000), ("chr1", "+", 2100))
    reads = [read("chr1", "+", 1000, 2000, CHAIN_A), read("chr1", "+", 1000, 2100, CHAIN_A)]
    m = judge.compute_chain_metrics(iso, None, peaks3, 50, read_ends=reads)
    assert m["3prime_altpair_51_200_truth_pairs"] == 1
    assert m["3prime_altpair_51_200_recall"] == 1.0
    assert m["3prime_altpair_51_200_precision"] == 1.0
    assert m["5prime_altpair_51_200_recall"] is None


# ------------------------------------------------------------------ --read-ends input
def dump_rows():
    return [
        # aligned ends differ from the corrected ones, so reading the wrong columns fails
        ["r1", "chr1", "+", "1000", "1600", "1012", "1590", "1000", "1600", "1012", "1590", "1100-1500", "True"],
        ["r2", "chr1", "-", "3800", "5300", "3809", "5288", "5300", "3800", "5288", "3809", "4000-4500",
         "False"],
        ["r3", "chr1", ".", "10", "90", "10", "90", "", "", "", "", "20-40", "False"],
    ]


EXPECTED = [judge.ReadEnd("chr1", "+", ((1100, 1500),), 1000, 1600),
            judge.ReadEnd("chr1", "-", ((4000, 4500),), 5300, 3800)]


def write_tsv(path, header, rows, opener=open):
    with opener(path, "wt") as fh:
        for row in [header, *rows]:
            fh.write("\t".join(row) + "\n")


def test_read_ends_from_a_dump_directory(tmp_path):
    d = tmp_path / "genome_record_dump"
    d.mkdir()
    rows = dump_rows()
    write_tsv(d / "chr1_0_2000.tsv", DUMP_COLUMNS, rows[:1])
    write_tsv(d / "chr1_3000_6000.tsv", DUMP_COLUMNS, rows[1:])
    assert sorted(judge.read_read_ends(str(d))) == sorted(EXPECTED)


def test_read_ends_from_a_dump_tgz(tmp_path):
    d = tmp_path / "genome_record_dump"
    d.mkdir()
    write_tsv(d / "chr1_0_6000.tsv", DUMP_COLUMNS, dump_rows())
    tgz = tmp_path / "genome_record_dump.tgz"
    with tarfile.open(tgz, "w:gz") as tar:
        tar.add(d, arcname="genome_record_dump")
    assert sorted(judge.read_read_ends(str(tgz))) == sorted(EXPECTED)


def test_read_ends_from_a_joined_tsv_gz(tmp_path):
    path = tmp_path / "lib.genome_records.tsv.gz"
    write_tsv(path, DUMP_COLUMNS + ["isoform"], [r + ["FLT1"] for r in dump_rows()], opener=gzip.open)
    assert sorted(judge.read_read_ends(str(path))) == sorted(EXPECTED)


def test_read_ends_reject_a_dump_without_reads(tmp_path):
    """An empty dump, or a directory with no region TSVs, is an error, not zero reads."""
    d = tmp_path / "not_a_dump"
    d.mkdir()
    with pytest.raises(ValueError, match="no reads"):
        judge.read_read_ends(str(d))
    path = tmp_path / "empty.tsv"
    write_tsv(path, DUMP_COLUMNS, [])
    with pytest.raises(ValueError, match="no reads"):
        judge.read_read_ends(str(path))


def test_read_ends_reject_a_file_without_the_dump_columns(tmp_path):
    path = tmp_path / "x.tsv"
    write_tsv(path, ["chrom", "strand"], [["chr1", "+"]])
    with pytest.raises(ValueError, match="five_prime"):
        judge.read_read_ends(str(path))


# ------------------------------------------------------------------ the CLI writes the new columns
def test_cli_appends_the_new_columns_after_the_old_ones(tmp_path):
    import csv
    import subprocess

    (tmp_path / "iso.bed").write_text(bed12("chr1", "+", exons_a(1000), score=3) + "\n"
                                      + bed12("chr1", "+", exons_a(1120), score=2) + "\n")
    (tmp_path / "ann.gtf").write_text(
        'chr1\tt\ttranscript\t1001\t1600\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
    (tmp_path / "cage.bed").write_text("#h\nchr1\t1000\t1001\tc1\t1\t+\nchr1\t1120\t1121\tc2\t1\t+\n")
    d = tmp_path / "dump"
    d.mkdir()
    write_tsv(d / "chr1_0_2000.tsv", DUMP_COLUMNS,
              [["r1", "chr1", "+", "1000", "1600", "1000", "1600", "1000", "1600", "1000", "1600",
                "1100-1500", "True"],
               ["r2", "chr1", "+", "1120", "1600", "1120", "1600", "1120", "1600", "1120", "1600",
                "1100-1500", "True"]])
    out = tmp_path / "out"
    r = subprocess.run([sys.executable, os.path.join(VENDORED, "ted_end_precision.py"),
                        "--isoforms-bed", str(tmp_path / "iso.bed"), "--gtf", str(tmp_path / "ann.gtf"),
                        "--peaks-5prime", str(tmp_path / "cage.bed"), "--read-ends", str(d),
                        "--outdir", str(out)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    with open(out / "precision_recall_summary.tsv") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        fh.seek(0)
        row = next(csv.DictReader(fh, delimiter="\t"))
    assert header[:len(judge.SUMMARY_FIELDS)] == judge.SUMMARY_FIELDS
    assert header[len(judge.SUMMARY_FIELDS):] == judge.CHAIN_FIELDS
    assert row["n_read_ends"] == "2"
    assert row["5prime_chain_recall_read_supported"] == "1.000000"
    assert row["5prime_altpair_51_200_recall"] == "1.000000"
    assert row["3prime_chain_precision"] == ""
