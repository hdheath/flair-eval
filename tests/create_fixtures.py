#!/usr/bin/env python3
"""
Generate synthetic test fixtures for nf-test module tests.

Creates a minimal set of genomics files on a tiny synthetic chromosome
(chr_test, 50kb) with 5 genes / 10 transcripts / 200 reads.

All files are self-consistent: reads align to the genome, isoforms match
annotation coordinates, peaks match TSS/TTS positions, etc.

Usage:
    python tests/create_fixtures.py

Output (written to tests/data/):
    tiny_genome.fa          - 50kb synthetic FASTA
    tiny_genome.fa.fai      - samtools faidx index
    tiny_annotation.gtf     - 5 genes, 10 transcripts
    tiny_reads.bam          - 200 synthetic reads (spliced)
    tiny_reads.bam.bai      - BAM index
    tiny_reads.bed          - BED12 from reads
    tiny_cage_peaks.bed     - 8 CAGE peaks at TSS positions
    tiny_quantseq_peaks.bed - 8 QuantSeq peaks at TTS positions
    tiny_isoforms.bed       - 10 BED12 isoforms (FLAIR-format)
    tiny_isoforms.gtf       - Same 10 isoforms in GTF format
    tiny_read_map.txt       - Read-to-isoform assignments
    tiny_ref_tss.bed        - Reference TSS (for Evaluation input)
    tiny_ref_tts.bed        - Reference TTS (for Evaluation input)
    tiny_junctions.tab      - STAR-format SJ.out.tab
    tiny_evaluation.tsv     - Pre-computed evaluation TSV (for SummaryPlots)
    tiny_samplesheet.csv    - Samplesheet pointing to fixture files
    tiny_params.json        - Params JSON for test runs
"""

import os
import random
import json
import subprocess
import tempfile
from pathlib import Path

random.seed(42)

CHROM = "chr_test"
CHROM_LEN = 50000
OUT_DIR = Path(__file__).parent / "data"

# ─────────────────────────────────────────────────────────────────────
# Gene definitions: (gene_id, strand, transcripts)
# Each transcript: (tx_id, exon_list as [(start, end), ...])
# Coordinates are 0-based, half-open (BED convention)
# ─────────────────────────────────────────────────────────────────────

GENES = [
    ("GENE1", "+", [
        ("TX1a", [(1000, 1200), (2000, 2300), (3000, 3500)]),
        ("TX1b", [(1000, 1200), (3000, 3500)]),           # skipped exon
    ]),
    ("GENE2", "-", [
        ("TX2a", [(5000, 5400), (6000, 6200), (7000, 7600)]),
        ("TX2b", [(5000, 5400), (6000, 6200), (7000, 7300)]),  # alt 3' end
    ]),
    ("GENE3", "+", [
        ("TX3a", [(10000, 10300), (11000, 11500), (12000, 12800)]),
        ("TX3b", [(10000, 10300), (11000, 11500)]),
    ]),
    ("GENE4", "-", [
        ("TX4a", [(20000, 20500), (21000, 21400), (22000, 22700)]),
        ("TX4b", [(20000, 20500), (22000, 22700)]),
    ]),
    ("GENE5", "+", [
        ("TX5a", [(30000, 30200), (31000, 31600), (32000, 32400), (33000, 33800)]),
        ("TX5b", [(30000, 30200), (31000, 31600), (33000, 33800)]),
    ]),
]


def make_genome():
    """Create a synthetic genome FASTA with embedded biological motifs.

    Embeds:
    - TATAAA (TATA box) ~30bp upstream of each + strand TSS
    - AATAAA (canonical polyA signal) ~20bp upstream of each + strand TTS
    - A-rich stretch (AAAAAAAAAAAAAAAAAAAA) at internal positions as
      internal priming traps (positions 15000 and 25000)
    """
    seq = list(''.join(random.choice("ACGT") for _ in range(CHROM_LEN)))

    # Embed TATA box upstream of + strand TSS positions
    # and AATAAA upstream of + strand TTS positions
    for gene_id, strand, transcripts in GENES:
        for tx_id, exons in transcripts:
            if strand == '+':
                tss = exons[0][0]
                tts = exons[-1][1]
                # TATA box ~30bp upstream of TSS
                tata_pos = max(0, tss - 30)
                for i, base in enumerate("TATAAA"):
                    if tata_pos + i < CHROM_LEN:
                        seq[tata_pos + i] = base
                # Canonical polyA signal ~20bp upstream of TTS
                polya_pos = max(0, tts - 20)
                for i, base in enumerate("AATAAA"):
                    if polya_pos + i < CHROM_LEN:
                        seq[polya_pos + i] = base
            else:
                # For - strand: TATA on reverse strand near the 3' end of gene
                # and polyA signal near the 5' end
                tss = exons[-1][1]  # TSS is at the high-coord end
                tts = exons[0][0]   # TTS is at the low-coord end
                # Reverse complement of TATAAA = TTTATA
                tata_pos = min(CHROM_LEN - 6, tss + 25)
                for i, base in enumerate("TTTATA"):
                    if tata_pos + i < CHROM_LEN:
                        seq[tata_pos + i] = base
                # Reverse complement of AATAAA = TTTATT
                polya_pos = min(CHROM_LEN - 6, tts + 15)
                for i, base in enumerate("TTTATT"):
                    if polya_pos + i < CHROM_LEN:
                        seq[polya_pos + i] = base

    # Internal priming traps: A-rich stretches at positions away from TTS
    for trap_pos in [15000, 25000, 35000]:
        a_stretch = "A" * 20
        for i, base in enumerate(a_stretch):
            if trap_pos + i < CHROM_LEN:
                seq[trap_pos + i] = base

    seq = ''.join(seq)
    fa_path = OUT_DIR / "tiny_genome.fa"
    with open(fa_path, 'w') as f:
        f.write(f">{CHROM}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")
    # Index
    subprocess.run(["samtools", "faidx", str(fa_path)], check=True)
    return fa_path, seq


def make_annotation():
    """Create GTF annotation for all genes/transcripts."""
    gtf_path = OUT_DIR / "tiny_annotation.gtf"
    lines = []
    for gene_id, strand, transcripts in GENES:
        # Gene line
        all_starts = [e[0] for tx_id, exons in transcripts for e in exons]
        all_ends = [e[1] for tx_id, exons in transcripts for e in exons]
        gene_start = min(all_starts)
        gene_end = max(all_ends)
        lines.append(f'{CHROM}\ttest\tgene\t{gene_start+1}\t{gene_end}\t.\t{strand}\t.\t'
                      f'gene_id "{gene_id}"; gene_name "{gene_id}";')

        for tx_id, exons in transcripts:
            tx_start = exons[0][0]
            tx_end = exons[-1][1]
            lines.append(f'{CHROM}\ttest\ttranscript\t{tx_start+1}\t{tx_end}\t.\t{strand}\t.\t'
                          f'gene_id "{gene_id}"; transcript_id "{tx_id}"; gene_name "{gene_id}";')
            for i, (es, ee) in enumerate(exons):
                lines.append(f'{CHROM}\ttest\texon\t{es+1}\t{ee}\t.\t{strand}\t.\t'
                              f'gene_id "{gene_id}"; transcript_id "{tx_id}"; '
                              f'exon_number "{i+1}"; gene_name "{gene_id}";')
    with open(gtf_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return gtf_path


def make_isoforms_bed():
    """Create FLAIR-format BED12 for isoforms."""
    bed_path = OUT_DIR / "tiny_isoforms.bed"
    lines = []
    for gene_id, strand, transcripts in GENES:
        for tx_id, exons in transcripts:
            chrom_start = exons[0][0]
            chrom_end = exons[-1][1]
            block_count = len(exons)
            block_sizes = ','.join(str(e[1] - e[0]) for e in exons)
            block_starts = ','.join(str(e[0] - chrom_start) for e in exons)
            name = f"{tx_id}_{gene_id}"
            score = random.randint(10, 500)
            lines.append(f"{CHROM}\t{chrom_start}\t{chrom_end}\t{name}\t{score}\t{strand}\t"
                          f"{chrom_start}\t{chrom_end}\t0,0,0\t{block_count}\t{block_sizes}\t{block_starts}")
    with open(bed_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return bed_path


def make_isoforms_gtf():
    """Create GTF version of isoforms (for Bambu/IsoQuant evaluation testing)."""
    gtf_path = OUT_DIR / "tiny_isoforms.gtf"
    lines = []
    for gene_id, strand, transcripts in GENES:
        for tx_id, exons in transcripts:
            tx_start = exons[0][0]
            tx_end = exons[-1][1]
            lines.append(f'{CHROM}\tflair\ttranscript\t{tx_start+1}\t{tx_end}\t.\t{strand}\t.\t'
                          f'gene_id "{gene_id}"; transcript_id "{tx_id}";')
            for i, (es, ee) in enumerate(exons):
                lines.append(f'{CHROM}\tflair\texon\t{es+1}\t{ee}\t.\t{strand}\t.\t'
                              f'gene_id "{gene_id}"; transcript_id "{tx_id}"; exon_number "{i+1}";')
    with open(gtf_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return gtf_path


def make_read_map():
    """Create read-to-isoform map file."""
    map_path = OUT_DIR / "tiny_read_map.txt"
    lines = []
    read_id = 0
    for gene_id, strand, transcripts in GENES:
        for tx_id, exons in transcripts:
            # Assign 10 reads per isoform
            reads = [f"read_{read_id + i}" for i in range(10)]
            read_id += 10
            lines.append(f"{tx_id}_{gene_id}\t" + ",".join(reads))
    with open(map_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return map_path


def make_reads_sam(genome_seq):
    """Create synthetic SAM reads that align to isoform exons."""
    sam_path = OUT_DIR / "tiny_reads.sam"
    bam_path = OUT_DIR / "tiny_reads.bam"

    sam_lines = []
    sam_lines.append(f"@HD\tVN:1.6\tSO:coordinate")
    sam_lines.append(f"@SQ\tSN:{CHROM}\tLN:{CHROM_LEN}")

    read_id = 0
    for gene_id, strand, transcripts in GENES:
        for tx_id, exons in transcripts:
            for r in range(20):  # 20 reads per isoform = 200 total
                name = f"read_{read_id}"
                read_id += 1
                # Build CIGAR from exon structure
                flag = 0 if strand == '+' else 16
                pos = exons[0][0] + 1  # 1-based SAM

                cigar_parts = []
                for i, (es, ee) in enumerate(exons):
                    match_len = ee - es
                    cigar_parts.append(f"{match_len}M")
                    if i < len(exons) - 1:
                        intron_len = exons[i+1][0] - ee
                        cigar_parts.append(f"{intron_len}N")
                cigar = ''.join(cigar_parts)

                # Extract sequence from genome
                seq_parts = []
                for es, ee in exons:
                    seq_parts.append(genome_seq[es:ee])
                seq = ''.join(seq_parts)
                qual = 'I' * len(seq)

                sam_lines.append(f"{name}\t{flag}\t{CHROM}\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}")

    with open(sam_path, 'w') as f:
        f.write('\n'.join(sam_lines) + '\n')

    # Convert to sorted BAM + index
    subprocess.run(["samtools", "sort", "-o", str(bam_path), str(sam_path)], check=True)
    subprocess.run(["samtools", "index", str(bam_path)], check=True)
    os.unlink(sam_path)  # remove temp SAM
    return bam_path


def make_reads_bed():
    """Create BED12 from reads (mimics flair align output)."""
    bed_path = OUT_DIR / "tiny_reads.bed"
    lines = []
    read_id = 0
    for gene_id, strand, transcripts in GENES:
        for tx_id, exons in transcripts:
            for r in range(20):
                name = f"read_{read_id}"
                read_id += 1
                chrom_start = exons[0][0]
                chrom_end = exons[-1][1]
                block_count = len(exons)
                block_sizes = ','.join(str(e[1] - e[0]) for e in exons)
                block_starts = ','.join(str(e[0] - chrom_start) for e in exons)
                lines.append(f"{CHROM}\t{chrom_start}\t{chrom_end}\t{name}\t60\t{strand}\t"
                              f"{chrom_start}\t{chrom_end}\t0,0,0\t{block_count}\t{block_sizes}\t{block_starts}")
    with open(bed_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return bed_path


def make_peaks():
    """Create CAGE (TSS) and QuantSeq (TTS) peak BED files."""
    cage_path = OUT_DIR / "tiny_cage_peaks.bed"
    quantseq_path = OUT_DIR / "tiny_quantseq_peaks.bed"

    cage_lines = []
    quantseq_lines = []
    peak_id = 0

    for gene_id, strand, transcripts in GENES:
        # Collect unique TSS and TTS positions
        tss_positions = set()
        tts_positions = set()
        for tx_id, exons in transcripts:
            if strand == '+':
                tss_positions.add(exons[0][0])
                tts_positions.add(exons[-1][1])
            else:
                tss_positions.add(exons[-1][1])
                tts_positions.add(exons[0][0])

        for tss in sorted(tss_positions):
            peak_id += 1
            # Peaks are narrow — ±50bp around the position
            start = max(0, tss - 50)
            end = min(CHROM_LEN, tss + 50)
            cage_lines.append(f"{CHROM}\t{start}\t{end}\tcage_peak_{peak_id}\t100\t{strand}")

        for tts in sorted(tts_positions):
            peak_id += 1
            start = max(0, tts - 50)
            end = min(CHROM_LEN, tts + 50)
            quantseq_lines.append(f"{CHROM}\t{start}\t{end}\tquantseq_peak_{peak_id}\t100\t{strand}")

    with open(cage_path, 'w') as f:
        f.write('\n'.join(cage_lines) + '\n')
    with open(quantseq_path, 'w') as f:
        f.write('\n'.join(quantseq_lines) + '\n')

    return cage_path, quantseq_path


def make_ref_peaks():
    """Create reference TSS/TTS BED files (ground truth for TED)."""
    tss_path = OUT_DIR / "tiny_ref_tss.bed"
    tts_path = OUT_DIR / "tiny_ref_tts.bed"

    tss_lines = []
    tts_lines = []

    for gene_id, strand, transcripts in GENES:
        for tx_id, exons in transcripts:
            if strand == '+':
                tss = exons[0][0]
                tts = exons[-1][1]
            else:
                tss = exons[-1][1]
                tts = exons[0][0]
            tss_lines.append(f"{CHROM}\t{tss}\t{tss+1}\t{tx_id}_tss\t0\t{strand}")
            tts_lines.append(f"{CHROM}\t{tts}\t{tts+1}\t{tx_id}_tts\t0\t{strand}")

    with open(tss_path, 'w') as f:
        f.write('\n'.join(tss_lines) + '\n')
    with open(tts_path, 'w') as f:
        f.write('\n'.join(tts_lines) + '\n')

    return tss_path, tts_path


def make_junctions():
    """Create STAR-format SJ.out.tab from transcript splice junctions."""
    sj_path = OUT_DIR / "tiny_junctions.tab"
    lines = []

    for gene_id, strand, transcripts in GENES:
        for tx_id, exons in transcripts:
            for i in range(len(exons) - 1):
                intron_start = exons[i][1] + 1  # 1-based
                intron_end = exons[i+1][0]        # 1-based inclusive
                strand_code = 1 if strand == '+' else 2
                # STAR SJ.out.tab: chr, intron_start, intron_end, strand, intron_motif, annotated, unique_reads, multi_reads, max_overhang
                lines.append(f"{CHROM}\t{intron_start}\t{intron_end}\t{strand_code}\t1\t1\t20\t0\t50")

    # Deduplicate (same junction from multiple transcripts)
    unique_lines = list(dict.fromkeys(lines))
    with open(sj_path, 'w') as f:
        f.write('\n'.join(unique_lines) + '\n')
    return sj_path


def make_evaluation_tsv():
    """Create a pre-computed evaluation TSV for SummaryPlots testing.

    Uses the real column names from synthesize_evaluations.py output:
    test_name, dataset, align_mode, partition_mode, transcriptome_mode,
    5prime_precision, 5prime_recall, 5prime_f1, 3prime_precision, 3prime_recall, 3prime_f1,
    ref5prime_precision, ref5prime_recall, ref5prime_f1, ref3prime_precision, ref3prime_recall, ref3prime_f1,
    isoforms_observed, genes_observed, FSM, ISM, NIC, NNC, SEM, SEN,
    total_sjc, supported_sjc, total_se, supported_se, ...
    """
    tsv_path = OUT_DIR / "tiny_evaluation.tsv"
    header = [
        "test_name", "dataset", "align_mode", "partition_mode", "transcriptome_mode",
        "isoforms_observed", "genes_observed",
        "assigned_unique_read_ids", "assigned_primary_alignments",
        "assigned_supplementary_alignments", "assigned_total_alignments",
        "reads_per_isoform_mean", "reads_per_isoform_median",
        "reads_per_isoform_min", "reads_per_isoform_max",
        "input_primary_alignments", "input_supplementary_alignments", "input_total_alignments",
        "assignment_rate", "primary_alignment_utilization", "total_alignment_utilization",
        "5prime_precision", "5prime_recall", "5prime_f1",
        "3prime_precision", "3prime_recall", "3prime_f1",
        "ref5prime_precision", "ref5prime_recall", "ref5prime_f1",
        "ref3prime_precision", "ref3prime_recall", "ref3prime_f1",
        "total_read_regions", "found_regions", "genic_reads",
        "total_sjc", "supported_sjc", "subset_sjc",
        "total_se", "supported_se",
        "FSM", "ISM", "NIC", "NNC", "SEM", "SEN",
    ]
    rows = [
        ["test_set", "sample1", "pre-aligned", "chr-test", "default",
         "10", "5",
         "180", "180", "0", "180",
         "18.0", "15", "2", "50",
         "200", "0", "200",
         "0.9", "0.9", "0.9",
         "0.8", "0.7", "0.747",
         "0.75", "0.65", "0.696",
         "0.6", "0.5", "0.545",
         "0.55", "0.45", "0.495",
         "10", "8", "170",
         "15", "12", "3",
         "20", "18",
         "4", "3", "2", "1", "0", "0"],
        ["test_set", "sample1", "pre-aligned", "chr-test", "bambu_default",
         "12", "5",
         "175", "175", "0", "175",
         "14.5", "12", "1", "40",
         "200", "0", "200",
         "0.875", "0.875", "0.875",
         "0.75", "0.65", "0.696",
         "0.7", "0.6", "0.646",
         "0.55", "0.45", "0.495",
         "0.5", "0.4", "0.444",
         "10", "7", "165",
         "15", "11", "4",
         "20", "16",
         "5", "3", "2", "1", "1", "0"],
    ]
    with open(tsv_path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for row in rows:
            f.write('\t'.join(row) + '\n')
    return tsv_path


def make_samplesheet():
    """Create a test samplesheet pointing to fixture files."""
    csv_path = OUT_DIR / "tiny_samplesheet.csv"
    data_dir = "${projectDir}/tests/data"
    row = ','.join([
        "sample1",
        f"{data_dir}/tiny_genome.fa",
        f"{data_dir}/tiny_annotation.gtf",
        "",  # reads (empty - using BAM)
        f"{data_dir}/tiny_reads.bam",
        f"{data_dir}/tiny_cage_peaks.bed",
        f"{data_dir}/tiny_quantseq_peaks.bed",
        f"{data_dir}/tiny_junctions.tab",
        "pacbio_cDNA",  # library_type
        "", "", "", ""  # signal bedgraphs (empty)
    ])
    header = "sample_id,genome,gtf,reads,bam,cage,quantseq,junction_tab,library_type,cage_signal_plus,cage_signal_minus,quantseq_signal_plus,quantseq_signal_minus"
    with open(csv_path, 'w') as f:
        f.write(header + '\n')
        f.write(row + '\n')
    return csv_path


def make_params_json():
    """Create a minimal params JSON for test runs."""
    params = {
        "align": {"pre-aligned": ""},
        "partition": {"chr-test": f"--region {CHROM}:1000-40000"},
        "transcriptome": {"default": ""},
        "bambu": {},
        "isoquant": {}
    }
    json_path = OUT_DIR / "tiny_params.json"
    with open(json_path, 'w') as f:
        json.dump(params, f, indent=4)
    return json_path


def make_cage_peak_reason_tsv():
    """Create a sample CAGE peak reason TSV for PeakReasonHeatmap testing."""
    tsv_path = OUT_DIR / "tiny_cage_peak_reasons.tsv"
    header = ["peak_id", "chrom", "start", "end", "strand", "reason", "distance", "isoform_id"]
    rows = [
        ["cage_peak_1", CHROM, "950", "1050", "+", "matched", "5", "TX1a_GENE1"],
        ["cage_peak_2", CHROM, "4950", "5050", "-", "no_isoform_nearby", "999", ""],
        ["cage_peak_3", CHROM, "9950", "10050", "+", "matched", "10", "TX3a_GENE3"],
        ["cage_peak_4", CHROM, "19950", "20050", "-", "wrong_strand", "15", ""],
    ]
    with open(tsv_path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for row in rows:
            f.write('\t'.join(row) + '\n')
    return tsv_path


def make_quantseq_peak_reason_tsv():
    """Create a sample QuantSeq peak reason TSV for PeakReasonHeatmap testing."""
    tsv_path = OUT_DIR / "tiny_quantseq_peak_reasons.tsv"
    header = ["peak_id", "chrom", "start", "end", "strand", "reason", "distance", "isoform_id"]
    rows = [
        ["quantseq_peak_1", CHROM, "3450", "3550", "+", "matched", "8", "TX1a_GENE1"],
        ["quantseq_peak_2", CHROM, "7250", "7350", "-", "matched", "12", "TX2b_GENE2"],
        ["quantseq_peak_3", CHROM, "12750", "12850", "+", "no_isoform_nearby", "999", ""],
    ]
    with open(tsv_path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for row in rows:
            f.write('\t'.join(row) + '\n')
    return tsv_path


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Creating test fixtures in {OUT_DIR}/")

    fa_path, genome_seq = make_genome()
    print(f"  ✓ {fa_path.name} ({CHROM_LEN} bp)")

    gtf_path = make_annotation()
    n_tx = sum(len(txs) for _, _, txs in GENES)
    print(f"  ✓ {gtf_path.name} ({len(GENES)} genes, {n_tx} transcripts)")

    bed_path = make_isoforms_bed()
    print(f"  ✓ {bed_path.name}")

    gtf_iso_path = make_isoforms_gtf()
    print(f"  ✓ {gtf_iso_path.name}")

    map_path = make_read_map()
    print(f"  ✓ {map_path.name}")

    bam_path = make_reads_sam(genome_seq)
    print(f"  ✓ {bam_path.name} (200 reads)")

    reads_bed_path = make_reads_bed()
    print(f"  ✓ {reads_bed_path.name}")

    cage_path, quantseq_path = make_peaks()
    print(f"  ✓ {cage_path.name}, {quantseq_path.name}")

    tss_path, tts_path = make_ref_peaks()
    print(f"  ✓ {tss_path.name}, {tts_path.name}")

    sj_path = make_junctions()
    print(f"  ✓ {sj_path.name}")

    eval_path = make_evaluation_tsv()
    print(f"  ✓ {eval_path.name}")

    cage_reason_path = make_cage_peak_reason_tsv()
    quantseq_reason_path = make_quantseq_peak_reason_tsv()
    print(f"  ✓ {cage_reason_path.name}, {quantseq_reason_path.name}")

    csv_path = make_samplesheet()
    print(f"  ✓ {csv_path.name}")

    json_path = make_params_json()
    print(f"  ✓ {json_path.name}")

    print(f"\nDone! {len(list(OUT_DIR.iterdir()))} files created.")


if __name__ == '__main__':
    main()
