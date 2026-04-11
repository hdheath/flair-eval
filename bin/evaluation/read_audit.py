#!/usr/bin/env python3
"""TED Read Auditor — per-read classification through FLAIR transcriptome assembly.

Compares each read's original BAM alignment position to its assigned
isoform's position to classify the type of end adjustment applied by TED
(or any other collapse mode).

Classifications
---------------
kept              Both ends within tolerance of the assigned isoform ends.
reassigned_5p     5' end moved beyond tolerance; 3' end unchanged.
reassigned_3p     3' end moved beyond tolerance; 5' end unchanged.
reassigned_both   Both ends moved beyond tolerance.
unassigned        Read not assigned to any isoform.

Outputs
-------
{output}.read_audit.bed   BED9 with itemRgb for IGV (colour = classification).
{output}.read_audit.tsv   Per-read detail table for downstream summary plots.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path

import pysam

# ── Classification categories & IGV colours (R,G,B) ────────────────────────

CLASSIFICATIONS = {
    'kept':             '0,128,0',      # Green
    'reassigned_5p':    '0,0,255',      # Blue
    'reassigned_3p':    '255,165,0',    # Orange
    'reassigned_both':  '255,0,0',      # Red
    'unassigned':       '150,150,150',  # Grey
}

# Stable ordering for summaries / plots
CLASS_ORDER = ['kept', 'reassigned_5p', 'reassigned_3p', 'reassigned_both', 'unassigned']


# ── Helpers ─────────────────────────────────────────────────────────────────

def parse_read_map(map_path: str) -> dict[str, str]:
    """Parse ``isoform.read.map.txt`` → {read_name: isoform_name}."""
    read_to_iso: dict[str, str] = {}
    with open(map_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or '\t' not in line:
                continue
            iso, reads_str = line.split('\t', 1)
            for r in reads_str.split(','):
                r = r.strip()
                if r:
                    read_to_iso[r] = iso
    return read_to_iso


def parse_isoforms_bed(bed_path: str) -> dict[str, tuple]:
    """Parse ``isoforms.bed`` → {isoform_name: (chrom, start, end, strand)}."""
    iso_info: dict[str, tuple] = {}
    with open(bed_path) as fh:
        for line in fh:
            cols = line.strip().split('\t')
            if len(cols) < 6:
                continue
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            name = cols[3]
            strand = cols[5]
            iso_info[name] = (chrom, start, end, strand)
    return iso_info


def _cigar_to_blocks(ref_start: int, cigar_tuples: list) -> list[tuple[int, int]]:
    """Convert CIGAR tuples to exon blocks (list of (block_start, block_size)).

    CIGAR ops that consume reference: M(0), D(2), N(3), =(7), X(8).
    N (skip/intron) splits blocks; M/D/=/X extend the current block.
    """
    _REF_CONSUMERS = {0, 2, 7, 8}  # M, D, =, X
    blocks: list[tuple[int, int]] = []
    pos = ref_start
    block_start = ref_start
    block_size = 0

    for op, length in cigar_tuples:
        if op == 3:  # N — intron / skip
            if block_size > 0:
                blocks.append((block_start, block_size))
            pos += length
            block_start = pos
            block_size = 0
        elif op in _REF_CONSUMERS:
            block_size += length
            pos += length
        # S(4), I(1), H(5), P(6) don't consume reference

    if block_size > 0:
        blocks.append((block_start, block_size))

    # Fallback: if CIGAR produced no blocks, make one spanning block
    if not blocks:
        blocks = [(ref_start, 1)]

    return blocks


def extract_bam_reads(bam_path: str, region: str | None = None) -> dict[str, tuple]:
    """Extract primary-alignment positions and exon blocks from a BAM file.

    Returns {read_name: (chrom, start, end, strand, blocks)}
    where blocks is a list of (block_start, block_size) tuples.
    """
    reads: dict[str, tuple] = {}
    bam = pysam.AlignmentFile(bam_path, 'rb')

    iterator = bam.fetch(region=region) if region else bam.fetch()
    for aln in iterator:
        if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
            continue
        strand = '-' if aln.is_reverse else '+'
        blocks = _cigar_to_blocks(aln.reference_start, aln.cigartuples or [])
        reads[aln.query_name] = (
            aln.reference_name,
            aln.reference_start,
            aln.reference_end,
            strand,
            blocks,
        )

    bam.close()
    return reads


def classify_read(
    orig_start: int, orig_end: int,
    iso_start: int, iso_end: int, iso_strand: str,
    tolerance: int,
) -> tuple[str, int, int]:
    """Classify a read's end-adjustment type.

    *Biological* direction is derived from the isoform strand:
      + strand → left = 5', right = 3'
      − strand → left = 3', right = 5'

    Returns (classification, delta_5p, delta_3p).
    Delta signs: positive = read end was downstream of isoform end.
    """
    left_delta = orig_start - iso_start
    right_delta = orig_end - iso_end

    left_changed = abs(left_delta) > tolerance
    right_changed = abs(right_delta) > tolerance

    if iso_strand == '-':
        fivep_changed = right_changed
        threep_changed = left_changed
        delta_5p = -right_delta   # higher coord = 5' on minus strand
        delta_3p = -left_delta
    else:
        fivep_changed = left_changed
        threep_changed = right_changed
        delta_5p = left_delta
        delta_3p = right_delta

    if fivep_changed and threep_changed:
        return 'reassigned_both', delta_5p, delta_3p
    elif fivep_changed:
        return 'reassigned_5p', delta_5p, delta_3p
    elif threep_changed:
        return 'reassigned_3p', delta_5p, delta_3p
    else:
        return 'kept', delta_5p, delta_3p


# ── Main ────────────────────────────────────────────────────────────────────

def run_audit(
    bam_path: str,
    read_map_path: str,
    isoforms_bed_path: str,
    output_prefix: str,
    tolerance: int = 50,
    region: str | None = None,
) -> dict[str, int]:
    """Run the full read audit and write output files.

    Returns classification counts dict.
    """
    log = logging.getLogger(__name__)

    # 1. Load data
    log.info('Parsing read map: %s', read_map_path)
    read_to_iso = parse_read_map(read_map_path)
    log.info('  %d reads assigned to isoforms', len(read_to_iso))

    log.info('Parsing isoforms BED: %s', isoforms_bed_path)
    iso_info = parse_isoforms_bed(isoforms_bed_path)
    log.info('  %d isoforms', len(iso_info))

    log.info('Extracting BAM reads: %s', bam_path)
    bam_reads = extract_bam_reads(bam_path, region=region)
    log.info('  %d primary alignments', len(bam_reads))

    # 2. Classify & write
    counts: dict[str, int] = defaultdict(int)

    bed_path = output_prefix + '.read_audit.bed'
    tsv_path = output_prefix + '.read_audit.tsv'

    log.info('Classifying reads (tolerance=%d bp) ...', tolerance)

    with open(bed_path, 'w') as bed_fh, open(tsv_path, 'w') as tsv_fh:
        # BED header for IGV colour support
        bed_fh.write('track name="read_audit" itemRgb="On"\n')

        # TSV header
        tsv_fh.write('\t'.join([
            'read_name', 'chrom', 'orig_start', 'orig_end', 'strand',
            'isoform', 'iso_start', 'iso_end', 'iso_strand',
            'delta_5p', 'delta_3p', 'classification',
        ]) + '\n')

        for read_name in sorted(bam_reads):
            chrom, start, end, strand, blocks = bam_reads[read_name]

            if read_name in read_to_iso:
                iso_name = read_to_iso[read_name]
                if iso_name in iso_info:
                    iso_chrom, iso_start, iso_end, iso_strand = iso_info[iso_name]
                    cls, d5, d3 = classify_read(
                        start, end, iso_start, iso_end, iso_strand, tolerance)
                else:
                    cls, d5, d3 = 'unassigned', 0, 0
                    iso_name = iso_start = iso_end = iso_strand = '.'
            else:
                cls = 'unassigned'
                d5 = d3 = 0
                iso_name = '.'
                iso_start = iso_end = iso_strand = '.'

            counts[cls] += 1
            rgb = CLASSIFICATIONS[cls]

            # BED12: chrom start end name score strand thickStart thickEnd
            #        itemRgb blockCount blockSizes blockStarts
            n_blocks = len(blocks)
            block_sizes = ','.join(str(bs) for _, bs in blocks)
            block_starts = ','.join(str(bstart - start) for bstart, _ in blocks)
            bed_fh.write(
                f'{chrom}\t{start}\t{end}\t{read_name}\t0\t{strand}'
                f'\t{start}\t{end}\t{rgb}'
                f'\t{n_blocks}\t{block_sizes}\t{block_starts}\n')

            # TSV detail
            tsv_fh.write(
                f'{read_name}\t{chrom}\t{start}\t{end}\t{strand}\t'
                f'{iso_name}\t{iso_start}\t{iso_end}\t{iso_strand}\t'
                f'{d5}\t{d3}\t{cls}\n')

    # 3. Summary
    total = sum(counts.values())
    log.info('Read audit complete: %d reads classified', total)
    for cls in CLASS_ORDER:
        c = counts.get(cls, 0)
        pct = 100.0 * c / total if total else 0.0
        log.info('  %-20s %8d  (%5.1f%%)', cls, c, pct)
    log.info('BED → %s', bed_path)
    log.info('TSV → %s', tsv_path)

    return dict(counts)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--bam', required=True,
                        help='Aligned BAM file (indexed)')
    parser.add_argument('--read-map', required=True,
                        help='isoform.read.map.txt from FLAIR')
    parser.add_argument('--isoforms-bed', required=True,
                        help='isoforms.bed from FLAIR')
    parser.add_argument('--output', required=True,
                        help='Output prefix (produces .read_audit.bed and .read_audit.tsv)')
    parser.add_argument('--tolerance', type=int, default=50,
                        help='BP tolerance for kept vs reassigned (default: 50)')
    parser.add_argument('--region', default=None,
                        help='Restrict to genomic region (e.g. chr22)')
    parser.add_argument('--verbose', '-v', action='store_true')
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s  %(levelname)-8s  %(message)s',
    )

    run_audit(
        bam_path=args.bam,
        read_map_path=args.read_map,
        isoforms_bed_path=args.isoforms_bed,
        output_prefix=args.output,
        tolerance=args.tolerance,
        region=args.region,
    )


if __name__ == '__main__':
    main()
