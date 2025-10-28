# src/flair_test_suite/qc/qc_utils.py
# -------------------------------------------------
# Shared helpers used by QC collector modules.
# Requirements:
#  • Python standard library: subprocess, itertools, collections, multiprocessing
#  • pysam               # for reading BAM and FASTA files
# Summary:
#   Provides common utility functions for QC collectors, including counting lines,
#   computing percentages, iterating primary alignments, counting splice junctions,
#   and extracting splice junction motifs in parallel.
# Functions:
#   count_lines(path) -> int
#       Runs `wc -l` on a file to return its line count.
#   percent(numerator, denominator, digits=2) -> float
#       Computes (numerator/denominator)*100 rounded to `digits`.
#   iter_primary(bam, limit) -> Iterator[AlignedSegment]
#       Streams up to `limit` primary alignments from a BAM.
#   count_unique_junctions(bed_path, sample_limit=None) -> int
#       Counts unique (chrom, donor, acceptor, strand) junctions in a BED12.
#   count_splice_junction_motifs(...)
#       Counts 4-mer splice site motifs around junctions using multiprocessing.

from __future__ import annotations
from pathlib import Path
import subprocess
from typing import Dict, Iterator

import pysam
from collections import Counter, defaultdict
from multiprocessing import Pool, cpu_count
import itertools
import logging

__all__ = [
    "count_lines",
    "bed_is_empty",
    "percent",
    "parse_gtf_attributes",
    "iter_primary",
    "count_unique_junctions",
    "SAMPLE_LIMIT",
    "count_splice_junction_motifs",
]

# ------------------------------------------------------------------
# Simple helpers ----------------------------------------------------

def count_lines(path: Path) -> int:
    """
    Count the number of lines in a file using `wc -l`.
    Raises if the command fails.
    """
    return int(
        subprocess.check_output(["wc", "-l", str(path)], text=True)
        .split()[0]
    )


def bed_is_empty(path: Path) -> bool:
    """Return True if `path` is missing or has zero data lines."""
    try:
        if (not path.exists()) or (path.stat().st_size == 0):
            return True
        return count_lines(path) == 0
    except Exception as e:
        logging.warning(f"[qc_utils] Failed to probe BED {path}: {e}; treating as empty.")
        return True


def percent(numerator: int, denominator: int, digits: int = 2) -> float:
    """
    Compute percentage = 100 * numerator / denominator, rounded to `digits`.
    Returns 0.0 if denominator is zero.
    """
    if denominator == 0:
        return 0.0
    return round(100 * numerator / denominator, digits)


def parse_gtf_attributes(attr_col: str) -> Dict[str, str]:
    """Return a dictionary mapping GTF attribute keys to values.

    Parameters
    ----------
    attr_col : str
        Raw attribute column from a GTF line, e.g. 'gene_id "g1"; transcript_id "t1";'

    Returns
    -------
    Dict[str, str]
        Mapping of attribute names to unquoted string values.
    """
    attrs: Dict[str, str] = {}
    for chunk in attr_col.rstrip(";").split(";"):
        chunk = chunk.strip()
        if chunk and " " in chunk:
            k, v = chunk.split(" ", 1)
            attrs[k] = v.strip('"')
    return attrs

# ------------------------------------------------------------------
# Streaming helper for large BAMs ----------------------------------

# Cap on sampled alignments for QC metrics
SAMPLE_LIMIT: int = 2_000_000
# Flag mask to exclude secondary, supplementary, and unmapped reads
PRIMARY_ONLY_FLAGMASK = 0x904

def iter_primary(
    bam: pysam.AlignmentFile,
    limit: int | None = SAMPLE_LIMIT
) -> Iterator[pysam.AlignedSegment]:
    """
    Yield up to `limit` primary alignments (not secondary/supplementary/unmapped).
    """
    count = 0
    for aln in bam.fetch(until_eof=True):
        # Skip if any excluded flags present
        if aln.flag & PRIMARY_ONLY_FLAGMASK:
            continue
        yield aln
        count += 1
        if limit is not None and count >= limit:
            break

# ------------------------------------------------------------------
# Junction counting from BED12 -------------------------------------

def _parse_blocks(block_sizes: str, block_starts: str) -> tuple[list[int], list[int]]:
    """
    Parse BED12 blockSizes and blockStarts into integer lists.
    """
    sizes = list(map(int, block_sizes.rstrip(",").split(",")))
    starts = list(map(int, block_starts.rstrip(",").split(",")))
    return sizes, starts


def count_unique_junctions(
    bed_path: Path,
    sample_limit: int | None = None
) -> int:
    """
    Count unique splice junctions in a BED12 file.
    Junction defined by (chrom, donor, acceptor, strand).
    Optionally limit scanning to first `sample_limit` lines.
    """
    seen: set[tuple[str, int, int, str]] = set()
    with open(bed_path) as fh:
        for idx, line in enumerate(fh):
            if not line.strip() or line.startswith('#'):
                continue
            if sample_limit is not None and idx >= sample_limit:
                break
            cols = line.rstrip("\n").split("\t")
            # Ensure BED12 format
            if len(cols) < 12:
                continue
            chrom = cols[0]
            chrom_start = int(cols[1])
            strand = cols[5]
            block_count = int(cols[9])
            if block_count < 2:
                continue
            sizes, starts = _parse_blocks(cols[10], cols[11])
            # Extract donor and acceptor positions for each intron
            for i in range(block_count - 1):
                donor = chrom_start + starts[i] + sizes[i]
                acceptor = chrom_start + starts[i+1]
                seen.add((chrom, donor, acceptor, strand))
    return len(seen)

# ------------------------------------------------------------------
# Splice junction motif extraction --------------------------------

def revcomp(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    """
    trans = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(trans)[::-1]


def bed_reader(path: Path) -> Iterator[list[str]]:
    """
    Yield BED fields lines, skipping headers and empty lines.
    """
    with open(path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith(("track","browser","#")):
                continue
            yield ln.rstrip().split('\t')


def introns_from_bed(fields: list[str]) -> Iterator[tuple[str,int,int,str]]:
    """
    Yield (chrom, donor_end, acceptor_start, strand) for each intron in BED12 fields.
    """
    chrom = fields[0]
    chromStart = int(fields[1])
    strand = fields[5]
    sizes = list(map(int, fields[10].rstrip(',').split(',')))
    starts = list(map(int, fields[11].rstrip(',').split(',')))
    for i in range(len(sizes) - 1):
        exon1_end = chromStart + starts[i] + sizes[i]
        exon2_start = chromStart + starts[i+1]
        yield chrom, exon1_end, exon2_start, strand


def chunk(lst: list, n: int) -> list[list]:
    """
    Split list `lst` into `n` roughly equal chunks.
    """
    k, m = divmod(len(lst), n)
    return [
        lst[i*k + min(i, m):(i+1)*k + min(i+1, m)]
        for i in range(n)
    ]


def _sj_worker(args) -> Counter:
    """
    Worker function: extract 4-mer motifs around splice junctions for a chunk of BED entries.
    """
    bed_fields_list, fasta_path, don_ex, don_in, acc_in, acc_ex = args
    fasta = pysam.FastaFile(str(fasta_path))
    ctr = Counter()
    for fields in bed_fields_list:
        for chrom, s, e, strand in introns_from_bed(fields):
            if strand == '+':
                donor = fasta.fetch(chrom, s-don_ex, s+don_in)
                acceptor = fasta.fetch(chrom, e-acc_in, e+acc_ex)
            else:
                donor = revcomp(fasta.fetch(chrom, e-don_in, e+don_ex))
                acceptor = revcomp(fasta.fetch(chrom, s-acc_ex, s+acc_ex))
            motif = (donor[-2:] + acceptor[:2]).upper()
            ctr[(strand, motif)] += 1
    fasta.close()
    return ctr


def count_splice_junction_motifs(
    bed_path: Path,
    fasta_path: Path,
    max_workers: int | None = None,
    max_bed_lines: int = 100_000,
    don_ex: int = 2, don_in: int = 2,
    acc_in: int = 2, acc_ex: int = 2
) -> Counter:
    """
    Count 4-mer splice site motifs by reading up to `max_bed_lines` BED entries,
    grouping by chromosome, and using up to `max_workers` parallel processes.
    Returns Counter mapping (strand, motif) -> count.
    """
    if max_workers is None:
        max_workers = 4

    # Read BED lines
    if max_bed_lines is None:
        all_fields = list(bed_reader(bed_path))
    else:
        all_fields = list(itertools.islice(bed_reader(bed_path), max_bed_lines))
    if not all_fields:
        raise ValueError("No usable BED lines!")

    # Group entries by chromosome
    chrom_groups: dict[str, list] = defaultdict(list)
    for fields in all_fields:
        chrom_groups[fields[0]].append(fields)

    # Determine number of workers
    n_workers = min(max_workers, cpu_count())
    # Create chunks for parallel work
    if len(chrom_groups) > 1:
        # Split chromosomes across workers
        chroms = list(chrom_groups)
        chrom_chunks = chunk(chroms, min(n_workers, len(chroms)))
        bed_chunks = [
            sum((chrom_groups[c] for c in grp), [])
            for grp in chrom_chunks
        ]
    else:
        # Single chromosome: split entries evenly
        single_list = next(iter(chrom_groups.values()))
        bed_chunks = [c for c in chunk(single_list, n_workers) if c]

    # Prepare arguments and run workers
    args = [(chunk, fasta_path, don_ex, don_in, acc_in, acc_ex) for chunk in bed_chunks]
    with Pool(processes=len(bed_chunks)) as pool:
        results = pool.map(_sj_worker, args)

    # Merge partial results
    total = Counter()
    for res in results:
        total.update(res)
    return total
