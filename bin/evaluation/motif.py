"""
Sequence motif analysis functions.

Provides utilities for extracting sequences around transcript ends
and computing position frequency matrices for motif analysis.
"""

import math
import os
import random
import statistics
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .utils import run, which, get_logger

logger = get_logger()


def extract_sequences_batch(
    genome_path: Path,
    regions: List[Tuple[str, int, str, int, int]],
) -> Dict[int, str]:
    """
    Extract multiple sequences in a single samtools call (50-100x faster than individual calls).

    Args:
        genome_path: Path to genome FASTA file
        regions: List of (chrom, pos, strand, upstream, downstream) tuples

    Returns:
        Dict mapping index -> sequence string (strand-aware)
    """
    if not regions or not which("samtools"):
        return {}

    # Create temp file with all regions
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.regions') as f:
        regions_file = f.name
        for idx, (chrom, pos, strand, upstream, downstream) in enumerate(regions):
            # Calculate genomic coordinates
            if strand == '+':
                start = max(0, pos - upstream)
                end = pos + downstream
            else:
                start = max(0, pos - downstream)
                end = pos + upstream

            # Write region (1-based for samtools)
            f.write(f"{chrom}:{start + 1}-{end}\n")

    try:
        # Single samtools call for all regions
        res = run(["samtools", "faidx", str(genome_path), "-r", regions_file])

        # Parse sequences from FASTA output
        sequences = {}
        current_seq = []
        current_idx = -1

        for line in res.stdout.strip().split('\n'):
            if line.startswith('>'):
                # Save previous sequence
                if current_idx >= 0 and current_seq:
                    seq = ''.join(current_seq).upper()
                    _, _, strand, _, _ = regions[current_idx]

                    # Reverse complement for minus strand
                    if strand == '-':
                        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
                        seq = ''.join(complement.get(b, 'N') for b in reversed(seq))

                    sequences[current_idx] = seq

                # Start new sequence
                current_idx += 1
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_idx >= 0 and current_seq:
            seq = ''.join(current_seq).upper()
            _, _, strand, _, _ = regions[current_idx]

            if strand == '-':
                complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
                seq = ''.join(complement.get(b, 'N') for b in reversed(seq))

            sequences[current_idx] = seq

        return sequences

    except Exception as e:
        logger.warning(f"Batch sequence extraction failed: {e}")
        return {}
    finally:
        # Clean up temp file
        try:
            os.unlink(regions_file)
        except Exception:
            pass


def extract_sequence_context(
    genome_path: Path,
    chrom: str,
    pos: int,
    strand: str,
    upstream: int = 50,
    downstream: int = 50,
) -> Optional[str]:
    """
    Extract sequence context around a position using samtools faidx.

    Note: For multiple sequences, use extract_sequences_batch() for 50-100x speedup.

    Args:
        genome_path: Path to genome FASTA file
        chrom: Chromosome name
        pos: Position (0-based)
        strand: Strand ('+' or '-')
        upstream: bp upstream (5' direction on strand)
        downstream: bp downstream (3' direction on strand)

    Returns:
        Sequence string (strand-aware), or None if extraction fails
    """
    # Use batch extraction for single sequence (simpler than duplicating logic)
    result = extract_sequences_batch(genome_path, [(chrom, pos, strand, upstream, downstream)])
    return result.get(0)


def compute_position_frequency_matrix(
    sequences: List[str],
) -> Optional[List[Dict[str, float]]]:
    """
    Compute position frequency matrix from aligned sequences.

    Args:
        sequences: List of equal-length sequences

    Returns:
        List of dicts, one per position, with nucleotide frequencies
    """
    if not sequences:
        return None

    seq_len = len(sequences[0])
    if not all(len(s) == seq_len for s in sequences):
        logger.warning("Sequences have unequal lengths; cannot compute PFM")
        return None

    pfm = []
    for i in range(seq_len):
        counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
        for seq in sequences:
            base = seq[i] if i < len(seq) else 'N'
            counts[base] = counts.get(base, 0) + 1

        total = sum(counts.values())
        freqs = {base: count / total for base, count in counts.items()}
        pfm.append(freqs)

    return pfm


def compute_information_content(pfm: List[Dict[str, float]]) -> List[float]:
    """Compute information content (bits) per position."""
    ic = []
    for pos_freqs in pfm:
        # Shannon entropy
        entropy = 0
        for base in ['A', 'T', 'G', 'C']:
            p = pos_freqs.get(base, 0)
            if p > 0:
                entropy -= p * math.log2(p)
        # Information content = max entropy (2 bits) - observed entropy
        ic.append(2.0 - entropy)
    return ic


def analyze_motifs_at_ends(
    read_ends: Dict[str, dict],
    genome_path: Path,
    end_type: str,
    upstream: int = 50,
    downstream: int = 50,
    max_sequences: int = 5000,
) -> dict:
    """
    Analyze sequence motifs around read transcript ends using batch extraction.

    Args:
        read_ends: Dict mapping read_id -> {'chrom', 'start', 'end', 'strand', 'tss', 'tts'}
        genome_path: Path to genome FASTA file
        end_type: 'tss' for 5' ends, 'tts' for 3' ends
        upstream: bp upstream of end to extract
        downstream: bp downstream of end to extract
        max_sequences: Maximum sequences to analyze (for performance)

    Returns:
        Dict with PFM, information content, and summary statistics
    """
    if not genome_path or not genome_path.exists():
        logger.warning("Genome file not available for motif analysis")
        return {}

    # Sample reads if too many
    read_list = list(read_ends.items())
    if len(read_list) > max_sequences:
        random.seed(42)
        read_list = random.sample(read_list, max_sequences)

    # Build regions list for batch extraction
    regions = []
    valid_reads = []
    for read_id, read_info in read_list:
        # Note: parse_reads_bed_ends returns lowercase keys ('chrom', 'strand', 'tss', 'tts')
        chrom = read_info['chrom']
        strand = read_info['strand']
        pos = read_info.get(end_type)

        if pos is not None:
            regions.append((chrom, pos, strand, upstream, downstream))
            valid_reads.append((read_id, chrom, pos, strand))

    if not regions:
        logger.warning(f"No valid positions for {end_type} motif analysis")
        return {}

    # Batch extract all sequences (50-100x faster than loop)
    logger.debug(f"Extracting {len(regions)} sequences for {end_type} motif analysis (batch mode)")
    sequence_dict = extract_sequences_batch(genome_path, regions)

    # Filter to expected length
    sequences = []
    positions_used = []
    expected_len = upstream + downstream
    for idx, (read_id, chrom, pos, strand) in enumerate(valid_reads):
        seq = sequence_dict.get(idx)
        if seq and len(seq) == expected_len:
            sequences.append(seq)
            positions_used.append((chrom, pos, strand))

    if not sequences:
        logger.warning(f"No sequences extracted for {end_type} motif analysis")
        return {}

    # Compute position frequency matrix
    pfm = compute_position_frequency_matrix(sequences)
    if not pfm:
        return {}

    ic = compute_information_content(pfm)

    # Find positions with high information content
    high_ic_positions = [(i - upstream, ic[i]) for i in range(len(ic)) if ic[i] > 0.5]

    result = {
        "n_sequences": len(sequences),
        "sequence_length": len(sequences[0]) if sequences else 0,
        "pfm": pfm,
        "information_content": ic,
        "mean_ic": statistics.mean(ic) if ic else 0,
        "high_ic_positions": high_ic_positions,
    }

    return result
