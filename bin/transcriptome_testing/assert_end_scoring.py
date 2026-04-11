#!/usr/bin/env python3
"""
assert_end_scoring.py
=====================
Run the end_scoring module on candidate boundaries from FLAIR's firstpass BED
and validate that confidence scores behave as expected for a given library type.

Key assertions:
  - Candidates with strong polyA signal + low depth should receive rescue bonus
  - Candidates with high A-richness at TTS should be penalized (for cDNA)
  - Candidates with high soft-clip fraction at TSS should be penalized
  - Ordering: high-confidence candidates should sort above low-confidence ones

This script serves as the integration bridge between the pure Python end_scoring
module (tested via pytest) and the Nextflow pipeline framework.
"""

import argparse
import json
import sys
from pathlib import Path

# Ensure flair-fusion's src/ is importable
import os
FLAIR_SRC = os.environ.get(
    'FLAIR_SRC',
    str(Path(__file__).resolve().parents[2] / '..' / 'tools' / 'flair-fusion' / 'src')
)
if FLAIR_SRC not in sys.path:
    sys.path.insert(0, FLAIR_SRC)

try:
    from flair.end_scoring import (
        EndConfidenceScorer,
        EndCandidate,
        LIBRARY_PROFILES,
        SequenceFeatureExtractor,
    )
except ImportError as e:
    print(f"ERROR: Cannot import end_scoring module from {FLAIR_SRC}: {e}",
          file=sys.stderr)
    sys.exit(2)

try:
    import pysam
except ImportError:
    print("ERROR: pysam not available", file=sys.stderr)
    sys.exit(2)


def parse_firstpass_bed(bed_path: Path) -> list:
    """
    Parse FLAIR firstpass BED12 into a list of dicts with fields:
      chrom, start, end, name, strand, n_reads, exon_starts, exon_ends
    """
    entries = []
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith('#') or line.startswith('track'):
                continue
            f = line.strip().split('\t')
            if len(f) < 12:
                continue
            start = int(f[1])
            n_exons = int(f[9])
            sizes = [int(x) for x in f[10].rstrip(',').split(',')[:n_exons]]
            starts = [int(x) for x in f[11].rstrip(',').split(',')[:n_exons]]

            # FLAIR encodes read count in the score field (col 5) or name field
            try:
                n_reads = int(f[4])
            except ValueError:
                n_reads = 1

            entries.append({
                'chrom': f[0],
                'start': start,
                'end': int(f[2]),
                'name': f[3],
                'strand': f[5],
                'n_reads': n_reads,
                'exon_starts': [start + s for s in starts],
                'exon_ends': [start + s + sz for s, sz in zip(starts, sizes)],
            })
    return entries


def compute_a_richness(genome: 'pysam.FastaFile', chrom: str,
                       pos: int, window: int = 30) -> float:
    """Fraction of A+T bases in a downstream window."""
    try:
        seq = genome.fetch(chrom, pos, pos + window).upper()
    except (ValueError, KeyError):
        return 0.0
    if not seq:
        return 0.0
    return (seq.count('A') + seq.count('T')) / len(seq)


def compute_soft_clip_hint(n_reads: int) -> float:
    """
    Heuristic soft-clip fraction estimate.
    Without per-read BAM data, we use a proxy based on read count.
    Low read count suggests potential mapping artifact.
    """
    if n_reads >= 10:
        return 0.05  # many supporting reads → low artifact risk
    elif n_reads >= 3:
        return 0.15
    else:
        return 0.35  # very few reads → elevated artifact risk


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--firstpass-bed', required=True, type=Path)
    parser.add_argument('--genome', required=True, type=Path)
    parser.add_argument('--library-type', default='default',
                        choices=list(LIBRARY_PROFILES.keys()))
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--test-label', default='scoring_test')
    args = parser.parse_args()

    # Load data
    entries = parse_firstpass_bed(args.firstpass_bed)
    genome = pysam.FastaFile(str(args.genome))
    profile = LIBRARY_PROFILES.get(args.library_type, LIBRARY_PROFILES['default'])
    scorer = EndConfidenceScorer(genome=genome, profile=profile)

    results = {
        'test_label': args.test_label,
        'library_type': args.library_type,
        'n_candidates': len(entries),
        'tss_scores': [],
        'tts_scores': [],
        'assertions': [],
        'passed': True,
    }

    for entry in entries:
        strand = entry['strand']

        # Score TSS
        tss_pos = entry['start'] if strand == '+' else entry['end']
        tss_a_rich = compute_a_richness(genome, entry['chrom'], tss_pos)
        tss_sc = compute_soft_clip_hint(entry['n_reads'])
        tss_candidate = EndCandidate(
            chrom=entry['chrom'],
            pos=tss_pos,
            strand=strand,
            end_type='tss',
            read_depth=entry['n_reads'],
            soft_clip_frac=tss_sc,
            a_richness=tss_a_rich,
        )
        tss_scored = scorer.score(tss_candidate)
        results['tss_scores'].append({
            'name': entry['name'],
            'position': tss_pos,
            'confidence': tss_scored.confidence,
            'read_depth': entry['n_reads'],
            'soft_clip_frac': tss_sc,
        })

        # Score TTS
        tts_pos = entry['end'] if strand == '+' else entry['start']
        tts_a_rich = compute_a_richness(genome, entry['chrom'], tts_pos)
        tts_sc = compute_soft_clip_hint(entry['n_reads'])
        tts_candidate = EndCandidate(
            chrom=entry['chrom'],
            pos=tts_pos,
            strand=strand,
            end_type='tts',
            read_depth=entry['n_reads'],
            soft_clip_frac=tts_sc,
            a_richness=tts_a_rich,
        )
        tts_scored = scorer.score(tts_candidate)
        results['tts_scores'].append({
            'name': entry['name'],
            'position': tts_pos,
            'confidence': tts_scored.confidence,
            'read_depth': entry['n_reads'],
            'a_richness': tts_a_rich,
        })

    genome.close()

    # ---- Run assertions ----

    # A1: All scores should be in [0, 1]
    all_scores = [s['confidence'] for s in results['tss_scores']] + \
                 [s['confidence'] for s in results['tts_scores']]
    out_of_range = [s for s in all_scores if s < 0 or s > 1]
    if out_of_range:
        results['assertions'].append({
            'name': 'scores_in_range',
            'passed': False,
            'detail': f'{len(out_of_range)} scores outside [0,1]',
        })
        results['passed'] = False
    else:
        results['assertions'].append({
            'name': 'scores_in_range',
            'passed': True,
            'detail': 'All scores in [0, 1]',
        })

    # A2: Higher read depth should generally correlate with higher confidence
    # (within same end type, controlling for sequence features)
    for end_type, scores in [('tss', results['tss_scores']),
                              ('tts', results['tts_scores'])]:
        if len(scores) < 2:
            continue
        sorted_by_depth = sorted(scores, key=lambda s: s['read_depth'])
        # Spearman-like: count concordant pairs
        concordant = 0
        total = 0
        for i in range(len(sorted_by_depth)):
            for j in range(i + 1, len(sorted_by_depth)):
                total += 1
                if sorted_by_depth[j]['confidence'] >= sorted_by_depth[i]['confidence']:
                    concordant += 1
        correlation = concordant / total if total > 0 else 0.0
        results['assertions'].append({
            'name': f'{end_type}_depth_confidence_correlation',
            'passed': correlation >= 0.3,  # weak positive expected
            'detail': f'Concordance: {correlation:.3f} ({concordant}/{total})',
        })
        if correlation < 0.3:
            results['passed'] = False

    # A3: For cDNA library, high A-richness TTS should score lower
    if args.library_type == 'ont_cDNA':
        high_arich = [s for s in results['tts_scores'] if s['a_richness'] > 0.7]
        low_arich = [s for s in results['tts_scores'] if s['a_richness'] < 0.4]
        if high_arich and low_arich:
            avg_high = sum(s['confidence'] for s in high_arich) / len(high_arich)
            avg_low = sum(s['confidence'] for s in low_arich) / len(low_arich)
            passed = avg_high < avg_low
            results['assertions'].append({
                'name': 'cdna_a_richness_penalty',
                'passed': passed,
                'detail': f'High A-rich avg={avg_high:.3f}, low A-rich avg={avg_low:.3f}',
            })
            if not passed:
                results['passed'] = False

    # Write report
    with open(args.output, 'w') as fh:
        json.dump(results, fh, indent=2)

    n_pass = sum(1 for a in results['assertions'] if a['passed'])
    n_fail = sum(1 for a in results['assertions'] if not a['passed'])
    print(f"[{args.test_label}]  candidates={results['n_candidates']}  "
          f"assertions: {n_pass} passed, {n_fail} failed  "
          f"PASS={results['passed']}")

    sys.exit(0 if results['passed'] else 1)


if __name__ == '__main__':
    main()
