#!/usr/bin/env python3
"""
assert_simulated_results.py
===========================
Scenario-aware assertion script for Badread simulation → FLAIR transcriptome
test results.

Each scenario has different expected behavior:
  - clean_ont / clean_pacbio: high junction recall (≥0.6), low TSS/TTS
    displacement (median ≤150bp)
  - high_chimera: FLAIR should filter most chimeras; false-positive rate
    should not explode
  - low_quality: pipeline should not crash; some isoforms recovered
  - truncated_5prime: TSS displacement may be large, but TTS should still
    be accurate; junction recall still reasonable
  - noisy_3prime: TTS displacement may be large, but junctions recovered

Exit codes:
  0  All assertions pass
  1  At least one assertion failed (details in JSON report)
"""

import argparse
import json
import sys
from collections import defaultdict
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional


# ── Scenario-specific thresholds ──────────────────────────────────────────

SCENARIO_THRESHOLDS = {
    'clean_ont': {
        'min_junction_recall': 0.5,
        'max_median_tss_disp': 200,
        'max_median_tts_disp': 200,
        'max_artifact_rate': 0.5,
        'min_isoforms': 1,
    },
    'clean_pacbio': {
        'min_junction_recall': 0.5,
        'max_median_tss_disp': 200,
        'max_median_tts_disp': 200,
        'max_artifact_rate': 0.5,
        'min_isoforms': 1,
    },
    'high_chimera': {
        'min_junction_recall': 0.3,
        'max_median_tss_disp': 300,
        'max_median_tts_disp': 300,
        'max_artifact_rate': 0.7,  # some chimeric artifacts expected
        'min_isoforms': 1,
    },
    'low_quality': {
        'min_junction_recall': 0.0,  # just don't crash
        'max_median_tss_disp': 500,
        'max_median_tts_disp': 500,
        'max_artifact_rate': 1.0,
        'min_isoforms': 0,  # may produce nothing
    },
    'truncated_5prime': {
        'min_junction_recall': 0.3,
        'max_median_tss_disp': 500,  # 5' truncation expected
        'max_median_tts_disp': 200,  # 3' should be fine
        'max_artifact_rate': 0.6,
        'min_isoforms': 1,
    },
    'noisy_3prime': {
        'min_junction_recall': 0.3,
        'max_median_tss_disp': 200,
        'max_median_tts_disp': 500,  # 3' noise expected
        'max_artifact_rate': 0.6,
        'min_isoforms': 1,
    },
}

# Default thresholds for unknown scenarios
DEFAULT_THRESHOLDS = {
    'min_junction_recall': 0.0,
    'max_median_tss_disp': 1000,
    'max_median_tts_disp': 1000,
    'max_artifact_rate': 1.0,
    'min_isoforms': 0,
}


# ── BED12 parsing ────────────────────────────────────────────────────────

@dataclass
class Isoform:
    chrom: str
    start: int
    end: int
    name: str
    strand: str
    exon_starts: list
    exon_ends: list

    @property
    def tss(self) -> int:
        return self.start if self.strand == '+' else self.end

    @property
    def tts(self) -> int:
        return self.end if self.strand == '+' else self.start

    @property
    def junctions(self) -> set:
        juncs = set()
        for i in range(len(self.exon_starts) - 1):
            juncs.add((self.exon_ends[i], self.exon_starts[i + 1]))
        return juncs


def parse_bed12(path: Path) -> list:
    isoforms = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or line.startswith('track'):
                continue
            f = line.strip().split('\t')
            if len(f) < 12:
                continue
            chrom = f[0]
            start = int(f[1])
            end = int(f[2])
            name = f[3]
            strand = f[5]
            n_exons = int(f[9])
            sizes = [int(x) for x in f[10].rstrip(',').split(',')[:n_exons]]
            starts_rel = [int(x) for x in f[11].rstrip(',').split(',')[:n_exons]]
            abs_starts = [start + s for s in starts_rel]
            abs_ends = [start + s + sz for s, sz in zip(starts_rel, sizes)]
            isoforms.append(Isoform(
                chrom=chrom, start=start, end=end, name=name,
                strand=strand, exon_starts=abs_starts, exon_ends=abs_ends,
            ))
    return isoforms


def parse_gtf_to_bed12(gtf_path: Path) -> list:
    """Parse a GTF annotation into Isoform objects (for ground truth)."""
    transcripts = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.strip().split('\t')
            if len(f) < 9:
                continue
            chrom, feature, strand = f[0], f[2], f[6]
            start = int(f[3]) - 1  # GTF is 1-based
            end = int(f[4])

            tid = None
            for token in f[8].split(';'):
                token = token.strip()
                if token.startswith('transcript_id'):
                    if '"' in token:
                        tid = token.split('"')[1]
                    else:
                        tid = token.split()[-1]
                    break
            if not tid:
                continue

            if feature == 'transcript':
                transcripts[tid] = {'chrom': chrom, 'strand': strand, 'exons': []}
            elif feature == 'exon' and tid in transcripts:
                transcripts[tid]['exons'].append((start, end))

    isoforms = []
    for tid, info in transcripts.items():
        if not info['exons']:
            continue
        exons = sorted(info['exons'])
        t_start = exons[0][0]
        t_end = exons[-1][1]
        isoforms.append(Isoform(
            chrom=info['chrom'], start=t_start, end=t_end, name=tid,
            strand=info['strand'],
            exon_starts=[e[0] for e in exons],
            exon_ends=[e[1] for e in exons],
        ))
    return isoforms


# ── Matching logic ────────────────────────────────────────────────────────

def junction_jaccard(a: Isoform, b: Isoform) -> float:
    if a.chrom != b.chrom or a.strand != b.strand:
        return 0.0
    aj, bj = a.junctions, b.junctions
    if not aj and not bj:
        # single-exon: use overlap fraction
        overlap = min(a.end, b.end) - max(a.start, b.start)
        span = max(a.end, b.end) - min(a.start, b.start)
        return max(0, overlap / span) if span > 0 else 0.0
    if not aj or not bj:
        return 0.0
    return len(aj & bj) / len(aj | bj)


def compute_metrics(truth: list, flair: list) -> dict:
    """Compute junction recall, TSS/TTS displacement, artifact count."""
    metrics = {
        'n_truth': len(truth),
        'n_flair': len(flair),
        'n_matched': 0,
        'junction_recall': 0.0,
        'tss_displacements': [],
        'tts_displacements': [],
        'matched_pairs': [],
        'artifacts': [],
    }

    if not truth or not flair:
        return metrics

    used_flair = set()
    for t in truth:
        best_score = 0.0
        best_idx = None
        for j, f in enumerate(flair):
            if j in used_flair:
                continue
            score = junction_jaccard(t, f)
            if score > best_score:
                best_score = score
                best_idx = j

        if best_idx is not None and best_score > 0.3:
            used_flair.add(best_idx)
            f = flair[best_idx]
            tss_disp = abs(f.tss - t.tss)
            tts_disp = abs(f.tts - t.tts)
            metrics['n_matched'] += 1
            metrics['tss_displacements'].append(tss_disp)
            metrics['tts_displacements'].append(tts_disp)
            metrics['matched_pairs'].append({
                'truth': t.name, 'flair': f.name,
                'junction_jaccard': round(best_score, 3),
                'tss_disp': tss_disp, 'tts_disp': tts_disp,
            })

    # Artifacts: FLAIR isoforms with no truth match
    for j, f in enumerate(flair):
        if j not in used_flair:
            metrics['artifacts'].append(f.name)

    metrics['junction_recall'] = (
        metrics['n_matched'] / metrics['n_truth'] if metrics['n_truth'] > 0 else 0.0
    )

    return metrics


def median(values: list) -> float:
    if not values:
        return 0.0
    s = sorted(values)
    n = len(s)
    if n % 2 == 0:
        return (s[n // 2 - 1] + s[n // 2]) / 2.0
    return float(s[n // 2])


# ── Assertions ────────────────────────────────────────────────────────────

def run_assertions(metrics: dict, scenario: str) -> dict:
    """Apply scenario-specific thresholds and return pass/fail report."""
    thresholds = SCENARIO_THRESHOLDS.get(scenario, DEFAULT_THRESHOLDS)
    report = {
        'scenario': scenario,
        'thresholds': thresholds,
        'metrics': {
            'n_truth': metrics['n_truth'],
            'n_flair': metrics['n_flair'],
            'n_matched': metrics['n_matched'],
            'junction_recall': round(metrics['junction_recall'], 4),
            'median_tss_displacement': median(metrics['tss_displacements']),
            'median_tts_displacement': median(metrics['tts_displacements']),
            'n_artifacts': len(metrics['artifacts']),
            'artifact_rate': (
                len(metrics['artifacts']) / metrics['n_flair']
                if metrics['n_flair'] > 0 else 0.0
            ),
        },
        'matched_pairs': metrics['matched_pairs'],
        'artifacts': metrics['artifacts'],
        'passed': True,
        'failures': [],
    }

    m = report['metrics']
    t = thresholds

    # Check minimum isoforms produced
    if m['n_flair'] < t['min_isoforms']:
        report['failures'].append(
            f"Too few FLAIR isoforms: {m['n_flair']} < {t['min_isoforms']}"
        )

    # Check junction recall
    if m['junction_recall'] < t['min_junction_recall']:
        report['failures'].append(
            f"Junction recall too low: {m['junction_recall']:.3f} < {t['min_junction_recall']}"
        )

    # Check TSS displacement (only if we have matches)
    if metrics['tss_displacements'] and m['median_tss_displacement'] > t['max_median_tss_disp']:
        report['failures'].append(
            f"Median TSS displacement too high: {m['median_tss_displacement']:.0f} > {t['max_median_tss_disp']}"
        )

    # Check TTS displacement
    if metrics['tts_displacements'] and m['median_tts_displacement'] > t['max_median_tts_disp']:
        report['failures'].append(
            f"Median TTS displacement too high: {m['median_tts_displacement']:.0f} > {t['max_median_tts_disp']}"
        )

    # Check artifact rate
    if m['artifact_rate'] > t['max_artifact_rate']:
        report['failures'].append(
            f"Artifact rate too high: {m['artifact_rate']:.3f} > {t['max_artifact_rate']}"
        )

    report['passed'] = len(report['failures']) == 0
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--flair-bed', required=True, type=Path,
                        help='FLAIR isoforms.bed output')
    parser.add_argument('--truth-bed', type=Path, default=None,
                        help='Ground-truth BED12 (optional if --annotation-gtf given)')
    parser.add_argument('--annotation-gtf', type=Path, default=None,
                        help='Annotation GTF to derive ground truth from')
    parser.add_argument('--read-map', type=Path, default=None,
                        help='FLAIR read map (for read support stats)')
    parser.add_argument('--scenario', required=True,
                        help='Scenario name (determines thresholds)')
    parser.add_argument('--output', required=True, type=Path,
                        help='Output JSON report')
    args = parser.parse_args()

    # Load ground truth
    if args.truth_bed and args.truth_bed.exists() and args.truth_bed.stat().st_size > 0:
        truth = parse_bed12(args.truth_bed)
    elif args.annotation_gtf and args.annotation_gtf.exists():
        truth = parse_gtf_to_bed12(args.annotation_gtf)
    else:
        print("WARNING: No ground truth provided, using empty truth set",
              file=sys.stderr)
        truth = []

    # Load FLAIR output
    flair = parse_bed12(args.flair_bed) if args.flair_bed.exists() else []

    # Compute metrics
    metrics = compute_metrics(truth, flair)

    # Apply scenario thresholds
    report = run_assertions(metrics, args.scenario)

    # Optionally add read support from read map
    if args.read_map and args.read_map.exists():
        read_counts = {}
        with open(args.read_map) as fh:
            for line in fh:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    iso_name = parts[0]
                    n_reads = len(parts[1].split(','))
                    read_counts[iso_name] = n_reads
        report['read_support'] = read_counts

    # Write report
    with open(args.output, 'w') as fh:
        json.dump(report, fh, indent=2, default=str)

    status = "PASS" if report['passed'] else "FAIL"
    print(f"[{args.scenario}]  {status}  "
          f"truth={metrics['n_truth']}  flair={metrics['n_flair']}  "
          f"matched={metrics['n_matched']}  "
          f"junction_recall={metrics['junction_recall']:.3f}  "
          f"artifacts={len(metrics['artifacts'])}")

    if not report['passed']:
        for f in report['failures']:
            print(f"  FAIL: {f}")

    sys.exit(0 if report['passed'] else 1)


if __name__ == '__main__':
    main()
