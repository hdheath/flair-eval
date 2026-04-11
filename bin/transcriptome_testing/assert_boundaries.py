#!/usr/bin/env python3
"""
assert_boundaries.py
====================
Compare FLAIR transcriptome output isoforms against ground-truth BED12 boundaries.

For each ground-truth isoform, find the best-matching FLAIR isoform by junction
chain overlap and compare TSS / TTS positions.

Reports:
  - Boundary precision: fraction of FLAIR isoforms that match a ground-truth
    isoform within TSS/TTS tolerance
  - Boundary recall: fraction of ground-truth isoforms recovered
  - Per-isoform TSS/TTS displacement (signed, bp)
  - Artifact detection: FLAIR isoforms with no ground-truth match may be artifacts

Exit codes:
  0  All assertions pass
  1  At least one assertion failed (details in JSON report)
"""

import argparse
import json
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional


@dataclass
class BedIsoform:
    chrom: str
    start: int
    end: int
    name: str
    strand: str
    exon_starts: list  # absolute coordinates
    exon_ends: list

    @property
    def tss(self) -> int:
        return self.start if self.strand == '+' else self.end

    @property
    def tts(self) -> int:
        return self.end if self.strand == '+' else self.start

    @property
    def junctions(self) -> list:
        """Return splice junctions as (end_of_exon_i, start_of_exon_i+1)."""
        juncs = []
        for i in range(len(self.exon_starts) - 1):
            juncs.append((self.exon_ends[i], self.exon_starts[i + 1]))
        return juncs


@dataclass
class MatchResult:
    truth_name: str
    flair_name: Optional[str]
    tss_displacement: Optional[int]  # flair - truth
    tts_displacement: Optional[int]
    junction_match: bool
    matched: bool


@dataclass
class Report:
    test_label: str
    n_truth: int = 0
    n_flair: int = 0
    n_matched: int = 0
    n_unmatched_truth: int = 0
    n_artifact_flair: int = 0
    precision: float = 0.0
    recall: float = 0.0
    matches: list = field(default_factory=list)
    artifacts: list = field(default_factory=list)
    passed: bool = True
    failure_reasons: list = field(default_factory=list)


def parse_bed12(path: Path) -> list:
    """Parse a BED12 file into BedIsoform objects."""
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
            starts = [int(x) for x in f[11].rstrip(',').split(',')[:n_exons]]
            abs_starts = [start + s for s in starts]
            abs_ends = [start + s + sz for s, sz in zip(starts, sizes)]
            isoforms.append(BedIsoform(
                chrom=chrom, start=start, end=end, name=name,
                strand=strand, exon_starts=abs_starts, exon_ends=abs_ends,
            ))
    return isoforms


def junction_overlap(truth: BedIsoform, flair: BedIsoform) -> float:
    """Jaccard similarity of junction sets."""
    if truth.chrom != flair.chrom or truth.strand != flair.strand:
        return 0.0
    t_juncs = set(truth.junctions)
    f_juncs = set(flair.junctions)
    if not t_juncs and not f_juncs:
        # Single-exon: check overlap region
        overlap = min(truth.end, flair.end) - max(truth.start, flair.start)
        span = max(truth.end, flair.end) - min(truth.start, flair.start)
        return max(0, overlap / span) if span > 0 else 0.0
    if not t_juncs or not f_juncs:
        return 0.0
    intersection = t_juncs & f_juncs
    union = t_juncs | f_juncs
    return len(intersection) / len(union)


def match_isoforms(
    truth_isos: list,
    flair_isos: list,
    tss_tol: int,
    tts_tol: int,
) -> Report:
    """Match ground-truth isoforms to FLAIR output."""
    report = Report(test_label='', n_truth=len(truth_isos), n_flair=len(flair_isos))

    used_flair = set()
    matches = []

    for t in truth_isos:
        best_score = 0.0
        best_flair = None
        for j, f in enumerate(flair_isos):
            if j in used_flair:
                continue
            score = junction_overlap(t, f)
            if score > best_score:
                best_score = score
                best_flair = (j, f)

        if best_flair is not None and best_score > 0.0:
            j, f = best_flair
            used_flair.add(j)
            tss_disp = f.tss - t.tss
            tts_disp = f.tts - t.tts
            junc_match = best_score >= 1.0
            boundary_ok = abs(tss_disp) <= tss_tol and abs(tts_disp) <= tts_tol
            matches.append(MatchResult(
                truth_name=t.name, flair_name=f.name,
                tss_displacement=tss_disp, tts_displacement=tts_disp,
                junction_match=junc_match, matched=boundary_ok,
            ))
        else:
            matches.append(MatchResult(
                truth_name=t.name, flair_name=None,
                tss_displacement=None, tts_displacement=None,
                junction_match=False, matched=False,
            ))

    report.matches = [asdict(m) for m in matches]
    report.n_matched = sum(1 for m in matches if m.matched)
    report.n_unmatched_truth = sum(1 for m in matches if not m.matched)

    # Artifacts: FLAIR isoforms with no truth match
    artifact_indices = set(range(len(flair_isos))) - used_flair
    report.artifacts = [flair_isos[i].name for i in artifact_indices]
    report.n_artifact_flair = len(report.artifacts)

    # Precision / Recall
    report.precision = report.n_matched / report.n_flair if report.n_flair > 0 else 0.0
    report.recall = report.n_matched / report.n_truth if report.n_truth > 0 else 0.0

    return report


def run_assertions(report: Report, min_recall: float = 0.5) -> Report:
    """Apply pass/fail assertions to the report."""
    if report.n_truth == 0:
        report.failure_reasons.append('No ground-truth isoforms provided')
        report.passed = False
        return report

    if report.recall < min_recall:
        report.failure_reasons.append(
            f'Recall {report.recall:.3f} < minimum {min_recall}'
        )
        report.passed = False

    # Check for large TSS/TTS displacements in matched isoforms
    for m in report.matches:
        if m['matched'] and m['tss_displacement'] is not None:
            if abs(m['tss_displacement']) > 50:
                report.failure_reasons.append(
                    f"Large TSS displacement for {m['truth_name']}: "
                    f"{m['tss_displacement']} bp"
                )
            if abs(m['tts_displacement']) > 50:
                report.failure_reasons.append(
                    f"Large TTS displacement for {m['truth_name']}: "
                    f"{m['tts_displacement']} bp"
                )

    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--flair-bed', required=True, type=Path,
                        help='FLAIR isoforms.bed output')
    parser.add_argument('--truth-bed', required=True, type=Path,
                        help='Ground-truth BED12')
    parser.add_argument('--read-map', type=Path, default=None,
                        help='FLAIR read map (for read support stats)')
    parser.add_argument('--output', required=True, type=Path,
                        help='Output JSON report')
    parser.add_argument('--test-label', default='test',
                        help='Label for the test run')
    parser.add_argument('--tss-tolerance', type=int, default=100,
                        help='Max TSS displacement in bp to count as match')
    parser.add_argument('--tts-tolerance', type=int, default=100,
                        help='Max TTS displacement in bp to count as match')
    parser.add_argument('--min-recall', type=float, default=0.5,
                        help='Minimum recall to pass assertion')
    args = parser.parse_args()

    truth_isos = parse_bed12(args.truth_bed)
    flair_isos = parse_bed12(args.flair_bed)

    report = match_isoforms(truth_isos, flair_isos,
                            tss_tol=args.tss_tolerance,
                            tts_tol=args.tts_tolerance)
    report.test_label = args.test_label
    report = run_assertions(report, min_recall=args.min_recall)

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
        for m in report.matches:
            if m.get('flair_name') and m['flair_name'] in read_counts:
                m['read_support'] = read_counts[m['flair_name']]

    with open(args.output, 'w') as fh:
        json.dump(asdict(report) if hasattr(report, '__dataclass_fields__') else report.__dict__,
                  fh, indent=2, default=str)

    print(f"[{args.test_label}]  truth={report.n_truth}  flair={report.n_flair}  "
          f"matched={report.n_matched}  recall={report.recall:.3f}  "
          f"precision={report.precision:.3f}  "
          f"artifacts={report.n_artifact_flair}  "
          f"PASS={report.passed}")

    sys.exit(0 if report.passed else 1)


if __name__ == '__main__':
    main()
