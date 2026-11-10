#!/usr/bin/env python3
"""
evaluate_end_trust.py — Evaluate when to trust 5' and 3' ends.

For each FLAIR firstpass BED, scores every candidate boundary with
EndConfidenceScorer under multiple library profiles and alpha values.
Compares against annotation (GTF) to measure displacement.  Produces a
comprehensive TSV that downstream plotting can consume.

Metrics computed per candidate boundary:
  1. Annotation displacement: signed distance to nearest annotated TSS/TTS
  2. Confidence score: EndConfidenceScorer composite for each (profile, alpha)
  3. Sequence feature scores: seq_score, depth_score, tech_penalty
  4. Rescue flag: whether the boundary was rescued by sequence features
  5. Trust bucket: categorised as trusted / uncertain / untrusted

Output: tab-separated file with one row per boundary candidate.
"""

import argparse
import csv
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam

# Allow running from bin/ directory
sys.path.insert(0, str(Path(__file__).resolve().parent))

from flair.end_scoring import (
    EndCandidate,
    EndConfidenceScorer,
    LIBRARY_PROFILES,
)
from flair.end_scoring.features import SequenceFeatureExtractor

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# GTF parsing — extract annotated TSS/TTS positions
# ---------------------------------------------------------------------------

def parse_annotated_ends(gtf_path: str) -> Dict[str, Dict[str, List[int]]]:
    """
    Extract unique TSS and TTS positions from a GTF file.

    Returns
    -------
    dict  chrom -> {"tss": sorted list of positions, "tts": sorted list}
    """
    ends: Dict[str, Dict[str, set]] = defaultdict(lambda: {"tss": set(), "tts": set()})

    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "transcript":
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1   # GTF 1-based → 0-based
            end = int(parts[4])
            strand = parts[6]

            if strand == "+":
                ends[chrom]["tss"].add(start)
                ends[chrom]["tts"].add(end)
            else:
                ends[chrom]["tss"].add(end)
                ends[chrom]["tts"].add(start)

    # Sort for binary-search later
    result = {}
    for chrom in ends:
        result[chrom] = {
            "tss": sorted(ends[chrom]["tss"]),
            "tts": sorted(ends[chrom]["tts"]),
        }
    return result


def _nearest_distance(query: int, sorted_positions: List[int]) -> int:
    """Signed distance to nearest position in a sorted list."""
    if not sorted_positions:
        return 999999
    import bisect
    idx = bisect.bisect_left(sorted_positions, query)
    best = 999999
    for i in (idx - 1, idx):
        if 0 <= i < len(sorted_positions):
            d = query - sorted_positions[i]
            if abs(d) < abs(best):
                best = d
    return best


# ---------------------------------------------------------------------------
# BED parsing — extract candidate boundaries from firstpass BED
# ---------------------------------------------------------------------------

def parse_firstpass_bed(bed_path: str) -> List[dict]:
    """
    Parse a FLAIR firstpass BED (BED12-like).

    Returns list of dicts with fields:
      chrom, start, end, name, score (read support), strand
    """
    candidates = []
    with open(bed_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            candidates.append({
                "chrom": parts[0],
                "start": int(parts[1]),
                "end": int(parts[2]),
                "name": parts[3],
                "score": int(parts[4]) if parts[4].isdigit() else 0,
                "strand": parts[5] if len(parts) > 5 else "+",
            })
    return candidates


# ---------------------------------------------------------------------------
# Core evaluation
# ---------------------------------------------------------------------------

ALPHA_VALUES = [0.0, 0.25, 0.5, 0.75, 1.0]
PROFILE_NAMES = ["default", "ont_cDNA", "ont_dRNA", "pacbio_isoseq", "pacbio_masseq"]

TRUST_THRESHOLDS = {
    "trusted":   0.5,     # confidence >= 0.5
    "uncertain": 0.25,    # 0.25 <= confidence < 0.5
    # untrusted: < 0.25
}


def categorise_trust(confidence: float) -> str:
    if confidence >= TRUST_THRESHOLDS["trusted"]:
        return "trusted"
    elif confidence >= TRUST_THRESHOLDS["uncertain"]:
        return "uncertain"
    return "untrusted"


def evaluate_boundaries(
    bed_path: str,
    genome_path: str,
    gtf_path: str,
    output_path: str,
    profiles: Optional[List[str]] = None,
    alphas: Optional[List[float]] = None,
) -> str:
    """
    Score every candidate boundary under multiple profiles/alphas and
    measure displacement from annotation.  Writes per-boundary TSV.
    """
    if profiles is None:
        profiles = PROFILE_NAMES
    if alphas is None:
        alphas = ALPHA_VALUES

    log.info(f"Loading genome: {genome_path}")
    genome = pysam.FastaFile(genome_path)

    log.info(f"Parsing annotation: {gtf_path}")
    annotated_ends = parse_annotated_ends(gtf_path)

    log.info(f"Parsing firstpass BED: {bed_path}")
    candidates = parse_firstpass_bed(bed_path)
    log.info(f"  {len(candidates)} isoform candidates")

    # Build scorers for each profile
    scorers = {}
    for pname in profiles:
        profile = LIBRARY_PROFILES.get(pname, LIBRARY_PROFILES["default"])
        scorers[pname] = EndConfidenceScorer(genome=genome, profile=profile)

    # Output header
    fieldnames = [
        "chrom", "pos", "strand", "end_type",
        "read_depth", "isoform_name",
        "annot_displacement",
        "profile", "alpha",
        "seq_score", "depth_score", "tech_penalty", "confidence",
        "rescue_reason", "trust_category",
        "tss_trust", "tts_trust",
        "seq_weight", "depth_weight", "tech_weight",
    ]

    n_rows = 0
    with open(output_path, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for cand in candidates:
            chrom = cand["chrom"]
            strand = cand["strand"]
            depth = cand["score"]

            # Determine TSS and TTS positions
            if strand == "+":
                tss_pos, tts_pos = cand["start"], cand["end"]
            else:
                tss_pos, tts_pos = cand["end"], cand["start"]

            for end_type, pos in [("tss", tss_pos), ("tts", tts_pos)]:
                # Annotation displacement
                chrom_ends = annotated_ends.get(chrom, {}).get(end_type, [])
                displacement = _nearest_distance(pos, chrom_ends)

                for pname in profiles:
                    scorer = scorers[pname]
                    profile = LIBRARY_PROFILES[pname]

                    ec = EndCandidate(
                        chrom=chrom,
                        pos=pos,
                        strand=strand,
                        end_type=end_type,
                        read_depth=depth,
                    )

                    scored = scorer.score(ec)

                    for alpha in alphas:
                        # Recompute blended confidence with alpha
                        heuristic = depth  # simple read count
                        conf_blend = (
                            (1.0 - alpha) * (scored.depth_score)
                            + alpha * scored.confidence
                        )

                        row = {
                            "chrom": chrom,
                            "pos": pos,
                            "strand": strand,
                            "end_type": end_type,
                            "read_depth": depth,
                            "isoform_name": cand["name"],
                            "annot_displacement": displacement,
                            "profile": pname,
                            "alpha": f"{alpha:.2f}",
                            "seq_score": f"{scored.seq_score:.4f}",
                            "depth_score": f"{scored.depth_score:.4f}",
                            "tech_penalty": f"{scored.tech_penalty:.4f}",
                            "confidence": f"{scored.confidence:.4f}",
                            "rescue_reason": scored.rescue_reason or "",
                            "trust_category": categorise_trust(scored.confidence),
                            "tss_trust": f"{profile.tss_trust:.2f}",
                            "tts_trust": f"{profile.tts_trust:.2f}",
                            "seq_weight": f"{getattr(profile, f'seq_weight_{end_type}'):.2f}",
                            "depth_weight": f"{getattr(profile, f'depth_weight_{end_type}'):.2f}",
                            "tech_weight": f"{getattr(profile, f'tech_weight_{end_type}'):.2f}",
                        }
                        writer.writerow(row)
                        n_rows += 1

    log.info(f"Wrote {n_rows} rows to {output_path}")
    genome.close()
    return output_path


# ---------------------------------------------------------------------------
# Aggregate summary — one row per (profile, alpha, end_type)
# ---------------------------------------------------------------------------

def summarise(input_path: str, output_path: str) -> str:
    """
    Read the per-boundary TSV and produce an aggregate summary:
    one row per (profile, alpha, end_type) with statistics.
    """
    import statistics

    groups = defaultdict(lambda: {
        "displacements": [],
        "confidences": [],
        "rescued": 0,
        "n_trusted": 0,
        "n_uncertain": 0,
        "n_untrusted": 0,
        "total": 0,
    })

    with open(input_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key = (row["profile"], row["alpha"], row["end_type"])
            g = groups[key]
            g["displacements"].append(int(row["annot_displacement"]))
            g["confidences"].append(float(row["confidence"]))
            if row["rescue_reason"]:
                g["rescued"] += 1
            g[f"n_{row['trust_category']}"] += 1
            g["total"] += 1

    fieldnames = [
        "profile", "alpha", "end_type",
        "n_boundaries", "n_rescued",
        "mean_displacement", "median_displacement", "std_displacement",
        "mean_abs_displacement", "pct_within_50bp", "pct_within_100bp",
        "mean_confidence", "median_confidence",
        "pct_trusted", "pct_uncertain", "pct_untrusted",
    ]

    with open(output_path, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for (profile, alpha, end_type), g in sorted(groups.items()):
            disps = g["displacements"]
            confs = g["confidences"]
            abs_disps = [abs(d) for d in disps]
            n = g["total"]
            if n == 0:
                continue

            writer.writerow({
                "profile": profile,
                "alpha": alpha,
                "end_type": end_type,
                "n_boundaries": n,
                "n_rescued": g["rescued"],
                "mean_displacement": f"{statistics.mean(disps):.1f}",
                "median_displacement": f"{statistics.median(disps):.0f}",
                "std_displacement": f"{statistics.stdev(disps):.1f}" if n > 1 else "0",
                "mean_abs_displacement": f"{statistics.mean(abs_disps):.1f}",
                "pct_within_50bp": f"{100 * sum(1 for d in abs_disps if d <= 50) / n:.1f}",
                "pct_within_100bp": f"{100 * sum(1 for d in abs_disps if d <= 100) / n:.1f}",
                "mean_confidence": f"{statistics.mean(confs):.4f}",
                "median_confidence": f"{statistics.median(confs):.4f}",
                "pct_trusted": f"{100 * g['n_trusted'] / n:.1f}",
                "pct_uncertain": f"{100 * g['n_uncertain'] / n:.1f}",
                "pct_untrusted": f"{100 * g['n_untrusted'] / n:.1f}",
            })

    log.info(f"Wrote summary to {output_path}")
    return output_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate when to trust 5' and 3' ends under different "
                    "library profiles and scoring alphas.",
    )
    parser.add_argument("--bed", required=True,
                        help="FLAIR firstpass BED (or isoforms BED)")
    parser.add_argument("--genome", required=True,
                        help="Reference genome FASTA (indexed)")
    parser.add_argument("--gtf", required=True,
                        help="Annotation GTF for displacement measurement")
    parser.add_argument("--output", required=True,
                        help="Output TSV path (per-boundary)")
    parser.add_argument("--summary", default=None,
                        help="Optional summary TSV path (per profile/alpha/end_type)")
    parser.add_argument("--profiles", nargs="+", default=PROFILE_NAMES,
                        help="Library profiles to evaluate")
    parser.add_argument("--alphas", nargs="+", type=float, default=ALPHA_VALUES,
                        help="Alpha blend values to evaluate")

    args = parser.parse_args()

    evaluate_boundaries(
        bed_path=args.bed,
        genome_path=args.genome,
        gtf_path=args.gtf,
        output_path=args.output,
        profiles=args.profiles,
        alphas=args.alphas,
    )

    if args.summary:
        summarise(args.output, args.summary)


if __name__ == "__main__":
    main()
