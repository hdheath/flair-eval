#!/usr/bin/env python3
"""
firstpass_vs_final.py — Compare firstpass and final isoform BEDs.

Runs junction-chain-deduplicated end precision/recall on both the firstpass
(pre-TED) and final (post-TED) isoform BED files, then produces:
  1. A grouped bar chart comparing P/R/F1 between firstpass and final
  2. A summary TSV with all metrics side-by-side
"""

import argparse
import csv
import logging
import sys
from pathlib import Path

try:
    from signal_utils import parse_bed12, tss_tts
    from ted_end_precision import (
        parse_gtf_ends, parse_gtf_transcripts, parse_peaks_bed,
        compute_jc_deduplicated_precision_recall,
    )
    from pub_style import apply_rc, savefig as pub_savefig, style_ax, W2, PALETTE
except ImportError:
    from evaluation.signal_utils import parse_bed12, tss_tts
    from evaluation.ted_end_precision import (
        parse_gtf_ends, parse_gtf_transcripts, parse_peaks_bed,
        compute_jc_deduplicated_precision_recall,
    )
    from evaluation.pub_style import apply_rc, savefig as pub_savefig, style_ax, W2, PALETTE

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(levelname)-8s  %(message)s")
log = logging.getLogger(__name__)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    apply_rc()
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

GOLDEN_RATIO = 1.618


def run_precision_recall(bed_path, annotated_ends, annot_transcripts,
                         window, peaks_5prime, peaks_3prime):
    """Parse a BED file and compute JC-deduplicated P/R."""
    isoforms = parse_bed12(bed_path)
    log.info(f"  {len(isoforms)} isoforms from {Path(bed_path).name}")
    return compute_jc_deduplicated_precision_recall(
        isoforms, annotated_ends, annot_transcripts,
        window=window, peaks_5prime=peaks_5prime, peaks_3prime=peaks_3prime,
    )


def plot_comparison(firstpass_results, final_results, output_path, mode_label,
                    fig_width=None):
    """Grouped bar chart: firstpass vs final for precision, recall, F1."""
    if not HAS_MPL:
        log.warning("matplotlib not available, skipping plot")
        return

    if fig_width is None:
        fig_width = W2

    metrics = []
    labels = []
    for end in ("5prime", "3prime"):
        end_nice = "5'" if end == "5prime" else "3'"
        for metric, metric_nice in [("dedup_precision", "Precision"),
                                     ("recall", "Recall"),
                                     ("dedup_f1", "F1")]:
            key = f"{end}_{metric}"
            fp_val = firstpass_results.get(key)
            fn_val = final_results.get(key)
            metrics.append((fp_val or 0, fn_val or 0))
            labels.append(f"{end_nice} {metric_nice}")

    fp_vals = [m[0] for m in metrics]
    fn_vals = [m[1] for m in metrics]

    x = np.arange(len(labels))
    bar_width = 0.35

    fig, ax = plt.subplots(figsize=(fig_width, fig_width / GOLDEN_RATIO))
    bars_fp = ax.bar(x - bar_width / 2, fp_vals, bar_width,
                     label="Firstpass", color=PALETTE[1],  # sky blue
                     edgecolor="white", linewidth=0.35)
    bars_fn = ax.bar(x + bar_width / 2, fn_vals, bar_width,
                     label="Final", color=PALETTE[5],  # vermillion
                     edgecolor="white", linewidth=0.35)

    style_ax(ax, ylabel="Score",
             title=f"Firstpass vs Final — {mode_label}" if mode_label else
                   "Firstpass vs Final End Precision/Recall")
    ax.set_ylim(0, 1.05)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.legend()

    # Value labels on bars
    for bars in (bars_fp, bars_fn):
        for bar in bars:
            h = bar.get_height()
            if h > 0.02:
                ax.text(bar.get_x() + bar.get_width() / 2, h + 0.01,
                        f"{h:.2f}", ha="center", va="bottom", fontsize=5)

    pub_savefig(fig, output_path, close=True)
    log.info(f"  → {output_path}")


def write_comparison_tsv(firstpass_results, final_results, output_path,
                         mode_label=""):
    """Write side-by-side summary TSV."""
    metric_keys = [
        "n_isoforms_total", "n_jc_groups", "n_single_exon",
        "5prime_dedup_precision", "5prime_recall", "5prime_dedup_f1",
        "5prime_n_annot_matched", "5prime_n_annot_total",
        "3prime_dedup_precision", "3prime_recall", "3prime_dedup_f1",
        "3prime_n_annot_matched", "3prime_n_annot_total",
        "paired_dedup_precision", "paired_both_hit", "paired_unique_pairs",
    ]
    fields = ["mode", "stage"] + metric_keys
    with open(output_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for stage, res in [("firstpass", firstpass_results),
                           ("final", final_results)]:
            row = {"mode": mode_label, "stage": stage}
            for k in metric_keys:
                v = res.get(k)
                if isinstance(v, float):
                    row[k] = f"{v:.6f}"
                elif v is not None:
                    row[k] = v
                else:
                    row[k] = ""
            w.writerow(row)
    log.info(f"  → {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare firstpass vs final isoform end precision/recall."
    )
    parser.add_argument("--firstpass-bed", required=True,
                        help="Firstpass isoforms BED12")
    parser.add_argument("--final-bed", required=True,
                        help="Final isoforms BED12")
    parser.add_argument("--gtf", required=True,
                        help="Reference annotation GTF")
    parser.add_argument("--peaks-5prime", default=None,
                        help="BED6 CAGE peaks for TSS evaluation")
    parser.add_argument("--peaks-3prime", default=None,
                        help="BED6 dRNA/dRNA peaks for TTS evaluation")
    parser.add_argument("--window", type=int, default=50,
                        help="Max distance (bp) for end matching (default: 50)")
    parser.add_argument("--region", default=None,
                        help="Restrict annotation to region (e.g. chr22:16000000-26000000)")
    parser.add_argument("--mode", default="",
                        help="Label for the transcriptome mode")
    parser.add_argument("--outdir", required=True,
                        help="Output directory")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Parse region
    region_chrom, region_start, region_end = None, None, None
    if args.region:
        if ":" in args.region:
            region_chrom, coords = args.region.split(":", 1)
            region_start, region_end = [int(x) for x in coords.split("-")]
        else:
            region_chrom = args.region

    log.info(f"Parsing GTF ends (region={args.region or 'all'})...")
    annotated_ends = parse_gtf_ends(args.gtf, region_chrom, region_start, region_end)

    log.info("Parsing GTF transcripts for JC-matched recall...")
    annot_transcripts = parse_gtf_transcripts(args.gtf, region_chrom, region_start, region_end)
    log.info(f"  {len(annot_transcripts)} annotation transcripts")

    peaks_5prime = None
    peaks_3prime = None
    if args.peaks_5prime:
        log.info(f"Parsing 5' peaks: {args.peaks_5prime}")
        peaks_5prime = parse_peaks_bed(args.peaks_5prime, region_chrom, region_start, region_end)
    if args.peaks_3prime:
        log.info(f"Parsing 3' peaks: {args.peaks_3prime}")
        peaks_3prime = parse_peaks_bed(args.peaks_3prime, region_chrom, region_start, region_end)

    log.info("Computing firstpass P/R...")
    fp_results = run_precision_recall(
        args.firstpass_bed, annotated_ends, annot_transcripts,
        args.window, peaks_5prime, peaks_3prime)

    log.info("Computing final P/R...")
    fn_results = run_precision_recall(
        args.final_bed, annotated_ends, annot_transcripts,
        args.window, peaks_5prime, peaks_3prime)

    write_comparison_tsv(fp_results, fn_results,
                         outdir / "firstpass_vs_final.tsv",
                         mode_label=args.mode)

    plot_comparison(fp_results, fn_results,
                    outdir / "firstpass_vs_final.png",
                    mode_label=args.mode)

    # Log the deltas
    for end in ("5prime", "3prime"):
        for m in ("dedup_precision", "recall", "dedup_f1"):
            key = f"{end}_{m}"
            fp_v = fp_results.get(key, 0) or 0
            fn_v = fn_results.get(key, 0) or 0
            delta = fn_v - fp_v
            sign = "+" if delta >= 0 else ""
            log.info(f"  {end} {m}: {fp_v:.4f} → {fn_v:.4f} ({sign}{delta:.4f})")

    log.info("Done.")


if __name__ == "__main__":
    main()
