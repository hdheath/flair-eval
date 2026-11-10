#!/usr/bin/env python3
"""
run_assembly_comparison_figs.py — driver for the assembly-comparison FIGURE SUITE.

Reads a multimethod manifest (pb_manifest.csv / ont_manifest.csv) and dispatches
the three figure scripts (Figs A/A2/B + C + D) for ONE chemistry, with the
CORRECTED within-50-of-interval metric + canonical chroms, both ends.

Per tool it picks the isoform input (isoforms_bed if present, else isoforms_gtf)
and the read-support source (isoform_counts if present, else isoform_read_map).

Figure scripts invoked (all in this dir):
  stratified_pr.py              -> <outdir>/stratified/   (Figs A/A2/B; all tools overlaid)
  sqanti_category_precise_bar.py-> <outdir>/sqanti/figC_<tool>.png (Fig C; per FLAIR/TED tool)
  end_dist_boxplot.py           -> <outdir>/end_dist/     (Fig D; all tools)

Usage (genome-wide; the driver hdheath runs AFTER chr22 validation):
  python run_assembly_comparison_figs.py \
      --manifest data/17_06_26_a549_multimethod/pb_manifest.csv \
      --outdir /private/home/hdheath/projects/plots_VALT/pb \
      [--region chrN]            # OMIT for genome-wide (canonical chroms auto-applied)
      [--fig-c-tools TED-hybrid-me1 FLAIR-default]
      [--window 50] [--min-support-figd 1]

Truth/GTF are read from the manifest columns (cage_peaks / drna_peaks / gtf) —
all rows in a manifest share them.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable
NO_TOKENS = ("NO_ISOFORMS_BED", "NO_ISOFORMS_GTF", "NO_COUNTS", "NO_TED_LOG",
             "NO_CDS", "NO_READ_MAP", "NO_")


def _present(path: str) -> bool:
    return bool(path) and not any(tok in path for tok in NO_TOKENS)


def _iso_input(row: dict) -> str:
    bed = row.get("isoforms_bed", "")
    gtf = row.get("isoforms_gtf", "")
    if _present(bed):
        return bed
    return gtf


def _counts_source(row: dict) -> str:
    cnt = row.get("isoform_counts", "")
    if _present(cnt):
        return cnt
    return row.get("isoform_read_map", "")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--region", default=None,
                    help="optional chrom/region; OMIT for genome-wide")
    ap.add_argument("--window", type=int, default=50)
    ap.add_argument("--min-support-figd", type=int, default=1)
    ap.add_argument("--fig-c-tools", nargs="+",
                    default=["TED-hybrid-me1", "FLAIR-default"],
                    help="tools to make a per-tool Fig C bar for (need BED12 input)")
    ap.add_argument("--all-contigs", action="store_true")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.manifest)))
    if not rows:
        sys.exit("empty manifest")
    cage = rows[0]["cage_peaks"]; drna = rows[0]["drna_peaks"]; gtf = rows[0]["gtf"]
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    labels, beds, counts = [], [], []
    for r in rows:
        lab = r["transcriptome_mode"]
        inp = _iso_input(r)
        cnt = _counts_source(r)
        if not _present(inp):
            print(f"  skip {lab}: no isoform input", file=sys.stderr)
            continue
        labels.append(lab); beds.append(inp)
        if _present(cnt):
            counts.append(f"{lab}:{cnt}")

    common = []
    if a.region:
        common += ["--region", a.region]
    if a.all_contigs:
        common += ["--all-contigs"]

    # ── Figs A/A2/B: stratified_pr (all tools overlaid) ──
    cmd = [PY, str(HERE / "stratified_pr.py"),
           "--beds", *beds, "--labels", *labels]
    if counts:
        cmd += ["--counts", *counts]
    cmd += ["--cage", cage, "--drna", drna, "--gtf", gtf,
            "--window", str(a.window), "--outdir", str(out / "stratified")] + common
    print("[run] stratified_pr", flush=True)
    subprocess.run(cmd, check=True)

    # ── Fig D: end_dist_boxplot (all tools) ──
    cmd = [PY, str(HERE / "end_dist_boxplot.py"),
           "--beds", *[f"{l}:{b}" for l, b in zip(labels, beds)]]
    if counts:
        cmd += ["--counts", *counts]
    cmd += ["--cage", cage, "--drna", drna, "--min-support", str(a.min_support_figd),
            "--outdir", str(out / "end_dist")]
    if a.region:
        cmd += ["--region", a.region]
    if a.all_contigs:
        cmd += ["--all-contigs"]
    print("[run] end_dist_boxplot", flush=True)
    subprocess.run(cmd, check=True)

    # ── Fig C: one bar plot per requested tool (needs BED12) ──
    cnt_map = {}
    for r in rows:
        lab = r["transcriptome_mode"]
        cnt_map[lab] = _counts_source(r)
    for tool in a.fig_c_tools:
        if tool not in labels:
            print(f"  Fig C skip {tool}: not in manifest", file=sys.stderr)
            continue
        idx = labels.index(tool)
        bed = beds[idx]
        if not (bed.endswith(".bed")):
            print(f"  Fig C skip {tool}: not a BED12 input ({bed})", file=sys.stderr)
            continue
        cmd = [PY, str(HERE / "sqanti_category_precise_bar.py"),
               "--bed", bed]
        if _present(cnt_map.get(tool, "")):
            cmd += ["--counts", cnt_map[tool]]
        cmd += ["--cage", cage, "--drna", drna, "--gtf", gtf,
                "--label", tool, "--window", str(a.window),
                "--min-support", str(1),
                "--outpng", str(out / "sqanti" / f"figC_{tool}.png")]
        if a.region:
            cmd += ["--region", a.region]
        if a.all_contigs:
            cmd += ["--all-contigs"]
        print(f"[run] sqanti_category_precise_bar {tool}", flush=True)
        subprocess.run(cmd, check=True)

    print("ASSEMBLY_COMPARISON_FIGS_DONE", out)


if __name__ == "__main__":
    main()
