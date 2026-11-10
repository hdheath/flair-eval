#!/usr/bin/env python
"""
ted_pileup_tracks.py — genome-wide TSS/TTS read-end pileup tracks (IGV bigWig).

Reproduces TED's global peak-caller pileup (ted.py::_build_pileup_coverage +
ted_global_collect_approved_peaks): the per-base count of strand-aware read ends.
  + strand read: TSS = reference_start, TTS = reference_end-1
  - strand read: TSS = reference_end-1, TTS = reference_start
(TSS track = where transcripts START; TTS track = where they END.) Primary,
mapped, non-secondary/supplementary reads only — the population the caller bins.

Outputs per BAM:
  RAW    : <prefix>.{tss,tts}_{fwd,rev}.bw           (1-bp read-end counts)
  SMOOTH : <prefix>.{tss,tts}_{fwd,rev}.smooth.bw    (--smooth; the find_peaks
           input: gaussian σ * (σ*2.5), σ_tss=75 / σ_tts=15 by default, matching
           the A549 full-genome runs; RLE-merged + thresholded to bound size)
"""
import argparse
import os
import subprocess
import numpy as np
import pysam

BG2BW = "/private/groups/brookslab/bin/bedGraphToBigWig"
# Canonical human chromosomes: autosomes + sex. chrM and unplaced scaffolds
# (GL*/KI*) are intentionally EXCLUDED from all analyses (hdheath 2026-06-20).
MAIN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def write_raw_bedgraph(path, chrom, arr):
    nz = np.nonzero(arr)[0]
    if nz.size == 0:
        return
    vals = arr[nz]
    with open(path, "a") as fh:
        fh.write("".join(f"{chrom}\t{p}\t{p+1}\t{int(v)}\n"
                          for p, v in zip(nz.tolist(), vals.tolist())))


def write_smoothed_bedgraph(path, chrom, raw, sigma, thresh=0.1, dp=2):
    """Gaussian-smoothed pile-up (TED's find_peaks input), RLE-merged."""
    from scipy.ndimage import gaussian_filter1d
    sm = gaussian_filter1d(raw.astype(np.float64), sigma=float(sigma), mode="nearest") * (float(sigma) * 2.5)
    r = np.round(sm, dp)
    r[r < thresh] = 0.0
    if not (r > 0).any():
        return
    change = np.nonzero(np.diff(r) != 0)[0] + 1
    starts = np.concatenate([[0], change])
    ends = np.concatenate([change, [len(r)]])
    with open(path, "a") as fh:
        fh.write("".join(f"{chrom}\t{int(s)}\t{int(e)}\t{r[s]}\n"
                         for s, e in zip(starts.tolist(), ends.tolist()) if r[s] > 0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--chroms", nargs="+", default=MAIN_CHROMS)
    ap.add_argument("--smooth", action="store_true", help="also emit gaussian-smoothed (find_peaks-input) tracks")
    ap.add_argument("--sigma-tss", type=float, default=75.0)
    ap.add_argument("--sigma-tts", type=float, default=15.0)
    args = ap.parse_args()

    sizes = {}
    with open(args.fai) as fh:
        for line in fh:
            f = line.split("\t"); sizes[f[0]] = int(f[1])

    bam = pysam.AlignmentFile(args.bam, "rb")
    bam_chroms = set(bam.references)
    chroms = [c for c in args.chroms if c in sizes and c in bam_chroms]
    print(f"[{os.path.basename(args.prefix)}] chroms: {len(chroms)}  smooth={args.smooth}")

    raw_bg = {k: f"{args.prefix}.{k}.bedgraph" for k in ("tss_fwd", "tss_rev", "tts_fwd", "tts_rev")}
    sm_bg = {k: f"{args.prefix}.{k}.smooth.bedgraph" for k in raw_bg} if args.smooth else {}
    for p in list(raw_bg.values()) + list(sm_bg.values()):
        open(p, "w").close()
    sigma = {"tss_fwd": args.sigma_tss, "tss_rev": args.sigma_tss,
             "tts_fwd": args.sigma_tts, "tts_rev": args.sigma_tts}

    for chrom in chroms:
        n = sizes[chrom]
        a = {k: np.zeros(n, dtype=np.uint32) for k in raw_bg}
        nreads = 0
        for r in bam.fetch(chrom):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            if r.mapping_quality < args.min_mapq:
                continue
            s, e = r.reference_start, r.reference_end
            if s is None or e is None:
                continue
            e1 = min(e - 1, n - 1)
            if e1 < 0:
                continue
            if not r.is_reverse:
                a["tss_fwd"][s] += 1; a["tts_fwd"][e1] += 1
            else:
                a["tss_rev"][e1] += 1; a["tts_rev"][s] += 1
            nreads += 1
        for k in raw_bg:
            write_raw_bedgraph(raw_bg[k], chrom, a[k])
        if args.smooth:
            for k in raw_bg:
                write_smoothed_bedgraph(sm_bg[k], chrom, a[k], sigma[k])
        print(f"  {chrom}: {nreads} reads", flush=True)
        del a

    bam.close()
    cs = f"{args.prefix}.chrom.sizes"
    with open(cs, "w") as fh:
        for c in chroms:
            fh.write(f"{c}\t{sizes[c]}\n")

    for label, bgmap, suffix in [("raw", raw_bg, ""), ("smooth", sm_bg, ".smooth")]:
        for k, bg in bgmap.items():
            if not os.path.exists(bg) or os.path.getsize(bg) == 0:
                continue
            bw = f"{args.prefix}.{k}{suffix}.bw"
            sbg = bg + ".sorted"
            with open(sbg, "w") as out:
                subprocess.run(["sort", "-k1,1", "-k2,2n", bg], check=True, stdout=out,
                               env={**os.environ, "LC_ALL": "C"})
            subprocess.run([BG2BW, sbg, cs, bw], check=True)
            os.remove(bg); os.remove(sbg)
            print(f"  [bigWig] {bw}")


if __name__ == "__main__":
    main()
