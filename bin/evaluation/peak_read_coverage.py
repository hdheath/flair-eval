#!/usr/bin/env python3
"""How much long-read support does each orthogonal peak actually have?

Recall counts a peak as missed when no assembler placed an end there — but that
conflates two very different failures: (a) no long read terminates near the peak at
all, so no method could ever recover it, and (b) reads are present and the end call
was lost. This separates them, and asks whether the two assays agree on expression.

For every peak (A549 CAGE for 5', A549 dRNA 3'-seq for 3'):
  * n_read_ends  — long reads whose strand-aware end lies within --window bp
  * peak_score   — the peak's own expression value (BED column 5)
Emits a per-peak TSV plus the coverage distribution and the rank correlation
between the two expression measures.
"""
import argparse, bisect, collections, sys
import pysam


def peak_positions(path, chroms):
    out = []
    for line in open(path):
        if line.startswith(('#', 'track', 'browser')):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 6 or f[0] not in chroms:
            continue
        try:
            s, e, sc = int(f[1]), int(f[2]), float(f[4])
        except ValueError:
            continue
        out.append((f[0], s, e, f[5], sc, f[3]))
    return out


def read_ends(bam, chroms, mapq=0):
    """-> {(chrom,strand): sorted[pos]} for 5' and 3' separately"""
    five = collections.defaultdict(list)
    three = collections.defaultdict(list)
    af = pysam.AlignmentFile(bam, 'rb')
    n = 0
    for c in chroms:
        for r in af.fetch(c):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            if r.mapping_quality < mapq:
                continue
            n += 1
            st = '-' if r.is_reverse else '+'
            s, e = r.reference_start, r.reference_end
            five[(c, st)].append(s if st == '+' else e)
            three[(c, st)].append(e if st == '+' else s)
    for d in (five, three):
        for k in d:
            d[k].sort()
    return five, three, n


def count_near(arr, lo, hi):
    """Count integer end coordinates in the half-open interval [lo, hi)."""
    return bisect.bisect_left(arr, hi) - bisect.bisect_left(arr, lo)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bam', required=True)
    ap.add_argument('--cage-peaks', required=True)
    ap.add_argument('--drna-peaks', required=True)
    ap.add_argument('--chroms', nargs='+', required=True)
    ap.add_argument('--window', type=int, default=50)
    ap.add_argument('--out-prefix', required=True)
    a = ap.parse_args()

    chroms = set(a.chroms)
    print(f"reading {a.bam} over {sorted(chroms)}", file=sys.stderr)
    five, three, n_reads = read_ends(a.bam, a.chroms)
    print(f"  {n_reads:,} primary alignments", file=sys.stderr)

    BINS = [(0, 0), (1, 2), (3, 5), (6, 10), (11, 25), (26, 100), (101, 10**9)]
    LBL = ['0', '1-2', '3-5', '6-10', '11-25', '26-100', '>100']

    for kind, peaks_path, ends in (("5'", a.cage_peaks, five), ("3'", a.drna_peaks, three)):
        peaks = peak_positions(peaks_path, chroms)
        rows = []
        for c, s, e, st, sc, name in peaks:
            arr = ends.get((c, st), [])
            n = count_near(arr, s - a.window, e + a.window) if arr else 0
            rows.append((c, s, e, st, sc, n, name))
        out = f"{a.out_prefix}_{'5p' if kind == chr(53)+chr(39) else '3p'}_peak_coverage.tsv"
        out = f"{a.out_prefix}_{'5p' if kind.startswith(chr(53)) else '3p'}_peak_coverage.tsv"
        with open(out, 'w') as fh:
            fh.write("chrom\tstart\tend\tstrand\tpeak_score\tn_read_ends\tpeak_name\n")
            for r in rows:
                fh.write('\t'.join(map(str, r)) + '\n')
        dist = collections.Counter()
        for *_, n, _n2 in [(r[0], r[1], r[2], r[3], r[4], r[5], r[6]) for r in rows]:
            for (lo, hi), lbl in zip(BINS, LBL):
                if lo <= n <= hi:
                    dist[lbl] += 1
                    break
        tot = len(rows)
        print(f"\n=== {kind} peaks: long-read end coverage within ±{a.window}bp ===")
        print(f"  {tot:,} peaks; wrote {out}")
        cum = 0
        for lbl in LBL:
            cum += dist[lbl]
            print(f"    {lbl:>7} reads : {dist[lbl]:>6} peaks ({100*dist[lbl]/tot:>5.1f}%)   cumulative {100*cum/tot:>5.1f}%")
        zero = dist['0']
        print(f"  --> {zero} peaks ({100*zero/tot:.1f}%) have ZERO long-read ends: unreachable, "
              f"caps recall at {100*(tot-zero)/tot:.1f}%")
        # rank correlation between the two expression measures, on covered peaks
        pairs = [(r[4], r[5]) for r in rows if r[5] > 0]
        if len(pairs) > 10:
            def rank(v):
                order = sorted(range(len(v)), key=lambda i: v[i])
                rk = [0.0]*len(v)
                for i, o in enumerate(order):
                    rk[o] = i
                return rk
            xs = rank([p[0] for p in pairs]); ys = rank([p[1] for p in pairs])
            mx, my = sum(xs)/len(xs), sum(ys)/len(ys)
            num = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
            den = (sum((x-mx)**2 for x in xs) * sum((y-my)**2 for y in ys)) ** 0.5
            print(f"  Spearman(peak_score, n_read_ends) on {len(pairs)} covered peaks = {num/den if den else float('nan'):.3f}")


if __name__ == '__main__':
    main()
