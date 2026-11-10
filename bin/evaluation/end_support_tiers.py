#!/usr/bin/env python3
"""Do our 'imprecise' transcript ends actually have support?

The precision metric counts an end as wrong when it misses a called peak. But a
called peak is a thresholded object: an end can sit on real orthogonal signal that
simply fell below the peak caller's cutoff, or on a site annotated in an external
database, and still be scored as a false positive.

This classifies every emitted end into a support tier, per assembler:

  1 PEAK        matched a called truth peak (A549 CAGE 5' / A549 dRNA 3') -> scored TP
  2 SIGNAL      no peak, but raw orthogonal signal above --min-signal at the end
  3 DATABASE    no peak, no signal, but within the window of an external annotation
                (FANTOM5 CAGE peaks 5' / PolyA_DB 3')
  4 NONE        no support of any kind -> the genuine artifact candidate

Tiers are checked in order and are mutually exclusive, so they sum to the end count.
Ends are deduplicated per (chrom, strand, position): the question is about distinct
genomic positions, not how many isoforms share one.
"""
import argparse, bisect, collections, gzip, os, sys


def opener(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def load_intervals(path, stranded=True):
    d = collections.defaultdict(list)
    if not path or not os.path.exists(path):
        return d
    with opener(path) as fh:
        for line in fh:
            if line.startswith(('#', 'track', 'browser')):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 3:
                continue
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                continue
            st = f[5] if (stranded and len(f) > 5 and f[5] in '+-') else '.'
            d[(f[0], st)].append((s, e))
    for k in d:
        d[k].sort()
    return d


def load_signal(path, min_val):
    """bedGraph -> {(chrom,'.'): sorted[(start,end)]} keeping only intervals >= min_val."""
    d = collections.defaultdict(list)
    if not path or not os.path.exists(path):
        return d
    with opener(path) as fh:
        for line in fh:
            if line.startswith(('#', 'track', 'browser')):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 4:
                continue
            try:
                v = float(f[3])
            except ValueError:
                continue
            if abs(v) < min_val:
                continue
            d[(f[0], '.')].append((int(f[1]), int(f[2])))
    for k in d:
        d[k].sort()
    return d


def hit(db, chrom, strand, pos, window):
    for key in ((chrom, strand), (chrom, '.')):
        arr = db.get(key)
        if not arr:
            continue
        i = bisect.bisect_left(arr, (pos - window - 1_000_000,))
        for s, e in arr[i:]:
            if s > pos + window:
                break
            if s - window <= pos <= e + window:
                return True
    return False


def ends_from_bed(path):
    """-> set of (chrom, strand, pos5), set of (chrom, strand, pos3)"""
    f5, f3 = set(), set()
    with opener(path) as fh:
        for line in fh:
            if line.startswith(('#', 'track')):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 6:
                continue
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                continue
            c, st = f[0], f[5]
            f5.add((c, st, s if st != '-' else e))
            f3.add((c, st, e if st != '-' else s))
    return f5, f3


def ends_from_gtf(path):
    """Transcript bounds from exon records (robust to missing transcript lines)."""
    span = {}
    with opener(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] not in ('exon', 'transcript'):
                continue
            i = f[8].find('transcript_id "')
            if i < 0:
                continue
            tid = f[8][i + 15:f[8].find('"', i + 15)]
            s, e = int(f[3]) - 1, int(f[4])
            k = (f[0], f[6], tid)
            if k in span:
                a, b = span[k]
                span[k] = (min(a, s), max(b, e))
            else:
                span[k] = (s, e)
    f5, f3 = set(), set()
    for (c, st, _), (s, e) in span.items():
        f5.add((c, st, s if st != '-' else e))
        f3.add((c, st, e if st != '-' else s))
    return f5, f3


def classify(ends, peaks, signal, db, window):
    t = collections.Counter()
    for c, st, p in ends:
        if hit(peaks, c, st, p, window):
            t['1_PEAK'] += 1
        elif hit(signal, c, st, p, window):
            t['2_SIGNAL'] += 1
        elif hit(db, c, st, p, window):
            t['3_DATABASE'] += 1
        else:
            t['4_NONE'] += 1
    return t


def build_expanded_truth(peaks, db, signal, window, min_sig_flank=0):
    """Expanded truth = called peaks UNION database sites that carry local A549 signal.

    Adding every FANTOM/PolyA_DB entry would import pan-tissue sites not expressed in
    A549 and crater recall for everyone. Requiring local signal keeps only annotated
    sites with evidence of being active here — the population the peak caller plausibly
    missed by thresholding.
    """
    out = collections.defaultdict(list)
    for k, v in peaks.items():
        out[k].extend(v)
    added = 0
    for (c, st), arr in db.items():
        for s_, e_ in arr:
            mid = (s_ + e_) // 2
            if hit(signal, c, st, mid, window):
                out[(c, st)].append((s_, e_)); added += 1
    for k in out:
        out[k].sort()
    return out, added


def recall_of(truth, ends_by_cs, window):
    """fraction of truth intervals with a called end within +/-window"""
    found = tot = 0
    for (c, st), arr in truth.items():
        pos = ends_by_cs.get((c, st), [])
        if not pos:
            tot += len(arr); continue
        for s_, e_ in arr:
            tot += 1
            i = bisect.bisect_left(pos, s_ - window)
            if i < len(pos) and pos[i] <= e_ + window:
                found += 1
    return found, tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--arms', nargs='+', required=True, help='label=path pairs')
    ap.add_argument('--cage-peaks', required=True)
    ap.add_argument('--drna-peaks', required=True)
    ap.add_argument('--cage-signal', nargs='*', default=[])
    ap.add_argument('--drna-signal', nargs='*', default=[])
    ap.add_argument('--fantom', required=True)
    ap.add_argument('--polyadb', required=True)
    ap.add_argument('--window', type=int, default=50)
    ap.add_argument('--min-signal', type=float, default=0.1)
    ap.add_argument('--chroms', nargs='*', default=None)
    a = ap.parse_args()

    print(f"loading evidence (window ±{a.window}bp, signal floor {a.min_signal})", file=sys.stderr)
    cage_pk, drna_pk = load_intervals(a.cage_peaks), load_intervals(a.drna_peaks)
    cage_sg = collections.defaultdict(list)
    for p in a.cage_signal:
        for k, v in load_signal(p, a.min_signal).items():
            cage_sg[k].extend(v)
    drna_sg = collections.defaultdict(list)
    for p in a.drna_signal:
        for k, v in load_signal(p, a.min_signal).items():
            drna_sg[k].extend(v)
    for d in (cage_sg, drna_sg):
        for k in d:
            d[k].sort()
    fantom, polyadb = load_intervals(a.fantom), load_intervals(a.polyadb)

    keep = set(a.chroms) if a.chroms else None
    if keep:
        # The truth denominator MUST be restricted to the same chromosomes as the
        # called ends, or recall is divided by a genome-wide count and deflated.
        def _restrict(d):
            return collections.defaultdict(list,
                {k: v for k, v in d.items() if k[0] in keep})
        cage_pk, drna_pk = _restrict(cage_pk), _restrict(drna_pk)
        cage_sg, drna_sg = _restrict(cage_sg), _restrict(drna_sg)
        fantom, polyadb = _restrict(fantom), _restrict(polyadb)
    TIERS = ['1_PEAK', '2_SIGNAL', '3_DATABASE', '4_NONE']
    for end_kind, pk, sg, db in (("5'", cage_pk, cage_sg, fantom),
                                 ("3'", drna_pk, drna_sg, polyadb)):
        print(f"\n=== {end_kind} ends ===")
        exp_truth, n_added = build_expanded_truth(pk, db, sg, a.window)
        n_pk = sum(len(v) for v in pk.values()); n_ex = sum(len(v) for v in exp_truth.values())
        print(f"  truth: {n_pk} called peaks -> {n_ex} expanded "
              f"(+{n_added} database sites carrying A549 signal)")
        print(f"{'':<22}{'':<7}{'--- PEAK truth ---':^24}{'--- EXPANDED truth ---':^26}")
        print(f"{'arm':<22}{'ends':>7}{'P':>8}{'R':>8}{'F1':>8}   {'P':>8}{'R':>8}{'F1':>8}{'unsup':>9}")
        print('-' * 96)
        for spec in a.arms:
            lbl, path = spec.split('=', 1)
            if not os.path.exists(path):
                print(f"{lbl:<22}  (missing)"); continue
            f5, f3 = (ends_from_gtf(path) if path.endswith(('.gtf', '.gff', '.gtf.gz'))
                      else ends_from_bed(path))
            ends = f5 if end_kind == "5'" else f3
            if keep:
                ends = {e for e in ends if e[0] in keep}
            t = classify(ends, pk, sg, db, a.window)
            n = sum(t.values()) or 1
            sup = t['1_PEAK'] + t['2_SIGNAL'] + t['3_DATABASE']
            by_cs = collections.defaultdict(list)
            for c, st, p_ in ends:
                by_cs[(c, st)].append(p_)
            for k in by_cs:
                by_cs[k].sort()
            fp_, tp_ = recall_of(pk, by_cs, a.window)
            fe_, te_ = recall_of(exp_truth, by_cs, a.window)
            p_pk, r_pk = t['1_PEAK']/n, (fp_/tp_ if tp_ else 0)
            p_ex, r_ex = sup/n, (fe_/te_ if te_ else 0)
            f_pk = 2*p_pk*r_pk/(p_pk+r_pk) if (p_pk+r_pk) else 0
            f_ex = 2*p_ex*r_ex/(p_ex+r_ex) if (p_ex+r_ex) else 0
            print(f"{lbl:<22}{n:>7}{p_pk:>8.3f}{r_pk:>8.3f}{f_pk:>8.3f}   "
                  f"{p_ex:>8.3f}{r_ex:>8.3f}{f_ex:>8.3f}{t['4_NONE']:>9}")


if __name__ == '__main__':
    main()
