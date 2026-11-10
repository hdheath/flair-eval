#!/usr/bin/env python3
"""OR-rescue TED end verdicts using external TSS / PAS databases.

An end is accepted if it EITHER passed TED's own gate (column 31, written by
`--ted_tag_ends`) OR lies within `--window` bp of a same-strand annotated end in
a reference database — FANTOM5 CAGE peaks for the 5' side, PolyA_DB for the 3'.

Motivation: the hybrid 3' gate cannot clear its threshold on depth + annotation
alone (0.67 + 0.15 = 0.82 < 0.84), so a well-supported end sitting on an
annotated PAS can still be rejected when the sequence model dislikes it. This
asks what recall that costs, and what precision the rescue gives back.

Writes the same filter modes as filter_isoforms_by_ted_tag.py so the two are
directly comparable: 5prime / 3prime / paired / any / none.
"""
import argparse, bisect, collections, os, sys


def load_db(path, pos_col, want_summit=False):
    """-> {(chrom, strand): sorted[int positions]}"""
    d = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith(('#', 'track', 'browser')):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 6:
                continue
            # FANTOM BED9: cols 7/8 are the peak summit; prefer it over the span.
            if want_summit and len(f) >= 8 and f[6].isdigit():
                pos = int(f[6])
            else:
                pos = int(f[pos_col])
            d[(f[0], f[5])].append(pos)
    for k in d:
        d[k].sort()
    return d


def near(db, chrom, strand, pos, window):
    arr = db.get((chrom, strand))
    if not arr:
        return False
    i = bisect.bisect_left(arr, pos - window)
    return i < len(arr) and arr[i] <= pos + window


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in-bed', required=True, help='tagged .isoforms.bed (31 cols)')
    ap.add_argument('--out-bed', required=True)
    ap.add_argument('--mode', required=True,
                    choices=['5prime', '3prime', 'paired', 'any', 'none'])
    ap.add_argument('--tss-db', help='FANTOM5 CAGE peaks BED (5-prime rescue)')
    ap.add_argument('--pas-db', help='PolyA_DB BED (3-prime rescue)')
    ap.add_argument('--window', type=int, default=50)
    ap.add_argument('--stats', action='store_true', help='report rescue counts to stderr')
    a = ap.parse_args()

    tss_db = load_db(a.tss_db, 1, want_summit=True) if a.tss_db else {}
    pas_db = load_db(a.pas_db, 1) if a.pas_db else {}

    n = kept = r5 = r3 = 0
    with open(a.in_bed) as fh, open(a.out_bed, 'w') as out:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 12:
                continue
            n += 1
            # A 30-column (untagged) bed has no verdict to rescue: pass through,
            # matching filter_isoforms_by_ted_tag.py's behaviour exactly.
            if len(f) < 31:
                out.write(line); kept += 1; continue
            chrom, strand = f[0], f[5]
            start, end = int(f[1]), int(f[2])
            tag = f[30] or ''
            p5 = start if strand != '-' else end
            p3 = end if strand != '-' else start

            ok5 = '5_prime_TED' in tag
            ok3 = '3_prime_TED' in tag
            if not ok5 and tss_db and near(tss_db, chrom, strand, p5, a.window):
                ok5 = True; r5 += 1
            if not ok3 and pas_db and near(pas_db, chrom, strand, p3, a.window):
                ok3 = True; r3 += 1

            keep = {'none': True, '5prime': ok5, '3prime': ok3,
                    'paired': ok5 and ok3, 'any': ok5 or ok3}[a.mode]
            if keep:
                out.write(line); kept += 1

    print(f"rescue_ends_by_db: mode={a.mode} window=±{a.window} "
          f"kept {kept}/{n}; rescued 5'={r5} 3'={r3}", file=sys.stderr)


if __name__ == '__main__':
    main()
