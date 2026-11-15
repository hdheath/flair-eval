#!/usr/bin/env python3
"""Build one line's 5' truth union as a truth BED (ted-clustering map ticket 05, decision 3).

    build_truth_union.py --out A549.union5.truth.bed CAGE=A549.CAGE.truth.bed PROcap=A549.PROcap.truth.bed

The sets are given in tier order (CAGE, then PRO-cap, then RAMPAGE). Every call of the first
set is kept at its reference point. Each next set's calls pair one-to-one with the calls kept so
far, same strand, within --window bp (the greedy sweep of the ortho-truth agreement baselines:
calls in position order each take the leftmost free kept call in range, a maximum matching for
equal windows on a line). A paired call is the same site and is dropped, so a site seen by more
than one set counts once; an unpaired call is added. A union is then one more entry on the
--truth sheet, scored like any truth set.

Input calls are 1-bp reference points (the ortho-truth BEDs); '#' lines are headers. The output
is BED6 plus the call's source set, sorted by chrom, point and strand, with a '#' header.
"""

import argparse
import sys
from bisect import bisect_left
from collections import defaultdict
from typing import List, NamedTuple, Sequence, Tuple

WINDOW = 50


class Call(NamedTuple):
    chrom: str
    point: int
    strand: str
    call_id: str
    score: str
    set: str = ""


def read_calls(path: str, set_name: str = "") -> List[Call]:
    """The 1-bp calls of a truth BED (columns 1-6), labelled with set_name."""
    calls = []
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                raise ValueError(f"{path}:{n}: a truth BED needs 6 columns")
            start, end = int(f[1]), int(f[2])
            if end - start != 1:
                raise ValueError(f"{path}:{n}: expected a 1-bp reference point, got {f[0]}:{start}-{end}")
            calls.append(Call(f[0], start, f[5], f[3], f[4], set_name))
    return calls


def _unpaired(kept: Sequence[int], new: Sequence[Tuple[int, int]], window: int) -> List[int]:
    """Indices of `new` ((point, index), sorted) left unpaired by the sweep against sorted `kept`."""
    out, j = [], 0
    for p, i in new:
        j = max(j, bisect_left(kept, p - window))   # kept calls before j are taken or out of range
        if j < len(kept) and kept[j] <= p + window:
            j += 1
        else:
            out.append(i)
    return out


def build_union(tiers: Sequence[Tuple[str, Sequence[Call]]], window: int = WINDOW) -> List[Call]:
    """The union of (set name, calls) tiers, in tier order (see the module docstring)."""
    kept: dict = defaultdict(list)          # (chrom, strand) -> [Call]
    for name, calls in tiers:
        calls = [c._replace(set=name) for c in calls]
        if not kept:
            for c in calls:
                kept[(c.chrom, c.strand)].append(c)
            continue
        by_key = defaultdict(list)
        for i, c in enumerate(calls):
            by_key[(c.chrom, c.strand)].append((c.point, i))
        added = []
        for key, new in by_key.items():
            points = sorted(c.point for c in kept.get(key, []))
            added += [calls[i] for i in _unpaired(points, sorted(new), window)]
        for c in added:
            kept[(c.chrom, c.strand)].append(c)
    return sorted((c for cs in kept.values() for c in cs), key=lambda c: (c.chrom, c.point, c.strand))


def write_bed(path: str, calls: Sequence[Call]) -> None:
    with open(path, "w") as fh:
        fh.write("#chrom\tstart\tend\tcall_id\tscore\tstrand\tsource_set\n")
        for c in calls:
            fh.write(f"{c.chrom}\t{c.point}\t{c.point + 1}\t{c.call_id}\t{c.score}\t{c.strand}\t{c.set}\n")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True, help="the union truth BED to write")
    p.add_argument("--window", type=int, default=WINDOW, help="pairing window, bp (default %(default)s)")
    p.add_argument("sets", nargs="+", metavar="NAME=BED", help="truth sets in tier order")
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = parse_args(argv)
    if len(a.sets) < 2:
        print("build_truth_union.py: error: a union needs at least two truth sets", file=sys.stderr)
        return 2
    tiers = []
    for s in a.sets:
        name, sep, path = s.partition("=")
        if not sep or not name or not path:
            print(f"build_truth_union.py: error: {s!r}: expected NAME=BED", file=sys.stderr)
            return 2
        try:
            tiers.append((name, read_calls(path, name)))
        except (OSError, ValueError) as e:
            print(f"build_truth_union.py: error: {e}", file=sys.stderr)
            return 2
    union = build_union(tiers, a.window)
    write_bed(a.out, union)
    counts = {name: sum(c.set == name for c in union) for name, _ in tiers}
    print(f"build_truth_union.py: {len(union)} calls ("
          + ", ".join(f"{n} {k}" for n, k in counts.items()) + f") -> {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
