#!/usr/bin/env python3
"""Cut a genome FASTA and GTF down to the chromosomes named in a regions BED.

Usage:
  flair_eval_subset_reference.py --genome G.fa[.gz] --gtf A.gtf[.gz] --regions R.bed --out DIR

Writes DIR/genome.fa (+ genome.fa.fai) and DIR/annotation.gtf holding whole chromosomes only, so
coordinates are unchanged. With --regions, FLAIR then never scans the rest of the genome (on
shared storage a 3 GB genome re-read per region made a 2 Mb test run take hours). Stdlib only.
"""
import argparse
import gzip
import sys
from pathlib import Path

WIDTH = 60


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def bed_chroms(path):
    chroms = []
    with opener(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c = line.split("\t", 1)[0].strip()
            if c and c not in chroms:
                chroms.append(c)
    return chroms


def write_fasta(genome, keep, out_fa, out_fai):
    found = []
    offset = 0
    with opener(genome) as src, open(out_fa, "w") as fa, open(out_fai, "w") as fai:
        name, chunks = None, []

        def flush():
            nonlocal offset
            if name is None or name not in keep:
                return
            seq = "".join(chunks)
            header = f">{name}\n"
            fa.write(header)
            offset += len(header)
            fai.write(f"{name}\t{len(seq)}\t{offset}\t{WIDTH}\t{WIDTH + 1}\n")
            for i in range(0, len(seq), WIDTH):
                line = seq[i:i + WIDTH] + "\n"
                fa.write(line)
                offset += len(line)
            found.append(name)

        for line in src:
            if line.startswith(">"):
                flush()
                name, chunks = line[1:].split()[0], []
            elif name in keep:
                chunks.append(line.strip())
        flush()
    return found


def write_gtf(gtf, keep, out_gtf):
    n = 0
    with opener(gtf) as src, open(out_gtf, "w") as out:
        for line in src:
            if line.startswith("#") or line.split("\t", 1)[0] in keep:
                out.write(line)
                n += not line.startswith("#")
    return n


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    chroms = bed_chroms(a.regions)
    if not chroms:
        sys.exit(f"flair_eval_subset_reference.py: {a.regions} names no chromosomes")
    a.out.mkdir(parents=True, exist_ok=True)
    found = write_fasta(a.genome, set(chroms), a.out / "genome.fa", a.out / "genome.fa.fai")
    missing = [c for c in chroms if c not in found]
    if missing:
        sys.exit(f"flair_eval_subset_reference.py: regions chromosome(s) {', '.join(missing)} are not in {a.genome}")
    n = write_gtf(a.gtf, set(chroms), a.out / "annotation.gtf")
    if n == 0:
        sys.exit(f"flair_eval_subset_reference.py: {a.gtf} has no records on {', '.join(chroms)}")
    print(f"kept {len(found)} chromosome(s): {', '.join(found)}; {n} GTF records")
    return 0


if __name__ == "__main__":
    sys.exit(main())
