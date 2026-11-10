"""Random access to a FASTA file, for splice-site motifs.

Uses `<fasta>.fai` when present (samtools faidx format). Without one the file is scanned
once to build the same index in memory; nothing is written next to the input.
Plain (uncompressed) FASTA only.
"""

from __future__ import annotations

import os
from typing import Dict, Tuple

from .annotation import InputError

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# Intron motifs read on the + strand of the genome: (first 2 nt, last 2 nt).
CANONICAL_PLUS = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
CANONICAL_MINUS = {("CT", "AC"), ("CT", "GC"), ("GT", "AT")}


class Fasta:
    def __init__(self, path: str):
        if path.endswith(".gz"):
            raise InputError(f"{path}: gzipped FASTA is not supported; pass the plain FASTA")
        self.path = path
        fai = path + ".fai"
        self.index: Dict[str, Tuple[int, int, int, int]] = (
            self._read_fai(fai) if os.path.exists(fai) else self._scan())
        self._fh = open(path, "rb")

    @staticmethod
    def _read_fai(fai: str) -> Dict[str, Tuple[int, int, int, int]]:
        idx = {}
        with open(fai) as fh:
            for line in fh:
                p = line.split("\t")
                idx[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
        return idx

    def _scan(self) -> Dict[str, Tuple[int, int, int, int]]:
        idx: Dict[str, Tuple[int, int, int, int]] = {}
        name = None
        length = offset = lbases = lbytes = 0
        pos = 0
        with open(self.path, "rb") as fh:
            for raw in fh:
                if raw.startswith(b">"):
                    if name is not None:
                        idx[name] = (length, offset, lbases, lbytes)
                    name = raw[1:].split()[0].decode()
                    length = lbases = lbytes = 0
                    offset = pos + len(raw)
                else:
                    seq = raw.rstrip(b"\r\n")
                    if lbases == 0 and seq:
                        lbases, lbytes = len(seq), len(raw)
                    length += len(seq)
                pos += len(raw)
        if name is not None:
            idx[name] = (length, offset, lbases, lbytes)
        if not idx:
            raise InputError(f"{self.path}: no FASTA records")
        return idx

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Sequence of [start, end), 0-based half-open, upper case. '' if off the contig."""
        rec = self.index.get(chrom)
        if rec is None or start < 0 or end > rec[0] or end <= start:
            return ""
        length, offset, lbases, lbytes = rec
        first = offset + (start // lbases) * lbytes + start % lbases
        last = offset + ((end - 1) // lbases) * lbytes + (end - 1) % lbases
        self._fh.seek(first)
        raw = self._fh.read(last - first + 1)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode().upper()

    def is_canonical(self, chrom: str, strand: str, a: int, b: int) -> bool:
        """True for a GT-AG, GC-AG or AT-AC intron [a, b) on the given strand."""
        motif = (self.fetch(chrom, a, a + 2), self.fetch(chrom, b - 2, b))
        if strand == "+":
            return motif in CANONICAL_PLUS
        if strand == "-":
            return motif in CANONICAL_MINUS
        return motif in CANONICAL_PLUS or motif in CANONICAL_MINUS

    def close(self):
        self._fh.close()
