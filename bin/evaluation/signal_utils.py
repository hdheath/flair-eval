"""Shared utilities for BED12/GTF parsing, bedGraph signal tracks, and KDE density.

Provides:
  - parse_bed12(path)               Full BED12 parser → list of dicts
  - parse_bed12_by_name(path)       Same but returns {name: dict}
  - parse_gtf(path)                 GTF parser → list of dicts (same shape)
  - parse_isoforms(path)            Auto-detect BED12 or GTF and parse
  - tss_tts(start, end, strand)     Strand-aware TSS/TTS extraction
  - group_by_junction_chain(isos)   Group multi-exon isoforms by SJC
  - gene_from_name(name)            Extract ENSG gene ID from an isoform name
  - BedGraphTrack                   Numpy-backed bedGraph reader with query()
  - isoform_signal(...)             Strand-aware TSS/TTS signal
  - kde_density(x, y)               Log-space Gaussian KDE, normalised 0–1

All evaluation scripts should import shared parsing from here:
    from signal_utils import parse_bed12, tss_tts, group_by_junction_chain
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from scipy.stats import gaussian_kde


# ── Constants ───────────────────────────────────────────────────────────────

SIG_WINDOW = 50          # bp window around each end for signal lookup
END_BIN    = 50          # bp bin size for binning TSS/TTS coordinates
_ENSG_RE   = re.compile(r"(ENSG\d+(?:\.\d+)?)")


# ── TSS / TTS helpers ──────────────────────────────────────────────────────

def tss_tts(start: int, end: int, strand: str) -> Tuple[int, int]:
    """Return (tss, tts) given BED start, end and strand."""
    return (start, end) if strand == '+' else (end, start)


# ── BED12 parsing ──────────────────────────────────────────────────────────

def parse_bed12(path: str | Path) -> List[dict]:
    """Parse a BED12 file and return a list of isoform dicts.

    Each dict contains:
      chrom, start, end, name, score, strand, tss, tts,
      junctions (tuple of (donor, acceptor) 2-tuples), n_exons (int),
      spliced_len (int, sum of exon/block sizes).

    Splice junctions are derived from block starts/sizes.
    """
    isoforms: List[dict] = []
    p = Path(path)
    if not p.exists():
        return isoforms
    with open(p) as f:
        for line in f:
            if line.startswith(("#", "track")):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 12:
                continue
            chrom, start, end = c[0], int(c[1]), int(c[2])
            name, strand = c[3], c[5]
            try:
                score = int(c[4])
            except ValueError:
                score = 0
            bc = int(c[9])
            bsz = [int(x) for x in c[10].rstrip(",").split(",") if x]
            bst = [int(x) for x in c[11].rstrip(",").split(",") if x]
            juncs: List[Tuple[int, int]] = []
            for i in range(bc - 1):
                juncs.append((start + bst[i] + bsz[i], start + bst[i + 1]))
            tss, tts = tss_tts(start, end, strand)
            isoforms.append(dict(
                chrom=chrom, start=start, end=end, name=name,
                score=score, strand=strand, tss=tss, tts=tts,
                junctions=tuple(juncs), n_exons=bc,
                spliced_len=sum(bsz),
            ))
    return isoforms


def parse_bed12_by_name(path: str | Path) -> Dict[str, dict]:
    """Parse a BED12 file and return {isoform_name: dict}.

    Same fields as :func:`parse_bed12` but keyed by isoform name.
    """
    return {iso["name"]: iso for iso in parse_bed12(path)}


# ── Junction-chain grouping ────────────────────────────────────────────────

def group_by_junction_chain(
    isoforms,
    min_exons: int = 2,
) -> Dict[Tuple[str, str, tuple], List[dict]]:
    """Group multi-exon isoforms by (chrom, strand, junction_chain).

    Accepts either a list of dicts (from :func:`parse_bed12`) or a
    dict-of-dicts keyed by name (from :func:`parse_bed12_by_name` or
    legacy parsers).

    Returns {(chrom, strand, junctions_tuple): [isoform_dicts]}.
    Single-exon isoforms (n_exons < *min_exons*) are excluded.
    """
    groups: Dict[Tuple[str, str, tuple], List[dict]] = defaultdict(list)
    items = isoforms.values() if isinstance(isoforms, dict) else isoforms
    for iso in items:
        # Support both list-of-dicts and dict-of-dicts with name in value
        if isinstance(iso, str):
            # dict-of-dicts: iso is the key, isoforms[iso] is the value
            info = isoforms[iso]
            if "name" not in info:
                info = {"name": iso, **info}
        else:
            info = iso
        junctions = tuple(info.get("junctions", ()))
        n_exons = int(info.get("n_exons", 1))
        if n_exons < min_exons or not junctions:
            continue
        key = (info["chrom"], info["strand"], junctions)
        groups[key].append(info)
    return dict(groups)


def gene_from_name(name: str) -> str:
    """Extract ENSG gene ID from an isoform name (e.g. ENST…_ENSG…).

    Falls back to the substring after the last underscore if no ENSG
    pattern is found.
    """
    m = _ENSG_RE.search(name)
    return m.group(1) if m else name.rsplit("_", 1)[-1]


# ── GTF parsing ────────────────────────────────────────────────────────────

_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def parse_gtf(path: str | Path) -> List[dict]:
    """Parse a GTF file and return isoform dicts (same shape as parse_bed12).

    Groups exons by transcript_id, builds junctions from exon boundaries,
    and returns one dict per transcript with: chrom, start, end, name,
    score, strand, junctions, n_exons.
    """
    p = Path(path)
    if not p.exists():
        return []

    # Collect exons per transcript
    tx_exons: Dict[str, List[Tuple[str, int, int, str]]] = defaultdict(list)
    tx_gene: Dict[str, str] = {}

    with open(p) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            if cols[2] != "exon":
                continue
            chrom = cols[0]
            start = int(cols[3]) - 1  # GTF is 1-based → 0-based
            end = int(cols[4])
            strand = cols[6]
            attrs = dict(_ATTR_RE.findall(cols[8]))
            tid = attrs.get("transcript_id", "")
            gid = attrs.get("gene_id", "")
            if not tid:
                continue
            tx_exons[tid].append((chrom, start, end, strand))
            if gid:
                tx_gene[tid] = gid

    isoforms: List[dict] = []
    for tid, exons in tx_exons.items():
        if not exons:
            continue
        exons.sort(key=lambda x: x[1])  # sort by start
        chrom = exons[0][0]
        strand = exons[0][3]
        tx_start = exons[0][1]
        tx_end = exons[-1][2]
        # Build splice junctions from consecutive exon boundaries
        juncs: List[Tuple[int, int]] = []
        for i in range(len(exons) - 1):
            juncs.append((exons[i][2], exons[i + 1][1]))
        gid = tx_gene.get(tid, "")
        name = f"{tid}_{gid}" if gid else tid
        tx_tss, tx_tts = tss_tts(tx_start, tx_end, strand)
        isoforms.append(dict(
            chrom=chrom, start=tx_start, end=tx_end, name=name,
            score=0, strand=strand, tss=tx_tss, tts=tx_tts,
            junctions=tuple(juncs), n_exons=len(exons),
        ))
    return isoforms


def parse_isoforms(path: str | Path) -> List[dict]:
    """Auto-detect BED12 or GTF and parse accordingly."""
    p = Path(path)
    if p.suffix in (".gtf", ".gff", ".gff3"):
        return parse_gtf(p)
    return parse_bed12(p)


# ── BedGraph signal track ──────────────────────────────────────────────────

class BedGraphTrack:
    """Numpy-backed bedGraph interval store with fast overlap queries.

    Usage::

        track = BedGraphTrack.from_file("signal.bedgraph")
        mean_signal = track.query("chr22", 10000, 10100)
    """

    def __init__(self) -> None:
        self.data: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    @classmethod
    def from_file(
        cls,
        path: str | Path,
        chroms: Optional[Set[str]] = None,
    ) -> "BedGraphTrack":
        """Load a bedGraph file into a BedGraphTrack.

        Parameters
        ----------
        path : path to the bedGraph file
        chroms : optional set of chromosome names to restrict loading
        """
        track = cls()
        buf: Dict[str, list] = defaultdict(list)
        with open(path) as f:
            for line in f:
                if line[0] in ("#", "t", "b"):
                    w = line.split(None, 1)[0] if line.strip() else ""
                    if w in ("track", "browser") or line[0] == "#":
                        continue
                p = line.split()
                if len(p) < 4:
                    continue
                ch = p[0]
                if chroms and ch not in chroms:
                    continue
                buf[ch].append((int(p[1]), int(p[2]), float(p[3])))
        for ch, rows in buf.items():
            rows.sort()
            s, e, v = zip(*rows)
            track.data[ch] = (
                np.array(s, dtype=np.int64),
                np.array(e, dtype=np.int64),
                np.array(v, dtype=np.float64),
            )
        return track

    def query(self, chrom: str, start: int, end: int) -> float:
        """Return the mean signal value over *[start, end)*.

        Computes a coverage-weighted average.  Returns 0.0 when no
        intervals overlap the query region.
        """
        if end <= start or chrom not in self.data:
            return 0.0
        starts, ends, vals = self.data[chrom]
        n = len(starts)
        if n == 0:
            return 0.0
        idx = int(np.searchsorted(starts, start, side="right")) - 1
        if idx < 0:
            idx = 0
        while idx < n and ends[idx] <= start:
            idx += 1
        total = 0.0
        while idx < n and starts[idx] < end:
            ov_s = max(start, int(starts[idx]))
            ov_e = min(end, int(ends[idx]))
            if ov_s < ov_e:
                total += vals[idx] * (ov_e - ov_s)
            idx += 1
        span = end - start
        return total / span if span > 0 else 0.0


# ── Signal helpers ──────────────────────────────────────────────────────────

def isoform_signal(
    iso: dict,
    cage_p: BedGraphTrack,
    cage_m: BedGraphTrack,
    qs_p: BedGraphTrack,
    qs_m: BedGraphTrack,
    window: int = SIG_WINDOW,
) -> Tuple[float, float]:
    """Return (TSS_signal, TTS_signal) for one isoform.

    TSS signal comes from CAGE; TTS signal from QuantSeq.
    The correct strand track is chosen automatically.
    """
    ch = iso["chrom"]
    if iso["strand"] == "+":
        tss = cage_p.query(ch, iso["start"] - window, iso["start"] + window)
        tts = qs_p.query(ch,   iso["end"]   - window, iso["end"]   + window)
    else:
        tss = cage_m.query(ch, iso["end"]   - window, iso["end"]   + window)
        tts = qs_m.query(ch,   iso["start"] - window, iso["start"] + window)
    return tss, tts


def load_signal_tracks(
    cage_plus: str,
    cage_minus: str,
    qs_plus: str,
    qs_minus: str,
    chroms: Optional[Set[str]] = None,
) -> Tuple[BedGraphTrack, BedGraphTrack, BedGraphTrack, BedGraphTrack]:
    """Load four bedGraph signal tracks and return (cage_p, cage_m, qs_p, qs_m)."""
    return (
        BedGraphTrack.from_file(cage_plus,  chroms),
        BedGraphTrack.from_file(cage_minus, chroms),
        BedGraphTrack.from_file(qs_plus,    chroms),
        BedGraphTrack.from_file(qs_minus,   chroms),
    )


# ── KDE density estimation ─────────────────────────────────────────────────

def kde_density(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Gaussian KDE density at each point (log-space).

    Returns 0–1 normalised density values suitable for scatter colouring.
    """
    lx = np.log10(np.maximum(x, 1e-10))
    ly = np.log10(np.maximum(y, 1e-10))
    xy = np.vstack([lx, ly])
    try:
        kde = gaussian_kde(xy, bw_method="scott")
        d = kde(xy)
    except np.linalg.LinAlgError:
        d = np.ones(len(x))
    dmax = d.max()
    return d / dmax if dmax > 0 else d
