#!/usr/bin/env python3
"""AltEndCDSAnalysis: alt-TSS / alt-TTS sibling pairs and their CDS / UTR
biology divergence.

For each junction-chain (SJC) group with >=2 isoforms, this analysis
asks: when TED keeps multiple alt-end siblings on the same chain, do
the alternative TSSs / TTSs produce different protein biology?

Specifically per pair of sibling isoforms in an SJC group:
  delta_5utr_len      change in 5' UTR length (bp)
  delta_3utr_len      change in 3' UTR length (bp)
  delta_productivity  productivity-class transition (e.g. PRO->PTC)
  axis                which end differs: tss / tts / both

Outputs:
  {prefix}_alt_end_cds.tsv         per-isoform CDS / UTR features
  {prefix}_alt_end_cds_pairs.tsv   per-sibling-pair delta features
  {prefix}_alt_end_cds_panel.png   4-panel summary

Inputs:
  --bed label:path  ...  (30-column FlairBed *.isoforms.bed files; cols
                          7/8 are the ORF start/end (thickStart/thickEnd)
                          populated by predictProductivity; col 18 is the
                          dedicated productivity column PRO/PTC/NGO/NST;
                          col 13 is the FLG gene_id).
  --output-prefix TEXT
  --output-dir DIR

CDS coords and productivity are read INLINE from the regular
*.isoforms.bed (the separate *.isoforms.CDS.bed file no longer exists in
the new FLAIR/TED schema).  When thickStart/thickEnd equals start/end and
the productivity column is empty, the BED carries no CDS prediction and
the analysis logs a warning + skips the file.

Legacy fallback: if col 18 is absent (old plain-BED12 with productivity
embedded in the name as <iso>_<PRO|PTC|NST|NGO>_<gene>), the name regex
is still consulted so old files keep working.

Productivity classes (from predictProductivity.py):
  PRO  productive ORF
  PTC  premature termination codon (NMD-triggering)
  NST  no stop codon (read-through)
  NGO  no good ORF
"""
from __future__ import annotations

import argparse
import logging
import re
import sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Bring sibling helpers in via the existing bin/ path convention.
_THIS_DIR = Path(__file__).resolve().parent
_BIN_DIR = _THIS_DIR.parent
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))

from evaluation.signal_utils import group_by_junction_chain, gene_from_name  # noqa: E402

logging.basicConfig(format="%(levelname)s %(message)s", level=logging.INFO)
log = logging.getLogger(__name__)


# ── BED12 + CDS parser (self-contained — signal_utils.parse_bed12 doesn't expose thick fields) ──


_PRODUCTIVITY_RE = re.compile(r"_(PRO|PTC|NST|NGO)_")


def _parse_productivity(name: str) -> str | None:
    """Pull the productivity class from a name like 'iso_PRO_ENSG...'.
    Returns None if no class is embedded (i.e., predictProductivity
    didn't run, or this isn't a CDS BED)."""
    m = _PRODUCTIVITY_RE.search(name)
    return m.group(1) if m else None


def parse_bed12_cds(path: Path) -> list[dict]:
    """Read a 30-column FlairBed (or legacy BED12) with CDS columns.

    Returns one dict per isoform with:

        chrom, start, end, name, strand, gene_id,
        thick_start, thick_end,             (CDS genomic span, cols 7/8)
        productivity,                       (PRO/PTC/NST/NGO or None)
        tss, tts,                           (strand-aware end coords)
        block_starts, block_sizes,          (lists; needed for UTR length)
        n_exons.

    Productivity is read from the dedicated FlairBed column 18 (the new
    schema); gene_id from the FLG column 13.  For legacy plain BED12 files
    (no col 18) the productivity is recovered from the name regex and
    gene_id from :func:`gene_from_name`.

    UTR length isn't computed here — it requires walking the exon
    blocks to count only the portion of UTR contained in exons (not in
    introns). That's done in compute_utr_lengths().
    """
    isos: list[dict] = []
    if not path.exists():
        return isos
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "track")):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 12:
                continue
            chrom = c[0]
            start = int(c[1])
            end = int(c[2])
            name = c[3]
            strand = c[5]
            thick_start = int(c[6])
            thick_end = int(c[7])
            block_count = int(c[9])
            block_sizes = [int(x) for x in c[10].rstrip(",").split(",") if x]
            block_starts = [int(x) for x in c[11].rstrip(",").split(",") if x]
            tss = start if strand == "+" else end
            tts = end if strand == "+" else start

            # ── FlairBed inline fields: gene_id (col 13), productivity
            #    (col 18). Fall back to name-derived values for legacy BEDs.
            is_flair_bed = len(c) >= 18
            flg_gene_id = c[12].strip().rstrip(",").strip() if len(c) > 12 else ""
            prod_inline = c[17].strip().rstrip(",").strip() if len(c) > 17 else ""
            productivity = prod_inline or _parse_productivity(name)
            if productivity == "":
                productivity = None
            gene_id = flg_gene_id if (is_flair_bed and flg_gene_id) \
                else gene_from_name(name)

            # Build splice junctions for SJC grouping (matches signal_utils).
            juncs: list[tuple[int, int]] = []
            for i in range(block_count - 1):
                juncs.append(
                    (start + block_starts[i] + block_sizes[i], start + block_starts[i + 1])
                )
            isos.append(
                dict(
                    chrom=chrom,
                    start=start,
                    end=end,
                    name=name,
                    strand=strand,
                    gene_id=gene_id,
                    thick_start=thick_start,
                    thick_end=thick_end,
                    productivity=productivity,
                    tss=tss,
                    tts=tts,
                    block_starts=block_starts,
                    block_sizes=block_sizes,
                    n_exons=block_count,
                    junctions=tuple(juncs),
                    spliced_len=sum(block_sizes),
                )
            )
    return isos


def has_cds(iso: dict) -> bool:
    """True iff predictProductivity ran for this isoform.

    In the new FlairBed schema productivity is a dedicated column (18),
    populated with one of PRO/PTC/NGO/NST whenever predictProductivity
    ran — including NGO (no good ORF), where thickStart/thickEnd legitimately
    equal start/end.  So the presence of a productivity class is the
    authoritative "CDS was predicted" signal; we no longer require
    thick != start/end (that would wrongly drop NGO/NST isoforms).

    For the UTR-length computation, a meaningful coding span is additionally
    required — see :func:`has_coding_span`."""
    return iso["productivity"] is not None


def has_coding_span(iso: dict) -> bool:
    """True iff the isoform has a non-trivial coding span (thick != span).

    UTR lengths are only defined when there is an actual ORF inside the
    transcript; NGO isoforms (productivity set but thick==start/end) have
    no UTR partition."""
    if iso["productivity"] is None:
        return False
    return not (iso["thick_start"] == iso["start"] and iso["thick_end"] == iso["end"])


# ── UTR length on spliced transcript (sum exon-bp from TSS to CDS start) ─


def _spliced_distance_in_blocks(
    block_starts: list[int],
    block_sizes: list[int],
    start: int,
    pos_a: int,
    pos_b: int,
) -> int:
    """Spliced (exon-only) distance between two genomic positions
    pos_a < pos_b that lie within the BED12 blocks. Sums the bp of each
    block that falls in [pos_a, pos_b], skipping intronic gaps.

    Returns 0 if pos_a == pos_b. Returns -1 if either position is
    outside the BED's blocks (e.g. UTR computation failure).
    """
    if pos_a > pos_b:
        pos_a, pos_b = pos_b, pos_a
    total = 0
    found_a = found_b = False
    for bsz, bst in zip(block_sizes, block_starts):
        block_lo = start + bst
        block_hi = block_lo + bsz
        # Slice this block to the [pos_a, pos_b] window.
        seg_lo = max(block_lo, pos_a)
        seg_hi = min(block_hi, pos_b)
        if seg_hi > seg_lo:
            total += seg_hi - seg_lo
        if block_lo <= pos_a < block_hi or pos_a == block_hi:
            found_a = True
        if block_lo <= pos_b < block_hi or pos_b == block_hi:
            found_b = True
    if not (found_a and found_b):
        return -1
    return total


def compute_utr_lengths(iso: dict) -> tuple[int, int]:
    """Return (5utr_len, 3utr_len) on the spliced transcript.
    Returns (-1, -1) if there is no meaningful coding span (e.g. NGO)."""
    if not has_coding_span(iso):
        return -1, -1
    if iso["strand"] == "+":
        utr5_len = _spliced_distance_in_blocks(
            iso["block_starts"], iso["block_sizes"],
            iso["start"], iso["start"], iso["thick_start"],
        )
        utr3_len = _spliced_distance_in_blocks(
            iso["block_starts"], iso["block_sizes"],
            iso["start"], iso["thick_end"], iso["end"],
        )
    else:
        # Minus strand: 5'UTR is from end (TSS) to thick_end (CDS start
        # in genomic terms); 3'UTR from thick_start (CDS end genomic)
        # to start (TTS).
        utr5_len = _spliced_distance_in_blocks(
            iso["block_starts"], iso["block_sizes"],
            iso["start"], iso["thick_end"], iso["end"],
        )
        utr3_len = _spliced_distance_in_blocks(
            iso["block_starts"], iso["block_sizes"],
            iso["start"], iso["start"], iso["thick_start"],
        )
    return utr5_len, utr3_len


# ── Per-pair sibling features ──────────────────────────────────────


def _classify_pair_axis(a: dict, b: dict) -> str:
    """Which end differs between sibling isoforms a and b?
    Returns 'tss' | 'tts' | 'both' | 'identical'."""
    tss_diff = a["tss"] != b["tss"]
    tts_diff = a["tts"] != b["tts"]
    if tss_diff and tts_diff:
        return "both"
    if tss_diff:
        return "tss"
    if tts_diff:
        return "tts"
    return "identical"


def _productivity_transition(a: dict, b: dict) -> str:
    """Format a->b transition (or just a class if both same)."""
    pa = a["productivity"] or "?"
    pb = b["productivity"] or "?"
    if pa == pb:
        return pa
    return f"{pa}->{pb}"


# ── Driver ─────────────────────────────────────────────────────────


def analyse_one_label(label: str, isos: list[dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute per-isoform features and per-sibling-pair deltas for one
    method. Returns (per_iso_df, per_pair_df)."""
    n_total = len(isos)
    n_with_cds = sum(1 for i in isos if has_cds(i))
    log.info(
        "  %s: %d isoforms (%d with CDS predicted)",
        label,
        n_total,
        n_with_cds,
    )
    if n_with_cds == 0:
        log.warning(
            "  %s: no isoforms have CDS data; alt-end CDS analysis "
            "skipped for this label",
            label,
        )
        return pd.DataFrame(), pd.DataFrame()

    # Per-isoform table.
    per_iso_rows = []
    for iso in isos:
        utr5, utr3 = compute_utr_lengths(iso)
        per_iso_rows.append(
            dict(
                method=label,
                iso=iso["name"],
                gene_id=iso["gene_id"],
                chrom=iso["chrom"],
                strand=iso["strand"],
                tss=iso["tss"],
                tts=iso["tts"],
                cds_start=iso["thick_start"],
                cds_end=iso["thick_end"],
                productivity=iso["productivity"],
                utr5_len=utr5,
                utr3_len=utr3,
                n_exons=iso["n_exons"],
            )
        )
    per_iso_df = pd.DataFrame(per_iso_rows)

    # Per-pair table — only iterate SJC groups with >=2 multi-exon isoforms.
    groups = group_by_junction_chain(isos)
    per_pair_rows = []
    for jc, members in groups.items():
        if len(members) < 2:
            continue
        # Restrict to members with CDS data.
        members_cds = [m for m in members if has_cds(m)]
        if len(members_cds) < 2:
            continue
        for a, b in combinations(members_cds, 2):
            utr5_a, utr3_a = compute_utr_lengths(a)
            utr5_b, utr3_b = compute_utr_lengths(b)
            if utr5_a < 0 or utr3_a < 0 or utr5_b < 0 or utr3_b < 0:
                # block-walk failed for one of them; skip to avoid
                # poisoning the delta distribution.
                continue
            per_pair_rows.append(
                dict(
                    method=label,
                    gene_id=a["gene_id"],
                    iso_a=a["name"],
                    iso_b=b["name"],
                    axis=_classify_pair_axis(a, b),
                    delta_5utr_len=utr5_b - utr5_a,
                    delta_3utr_len=utr3_b - utr3_a,
                    productivity_a=a["productivity"],
                    productivity_b=b["productivity"],
                    productivity_transition=_productivity_transition(a, b),
                    productivity_changed=(a["productivity"] != b["productivity"]),
                )
            )
    per_pair_df = pd.DataFrame(per_pair_rows)
    return per_iso_df, per_pair_df


def plot_panel(
    per_iso: pd.DataFrame,
    per_pair: pd.DataFrame,
    path: Path,
) -> None:
    """4-panel summary:
        (1) Δ5'UTR-len histogram per method (alt-TSS pairs only)
        (2) Δ3'UTR-len histogram per method (alt-TTS pairs only)
        (3) productivity-transition stacked bar per method (% PRO->PRO,
            PRO->PTC, PRO->NGO, etc.)
        (4) by-axis pair count per method (tss / tts / both / identical)
    """
    if per_pair.empty:
        log.warning("per_pair is empty; skipping panel plot")
        return

    methods = sorted(per_pair["method"].unique())
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # (1) Δ5'UTR for alt-TSS pairs (axis in {tss, both}).
    ax = axes[0, 0]
    for method in methods:
        sub = per_pair[
            (per_pair["method"] == method) & per_pair["axis"].isin(["tss", "both"])
        ]
        if not len(sub):
            continue
        # Clip extreme tails for visual clarity.
        d = np.clip(sub["delta_5utr_len"], -1000, 1000)
        ax.hist(d, bins=40, alpha=0.5, label=f"{method} (n={len(sub)})", histtype="step", linewidth=1.5)
    ax.set_xlabel("Δ 5'UTR length (b - a, bp)")
    ax.set_ylabel("# sibling pairs")
    ax.set_title("Alt-TSS sibling pairs: Δ5'UTR length")
    ax.legend(fontsize=8)
    ax.axvline(0, color="grey", linewidth=0.5)

    # (2) Δ3'UTR for alt-TTS pairs (axis in {tts, both}).
    ax = axes[0, 1]
    for method in methods:
        sub = per_pair[
            (per_pair["method"] == method) & per_pair["axis"].isin(["tts", "both"])
        ]
        if not len(sub):
            continue
        d = np.clip(sub["delta_3utr_len"], -1000, 1000)
        ax.hist(d, bins=40, alpha=0.5, label=f"{method} (n={len(sub)})", histtype="step", linewidth=1.5)
    ax.set_xlabel("Δ 3'UTR length (b - a, bp)")
    ax.set_ylabel("# sibling pairs")
    ax.set_title("Alt-TTS sibling pairs: Δ3'UTR length")
    ax.legend(fontsize=8)
    ax.axvline(0, color="grey", linewidth=0.5)

    # (3) Productivity transitions stacked bar.
    ax = axes[1, 0]
    transitions_per_method = {}
    all_transitions = set()
    for method in methods:
        sub = per_pair[per_pair["method"] == method]
        if not len(sub):
            transitions_per_method[method] = Counter()
            continue
        c = Counter(sub["productivity_transition"])
        transitions_per_method[method] = c
        all_transitions.update(c.keys())
    # Order so PRO-only comes first, then transitions, then NGO/NST.
    ordered = sorted(
        all_transitions,
        key=lambda t: (
            0 if t == "PRO" else 1 if "->" in t else 2,
            t,
        ),
    )
    bottom = np.zeros(len(methods))
    for trans in ordered:
        vals = [transitions_per_method[m].get(trans, 0) for m in methods]
        ax.bar(methods, vals, bottom=bottom, label=trans, alpha=0.85)
        bottom += np.array(vals)
    ax.set_ylabel("# sibling pairs")
    ax.set_title("Productivity-class transitions in sibling pairs")
    ax.legend(fontsize=7, loc="upper right", ncol=2)
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right", fontsize=7)

    # (4) Pair counts by axis.
    ax = axes[1, 1]
    axis_order = ["tss", "tts", "both", "identical"]
    for axis_name in axis_order:
        vals = [
            len(per_pair[(per_pair["method"] == m) & (per_pair["axis"] == axis_name)])
            for m in methods
        ]
        ax.bar(
            np.arange(len(methods)),
            vals,
            label=axis_name,
            alpha=0.85,
            bottom=[
                sum(
                    len(
                        per_pair[
                            (per_pair["method"] == m) & (per_pair["axis"] == prev)
                        ]
                    )
                    for prev in axis_order[: axis_order.index(axis_name)]
                )
                for m in methods
            ],
        )
    ax.set_xticks(np.arange(len(methods)))
    ax.set_xticklabels(methods, rotation=20, ha="right", fontsize=7)
    ax.set_ylabel("# sibling pairs")
    ax.set_title("Pair count by which end differs")
    ax.legend(fontsize=7, loc="upper right")

    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ── CLI ────────────────────────────────────────────────────────────


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--bed", nargs="+", required=True)
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--output-dir", type=Path, default=Path("."))
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.verbose:
        log.setLevel(logging.DEBUG)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    per_iso_dfs: list[pd.DataFrame] = []
    per_pair_dfs: list[pd.DataFrame] = []

    for entry in args.bed:
        if ":" not in entry:
            log.warning("skipping malformed --bed entry %r", entry)
            continue
        label, path_str = entry.split(":", 1)
        path = Path(path_str)
        if not path.exists():
            log.warning("BED %s not found for label %s; skipping", path, label)
            continue
        isos = parse_bed12_cds(path)
        log.info("loaded %d isoforms for %s from %s", len(isos), label, path.name)
        per_iso, per_pair = analyse_one_label(label, isos)
        if not per_iso.empty:
            per_iso_dfs.append(per_iso)
        if not per_pair.empty:
            per_pair_dfs.append(per_pair)

    if not per_iso_dfs:
        log.warning("no labels produced CDS-enriched data; nothing to write")
        sys.exit(0)

    per_iso_all = pd.concat(per_iso_dfs, ignore_index=True)
    per_pair_all = (
        pd.concat(per_pair_dfs, ignore_index=True)
        if per_pair_dfs
        else pd.DataFrame()
    )

    prefix = args.output_prefix
    out = args.output_dir
    per_iso_all.to_csv(out / f"{prefix}_alt_end_cds.tsv", sep="\t", index=False)
    per_pair_all.to_csv(out / f"{prefix}_alt_end_cds_pairs.tsv", sep="\t", index=False)
    plot_panel(per_iso_all, per_pair_all, out / f"{prefix}_alt_end_cds_panel.png")

    log.info(
        "wrote: %d per-iso rows, %d per-pair rows across %d methods",
        len(per_iso_all),
        len(per_pair_all),
        per_iso_all["method"].nunique(),
    )


if __name__ == "__main__":
    main()
