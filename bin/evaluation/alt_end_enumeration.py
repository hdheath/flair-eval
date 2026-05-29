#!/usr/bin/env python3
"""
alt_end_enumeration.py — Per-gene #TSS / #PAS enumeration (Anvar-style An-F1).

For each gene, count the number of distinct transcription start sites (TSS) and
distinct polyadenylation sites (PAS/TTS) used by its isoforms, and report the
fraction of genes that are multi-TSS and/or multi-PAS. This is the per-gene
alternative-end *repertoire* — the foundation table the usage-fraction and
coupling analyses build on.

Distinct ends are counted after binning end coordinates to a tolerance window
(default 50 bp, the suite's ±50 bp convention) so that small coordinate wobble
within one true site does not inflate the count.

TED-confidence stratification (optional, --ted-log): when a flair_ted.ted_log
TSV is supplied, each gene's ends are split into high- vs low-confidence by the
per-end TED reality score (matched to isoform ends by position). The script then
reports the multi-TSS/multi-PAS fraction within the HIGH-confidence end set vs
the full set — the assay-free demonstration that high-confidence TED ends drive
a cleaner alternative-end repertoire. No competitor tool can draw this split.

Outputs (per method):
  <prefix>_alt_end_enumeration.tsv   per-gene n_tss / n_pas / multi flags
  <prefix>_alt_end_enumeration.png   histogram of #TSS and #PAS per gene
                                     + multi-end fraction bars (overall + by
                                     confidence bin when --ted-log given)

EMPTY INPUT: never exits non-zero. Writes a header-only TSV and skips the plot.

Usage:
    python alt_end_enumeration.py \\
        --bed flair:isoforms.bed ted:isoforms.bed ... \\
        [--ted-log flair:none ted:flair_ted.ted_log.tsv ...] \\
        --output-prefix run_name \\
        --output-dir out/ \\
        [--end-tol 50] [--min-exons 2] [--reality-bin 0.5]
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
    from signal_utils import parse_isoforms, gene_from_name
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
    from evaluation.signal_utils import parse_isoforms, gene_from_name

apply_rc()

DEFAULT_END_TOL   = 50    # bp; ends within this distance count as the same site
DEFAULT_MIN_EXONS = 2     # single-exon isoforms have no junction chain
DEFAULT_REALITY_BIN = 0.5  # TED reality >= this = "high-confidence" end

TSV_COLUMNS = [
    "method", "gene_id", "n_isoforms",
    "n_tss", "n_pas", "multi_tss", "multi_pas", "multi_end",
    "n_tss_highconf", "n_pas_highconf",
]


# ── Distinct-site counting ────────────────────────────────────────────────────

def count_distinct_sites(positions: List[int], tol: int) -> int:
    """Greedy single-linkage collapse: count clusters of positions within *tol*.

    Sorted positions; a new cluster starts whenever the gap to the previous
    position exceeds *tol*. Returns the number of clusters (= distinct sites).
    """
    if not positions:
        return 0
    sp = sorted(positions)
    n = 1
    last = sp[0]
    for p in sp[1:]:
        if p - last > tol:
            n += 1
        last = p
    return n


# ── ted_log per-end reality matching ──────────────────────────────────────────

def load_end_reality(ted_log_path: Path) -> Dict[Tuple[str, str], List[Tuple[int, float]]]:
    """Load per-end TED reality from a ted_log TSV.

    Returns {('tss'|'tts', chrom): [(pos, reality), ...]} for per-isoform rows
    only (global_peak partition-summary rows and -1 sentinels dropped).

    Imports the schema-aware loader so global-mode partition rows never leak in.
    Degrades to an empty dict when the file is missing/empty/unparseable.
    """
    out: Dict[Tuple[str, str], List[Tuple[int, float]]] = defaultdict(list)
    try:
        from ted_log_loader import load_ted_log
    except ImportError:
        try:
            from evaluation.ted_log_loader import load_ted_log
        except ImportError:
            return {}
    try:
        df = load_ted_log(ted_log_path, drop_partition_summary=True)
    except Exception:
        return {}
    if df is None or len(df) == 0:
        return {}
    for _, row in df.iterrows():
        junc = str(row.get("junc_id", ""))
        chrom = junc.split(":", 1)[0] if junc else ""
        for axis, pos_col, real_col in (
            ("tss", "tss_pos", "TED_tss_reality"),
            ("tts", "tts_pos", "TED_tts_reality"),
        ):
            pos = row.get(pos_col)
            real = row.get(real_col)
            try:
                pos = int(pos)
                real = float(real)
            except (TypeError, ValueError):
                continue
            if pos < 0 or real < 0:
                continue
            out[(axis, chrom)].append((pos, real))
    return dict(out)


def reality_at(axis: str, chrom: str, pos: int,
               reality_index: Dict[Tuple[str, str], List[Tuple[int, float]]],
               tol: int) -> Optional[float]:
    """Nearest-within-tol TED reality for an isoform end, or None if no match."""
    cands = reality_index.get((axis, chrom))
    if not cands:
        return None
    best_d, best_r = tol + 1, None
    for p, r in cands:
        d = abs(p - pos)
        if d <= tol and d < best_d:
            best_d, best_r = d, r
    return best_r


# ── Per-method enumeration ─────────────────────────────────────────────────────

def enumerate_method(
    isoforms: List[dict],
    end_tol: int,
    min_exons: int,
    reality_index: Optional[Dict[Tuple[str, str], List[Tuple[int, float]]]],
    reality_bin: float,
) -> List[dict]:
    """Return one row per gene with TSS/PAS counts (+ high-conf counts)."""
    genes: Dict[str, List[dict]] = defaultdict(list)
    for iso in isoforms:
        if int(iso.get("n_exons", 1)) < min_exons:
            continue
        gid = iso.get("gene_id") or gene_from_name(iso.get("name", ""))
        if not gid:
            continue
        genes[gid].append(iso)

    rows: List[dict] = []
    for gid, isos in genes.items():
        tss_pos = [int(i["tss"]) for i in isos]
        pas_pos = [int(i["tts"]) for i in isos]
        n_tss = count_distinct_sites(tss_pos, end_tol)
        n_pas = count_distinct_sites(pas_pos, end_tol)

        n_tss_hc = n_pas_hc = -1
        if reality_index:
            hc_tss, hc_pas = [], []
            for iso in isos:
                rt = reality_at("tss", iso["chrom"], int(iso["tss"]),
                                reality_index, end_tol)
                if rt is not None and rt >= reality_bin:
                    hc_tss.append(int(iso["tss"]))
                rp = reality_at("tts", iso["chrom"], int(iso["tts"]),
                                reality_index, end_tol)
                if rp is not None and rp >= reality_bin:
                    hc_pas.append(int(iso["tts"]))
            n_tss_hc = count_distinct_sites(hc_tss, end_tol)
            n_pas_hc = count_distinct_sites(hc_pas, end_tol)

        rows.append({
            "gene_id": gid,
            "n_isoforms": len(isos),
            "n_tss": n_tss,
            "n_pas": n_pas,
            "multi_tss": int(n_tss > 1),
            "multi_pas": int(n_pas > 1),
            "multi_end": int(n_tss > 1 or n_pas > 1),
            "n_tss_highconf": n_tss_hc,
            "n_pas_highconf": n_pas_hc,
        })
    return rows


# ── TSV ────────────────────────────────────────────────────────────────────────

def write_tsv(rows_by_method: Dict[str, List[dict]], path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(TSV_COLUMNS)
        for method, rows in rows_by_method.items():
            for r in sorted(rows, key=lambda x: x["gene_id"]):
                w.writerow([method] + [r[c] for c in TSV_COLUMNS[1:]])


# ── Plot ────────────────────────────────────────────────────────────────────────

def make_plot(rows_by_method: Dict[str, List[dict]], path: Path,
              has_reality: bool, reality_bin: float,
              title_prefix: str = "") -> None:
    methods = [m for m in rows_by_method if rows_by_method[m]]
    if not methods:
        return
    styler = ModeStyler(methods)

    fig, axes = plt.subplots(1, 3, figsize=(W2, W1 * 0.7))

    # Panel 1: #TSS per gene
    ax = axes[0]
    maxn = max((r["n_tss"] for m in methods for r in rows_by_method[m]), default=1)
    bins = np.arange(0.5, max(maxn, 1) + 1.5, 1)
    for m in methods:
        vals = [r["n_tss"] for r in rows_by_method[m]]
        ax.hist(vals, bins=bins, histtype="step", linewidth=1.3,
                color=styler.color(m), label=m)
    style_ax(ax, xlabel="# distinct TSS per gene", ylabel="genes")

    # Panel 2: #PAS per gene
    ax = axes[1]
    maxn = max((r["n_pas"] for m in methods for r in rows_by_method[m]), default=1)
    bins = np.arange(0.5, max(maxn, 1) + 1.5, 1)
    for m in methods:
        vals = [r["n_pas"] for r in rows_by_method[m]]
        ax.hist(vals, bins=bins, histtype="step", linewidth=1.3,
                color=styler.color(m), label=m)
    style_ax(ax, xlabel="# distinct PAS per gene", ylabel="genes")

    # Panel 3: multi-end fraction bars (overall; + high-conf when available)
    ax = axes[2]
    x = np.arange(len(methods))
    width = 0.35 if has_reality else 0.6
    frac_all = []
    frac_hc = []
    for m in methods:
        rows = rows_by_method[m]
        n = len(rows)
        frac_all.append(sum(r["multi_end"] for r in rows) / n if n else 0.0)
        if has_reality:
            hc = [(r["n_tss_highconf"], r["n_pas_highconf"]) for r in rows
                  if r["n_tss_highconf"] >= 0]
            if hc:
                frac_hc.append(
                    sum(1 for t, p in hc if t > 1 or p > 1) / len(hc))
            else:
                frac_hc.append(0.0)
    if has_reality:
        ax.bar(x - width / 2, frac_all, width,
               color=[styler.color(m) for m in methods], alpha=0.55,
               label="all ends")
        ax.bar(x + width / 2, frac_hc, width,
               color=[styler.color(m) for m in methods], alpha=1.0,
               edgecolor="black", linewidth=0.4,
               label=f"TED reality ≥ {reality_bin:g}")
        from matplotlib.patches import Patch
        ax.legend(handles=[
            Patch(facecolor="0.6", alpha=0.55, label="all ends"),
            Patch(facecolor="0.3", label=f"reality ≥ {reality_bin:g}")],
            fontsize=5, frameon=False, loc="upper right")
    else:
        ax.bar(x, frac_all, width, color=[styler.color(m) for m in methods])
    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=5, rotation=45, ha="right")
    ax.set_ylim(0, 1)
    style_ax(ax, ylabel="fraction multi-end genes")

    if title_prefix:
        fig.suptitle(f"{title_prefix}per-gene alternative-end repertoire",
                     fontsize=8)
    fig.tight_layout(pad=0.4)
    savefig(fig, path)


# ── CLI ────────────────────────────────────────────────────────────────────────

def _parse_pairs(entries: List[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for e in entries or []:
        if ":" not in e:
            print(f"WARNING: skipping malformed pair '{e}'", file=sys.stderr)
            continue
        label, path = e.split(":", 1)
        out[label] = path
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed", nargs="+", required=True,
                   help="label:path pairs for isoform BED12/GTF")
    p.add_argument("--ted-log", nargs="+", default=[],
                   help="label:path pairs for flair_ted.ted_log TSVs "
                        "(use 'none' or omit a label for non-TED methods)")
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--output-dir", default=".")
    p.add_argument("--end-tol", type=int, default=DEFAULT_END_TOL)
    p.add_argument("--min-exons", type=int, default=DEFAULT_MIN_EXONS)
    p.add_argument("--reality-bin", type=float, default=DEFAULT_REALITY_BIN)
    p.add_argument("--title-prefix", default="")
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args(argv)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = out_dir / f"{args.output_prefix}_alt_end_enumeration.tsv"
    png_path = out_dir / f"{args.output_prefix}_alt_end_enumeration.png"

    bed_pairs = _parse_pairs(args.bed)
    log_pairs = _parse_pairs(args.ted_log)

    rows_by_method: Dict[str, List[dict]] = {}
    has_reality = False
    for method, bed_path in bed_pairs.items():
        if not Path(bed_path).exists():
            print(f"WARNING: BED not found for {method}: {bed_path}",
                  file=sys.stderr)
            rows_by_method[method] = []
            continue
        isos = parse_isoforms(bed_path)

        reality_index = None
        lp = log_pairs.get(method)
        if lp and lp.lower() != "none" and Path(lp).exists():
            reality_index = load_end_reality(Path(lp))
            if reality_index:
                has_reality = True

        rows_by_method[method] = enumerate_method(
            isos, args.end_tol, args.min_exons, reality_index, args.reality_bin)
        if args.verbose:
            print(f"  {method}: {len(isos)} isoforms → "
                  f"{len(rows_by_method[method])} genes", file=sys.stderr)

    # Always write the TSV (header-only if everything empty) — never exit 1.
    write_tsv(rows_by_method, tsv_path)
    print(f"Wrote {tsv_path}")

    if any(rows_by_method.values()):
        make_plot(rows_by_method, png_path, has_reality, args.reality_bin,
                  args.title_prefix)
        print(f"Wrote {png_path}")
    else:
        print("No genes with multi-exon isoforms — skipping plot "
              "(empty input handled gracefully).", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
