#!/usr/bin/env python3
"""GeneRecoveryAnalysis: per-gene isoform-count diff vs a baseline.

Answers two precision/recall questions about TED's effect on gene-level
output:

  1. Drop-entirely: which genes does TED reject ALL isoforms for?
     (If FLAIR baseline assembled isoforms for gene G but no TED mode
     keeps any, G is "dropped" by that mode.)

  2. Alt-isoform collapse: which genes go from N>=2 FLAIR isoforms to
     exactly 1 TED isoform? (TED's peak-snapping effectively merged the
     alt-end siblings into one.)

For each non-baseline mode, computes:
    n_genes_dropped     genes with >=1 FLAIR iso, 0 in this mode
    n_genes_collapsed   genes with >=2 FLAIR isos, exactly 1 in this mode
    n_genes_expanded    this mode has > FLAIR's count (rare; recorded)
    n_genes_unchanged   counts match
    n_genes_partial     mode has 0<count<FLAIR (some isoforms lost,
                        but not all and not collapsed-to-1)

Outputs:
    {prefix}_gene_recovery_dropped.tsv     per-mode list of dropped genes
    {prefix}_gene_recovery_collapsed.tsv   per-mode list of collapsed genes
    {prefix}_gene_recovery_summary.tsv     per-mode count summary
    {prefix}_gene_recovery_panel.png       2-panel summary bar chart

Inputs:
    --bed label:path  ...  (BED12 isoform files, one per mode/method)
    --baseline-label TEXT  (default "FLAIR" — the no-TED baseline)
    --output-prefix TEXT
    --output-dir DIR

Conventions match UTRFeatures (modules/divergence/utr_features/main.nf).
Reuses signal_utils.parse_isoforms for BED parsing — avoids hand-rolling
the same BED12 logic.
"""
from __future__ import annotations

import argparse
import logging
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Resolve sibling helpers without each test/script hard-coding sys.path
# manipulations: walk up to repo root, add bin/ once.
_THIS_DIR = Path(__file__).resolve().parent
_BIN_DIR = _THIS_DIR.parent
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))

from evaluation.signal_utils import parse_isoforms  # noqa: E402

logging.basicConfig(format="%(levelname)s %(message)s", level=logging.INFO)
log = logging.getLogger(__name__)


# ── Per-mode iso → gene counts ─────────────────────────────────────


def count_isoforms_per_gene(isoforms: list[dict]) -> dict[str, int]:
    """For each gene_id, count how many isoforms reference it.

    Isoforms with no detectable gene_id are bucketed under "_NO_GENE_ID"
    so we don't silently drop them (that bucket should be small for
    flair output; non-trivial size means a parse drift).
    """
    counts: Counter[str] = Counter()
    no_gene = 0
    for iso in isoforms:
        gid = iso.get("gene_id")
        if not gid:
            no_gene += 1
            counts["_NO_GENE_ID"] += 1
            continue
        counts[gid] += 1
    if no_gene:
        log.info(
            "  %d isoforms had no detectable gene_id "
            "(bucketed under _NO_GENE_ID)",
            no_gene,
        )
    return dict(counts)


# ── Diff against baseline ──────────────────────────────────────────


def classify_gene(baseline_n: int, mode_n: int) -> str:
    """Return one of: dropped, collapsed, partial, unchanged, expanded,
    not_in_baseline. The classification is asymmetric — "FLAIR has N,
    this mode has M" — and pinned to specific definitions:

      not_in_baseline : baseline has 0 isoforms for this gene; the mode
                        invented one (TED finds a novel gene? more
                        likely a gene_id parse drift; flagged here).
      dropped         : baseline >= 1, mode == 0.
      collapsed       : baseline >= 2, mode == 1.
      partial         : baseline > mode > 0 (some lost, not collapsed).
      unchanged       : baseline == mode.
      expanded        : mode > baseline (TED found more isoforms than
                        FLAIR — surprising; flagged but rare).
    """
    if baseline_n == 0 and mode_n > 0:
        return "not_in_baseline"
    if baseline_n >= 1 and mode_n == 0:
        return "dropped"
    if baseline_n >= 2 and mode_n == 1:
        return "collapsed"
    if baseline_n == mode_n:
        return "unchanged"
    if mode_n > baseline_n:
        return "expanded"
    # baseline > mode > 0 and not collapsed
    return "partial"


def diff_modes(
    counts_by_mode: dict[str, dict[str, int]], baseline_label: str
) -> pd.DataFrame:
    """Build a long-form DataFrame: one row per (mode, gene_id) showing
    the baseline count, mode count, and classification. Modes are
    everything in counts_by_mode except `baseline_label`."""
    if baseline_label not in counts_by_mode:
        raise ValueError(
            f"baseline label {baseline_label!r} not present in BED inputs; "
            f"available: {list(counts_by_mode.keys())}"
        )
    baseline = counts_by_mode[baseline_label]
    rows = []
    for mode, mode_counts in counts_by_mode.items():
        if mode == baseline_label:
            continue
        all_genes = set(baseline) | set(mode_counts)
        for gid in sorted(all_genes):
            b_n = baseline.get(gid, 0)
            m_n = mode_counts.get(gid, 0)
            rows.append(
                {
                    "mode": mode,
                    "gene_id": gid,
                    "baseline_n": b_n,
                    "mode_n": m_n,
                    "classification": classify_gene(b_n, m_n),
                }
            )
    return pd.DataFrame(rows)


# ── Outputs ────────────────────────────────────────────────────────


def write_dropped_tsv(df: pd.DataFrame, path: Path) -> int:
    sub = df[df["classification"] == "dropped"]
    sub = sub[["mode", "gene_id", "baseline_n", "mode_n"]].sort_values(
        ["mode", "baseline_n"], ascending=[True, False]
    )
    sub.to_csv(path, sep="\t", index=False)
    return len(sub)


def write_collapsed_tsv(df: pd.DataFrame, path: Path) -> int:
    sub = df[df["classification"] == "collapsed"]
    sub = sub[["mode", "gene_id", "baseline_n", "mode_n"]].sort_values(
        ["mode", "baseline_n"], ascending=[True, False]
    )
    sub.to_csv(path, sep="\t", index=False)
    return len(sub)


def write_summary_tsv(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    """Per-mode summary: count of each classification."""
    summary = (
        df.groupby(["mode", "classification"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    # Order columns consistently
    cols_order = [
        "mode",
        "dropped",
        "collapsed",
        "partial",
        "unchanged",
        "expanded",
        "not_in_baseline",
    ]
    cols_present = [c for c in cols_order if c in summary.columns]
    summary = summary[cols_present]
    summary.to_csv(path, sep="\t", index=False)
    return summary


def plot_panel(summary: pd.DataFrame, baseline_label: str, path: Path) -> None:
    """2-panel horizontal bar chart: dropped count and collapsed count
    per mode. Mode labels go on y-axis, count on x-axis. Modes are
    ordered by dropped count ascending (least-disruptive at top)."""
    if summary.empty:
        log.warning("summary is empty; skipping panel plot")
        return

    summary_sorted = summary.sort_values("dropped").reset_index(drop=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, max(2, 0.3 * len(summary_sorted))))

    y = np.arange(len(summary_sorted))

    axes[0].barh(y, summary_sorted["dropped"], color="#D55E00", alpha=0.85)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(summary_sorted["mode"], fontsize=7)
    axes[0].set_xlabel(f"# genes dropped (vs {baseline_label})")
    axes[0].set_title("Drop-entirely: 0 isoforms passed")
    axes[0].invert_yaxis()

    axes[1].barh(y, summary_sorted["collapsed"], color="#0072B2", alpha=0.85)
    axes[1].set_yticks(y)
    axes[1].set_yticklabels([])  # share with left
    axes[1].set_xlabel(f"# genes collapsed (vs {baseline_label})")
    axes[1].set_title("Alt-isoform collapse: N>=2 -> 1")

    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ── CLI ────────────────────────────────────────────────────────────


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--bed",
        nargs="+",
        required=True,
        help="One or more 'label:path' pairs, e.g. FLAIR:flair.bed TED:ted.bed",
    )
    p.add_argument(
        "--baseline-label",
        default="FLAIR",
        help="Which label to treat as baseline (default: FLAIR)",
    )
    p.add_argument(
        "--output-prefix",
        required=True,
        help="Output filename prefix",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Output directory (default: current)",
    )
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.verbose:
        log.setLevel(logging.DEBUG)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    counts_by_mode: dict[str, dict[str, int]] = {}
    for entry in args.bed:
        if ":" not in entry:
            log.warning("skipping malformed --bed entry %r (expected label:path)", entry)
            continue
        label, path_str = entry.split(":", 1)
        path = Path(path_str)
        if not path.exists():
            log.warning("BED %s not found for label %s; skipping", path, label)
            continue
        isos = parse_isoforms(path)
        log.info("loaded %d isoforms for %s", len(isos), label)
        counts_by_mode[label] = count_isoforms_per_gene(isos)

    if args.baseline_label not in counts_by_mode:
        log.error(
            "baseline label %r not present; available: %s. Skipping analysis.",
            args.baseline_label,
            list(counts_by_mode.keys()),
        )
        sys.exit(0)

    if len(counts_by_mode) < 2:
        log.warning("only one method present; nothing to diff against baseline")
        sys.exit(0)

    df = diff_modes(counts_by_mode, args.baseline_label)
    log.info("diff: %d (mode, gene) rows across %d modes", len(df), df["mode"].nunique())

    prefix = args.output_prefix
    out = args.output_dir
    n_dropped = write_dropped_tsv(df, out / f"{prefix}_gene_recovery_dropped.tsv")
    n_collapsed = write_collapsed_tsv(df, out / f"{prefix}_gene_recovery_collapsed.tsv")
    summary = write_summary_tsv(df, out / f"{prefix}_gene_recovery_summary.tsv")
    plot_panel(summary, args.baseline_label, out / f"{prefix}_gene_recovery_panel.png")

    log.info(
        "wrote: %d dropped, %d collapsed rows; summary covers %d modes",
        n_dropped,
        n_collapsed,
        len(summary),
    )


if __name__ == "__main__":
    main()
