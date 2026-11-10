#!/usr/bin/env python3
"""CDSCoverageRecovery: are TED-retained isoforms keeping their productive ORFs?

For each gene that survives in BOTH the baseline (FLAIR) and a TED mode,
compares the productivity-class composition of the gene's isoforms:

  retained_PRO    gene has >=1 PRO isoform in both baseline and mode
  lost_PRO        baseline had >=1 PRO; mode has 0 PRO (productivity gone)
  gained_PRO      baseline had 0 PRO; mode has >=1 PRO (rare; flagged)
  PRO_to_PTC      baseline-PRO genes where mode's only ORF is PTC
                  (premature termination, NMD-triggering)
  PRO_to_NGO      baseline-PRO genes where mode has no good ORF at all

Outputs:
  {prefix}_cds_coverage_recovery.tsv         per-(mode, gene) class
  {prefix}_cds_coverage_recovery_summary.tsv per-mode count summary
  {prefix}_cds_coverage_recovery_panel.png   stacked-bar summary

Inputs:
  --bed label:path  ...  (30-column FlairBed *.isoforms.bed files;
                          productivity read inline from col 18, gene_id
                          from col 13 — no separate *.isoforms.CDS.bed)
  --baseline-label TEXT  (default "FLAIR")
  --output-prefix TEXT
  --output-dir DIR

Productivity classes (from predictProductivity.py):
  PRO  productive ORF
  PTC  premature termination codon (NMD-triggering)
  NST  no stop codon (read-through)
  NGO  no good ORF
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

_THIS_DIR = Path(__file__).resolve().parent
_BIN_DIR = _THIS_DIR.parent
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))

# Reuse the CDS BED parser from alt_end_cds_analysis to keep behavior
# consistent (same productivity regex, same CDS-coord parse).
#
# NOTE on the CDS gate (hdheath decision 2026-05-29): coverage-recovery uses
# the CONSERVATIVE definition of "CDS predicted" — only isoforms with an
# ACTUAL CODING SPAN (thickStart != thickEnd, i.e. PRO/PTC) count.  NGO (no
# good ORF) and NST (start but no stop) have no real coding span (thick ==
# span) and are treated as NO-CDS, so they are dropped from the
# coverage-recovery diff.  That is what alt_end_cds_analysis.has_coding_span
# already means, so we import it under the local name ``has_cds`` here.  The
# productivity-present ``has_cds`` (which counts NGO/NST) is used only by the
# OTHER analyses (alt_end_cds productivity-transition classification).
from evaluation.alt_end_cds_analysis import (  # noqa: E402
    parse_bed12_cds,
    has_coding_span as has_cds,
)

logging.basicConfig(format="%(levelname)s %(message)s", level=logging.INFO)
log = logging.getLogger(__name__)


# ── Per-gene productivity composition ──────────────────────────────


def gene_productivity_classes(isos: list[dict]) -> dict[str, Counter]:
    """For each gene_id, count productivity classes across its isoforms.
    Isoforms without CDS data are skipped (counter stays zero on those)."""
    by_gene: dict[str, Counter] = defaultdict(Counter)
    for iso in isos:
        if not has_cds(iso):
            continue
        gid = iso.get("gene_id") or "_NO_GENE"
        cls = iso["productivity"]
        if cls is None:
            continue
        by_gene[gid][cls] += 1
    return dict(by_gene)


def classify_gene(
    base_counts: Counter, mode_counts: Counter
) -> str:
    """Classify how productivity composition changed for a gene."""
    base_pro = base_counts.get("PRO", 0)
    mode_pro = mode_counts.get("PRO", 0)
    base_total = sum(base_counts.values())
    mode_total = sum(mode_counts.values())

    if base_total == 0 and mode_total == 0:
        return "no_cds_either"
    if base_total == 0:
        # Mode has CDS data, baseline didn't. Mostly a parse anomaly;
        # bucket separately.
        return "no_baseline_cds"
    if mode_total == 0:
        # Baseline had CDS data, mode kept the gene but lost CDS data.
        # This means TED retained isoforms but predictProductivity
        # couldn't find an ORF — it's worth flagging.
        return "mode_lost_cds"

    if base_pro > 0 and mode_pro > 0:
        return "retained_PRO"
    if base_pro == 0 and mode_pro > 0:
        return "gained_PRO"
    if base_pro > 0 and mode_pro == 0:
        # PRO disappeared. What replaced it?
        if mode_counts.get("PTC", 0) > 0:
            return "PRO_to_PTC"
        if mode_counts.get("NGO", 0) > 0:
            return "PRO_to_NGO"
        if mode_counts.get("NST", 0) > 0:
            return "PRO_to_NST"
        return "lost_PRO"
    return "no_PRO_either"


# ── Driver ─────────────────────────────────────────────────────────


def diff_modes(
    counts_by_mode: dict[str, dict[str, Counter]], baseline_label: str
) -> pd.DataFrame:
    """Build per-(mode, gene) classification rows."""
    if baseline_label not in counts_by_mode:
        raise ValueError(
            f"baseline label {baseline_label!r} not in inputs; "
            f"available: {list(counts_by_mode.keys())}"
        )
    baseline = counts_by_mode[baseline_label]
    rows = []
    for mode, mode_counts in counts_by_mode.items():
        if mode == baseline_label:
            continue
        all_genes = set(baseline) | set(mode_counts)
        for gid in sorted(all_genes):
            b = baseline.get(gid, Counter())
            m = mode_counts.get(gid, Counter())
            rows.append(
                dict(
                    mode=mode,
                    gene_id=gid,
                    baseline_PRO=b.get("PRO", 0),
                    baseline_PTC=b.get("PTC", 0),
                    baseline_NGO=b.get("NGO", 0),
                    baseline_NST=b.get("NST", 0),
                    mode_PRO=m.get("PRO", 0),
                    mode_PTC=m.get("PTC", 0),
                    mode_NGO=m.get("NGO", 0),
                    mode_NST=m.get("NST", 0),
                    classification=classify_gene(b, m),
                )
            )
    return pd.DataFrame(rows)


def write_summary(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    summary = (
        df.groupby(["mode", "classification"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    cols_order = [
        "mode",
        "retained_PRO",
        "PRO_to_PTC",
        "PRO_to_NGO",
        "PRO_to_NST",
        "lost_PRO",
        "gained_PRO",
        "no_PRO_either",
        "mode_lost_cds",
        "no_baseline_cds",
        "no_cds_either",
    ]
    cols_present = [c for c in cols_order if c in summary.columns]
    summary = summary[cols_present]
    summary.to_csv(path, sep="\t", index=False)
    return summary


def plot_panel(
    summary: pd.DataFrame, baseline_label: str, path: Path
) -> None:
    """Stacked-bar: per mode, fraction of (mode, gene) rows in each
    classification bucket. Most informative buckets are colored
    distinctly; "no PRO either" / "no_cds_either" are greyed out."""
    if summary.empty:
        log.warning("summary empty; skipping panel plot")
        return

    # Order modes alphabetically (caller may have a preference; this
    # keeps the plot deterministic).
    summary_sorted = summary.sort_values("mode").reset_index(drop=True)
    methods = summary_sorted["mode"].tolist()
    n = len(methods)

    # Buckets and their colors. The order is meaningful — left of the
    # bar is "good" (retained_PRO), right is "lost" (PRO_to_*).
    bucket_specs = [
        ("retained_PRO",     "#2CA02C"),  # green
        ("gained_PRO",       "#98DF8A"),  # light green
        ("PRO_to_PTC",       "#D62728"),  # red
        ("PRO_to_NGO",       "#FF7F0E"),  # orange
        ("PRO_to_NST",       "#FFBB78"),  # peach
        ("lost_PRO",         "#9467BD"),  # purple
        ("no_PRO_either",    "#7F7F7F"),  # mid-grey
        ("mode_lost_cds",    "#C5B0D5"),  # pale purple
        ("no_baseline_cds",  "#C7C7C7"),  # light grey
        ("no_cds_either",    "#E0E0E0"),  # very light grey
    ]
    fig, ax = plt.subplots(figsize=(11, max(2.5, 0.35 * n + 1)))
    bottom = np.zeros(n)
    y = np.arange(n)
    for bucket, color in bucket_specs:
        if bucket not in summary_sorted.columns:
            continue
        vals = summary_sorted[bucket].astype(float).values
        if vals.sum() == 0:
            continue
        ax.barh(y, vals, left=bottom, label=bucket, color=color, alpha=0.92)
        bottom += vals

    ax.set_yticks(y)
    ax.set_yticklabels(methods, fontsize=7)
    ax.set_xlabel(f"# genes (vs {baseline_label})")
    ax.set_title(
        f"CDS Coverage Recovery: per-mode productivity-class transitions\n"
        f"(baseline = {baseline_label}; one row per (mode, gene))"
    )
    ax.invert_yaxis()
    ax.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.0, 0.5))
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
    p.add_argument(
        "--baseline-label",
        default="FLAIR",
        help="Label to treat as baseline (default: FLAIR)",
    )
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--output-dir", type=Path, default=Path("."))
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.verbose:
        log.setLevel(logging.DEBUG)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    counts_by_mode: dict[str, dict[str, Counter]] = {}
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
        n_with_cds = sum(1 for i in isos if has_cds(i))
        log.info(
            "  %s: %d isoforms (%d with CDS predicted)",
            label,
            len(isos),
            n_with_cds,
        )
        if n_with_cds == 0:
            log.warning(
                "  %s: no CDS data; skipping (was --predict_cds run?)",
                label,
            )
            continue
        counts_by_mode[label] = gene_productivity_classes(isos)

    if args.baseline_label not in counts_by_mode:
        log.warning(
            "baseline label %r not present (or has no CDS data); "
            "skipping analysis",
            args.baseline_label,
        )
        sys.exit(0)
    if len(counts_by_mode) < 2:
        log.warning("need at least 2 labels with CDS data to diff")
        sys.exit(0)

    df = diff_modes(counts_by_mode, args.baseline_label)
    log.info("diff: %d (mode, gene) rows across %d modes",
             len(df), df["mode"].nunique())

    prefix = args.output_prefix
    out = args.output_dir
    df.to_csv(out / f"{prefix}_cds_coverage_recovery.tsv", sep="\t", index=False)
    summary = write_summary(df, out / f"{prefix}_cds_coverage_recovery_summary.tsv")
    plot_panel(summary, args.baseline_label, out / f"{prefix}_cds_coverage_recovery_panel.png")

    log.info(
        "wrote: %d (mode, gene) rows; summary covers %d modes",
        len(df),
        len(summary),
    )


if __name__ == "__main__":
    main()
