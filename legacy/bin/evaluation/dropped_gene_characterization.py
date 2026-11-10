#!/usr/bin/env python3
"""DroppedGeneCharacterization: characterize the gene-level cost of TED.

Companion to GeneRecoveryAnalysis. Where GeneRecoveryAnalysis answers
"how many genes did each TED mode drop?", this script answers "what
KIND of genes are getting dropped?" — biotype, annotated transcript
length, isoform-repertoire complexity.

Driving question:
    Are TED-dropped genes systematically (a) too long for reads to span
    [length-driven], (b) sparsely annotated and likely non-coding
    [biotype-driven], or (c) something else?

Outputs per sample:
    {prefix}_dropped_gene_chars.tsv         per-(mode, gene) row with
                                            biotype + length + n_tx
    {prefix}_dropped_gene_summary.tsv       per-mode aggregate (drop rate
                                            per biotype bucket, length
                                            stats for dropped vs kept)
    {prefix}_dropped_gene_panel.png         2-panel:
                                              left  = drop rate by biotype
                                              right = length distribution
                                                      of dropped vs kept

Inputs:
    --bed label:path  ...  (BED12 isoform files; first label is treated
                            as baseline. Default expects 'FLAIR' first
                            but configurable via --baseline-label.)
    --gtf PATH             reference annotation (chrom-restricted ok)
    --output-prefix TEXT
    --output-dir DIR
    --baseline-label TEXT  (default 'FLAIR')

Caveats:
    Annotated transcript lengths are computed from the GTF the pipeline
    consumed. If the reference GTF doesn't have an entry for a gene, that
    gene's `n_annotated_tx=0`, which is itself an interesting signal
    (novel / un-annotated locus).
"""
from __future__ import annotations

import argparse
import logging
import re
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

from evaluation.signal_utils import parse_isoforms  # noqa: E402

logging.basicConfig(format="%(levelname)s %(message)s", level=logging.INFO)
log = logging.getLogger(__name__)


# ── GTF parsing — collect biotype + transcript lengths per gene ────


def _attr(attrs: str, key: str) -> str | None:
    m = re.search(rf'{key} "([^"]+)"', attrs)
    return m.group(1) if m else None


def parse_gtf_gene_features(gtf_path: Path) -> tuple[dict[str, str], dict[str, list[int]]]:
    """Return (gene_to_biotype, gene_to_transcript_lengths) by walking
    the GTF once. Transcript lengths are the sum of exon sizes per
    transcript_id; gene_to_transcript_lengths[gene] is a list of those
    sums for every transcript in that gene.

    Skips comment / malformed lines silently. Restricts to whatever
    chromosomes appear in the GTF — caller is expected to feed a
    chrom-filtered GTF if they want chr22-only stats."""
    tx_to_gene: dict[str, str] = {}
    tx_to_len: dict[str, int] = defaultdict(int)
    gene_to_biotype: dict[str, str] = {}
    if not gtf_path.exists():
        log.warning("GTF not found: %s — biotype/length features will be empty", gtf_path)
        return gene_to_biotype, defaultdict(list)
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            feature = cols[2]
            attrs = cols[8]
            if feature == "exon":
                tid = _attr(attrs, "transcript_id")
                gid = _attr(attrs, "gene_id")
                if not (tid and gid):
                    continue
                tx_to_gene[tid] = gid
                # GTF coords are 1-based inclusive
                tx_to_len[tid] += int(cols[4]) - (int(cols[3]) - 1)
            elif feature == "gene":
                gid = _attr(attrs, "gene_id")
                # gencode uses gene_type; ensembl/refseq use gene_biotype.
                # Fall back gracefully.
                bio = _attr(attrs, "gene_type") or _attr(attrs, "gene_biotype")
                if gid and bio:
                    gene_to_biotype[gid] = bio
    gene_to_lens: dict[str, list[int]] = defaultdict(list)
    for tid, gid in tx_to_gene.items():
        gene_to_lens[gid].append(tx_to_len[tid])
    return gene_to_biotype, gene_to_lens


# ── Per-(mode, gene) characterization ──────────────────────────────


def _strip_ensg_version(gid: str | None) -> str | None:
    """Strip an Ensembl version suffix (e.g. ``ENSG00000123456.7`` ->
    ``ENSG00000123456``). Leaves non-versioned / non-ENSG ids untouched
    except for trimming a trailing ``.<digits>`` so GTF and BED keys
    match regardless of which side carries the version."""
    if not gid:
        return gid
    return re.sub(r"\.\d+$", "", gid)


def _ref_gene_lookup_id(iso: dict) -> str:
    """Resolve the reference (ENSG) gene id to use for biotype/length
    lookups.

    For a 30-column FlairBed, ``gene_id`` is the FLAIR-internal ``FLG…``
    id (used for grouping) which never matches the GTF's ENSG keys; the
    real ENSG lives in ``ref_gene_mapping`` (col 15). Prefer that. For a
    legacy 12-column BED there is no ``ref_gene_mapping`` and ``gene_id``
    already holds the ENSG (or fixture gene id), so fall back to it."""
    ref = iso.get("ref_gene_mapping") or ""
    return ref if ref else iso.get("gene_id", "")


def classify_per_mode(
    isoforms_per_mode: dict[str, list[dict]],
    baseline_label: str,
    gene_to_biotype: dict[str, str],
    gene_to_lens: dict[str, list[int]],
) -> pd.DataFrame:
    """For each non-baseline mode, mark each gene from the baseline as
    dropped/kept and attach biotype + length features.

    Genes are grouped by ``gene_id`` (the FLAIR-internal ``FLG…`` id on a
    FlairBed, or the ENSG-from-name on a legacy BED). Biotype + length
    features, however, are keyed in the GTF by ENSG, so the lookup is
    performed against the reference gene mapping (``ref_gene_mapping``,
    col 15) with Ensembl version suffixes stripped on both sides.
    """
    if baseline_label not in isoforms_per_mode:
        raise ValueError(
            f"baseline label {baseline_label!r} not in inputs: "
            f"{list(isoforms_per_mode.keys())}"
        )

    # Version-strip the GTF-derived dicts once so the ENSG join is
    # robust to either side carrying a ``.<version>`` suffix.
    biotype_by_ensg = {
        _strip_ensg_version(g): b for g, b in gene_to_biotype.items()
    }
    lens_by_ensg = {
        _strip_ensg_version(g): l for g, l in gene_to_lens.items()
    }

    # Map each grouping gene_id -> the (version-stripped) ENSG to look up
    # in the GTF dicts. Built from the baseline isoforms.
    geneid_to_ensg: dict[str, str] = {}
    baseline_genes: dict[str, int] = Counter()
    for iso in isoforms_per_mode[baseline_label]:
        gid = iso["gene_id"]
        baseline_genes[gid] += 1
        geneid_to_ensg.setdefault(gid, _strip_ensg_version(_ref_gene_lookup_id(iso)))

    rows: list[dict] = []
    for mode, isos in isoforms_per_mode.items():
        if mode == baseline_label:
            continue
        mode_genes: dict[str, int] = Counter()
        for iso in isos:
            mode_genes[iso["gene_id"]] += 1
        # Iterate only genes present in baseline; the inverse direction
        # (genes the mode added that baseline lacked) is handled by
        # GeneRecoveryAnalysis already.
        for gene_id in baseline_genes:
            n_baseline = baseline_genes[gene_id]
            n_mode = mode_genes.get(gene_id, 0)
            status = "dropped" if n_mode == 0 else "kept"
            ensg = geneid_to_ensg.get(gene_id) or gene_id
            biotype = biotype_by_ensg.get(ensg, "unknown")
            lens = lens_by_ensg.get(ensg, [])
            rows.append(
                {
                    "mode": mode,
                    "gene_id": gene_id,
                    "status": status,
                    "biotype": biotype,
                    "n_baseline_iso": n_baseline,
                    "n_mode_iso": n_mode,
                    "n_annotated_tx": len(lens),
                    "annotated_mean_len_bp": (
                        float(np.mean(lens)) if lens else float("nan")
                    ),
                    "annotated_max_len_bp": (
                        int(max(lens)) if lens else 0
                    ),
                }
            )
    return pd.DataFrame(rows)


# ── Per-mode aggregate summary ─────────────────────────────────────


def _drop_rate(df: pd.DataFrame, where: pd.Series) -> float:
    sub = df[where]
    if len(sub) == 0:
        return 0.0
    return 100 * (sub["status"] == "dropped").sum() / len(sub)


def build_summary(detail: pd.DataFrame) -> pd.DataFrame:
    """One row per mode; columns aggregate biotype-level drop rates and
    length stats."""
    summaries: list[dict] = []
    for mode, sub in detail.groupby("mode", sort=False):
        n_total = len(sub)
        n_dropped = (sub["status"] == "dropped").sum()
        n_kept = n_total - n_dropped

        row = {
            "mode": mode,
            "n_baseline_genes": n_total,
            "n_dropped": int(n_dropped),
            "n_kept": int(n_kept),
            "drop_rate_overall_pct": (
                round(100 * n_dropped / n_total, 1) if n_total else 0.0
            ),
        }

        # Per-biotype drop rate (top biotypes by gene count).
        for biotype in (
            "protein_coding",
            "lncRNA",
            "transcribed_unprocessed_pseudogene",
            "processed_pseudogene",
            "unknown",
        ):
            mask = sub["biotype"] == biotype
            row[f"drop_rate_{biotype}_pct"] = round(
                _drop_rate(sub, mask), 1
            )
            row[f"n_{biotype}"] = int(mask.sum())

        # Length stats — dropped vs kept (mean + median bp).
        d_lens = sub.loc[sub["status"] == "dropped", "annotated_mean_len_bp"].dropna()
        k_lens = sub.loc[sub["status"] == "kept", "annotated_mean_len_bp"].dropna()
        row["dropped_median_len_bp"] = (
            int(d_lens.median()) if len(d_lens) else 0
        )
        row["kept_median_len_bp"] = (
            int(k_lens.median()) if len(k_lens) else 0
        )
        row["dropped_p90_len_bp"] = (
            int(d_lens.quantile(0.9)) if len(d_lens) else 0
        )
        row["kept_p90_len_bp"] = (
            int(k_lens.quantile(0.9)) if len(k_lens) else 0
        )

        # Iso-repertoire stats — n annotated tx per gene.
        d_ntx = sub.loc[sub["status"] == "dropped", "n_annotated_tx"]
        k_ntx = sub.loc[sub["status"] == "kept", "n_annotated_tx"]
        row["dropped_median_n_tx"] = (
            int(d_ntx.median()) if len(d_ntx) else 0
        )
        row["kept_median_n_tx"] = (
            int(k_ntx.median()) if len(k_ntx) else 0
        )

        summaries.append(row)
    return pd.DataFrame(summaries)


# ── Plot ───────────────────────────────────────────────────────────


def plot_panel(detail: pd.DataFrame, out_path: Path) -> None:
    """Two-panel summary:
      Left  — per-biotype drop rate, grouped bar (one bar per mode per biotype).
      Right — annotated mean length distribution, dropped vs kept (log-x ECDF).
    """
    if detail.empty:
        log.warning("detail df empty; skipping panel")
        return

    modes = sorted(detail["mode"].unique())
    biotypes = (
        detail.groupby("biotype").size().sort_values(ascending=False).index.tolist()
    )
    # Cap to top 6 biotypes for readability.
    if len(biotypes) > 6:
        biotypes = biotypes[:6]

    fig, axes = plt.subplots(1, 2, figsize=(14, max(4, 0.4 * len(biotypes) + 1.5)))

    # --- Panel 1: drop rate by biotype × mode (horizontal grouped bar) ---
    ax = axes[0]
    bar_h = 0.8 / max(1, len(modes))
    y_centers = np.arange(len(biotypes))
    for i, mode in enumerate(modes):
        sub = detail[detail["mode"] == mode]
        rates = []
        for b in biotypes:
            mask = sub["biotype"] == b
            n = int(mask.sum())
            n_drop = int((sub.loc[mask, "status"] == "dropped").sum())
            rates.append(100 * n_drop / n if n else 0)
        offset = (i - (len(modes) - 1) / 2) * bar_h
        bars = ax.barh(
            y_centers + offset,
            rates,
            height=bar_h * 0.9,
            label=mode,
            alpha=0.88,
        )
        # Annotate counts inside bars to preview the n behind each %.
        for bar, b in zip(bars, biotypes):
            mask_full = sub["biotype"] == b
            n_total = int(mask_full.sum())
            n_drop = int((sub.loc[mask_full, "status"] == "dropped").sum())
            if n_total > 0:
                ax.text(
                    bar.get_width() + 1,
                    bar.get_y() + bar.get_height() / 2,
                    f"{n_drop}/{n_total}",
                    va="center",
                    fontsize=6,
                    color="dimgrey",
                )
    ax.set_yticks(y_centers)
    ax.set_yticklabels(biotypes, fontsize=7)
    ax.set_xlabel("drop rate (%)")
    ax.set_title("TED-dropped fraction by gene biotype")
    ax.set_xlim(0, max(105, ax.get_xlim()[1] + 5))
    ax.legend(fontsize=6, loc="lower right", framealpha=0.85)
    ax.invert_yaxis()

    # --- Panel 2: length ECDF, dropped vs kept (per mode) ---
    ax = axes[1]
    for mode in modes:
        sub = detail[detail["mode"] == mode]
        d_lens = sub.loc[sub["status"] == "dropped", "annotated_mean_len_bp"].dropna().values
        k_lens = sub.loc[sub["status"] == "kept", "annotated_mean_len_bp"].dropna().values
        if len(d_lens):
            d_sorted = np.sort(d_lens)
            ax.plot(
                d_sorted,
                np.arange(1, len(d_sorted) + 1) / len(d_sorted),
                linestyle="--",
                label=f"{mode} (dropped, n={len(d_lens)})",
                alpha=0.85,
            )
        if len(k_lens):
            k_sorted = np.sort(k_lens)
            ax.plot(
                k_sorted,
                np.arange(1, len(k_sorted) + 1) / len(k_sorted),
                linestyle="-",
                label=f"{mode} (kept, n={len(k_lens)})",
                alpha=0.85,
            )
    ax.set_xscale("log")
    ax.set_xlabel("annotated mean transcript length (bp, log scale)")
    ax.set_ylabel("ECDF")
    ax.set_title("Length distribution of dropped vs kept genes")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6, loc="lower right", framealpha=0.85)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
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
        "--gtf",
        type=Path,
        required=True,
        help="Reference GTF (with 'gene' rows carrying gene_type/gene_biotype, "
        "and 'exon' rows for transcript length aggregation).",
    )
    p.add_argument(
        "--baseline-label",
        default="FLAIR",
        help="Label whose genes serve as the baseline (default: FLAIR)",
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

    log.info("Parsing GTF for biotype + transcript-length features...")
    gene_to_biotype, gene_to_lens = parse_gtf_gene_features(args.gtf)
    log.info(
        "  loaded biotype for %d genes; transcripts for %d genes",
        len(gene_to_biotype),
        len(gene_to_lens),
    )

    isoforms_per_mode: dict[str, list[dict]] = {}
    for entry in args.bed:
        if ":" not in entry:
            log.warning("skipping malformed --bed entry %r", entry)
            continue
        label, path_str = entry.split(":", 1)
        path = Path(path_str)
        if not path.exists():
            log.warning("BED %s not found for label %s; skipping", path, label)
            continue
        isos = parse_isoforms(path)
        log.info("  %s: %d isoforms", label, len(isos))
        isoforms_per_mode[label] = isos

    if args.baseline_label not in isoforms_per_mode:
        log.error(
            "baseline label %r not present; available: %s",
            args.baseline_label,
            list(isoforms_per_mode.keys()),
        )
        sys.exit(0)
    if len(isoforms_per_mode) < 2:
        log.warning("need at least 2 labels (baseline + 1 mode) — skipping")
        sys.exit(0)

    detail = classify_per_mode(
        isoforms_per_mode, args.baseline_label, gene_to_biotype, gene_to_lens
    )
    if detail.empty:
        log.warning("no rows produced — exiting cleanly")
        sys.exit(0)

    summary = build_summary(detail)

    prefix = args.output_prefix
    out = args.output_dir
    detail.to_csv(out / f"{prefix}_dropped_gene_chars.tsv", sep="\t", index=False)
    summary.to_csv(out / f"{prefix}_dropped_gene_summary.tsv", sep="\t", index=False)
    plot_panel(detail, out / f"{prefix}_dropped_gene_panel.png")

    log.info(
        "wrote: %d (mode, gene) detail rows; summary covers %d modes",
        len(detail),
        len(summary),
    )


if __name__ == "__main__":
    main()
