#!/usr/bin/env python3
"""Cell-line precision support panel.

Builds one three-column figure per cell line:

  1. 5' GTF-reference precision vs CAGE precision
  2. 3' GTF-reference precision vs dRNA precision
  3. GTF splice-junction precision vs STAR splice-junction precision

Each point is one sample x transcriptome method.  STAR SJ support is an exact
splice-junction match to SJ.out.tab entries with unique_reads >= 1.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pub_style import ModeStyler, legend_outside, savefig, style_ax
from signal_utils import parse_isoforms


STAR_STRAND = {"1": "+", "2": "-"}
_GTF_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def _split_label_path(entry: str, expected_parts: int) -> list[str] | None:
    parts = entry.split(":", expected_parts - 1)
    if len(parts) != expected_parts:
        print(f"WARNING: malformed label:path entry skipped: {entry}", file=sys.stderr)
        return None
    return parts


def _sample_from_pr_path(path: Path, test_name: str | None) -> str:
    stem = path.name
    for suffix in (
        "_precision_recall_ortho.tsv",
        "_precision_recall_gtf.tsv",
        "_precision_recall.tsv",
    ):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    if test_name and stem.startswith(f"{test_name}_"):
        return stem[len(test_name) + 1 :]
    return stem


def _cell_line(sample_id: str, star_labels: Iterable[str] = ()) -> str:
    upper = sample_id.upper()
    for label in star_labels:
        if upper.startswith(str(label).upper()):
            return str(label)
    return sample_id.split("_", 1)[0]


def _sample_short(sample_id: str, cell_line: str) -> str:
    label = sample_id
    prefix = f"{cell_line}_"
    if label.startswith(prefix):
        label = label[len(prefix) :]
    return label.replace("_", "-")


def _short_mode(name: str) -> str:
    return (
        str(name)
        .replace("lib-", "")
        .replace("-TED-global-logreg-05-14", "-global")
        .replace("TED-per-SJ", "TED-SJ")
        .replace("isoquant_auto", "isoquant")
        .replace("_", "-")
    )


def _mode_order(modes: Iterable[str]) -> list[str]:
    preferred = [
        "FLAIR",
        "TED-per-SJ",
        "isoquant_auto",
        "lib-pacbio-cDNA-TED-global-logreg-05-14",
        "lib-ont-cDNA-TED-global-logreg-05-14",
        "lib-ont-dRNA-TED-global-logreg-05-14",
    ]
    present = [str(m) for m in modes if pd.notna(m)]
    seen = set()
    present = [m for m in present if not (m in seen or seen.add(m))]
    out = [m for m in preferred if m in seen]
    out.extend(sorted(m for m in seen if m not in out))
    return out


def load_pr_tables(paths: list[str], *, kind: str, test_name: str | None) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for raw in paths:
        path = Path(raw)
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception as exc:
            print(f"WARNING: could not read {path}: {exc}", file=sys.stderr)
            continue
        if df.empty:
            continue
        sample_id = (
            str(df["dataset_name"].dropna().iloc[0])
            if "dataset_name" in df.columns and df["dataset_name"].notna().any()
            else _sample_from_pr_path(path, test_name)
        )
        df["sample_id"] = sample_id
        df["metric_source"] = kind
        df["source_file"] = str(path)
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _junction_set_from_isoforms(path: str | Path) -> set[tuple[str, int, int, str]]:
    junctions: set[tuple[str, int, int, str]] = set()
    for iso in parse_isoforms(path):
        chrom = iso.get("chrom")
        strand = iso.get("strand")
        if not chrom or strand not in {"+", "-"}:
            continue
        for donor, acceptor in iso.get("junctions", ()):
            junctions.add((chrom, int(donor), int(acceptor), strand))
    return junctions


def _transcript_id(attrs: str) -> str | None:
    parsed = dict(_GTF_ATTR_RE.findall(attrs))
    tid = parsed.get("transcript_id")
    if tid:
        return tid
    for part in attrs.strip().rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            if key in {"transcript_id", "Parent"}:
                return value.split(",", 1)[0].strip().strip('"')
    return None


def _junction_set_from_gtf(
    path: str | Path,
    chroms: set[str] | None = None,
) -> set[tuple[str, int, int, str]]:
    tx_exons: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            chrom = cols[0]
            if chroms and chrom not in chroms:
                continue
            strand = cols[6]
            if strand not in {"+", "-"}:
                continue
            tid = _transcript_id(cols[8])
            if not tid:
                continue
            try:
                start = int(cols[3]) - 1
                end = int(cols[4])
            except ValueError:
                continue
            tx_exons[(chrom, strand, tid)].append((start, end))

    junctions: set[tuple[str, int, int, str]] = set()
    for (chrom, strand, _tid), exons in tx_exons.items():
        if len(exons) < 2:
            continue
        exons.sort()
        for (_, donor), (acceptor, _) in zip(exons, exons[1:]):
            junctions.add((chrom, int(donor), int(acceptor), strand))
    return junctions


def _junction_set_from_star(
    path: str | Path,
    min_unique_reads: int,
    chroms: set[str] | None = None,
) -> set[tuple[str, int, int, str]]:
    junctions: set[tuple[str, int, int, str]] = set()
    p = Path(path)
    if not p.exists():
        print(f"WARNING: STAR SJ file not found: {path}", file=sys.stderr)
        return junctions
    with open(p) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            if chroms and cols[0] not in chroms:
                continue
            strand = STAR_STRAND.get(cols[3])
            if strand is None:
                continue
            try:
                unique_reads = int(cols[6])
                if unique_reads < min_unique_reads:
                    continue
                donor = int(cols[1]) - 1
                acceptor = int(cols[2])
            except ValueError:
                continue
            junctions.add((cols[0], donor, acceptor, strand))
    return junctions


def compute_sj_precision(
    isoform_entries: list[str],
    gtf_entries: list[str],
    star_entries: list[str],
    *,
    min_unique_reads: int,
) -> pd.DataFrame:
    gtf_by_sample: dict[str, str] = {}
    for entry in gtf_entries:
        parts = _split_label_path(entry, 2)
        if not parts:
            continue
        sample_id, path = parts
        gtf_by_sample[sample_id] = path

    star_paths: dict[str, str] = {}
    for entry in star_entries:
        parts = _split_label_path(entry, 2)
        if parts:
            star_paths[parts[0]] = parts[1]

    parsed_entries = []
    gtf_chroms: dict[str, set[str]] = defaultdict(set)
    star_chroms: dict[str, set[str]] = defaultdict(set)
    for entry in isoform_entries:
        parts = _split_label_path(entry, 3)
        if not parts:
            continue
        sample_id, mode, isoform_path = parts
        cell_line = _cell_line(sample_id, star_paths.keys())
        ref_gtf = gtf_by_sample.get(sample_id)
        if not ref_gtf:
            print(f"WARNING: no GTF supplied for sample {sample_id}", file=sys.stderr)
            continue

        called = _junction_set_from_isoforms(isoform_path)
        called_chroms = {j[0] for j in called}
        gtf_chroms[ref_gtf].update(called_chroms)
        star_chroms[cell_line].update(called_chroms)
        parsed_entries.append((sample_id, mode, cell_line, ref_gtf, called))

    gtf_cache = {
        ref_gtf: _junction_set_from_gtf(ref_gtf, chroms=chroms)
        for ref_gtf, chroms in gtf_chroms.items()
    }
    star_cache = {
        cell_line: _junction_set_from_star(
            path,
            min_unique_reads,
            chroms=star_chroms.get(cell_line) or None,
        )
        for cell_line, path in star_paths.items()
    }

    rows = []
    for sample_id, mode, cell_line, ref_gtf, called in parsed_entries:
        ref = gtf_cache.get(ref_gtf, set())
        star = star_cache.get(cell_line, set())
        n_called = len(called)
        n_ref = len(called & ref)
        n_star = len(called & star)
        rows.append({
            "sample_id": sample_id,
            "cell_line": cell_line,
            "transcriptome_mode": mode,
            "sj_n_called": n_called,
            "sj_n_gtf_supported": n_ref,
            "sj_n_star_supported": n_star,
            "sj_gtf_precision": n_ref / n_called if n_called else np.nan,
            "sj_star_precision": n_star / n_called if n_called else np.nan,
        })
    return pd.DataFrame(rows)


def build_summary(
    ortho_df: pd.DataFrame,
    gtf_df: pd.DataFrame,
    sj_df: pd.DataFrame,
    *,
    star_labels: Iterable[str],
) -> pd.DataFrame:
    metric_cols = [
        "sample_id", "transcriptome_mode",
        "5prime_precision", "3prime_precision",
    ]
    ortho = ortho_df[[c for c in metric_cols if c in ortho_df.columns]].copy()
    ref = gtf_df[[c for c in metric_cols if c in gtf_df.columns]].copy()
    ortho = ortho.rename(columns={
        "5prime_precision": "5prime_cage_precision",
        "3prime_precision": "3prime_drna_precision",
    })
    ref = ref.rename(columns={
        "5prime_precision": "5prime_ref_precision",
        "3prime_precision": "3prime_ref_precision",
    })
    merged = pd.merge(
        ref,
        ortho,
        on=["sample_id", "transcriptome_mode"],
        how="outer",
    )
    if not sj_df.empty:
        merged = pd.merge(
            merged,
            sj_df,
            on=["sample_id", "transcriptome_mode"],
            how="outer",
        )
    if "cell_line" not in merged.columns:
        merged["cell_line"] = [
            _cell_line(str(sample), star_labels) for sample in merged["sample_id"]
        ]
    else:
        merged["cell_line"] = merged["cell_line"].fillna(
            merged["sample_id"].map(lambda s: _cell_line(str(s), star_labels))
        )

    for col in [
        "5prime_ref_precision", "5prime_cage_precision",
        "3prime_ref_precision", "3prime_drna_precision",
        "sj_gtf_precision", "sj_star_precision",
    ]:
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce")
    return merged


def _plot_panel(
    ax,
    df: pd.DataFrame,
    *,
    xcol: str,
    ycol: str,
    xlabel: str,
    ylabel: str,
    title: str,
    mode_order: list[str],
    styler: ModeStyler,
    cell_line: str,
    show_sample_labels: bool,
) -> None:
    sub = df.dropna(subset=[xcol, ycol]).copy()
    for mode in mode_order:
        rows = sub[sub["transcriptome_mode"] == mode]
        for _, row in rows.iterrows():
            x = float(row[xcol]) * 100.0
            y = float(row[ycol]) * 100.0
            ax.scatter(
                x, y,
                s=52,
                c=styler.color(mode),
                marker=styler.marker(mode),
                edgecolors="white",
                linewidth=0.45,
                alpha=0.90,
                zorder=2,
            )
            if show_sample_labels:
                ax.annotate(
                    _sample_short(str(row["sample_id"]), cell_line),
                    (x, y),
                    textcoords="offset points",
                    xytext=(3, 2),
                    fontsize=4.8,
                    color="#555555",
                    alpha=0.9,
                    zorder=3,
                )
    ax.plot([0, 100], [0, 100], "--", color="#cccccc", linewidth=0.8, zorder=0)
    ax.set_xlim(0, 105)
    ax.set_ylim(0, 105)
    style_ax(ax, xlabel=xlabel, ylabel=ylabel, title=title)
    ax.grid(axis="y", color="#dddddd", linestyle=":", linewidth=0.55, alpha=0.8)
    ax.set_axisbelow(True)


def write_plots(summary: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_dir / "cell_line_precision_support.tsv", sep="\t", index=False)

    mode_order = _mode_order(summary["transcriptome_mode"].dropna().unique())
    styler = ModeStyler(mode_order)
    handles = [styler.legend_handle(m, label=_short_mode(m), markersize=6) for m in mode_order]

    for cell_line, cell_df in summary.groupby("cell_line", sort=True):
        if cell_df.empty:
            continue
        fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.25), sharex=False, sharey=False)
        specs = [
            (
                "5prime_ref_precision",
                "5prime_cage_precision",
                "5′ reference precision (%)",
                "5′ CAGE precision (%)",
                "5′ ends",
            ),
            (
                "3prime_ref_precision",
                "3prime_drna_precision",
                "3′ reference precision (%)",
                "3′ dRNA precision (%)",
                "3′ ends",
            ),
            (
                "sj_gtf_precision",
                "sj_star_precision",
                "GTF SJ precision (%)",
                "STAR SJ precision (%)",
                "splice junctions",
            ),
        ]
        for ax, (xcol, ycol, xlabel, ylabel, title) in zip(axes, specs):
            show_sample_labels = cell_df["sample_id"].nunique() > 1
            _plot_panel(
                ax,
                cell_df,
                xcol=xcol,
                ycol=ycol,
                xlabel=xlabel,
                ylabel=ylabel,
                title=title,
                mode_order=mode_order,
                styler=styler,
                cell_line=str(cell_line),
                show_sample_labels=show_sample_labels,
            )
        fig.suptitle(str(cell_line), fontsize=10, y=0.98)
        legend_outside(
            fig,
            handles=handles,
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            ncol=1,
            fontsize=6,
        )
        fig.tight_layout(rect=(0, 0, 0.86, 0.94))
        safe_name = str(cell_line).replace("/", "_").replace(" ", "_")
        savefig(fig, output_dir / f"{safe_name}_precision_support_three_panel.png", dpi=300)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test-name", default=None)
    parser.add_argument("--ortho-pr", nargs="+", required=True)
    parser.add_argument("--gtf-pr", nargs="+", required=True)
    parser.add_argument("--isoform", nargs="+", required=True,
                        help="sample_id:transcriptome_mode:path entries")
    parser.add_argument("--gtf", nargs="+", required=True,
                        help="sample_id:path entries")
    parser.add_argument("--star-sj", nargs="+", required=True,
                        help="cell_line:path STAR SJ.out.tab entries")
    parser.add_argument("--min-star-unique", type=int, default=1)
    parser.add_argument("--output", required=True)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    ortho_df = load_pr_tables(args.ortho_pr, kind="ortho", test_name=args.test_name)
    gtf_df = load_pr_tables(args.gtf_pr, kind="gtf", test_name=args.test_name)
    if ortho_df.empty or gtf_df.empty:
        print("ERROR: missing orthogonal or GTF precision tables", file=sys.stderr)
        sys.exit(1)

    star_labels = []
    for entry in args.star_sj:
        parts = _split_label_path(entry, 2)
        if parts:
            star_labels.append(parts[0])

    sj_df = compute_sj_precision(
        args.isoform,
        args.gtf,
        args.star_sj,
        min_unique_reads=args.min_star_unique,
    )
    summary = build_summary(ortho_df, gtf_df, sj_df, star_labels=star_labels)
    if summary.empty:
        print("ERROR: no rows available for plotting", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print(
            f"Loaded {len(summary)} sample x method rows across "
            f"{summary['cell_line'].nunique()} cell lines",
            file=sys.stderr,
        )

    write_plots(summary, Path(args.output))


if __name__ == "__main__":
    main()
