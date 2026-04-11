#!/usr/bin/env python3
"""
isoform_diversity_plot.py
─────────────────────────
Parallel-coordinates visualisation of isoform diversity across genomic
regions, computed directly from GTF files.

For each (GTF × region) combination the script computes:
    • gene_count           – number of distinct gene_id values
    • transcript_count     – number of distinct transcript_id values
    • mean_tx_per_gene     – mean transcripts per gene
    • median_tx_per_gene   – median transcripts per gene
    • mean_exons_per_tx    – mean exon count across transcripts
    • biotype composition  – top-5 transcript_biotype fractions + other
    • shannon_entropy      – transcript-per-gene entropy (diversity index)

Usage:
    # Compare reference GTF diversity across regions
    python isoform_diversity_plot.py \
        --gtfs ref.gtf \
        --labels "Reference" \
        --regions chr3:48000000-53000000 chr11:64000000-69000000 \
        --output diversity.png

    # Compare reference vs FLAIR output
    python isoform_diversity_plot.py \
        --gtfs ref.gtf flair.gtf \
        --labels Reference FLAIR_default \
        --regions chr3:48000000-53000000 chr11:64000000-69000000 \
        --output diversity.png

Requires: plotly, scipy, python-kaleido (for PNG export).
"""

import argparse
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import entropy

from .pub_style import PLOTLY_LAYOUT, PALETTE as _OI_PUB_COLORS

# ---------------------------------------------------------------------------
# GTF parsing
# ---------------------------------------------------------------------------

_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def parse_gtf_regions(gtf_path, regions):
    """
    Parse a GTF and return per-region gene/transcript structures.

    Parameters
    ----------
    gtf_path : str
        Path to GTF file.
    regions : list[tuple[str, int, int]]
        Parsed regions as (chrom, start, end).

    Returns
    -------
    dict[str, dict]
        Keyed by region string, value is dict with:
            'genes': {gene_id: {'biotype': str, 'transcripts': {tx_id: {'biotype': str, 'exon_count': int}}}}
    """
    region_data = {}
    for chrom, start, end in regions:
        region_data[f"{chrom}:{start}-{end}"] = {"genes": defaultdict(lambda: {"biotype": "", "transcripts": defaultdict(lambda: {"biotype": "", "exon_count": 0, "genomic_span": 0, "spliced_length": 0})})}

    # Build chrom -> [(start, end, region_key)] index for fast look-up
    chrom_index = defaultdict(list)
    for chrom, start, end in regions:
        key = f"{chrom}:{start}-{end}"
        chrom_index[chrom].append((start, end, key))

    with open(gtf_path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            # Fast chromosome check before full parsing
            tab1 = line.index('\t')
            chrom = line[:tab1]
            if chrom not in chrom_index:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature not in ("gene", "transcript", "exon"):
                continue
            fstart = int(fields[3])
            fend = int(fields[4])
            attrs = dict(_ATTR_RE.findall(fields[8]))

            # Check which region(s) this feature falls in
            for rstart, rend, rkey in chrom_index[chrom]:
                if fstart >= rstart and fend <= rend:
                    gene_id = attrs.get("gene_id", "")
                    tx_id = attrs.get("transcript_id", "")
                    if not gene_id:
                        continue

                    rd = region_data[rkey]["genes"]
                    if feature == "gene":
                        rd[gene_id]["biotype"] = attrs.get("gene_type", attrs.get("gene_biotype", ""))
                    elif feature == "transcript":
                        if tx_id:
                            rd[gene_id]["transcripts"][tx_id]["biotype"] = attrs.get(
                                "transcript_type", attrs.get("transcript_biotype", "")
                            )
                            rd[gene_id]["transcripts"][tx_id]["genomic_span"] = fend - fstart
                    elif feature == "exon":
                        if tx_id:
                            rd[gene_id]["transcripts"][tx_id]["exon_count"] += 1
                            rd[gene_id]["transcripts"][tx_id]["spliced_length"] += (fend - fstart)

    return region_data


def parse_region_string(region_str):
    """Parse 'chr1:100-200' into ('chr1', 100, 200)."""
    m = re.match(r"^(\w+):(\d+)-(\d+)$", region_str)
    if not m:
        raise ValueError(f"Cannot parse region: {region_str}")
    return m.group(1), int(m.group(2)), int(m.group(3))


# ---------------------------------------------------------------------------
# Metrics computation
# ---------------------------------------------------------------------------


def compute_metrics(region_data, top_biotypes=None):
    """
    Compute diversity metrics for one GTF's region data.

    Returns
    -------
    list[dict]
        One dict per region with computed metrics.
    top_biotypes : list[str] or None
        If None, determined from data.
    """
    # First pass: determine top biotypes across all regions if not given
    if top_biotypes is None:
        all_bt = Counter()
        for rkey, rinfo in region_data.items():
            for gene_id, ginfo in rinfo["genes"].items():
                for tx_id, tinfo in ginfo["transcripts"].items():
                    bt = tinfo["biotype"].strip().lower() if tinfo["biotype"] else "unknown"
                    all_bt[bt] += 1
        top_biotypes = [bt for bt, _ in all_bt.most_common(5)]

    rows = []
    for rkey in sorted(region_data.keys()):
        genes = region_data[rkey]["genes"]
        gene_count = len(genes)
        if gene_count == 0:
            rows.append(
                {
                    "region": rkey,
                    "gene_count": 0,
                    "transcript_count": 0,
                    "mean_tx_per_gene": 0,
                    "median_tx_per_gene": 0,
                    "mean_exons_per_tx": 0,
                    "shannon_entropy": 0,
                    **{f"pct_{bt}": 0 for bt in top_biotypes},
                    "pct_other": 0,
                }
            )
            continue

        tx_per_gene = []
        exon_counts = []
        bt_counter = Counter()
        spliced_lengths = []
        for gene_id, ginfo in genes.items():
            n_tx = len(ginfo["transcripts"])
            tx_per_gene.append(n_tx)
            for tx_id, tinfo in ginfo["transcripts"].items():
                exon_counts.append(tinfo["exon_count"])
                bt = tinfo["biotype"].strip().lower() if tinfo["biotype"] else "unknown"
                bt_counter[bt] += 1
                sl = tinfo.get("spliced_length", 0)
                if sl > 0:
                    spliced_lengths.append(sl)

        tx_per_gene = np.array(tx_per_gene, dtype=float)
        total_tx = int(tx_per_gene.sum())

        # Shannon entropy of transcript distribution across genes
        if tx_per_gene.sum() > 0:
            probs = tx_per_gene / tx_per_gene.sum()
            h = float(entropy(probs, base=2))
        else:
            h = 0.0

        # Biotype fractions
        total_bt = sum(bt_counter.values()) or 1
        pct = {}
        for bt in top_biotypes:
            pct[f"pct_{bt}"] = 100.0 * bt_counter.get(bt, 0) / total_bt
        covered = sum(bt_counter.get(bt, 0) for bt in top_biotypes)
        pct["pct_other"] = 100.0 * (total_bt - covered) / total_bt

        # Truncation-relevant metrics
        if spliced_lengths:
            median_spliced_len = float(np.median(spliced_lengths))
            pct_long = 100.0 * sum(1 for sl in spliced_lengths if sl > 2500) / len(spliced_lengths)
            pct_very_long = 100.0 * sum(1 for sl in spliced_lengths if sl > 5000) / len(spliced_lengths)
        else:
            median_spliced_len = 0.0
            pct_long = 0.0
            pct_very_long = 0.0

        rows.append(
            {
                "region": rkey,
                "gene_count": gene_count,
                "transcript_count": total_tx,
                "mean_tx_per_gene": float(np.mean(tx_per_gene)),
                "median_tx_per_gene": float(np.median(tx_per_gene)),
                "mean_exons_per_tx": float(np.mean(exon_counts)) if exon_counts else 0.0,
                "shannon_entropy": h,
                "median_spliced_len": median_spliced_len,
                "pct_long_tx": pct_long,
                "pct_very_long_tx": pct_very_long,
                **pct,
            }
        )

    return rows, top_biotypes


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

# Okabe-Ito palette (colorblind-safe)
_OI_COLORS = [
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # pink
    "#000000",  # black
]


def plot_parallel_coordinates(all_rows, labels, top_biotypes, output, title_prefix=""):
    """
    Create parallel coordinates plot.  One trace per GTF label, one line
    per region.
    """
    metrics = [
        "gene_count",
        "transcript_count",
        "mean_tx_per_gene",
        "median_tx_per_gene",
        "mean_exons_per_tx",
        "shannon_entropy",
        "median_spliced_len",
        "pct_long_tx",
        "pct_very_long_tx",
    ] + [f"pct_{bt}" for bt in top_biotypes] + ["pct_other"]

    pretty = {
        "gene_count": "Genes",
        "transcript_count": "Transcripts",
        "mean_tx_per_gene": "Mean Tx/Gene",
        "median_tx_per_gene": "Med Tx/Gene",
        "mean_exons_per_tx": "Mean Exons/Tx",
        "shannon_entropy": "Shannon H",
        "median_spliced_len": "Med Spliced Len",
        "pct_long_tx": "% Tx >2.5kb",
        "pct_very_long_tx": "% Tx >5kb",
        "pct_other": "% other",
    }
    for bt in top_biotypes:
        pretty[f"pct_{bt}"] = f"% {bt[:20]}"

    n_labels = len(labels)

    if n_labels == 1:
        # Single GTF — colour lines by transcript count
        rows = all_rows[0]
        dims = _build_dims(rows, metrics, pretty)
        tx_vals = [r["transcript_count"] for r in rows]
        trace = go.Parcoords(
            line=dict(
                color=tx_vals,
                colorscale="Viridis",
                showscale=True,
                colorbar=dict(title="Transcripts"),
            ),
            dimensions=dims,
        )
        fig = go.Figure(data=[trace])
        title = f"{title_prefix}Isoform Diversity — {labels[0]}"
    else:
        # Multiple GTFs — one trace per label, each a distinct colour
        fig = go.Figure()
        for idx, (rows, label) in enumerate(zip(all_rows, labels)):
            color = _OI_COLORS[idx % len(_OI_COLORS)]
            dims = _build_dims(rows, metrics, pretty)
            # Encode constant colour as numeric to satisfy Parcoords API
            fig.add_trace(
                go.Parcoords(
                    line=dict(color=[idx] * len(rows), colorscale=[[0, color], [1, color]]),
                    dimensions=dims,
                    name=label,
                )
            )
        title = f"{title_prefix}Isoform Diversity Comparison"

    fig.update_layout(
        title=title,
        **PLOTLY_LAYOUT,
    )

    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    # Save HTML (interactive) alongside PNG
    html_path = output.rsplit(".", 1)[0] + ".html"
    fig.write_html(html_path)
    try:
        fig.write_image(output, width=1400, height=550, scale=2)
        print(f"Saved PNG: {output}")
    except Exception as e:
        print(f"Warning: PNG export failed ({e}). HTML saved at {html_path}", file=sys.stderr)
    print(f"Saved HTML: {html_path}")


def _build_dims(rows, metrics, pretty, num_ticks=5):
    """Build plotly dimension dicts for a set of rows."""
    dims = []
    for m in metrics:
        vals = [r[m] for r in rows]
        mn, mx = min(vals), max(vals)
        if mn == mx:
            mn -= 1
            mx += 1
        tick_vals = np.linspace(mn, mx, num_ticks)
        if m.startswith("pct_") or m == "shannon_entropy":
            tick_txt = [f"{v:.1f}" for v in tick_vals]
        else:
            tick_vals = np.round(tick_vals).astype(int)
            tick_txt = [str(v) for v in tick_vals]

        dims.append(
            dict(
                range=[mn, mx],
                tickvals=tick_vals.tolist(),
                ticktext=tick_txt,
                label=pretty.get(m, m),
                values=vals,
            )
        )
    # Add region as categorical axis (first dimension)
    region_labels = [r["region"] for r in rows]
    unique_regions = sorted(set(region_labels))
    region_to_idx = {r: i for i, r in enumerate(unique_regions)}
    dims.insert(
        0,
        dict(
            range=[0, len(unique_regions) - 1],
            tickvals=list(range(len(unique_regions))),
            ticktext=unique_regions,
            label="Region",
            values=[region_to_idx[r] for r in region_labels],
        ),
    )
    return dims


# ---------------------------------------------------------------------------
# Bar chart fallback (static matplotlib) for 1-GTF, few-region use
# ---------------------------------------------------------------------------


def plot_diversity_bars(rows, output, title_prefix=""):
    """
    Simple matplotlib grouped-bar chart when plotly is overkill / for
    quick publication-style snapshot.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping bar chart", file=sys.stderr)
        return

    regions = [r["region"] for r in rows]
    fig, axes = plt.subplots(3, 3, figsize=(7.2, 5.5))
    bar_metrics = [
        ("gene_count", "Genes"),
        ("transcript_count", "Transcripts"),
        ("mean_tx_per_gene", "Mean Tx/Gene"),
        ("median_tx_per_gene", "Med Tx/Gene"),
        ("mean_exons_per_tx", "Mean Exons/Tx"),
        ("shannon_entropy", "Shannon Entropy"),
        ("median_spliced_len", "Med Spliced Len (bp)"),
        ("pct_long_tx", "% Tx >2.5kb"),
        ("pct_very_long_tx", "% Tx >5kb"),
    ]
    for ax, (metric, ylabel) in zip(axes.flat, bar_metrics):
        vals = [r[metric] for r in rows]
        x = np.arange(len(regions))
        ax.bar(x, vals, color=_OI_COLORS[: len(regions)])
        ax.set_xticks(x)
        ax.set_xticklabels([r.split(":")[0] for r in regions], rotation=45, ha="right", fontsize=8)
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel)

    fig.suptitle(f"{title_prefix}Per-Region Isoform Diversity", fontsize=8)
    fig.tight_layout()
    bar_path = output.rsplit(".", 1)[0] + "_bars.png"
    fig.savefig(bar_path, dpi=300)
    plt.close(fig)
    print(f"Saved bar chart: {bar_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    p = argparse.ArgumentParser(
        description="Compute and plot isoform diversity across genomic regions from GTF files."
    )
    p.add_argument(
        "--gtfs",
        nargs="+",
        required=True,
        help="One or more GTF files to analyse.",
    )
    p.add_argument(
        "--labels",
        nargs="+",
        help="Labels for each GTF (default: filename stems).",
    )
    p.add_argument(
        "--regions",
        nargs="+",
        required=True,
        help="Regions in chr:start-end format.",
    )
    p.add_argument(
        "--output",
        required=True,
        help="Output path (PNG). An HTML will also be saved alongside.",
    )
    p.add_argument("--title-prefix", default="", help="Prefix for plot title.")
    p.add_argument("--bars", action="store_true", help="Also produce a bar chart (matplotlib).")
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args()

    parsed_regions = [parse_region_string(r) for r in args.regions]
    labels = args.labels or [os.path.splitext(os.path.basename(g))[0] for g in args.gtfs]
    if len(labels) != len(args.gtfs):
        p.error("Number of --labels must match number of --gtfs")

    all_rows = []
    top_biotypes = None  # share across GTFs for consistent axes
    for gtf_path, label in zip(args.gtfs, labels):
        if args.verbose:
            print(f"Parsing {gtf_path} as '{label}' …")
        region_data = parse_gtf_regions(gtf_path, parsed_regions)
        rows, top_biotypes = compute_metrics(region_data, top_biotypes)
        if args.verbose:
            for row in rows:
                print(f"  {row['region']}: {row['gene_count']} genes, {row['transcript_count']} tx")
        all_rows.append(rows)

    plot_parallel_coordinates(all_rows, labels, top_biotypes, args.output, title_prefix=args.title_prefix)

    if args.bars and len(all_rows) == 1:
        plot_diversity_bars(all_rows[0], args.output, title_prefix=args.title_prefix)

    print("Done.")


if __name__ == "__main__":
    main()
