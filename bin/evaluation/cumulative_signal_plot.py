#!/usr/bin/env python3
"""
Cumulative orthogonal signal plot — isoforms ranked by read support.

For each method, isoforms are sorted by read-support score (descending),
and the cumulative fraction of total TSS + TTS signal is plotted against
the fraction of isoforms.  Requires BED12 isoform files plus four bedGraph
signal tracks.

Read counts come from isoform read-map files (--read-map label:path pairs).
Methods without a matching --read-map entry are skipped.

Usage:
    python cumulative_signal_plot.py \\
        --bed label1:bed1.bed label2:bed2.bed ... \\
        --read-map label1:map1.txt label2:map2.txt ... \\
        --cage-plus cage_plus.bg --cage-minus cage_minus.bg \\
        --qs-plus qs_plus.bg --qs-minus qs_minus.bg \\
        --output output_dir/
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pub_style import ModeStyler, style_ax, legend_outside, savefig, W1, apply_rc
from signal_utils import parse_isoforms, load_signal_tracks, isoform_signal


# ── Read-map helpers ────────────────────────────────────────────────────────

def load_read_map(path: str | Path) -> dict[str, int]:
    """Load an isoform read-map file → {isoform_id: read_count}.

    Self-referencing maps (e.g. StringTie2 where each 'read' is the
    isoform name itself) are detected and return an empty dict.
    """
    counts: dict[str, int] = {}
    n_self = 0
    n_total = 0
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            iso_id = parts[0]
            reads = parts[1].split(",")
            n_reads = len(reads)
            counts[iso_id] = n_reads
            n_total += 1
            if n_reads == 1 and reads[0] == iso_id:
                n_self += 1
    # Skip self-referencing read maps (e.g., StringTie2)
    if n_total > 0 and n_self / n_total > 0.9:
        return {}
    return counts


_ENSG_RE = re.compile(r"_ENSG\d")


def _lookup_read_count(name: str, read_counts: dict[str, int]) -> int | None:
    """Look up read count for an isoform name, handling tid_gid naming."""
    if name in read_counts:
        return read_counts[name]
    # GTF-parsed names are "tid_gid" — try stripping the gene_id suffix.
    # First try splitting at _ENSG (covers GENCODE gene IDs).
    m = _ENSG_RE.search(name)
    if m:
        tid = name[: m.start()]
        if tid in read_counts:
            return read_counts[tid]
    # Fallback: split on last underscore (covers BambuGene, MSTRG, etc.)
    idx = name.rfind("_")
    if idx > 0:
        tid = name[:idx]
        if tid in read_counts:
            return read_counts[tid]
    return None


# ── Data preparation ────────────────────────────────────────────────────────

def _build_iso_data(
    beds_by_method: dict,
    read_maps_by_method: dict,
    cage_p, cage_m, qs_p, qs_m,
) -> dict[str, list[tuple[int, float, float, float]]]:
    """Build per-isoform (read_count, tss_signal, tts_signal, total_signal)
    for each method that has a read-map.

    Returns dict[method_label, list[(count, tss, tts, total)]].
    """
    data: dict[str, list[tuple[int, float, float, float]]] = {}
    for m in beds_by_method:
        if m not in read_maps_by_method:
            continue
        rc = read_maps_by_method[m]
        rows = []
        for iso in beds_by_method[m]:
            count = _lookup_read_count(iso["name"], rc)
            if count is None:
                count = 0
            tss_sig, tts_sig = isoform_signal(iso, cage_p, cage_m, qs_p, qs_m)
            rows.append((count, tss_sig, tts_sig, tss_sig + tts_sig))
        data[m] = rows
    return data


# ── Plotting ────────────────────────────────────────────────────────────────

def plot_cumulative_signal(
    iso_data_by_method: dict[str, list],
    output_dir: Path,
):
    """Cumulative orthogonal signal curves per method."""
    methods = list(iso_data_by_method.keys())
    if not methods:
        return

    apply_rc()
    styler = ModeStyler(methods)

    # Wide figure: plot on left, legend on right
    fig, ax = plt.subplots(figsize=(W1 * 1.8, W1 * 0.82))
    handles = []

    for m in methods:
        rows = sorted(iso_data_by_method[m], key=lambda x: -x[0])
        sigs = np.array([r[3] for r in rows])
        cum = np.cumsum(sigs)
        frac = np.arange(1, len(cum) + 1) / len(cum)
        cum_n = cum / cum[-1] if cum[-1] > 0 else cum

        line, = ax.plot(
            frac, cum_n,
            color=styler.color(m),
            linewidth=1.0,
            alpha=0.85,
        )
        line.set_dashes(styler.dash(m))
        # Marker at ~8 evenly-spaced positions so lines stay distinguishable
        step = max(1, len(frac) // 8)
        ax.plot(
            frac[::step], cum_n[::step],
            color=styler.color(m),
            marker=styler.marker(m),
            markersize=3.5,
            linestyle="",
            alpha=0.9,
        )
        handles.append(styler.legend_handle(m, label=m, markersize=5))

    style_ax(ax,
             xlabel="Fraction of isoforms (by read support)",
             ylabel="Cumulative signal (norm.)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_axisbelow(True)

    ncol = max(1, len(methods) // 20)
    legend_outside(fig, handles=handles, loc="upper left",
                   bbox_to_anchor=(1.02, 1.0), ncol=ncol, fontsize=6)
    fig.tight_layout(pad=0.3)
    savefig(fig, output_dir / "cumulative_signal.png")


def plot_gini_barplot(
    iso_data_by_method: dict[str, list],
    output_dir: Path,
):
    """Gini coefficient bar plot — higher = signal concentrated in top isoforms."""
    methods = list(iso_data_by_method.keys())
    if not methods:
        return

    gini_vals = {}
    for m in methods:
        rows = sorted(iso_data_by_method[m], key=lambda x: x[3])  # ascending signal
        sigs = np.array([r[3] for r in rows])
        n = len(sigs)
        if n < 2 or sigs.sum() == 0:
            continue
        cumw = np.cumsum(sigs)
        gini = 1 - 2 * cumw.sum() / (n * sigs.sum()) + 1 / n
        gini_vals[m] = gini

    if not gini_vals:
        return

    labels = list(gini_vals.keys())
    vals = [gini_vals[l] for l in labels]
    colors = [ModeStyler(labels).color(l) for l in labels]

    fig, ax = plt.subplots(figsize=(W1, W1 * 0.6))
    ax.barh(range(len(labels)), vals, color=colors, height=0.6, edgecolor="none")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=6)
    ax.invert_yaxis()
    style_ax(ax, xlabel="Gini coefficient (signal concentration)")
    ax.set_xlim(0, 1)
    fig.tight_layout(pad=0.3)
    savefig(fig, output_dir / "signal_gini.png")


def plot_zero_signal_fraction(
    iso_data_by_method: dict[str, list],
    output_dir: Path,
):
    """Bar plot of fraction of isoforms with zero signal at both TSS and TTS."""
    methods = list(iso_data_by_method.keys())
    if not methods:
        return

    fracs_tss = []
    fracs_tts = []
    fracs_both = []
    labels = []
    for m in methods:
        rows = iso_data_by_method[m]
        n = len(rows)
        if n == 0:
            continue
        n_zero_tss = sum(1 for r in rows if r[1] == 0)
        n_zero_tts = sum(1 for r in rows if r[2] == 0)
        n_zero_both = sum(1 for r in rows if r[1] == 0 and r[2] == 0)
        fracs_tss.append(n_zero_tss / n)
        fracs_tts.append(n_zero_tts / n)
        fracs_both.append(n_zero_both / n)
        labels.append(m)

    if not labels:
        return

    x = np.arange(len(labels))
    w = 0.25
    fig, ax = plt.subplots(figsize=(W1, W1 * 0.65))
    ax.bar(x - w, fracs_tss, w, label="Zero TSS (CAGE)", color="#4C72B0", edgecolor="none")
    ax.bar(x, fracs_tts, w, label="Zero TTS (dRNA)", color="#DD8452", edgecolor="none")
    ax.bar(x + w, fracs_both, w, label="Zero both", color="#C44E52", edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=5, rotation=45, ha="right")
    ax.set_ylim(0, 1)
    style_ax(ax, ylabel="Fraction of isoforms")
    ax.legend(frameon=False, fontsize=5, loc="upper right")
    fig.tight_layout(pad=0.3)
    savefig(fig, output_dir / "zero_signal_fraction.png")


def plot_zero_signal_counts(
    iso_data_by_method: dict[str, list],
    output_dir: Path,
):
    """Bar plot of raw isoform counts with zero signal at TSS, TTS, or both."""
    methods = list(iso_data_by_method.keys())
    if not methods:
        return

    counts_tss = []
    counts_tts = []
    counts_both = []
    labels = []
    for m in methods:
        rows = iso_data_by_method[m]
        if not rows:
            continue
        counts_tss.append(sum(1 for r in rows if r[1] == 0))
        counts_tts.append(sum(1 for r in rows if r[2] == 0))
        counts_both.append(sum(1 for r in rows if r[1] == 0 and r[2] == 0))
        labels.append(m)

    if not labels:
        return

    x = np.arange(len(labels))
    w = 0.25
    fig, ax = plt.subplots(figsize=(W1, W1 * 0.65))
    ax.bar(x - w, counts_tss, w, label="Zero TSS (CAGE)", color="#4C72B0", edgecolor="none")
    ax.bar(x, counts_tts, w, label="Zero TTS (dRNA)", color="#DD8452", edgecolor="none")
    ax.bar(x + w, counts_both, w, label="Zero both", color="#C44E52", edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=5, rotation=45, ha="right")
    style_ax(ax, ylabel="Number of isoforms")
    ax.legend(frameon=False, fontsize=5, loc="upper right")
    fig.tight_layout(pad=0.3)
    savefig(fig, output_dir / "zero_signal_counts.png")


def _signal_efficiency_scatter(
    iso_data_by_method: dict[str, list],
    signal_extractor,
    ylabel: str,
    output_path: Path,
):
    """Scatter: X = total isoforms, Y = normalised captured signal (one point per method).

    signal_extractor(rows) → float  — compute the raw signal total for one method.
    Signal is normalised by the global max across all methods before plotting.
    Upper-left = efficient (few isoforms, high signal capture).
    """
    methods = list(iso_data_by_method.keys())
    if not methods:
        return

    raw = {m: signal_extractor(iso_data_by_method[m]) for m in methods}
    global_max = max(raw.values()) or 1.0

    styler = ModeStyler(methods)
    xs = [len(iso_data_by_method[m]) for m in methods]
    ys = [raw[m] / global_max for m in methods]
    colors = [styler.color(m) for m in methods]
    markers = [styler.marker(m) for m in methods]

    fig, ax = plt.subplots(figsize=(W1, W1 * 0.82))
    for i, m in enumerate(methods):
        ax.scatter(xs[i], ys[i], c=colors[i], marker=markers[i],
                   s=35, zorder=3, edgecolors="white", linewidths=0.3)
        ax.annotate(m, (xs[i], ys[i]), fontsize=5, ha="left", va="bottom",
                    xytext=(3, 3), textcoords="offset points")
    style_ax(ax, xlabel="Total isoforms", ylabel=ylabel)
    ax.set_axisbelow(True)
    fig.tight_layout(pad=0.3)
    savefig(fig, output_path)


def plot_signal_efficiency(
    iso_data_by_method: dict[str, list],
    output_dir: Path,
):
    """Three scatter plots: combined, 5'-only (CAGE), 3'-only (dRNA)."""
    # Combined: TSS + TTS each normalised independently, then summed → max 2.0.
    # Pre-compute the per-end maxima so the extractor closure sees fixed values.
    methods = list(iso_data_by_method.keys())
    tss_max = max(sum(r[1] for r in iso_data_by_method[m]) for m in methods) or 1.0
    tts_max = max(sum(r[2] for r in iso_data_by_method[m]) for m in methods) or 1.0

    def _combined(rows):
        return sum(r[1] for r in rows) / tss_max + sum(r[2] for r in rows) / tts_max

    _signal_efficiency_scatter(
        iso_data_by_method,
        signal_extractor=_combined,
        ylabel="Normalised captured signal\n(CAGE + dRNA, each 0–1)",
        output_path=output_dir / "signal_efficiency.png",
    )
    _signal_efficiency_scatter(
        iso_data_by_method,
        signal_extractor=lambda rows: sum(r[1] for r in rows),
        ylabel="Total 5′ signal captured (CAGE, raw)",
        output_path=output_dir / "signal_efficiency_5prime.png",
    )
    _signal_efficiency_scatter(
        iso_data_by_method,
        signal_extractor=lambda rows: sum(r[2] for r in rows),
        ylabel="Total 3′ signal captured (dRNA, raw)",
        output_path=output_dir / "signal_efficiency_3prime.png",
    )


def plot_read_count_vs_signal(
    iso_data_by_method: dict[str, list],
    output_dir: Path,
):
    """Per-method scatter: isoform read count vs boundary signal.

    Reports Spearman rho in each panel title.
    """
    from scipy.stats import spearmanr

    methods = list(iso_data_by_method.keys())
    if not methods:
        return

    styler = ModeStyler(methods)
    n = len(methods)
    ncols = min(n, 3)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(W1 * ncols / 3, W1 * 0.82 * nrows / 3),
                             squeeze=False)

    for idx, m in enumerate(methods):
        ax = axes[idx // ncols][idx % ncols]
        rows = iso_data_by_method[m]
        counts = np.array([r[0] for r in rows], dtype=float)
        sigs = np.array([r[3] for r in rows])
        mask = counts > 0
        counts_f, sigs_f = counts[mask], sigs[mask]
        if len(counts_f) > 1:
            rho, _ = spearmanr(counts_f, sigs_f)
        else:
            rho = float("nan")
        eps = 1e-3
        ax.scatter(counts_f + eps, sigs_f + eps, s=1.5, alpha=0.3,
                   color=styler.color(m), edgecolors="none", rasterized=True)
        ax.set_xscale("log")
        ax.set_yscale("log")
        style_ax(ax)
        ax.text(0.04, 0.96, f"{m}\n\u03C1 = {rho:.2f}", transform=ax.transAxes,
                ha="left", va="top", fontsize=5, fontweight="bold")
        if idx >= n - ncols:
            ax.set_xlabel("Read count", fontsize=6)
        if idx % ncols == 0:
            ax.set_ylabel("Boundary signal", fontsize=6)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.tight_layout(pad=0.3)
    savefig(fig, output_dir / "read_count_vs_signal.png")


# ── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--bed", nargs="+", required=True,
        help="label:path pairs for BED12 isoform files",
    )
    parser.add_argument(
        "--read-map", nargs="+", default=[],
        help="label:path pairs for isoform read-map files",
    )
    parser.add_argument("--cage-plus",  required=True, help="CAGE bedGraph (+ strand)")
    parser.add_argument("--cage-minus", required=True, help="CAGE bedGraph (- strand)")
    parser.add_argument("--qs-plus",    required=True, help="dRNA bedGraph (+ strand)")
    parser.add_argument("--qs-minus",   required=True, help="dRNA bedGraph (- strand)")
    parser.add_argument("--output",     required=True, help="Output directory")
    parser.add_argument("--verbose",    action="store_true")
    args = parser.parse_args()

    beds_by_method = {}
    for entry in args.bed:
        if ":" not in entry:
            print(f"WARNING: skipping malformed entry '{entry}'", file=sys.stderr)
            continue
        label, path = entry.split(":", 1)
        if not Path(path).exists():
            print(f"WARNING: file not found: {path}", file=sys.stderr)
            continue
        isoforms = parse_isoforms(path)
        if isoforms:
            beds_by_method[label] = isoforms
            if args.verbose:
                print(f"  {label}: {len(isoforms)} isoforms")

    if not beds_by_method:
        print("No BED data loaded — skipping", file=sys.stderr)
        sys.exit(1)

    # Load read-map files
    read_maps_by_method: dict[str, dict[str, int]] = {}
    for entry in args.read_map:
        if ":" not in entry:
            print(f"WARNING: skipping malformed read-map entry '{entry}'", file=sys.stderr)
            continue
        label, path = entry.split(":", 1)
        if not Path(path).exists():
            print(f"WARNING: read-map not found: {path}", file=sys.stderr)
            continue
        if label not in beds_by_method:
            if args.verbose:
                print(f"  Skipping read-map for '{label}' (no matching BED)", file=sys.stderr)
            continue
        rc = load_read_map(path)
        if rc:
            read_maps_by_method[label] = rc
            if args.verbose:
                print(f"  {label}: read-map has {len(rc)} isoforms")

    if not read_maps_by_method:
        print("No read-map data loaded — skipping cumulative signal plot", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print("  Loading signal tracks...")
    cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
        args.cage_plus, args.cage_minus, args.qs_plus, args.qs_minus,
    )

    if args.verbose:
        print("  Computing per-isoform signal...")
    iso_data = _build_iso_data(
        beds_by_method, read_maps_by_method,
        cage_p, cage_m, qs_p, qs_m,
    )
    if not iso_data:
        print("No methods with matched read-map data — skipping", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_cumulative_signal(iso_data, output_dir)
    plot_gini_barplot(iso_data, output_dir)
    plot_zero_signal_fraction(iso_data, output_dir)
    plot_zero_signal_counts(iso_data, output_dir)
    plot_signal_efficiency(iso_data, output_dir)
    plot_read_count_vs_signal(iso_data, output_dir)
    print(f"Saved cumulative signal plots to {args.output}")


if __name__ == "__main__":
    main()
