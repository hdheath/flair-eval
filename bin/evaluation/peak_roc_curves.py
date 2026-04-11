#!/usr/bin/env python3
"""
peak_roc_curves.py — Signal-stratified peak recovery curves per tool.

Each peak_reason file has columns:
    peak_id  score  read_count  reason

The `score` is the orthogonal signal intensity of the REFERENCE peak —
identical for every tool on any given peak.  Traditional ROC (sweep tool
confidence → TP/FP rate) is not meaningful here.

What IS meaningful: signal-stratified recall.
  "Among all reference peaks with signal >= t, what fraction does each
   tool recover?"

High-signal peaks should be easy to recover; curves diverge at low signal
where tools make different tradeoffs.

Outputs
-------
  <outdir>/signal_recall_cage.png
  <outdir>/signal_recall_drna.png
  <outdir>/signal_recall_cage_cross_tech.png   (mean ± SD across datasets)
  <outdir>/signal_recall_drna_cross_tech.png
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pub_style import ModeStyler, style_ax, legend_outside, savefig, apply_rc

# ── helpers ──────────────────────────────────────────────────────────────────

_RECOVERED = "recovered"
_TOOL_RE = re.compile(
    r'(?:^|/)([^/]+?)_transcriptome_(?:cage|drna)_peak_reasons\.tsv$'
)
_SIG_RE = re.compile(r'_transcriptome_(cage|drna)_peak_reasons\.tsv$')


def _tool_name(path: Path) -> str:
    m = _TOOL_RE.search(str(path))
    if m:
        full = m.group(1)
        parts = full.rsplit('_', 1)
        return parts[-1] if len(parts) > 1 else full
    return path.stem


def _sig_type(path: Path) -> str:
    m = _SIG_RE.search(str(path))
    return m.group(1) if m else 'unknown'


def _load(path: Path):
    """Return list of (score, is_recovered) for non-header rows."""
    rows = []
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            try:
                score = float(parts[1])
            except ValueError:
                score = 0.0
            rows.append((score, parts[3].strip() == _RECOVERED))
    return rows


def _signal_recall_curve(rows):
    """
    Sweep signal threshold high → low.
    At each threshold t, consider only peaks with score >= t.
    recall(t) = recovered_at_or_above_t / total_at_or_above_t
    """
    if not rows:
        return np.array([0.0]), np.array([0.0])
    scores = np.array([r[0] for r in rows])
    is_rec = np.array([r[1] for r in rows], dtype=bool)
    thresholds = np.unique(scores)[::-1]
    recalls = []
    for t in thresholds:
        mask = scores >= t
        total = mask.sum()
        recalls.append((mask & is_rec).sum() / total if total > 0 else 0.0)
    return thresholds, np.array(recalls)


def _auc_threshold(thresholds, recalls):
    """AUC over normalised threshold axis (0→1 after min-max scaling)."""
    if len(thresholds) < 2:
        return float(recalls[0]) if len(recalls) else 0.0
    t_min, t_max = thresholds.min(), thresholds.max()
    if t_max == t_min:
        return float(recalls.mean())
    t_norm = (thresholds - t_min) / (t_max - t_min)
    order = np.argsort(t_norm)
    return float(np.trapz(recalls[order], t_norm[order]))


# ── plotting ─────────────────────────────────────────────────────────────────

def _plot_signal_recall(curves_by_tool, output_path, title):
    tools = sorted(curves_by_tool)
    styler = ModeStyler(tools)
    apply_rc()
    fig, ax = plt.subplots(figsize=(3.8, 3.5))

    for tool in tools:
        thresholds, recalls = curves_by_tool[tool]
        line, = ax.plot(thresholds, recalls * 100,
                        color=styler.color(tool), linewidth=1.2, alpha=0.85)
        line.set_dashes(styler.dash(tool))
        step = max(1, len(thresholds) // 8)
        ax.plot(thresholds[::step], (recalls * 100)[::step],
                color=styler.color(tool), marker=styler.marker(tool),
                markersize=4, linestyle='', alpha=0.9)

    ax.set_xscale('log')
    ax.invert_xaxis()
    style_ax(ax,
             xlabel="Signal score threshold (log, high→low)",
             ylabel="Recall at threshold (%)",
             title=title)
    ax.set_ylim(-2, 102)

    handles = [
        styler.legend_handle(t,
            label=f"{t}  (AUC={_auc_threshold(*curves_by_tool[t]):.3f})",
            markersize=6)
        for t in tools
    ]
    legend_outside(fig, handles=handles, loc='upper left',
                   bbox_to_anchor=(1.02, 1.0), ncol=1, fontsize=6.5)
    fig.tight_layout()
    savefig(fig, output_path, dpi=300)
    plt.close(fig)


def _interpolate_recall(thresholds, recalls, grid):
    order = np.argsort(thresholds)
    return np.interp(grid, thresholds[order], recalls[order],
                     left=recalls[order][0], right=recalls[order][-1])


def _plot_cross_tech_signal_recall(per_tool_per_dataset, output_path, title,
                                    n_grid=300):
    all_thresholds = []
    for curves in per_tool_per_dataset.values():
        for t, _ in curves:
            all_thresholds.extend(t.tolist())
    if not all_thresholds:
        return
    t_min = max(1.0, min(all_thresholds))
    t_max = max(all_thresholds)
    grid = np.logspace(np.log10(t_min), np.log10(t_max), n_grid)

    tools = sorted(per_tool_per_dataset)
    styler = ModeStyler(tools)
    apply_rc()
    fig, ax = plt.subplots(figsize=(3.8, 3.5))

    for tool in tools:
        curves = per_tool_per_dataset[tool]
        if not curves:
            continue
        interped = np.array([_interpolate_recall(t, r, grid) for t, r in curves])
        mean_r = interped.mean(axis=0)
        std_r  = interped.std(axis=0)

        line, = ax.plot(grid, mean_r * 100,
                        color=styler.color(tool), linewidth=1.2, alpha=0.9)
        line.set_dashes(styler.dash(tool))
        if len(curves) > 1:
            ax.fill_between(grid,
                            np.clip((mean_r - std_r) * 100, 0, 100),
                            np.clip((mean_r + std_r) * 100, 0, 100),
                            color=styler.color(tool), alpha=0.15)

    ax.set_xscale('log')
    ax.invert_xaxis()
    style_ax(ax,
             xlabel="Signal score threshold (log, high→low)",
             ylabel="Recall at threshold (%)",
             title=title)
    ax.set_ylim(-2, 102)

    handles = []
    for tool in tools:
        curves = per_tool_per_dataset[tool]
        if not curves:
            continue
        interped = np.array([_interpolate_recall(t, r, grid) for t, r in curves])
        auc = _auc_threshold(grid, interped.mean(axis=0))
        handles.append(styler.legend_handle(tool,
            label=f"{tool}  (AUC={auc:.3f})", markersize=6))
    legend_outside(fig, handles=handles, loc='upper left',
                   bbox_to_anchor=(1.02, 1.0), ncol=1, fontsize=6.5)
    fig.tight_layout()
    savefig(fig, output_path, dpi=300)
    plt.close(fig)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Signal-stratified peak recovery curves per tool"
    )
    ap.add_argument('--input', '-i', nargs='+', required=True,
                    help='peak_reason TSV files (cage and/or drna)')
    ap.add_argument('--output', '-o', required=True,
                    help='Output directory')
    ap.add_argument('--dataset', default=None,
                    help='Dataset label (used in plot titles)')
    ap.add_argument('--cross-tech-only', action='store_true',
                    help='Only write cross-tech aggregate curves')
    ap.add_argument('--verbose', '-v', action='store_true')
    args = ap.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    by_sig = defaultdict(dict)
    cross_tech = defaultdict(lambda: defaultdict(list))

    for f in args.input:
        p = Path(f)
        if not p.exists():
            print(f"Warning: {f} not found, skipping", file=sys.stderr)
            continue
        sig  = _sig_type(p)
        tool = _tool_name(p)
        rows = _load(p)
        if not rows:
            continue
        thresholds, recalls = _signal_recall_curve(rows)
        by_sig[sig][tool] = (thresholds, recalls)
        cross_tech[sig][tool].append((thresholds, recalls))
        if args.verbose:
            auc = _auc_threshold(thresholds, recalls)
            n_rec = sum(1 for _, r in rows if r)
            print(f"  {tool} ({sig}): {len(rows)} peaks, "
                  f"{n_rec} recovered, AUC={auc:.3f}")

    dataset_label = args.dataset or ''

    for sig in ('cage', 'drna'):
        sig_label = "CAGE (5')" if sig == 'cage' else "3' dRNA"

        if not args.cross_tech_only and sig in by_sig and by_sig[sig]:
            title = f"{sig_label} signal-stratified recall"
            if dataset_label:
                title = f"{dataset_label} — {title}"
            _plot_signal_recall(
                by_sig[sig],
                out / f"signal_recall_{sig}.png",
                title=title,
            )
            if args.verbose:
                print(f"Wrote signal_recall_{sig}.png")

        if sig in cross_tech and cross_tech[sig]:
            n_datasets = max(len(v) for v in cross_tech[sig].values())
            cross_title = (
                f"{sig_label} signal-stratified recall "
                f"(cross-tech mean ± SD, n={n_datasets})"
                if n_datasets > 1
                else f"{sig_label} signal-stratified recall"
            )
            _plot_cross_tech_signal_recall(
                cross_tech[sig],
                out / f"signal_recall_{sig}_cross_tech.png",
                title=cross_title,
            )
            if args.verbose:
                print(f"Wrote signal_recall_{sig}_cross_tech.png")


if __name__ == '__main__':
    main()
