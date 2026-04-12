#!/usr/bin/env python3
"""
peak_roc_curves.py — Signal-stratified peak recovery curves per tool.

Each peak_reason file has columns:
    peak_id  score  read_count  reason

The `score` is the orthogonal signal intensity of the REFERENCE peak —
identical for every tool on any given peak.

Two complementary views are produced for each signal type (CAGE=5', dRNA=3'):

1. Signal-recall curves (existing)
   x = signal score threshold (log scale, high→low)
   y = recall at that threshold
   AUC over normalised threshold axis

2. Normalised signal-recall curves (new)
   Same as (1) but x is normalised to [0,1] across the shared score range
   for that end type, so CAGE and dRNA are directly comparable.

3. TPR-FPR ROC curves (new)
   Treat peak recovery as a binary classifier:
     "Would this tool recover a peak given that it has enough reads?"
   Sweep signal threshold high→low; peaks above threshold are "called":
     TPR(t) = recovered_above_t  / total_recoverable
     FPR(t) = missed_above_t     / total_recoverable
   Diagonal = random classifier. AUC measured vs FPR axis.

4. 5'-vs-3' normalised overlay (new)
   Single plot combining CAGE (5') and dRNA (3') normalised recall
   on the same normalised axis — shows whether a tool is stronger at 5'
   or 3' end recovery.

Outputs (per dataset)
----------------------
  signal_recall_5prime.png               (cage raw signal x-axis)
  signal_recall_3prime.png               (drna raw signal x-axis)
  signal_recall_5prime_norm.png          (normalised x-axis)
  signal_recall_3prime_norm.png          (normalised x-axis)
  signal_roc_5prime.png                  (TPR vs FPR, cage)
  signal_roc_3prime.png                  (TPR vs FPR, drna)
  signal_recall_5v3_overlay.png          (normalised overlay, one panel per tool)

Cross-tech aggregates (mean ± SD when >1 dataset)
  signal_recall_5prime_cross_tech.png
  signal_recall_3prime_cross_tech.png
  signal_recall_5prime_norm_cross_tech.png
  signal_recall_3prime_norm_cross_tech.png
  signal_roc_5prime_cross_tech.png
  signal_roc_3prime_cross_tech.png
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

# ── constants ─────────────────────────────────────────────────────────────────

_RECOVERED = "recovered"
_TOOL_RE = re.compile(
    r'(?:^|/)([^/]+?)_transcriptome_(?:cage|drna)_peak_reasons\.tsv$'
)
_SIG_RE = re.compile(r'_transcriptome_(cage|drna)_peak_reasons\.tsv$')

# Canonical end-type names for both raw signal type labels and file suffixes
_SIG_END = {'cage': '5prime', 'drna': '3prime'}
_END_LABEL = {'5prime': "5′ CAGE", '3prime': "3′ dRNA"}


# ── file parsing ──────────────────────────────────────────────────────────────

def _tool_name(path: Path) -> str:
    m = _TOOL_RE.search(str(path))
    if m:
        full = m.group(1)
        # Strip dataset/align/partition prefix — keep only mode (last underscore segment)
        parts = full.rsplit('_', 1)
        return parts[-1] if len(parts) > 1 else full
    return path.stem


def _sig_type(path: Path) -> str:
    """Return canonical signal type: 'cage' or 'drna'."""
    m = _SIG_RE.search(str(path))
    return m.group(1) if m else 'unknown'


def _load(path: Path):
    """Return list of (score, is_recovered) for all non-header rows."""
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


# ── curve computation ─────────────────────────────────────────────────────────

def _signal_recall_curve(rows):
    """Signal-threshold recall curve (x = raw score threshold, y = recall)."""
    if not rows:
        return np.array([0.0]), np.array([0.0])
    scores = np.array([r[0] for r in rows])
    is_rec = np.array([r[1] for r in rows], dtype=bool)
    thresholds = np.unique(scores)[::-1]   # high → low
    recalls = []
    for t in thresholds:
        mask = scores >= t
        total = mask.sum()
        recalls.append((mask & is_rec).sum() / total if total > 0 else 0.0)
    return thresholds, np.array(recalls)


def _signal_recall_curve_norm(rows, t_min=None, t_max=None):
    """Same curve but x-axis normalised to [0, 1] over the provided score range."""
    thresholds, recalls = _signal_recall_curve(rows)
    if t_min is None:
        t_min = thresholds.min()
    if t_max is None:
        t_max = thresholds.max()
    if t_max > t_min:
        norm = (thresholds - t_min) / (t_max - t_min)
    else:
        norm = np.zeros_like(thresholds)
    return norm, recalls


def _tpr_fpr_roc_curve(rows):
    """ROC curve: TPR vs FPR by sweeping signal threshold.

    At each threshold t (high→low), peaks with score >= t are "called positive".
    Among recoverable peaks:
        TPR(t) = recovered_at_or_above_t / total_recoverable
        FPR(t) = missed_at_or_above_t   / total_recoverable

    FPR here measures how many recoverable peaks are *missed* at each
    threshold — i.e. called positive but not actually recovered.
    This is more interpretable than a traditional FPR because every peak
    is a potential TP (there are no true negatives in the classical sense).
    """
    if not rows:
        return np.array([0.0, 1.0]), np.array([0.0, 0.0]), 0.0

    scores = np.array([r[0] for r in rows])
    is_rec = np.array([r[1] for r in rows], dtype=bool)
    total_recoverable = is_rec.sum()
    if total_recoverable == 0:
        return np.array([0.0, 1.0]), np.array([0.0, 0.0]), 0.0

    thresholds = np.unique(scores)[::-1]
    tprs, fprs = [], []
    for t in thresholds:
        above = scores >= t
        tpr = (above & is_rec).sum() / total_recoverable
        fpr = (above & ~is_rec).sum() / max((~is_rec).sum(), 1)
        tprs.append(tpr)
        fprs.append(fpr)

    fprs = np.array(fprs)
    tprs = np.array(tprs)
    # Sort by FPR for AUC integration
    order = np.argsort(fprs)
    auc = float(np.trapz(tprs[order], fprs[order]))
    return fprs, tprs, auc


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


def _interpolate(x_src, y_src, x_grid):
    """Interpolate y_src at x_grid points."""
    order = np.argsort(x_src)
    return np.interp(x_grid, x_src[order], y_src[order],
                     left=y_src[order][0], right=y_src[order][-1])


# ── plotting helpers ──────────────────────────────────────────────────────────

def _make_styler_and_fig(tools, figsize=(3.8, 3.5)):
    apply_rc()
    styler = ModeStyler(tools)
    fig, ax = plt.subplots(figsize=figsize)
    return styler, fig, ax


def _finish(fig, ax, handles, output_path, xlabel, ylabel, title,
            xlim=None, ylim=None, xscale=None, invert_x=False):
    style_ax(ax, xlabel=xlabel, ylabel=ylabel, title=title)
    if xscale:
        ax.set_xscale(xscale)
    if invert_x:
        ax.invert_xaxis()
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    legend_outside(fig, handles=handles, loc='upper left',
                   bbox_to_anchor=(1.02, 1.0), ncol=1, fontsize=6.5)
    fig.tight_layout()
    savefig(fig, output_path, dpi=300)


# ── per-dataset plots ─────────────────────────────────────────────────────────

def plot_signal_recall(curves_by_tool, output_path, title, normalised=False,
                       shared_range=None):
    """Signal-threshold recall curve.

    Args:
        curves_by_tool: dict[tool] -> (thresholds, recalls)
        normalised: if True, normalise x-axis to [0,1] using shared_range
        shared_range: (t_min, t_max) for normalisation
    """
    tools = sorted(curves_by_tool)
    if not tools:
        return
    styler, fig, ax = _make_styler_and_fig(tools)
    handles = []

    for tool in tools:
        thresholds, recalls = curves_by_tool[tool]
        if normalised and shared_range:
            t_min, t_max = shared_range
            if t_max > t_min:
                x = (thresholds - t_min) / (t_max - t_min)
            else:
                x = np.zeros_like(thresholds)
        else:
            x = thresholds

        line, = ax.plot(x, recalls * 100,
                        color=styler.color(tool), linewidth=1.2, alpha=0.85)
        line.set_dashes(styler.dash(tool))
        step = max(1, len(x) // 8)
        ax.plot(x[::step], (recalls * 100)[::step],
                color=styler.color(tool), marker=styler.marker(tool),
                markersize=4, linestyle='', alpha=0.9)

        if normalised:
            auc = _auc_threshold(x, recalls)
        else:
            auc = _auc_threshold(thresholds, recalls)
        handles.append(styler.legend_handle(
            tool, label=f"{tool}  (AUC={auc:.3f})", markersize=6))

    if normalised:
        xlabel = "Normalised signal threshold (high→low)"
        xscale = None
        invert_x = True
        xlim = (-0.02, 1.02)
    else:
        xlabel = "Signal score threshold (log, high→low)"
        xscale = 'log'
        invert_x = True
        xlim = None

    _finish(fig, ax, handles, output_path,
            xlabel=xlabel,
            ylabel="Recall at threshold (%)",
            title=title,
            xlim=xlim, ylim=(-2, 102),
            xscale=xscale, invert_x=invert_x)


def plot_tpr_fpr_roc(roc_by_tool, output_path, title):
    """TPR-vs-FPR ROC curve.

    Args:
        roc_by_tool: dict[tool] -> (fprs, tprs, auc)
    """
    tools = sorted(roc_by_tool)
    if not tools:
        return
    styler, fig, ax = _make_styler_and_fig(tools)

    # Diagonal reference line
    ax.plot([0, 1], [0, 1], '--', color='#cccccc', linewidth=0.8, zorder=0)

    handles = []
    for tool in tools:
        fprs, tprs, auc = roc_by_tool[tool]
        line, = ax.plot(fprs, tprs * 100,
                        color=styler.color(tool), linewidth=1.2, alpha=0.85)
        line.set_dashes(styler.dash(tool))
        step = max(1, len(fprs) // 8)
        ax.plot(fprs[::step], (tprs * 100)[::step],
                color=styler.color(tool), marker=styler.marker(tool),
                markersize=4, linestyle='', alpha=0.9)
        handles.append(styler.legend_handle(
            tool, label=f"{tool}  (AUC={auc:.3f})", markersize=6))

    _finish(fig, ax, handles, output_path,
            xlabel="FPR — Missed recoverable peaks (%)",
            ylabel="TPR — Recovered peaks (%)",
            title=title,
            xlim=(-2, 102), ylim=(-2, 102))


def plot_5v3_overlay(curves_5, curves_3, output_path, title, shared_range_5,
                     shared_range_3):
    """Overlay 5' (CAGE) and 3' (dRNA) normalised recall per tool.

    One subplot per tool, all on normalised [0,1] x-axis.
    5' shown solid, 3' shown dashed.  Lets you see which end a tool
    does better at across the signal range.
    """
    all_tools = sorted(set(curves_5) | set(curves_3))
    if not all_tools:
        return

    n = len(all_tools)
    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    apply_rc()
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols * 2.8, nrows * 2.6),
                              squeeze=False)

    for idx, tool in enumerate(all_tools):
        ax = axes[idx // ncols][idx % ncols]
        color = ModeStyler(all_tools).color(tool)

        for curves, label, ls, end, sr in [
            (curves_5, "5′ CAGE",        (1, 0), '5prime', shared_range_5),
            (curves_3, "3′ dRNA",    (4, 2), '3prime', shared_range_3),
        ]:
            if tool not in curves:
                continue
            thresholds, recalls = curves[tool]
            t_min, t_max = sr if sr else (thresholds.min(), thresholds.max())
            if t_max > t_min:
                x = (thresholds - t_min) / (t_max - t_min)
            else:
                x = np.linspace(0, 1, len(thresholds))
            line, = ax.plot(x, recalls * 100, color=color, linewidth=1.3, alpha=0.9,
                            label=label)
            line.set_dashes(ls)

        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-2, 102)
        ax.invert_xaxis()
        style_ax(ax, title=tool, xlabel="Norm. signal (high→low)", ylabel="Recall (%)")
        ax.legend(fontsize=6, loc='lower left')

    # Hide unused axes
    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle(title, fontsize=8, y=1.01)
    fig.tight_layout()
    savefig(fig, output_path, dpi=300)


# ── cross-tech aggregates ─────────────────────────────────────────────────────

def plot_cross_tech_signal_recall(per_tool_per_dataset, output_path, title,
                                  normalised=False, shared_range=None,
                                  n_grid=300):
    """Mean ± SD signal recall across multiple datasets."""
    all_thresholds = []
    for curves in per_tool_per_dataset.values():
        for t, _ in curves:
            all_thresholds.extend(t.tolist())
    if not all_thresholds:
        return

    t_min = max(1e-6, min(all_thresholds))
    t_max = max(all_thresholds)

    if normalised and shared_range:
        sr_min, sr_max = shared_range
        grid = np.linspace(0.0, 1.0, n_grid)
        to_x = lambda t: (t - sr_min) / (sr_max - sr_min) if sr_max > sr_min else np.zeros_like(t)
    else:
        grid = np.logspace(np.log10(t_min), np.log10(t_max), n_grid)
        to_x = lambda t: t

    tools = sorted(per_tool_per_dataset)
    apply_rc()
    styler = ModeStyler(tools)
    fig, ax = plt.subplots(figsize=(3.8, 3.5))
    handles = []

    for tool in tools:
        curves = per_tool_per_dataset[tool]
        if not curves:
            continue
        interped = np.array([
            _interpolate(to_x(t), r, grid)
            for t, r in curves
        ])
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

        auc = _auc_threshold(grid, mean_r)
        handles.append(styler.legend_handle(
            tool, label=f"{tool}  (AUC={auc:.3f})", markersize=6))

    if normalised:
        xlabel = "Normalised signal threshold (high→low)"
        xscale = None
        invert_x = True
        xlim = (-0.02, 1.02)
    else:
        xlabel = "Signal score threshold (log, high→low)"
        xscale = 'log'
        invert_x = True
        xlim = None

    _finish(fig, ax, handles, output_path,
            xlabel=xlabel,
            ylabel="Recall at threshold (%)",
            title=title,
            xlim=xlim, ylim=(-2, 102),
            xscale=xscale, invert_x=invert_x)


def plot_cross_tech_roc(per_tool_per_dataset, output_path, title, n_grid=200):
    """Mean ± SD TPR-FPR ROC across multiple datasets."""
    tools = sorted(per_tool_per_dataset)
    apply_rc()
    styler = ModeStyler(tools)
    fig, ax = plt.subplots(figsize=(3.8, 3.5))
    ax.plot([0, 1], [0, 1], '--', color='#cccccc', linewidth=0.8, zorder=0)
    handles = []
    fpr_grid = np.linspace(0, 1, n_grid)

    for tool in tools:
        curves = per_tool_per_dataset[tool]
        if not curves:
            continue
        # curves is list of (fprs, tprs, auc)
        interped = np.array([
            _interpolate(fprs, tprs, fpr_grid)
            for fprs, tprs, _ in curves
        ])
        mean_tpr = interped.mean(axis=0)
        std_tpr  = interped.std(axis=0)
        auc = float(np.trapz(mean_tpr, fpr_grid))

        line, = ax.plot(fpr_grid * 100, mean_tpr * 100,
                        color=styler.color(tool), linewidth=1.2, alpha=0.9)
        line.set_dashes(styler.dash(tool))
        if len(curves) > 1:
            ax.fill_between(fpr_grid * 100,
                            np.clip((mean_tpr - std_tpr) * 100, 0, 100),
                            np.clip((mean_tpr + std_tpr) * 100, 0, 100),
                            color=styler.color(tool), alpha=0.15)
        handles.append(styler.legend_handle(
            tool, label=f"{tool}  (AUC={auc:.3f})", markersize=6))

    _finish(fig, ax, handles, output_path,
            xlabel="FPR — Missed recoverable peaks (%)",
            ylabel="TPR — Recovered peaks (%)",
            title=title,
            xlim=(-2, 102), ylim=(-2, 102))


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Signal-stratified peak recovery + ROC curves per tool"
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

    # Group raw data by end type (5prime / 3prime) and tool
    # raw_curves[end_type][tool] = (thresholds, recalls)
    raw_curves:  dict = defaultdict(dict)
    roc_curves:  dict = defaultdict(dict)
    # For cross-tech aggregation:
    # ct_recall[end_type][tool] = list of (thresholds, recalls)
    # ct_roc[end_type][tool]    = list of (fprs, tprs, auc)
    ct_recall: dict = defaultdict(lambda: defaultdict(list))
    ct_roc:    dict = defaultdict(lambda: defaultdict(list))

    for f in args.input:
        p = Path(f)
        if not p.exists():
            print(f"Warning: {f} not found, skipping", file=sys.stderr)
            continue
        sig_raw = _sig_type(p)          # 'cage' or 'drna'
        end     = _SIG_END.get(sig_raw, sig_raw)   # '5prime' or '3prime'
        tool    = _tool_name(p)
        rows    = _load(p)
        if not rows:
            continue

        thresholds, recalls = _signal_recall_curve(rows)
        fprs, tprs, auc     = _tpr_fpr_roc_curve(rows)

        raw_curves[end][tool] = (thresholds, recalls)
        roc_curves[end][tool] = (fprs, tprs, auc)
        ct_recall[end][tool].append((thresholds, recalls))
        ct_roc[end][tool].append((fprs, tprs, auc))

        if args.verbose:
            n_rec = sum(1 for _, r in rows if r)
            print(f"  {tool} ({sig_raw}/{end}): {len(rows)} peaks, "
                  f"{n_rec} recovered, recall-AUC={_auc_threshold(thresholds, recalls):.3f}, "
                  f"ROC-AUC={auc:.3f}")

    dataset_label = args.dataset or ''

    # Compute shared score range per end type for normalisation
    shared_range: dict = {}
    for end in ('5prime', '3prime'):
        all_t = [t for tool in raw_curves[end]
                   for t in raw_curves[end][tool][0].tolist()]
        if all_t:
            shared_range[end] = (min(all_t), max(all_t))
        else:
            shared_range[end] = (0.0, 1.0)

    end_label_map = {
        '5prime': "5′ CAGE",
        '3prime': "3′ dRNA",
    }

    for end in ('5prime', '3prime'):
        end_label = end_label_map.get(end, end)
        sr = shared_range.get(end)

        def _title(base):
            return f"{dataset_label} — {base}" if dataset_label else base

        # ── per-dataset plots ──────────────────────────────────────────────

        if not args.cross_tech_only and raw_curves[end]:

            # 1. Raw-signal recall
            plot_signal_recall(
                raw_curves[end],
                out / f"signal_recall_{end}.png",
                title=_title(f"{end_label} signal-stratified recall"),
                normalised=False,
            )

            # 2. Normalised-signal recall
            plot_signal_recall(
                raw_curves[end],
                out / f"signal_recall_{end}_norm.png",
                title=_title(f"{end_label} normalised signal recall"),
                normalised=True,
                shared_range=sr,
            )

            # 3. TPR-FPR ROC
            plot_tpr_fpr_roc(
                roc_curves[end],
                out / f"signal_roc_{end}.png",
                title=_title(f"{end_label} peak recovery ROC"),
            )

            if args.verbose:
                print(f"Wrote 3 plots for {end}")

        # ── cross-tech aggregates ──────────────────────────────────────────

        if ct_recall[end]:
            n_ds = max(len(v) for v in ct_recall[end].values())
            suffix = f" (n={n_ds})" if n_ds > 1 else ""

            plot_cross_tech_signal_recall(
                ct_recall[end],
                out / f"signal_recall_{end}_cross_tech.png",
                title=_title(f"{end_label} signal recall{suffix}"),
                normalised=False,
            )
            plot_cross_tech_signal_recall(
                ct_recall[end],
                out / f"signal_recall_{end}_norm_cross_tech.png",
                title=_title(f"{end_label} normalised signal recall{suffix}"),
                normalised=True,
                shared_range=sr,
            )
            plot_cross_tech_roc(
                ct_roc[end],
                out / f"signal_roc_{end}_cross_tech.png",
                title=_title(f"{end_label} ROC{suffix}"),
            )

            if args.verbose:
                print(f"Wrote cross-tech plots for {end}")

    # ── 5'-vs-3' overlay ──────────────────────────────────────────────────────
    if (not args.cross_tech_only
            and raw_curves.get('5prime') and raw_curves.get('3prime')):
        plot_5v3_overlay(
            raw_curves['5prime'],
            raw_curves['3prime'],
            out / "signal_recall_5v3_overlay.png",
            title=_title("5′ vs 3′ normalised recall per tool"),
            shared_range_5=shared_range.get('5prime'),
            shared_range_3=shared_range.get('3prime'),
        )
        if args.verbose:
            print("Wrote signal_recall_5v3_overlay.png")


if __name__ == '__main__':
    main()
