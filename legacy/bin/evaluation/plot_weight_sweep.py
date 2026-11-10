#!/usr/bin/env python3
"""
plot_weight_sweep.py — Pub-quality summary plots for weight-sweep results.

Reads the combined weight_sweep_summary.tsv (produced by CollectSweepResults)
and generates:
  1. weight_sweep_heatmap_precision.png — alpha × profile heatmaps for 5'/3'
     precision, averaged over SQANTI categories (weighted by isoform count).
  2. weight_sweep_pareto.png — Pareto frontier: mean precision vs end redundancy.
  3. weight_sweep_category_sensitivity.png — Per-SQANTI-category precision
     across alpha values, one line per profile.
  4. weight_sweep_redundancy.png — End redundancy by (profile, alpha).

Usage:
  python plot_weight_sweep.py --summary weight_sweep_summary.tsv --outdir .
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

import csv

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

# ─── Pub style defaults ───
DPI = 300
SPINE_VISIBLE = {'top': False, 'right': False}

# Okabe-Ito colorblind-safe palette
PALETTE = [
    '#E69F00',  # orange
    '#56B4E9',  # sky blue
    '#009E73',  # bluish green
    '#F0E442',  # yellow
    '#0072B2',  # blue
    '#D55E00',  # vermillion
    '#CC79A7',  # reddish purple
    '#000000',  # black
]

CATEGORY_ORDER = ['FSM', 'ISM', 'NIC', 'NNC', 'SEM', 'SEN']
CATEGORY_COLORS = {
    'FSM': '#2ecc71', 'ISM': '#3498db', 'NIC': '#f39c12',
    'NNC': '#e74c3c', 'SEM': '#9b59b6', 'SEN': '#95a5a6',
}


def _apply_style(ax):
    """Apply pub-style formatting to an axis."""
    for spine, visible in SPINE_VISIBLE.items():
        ax.spines[spine].set_visible(visible)


def _save(fig, outpath):
    """Save figure with tight layout."""
    fig.savefig(outpath, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: {outpath}")


def parse_summary(tsv_path: Path) -> List[Dict]:
    """Read the combined sweep summary TSV."""
    rows = []
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rows.append(row)
    return rows


def _float(val, default=0.0):
    """Safe float conversion."""
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def _int(val, default=0):
    """Safe int conversion."""
    try:
        return int(val)
    except (ValueError, TypeError):
        return default


# ═══════════════════════════════════════════════════════════════════════════
# Plot 1: Heatmap — alpha × profile → weighted-mean precision
# ═══════════════════════════════════════════════════════════════════════════

def plot_heatmap_precision(rows: List[Dict], outdir: Path):
    """Create heatmaps of weighted-mean 5'/3' precision by (profile, alpha)."""
    # Aggregate: for each (profile, alpha), compute count-weighted mean precision
    agg = defaultdict(lambda: {'w5p': 0.0, 'w3p': 0.0, 'n': 0})
    for row in rows:
        key = (row.get('profile', ''), row.get('alpha', ''))
        count = _int(row.get('count', 1), 1)
        p5 = _float(row.get('5prime_precision', 0))
        p3 = _float(row.get('3prime_precision', 0))
        agg[key]['w5p'] += p5 * count
        agg[key]['w3p'] += p3 * count
        agg[key]['n'] += count

    if not agg:
        return

    profiles = sorted(set(k[0] for k in agg.keys()))
    alphas = sorted(set(float(k[1]) for k in agg.keys()))
    alpha_strs = [str(a) for a in alphas]

    # Build matrices
    mat5 = np.full((len(profiles), len(alphas)), np.nan)
    mat3 = np.full((len(profiles), len(alphas)), np.nan)

    for (prof, alph), vals in agg.items():
        pi = profiles.index(prof)
        ai = alpha_strs.index(alph) if alph in alpha_strs else -1
        if ai < 0:
            continue
        if vals['n'] > 0:
            mat5[pi, ai] = vals['w5p'] / vals['n']
            mat3[pi, ai] = vals['w3p'] / vals['n']

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    for ax, mat, title, cmap in [
        (axes[0], mat5, "5' Weighted-Mean Precision", 'YlGnBu'),
        (axes[1], mat3, "3' Weighted-Mean Precision", 'YlOrRd'),
    ]:
        im = ax.imshow(mat, aspect='auto', cmap=cmap, vmin=0, vmax=1)
        ax.set_xticks(range(len(alphas)))
        ax.set_xticklabels([f'{a:.2f}' for a in alphas], fontsize=7)
        ax.set_yticks(range(len(profiles)))
        ax.set_yticklabels(profiles, fontsize=7)
        ax.set_xlabel('Alpha', fontsize=7)
        ax.set_ylabel('Profile', fontsize=7)
        ax.set_title(title, fontsize=8, fontweight='normal')

        # Annotate cells
        for i in range(len(profiles)):
            for j in range(len(alphas)):
                val = mat[i, j]
                if not np.isnan(val):
                    color = 'white' if val > 0.6 else 'black'
                    ax.text(j, i, f'{val:.3f}', ha='center', va='center',
                            fontsize=8, color=color)

        plt.colorbar(im, ax=ax, shrink=0.8)

    fig.suptitle('Weight Sweep: Precision by Profile × Alpha', fontsize=8, fontweight='normal')
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    _save(fig, outdir / 'weight_sweep_heatmap_precision.png')


# ═══════════════════════════════════════════════════════════════════════════
# Plot 2: Pareto — mean precision vs end redundancy
# ═══════════════════════════════════════════════════════════════════════════

def plot_pareto(rows: List[Dict], outdir: Path):
    """Pareto frontier: mean precision vs redundancy rate."""
    # Aggregate per (profile, alpha)
    agg = defaultdict(lambda: {'w5p': 0.0, 'w3p': 0.0, 'n': 0,
                                'redund_tss': 0, 'redund_tts': 0,
                                'total_groups': 0})
    for row in rows:
        key = (row.get('profile', ''), row.get('alpha', ''))
        count = _int(row.get('count', 1), 1)
        agg[key]['w5p'] += _float(row.get('5prime_precision', 0)) * count
        agg[key]['w3p'] += _float(row.get('3prime_precision', 0)) * count
        agg[key]['n'] += count
        agg[key]['redund_tss'] = _int(row.get('redundant_tss', 0))
        agg[key]['redund_tts'] = _int(row.get('redundant_tts', 0))
        agg[key]['total_groups'] = _int(row.get('total_groups', 0))

    if not agg:
        return

    fig, ax = plt.subplots(figsize=(7, 5.5))
    _apply_style(ax)

    profiles = sorted(set(k[0] for k in agg.keys()))
    profile_colors = {p: PALETTE[i % len(PALETTE)] for i, p in enumerate(profiles)}

    for (prof, alph), vals in agg.items():
        if vals['n'] == 0:
            continue
        mean_prec = (vals['w5p'] + vals['w3p']) / (2 * vals['n'])
        redund_total = vals['redund_tss'] + vals['redund_tts']
        ax.scatter(mean_prec, redund_total,
                   color=profile_colors[prof], s=80, zorder=3,
                   edgecolors='black', linewidths=0.5)
        ax.annotate(f'α={alph}', (mean_prec, redund_total),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=7, color=profile_colors[prof])

    # Legend for profiles
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker='o', color='w',
                       markerfacecolor=profile_colors[p],
                       markeredgecolor='black', markersize=8,
                       label=p)
               for p in profiles]
    ax.legend(handles=handles, fontsize=8, loc='upper right',
              frameon=True, framealpha=0.9)

    ax.set_xlabel('Mean Precision (5\' + 3\') / 2', fontsize=7)
    ax.set_ylabel('Total Redundant End Positions (TSS + TTS)', fontsize=7)
    ax.set_title('Precision vs End Redundancy\n(lower-right = better)',
                 fontsize=8, fontweight='normal')
    ax.grid(True, alpha=0.2)
    _save(fig, outdir / 'weight_sweep_pareto.png')


# ═══════════════════════════════════════════════════════════════════════════
# Plot 3: Category sensitivity — precision per category across alphas
# ═══════════════════════════════════════════════════════════════════════════

def plot_category_sensitivity(rows: List[Dict], outdir: Path):
    """Per-SQANTI-category 5'/3' precision across alpha, one subplot per profile."""
    # Group by (profile, alpha, category)
    data = defaultdict(dict)
    for row in rows:
        prof = row.get('profile', '')
        alpha = _float(row.get('alpha', 0))
        cat = row.get('category', '')
        if cat not in CATEGORY_ORDER:
            continue
        data[(prof, alpha)][cat] = {
            'p5': _float(row.get('5prime_precision', 0)),
            'p3': _float(row.get('3prime_precision', 0)),
        }

    if not data:
        return

    profiles = sorted(set(k[0] for k in data.keys()))
    alphas = sorted(set(k[1] for k in data.keys()))

    n_profiles = len(profiles)
    fig, axes = plt.subplots(n_profiles, 2, figsize=(12, 3 * n_profiles + 1),
                              squeeze=False)

    for pi, prof in enumerate(profiles):
        for col, (end_label, key) in enumerate([("5'", 'p5'), ("3'", 'p3')]):
            ax = axes[pi, col]
            _apply_style(ax)

            for cat in CATEGORY_ORDER:
                vals = []
                for alpha in alphas:
                    cat_data = data.get((prof, alpha), {}).get(cat, {})
                    vals.append(cat_data.get(key, np.nan))

                ax.plot(alphas, vals, marker='o', linewidth=1.5, markersize=5,
                        color=CATEGORY_COLORS.get(cat, '#999'),
                        label=cat)

            ax.set_ylim(-0.05, 1.05)
            ax.set_xlabel('Alpha', fontsize=7)
            ax.set_ylabel('Precision', fontsize=7)
            ax.set_title(f'{prof} — {end_label} Precision', fontsize=7, fontweight='normal')
            ax.grid(True, alpha=0.2)
            if pi == 0 and col == 1:
                ax.legend(fontsize=7, loc='lower left', ncol=2)

    fig.suptitle('Category-Level Precision Sensitivity to Alpha',
                 fontsize=8, fontweight='normal')
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _save(fig, outdir / 'weight_sweep_category_sensitivity.png')


# ═══════════════════════════════════════════════════════════════════════════
# Plot 4: End redundancy by (profile, alpha)
# ═══════════════════════════════════════════════════════════════════════════

def plot_redundancy_grid(rows: List[Dict], outdir: Path):
    """Grouped bar chart showing TSS/TTS redundancy by (profile, alpha)."""
    # Extract unique (profile, alpha) combos
    seen = {}
    for row in rows:
        key = (row.get('profile', ''), row.get('alpha', ''))
        if key not in seen:
            seen[key] = {
                'redund_tss': _int(row.get('redundant_tss', 0)),
                'redund_tts': _int(row.get('redundant_tts', 0)),
                'total_groups': _int(row.get('total_groups', 0)),
            }

    if not seen:
        return

    profiles = sorted(set(k[0] for k in seen.keys()))
    alphas = sorted(set(float(k[1]) for k in seen.keys()))

    fig, ax = plt.subplots(figsize=(max(8, len(seen) * 0.6), 5))
    _apply_style(ax)

    x_labels = []
    tss_vals = []
    tts_vals = []

    for prof in profiles:
        for alpha in alphas:
            key = (prof, str(alpha))
            vals = seen.get(key, {'redund_tss': 0, 'redund_tts': 0})
            x_labels.append(f'{prof}\nα={alpha}')
            tss_vals.append(vals['redund_tss'])
            tts_vals.append(vals['redund_tts'])

    x = np.arange(len(x_labels))
    width = 0.35

    ax.bar(x - width / 2, tss_vals, width, label="TSS redundant",
           color='#3498db', edgecolor='white')
    ax.bar(x + width / 2, tts_vals, width, label="TTS redundant",
           color='#e74c3c', edgecolor='white')

    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, fontsize=7, rotation=45, ha='right')
    ax.set_ylabel('Redundant End Positions', fontsize=7)
    ax.set_title('End Redundancy Across Weight Configurations',
                 fontsize=8, fontweight='normal')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.2, axis='y')
    fig.tight_layout()
    _save(fig, outdir / 'weight_sweep_redundancy.png')


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Generate weight-sweep summary plots from combined TSV")
    parser.add_argument("--summary", required=True, type=Path,
                        help="Combined weight_sweep_summary.tsv")
    parser.add_argument("--outdir", required=True, type=Path,
                        help="Output directory for plots")
    args = parser.parse_args()

    if not HAS_MPL:
        print("ERROR: matplotlib is required for plotting", file=sys.stderr)
        sys.exit(1)

    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = parse_summary(args.summary)
    print(f"Loaded {len(rows)} rows from {args.summary}")

    if not rows:
        print("No data to plot.")
        sys.exit(0)

    # Check if we have the expected columns
    sample_row = rows[0]
    has_redund = 'redundant_tss' in sample_row
    has_category = 'category' in sample_row

    plot_heatmap_precision(rows, args.outdir)
    if has_redund:
        plot_pareto(rows, args.outdir)
        plot_redundancy_grid(rows, args.outdir)
    if has_category:
        plot_category_sensitivity(rows, args.outdir)

    print("Done.")


if __name__ == "__main__":
    main()
