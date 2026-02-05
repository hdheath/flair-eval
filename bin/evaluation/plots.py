"""
Visualization functions for evaluation plots.

Provides functions for generating distance histograms, entropy distributions,
read classification summaries, truncation pattern plots, and sequence logos.
"""

import statistics
from pathlib import Path
from typing import Dict, List, Optional

from .utils import get_logger
from .motif import compute_information_content

logger = get_logger()

# Optional matplotlib/numpy import for plotting
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend for server environments
    import matplotlib.pyplot as plt
    import numpy as np
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    np = None

# Golden ratio for aesthetically proportioned figures
GOLDEN_RATIO = 1.618
PLOT_DPI = 300


def calculate_histogram_max_count(
    distances: List[int],
    bin_size: int = 50,
    min_dist: int = -1000,
    max_dist: int = 1000,
) -> int:
    """Calculate the maximum bin count for a histogram.
    
    Args:
        distances: List of signed distances
        bin_size: Size of each bin in bp
        min_dist: Minimum distance for binning
        max_dist: Maximum distance for binning
        
    Returns:
        Maximum count in any bin
    """
    if not distances:
        return 0
    
    # Clamp distances to range
    clamped = [max(min_dist, min(max_dist, d)) for d in distances]
    
    # Create bins
    bins = list(range(min_dist, max_dist + bin_size, bin_size))
    
    # Count values in each bin
    counts, _ = np.histogram(clamped, bins=bins)
    
    return int(max(counts)) if len(counts) > 0 else 0


def plot_distance_histogram(
    distances: List[int],
    output_path: Path,
    title: str,
    bin_size: int = 50,
    min_dist: int = -1000,
    max_dist: int = 1000,
    max_count: Optional[int] = None,
) -> bool:
    """
    Create a histogram of signed distances to experimental peaks.

    Args:
        distances: List of signed distances
        output_path: Path to save the plot
        title: Plot title
        bin_size: Size of each bin in bp (default: 50)
        min_dist: Minimum distance for x-axis (default: -1000)
        max_dist: Maximum distance for x-axis (default: 1000)
        max_count: Fixed y-axis limit for consistent comparison across runs (optional)

    Returns:
        True if plot was created successfully, False otherwise
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available; cannot create distance histogram")
        return False

    if not distances:
        logger.warning(f"No distances to plot for {title}")
        return False

    # Clamp distances into overflow bins at the edges
    # Distances beyond min_dist/max_dist are collected into the outermost bins
    clamped = [max(min_dist, min(max_dist, d)) for d in distances]
    n_clamped_low = sum(1 for d in distances if d < min_dist)
    n_clamped_high = sum(1 for d in distances if d > max_dist)

    # Create figure with golden ratio proportions
    fig_width = 10
    fig, ax = plt.subplots(figsize=(fig_width, fig_width / GOLDEN_RATIO))

    # Calculate bins
    bins = list(range(min_dist, max_dist + bin_size, bin_size))

    # Create histogram with clamped distances so outliers appear in edge bins
    ax.hist(clamped, bins=bins, edgecolor='black', alpha=0.7, color='steelblue')

    # Add vertical line at 0
    ax.axvline(x=0, color='red', linestyle='--', linewidth=1.5, label='Perfect alignment')

    # Labels and title
    ax.set_xlabel('Distance to Nearest Peak (bp)', fontsize=12)
    ax.set_ylabel('Number of Transcripts', fontsize=12)
    ax.set_title(title, fontsize=14)

    # Add statistics annotation
    mean_dist = statistics.mean(distances)
    median_dist = statistics.median(distances)
    std_dist = statistics.stdev(distances) if len(distances) > 1 else 0

    stats_text = f'n = {len(distances)}\nMean = {mean_dist:.1f} bp\nMedian = {median_dist:.1f} bp\nStd = {std_dist:.1f} bp'
    if n_clamped_low or n_clamped_high:
        stats_text += f'\n<{min_dist}: {n_clamped_low}  >{max_dist}: {n_clamped_high}'
    ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Set x-axis limits
    ax.set_xlim(min_dist, max_dist)
    
    # Set y-axis limit if provided (for consistent comparison across runs)
    if max_count is not None:
        ax.set_ylim(0, max_count)

    # Add legend
    ax.legend(loc='upper left')

    # Tight layout
    plt.tight_layout()

    # Save figure
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=PLOT_DPI, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved distance histogram to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to save histogram to {output_path}: {e}")
        plt.close(fig)
        return False


def plot_distance_histogram_colored(
    distances: List[int],
    recoverable_mask: List[bool],
    output_path: Path,
    title: str,
    bin_size: int = 50,
    min_dist: int = -1000,
    max_dist: int = 1000,
    max_count: Optional[int] = None,
) -> bool:
    """Create a histogram with bars colored by peak recoverability.

    Recoverable = at least one long read end is within the window of the peak.

    Args:
        distances: List of signed distances
        recoverable_mask: List of boolean flags indicating recoverability
        output_path: Path to save the plot
        title: Plot title
        bin_size: Size of each bin in bp (default: 50)
        min_dist: Minimum distance for x-axis (default: -1000)
        max_dist: Maximum distance for x-axis (default: 1000)
        max_count: Fixed y-axis limit for consistent comparison across runs (optional)

    Returns:
        True if plot was created successfully, False otherwise
    """
    if not HAS_MATPLOTLIB:
        return False
    if not distances:
        return False
    if len(distances) != len(recoverable_mask):
        logger.warning(f"Distance/mask length mismatch: {len(distances)} vs {len(recoverable_mask)}, "
                       f"falling back to uncolored plot")
        return plot_distance_histogram(distances, output_path, title, bin_size, min_dist, max_dist, max_count)

    # Split distances by recoverability
    dists_recoverable = [d for d, r in zip(distances, recoverable_mask) if r]
    dists_unrecoverable = [d for d, r in zip(distances, recoverable_mask) if not r]

    # Clamp
    clamp = lambda d: max(min_dist, min(max_dist, d))
    clamped_rec = [clamp(d) for d in dists_recoverable]
    clamped_unrec = [clamp(d) for d in dists_unrecoverable]

    n_clamped_low = sum(1 for d in distances if d < min_dist)
    n_clamped_high = sum(1 for d in distances if d > max_dist)

    fig_width = 10
    fig, ax = plt.subplots(figsize=(fig_width, fig_width / GOLDEN_RATIO))
    bins = list(range(min_dist, max_dist + bin_size, bin_size))

    # Stacked histogram: recoverable on bottom, unrecoverable on top
    ax.hist([clamped_rec, clamped_unrec], bins=bins, stacked=True,
            edgecolor='black', alpha=0.7,
            color=['steelblue', 'lightcoral'],
            label=[f'Read-supported peak ({len(dists_recoverable)})',
                   f'No read support ({len(dists_unrecoverable)})'])

    ax.axvline(x=0, color='red', linestyle='--', linewidth=1.5, label='Perfect alignment')
    ax.set_xlabel('Distance to Nearest Peak (bp)', fontsize=12)
    ax.set_ylabel('Number of Transcripts', fontsize=12)
    ax.set_title(title, fontsize=14)

    # Stats on full distribution
    mean_dist = statistics.mean(distances)
    median_dist = statistics.median(distances)
    std_dist = statistics.stdev(distances) if len(distances) > 1 else 0
    stats_text = (f'n = {len(distances)}\nMean = {mean_dist:.1f} bp\n'
                  f'Median = {median_dist:.1f} bp\nStd = {std_dist:.1f} bp')
    if n_clamped_low or n_clamped_high:
        stats_text += f'\n<{min_dist}: {n_clamped_low}  >{max_dist}: {n_clamped_high}'
    ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.set_xlim(min_dist, max_dist)
    
    # Set y-axis limit if provided (for consistent comparison across runs)
    if max_count is not None:
        ax.set_ylim(0, max_count)
    
    ax.legend(loc='upper left')
    plt.tight_layout()

    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=PLOT_DPI, bbox_inches='tight')
        plt.close(fig)
        return True
    except Exception as e:
        logger.error(f"Failed to save colored histogram: {e}")
        plt.close(fig)
        return False


def plot_read_end_entropy(
    entropy_data: dict,
    plot_output_dir: Path,
    plot_prefix: str,
) -> None:
    """Generate read-end entropy plots:
    1. Aggregate histogram of read-to-model offsets (TSS and TTS)
    2. Distribution of per-isoform entropies (TSS and TTS)
    """
    if not HAS_MATPLOTLIB:
        return

    plot_output_dir.mkdir(parents=True, exist_ok=True)

    # Plot 1a: Aggregate TSS offset histogram
    tss_offsets = entropy_data["all_tss_offsets"]
    if tss_offsets:
        plot_distance_histogram(
            distances=tss_offsets,
            output_path=plot_output_dir / f"{plot_prefix}_read_tss_offset_histogram.png",
            title="Read 5' End Offset from Isoform TSS",
        )

    # Plot 1b: Aggregate TTS offset histogram
    tts_offsets = entropy_data["all_tts_offsets"]
    if tts_offsets:
        plot_distance_histogram(
            distances=tts_offsets,
            output_path=plot_output_dir / f"{plot_prefix}_read_tts_offset_histogram.png",
            title="Read 3' End Offset from Isoform TTS",
        )

    # Plot 2a: Per-isoform TSS entropy distribution
    tss_entropies = [e for _, e, _ in entropy_data["tss_entropy_per_isoform"]]
    if tss_entropies:
        fig_width = 10
        fig, ax = plt.subplots(figsize=(fig_width, fig_width / GOLDEN_RATIO))
        ax.hist(tss_entropies, bins=50, edgecolor='black', alpha=0.7, color='coral')
        ax.set_xlabel("Shannon Entropy (bits, 10bp bins)", fontsize=12)
        ax.set_ylabel("Number of Isoforms", fontsize=12)
        ax.set_title("Per-Isoform 5' Read-End Entropy", fontsize=14)
        mean_e = statistics.mean(tss_entropies)
        median_e = statistics.median(tss_entropies)
        stats_text = f'n = {len(tss_entropies)}\nMean = {mean_e:.2f} bits\nMedian = {median_e:.2f} bits'
        ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.tight_layout()
        plt.savefig(plot_output_dir / f"{plot_prefix}_tss_entropy_distribution.png",
                    dpi=PLOT_DPI, bbox_inches='tight')
        plt.close(fig)

    # Plot 2b: Per-isoform TTS entropy distribution
    tts_entropies = [e for _, e, _ in entropy_data["tts_entropy_per_isoform"]]
    if tts_entropies:
        fig_width = 10
        fig, ax = plt.subplots(figsize=(fig_width, fig_width / GOLDEN_RATIO))
        ax.hist(tts_entropies, bins=50, edgecolor='black', alpha=0.7, color='mediumpurple')
        ax.set_xlabel("Shannon Entropy (bits, 10bp bins)", fontsize=12)
        ax.set_ylabel("Number of Isoforms", fontsize=12)
        ax.set_title("Per-Isoform 3' Read-End Entropy", fontsize=14)
        mean_e = statistics.mean(tts_entropies)
        median_e = statistics.median(tts_entropies)
        stats_text = f'n = {len(tts_entropies)}\nMean = {mean_e:.2f} bits\nMedian = {median_e:.2f} bits'
        ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.tight_layout()
        plt.savefig(plot_output_dir / f"{plot_prefix}_tts_entropy_distribution.png",
                    dpi=PLOT_DPI, bbox_inches='tight')
        plt.close(fig)


def plot_read_support_distribution(
    recoverable_counts: Dict[str, int],
    output_path: Path,
    title: str,
) -> bool:
    """Plot read support counts for recoverable peaks.

    Displays:
      - 0-100: unit-width bins (one bar per integer)
      - 100+: a single overflow bar whose height equals the number of peaks with support >100
    """
    if not HAS_MATPLOTLIB:
        return False
    if not recoverable_counts:
        return False

    counts = list(recoverable_counts.values())

    fig_width = 12
    fig, ax = plt.subplots(figsize=(fig_width, fig_width / GOLDEN_RATIO))

    # Build frequency table for 0..100 and overflow (100+)
    max_bin = 100
    freq = [0] * (max_bin + 1)  # index i = count i
    overflow = 0

    for c in counts:
        if c < 0:
            continue
        if c <= max_bin:
            freq[c] += 1
        else:
            overflow += 1

    # X positions: 0..100, plus overflow bar with a visible gap from the main histogram
    x_main = list(range(0, max_bin + 1))
    y_main = freq

    overflow_x = max_bin + 5  # 105 — gap separates it from the 100 bar
    ax.bar(x_main, y_main, width=1.0, edgecolor='black', alpha=0.7)
    ax.bar([overflow_x], [overflow], width=1.0, edgecolor='black', alpha=0.7)

    ax.set_xlabel('Number of Supporting Long Reads', fontsize=12)
    ax.set_ylabel('Number of Peaks', fontsize=12)
    ax.set_title(title, fontsize=14)

    # Ticks: keep sparse for readability, and label overflow as 100+
    tick_positions = [0, 1, 2, 3, 4, 5, 10, 20, 30, 50, 75, 100, overflow_x]
    tick_labels = [str(t) for t in tick_positions[:-1]] + ['100+']
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)

    # Limits so the overflow bar is fully visible
    ax.set_xlim(-0.5, overflow_x + 1.5)

    # Summary stats
    mean_c = statistics.mean(counts)
    median_c = statistics.median(counts)
    std_c = statistics.stdev(counts) if len(counts) > 1 else 0
    stats_text = (
        f'n = {len(counts)} peaks\n'
        f'Mean = {mean_c:.1f} reads\n'
        f'Median = {median_c:.1f} reads\n'
        f'Std = {std_c:.1f}\n'
        f'>100: {overflow}'
    )
    ax.text(
        0.98, 0.98, stats_text,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment='top',
        horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
    )

    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=PLOT_DPI, bbox_inches='tight')
        plt.close(fig)
        return True
    except Exception as e:
        logger.error(f"Failed to save read support histogram: {e}")
        plt.close(fig)
        return False


def plot_read_classification_summary(
    classification_summary: Dict[str, int],
    output_path: Path,
    title: str,
) -> bool:
    """Plot bar chart of read classification categories for missed peaks."""
    if not HAS_MATPLOTLIB:
        return False

    if not classification_summary or sum(classification_summary.values()) == 0:
        return False

    fig_width = 8
    fig, ax = plt.subplots(figsize=(fig_width, fig_width / GOLDEN_RATIO))

    category_order = ['unassigned', 'assigned_nearby', 'assigned_distant', 'assigned_wrong_strand']
    category_labels = {
        'unassigned': 'Unassigned',
        'assigned_nearby': 'Nearby',
        'assigned_distant': 'Distant',
        'assigned_wrong_strand': 'Wrong Strand',
    }
    colors = ['#E57373', '#64B5F6', '#81C784', '#FFB74D']

    bar_labels = []
    counts = []
    bar_colors = []
    for i, cat in enumerate(category_order):
        count = classification_summary.get(cat, 0)
        if count > 0:
            bar_labels.append(category_labels[cat])
            counts.append(count)
            bar_colors.append(colors[i])

    if not counts:
        plt.close(fig)
        return False

    total = sum(counts)
    bars = ax.bar(bar_labels, counts, color=bar_colors, edgecolor='black', alpha=0.85)
    for bar, count in zip(bars, counts):
        pct = 100.0 * count / total
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + total * 0.01,
                f"{count}\n({pct:.1f}%)", ha='center', va='bottom', fontsize=10)

    ax.set_xlabel('Read Classification', fontsize=11)
    ax.set_ylabel('Number of Peaks', fontsize=11)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax.set_axisbelow(True)

    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=PLOT_DPI, bbox_inches='tight')
        plt.close(fig)
        return True
    except Exception as e:
        logger.error(f"Failed to save classification bar chart: {e}")
        plt.close(fig)
        return False


def plot_truncation_patterns(
    truncation_patterns: Dict[str, int],
    output_path: Path,
    title: str,
) -> bool:
    """Plot bar chart of truncation pattern categories."""
    if not HAS_MATPLOTLIB:
        return False

    if not truncation_patterns or sum(truncation_patterns.values()) == 0:
        return False

    fig_width = 10
    fig, ax = plt.subplots(figsize=(fig_width, fig_width / GOLDEN_RATIO))

    pattern_order = ['sharp', 'trailing', 'bimodal', 'dispersed', 'sparse']
    colors = ['#2ecc71', '#e74c3c', '#9b59b6', '#3498db', '#95a5a6']

    patterns = []
    counts = []
    bar_colors = []

    for i, pattern in enumerate(pattern_order):
        if pattern in truncation_patterns:
            patterns.append(pattern.capitalize())
            counts.append(truncation_patterns[pattern])
            bar_colors.append(colors[i])

    if not patterns:
        plt.close(fig)
        return False

    bars = ax.bar(patterns, counts, color=bar_colors, edgecolor='black', alpha=0.8)

    # Add count labels on bars
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(count), ha='center', va='bottom', fontsize=10)

    ax.set_xlabel('Truncation Pattern', fontsize=12)
    ax.set_ylabel('Number of Peaks', fontsize=12)
    ax.set_title(title, fontsize=14)

    # Add descriptions
    descriptions = {
        'Sharp': 'Reads cluster at peak',
        'Trailing': 'RT truncation gradient',
        'Bimodal': 'Two clusters (alt TSS?)',
        'Dispersed': 'Spread out reads',
        'Sparse': '<5 reads',
    }
    desc_text = '\n'.join([f"{k}: {v}" for k, v in descriptions.items() if k in patterns])
    ax.text(0.98, 0.98, desc_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=PLOT_DPI, bbox_inches='tight')
        plt.close(fig)
        return True
    except Exception as e:
        logger.error(f"Failed to save truncation patterns plot: {e}")
        plt.close(fig)
        return False


def plot_all_truncation_patterns(
    peak_patterns: Dict[str, dict],
    output_path: Path,
    title: str,
) -> bool:
    """Plot bar chart of truncation patterns for ALL recoverable peaks, split by captured/missed."""
    if not HAS_MATPLOTLIB:
        return False

    if not peak_patterns:
        return False

    from collections import defaultdict

    # Count patterns separately for captured vs missed
    captured_patterns = defaultdict(int)
    missed_patterns = defaultdict(int)

    for peak_id, info in peak_patterns.items():
        pattern = info.get('pattern', 'unknown')
        if info.get('is_captured'):
            captured_patterns[pattern] += 1
        else:
            missed_patterns[pattern] += 1

    pattern_order = ['sharp', 'trailing', 'bimodal', 'dispersed', 'sparse']
    colors_captured = '#2ecc71'  # Green
    colors_missed = '#e74c3c'    # Red

    fig_width = 12
    fig, ax = plt.subplots(figsize=(fig_width, fig_width / GOLDEN_RATIO))

    x = np.arange(len(pattern_order))
    width = 0.35

    captured_counts = [captured_patterns.get(p, 0) for p in pattern_order]
    missed_counts = [missed_patterns.get(p, 0) for p in pattern_order]

    bars1 = ax.bar(x - width/2, captured_counts, width, label='Captured', color=colors_captured, alpha=0.8)
    bars2 = ax.bar(x + width/2, missed_counts, width, label='Missed', color=colors_missed, alpha=0.8)

    # Add count labels
    for bar, count in zip(bars1, captured_counts):
        if count > 0:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    str(count), ha='center', va='bottom', fontsize=9)
    for bar, count in zip(bars2, missed_counts):
        if count > 0:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    str(count), ha='center', va='bottom', fontsize=9)

    ax.set_xlabel('Truncation Pattern', fontsize=12)
    ax.set_ylabel('Number of Peaks', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels([p.capitalize() for p in pattern_order])
    ax.legend()

    # Add descriptions
    descriptions = {
        'Sharp': 'Reads cluster at peak',
        'Trailing': 'RT truncation gradient',
        'Bimodal': 'Two clusters (alt TSS?)',
        'Dispersed': 'Spread out reads',
        'Sparse': '<5 reads',
    }
    desc_text = '\n'.join([f"{k}: {v}" for k, v in descriptions.items()])
    ax.text(0.98, 0.98, desc_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=PLOT_DPI, bbox_inches='tight')
        plt.close(fig)
        return True
    except Exception as e:
        logger.error(f"Failed to save truncation patterns plot: {e}")
        plt.close(fig)
        return False


def plot_sequence_logo(
    pfm: List[Dict[str, float]],
    output_path: Path,
    title: str,
    upstream: int = 50,
) -> bool:
    """
    Generate sequence logo using ggseqlogo R package.
    
    Requires ggseqlogo to be installed in R:
      install.packages('ggseqlogo')
    """
    if not pfm:
        return False

    import subprocess
    import json
    import tempfile
    import shutil

    # Require Rscript to be available in PATH from the activated environment
    if shutil.which('Rscript') is None:
        logger.error('Rscript not found in PATH; please activate an R-enabled conda env')
        raise RuntimeError('Rscript not found')

    # Prepare JSON input and R script (no auto-install, fail if ggseqlogo missing)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pfm_data = {
        'pfm': pfm,
        'upstream': upstream,
        'positions': list(range(-upstream, len(pfm) - upstream))
    }

    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(pfm_data, f)
        pfm_file = f.name

    # Escape double quotes in title for R string (we'll use double quotes in R)
    r_title = title.replace('"', '\\"')
    
    r_script = f"""
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

if (!require("ggseqlogo", quietly = TRUE)) {{
    stop("ggseqlogo R package not installed; please install in the R environment")
}}

library(jsonlite)
library(ggplot2)

# Read PFM data
pfm_data <- fromJSON('{pfm_file}', simplifyVector = FALSE)
pfm <- pfm_data$pfm
positions <- pfm_data$positions
upstream <- pfm_data$upstream

# Convert PFM to matrix format (A, C, G, T)
if (is.data.frame(pfm)) {{
    pfm_matrix <- as.matrix(pfm)
    # pfm_matrix currently has positions as rows and bases as columns; transpose
    pfm_matrix <- t(pfm_matrix)
    colnames(pfm_matrix) <- positions
}} else {{
    pfm_matrix <- do.call(rbind, lapply(pfm, function(x) {{
        ca <- if (is.null(x[['A']])) 0 else as.numeric(x[['A']])
        cc <- if (is.null(x[['C']])) 0 else as.numeric(x[['C']])
        cg <- if (is.null(x[['G']])) 0 else as.numeric(x[['G']])
        ct <- if (is.null(x[['T']])) 0 else as.numeric(x[['T']])
        c(A = ca, C = cc, G = cg, T = ct)
    }}))
    # currently rows are positions; transpose to have rows=A,C,G,T and columns positions
    pfm_matrix <- t(pfm_matrix)
    colnames(pfm_matrix) <- positions
}}

# Create sequence logo using ggseqlogo (specify dna PFM and bits method)
p <- ggseqlogo(pfm_matrix, seq_type = 'dna', method = 'bits', col_scheme = 'nucleotide') +
    geom_vline(xintercept = 50, linetype = 'dotted', color = 'red', linewidth = 1) +
    labs(title = "{r_title}",
         x = "Position relative to transcript end (bp)",
         y = "Bits") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave('{output_path}', p, width = 16, height = 4, dpi = 150)
cat("Saved sequence logo to {output_path}\\n")
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as f:
        f.write(r_script)
        r_file = f.name

    # Run R script and fail loudly if it errors
    result = subprocess.run(['Rscript', r_file], capture_output=True, text=True, timeout=120)
    if result.returncode != 0:
        logger.error(f"ggseqlogo R script failed: {result.stderr}")
        raise RuntimeError(f"ggseqlogo R script failed: {result.stderr}")

    logger.info(f"Saved sequence logo to {output_path}")
    return True


def plot_transcript_classification(
    classification_counts: Dict[str, int],
    output_path: Path,
    title: str = "Transcript Structural Classification",
) -> bool:
    """Plot bar chart of transcript structural classification categories.

    Categories: FSM, ISM, NIC, NNC, SEM, SEN
    (Full Splice Match, Incomplete Splice Match, Novel In Catalog,
     Novel Not in Catalog, Single-Exon Match, Single-Exon Novel)

    Args:
        classification_counts: Dict mapping category name to count
        output_path: Path to save the plot
        title: Plot title

    Returns:
        True if plot was created successfully, False otherwise
    """
    if not HAS_MATPLOTLIB:
        return False

    if not classification_counts or sum(classification_counts.values()) == 0:
        return False

    # Canonical category order and colors
    category_order = ['FSM', 'ISM', 'NIC', 'NNC', 'SEM', 'SEN']
    category_colors = {
        'FSM': '#2ecc71',   # green - full splice match
        'ISM': '#3498db',   # blue - incomplete splice match
        'NIC': '#f39c12',   # yellow/orange - novel in catalog
        'NNC': '#e74c3c',   # red - novel not in catalog
        'SEM': '#9b59b6',   # purple - single-exon match
        'SEN': '#95a5a6',   # gray - single-exon novel
    }

    categories = []
    counts = []
    colors = []
    for cat in category_order:
        count = classification_counts.get(cat, 0)
        if count > 0:
            categories.append(cat)
            counts.append(count)
            colors.append(category_colors.get(cat, '#7f8c8d'))

    if not categories:
        return False

    fig_width = 10
    fig, ax = plt.subplots(figsize=(fig_width, fig_width / GOLDEN_RATIO))

    bars = ax.bar(categories, counts, color=colors, edgecolor='black', alpha=0.8)

    # Add count labels on bars
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                str(count), ha='center', va='bottom', fontsize=10)

    ax.set_xlabel('Classification', fontsize=12)
    ax.set_ylabel('Number of Transcripts', fontsize=12)
    ax.set_title(title, fontsize=14)

    total = sum(counts)
    pct_text = '\n'.join(
        f"{cat}: {cnt} ({100 * cnt / total:.1f}%)" for cat, cnt in zip(categories, counts)
    )
    ax.text(0.98, 0.98, f'Total: {total}\n{pct_text}',
            transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax.set_axisbelow(True)

    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=PLOT_DPI, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved transcript classification plot to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to save transcript classification plot: {e}")
        plt.close(fig)
        return False


def plot_splice_junction_support(
    supported_sjc: int,
    subset_sjc: int,
    unsupported_sjc: int,
    supported_se: int,
    unsupported_se: int,
    output_path: Path,
    title: str = "Splice Junction Support",
) -> bool:
    """Plot stacked bar chart of splice junction support categories.

    Shows multi-exon (supported / subset / unsupported) and single-exon
    (supported / unsupported) splice junction chain counts.

    Args:
        supported_sjc: Fully supported splice junction chains
        subset_sjc: Subset-supported splice junction chains
        unsupported_sjc: Unsupported splice junction chains
        supported_se: Supported single-exon transcripts
        unsupported_se: Unsupported single-exon transcripts
        output_path: Path to save the plot
        title: Plot title

    Returns:
        True if plot was created successfully, False otherwise
    """
    if not HAS_MATPLOTLIB:
        return False

    total_sjc = supported_sjc + subset_sjc + unsupported_sjc
    total_se = supported_se + unsupported_se
    if total_sjc == 0 and total_se == 0:
        return False

    fig_width = 8
    fig, ax = plt.subplots(figsize=(fig_width, fig_width / GOLDEN_RATIO))

    bar_labels = []
    bar_data = []  # list of (bottom_vals, heights, colors, labels)

    if total_sjc > 0:
        bar_labels.append('Multi-Exon')
        bar_data.append([
            (0, supported_sjc, '#2ecc71', 'Supported'),
            (supported_sjc, subset_sjc, '#f39c12', 'Subset'),
            (supported_sjc + subset_sjc, unsupported_sjc, '#e74c3c', 'Unsupported'),
        ])
    if total_se > 0:
        bar_labels.append('Single-Exon')
        bar_data.append([
            (0, supported_se, '#2ecc71', 'Supported'),
            (supported_se, unsupported_se, '#e74c3c', 'Unsupported'),
        ])

    x = np.arange(len(bar_labels))
    width = 0.5

    # Track legend entries (avoid duplicates)
    legend_entries = {}
    for i, segments in enumerate(bar_data):
        for bottom, height, color, label in segments:
            if height == 0:
                continue
            bar = ax.bar(x[i], height, width, bottom=bottom, color=color, edgecolor='black', alpha=0.8)
            if label not in legend_entries:
                legend_entries[label] = bar[0]
            # Add count label
            if height > 0:
                ax.text(x[i], bottom + height / 2, str(height),
                        ha='center', va='center', fontsize=10, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(bar_labels, fontsize=11)
    ax.set_ylabel('Number of Read Splice Junction Chains', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.legend(legend_entries.values(), legend_entries.keys(), loc='upper right')

    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax.set_axisbelow(True)

    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=PLOT_DPI, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved splice junction support plot to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to save splice junction support plot: {e}")
        plt.close(fig)
        return False


def plot_missed_peak_sj_support(
    sj_support_summary: Dict[str, int],
    output_path: Path,
    title: str = "Splice Junction Support at Missed Recoverable Peaks",
) -> bool:
    """Plot bar chart showing SJ support levels for missed recoverable peaks.

    Shows how many missed peaks have reads with full splice junction match,
    subset match, unsupported junctions, or single-exon reads as their best
    SJ evidence.

    Args:
        sj_support_summary: Dict from analyze_missed_peaks_comprehensive with
            keys like 'full_match', 'subset_match', 'unsupported', 'single_exon'
            (peak-level best SJ support) and 'reads_full_match', etc. (read counts)
        output_path: Path to save the plot
        title: Plot title

    Returns:
        True if plot was created successfully, False otherwise
    """
    if not HAS_MATPLOTLIB:
        return False

    # Peak-level best SJ support counts
    categories = ['full_match', 'subset_match', 'unsupported', 'single_exon']
    labels = ['Full Match', 'Subset Match', 'Unsupported', 'Single-Exon']
    colors = ['#2ecc71', '#f39c12', '#e74c3c', '#95a5a6']

    peak_counts = [sj_support_summary.get(cat, 0) for cat in categories]
    total_peaks = sum(peak_counts)
    if total_peaks == 0:
        return False

    fig_width = 10
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fig_width, fig_width / GOLDEN_RATIO))

    # Left panel: peak-level best SJ support
    present = [(l, c, col) for l, c, col in zip(labels, peak_counts, colors) if c > 0]
    if present:
        bar_labels, bar_counts, bar_colors = zip(*present)
        bars = ax1.bar(bar_labels, bar_counts, color=bar_colors, edgecolor='black', alpha=0.8)
        for bar, cnt in zip(bars, bar_counts):
            ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                     str(cnt), ha='center', va='bottom', fontsize=10)
    ax1.set_xlabel('Best SJ Support Level', fontsize=11)
    ax1.set_ylabel('Number of Missed Peaks', fontsize=11)
    ax1.set_title('Best Read SJ Support per Peak', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax1.set_axisbelow(True)
    ax1.tick_params(axis='x', rotation=20)

    # Right panel: read-level SJ support counts
    read_counts = [sj_support_summary.get(f'reads_{cat}', 0) for cat in categories]
    total_reads = sum(read_counts)
    present_r = [(l, c, col) for l, c, col in zip(labels, read_counts, colors) if c > 0]
    if present_r:
        bar_labels_r, bar_counts_r, bar_colors_r = zip(*present_r)
        bars_r = ax2.bar(bar_labels_r, bar_counts_r, color=bar_colors_r, edgecolor='black', alpha=0.8)
        for bar, cnt in zip(bars_r, bar_counts_r):
            ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                     str(cnt), ha='center', va='bottom', fontsize=10)
    ax2.set_xlabel('SJ Support Level', fontsize=11)
    ax2.set_ylabel('Number of Supporting Reads', fontsize=11)
    ax2.set_title('All Supporting Reads SJ Support', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax2.set_axisbelow(True)
    ax2.tick_params(axis='x', rotation=20)

    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=PLOT_DPI, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved missed peak SJ support plot to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to save missed peak SJ support plot: {e}")
        plt.close(fig)
        return False


def plot_peak_recovery_by_expression(
    peak_metadata: List[dict],
    output_path: Path,
    title: str = "Peak Recovery Rate by Expression Level",
) -> bool:
    """Plot peak recovery rate as a function of TPM expression.

    Two-panel plot:
    - Left: recovery rate per TPM bin with bar chart and count annotations
    - Right: cumulative recovery curve (peaks sorted by TPM descending)

    Args:
        peak_metadata: list of dicts with keys 'score' (TPM), 'status' ('captured'|'missed'|'no_reads')
        output_path: where to save the plot
        title: figure title
    """
    if not HAS_MATPLOTLIB or not peak_metadata:
        return False

    scores = [p['score'] for p in peak_metadata]
    statuses = [p['status'] for p in peak_metadata]

    if not scores or max(scores) == 0:
        logger.warning("No TPM data available for recovery-by-expression plot")
        return False

    # Define TPM bins: 0-10, 10-25, 25-50, 50-100, 100-250, 250-500, 500+
    bin_edges = [0, 10, 25, 50, 100, 250, 500, float('inf')]
    bin_labels = ['0-10', '10-25', '25-50', '50-100', '100-250', '250-500', '500+']

    bin_captured = [0] * len(bin_labels)
    bin_total = [0] * len(bin_labels)

    for score, status in zip(scores, statuses):
        for i in range(len(bin_edges) - 1):
            if bin_edges[i] <= score < bin_edges[i + 1]:
                bin_total[i] += 1
                if status == 'captured':
                    bin_captured[i] += 1
                break

    bin_rates = []
    for c, t in zip(bin_captured, bin_total):
        bin_rates.append(c / t if t > 0 else 0)

    fig_width = 10
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fig_width, fig_width / GOLDEN_RATIO))

    # Left panel: recovery rate per TPM bin
    non_empty = [(lbl, rate, cap, tot) for lbl, rate, cap, tot in
                 zip(bin_labels, bin_rates, bin_captured, bin_total) if tot > 0]
    if non_empty:
        x_labels, x_rates, x_cap, x_tot = zip(*non_empty)
        colors = [plt.cm.RdYlGn(r) for r in x_rates]
        bars = ax1.bar(range(len(x_labels)), x_rates, color=colors, edgecolor='black', alpha=0.85)
        for i, (bar, cap, tot) in enumerate(zip(bars, x_cap, x_tot)):
            ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                     f"{cap}/{tot}", ha='center', va='bottom', fontsize=9)
        ax1.set_xticks(range(len(x_labels)))
        ax1.set_xticklabels(x_labels, rotation=30, ha='right')
    ax1.set_xlabel('Peak Expression (TPM)', fontsize=11)
    ax1.set_ylabel('Recovery Rate', fontsize=11)
    ax1.set_title('Recovery Rate by TPM Bin', fontsize=12, fontweight='bold')
    ax1.set_ylim(0, 1.15)
    ax1.axhline(y=1.0, color='gray', linestyle='--', alpha=0.4)
    ax1.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax1.set_axisbelow(True)

    # Right panel: cumulative recovery curve sorted by TPM descending
    paired = sorted(zip(scores, statuses), key=lambda x: -x[0])
    cum_total = 0
    cum_captured = 0
    cum_totals = []
    cum_rates = []
    tpm_thresholds = []
    for score, status in paired:
        cum_total += 1
        if status == 'captured':
            cum_captured += 1
        cum_totals.append(cum_total)
        cum_rates.append(cum_captured / cum_total)
        tpm_thresholds.append(score)

    ax2.plot(range(len(cum_rates)), cum_rates, color='#2166AC', linewidth=1.5)
    ax2.fill_between(range(len(cum_rates)), cum_rates, alpha=0.15, color='#2166AC')

    # Add TPM threshold markers
    marker_tpms = [10, 25, 50, 100, 250, 500]
    for mtpm in marker_tpms:
        idx = None
        for j, t in enumerate(tpm_thresholds):
            if t <= mtpm:
                idx = j
                break
        if idx is not None and idx < len(cum_rates):
            ax2.axvline(x=idx, color='gray', linestyle=':', alpha=0.4)
            ax2.text(idx, 0.02, f'TPM\u2264{mtpm}', fontsize=7, rotation=90,
                     va='bottom', ha='right', alpha=0.6)

    ax2.set_xlabel(f'Peaks (sorted by TPM, n={len(cum_rates)})', fontsize=11)
    ax2.set_ylabel('Cumulative Recovery Rate', fontsize=11)
    ax2.set_title('Cumulative Recovery (High TPM First)', fontsize=12, fontweight='bold')
    ax2.set_ylim(0, 1.05)
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.set_axisbelow(True)

    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=PLOT_DPI, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved peak recovery by expression plot to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to save peak recovery by expression plot: {e}")
        plt.close(fig)
        return False


def plot_peak_recovery_by_width(
    peak_metadata: List[dict],
    output_path: Path,
    title: str = "Peak Recovery Rate by Peak Width",
) -> bool:
    """Plot peak recovery rate as a function of peak width (CAGE peaks only).

    Two-panel plot:
    - Left: recovery rate per width bin with bar chart
    - Right: scatter of peak width vs TPM, colored by recovery status

    Args:
        peak_metadata: list of dicts with keys 'width', 'score' (TPM), 'status'
        output_path: where to save the plot
        title: figure title
    """
    if not HAS_MATPLOTLIB or not peak_metadata:
        return False

    widths = [p['width'] for p in peak_metadata]
    if not widths or max(widths) <= 1:
        logger.warning("No variable-width peaks for recovery-by-width plot")
        return False

    scores = [p['score'] for p in peak_metadata]
    statuses = [p['status'] for p in peak_metadata]

    # Define width bins: 1, 2-5, 6-10, 11-25, 26-50, 51-100, 100+
    bin_edges = [0, 2, 6, 11, 26, 51, 101, float('inf')]
    bin_labels = ['1', '2-5', '6-10', '11-25', '26-50', '51-100', '100+']

    bin_captured = [0] * len(bin_labels)
    bin_total = [0] * len(bin_labels)

    for width, status in zip(widths, statuses):
        for i in range(len(bin_edges) - 1):
            if bin_edges[i] <= width < bin_edges[i + 1]:
                bin_total[i] += 1
                if status == 'captured':
                    bin_captured[i] += 1
                break

    bin_rates = []
    for c, t in zip(bin_captured, bin_total):
        bin_rates.append(c / t if t > 0 else 0)

    fig_width = 10
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fig_width, fig_width / GOLDEN_RATIO))

    # Left panel: recovery rate per width bin
    non_empty = [(lbl, rate, cap, tot) for lbl, rate, cap, tot in
                 zip(bin_labels, bin_rates, bin_captured, bin_total) if tot > 0]
    if non_empty:
        x_labels, x_rates, x_cap, x_tot = zip(*non_empty)
        colors = [plt.cm.RdYlGn(r) for r in x_rates]
        bars = ax1.bar(range(len(x_labels)), x_rates, color=colors, edgecolor='black', alpha=0.85)
        for i, (bar, cap, tot) in enumerate(zip(bars, x_cap, x_tot)):
            ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                     f"{cap}/{tot}", ha='center', va='bottom', fontsize=9)
        ax1.set_xticks(range(len(x_labels)))
        ax1.set_xticklabels(x_labels, rotation=30, ha='right')
    ax1.set_xlabel('Peak Width (bp)', fontsize=11)
    ax1.set_ylabel('Recovery Rate', fontsize=11)
    ax1.set_title('Recovery Rate by Peak Width', fontsize=12, fontweight='bold')
    ax1.set_ylim(0, 1.15)
    ax1.axhline(y=1.0, color='gray', linestyle='--', alpha=0.4)
    ax1.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax1.set_axisbelow(True)

    # Right panel: scatter of width vs TPM colored by status
    captured_w = [w for w, s in zip(widths, statuses) if s == 'captured']
    captured_s = [sc for sc, s in zip(scores, statuses) if s == 'captured']
    missed_w = [w for w, s in zip(widths, statuses) if s == 'missed']
    missed_s = [sc for sc, s in zip(scores, statuses) if s == 'missed']
    noread_w = [w for w, s in zip(widths, statuses) if s == 'no_reads']
    noread_s = [sc for sc, s in zip(scores, statuses) if s == 'no_reads']

    if noread_w:
        ax2.scatter(noread_w, noread_s, c='#BDBDBD', alpha=0.4, s=15, label=f'No reads ({len(noread_w)})', zorder=1)
    if missed_w:
        ax2.scatter(missed_w, missed_s, c='#D32F2F', alpha=0.5, s=20, label=f'Missed ({len(missed_w)})', zorder=2)
    if captured_w:
        ax2.scatter(captured_w, captured_s, c='#388E3C', alpha=0.5, s=20, label=f'Captured ({len(captured_w)})', zorder=3)

    ax2.set_xlabel('Peak Width (bp)', fontsize=11)
    ax2.set_ylabel('Peak Expression (TPM)', fontsize=11)
    ax2.set_title('Width vs Expression by Recovery', fontsize=12, fontweight='bold')
    if scores and max(scores) > 100:
        ax2.set_yscale('log')
    ax2.legend(fontsize=9, loc='lower right', framealpha=0.9)
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.set_axisbelow(True)

    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=PLOT_DPI, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved peak recovery by width plot to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to save peak recovery by width plot: {e}")
        plt.close(fig)
        return False


def plot_read_end_frequency_at_peaks(
    read_positions: List[dict],
    peaks: List[dict],
    output_path: Path,
    window: int = 100,
    title: str = "Read End Frequency Around Peaks",
    end_type: str = "tts",
) -> bool:
    """Plot read end frequency distribution relative to peak positions.

    For 1bp peaks (QuantSeq): single histogram of read end offsets from peak center.
    Shows where reads terminate relative to the orthogonal peak signal.

    Args:
        read_positions: List of dicts with 'Chrom', 'Start', 'End', 'Strand'
        peaks: List of dicts with 'Chrom', 'Start', 'End', 'Strand', 'Score'
        output_path: Where to save the plot
        window: bp upstream/downstream to consider
        title: Figure title
        end_type: 'tss' or 'tts' for labeling
    """
    if not HAS_MATPLOTLIB or not read_positions or not peaks:
        return False

    # Build lookup by chrom/strand for efficient matching
    from collections import defaultdict
    reads_by_cs = defaultdict(list)
    for r in read_positions:
        key = (r['Chrom'], r['Strand'])
        pos = (r['Start'] + r['End']) // 2  # midpoint of 1bp interval
        reads_by_cs[key].append(pos)

    # Compute offsets from peak centers
    offsets = []
    for peak in peaks:
        peak_center = (peak['Start'] + peak['End']) // 2
        strands = [peak['Strand']] if peak['Strand'] != '.' else ['+', '-']
        for strand in strands:
            key = (peak['Chrom'], strand)
            for read_pos in reads_by_cs.get(key, []):
                offset = read_pos - peak_center
                # Flip offset for minus strand so upstream is always negative
                if strand == '-':
                    offset = -offset
                if -window <= offset <= window:
                    offsets.append(offset)

    if not offsets:
        logger.warning("No read ends found within window of peaks")
        return False

    fig_width = 10
    fig, ax = plt.subplots(figsize=(fig_width, fig_width / GOLDEN_RATIO))

    bins = np.arange(-window, window + 2, 2)  # 2bp resolution
    ax.hist(offsets, bins=bins, color='#2196F3', edgecolor='black', alpha=0.7)
    ax.axvline(x=0, color='red', linestyle='--', linewidth=2, label='Peak center')
    ax.set_xlabel(f'Distance from peak center (bp)\n← Upstream | Downstream →', fontsize=11)
    ax.set_ylabel('Read end count', fontsize=11)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax.set_axisbelow(True)

    # Add summary stats
    median_offset = np.median(offsets)
    ax.text(0.98, 0.95, f'n={len(offsets):,}\nmedian={median_offset:.1f}bp',
            transform=ax.transAxes, ha='right', va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=PLOT_DPI, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved read end frequency plot to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to save read end frequency plot: {e}")
        plt.close(fig)
        return False


def plot_read_end_frequency_stratified_by_width(
    read_positions: List[dict],
    peaks: List[dict],
    output_path: Path,
    window: int = 100,
    title: str = "Read End Frequency Around CAGE Peaks by Width",
) -> bool:
    """Plot read end frequency stratified by peak width quartiles.

    For variable-width peaks (CAGE): 2x2 grid showing read end distributions
    for narrow, medium, wide, and very wide peaks separately.

    Args:
        read_positions: List of dicts with 'Chrom', 'Start', 'End', 'Strand'
        peaks: List of dicts with 'Chrom', 'Start', 'End', 'Strand', 'Score'
        output_path: Where to save the plot
        window: bp upstream/downstream to consider (scaled by peak width for wider peaks)
        title: Figure title
    """
    if not HAS_MATPLOTLIB or not read_positions or not peaks:
        return False

    # Calculate peak widths and stratify
    widths = [p['End'] - p['Start'] for p in peaks]
    if max(widths) <= 1:
        logger.warning("All peaks are 1bp, use plot_read_end_frequency_at_peaks instead")
        return False

    # Define width bins: 1-10bp, 11-25bp, 26-50bp, 51+bp
    width_bins = [
        ('1-10bp', 1, 10),
        ('11-25bp', 11, 25),
        ('26-50bp', 26, 50),
        ('51+bp', 51, float('inf')),
    ]

    # Build read lookup
    from collections import defaultdict
    reads_by_cs = defaultdict(list)
    for r in read_positions:
        key = (r['Chrom'], r['Strand'])
        pos = (r['Start'] + r['End']) // 2
        reads_by_cs[key].append(pos)

    # Compute offsets for each width bin
    offsets_by_bin = {label: [] for label, _, _ in width_bins}
    peak_counts = {label: 0 for label, _, _ in width_bins}

    for peak in peaks:
        peak_width = peak['End'] - peak['Start']
        peak_center = (peak['Start'] + peak['End']) // 2

        # Find which bin this peak belongs to
        bin_label = None
        for label, lo, hi in width_bins:
            if lo <= peak_width <= hi:
                bin_label = label
                break
        if bin_label is None:
            continue

        peak_counts[bin_label] += 1
        strands = [peak['Strand']] if peak['Strand'] != '.' else ['+', '-']

        # Use window scaled to peak width (at least 'window', at most 2x window)
        effective_window = max(window, min(peak_width * 2, window * 2))

        for strand in strands:
            key = (peak['Chrom'], strand)
            for read_pos in reads_by_cs.get(key, []):
                offset = read_pos - peak_center
                if strand == '-':
                    offset = -offset
                if -effective_window <= offset <= effective_window:
                    # Normalize offset relative to peak half-width for comparison
                    half_width = max(peak_width // 2, 1)
                    normalized_offset = offset / half_width
                    offsets_by_bin[bin_label].append((offset, normalized_offset))

    fig_width = 12
    fig, axes = plt.subplots(2, 2, figsize=(fig_width, fig_width / GOLDEN_RATIO * 1.2))
    axes = axes.flatten()

    colors = ['#4CAF50', '#2196F3', '#FF9800', '#9C27B0']

    for idx, (label, lo, hi) in enumerate(width_bins):
        ax = axes[idx]
        data = offsets_by_bin[label]
        n_peaks = peak_counts[label]

        if data:
            raw_offsets = [d[0] for d in data]
            # Use appropriate binning based on peak width range
            if hi <= 10:
                bins = np.arange(-window, window + 2, 2)
            elif hi <= 25:
                bins = np.arange(-window * 1.5, window * 1.5 + 3, 3)
            else:
                bins = np.arange(-window * 2, window * 2 + 5, 5)

            ax.hist(raw_offsets, bins=bins, color=colors[idx], edgecolor='black', alpha=0.7)
            median_off = np.median(raw_offsets)
            ax.text(0.98, 0.95, f'n={len(raw_offsets):,}\npeaks={n_peaks}\nmed={median_off:.0f}bp',
                    transform=ax.transAxes, ha='right', va='top', fontsize=9,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        else:
            ax.text(0.5, 0.5, f'No data\n({n_peaks} peaks)', transform=ax.transAxes,
                    ha='center', va='center', fontsize=12)

        ax.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
        ax.set_title(f'Peak width: {label}', fontsize=11, fontweight='bold')
        ax.set_xlabel('Distance from peak center (bp)', fontsize=10)
        ax.set_ylabel('Read end count', fontsize=10)
        ax.grid(True, alpha=0.3, linestyle='--', axis='y')
        ax.set_axisbelow(True)

    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=PLOT_DPI, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved stratified read end frequency plot to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to save stratified read end frequency plot: {e}")
        plt.close(fig)
        return False


def plot_peak_width_histogram(
    peak_metadata: List[dict],
    output_path: Path,
    title: str = "CAGE Peak Width Distribution by Recovery Status",
) -> bool:
    """Plot histogram of peak widths with stacked bars colored by recovery status.

    Args:
        peak_metadata: List of dicts with 'width' and 'status' ('captured'/'missed'/'no_reads')
        output_path: Where to save the plot
        title: Figure title
    """
    if not HAS_MATPLOTLIB or not peak_metadata:
        return False

    widths = [p['width'] for p in peak_metadata]
    if max(widths) <= 1:
        logger.warning("All peaks are 1bp, width histogram not meaningful")
        return False

    # Separate by status
    captured_widths = [p['width'] for p in peak_metadata if p['status'] == 'captured']
    missed_widths = [p['width'] for p in peak_metadata if p['status'] == 'missed']
    no_reads_widths = [p['width'] for p in peak_metadata if p['status'] == 'no_reads']

    fig_width = 10
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fig_width, fig_width / GOLDEN_RATIO))

    # Left panel: stacked histogram
    max_width = max(widths)
    if max_width <= 50:
        bins = np.arange(0, max_width + 2, 1)
    elif max_width <= 200:
        bins = np.arange(0, max_width + 5, 5)
    else:
        bins = np.arange(0, min(max_width + 10, 500), 10)

    ax1.hist([captured_widths, missed_widths, no_reads_widths], bins=bins, stacked=True,
             color=['#388E3C', '#D32F2F', '#BDBDBD'], edgecolor='black', alpha=0.85,
             label=[f'Captured ({len(captured_widths)})',
                    f'Missed ({len(missed_widths)})',
                    f'No reads ({len(no_reads_widths)})'])
    ax1.set_xlabel('Peak Width (bp)', fontsize=11)
    ax1.set_ylabel('Number of Peaks', fontsize=11)
    ax1.set_title('Peak Count by Width', fontsize=12, fontweight='bold')
    ax1.legend(fontsize=9, loc='upper right')
    ax1.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax1.set_axisbelow(True)

    # Right panel: recovery rate by width bin
    bin_edges = [0, 5, 10, 20, 30, 50, 100, float('inf')]
    bin_labels = ['1-5', '6-10', '11-20', '21-30', '31-50', '51-100', '100+']
    bin_captured = [0] * len(bin_labels)
    bin_total = [0] * len(bin_labels)

    for p in peak_metadata:
        w = p['width']
        for i in range(len(bin_edges) - 1):
            if bin_edges[i] < w <= bin_edges[i + 1]:
                bin_total[i] += 1
                if p['status'] == 'captured':
                    bin_captured[i] += 1
                break

    bin_rates = [c / t if t > 0 else 0 for c, t in zip(bin_captured, bin_total)]
    non_empty = [(lbl, rate, cap, tot) for lbl, rate, cap, tot in
                 zip(bin_labels, bin_rates, bin_captured, bin_total) if tot > 0]

    if non_empty:
        x_labels, x_rates, x_cap, x_tot = zip(*non_empty)
        colors = [plt.cm.RdYlGn(r) for r in x_rates]
        bars = ax2.bar(range(len(x_labels)), x_rates, color=colors, edgecolor='black', alpha=0.85)
        for i, (bar, cap, tot) in enumerate(zip(bars, x_cap, x_tot)):
            ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                     f"{cap}/{tot}", ha='center', va='bottom', fontsize=9)
        ax2.set_xticks(range(len(x_labels)))
        ax2.set_xticklabels(x_labels, rotation=30, ha='right')

    ax2.set_xlabel('Peak Width (bp)', fontsize=11)
    ax2.set_ylabel('Recovery Rate', fontsize=11)
    ax2.set_title('Recovery Rate by Width', fontsize=12, fontweight='bold')
    ax2.set_ylim(0, 1.15)
    ax2.axhline(y=1.0, color='gray', linestyle='--', alpha=0.4)
    ax2.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax2.set_axisbelow(True)

    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=PLOT_DPI, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved peak width histogram to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to save peak width histogram: {e}")
        plt.close(fig)
        return False
