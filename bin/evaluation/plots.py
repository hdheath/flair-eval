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

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

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
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
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

    fig, ax = plt.subplots(figsize=(10, 6))
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
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
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
            title=f"Read 5' End Offset from Isoform Model TSS\n{plot_prefix}",
        )

    # Plot 1b: Aggregate TTS offset histogram
    tts_offsets = entropy_data["all_tts_offsets"]
    if tts_offsets:
        plot_distance_histogram(
            distances=tts_offsets,
            output_path=plot_output_dir / f"{plot_prefix}_read_tts_offset_histogram.png",
            title=f"Read 3' End Offset from Isoform Model TTS\n{plot_prefix}",
        )

    # Plot 2a: Per-isoform TSS entropy distribution
    tss_entropies = [e for _, e, _ in entropy_data["tss_entropy_per_isoform"]]
    if tss_entropies:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(tss_entropies, bins=50, edgecolor='black', alpha=0.7, color='coral')
        ax.set_xlabel("Shannon Entropy (bits, 10bp bins)", fontsize=12)
        ax.set_ylabel("Number of Isoforms", fontsize=12)
        ax.set_title(f"Per-Isoform TSS Read-End Entropy (5' End)\n{plot_prefix}", fontsize=14)
        mean_e = statistics.mean(tss_entropies)
        median_e = statistics.median(tss_entropies)
        stats_text = f'n = {len(tss_entropies)}\nMean = {mean_e:.2f} bits\nMedian = {median_e:.2f} bits'
        ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.tight_layout()
        plt.savefig(plot_output_dir / f"{plot_prefix}_tss_entropy_distribution.png",
                    dpi=150, bbox_inches='tight')
        plt.close(fig)

    # Plot 2b: Per-isoform TTS entropy distribution
    tts_entropies = [e for _, e, _ in entropy_data["tts_entropy_per_isoform"]]
    if tts_entropies:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(tts_entropies, bins=50, edgecolor='black', alpha=0.7, color='mediumpurple')
        ax.set_xlabel("Shannon Entropy (bits, 10bp bins)", fontsize=12)
        ax.set_ylabel("Number of Isoforms", fontsize=12)
        ax.set_title(f"Per-Isoform TTS Read-End Entropy (3' End)\n{plot_prefix}", fontsize=14)
        mean_e = statistics.mean(tts_entropies)
        median_e = statistics.median(tts_entropies)
        stats_text = f'n = {len(tts_entropies)}\nMean = {mean_e:.2f} bits\nMedian = {median_e:.2f} bits'
        ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.tight_layout()
        plt.savefig(plot_output_dir / f"{plot_prefix}_tts_entropy_distribution.png",
                    dpi=150, bbox_inches='tight')
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

    fig, ax = plt.subplots(figsize=(12, 6))  # wider helps readability

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

    # X positions: 0..100, plus overflow at 101 labeled "100+"
    x_main = list(range(0, max_bin + 1))
    y_main = freq

    overflow_x = max_bin + 1  # 101
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
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
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
    """Plot pie chart of read classification categories for missed peaks."""
    if not HAS_MATPLOTLIB:
        return False

    if not classification_summary or sum(classification_summary.values()) == 0:
        return False

    fig, ax = plt.subplots(figsize=(8, 8))

    labels = []
    sizes = []
    colors = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99']
    category_labels = {
        'unassigned': 'Unassigned',
        'assigned_nearby': 'Assigned Nearby',
        'assigned_distant': 'Assigned Distant',
        'assigned_wrong_strand': 'Wrong Strand',
    }

    for cat in ['unassigned', 'assigned_nearby', 'assigned_distant', 'assigned_wrong_strand']:
        count = classification_summary.get(cat, 0)
        if count > 0:
            labels.append(f"{category_labels[cat]}\n({count})")
            sizes.append(count)

    if not sizes:
        plt.close(fig)
        return False

    ax.pie(sizes, labels=labels, colors=colors[:len(sizes)], autopct='%1.1f%%',
           startangle=90, textprops={'fontsize': 10})
    ax.set_title(title, fontsize=14)

    plt.tight_layout()
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return True
    except Exception as e:
        logger.error(f"Failed to save classification pie chart: {e}")
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

    fig, ax = plt.subplots(figsize=(10, 6))

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
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
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

    fig, ax = plt.subplots(figsize=(12, 6))

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
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
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
