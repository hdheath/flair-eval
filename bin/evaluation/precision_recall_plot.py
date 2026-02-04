#!/usr/bin/env python3
"""
Precision-Recall scatter plots for comparing assemblers.

Creates a 3-panel figure:
  1. CAGE (5') Precision vs Recall
  2. QuantSeq (3') Precision vs Recall
  3. Reference 5' Precision vs Reference 3' Precision

Points are colored by assembler (FLAIR, Bambu, IsoQuant).
"""

import argparse
import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Optional matplotlib import
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# Assembler color scheme
ASSEMBLER_COLORS = {
    'flair': '#1f77b4',      # blue
    'bambu': '#ff7f0e',      # orange
    'isoquant': '#2ca02c',   # green
    'unknown': '#7f7f7f',    # gray
}

# Assembler display names
ASSEMBLER_NAMES = {
    'flair': 'FLAIR',
    'bambu': 'Bambu',
    'isoquant': 'IsoQuant',
    'unknown': 'Unknown',
}


def identify_assembler(row):
    """Identify which assembler produced this result based on transcriptome_mode."""
    mode = str(row.get('transcriptome_mode', '')).lower()

    if 'bambu' in mode:
        return 'bambu'
    elif 'isoquant' in mode:
        return 'isoquant'
    elif mode and mode not in ['', 'nan', 'none']:
        # Assume FLAIR for other modes (default, max-ends, etc.)
        return 'flair'
    return 'unknown'


def load_evaluation_files(input_files):
    """Load and concatenate multiple evaluation TSV files."""
    dfs = []
    for f in input_files:
        try:
            df = pd.read_csv(f, sep='\t')
            df['source_file'] = str(f)
            dfs.append(df)
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}", file=sys.stderr)

    if not dfs:
        return None

    combined = pd.concat(dfs, ignore_index=True)
    combined['assembler'] = combined.apply(identify_assembler, axis=1)
    return combined


def create_precision_recall_plot(df, output_path, title_prefix=""):
    """Create the 3-panel precision-recall plot."""
    if not HAS_MATPLOTLIB:
        print("Warning: matplotlib not available, skipping plot generation", file=sys.stderr)
        return False

    # Filter to rows with valid data
    df = df.copy()

    # Set up the figure with 3 subplots (golden ratio proportions)
    fig_width = 15
    golden_ratio = 1.618
    fig, axes = plt.subplots(1, 3, figsize=(fig_width, fig_width / (3 * golden_ratio / 1.2)))
    fig.suptitle(f'{title_prefix}Precision vs Recall by Assembler', fontsize=14, fontweight='bold')

    # Common plot settings
    plot_configs = [
        {
            'ax': axes[0],
            'x_col': '5prime_recall',
            'y_col': '5prime_precision',
            'title': "CAGE TSS",
            'xlabel': "5' Recall",
            'ylabel': "5' Precision",
        },
        {
            'ax': axes[1],
            'x_col': '3prime_recall',
            'y_col': '3prime_precision',
            'title': "QuantSeq TTS",
            'xlabel': "3' Recall",
            'ylabel': "3' Precision",
        },
        {
            'ax': axes[2],
            'x_col': 'ref3prime_precision',
            'y_col': 'ref5prime_precision',
            'title': "Reference End Precision",
            'xlabel': "Ref 3' Precision",
            'ylabel': "Ref 5' Precision",
        },
    ]

    # Track which assemblers have data for legend
    assemblers_with_data = set()

    for config in plot_configs:
        ax = config['ax']
        x_col = config['x_col']
        y_col = config['y_col']

        # Filter to valid data points for this plot
        valid_mask = df[x_col].notna() & df[y_col].notna()
        plot_df = df[valid_mask]

        # Plot each assembler
        for assembler in ['flair', 'bambu', 'isoquant', 'unknown']:
            assembler_df = plot_df[plot_df['assembler'] == assembler]
            if len(assembler_df) == 0:
                continue

            assemblers_with_data.add(assembler)

            x_vals = assembler_df[x_col].values * 100  # Convert to percentage
            y_vals = assembler_df[y_col].values * 100

            ax.scatter(x_vals, y_vals,
                      c=ASSEMBLER_COLORS[assembler],
                      label=ASSEMBLER_NAMES[assembler],
                      s=80, alpha=0.7, edgecolors='black', linewidths=0.5,
                      zorder=3)

        # Configure axes
        ax.set_xlabel(config['xlabel'], fontsize=11)
        ax.set_ylabel(config['ylabel'], fontsize=11)
        ax.set_title(config['title'], fontsize=12, fontweight='bold')

        # Set axis limits (0-100 for percentages)
        ax.set_xlim(-5, 105)
        ax.set_ylim(-5, 105)

        # Add grid
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)

        # Add diagonal reference line for P/R plots
        if 'recall' in x_col.lower() and 'precision' in y_col.lower():
            ax.plot([0, 100], [0, 100], 'k--', alpha=0.3, linewidth=1, zorder=1)

        # Tick formatting
        ax.set_xticks(range(0, 101, 20))
        ax.set_yticks(range(0, 101, 20))

    # Create legend
    legend_handles = [
        Line2D([0], [0], marker='o', linestyle='',
               markerfacecolor=ASSEMBLER_COLORS[a],
               markeredgecolor='black', markersize=10,
               label=ASSEMBLER_NAMES[a], alpha=0.7)
        for a in ['flair', 'bambu', 'isoquant']
        if a in assemblers_with_data
    ]

    if legend_handles:
        fig.legend(handles=legend_handles,
                  loc='upper center',
                  bbox_to_anchor=(0.5, 0.02),
                  ncol=len(legend_handles),
                  frameon=True,
                  fontsize=10)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)  # Make room for legend

    # Save figure
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Create precision-recall scatter plots comparing assemblers"
    )
    parser.add_argument(
        '--input', '-i',
        nargs='+',
        required=True,
        help="Input evaluation TSV file(s)"
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help="Output PNG file path"
    )
    parser.add_argument(
        '--title-prefix',
        default="",
        help="Optional prefix for the plot title"
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help="Print verbose output"
    )

    args = parser.parse_args()

    if not HAS_MATPLOTLIB:
        print("Error: matplotlib is required for plotting", file=sys.stderr)
        sys.exit(1)

    # Load data
    if args.verbose:
        print(f"Loading {len(args.input)} evaluation file(s)...")

    df = load_evaluation_files(args.input)

    if df is None or len(df) == 0:
        print("Error: No valid data found in input files", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print(f"Loaded {len(df)} evaluation records")
        print(f"Assemblers: {df['assembler'].value_counts().to_dict()}")

    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Create plot
    success = create_precision_recall_plot(df, output_path, args.title_prefix)

    if success:
        print(f"Saved precision-recall plot to {output_path}")
    else:
        print("Failed to create plot", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
