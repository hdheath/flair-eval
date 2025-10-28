#!/usr/bin/env python3
"""
sqanti_plot.py

Find every sqanti_results.tsv under results/, plot a 3‐panel figure,
and write sqanti.png under plots/region/run_name/, skipping if present.
"""

import os, sys, glob, argparse, logging
import pandas as pd
import matplotlib.pyplot as plt

# SQANTI isoform classes in canonical order
ISOCLASSES = [
    'FSM', 'ISM', 'NIC', 'NNC',
    'Genic_Genomic', 'Genic_Intron',
    'Antisense', 'Fusion', 'Intergenic'
]

# ─── same hack ───────────────────────────────────────────────────────────────
THIS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(THIS, "..", ".."))
SRC  = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

logger = logging.getLogger(__name__)


def plot_summary(summary_tsv_or_df, plot_dir):
    """Plot a summary from a TSV filepath or a pandas.DataFrame.

    summary_tsv_or_df may be either a string path to a TSV (tab-separated)
    or a pandas.DataFrame already in memory (useful for piping results
    from another script).
    """
    out_png = os.path.join(plot_dir, "sqanti.png")
    if os.path.exists(out_png):
        logger.debug(f"[SKIP] plot already exists: {out_png}")
        return

    os.makedirs(plot_dir, exist_ok=True)
    if isinstance(summary_tsv_or_df, pd.DataFrame):
        df = summary_tsv_or_df
    else:
        df = pd.read_csv(summary_tsv_or_df, sep='\t')
    samples = df['sample'].tolist()
    x       = range(len(df))

    # panel 1
    fig, axs = plt.subplots(3,1, figsize=(max(8,len(df)*0.4), 10), constrained_layout=True)
    axs[0].bar(x, df['fracSRTM'], label='SRTM')
    axs[0].bar(x, df['fracSNTM'], bottom=df['fracSRTM'], label='SNTM')
    axs[0].set_ylabel("SRTM/SNTM"); axs[0].legend(loc='center left', bbox_to_anchor=(1,0.5)); axs[0].set_xticks([])

    # panel 2
    axs[1].bar([i-0.2 for i in x], df["5'endsupport"], width=0.4, label="5'")
    axs[1].bar([i+0.2 for i in x], df["3'endsupport"], width=0.4, label="3'")
    axs[1].set_ylabel("end support"); axs[1].legend(loc='center left', bbox_to_anchor=(1,0.5)); axs[1].set_xticks([])

    # panel 3: isoform counts
    class_cols = [c for c in df.columns if c in ISOCLASSES]
    bottom = [0]*len(df)
    for cls in class_cols:
        axs[2].bar(x, df[cls], bottom=bottom, label=cls)
        bottom = [b+v for b,v in zip(bottom, df[cls])]
    axs[2].set_xticks(x)
    axs[2].set_xticklabels(samples, rotation=90, fontsize=6)
    axs[2].set_ylabel("isoform count")
    axs[2].legend(loc='center left', bbox_to_anchor=(1,0.5))

    plt.savefig(out_png, dpi=300)
    plt.close()
    logger.debug(f"[✓] Wrote plot {out_png}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument('-o','--outdir', required=True)
    p.add_argument('-s', '--summary', default=None,
                   help='Optional path to a single summary TSV to plot (skips glob).')
    args = p.parse_args()

    if args.summary:
        # User provided a single summary TSV to plot
        tsv = args.summary
        # derive a reasonable plot directory under outdir/plots/<basename>
        base = os.path.splitext(os.path.basename(tsv))[0]
        plot_dir = os.path.join(args.outdir, 'plots', base)
        logger.debug(f"Plotting {tsv} → {plot_dir}")
        plot_summary(tsv, plot_dir)
    else:
        pattern = os.path.join(args.outdir,'results','*','*','qc','sqanti','sqanti_results.tsv')
        for tsv in sorted(glob.glob(pattern)):
            parts = os.path.normpath(tsv).split(os.sep)
            region = parts[-5]
            runn   = parts[-4]
            plot_dir = os.path.join(args.outdir,'plots',region,runn)
            logger.debug(f"Plotting {tsv} → {plot_dir}")
            plot_summary(tsv, plot_dir)

if __name__ == '__main__':
    main()
