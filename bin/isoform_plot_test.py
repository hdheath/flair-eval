import pysam
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.colors as mcolors
import argparse


def plot_read_with_features(ax, read_alignment, ypos, base_color, alpha=1):
    """
    Plot read blocks with additional BAM features highlighted
    - Soft clipping shown in orange
    - Insertions shown as vertical lines in green
    - Deletions shown as gaps in red
    - Mismatches shown as red marks
    """
    cigar = read_alignment.cigartuples
    pos = read_alignment.reference_start
    query_pos = 0

    blocks = []
    soft_clips = []
    insertions = []
    deletions = []

    # Parse CIGAR to identify features
    for op, length in cigar:
        if op == 0:  # M - alignment match/mismatch
            blocks.append((pos, pos + length))
            query_pos += length
            pos += length
        elif op == 1:  # I - insertion
            insertions.append((pos, ypos))
            query_pos += length
        elif op == 2:  # D - deletion
            deletions.append((pos, pos + length))
            pos += length
        elif op == 3:  # N - intron/skipped region
            pos += length
        elif op == 4:  # S - soft clipping
            soft_clips.append((pos if query_pos == 0 else pos, length))
            query_pos += length
        elif op in {7, 8}:  # = or X - sequence match/mismatch
            blocks.append((pos, pos + length))
            query_pos += length
            pos += length

    if not blocks:
        return ax

    # Draw intron line (thin line connecting all blocks)
    r = mplpatches.Rectangle((blocks[0][0], ypos - .65),
                             blocks[-1][1] - blocks[0][0], 0.1,
                             facecolor=base_color, linewidth=0, alpha=alpha)
    ax.add_patch(r)

    # Draw exon blocks
    for bs, be in blocks:
        r = mplpatches.Rectangle((bs, ypos-1), be-bs, 0.8,
                                facecolor=base_color, linewidth=0, alpha=alpha)
        ax.add_patch(r)

    # Draw soft clips in bright yellow (highly visible, won't clash)
    for clip_pos, clip_len in soft_clips:
        r = mplpatches.Rectangle((clip_pos, ypos-1), clip_len, 0.8,
                                facecolor='#FFD700', linewidth=0, alpha=0.9)
        ax.add_patch(r)

    # Draw insertions as vertical lime green lines (very distinct from pink)
    for ins_pos, ins_y in insertions:
        ax.plot([ins_pos, ins_pos], [ins_y-1, ins_y-0.2],
               color='#00FF00', linewidth=2, alpha=0.9)

    # Draw deletions as bright red gaps (visible against grey and all isoform colors)
    for del_start, del_end in deletions:
        r = mplpatches.Rectangle((del_start, ypos-1), del_end-del_start, 0.8,
                                facecolor='#FF0000', linewidth=0, alpha=0.7)
        ax.add_patch(r)

    # Get mismatches from MD tag if available - use deep purple (distinct from all other colors)
    if read_alignment.has_tag('MD'):
        try:
            md_tag = read_alignment.get_tag('MD')
            ref_pos = read_alignment.reference_start
            aligned_pairs = read_alignment.get_aligned_pairs(with_seq=True)

            for qpos, rpos, refbase in aligned_pairs:
                if qpos is not None and rpos is not None and refbase:
                    if refbase.islower():  # Mismatch indicated by lowercase in MD
                        r = mplpatches.Rectangle((rpos, ypos-1), 1, 0.8,
                                                facecolor='#9400D3', linewidth=0, alpha=0.9)
                        ax.add_patch(r)
        except:
            pass  # Skip if MD tag parsing fails

    return ax


def plot_isoform_blocks(ax, blocks, ypos, color, alpha):
    """Plot isoform structure as thicker blocks"""
    r = mplpatches.Rectangle((blocks[0][0], ypos - .65),
                             blocks[-1][1]-blocks[0][0], 0.15,
                             facecolor=color, linewidth=0, alpha=alpha)
    ax.add_patch(r)
    for bs, be in blocks:
        r = mplpatches.Rectangle((bs, ypos-1), be-bs, 1.0,
                                facecolor=color, linewidth=0, alpha=alpha)
        ax.add_patch(r)
    return ax


def main():
    parser = argparse.ArgumentParser(description='Plot isoforms and reads from FLAIR output')
    parser.add_argument('--bam', required=True, help='BAM file')
    parser.add_argument('--readmap', required=True, help='Isoform read map file')
    parser.add_argument('--isoforms', required=True, help='Isoforms BED file')
    parser.add_argument('--output', required=True, help='Output prefix for plots')
    parser.add_argument('--cage', default=None, help='CAGE file (optional)')
    parser.add_argument('--quantseq', default=None, help='QuantSeq file (optional)')

    args = parser.parse_args()

    bamfile = args.bam
    readmapfile = args.readmap
    isoformsfile = args.isoforms
    outname = args.output

    # Create output directory if needed
    outdir = os.path.dirname(outname)
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)

    print('Loading data...')

    # Load isoform to reads mapping
    print('  - Reading read map')
    isotoreads = {}
    readtoiso = {}
    for line in open(readmapfile):
        iso, reads = line.rstrip().split('\t', 1)
        reads = reads.split(',')
        isotoreads[iso] = reads
        for r in reads:
            readtoiso[r] = iso

    # Load isoform structures
    print('  - Reading isoforms BED')
    isotoblocks = {}
    for line in open(isoformsfile):
        line = line.rstrip().split('\t')
        chrom, start, end = line[0], int(line[1]), int(line[2])
        isoname, strand = line[3], line[5]
        esizes, estarts = [int(x) for x in line[-2].rstrip(',').split(',')],  [int(x) for x in line[-1].rstrip(',').split(',')]
        exonblocks = [(start + estarts[i], start + estarts[i] + esizes[i]) for i in range(len(esizes))]
        isotoblocks[isoname] = [chrom, strand, start, end, exonblocks]

    # Load read alignments with full alignment objects first
    print('  - Reading BAM file')
    readtoblocks = {}
    readalignments = {}
    samfile = pysam.AlignmentFile(bamfile, 'rb')
    for s in samfile:
        if not s.is_supplementary and not s.is_secondary and s.is_mapped:
            readtoblocks[s.query_name] = s.is_reverse
            readalignments[s.query_name] = s
    samfile.close()

    print(f'Found {len(isotoblocks)} isoforms and {len(readtoblocks)} reads')

    # Determine plot region from isoforms or reads
    if isotoblocks:
        all_chroms = set(v[0] for v in isotoblocks.values())
        if len(all_chroms) > 1:
            print(f"Warning: Multiple chromosomes found: {all_chroms}")
        rchrom = list(all_chroms)[0]
        rstart = min(v[2] for v in isotoblocks.values())
        rend = max(v[3] for v in isotoblocks.values())
        print(f'  - Region: {rchrom}:{rstart}-{rend}')
    elif readtoblocks:
        # No isoforms but we have reads - use read alignment region
        print("No isoforms found - plotting unassigned reads only")
        rchrom = readalignments[list(readalignments.keys())[0]].reference_name
        rstart = min(s.reference_start for s in readalignments.values())
        rend = max(s.reference_end for s in readalignments.values())
        print(f'  - Region: {rchrom}:{rstart}-{rend}')
    else:
        print("Error: No isoforms or reads found")
        sys.exit(1)

    # Add unassigned reads
    noisoreads = [r for r in readtoblocks if r not in readtoiso]
    if noisoreads:
        isotoreads['unassigned'] = noisoreads

    # Sort isoforms by start position (leftmost first)
    isowithreadcount = []
    for i in isotoreads:
        if i in isotoblocks:
            start_pos = isotoblocks[i][2]  # Get start position
            isowithreadcount.append((i, len(isotoreads[i]), start_pos))
        else:
            # Unassigned reads - put at the end with max value
            isowithreadcount.append((i, len(isotoreads[i]), float('inf')))

    isowithreadcount.sort(key=lambda x:x[2])  # Sort by start position

    # Sort reads within each isoform by start position (leftmost first)
    for iso in isotoreads:
        myreads = [x for x in isotoreads[iso] if x in readtoblocks]
        readposinfo = [(readalignments[x].reference_start, x) for x in myreads]
        readposinfo.sort()  # Sort by start position
        isotoreads[iso] = [x[-1] for x in readposinfo]

    print('Plotting...')

    # Rich color palette for isoforms (distinct, saturated colors)
    isoform_colors = [
        '#E63946',  # Red
        '#1D3557',  # Navy
        '#F77F00',  # Orange
        '#06A77D',  # Teal
        '#A53860',  # Burgundy
        '#457B9D',  # Blue
        '#2A9D8F',  # Turquoise
        '#E76F51',  # Coral
        '#264653',  # Dark teal
        '#E9C46A',  # Gold
    ]

    # Calculate figure height
    mypos = len(readtoblocks) + (len(isotoblocks) * 2) + 2
    plt.figure(figsize=(12, max(6, mypos/10)))
    ax = plt.axes(frameon=False)

    # Assign colors to isoforms
    isoform_color_map = {}
    for idx, (iso, _, _) in enumerate(isowithreadcount):
        if iso != 'unassigned':
            base_color = isoform_colors[idx % len(isoform_colors)]
            cycle_num = idx // len(isoform_colors)

            if cycle_num == 0:
                # First cycle - use original colors
                isoform_color_map[iso] = base_color
            else:
                # Subsequent cycles - modify brightness
                # Convert hex to RGB
                rgb = mcolors.hex2color(base_color)
                # Darken by 20% for each cycle
                darken_factor = 0.8 ** cycle_num
                darker_rgb = tuple(c * darken_factor for c in rgb)
                isoform_color_map[iso] = mcolors.rgb2hex(darker_rgb)

            # Warn if too many cycles
            if cycle_num > 0 and idx == len(isoform_colors) * cycle_num:
                print(f'Warning: More than {len(isoform_colors) * cycle_num} isoforms detected. Colors will repeat with modified brightness.')

    # Plot each isoform and its reads
    for iso, isoreadcount, _ in isowithreadcount:
        ax.text(rstart, mypos + 0.5, f'{iso} (n={isoreadcount})', fontsize=7, ha='left', va='top')
        mypos -= 1.5  # More spacing between name and model

        if iso != 'unassigned':
            # Plot isoform structure
            iso_color = isoform_color_map[iso]
            ax = plot_isoform_blocks(ax, isotoblocks[iso][4], mypos, iso_color, 1)
        mypos -= 1

        # Plot reads
        for r in isotoreads[iso]:
            if r in readalignments:
                if iso != 'unassigned':
                    # Use lighter version of isoform color for reads
                    base_color = isoform_color_map[iso]
                    # Lighten the color
                    rgb = mcolors.hex2color(base_color)
                    lighter_rgb = tuple(min(1.0, c + 0.3) for c in rgb)
                    read_color = mcolors.rgb2hex(lighter_rgb)
                else:
                    # Grey for unassigned reads
                    read_color = 'grey'

                ax = plot_read_with_features(ax, readalignments[r], mypos, read_color, 1)
                mypos -= 1
        mypos -= 1

    ax.set_ylim((mypos-1, len(readtoblocks) + (len(isotoblocks) * 2) + 2))
    ax.set_yticks([])
    ax.set_xlim((rstart, rend))
    ax.set_xlabel(f'{rchrom}:{rstart:,}-{rend:,}')

    # Add legend for features - place outside plot area to avoid overlap
    legend_elements = [
        mplpatches.Patch(facecolor='#FFD700', label='Soft clipping'),
        mplpatches.Patch(facecolor='#FF0000', alpha=0.7, label='Deletion'),
        mplpatches.Patch(facecolor='#9400D3', label='Mismatch'),
        plt.Line2D([0], [0], color='#00FF00', linewidth=2, label='Insertion'),
        mplpatches.Patch(facecolor='grey', label='Unassigned reads'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.01, 1),
             fontsize=8, framealpha=0.95, borderpad=1)

    print(f'Saving to {outname}.png')
    plt.savefig(outname + '.png', dpi=300, bbox_inches='tight')
    print('Done!')


if __name__ == '__main__':
    main()
