#!/usr/bin/env python3
"""
plot_genome_browser_v13.py (reads + isoform bars only)

- Keeps original row-packing logic:
  * per-isoform relative packing + global base scan with adjacency guards
  * reserved rows scaffold (space1, collapsed_bar, ...)
  * flip rows at the end so larger stacks appear on top
  * unassigned reads packed with 3-row guard
- Removes ALL GTF and ALL TSS/TTS code/plots
- Fixes isoform bars: draw at reserved 'collapsed_bar' row (br)

JSON config keys (unchanged):
  bam, gtf, mapping, collapsed_isoforms, genome, outdir,
  gene_height, gene_row_height, read_row_height, fig_width, iso_kde
"""

import json
import argparse
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple
from collections import defaultdict
from itertools import chain

import pandas as pd
import pysam
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from intervaltree import IntervalTree
from matplotlib.collections import PatchCollection, LineCollection


@dataclass(frozen=True)
class Config:
    gtf: Path
    bam: Optional[Path] = None
    genome: str = ""
    outdir: Path = Path(".")
    mapping: Optional[Path] = None
    collapsed_isoforms: Optional[Path] = None
    gene_height: float = 1.5
    gene_row_height: float = 3.0      # kept for backward compat; not used now
    read_row_height: float = 0.02
    fig_width: float = 12.0
    iso_kde: bool = True              # kept for compat; not used now

    @staticmethod
    def from_json(p: Path) -> "Config":
        cfg = json.loads(p.read_text())
        return Config(
            bam=Path(cfg["bam"]) if cfg.get("bam") else None,
            gtf=Path(cfg["gtf"]),
            genome=cfg.get("genome", ""),
            outdir=Path(cfg.get("outdir", ".")),
            mapping=Path(cfg["mapping"]) if cfg.get("mapping") else None,
            collapsed_isoforms=Path(cfg["collapsed_isoforms"]) if cfg.get("collapsed_isoforms") else None,
            gene_height=cfg.get("gene_height", 1.5),
            gene_row_height=cfg.get("gene_row_height", 3.0),
            read_row_height=cfg.get("read_row_height", 0.02),
            fig_width=cfg.get("fig_width", 12.0),
            iso_kde=cfg.get("iso_kde", True),
        )


def parse_args():
    p = argparse.ArgumentParser(
        description="Reads + collapsed isoform bars (no genes, no TSS/TTS).",
    )
    p.add_argument(
        "-c", "--config", type=Path, required=True,
        help=("JSON config keys: (bam), gtf, mapping, collapsed_isoforms, genome, outdir, "
              "gene_height, gene_row_height, read_row_height, fig_width, iso_kde"),
    )
    p.add_argument("--region", help="Limit plotting region (e.g. chr1:100-200)")
    return p.parse_args()


def pleasant_colors(seq):
    def ok(rgb):
        mx, mn = max(rgb), min(rgb)
        return (mx < 0.9) and not (mx - mn < 0.2 and mx > 0.7)
    return [c for c in seq if ok(c)]


def read_strand(rd: pysam.AlignedSegment) -> str:
    xs = rd.get_tag("XS") if rd.has_tag("XS") else None
    if xs in ("+", "-"):
        return xs
    return "-" if rd.is_reverse else "+"


def overlaps(tree: IntervalTree, blks) -> bool:
    return any(tree.overlap(s, e) for s, e in blks)


def add_blocks(tree: IntervalTree, blks) -> None:
    for s, e in blks:
        tree.addi(s, e, True)


def ensure_rows(occupancy, upto):
    missing = upto - (len(occupancy) - 1)
    for _ in range(max(0, missing)):
        occupancy.append(IntervalTree())


def reserve_iso_rows(base: int, max_rel: int) -> dict[str, int]:
    return {
        "space1": base + max_rel + 1,
        "collapsed_bar": base + max_rel + 2,  # <— isoform bar row
        "space2": base + max_rel + 4,
        "empty1": base + max_rel + 6,
        "space3": base + max_rel + 9,
        "empty2": base + max_rel + 11,
        "space4": base + max_rel + 14,
    }


def load_mapping(map_path: Optional[Path]):
    m: dict[str, str] = {}
    if not map_path or not map_path.exists():
        logging.warning(f"Mapping not found: {map_path}")
        return m
    with map_path.open() as f:
        for line in f:
            iso, reads = line.strip().split('\t')
            for r in reads.split(','):
                r = r.strip()
                if r:
                    m[r] = iso
    if not m:
        logging.warning("Mapping is empty.")
    return m


def assign_read_rows(blocks_list):
    """
    Row assignment using interval trees for fast overlap checks.
    Returns a list of row indices and total number of rows.
    """
    rows = []  # list of IntervalTrees
    assign = []
    for blks in blocks_list:
        placed = False
        for ridx, tree in enumerate(rows):
            # ensure none of these blocks overlap existing intervals
            if all(len(tree.overlap(b[0], b[1])) == 0 for b in blks):
                # insert blocks into this row's tree
                for b in blks:
                    tree[b[0]:b[1]] = True
                assign.append(ridx)
                placed = True
                break
        if not placed:
            tree = IntervalTree()
            for b in blks:
                tree[b[0]:b[1]] = True
            rows.append(tree)
            assign.append(len(rows) - 1)
    return assign, len(rows)


def _parse_region(region: Optional[str]) -> Tuple[Optional[Tuple[str, int, int]], bool]:
    """
    Parse a region string like ``"chr1:100-200"`` into a tuple and flag
    whether the span is too large for plotting.

    Returns a tuple ``(region_tuple, skip_plot)`` where ``region_tuple`` is
    ``(chrom, start, end)`` if parsing succeeds or ``None`` otherwise, and
    ``skip_plot`` is ``True`` when the span exceeds the 20kb plotting limit.
    """
    if not region:
        return None, False
    m = re.match(r"^([^:]+):(\d+)-(\d+)$", region.replace(",", ""))
    if not m:
        logging.warning(f"Could not parse region string: {region}")
        return None, False
    chrom, start_s, end_s = m.groups()
    start_i, end_i = int(start_s), int(end_s)
    # Do not decide to skip plotting here; leave span-based gating to the
    # caller (e.g., TED) so that policy is centralized. Always return the
    # parsed region tuple and let the caller decide whether to invoke the
    # plotting routine for very large regions.
    return (chrom, start_i, end_i), False


def generate(cfg: Config, region: Optional[str] = None) -> Optional[Path]:
    """Generate plot from configuration (reads + isoform bars only).

    Returns the Path to the saved PNG on success, or ``None`` if the plot was
    skipped (e.g., due to missing inputs).
    """
    region_tuple, skip_plot = _parse_region(region)
    # If region string was supplied but couldn't be parsed, bail out
    if region and region_tuple is None:
        return
    chrom = None
    r0 = r1 = None
    if region_tuple:
        chrom, r0, r1 = region_tuple

    # Resolve BAM; gtf used only to locate regionalized BAM if --region specified
    gtf = cfg.gtf
    bam = None
    if region_tuple:
        tag = f"{chrom}_{r0}_{r1}"
        run_root = gtf.parent.parent.parent
        reg_root = run_root / "regionalize"
        if reg_root.exists():
            for d in reg_root.iterdir():
                cand = d / f"{tag}.bam"
                if cand.exists():
                    bam = cand
                    break
    if bam is None:
        bam = cfg.bam

    if not bam or not bam.exists():
        logging.warning(f"BAM not found: {bam}")
        return
    if not gtf.exists():
        # Only warn (still allowed to plot); we keep this because gtf path may be used for regionalized BAM resolution
        logging.warning(f"GTF not found (only used for regionalized BAM lookup): {gtf}")

    mapping_file = cfg.mapping if cfg.mapping and cfg.mapping.exists() else None
    if cfg.mapping and not mapping_file:
        logging.warning(f"Mapping not found: {cfg.mapping}")
    cis_bed = cfg.collapsed_isoforms if cfg.collapsed_isoforms and cfg.collapsed_isoforms.exists() else None
    if cfg.collapsed_isoforms and not cis_bed:
        logging.warning(f"Collapsed isoforms BED not found: {cfg.collapsed_isoforms}")

    genome = cfg.genome
    outdir = cfg.outdir
    gene_h = cfg.gene_height
    read_row_h = cfg.read_row_height
    fig_w = cfg.fig_width
    outdir.mkdir(parents=True, exist_ok=True)

    # load reads
    reads_blocks, reads_introns, read_strands, read_names = [], [], [], []
    contig_for_title = None
    with pysam.AlignmentFile(bam, "rb") as bf:
        iterator = bf.fetch(*region_tuple) if region_tuple else bf.fetch()
        for rd in iterator:
            # Only consider primary alignments.
            if rd.is_unmapped or rd.is_secondary or rd.is_supplementary:
                continue
            if contig_for_title is None:
                contig_for_title = bf.get_reference_name(rd.reference_id)
            blks = rd.get_blocks()
            if not blks:
                continue
            # Clip blocks to region bounds if provided
            if r0 is not None:
                clipped = []
                for s, e in blks:
                    s2 = max(s, r0)
                    e2 = min(e, r1)
                    if e2 > s2:
                        clipped.append((s2, e2))
                blks = clipped
                if not blks:
                    continue
            reads_blocks.append(blks)
            reads_introns.append([(blks[i][1], blks[i + 1][0]) for i in range(len(blks) - 1)])
            read_strands.append(read_strand(rd))
            read_names.append(rd.query_name)
    if not reads_blocks:
        logging.warning("No reads loaded.")
        return

    all_starts = [b[0] for blks in reads_blocks for b in blks]
    all_ends = [b[1] for blks in reads_blocks for b in blks]
    read_min, read_max = min(all_starts), max(all_ends)

    # load mapping & collapsed isoforms (keep 'chr', then later filter per contig/region)
    mapping = load_mapping(mapping_file)

    cis_df = pd.DataFrame(columns=[
        "chr","iso_id","start","end","strand","blockCount","blockSizes","blockStarts"
    ])
    iso_blocks: dict[str, list[tuple[int, int]]] = {}

    if cis_bed and cis_bed.exists():
        tmp = pd.read_csv(cis_bed, sep="\t", header=None, comment="#")
        if tmp.shape[1] < 12:
            logging.warning("Collapsed isoforms BED missing BED12 fields.")
            # Map what we can (BED6-style)
            tmp = tmp.iloc[:, :min(tmp.shape[1], 6)]
            # Pad/truncate columns to expected names
            cols = ["chr","start","end","iso_id","score","strand"][:tmp.shape[1]]
            tmp.columns = cols
            tmp["blockCount"] = np.nan
            tmp["blockSizes"] = ""
            tmp["blockStarts"] = ""
        else:
            tmp = tmp.iloc[:, :12]
            tmp.columns = [
                "chr","start","end","iso_id","score","strand",
                "thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"
            ]
        cis_df = tmp[["chr","iso_id","start","end","strand","blockCount","blockSizes","blockStarts"]]

        # Build iso_blocks (exon rectangles) from BED12 if present
        for _, r in cis_df.iterrows():
            sizes = [int(x) for x in str(r.blockSizes).rstrip(",").split(",") if x]
            starts = [int(x) for x in str(r.blockStarts).rstrip(",").split(",") if x]
            if sizes and starts:
                blks = [(int(r.start) + st, int(r.start) + st + sz) for st, sz in zip(starts, sizes)]
            else:
                blks = [(int(r.start), int(r.end))]
            iso_blocks[r.iso_id] = blks

    # group by isoform
    iso_to_idxs, unassigned = defaultdict(list), []
    for i, name in enumerate(read_names):
        iso = mapping.get(name)
        if iso:
            iso_to_idxs[iso].append(i)
        else:
            unassigned.append(i)

    base_colors = list(chain(plt.cm.tab20.colors, plt.cm.Set3.colors, plt.cm.Dark2.colors))
    pal = pleasant_colors(base_colors)
    isos = sorted(iso_to_idxs, key=lambda i: len(iso_to_idxs[i]), reverse=True)
    iso_colors = {iso: pal[i % len(pal)] for i, iso in enumerate(isos)}
    ua_color = 'lightgrey'

    # layout
    occupancy, row_assign, collapsed = [], [None] * len(reads_blocks), {}
    ih, rh = gene_h * 0.6, gene_h * 0.6

    for iso in isos:
        idxs = iso_to_idxs[iso]
        sorted_idxs = sorted(idxs, key=lambda i: reads_blocks[i][0][0])
        blks_list, _ = zip(*[(reads_blocks[i], None) for i in sorted_idxs])
        rel, _ = assign_read_rows(blks_list)

        # Collapsed isoform span from BED; if missing, treat reads as unassigned (preserve OG behavior)
        st, en = (None, None)
        if not cis_df.empty and iso in set(cis_df["iso_id"]):
            row = cis_df.loc[cis_df["iso_id"] == iso].iloc[0]
            st, en = int(row.start), int(row.end)
        if st is None:
            logging.warning(f"No collapsed bar for {iso}")
            unassigned.extend(idxs)
            continue

        left = min(reads_blocks[j][0][0] for j in idxs)
        right = max(reads_blocks[j][-1][1] for j in idxs)
        max_rel = max(rel)

        # Find a base where this iso's stack + reserved rows fit, with adjacency guard
        for base in range(len(occupancy) + 1):
            ok = True
            for ri, blks in zip(rel, blks_list):
                for off in (0, 1):
                    tgt = base + ri + off
                    if tgt < len(occupancy) and overlaps(occupancy[tgt], blks):
                        ok = False
                        break
                if not ok:
                    break
            if not ok:
                continue
            spec_rows = reserve_iso_rows(base, max_rel)
            for key, pr in spec_rows.items():
                seg = (left, right) if key not in ("space1", "collapsed_bar") else (st, en)
                if pr < len(occupancy) and overlaps(occupancy[pr], [seg]):
                    ok = False
                    break
            if ok:
                break
        else:
            base = len(occupancy)
            spec_rows = reserve_iso_rows(base, max_rel)

        ensure_rows(occupancy, spec_rows["space4"])
        for ri, idx in zip(rel, sorted_idxs):
            r = base + ri
            row_assign[idx] = r
            add_blocks(occupancy[r], reads_blocks[idx])
        for key in ("space1", "collapsed_bar", "space2", "empty1", "space3", "empty2", "space4"):
            r = spec_rows[key]
            seg = (left, right) if key not in ("space1", "collapsed_bar") else (st, en)
            add_blocks(occupancy[r], [seg])
        collapsed[iso] = (
            spec_rows["space1"],
            spec_rows["collapsed_bar"],  # <- use this for bar y
            spec_rows["space2"],
            spec_rows["empty1"],
            spec_rows["space3"],
            spec_rows["empty2"],
            spec_rows["space4"],
        )

    # place unassigned reads
    for idx in unassigned:
        blks = reads_blocks[idx]
        placed = False
        for r in range(len(occupancy)):
            if overlaps(occupancy[r], blks): continue
            if r > 0 and overlaps(occupancy[r - 1], blks): continue
            if r < len(occupancy) - 1 and overlaps(occupancy[r + 1], blks): continue
            row_assign[idx] = r
            add_blocks(occupancy[r], blks)
            placed = True
            break
        if not placed:
            t = IntervalTree()
            add_blocks(t, blks)
            occupancy.append(t)
            row_assign[idx] = len(occupancy) - 1

    # Flip rows so stacks are at top
    tot = len(occupancy)
    row_assign = [tot - 1 - r for r in row_assign]
    collapsed = {iso: tuple(tot - 1 - r for r in rows) for iso, rows in collapsed.items()}

    # Build iso span map for drawing bars/exons (BED12 blocks if present)
    ci_map = {r.iso_id: (int(r.start), int(r.end)) for _, r in cis_df.iterrows()} if not cis_df.empty else {}

    # ---- Plot (reads + iso bars) ----
    # Figure height: scale by total rows; add small margin
    fig_h = max(3.5, read_row_h * (tot + 12))
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=600)
    ax = fig.add_subplot(1, 1, 1)
    # Remove plot border/spines
    for side in ['left', 'right', 'bottom', 'top']:
        ax.spines[side].set_visible(False)

    # Draw reads (x-ordered) with iso colors
    order, drawn = (sorted(range(len(reads_blocks)), key=lambda i: reads_blocks[i][0][0]), set())
    rects, lines, arrows = [], [], []
    for i in order:
        blks = reads_blocks[i]; intr = reads_introns[i]
        strand = read_strands[i]; rn = read_names[i]; y = row_assign[i]
        grp = mapping.get(rn, '__UNASSIGNED__'); col = iso_colors.get(grp, ua_color)
        for s, e in blks: rects.append(Rectangle((s, y - rh/2), e - s, rh, facecolor=col, linewidth=0))
        for a, b in intr: lines.append([(a, y), (b, y)])
        arr = '>' if strand == '+' else '<'
        arrows += [(blks[0][0], y, 'left', arr), (blks[-1][1], y, 'right', arr)]
    if rects: ax.add_collection(PatchCollection(rects, match_original=True))
    if lines: ax.add_collection(LineCollection(lines, colors='black', linewidths=0.5))
    for x, y, ha, arr in arrows: ax.text(x, y, arr, ha=ha, va='center', fontsize=2, color='white')

    # Draw isoform bars + labels above grouped isoform reads (single block, no duplicate)
    for iso in isos:
        if iso not in collapsed:
            continue
        s1, br, s2, er, s3, e2, s4 = collapsed[iso]
        if iso in ci_map:
            cs, ce = ci_map[iso]
        else:
            idxs = iso_to_idxs[iso]
            cs = min(reads_blocks[j][0][0] for j in idxs)
            ce = max(reads_blocks[j][-1][1] for j in idxs)
        if r0 is not None:
            cs = max(cs, r0); ce = min(ce, r1)
            if ce <= cs:
                continue
        blks_iso = iso_blocks.get(iso, [(cs, ce)])
        if r0 is not None:
            blks_iso = [(max(bs, r0), min(be, r1)) for (bs, be) in blks_iso if be > r0 and bs < r1]
            if not blks_iso:
                continue
        color = iso_colors.get(iso, 'black')
        idxs = iso_to_idxs[iso]
        if idxs:
            max_row = max(row_assign[j] for j in idxs)
            min_row = min(row_assign[j] for j in idxs)
            # Add extra spacing above the isoform bar and below the label to avoid overlap
            bar_y = max_row + rh * 1.2 + ih * 0.5
            label_y = bar_y + ih * 1.5
            # Also, ensure the first assigned read is not too close to the bar
            # (by increasing bar_y and label_y as above)
        else:
            # If no assigned reads, place bar at reserved row
            bar_y = br
            label_y = br + ih * 1.5
        span_left = min(bs for bs, _ in blks_iso); span_right = max(be for _, be in blks_iso)
        for bs, be in blks_iso:
            ax.add_patch(Rectangle((bs, bar_y - ih/2), be - bs, ih,
                                   facecolor=color, edgecolor='black', linewidth=0.1, zorder=3))
        intr = [(blks_iso[i][1], blks_iso[i+1][0]) for i in range(len(blks_iso) - 1)]
        if intr:
            ax.add_collection(LineCollection([((a, bar_y), (b, bar_y)) for a, b in intr],
                                             colors='black', linewidths=0.5, zorder=3))
        iso_label = iso.split('_')[-1] if '_' in iso else iso
        ax.add_patch(Rectangle((span_left, label_y - 0.8), span_right - span_left, 1.6,
                               facecolor='white', edgecolor='none', alpha=0.7, zorder=9))
        ax.text(span_left, label_y, iso_label, fontsize=8, ha='left', va='bottom', color=color, zorder=10)

    # Label unassigned reads as 'Unassigned'
    if unassigned:
        ua_rows = [row_assign[idx] for idx in unassigned if row_assign[idx] is not None]
        if ua_rows:
            max_ua_row = max(ua_rows)
            label_y = max_ua_row + ih * 1.2
            left = min(reads_blocks[idx][0][0] for idx in unassigned)
            right = max(reads_blocks[idx][-1][1] for idx in unassigned)
            ax.add_patch(Rectangle((left, label_y - 0.8), right - left, 1.6,
                                   facecolor='white', edgecolor='none', alpha=0.7, zorder=9))
            ax.text(left, label_y, 'Unassigned', fontsize=8, ha='left', va='bottom', color='gray', zorder=10)

    # Labels & limits
    ax.set_yticks([])
    contig = chrom if chrom else (contig_for_title or "unknown")
    # X limits prefer region, else reads span with padding
    if r0 is not None:
        xmin, xmax = max(r0, read_min), min(r1, read_max)
        if xmax <= xmin:
            xmin, xmax = r0, r1
    else:
        xmin, xmax = read_min, read_max
    xpad = max((xmax - xmin) * 0.05, 50)
    ax.set_xlim(xmin - xpad, xmax + xpad)

    all_rows = []
    all_rows.extend([r for r in row_assign if r is not None])
    all_rows.extend([r for rows in collapsed.values() for r in rows])
    y_lo = (min(all_rows) if all_rows else 0) - 2
    y_hi = (max(all_rows) if all_rows else 10) + 6
    ax.set_ylim(y_lo, y_hi)

    # Set x-axis label and scientific notation for large values
    ax.set_xlabel(f"{contig}")
    ax.ticklabel_format(style='sci', axis='x', scilimits=(6,6))
    fig.subplots_adjust(bottom=0.18)  # Add more space at the bottom for x-label

    # save
    if skip_plot:
        plt.close(fig)
        return None
    out_png = outdir / f"{genome}_{contig}_{xmin}-{xmax}.png"
    fig.savefig(out_png, dpi=fig.dpi)
    plt.close(fig)
    width_px, height_px = int(fig_w * fig.dpi), int(fig_h * fig.dpi)
    logger = logging.getLogger(__name__)
    logger.debug(f"Figure dimensions: {width_px} x {height_px} pixels")
    logger.debug(f"Saved: {out_png}")
    return out_png


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    cfg = Config.from_json(args.config)
    generate(cfg, args.region)


if __name__ == '__main__':
    main()

