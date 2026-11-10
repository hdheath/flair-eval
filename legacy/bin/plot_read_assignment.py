#!/usr/bin/env python3
"""
Publication figure: TED isoforms and the reads the read map assigns to them, on a compressed x axis.

Tracks (top to bottom, one shared x axis): header line; CAGE 5' signal (gene strand); isoforms (one row each,
left label = short ID, assigned reads, MANE); reads grouped in the same order and colour, each group with dashed
guides at its isoform's 5' and 3' end and a right label with the share of reads whose OWN end stops > --off-bp
short of the isoform's end; a grey group of gene-strand reads that are in no isoform of the read map; dRNA 3'
signal (gene strand); compressed axis (break marks + true length, scale bar, window coordinates). Compression:
gaps covered by no exon of a shown isoform and no block of a drawn ASSIGNED read, longer than --intron-display-bp,
are drawn that wide (unassigned reads are drawn on the same transform, so intronic reads stay in their intron).

Isoforms shown: gene strand, >= --min-isoform-reads assigned reads, --keep first, then by read count up to
--max-isoforms (6 = the validated colour-blind-safe Okabe-Ito order). With --filter-bed (a filtered isoform set, same
IDs), shown isoforms missing from it are hatched by the rule that drops them (labels 5′ and 3′: / 5′ only, \\ 3′ only,
x both; other labels: dots) and marked † + label in the isoform column and the read-group label. --show-mismatches
draws black ticks at mismatched bases (=/X CIGAR or MD tag); --show-softclips draws soft clips as thin dark bars beyond
the aligned read end, at true length. A legend row at the bottom explains the marks used. With --annotation (BED12,
name = transcript_name|transcript_id|MANE), a reference track above the TED isoforms shows, for each shown isoform, the
annotated transcript with the most similar intron chain (junction Jaccard, then nearest ends), plus MANE; CDS thick,
UTR / non-coding thin; the right margin says which isoforms share its chain (=) or are closest to it (≈). Reads per group: at most
--reads-per-group,
sampled by a fixed hash so that the share of off-end reads is kept. Reuses the loaders, the intron-compression
transform and the signal binning of isoform_plot_test.py.

--compact: streamlined layout. Each isoform's model sits directly above its own reads (no separate model panel, no
right-hand text, no dashed end guides: the model directly above shows where reads should end); the left gutter holds one-word track names and "iN · reads" († = not in the filtered set); GENCODE
rows are unlabelled except MANE; --procap-signal adds a PRO-cap track under CAGE; unassigned reads are one grey coverage strip; the axis keeps only // breaks and a scale
bar; the only legend is the filter-hatch key. Coordinates, intron lengths, signal ranges and end-agreement percentages
belong in the caption (they are printed to stdout).
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import isoform_plot_test as IPT  # noqa: E402  (also applies pub_style rcParams)

import matplotlib.font_manager as fm  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import pysam  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.transforms import blended_transform_factory as blend  # noqa: E402

MM = 1 / 25.4
COLORS = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7"]  # Okabe-Ito, validated order
GREY = "#9A9A9A"
INK, INK2, INK3 = "#222222", "#555555", "#8A8A8A"
SIGNAL = "#6B7B8C"   # muted slate for both orthogonal tracks
MARK = "#111111"     # mismatches and soft clips on reads
ANNOT = "#3A3A3A"    # reference annotation track


def aln_marks(aln, want_mm=True):
    """(left soft clip, right soft clip, mismatch reference positions) of one alignment. Mismatches come from =/X
    CIGAR operations when present, else from the MD tag; without either none are reported."""
    ct = aln.cigartuples or []
    cl = ct[0][1] if ct and ct[0][0] == 4 else 0
    cr = ct[-1][1] if ct and ct[-1][0] == 4 else 0
    mm = []
    if want_mm:
        if any(op == 8 for op, _ in ct):
            r = aln.reference_start
            for op, n in ct:
                if op == 8:
                    mm.extend(range(r, r + n))
                if op in (0, 2, 3, 7, 8):
                    r += n
        elif aln.has_tag("MD"):
            mm = [rp for _, rp, base in aln.get_aligned_pairs(matches_only=True, with_seq=True) if base.islower()]
    return cl, cr, mm


def stream_readmap(path, isoform_ids, read_names):
    """iso -> reads for isoform_ids; read -> iso for any read in read_names (whole map)."""
    iso_reads, read_iso = {}, {}
    with open(path) as fh:
        for line in fh:
            i = line.find("\t")
            iso = line[:i]
            reads = line[i + 1:].rstrip("\n").split(",")
            if iso in isoform_ids:
                iso_reads[iso] = reads
            for r in reads:
                if r in read_names:
                    read_iso[r] = iso
    return iso_reads, read_iso


def chain(blocks):
    return tuple((blocks[i][1], blocks[i + 1][0]) for i in range(len(blocks) - 1))


def hash_rank(name):
    return hashlib.blake2b(name.encode(), digest_size=8, key=b"plot_read_assignment").hexdigest()


def sample_group(reads, k):
    """Keep at most k reads; keep the share of off-end reads (fixed hash order within each stratum)."""
    if len(reads) <= k:
        return list(reads)
    off = sorted([r for r in reads if r["off"]], key=lambda r: hash_rank(r["name"]))
    ok = sorted([r for r in reads if not r["off"]], key=lambda r: hash_rank(r["name"]))
    k_off = round(k * len(off) / len(reads))
    if off and k_off == 0:
        k_off = 1
    k_off = min(k_off, len(off))
    return off[:k_off] + ok[:k - k_off]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--readmap", required=True)
    ap.add_argument("--isoforms", required=True)
    ap.add_argument("--region", required=True, help="chr:start-end")
    ap.add_argument("--strand", required=True, choices=["+", "-"], help="gene strand")
    ap.add_argument("--gene", required=True)
    ap.add_argument("--sample-label", required=True, help="e.g. 'K562 · PacBio MAS-Iso-Seq · TED'")
    ap.add_argument("--cage-signal")
    ap.add_argument("--drna-signal")
    ap.add_argument("--procap-signal", help="gene-strand PRO-cap 5'-end bedgraph (compact layout: track under CAGE)")
    ap.add_argument("--cage-peaks")
    ap.add_argument("--drna-peaks")
    ap.add_argument("--cage-unit", default="TPM")
    ap.add_argument("--drna-unit", default="norm.")
    ap.add_argument("--mane", help="MANE Select BED12 (isoforms with the same intron chain are marked MANE)")
    ap.add_argument("--min-isoform-reads", type=int, default=10)
    ap.add_argument("--max-isoforms", type=int, default=len(COLORS))
    ap.add_argument("--keep", default="", help="comma list of isoform IDs always shown")
    ap.add_argument("--reads-per-group", type=int, default=22)
    ap.add_argument("--off-bp", type=int, default=50)
    ap.add_argument("--intron-display-bp", type=int, default=300)
    ap.add_argument("--filter-bed", action="append", default=[], metavar="LABEL=BED",
                    help="filtered isoform set (same IDs, column 4); shown isoforms missing from it are marked "
                         "† LABEL. Repeatable (e.g. 5′=… 3′=…); LABEL may be empty ('=path')")
    ap.add_argument("--filter-note", default="", help="legend text naming the filter, e.g. 'ladder L4'")
    ap.add_argument("--annotation", help="reference BED12 (name = transcript_name|transcript_id|MANE) of the gene")
    ap.add_argument("--annotation-label", default="GENCODE basic", help="track title")
    ap.add_argument("--show-mismatches", action="store_true", help="black ticks at mismatches (=/X CIGAR or MD tag)")
    ap.add_argument("--show-softclips", action="store_true", help="soft clips as thin bars beyond the read end")
    ap.add_argument("--compact", action="store_true", help="streamlined layout (see docstring)")
    ap.add_argument("--font-dir", help="directory with Arial*.ttf to register (compute nodes lack Arial)")
    ap.add_argument("--output", required=True, help="prefix; writes .pdf and .png (600 dpi)")
    a = ap.parse_args()

    if a.font_dir:
        for f in Path(a.font_dir).glob("*.ttf"):
            fm.fontManager.addfont(str(f))
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                         "savefig.bbox": "standard", "hatch.linewidth": 0.5})   # keep the exact 183-mm page width (pub_style sets "tight")

    chrom, rs, re_ = IPT.parse_region(a.region)
    st = a.strand
    max_iso = min(a.max_isoforms, len(COLORS))

    # ---- isoforms, reads, assignment -------------------------------------------------------------------------
    isoforms = IPT.filter_isoforms_to_region(IPT.load_isoforms_bed(Path(a.isoforms)), chrom, rs, re_)
    isoforms = {k: v for k, v in isoforms.items() if v["strand"] == st}
    alns = {}
    with pysam.AlignmentFile(a.bam) as bam:
        for aln in bam.fetch(chrom, rs, re_):
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.query_name in alns:
                continue
            alns[aln.query_name] = aln
    iso_reads, read_iso = stream_readmap(a.readmap, set(isoforms), set(alns))
    n_assigned = {k: len(iso_reads.get(k, [])) for k in isoforms}
    keep = [k for k in a.keep.split(",") if k in isoforms]
    ranked = sorted((k for k in isoforms if n_assigned[k] >= a.min_isoform_reads and k not in keep),
                    key=lambda k: -n_assigned[k])
    shown = (keep + ranked)[:max_iso]
    hidden = [k for k in isoforms if k not in shown and n_assigned[k] > 0]
    hidden_reads = sum(n_assigned[k] for k in hidden)

    def five(b):
        return b[0][0] if st == "+" else b[-1][1]

    def three(b):
        return b[-1][1] if st == "+" else b[0][0]

    sgn = 1 if st == "+" else -1   # transcript orientation
    shown.sort(key=lambda k: (sgn * five(isoforms[k]["blocks"]), sgn * three(isoforms[k]["blocks"]), -n_assigned[k]))
    label = {k: "i%d" % (i + 1) for i, k in enumerate(shown)}
    color = {k: COLORS[i] for i, k in enumerate(shown)}
    failed = {k: [] for k in shown}
    for spec in a.filter_bed:
        flab, _, fpath = spec.partition("=")
        kept = set()
        with open(fpath) as fh:
            for line in fh:
                f = line.split("\t", 4)
                if len(f) > 3 and f[3] in failed:
                    kept.add(f[3])
        for k in shown:
            if k not in kept:
                failed[k].append(flab)
    mark = {k: ("†" + " ".join(x for x in failed[k] if x)) if failed[k] else "" for k in shown}

    def reason(k):
        s = set(failed[k])
        if "5′" in s and "3′" in s:
            return "both ends", "xxxxxx"
        if "5′" in s:
            return "5′ end", "//////"
        if "3′" in s:
            return "3′ end", "\\\\\\"
        return "filtered", "......"
    reasons = {k: reason(k) for k in shown if failed[k]}

    mane_chains = set()
    if a.mane:
        for line in open(a.mane):
            f = line.rstrip("\n").split("\t")
            if f[0] != chrom or f[5] != st or int(f[2]) <= rs or int(f[1]) >= re_:
                continue
            s0 = int(f[1])
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            starts = [int(x) for x in f[11].rstrip(",").split(",")]
            mane_chains.add(chain([(s0 + b, s0 + b + z) for b, z in zip(starts, sizes)]))
    is_mane = {k: bool(chain(isoforms[k]["blocks"])) and chain(isoforms[k]["blocks"]) in mane_chains for k in shown}

    # ---- reference annotation: closest transcript per shown isoform, plus MANE -------------------------------------
    ann, ann_rows, match = {}, [], {}
    if a.annotation:
        for line in open(a.annotation):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or f[0] != chrom or f[5] != st or int(f[2]) <= rs or int(f[1]) >= re_:
                continue
            s0 = int(f[1])
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            starts = [int(x) for x in f[11].rstrip(",").split(",")]
            nm, tid, tag = (f[3].split("|") + ["", ""])[:3]
            ann[tid or nm] = dict(name=nm, mane=(tag == "MANE"), thick=(int(f[6]), int(f[7])),
                                  blocks=[(s0 + b, s0 + b + z) for b, z in zip(starts, sizes)])

        def sim(b1, b2):
            j1, j2 = set(chain(b1)), set(chain(b2))
            if j1 or j2:
                return len(j1 & j2) / len(j1 | j2)
            ov = max(0, min(b1[-1][1], b2[-1][1]) - max(b1[0][0], b2[0][0]))
            return ov / (max(b1[-1][1], b2[-1][1]) - min(b1[0][0], b2[0][0]))

        for k in shown:
            ib = isoforms[k]["blocks"]
            best = max(ann, default=None, key=lambda t: (sim(ib, ann[t]["blocks"]),
                       -abs(five(ib) - five(ann[t]["blocks"])) - abs(three(ib) - three(ann[t]["blocks"]))))
            if best is not None and sim(ib, ann[best]["blocks"]) > 0:
                match.setdefault(best, []).append((label[k], chain(ib) == chain(ann[best]["blocks"])))
        ann_rows = sorted(set(match) | {t for t in ann if ann[t]["mane"]},
                          key=lambda t: (sgn * five(ann[t]["blocks"]), sgn * three(ann[t]["blocks"])))

    def read_rec(name, iso):
        aln = alns[name]
        blocks = aln.get_blocks()
        rec = dict(name=name, blocks=blocks, off=False, d5=0, d3=0)
        if iso is not None:
            ib = isoforms[iso]["blocks"]
            rec["d5"] = sgn * (five(blocks) - five(ib))    # > 0: read starts downstream (inside) the isoform 5' end
            rec["d3"] = sgn * (three(ib) - three(blocks))  # > 0: read ends upstream (inside) the isoform 3' end
            rec["off"] = rec["d5"] > a.off_bp or rec["d3"] > a.off_bp
        return rec

    groups = []
    for k in shown:
        recs = [read_rec(r, k) for r in iso_reads.get(k, []) if r in alns]
        n = len(recs)
        p5 = 100.0 * sum(r["d5"] > a.off_bp for r in recs) / n if n else 0
        p3 = 100.0 * sum(r["d3"] > a.off_bp for r in recs) / n if n else 0
        pick = sample_group(recs, a.reads_per_group)
        pick.sort(key=lambda r: (r["d3"] > a.off_bp, r["d5"] > a.off_bp, r["d3"], r["d5"]))
        groups.append(dict(iso=k, n=n_assigned[k], found=n, p5=p5, p3=p3, reads=pick))
    unassigned = [n for n, aln in alns.items() if n not in read_iso and (aln.is_reverse == (st == "-"))]
    un_pick = sorted(unassigned, key=hash_rank)[:a.reads_per_group]
    un_recs = [read_rec(r, None) for r in un_pick]
    un_recs.sort(key=lambda r: (r["blocks"][0][0], r["blocks"][-1][1]))
    groups.append(dict(iso=None, n=len(unassigned), found=len(unassigned), reads=un_recs))

    # ---- x transform: exon scale, gaps > N without exons/reads drawn N wide --------------------------------------
    covered = [b for k in shown for b in isoforms[k]["blocks"]]
    covered += [b for g in groups if g["iso"] is not None for r in g["reads"] for b in r["blocks"]]
    covered += [(max(s, rs), min(e, re_)) for t in ann_rows for s, e in ann[t]["blocks"] if e > rs and s < re_]
    tx, segments, gaps = IPT.build_intron_compression(covered, rs, re_, a.intron_display_bp)
    x0, x1 = float(rs), tx(re_)

    if a.compact:
        draw_compact(a, locals())
        return

    # ---- layout (mm) --------------------------------------------------------------------------------------------
    W = 183.0
    left, right = 27.0, 38.0
    track_h, model_row, read_row, group_gap = 7.0, 2.7, 0.42, 1.6
    n_read_rows = sum(len(g["reads"]) for g in groups)
    reads_h = n_read_rows * read_row + group_gap * (len(groups) - 1)
    if reads_h > 92:
        read_row *= 92 / reads_h
        reads_h = 92
    models_h = max(1, len(shown)) * model_row + 1.0
    ann_h = len(ann_rows) * 2.3 + 1.0 if ann_rows else 0
    heights = [("header", 5.0), ("cage", track_h if a.cage_signal else 0), ("pad", 4.0), ("gencode", ann_h),
               ("pad", 3.4 if ann_h else 0), ("models", models_h),
               ("pad", 2.2), ("reads", reads_h), ("pad", 2.2), ("drna", track_h if a.drna_signal else 0),
               ("axis", 9.0)]
    want_legend = bool(reasons) or a.show_mismatches or a.show_softclips
    if want_legend:
        heights.append(("legend", 4.5))
    H = sum(h for _, h in heights)
    fig = plt.figure(figsize=(W * MM, H * MM))
    axes, top = {}, H
    for name, h in heights:
        top -= h
        if h <= 0 or name in ("pad", "header"):
            continue
        ax = fig.add_axes([left / W, top / H, (W - left - right) / W, h / H])
        ax.set_xlim(x0, x1)
        ax.axis("off")
        axes[name] = ax
    fs, fs_s = 7, 6

    def mm_x(mm):  # x position (figure fraction) from the left edge, in mm
        return mm / W

    # ---- header ---------------------------------------------------------------------------------------------
    y_head = (H - 2.8) / H
    t = fig.text(mm_x(left), y_head, a.gene, fontsize=8, fontstyle="italic", fontweight="bold", color=INK,
                 va="center", ha="left")
    fig.canvas.draw()
    bb = t.get_window_extent(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    arrow = "→" if st == "+" else "←"
    fig.text(bb.x1 + 1.2 / W, y_head,
             f"{a.sample_label}   {'+' if st == '+' else chr(0x2212)} strand {arrow}   "
             f"{chrom}:{rs + 1:,}–{re_:,}",
             fontsize=fs, color=INK2, va="center", ha="left")

    # ---- signal tracks ----------------------------------------------------------------------------------------
    def signal_track(ax, path, peaks_path, name, unit):
        span = re_ - rs
        xs, ys, _ = IPT.compute_binned_signal_from_bedgraphs([Path(path)], chrom, rs, re_,
                                                             max(300, min(40000, span // 10)))
        ymax = max(ys) if ys and max(ys) > 0 else 1.0
        ax.fill_between([tx(x) for x in xs], ys, color=SIGNAL, alpha=0.85, linewidth=0)
        ax.set_ylim(-0.30 * ymax, ymax * 1.05)
        ax.plot([x0, x1], [0, 0], color="#C8C8C8", lw=0.4)
        if peaks_path:
            for s, e in IPT.load_peaks(Path(peaks_path), chrom, rs, re_):
                ax.add_patch(Rectangle((tx(s), -0.24 * ymax), max(tx(e) - tx(s), (x1 - x0) * 0.004), 0.14 * ymax,
                                       facecolor=INK, edgecolor="none"))
        tf = blend(fig.transFigure, ax.transAxes)
        ax.text(mm_x(left - 2), 0.62, name, transform=tf, fontsize=fs, color=INK, ha="right", va="center")
        ax.text(mm_x(left - 2), 0.22, f"0–{ymax:.3g} {unit}", transform=tf, fontsize=fs_s, color=INK3,
                ha="right", va="center")

    if "cage" in axes:
        signal_track(axes["cage"], a.cage_signal, a.cage_peaks, "CAGE 5′ ends", a.cage_unit)
    if "drna" in axes:
        signal_track(axes["drna"], a.drna_signal, a.drna_peaks, "dRNA 3′ ends", a.drna_unit)

    # ---- reference annotation track ---------------------------------------------------------------------------
    if ann_rows:
        ax = axes["gencode"]
        ax.set_ylim(len(ann_rows) - 0.4, -0.6)
        tf = blend(fig.transFigure, ax.transData)
        wmin = (x1 - x0) * 0.0008
        for i, t in enumerate(ann_rows):
            b, (ts, te) = ann[t]["blocks"], ann[t]["thick"]
            ax.plot([tx(max(b[0][0], rs)), tx(min(b[-1][1], re_))], [i, i], color=ANNOT, lw=0.5,
                    solid_capstyle="butt")
            thin, thick = [], []
            for s, e in b:
                cs, ce = max(s, ts), min(e, te)
                if ce > cs:
                    thick.append((cs, ce))
                    thin += [(s, cs)] if cs > s else []
                    thin += [(ce, e)] if e > ce else []
                else:
                    thin.append((s, e))
            for parts, hh in ((thin, 0.32), (thick, 0.62)):
                if parts:
                    ax.broken_barh([(tx(s), max(tx(e) - tx(s), wmin)) for s, e in parts], (i - hh / 2, hh),
                                   facecolors=ANNOT, edgecolor="none")
            ax.text(mm_x(left - 25), i, ann[t]["name"], transform=tf, fontsize=fs_s - 0.5, color=INK, ha="left",
                    va="center")
            if ann[t]["mane"]:
                ax.text(mm_x(left - 9.3), i, "MANE", transform=tf, fontsize=fs_s - 0.5, color=INK2, ha="left",
                        va="center")
            same = [l for l, ex in match.get(t, []) if ex]
            near = [l for l, ex in match.get(t, []) if not ex]
            side_lab = "   ".join(x for x in (("= " + " ".join(same)) if same else "",
                                              ("≈ " + " ".join(near)) if near else "") if x)
            if side_lab:
                ax.text(mm_x(W - right + 2), i, side_lab, transform=tf, fontsize=fs_s - 0.5, color=INK2,
                        ha="left", va="center")
        ax.text(mm_x(left - 25), -1.05, f"{a.annotation_label} ({len(ann_rows)} of {len(ann)})", transform=tf,
                fontsize=fs_s - 0.5, color=INK3, ha="left", va="bottom")
        ax.text(mm_x(W - right + 2), -1.05, "= same intron chain   ≈ closest", transform=tf,
                fontsize=fs_s - 0.5, color=INK3, ha="left", va="bottom")

    # ---- isoform models -------------------------------------------------------------------------------------
    ax = axes["models"]
    n_rows = max(1, len(shown))
    ax.set_ylim(n_rows - 0.4, -0.6)
    pos = ax.get_position()
    tf = blend(fig.transFigure, ax.transData)
    for i, k in enumerate(shown):
        b = isoforms[k]["blocks"]
        ax.plot([tx(b[0][0]), tx(b[-1][1])], [i, i], color=color[k], lw=0.6, solid_capstyle="butt")
        if k in reasons:   # not in the filtered set: hatched by the rule that drops it
            ax.broken_barh([(tx(s), tx(e) - tx(s)) for s, e in b], (i - 0.30, 0.60), facecolors=color[k],
                           edgecolor="white", linewidth=0, hatch=reasons[k][1])
        else:
            ax.broken_barh([(tx(s), tx(e) - tx(s)) for s, e in b], (i - 0.30, 0.60), facecolors=color[k],
                           edgecolor="none")
        ax.text(mm_x(left - 25), i, f"{label[k]} {mark[k]}" if mark[k] else label[k], transform=tf,
                fontsize=fs_s, color=INK, ha="left", va="center")
        ax.text(mm_x(left - 10.5), i, f"{n_assigned[k]:,}", transform=tf, fontsize=fs_s, color=INK, ha="right",
                va="center")
        if is_mane[k]:
            ax.text(mm_x(left - 9.3), i, "MANE", transform=tf, fontsize=fs_s - 0.5, color=INK2, ha="left",
                    va="center")
    ax.text(mm_x(left - 25), -1.05, "TED", transform=tf, fontsize=fs_s - 0.5, color=INK3, ha="left",
            va="bottom")
    ax.text(mm_x(left - 10.5), -1.05, "reads", transform=tf, fontsize=fs_s - 0.5, color=INK3, ha="right",
            va="bottom")
    if not shown:
        ax.text((x0 + x1) / 2, 0, "no isoform passes the read threshold", fontsize=fs_s, color=INK3, ha="center")
    side = []
    if hidden:
        side.append(f"{len(hidden)} more isoform{'s' if len(hidden) > 1 else ''} not shown\n"
                    f"(< {a.min_isoform_reads} reads or beyond {max_iso}; {hidden_reads:,} reads)")
    if side:
        fig.text(mm_x(W - right + 2), pos.y0 + pos.height / 2, "\n\n".join(side),
                 fontsize=fs_s, color=INK3, ha="left", va="center", linespacing=1.2)

    # ---- reads -------------------------------------------------------------------------------------------------
    ax = axes["reads"]
    y = 0.0
    gap_rows = group_gap / read_row
    spans = []
    for g in groups:
        col = color[g["iso"]] if g["iso"] else GREY
        y_start = y
        for r in g["reads"]:
            b = r["blocks"]
            ax.plot([tx(b[0][0]), tx(b[-1][1])], [y, y], color=col, lw=0.25, alpha=0.45, solid_capstyle="butt")
            ax.broken_barh([(tx(s), max(tx(e) - tx(s), (x1 - x0) * 0.0008)) for s, e in b], (y - 0.36, 0.72),
                           facecolors=col, edgecolor="none")
            if a.show_mismatches or a.show_softclips:
                cl, cr, mm = aln_marks(alns[r["name"]], a.show_mismatches)
                if mm:
                    ax.broken_barh([(tx(q), max(tx(q + 1) - tx(q), (x1 - x0) * 0.0006)) for q in mm],
                                   (y - 0.36, 0.72), facecolors=MARK, edgecolor="none")
                if a.show_softclips:
                    xa, xb = tx(b[0][0]), tx(b[-1][1])
                    wmin = (x1 - x0) * 0.0008
                    clips = ([(xa - max(cl, wmin), max(cl, wmin))] if cl else []) + ([(xb, max(cr, wmin))] if cr else [])
                    if clips:
                        ax.broken_barh(clips, (y - 0.2, 0.4), facecolors=MARK, edgecolor="none")
            y += 1
        spans.append((g, y_start, y - 1))
        y += gap_rows
    ax.set_ylim(y - gap_rows + 0.2, -0.8)
    tf = blend(fig.transFigure, ax.transData)
    for g, ya, yb in spans:
        if g["iso"] is not None and g["reads"]:
            ib = isoforms[g["iso"]]["blocks"]
            for xg in (tx(five(ib)), tx(three(ib))):
                ax.plot([xg, xg], [ya - 0.7, yb + 0.7], color=INK, lw=0.45, ls=(0, (2.2, 1.6)), alpha=0.85)
        if g["iso"] is None:
            line1 = "not in read map"
            line2 = f"n = {g['n']:,} ({len(g['reads'])} shown)"
        else:
            lab = f"{label[g['iso']]} {mark[g['iso']]}" if mark[g["iso"]] else label[g["iso"]]
            line1 = f"{lab}   n = {g['n']:,} ({len(g['reads'])} shown)"
            if g["p3"] >= g["p5"] and g["p3"] >= 10:
                line2 = f"{g['p3']:.0f}% end > {a.off_bp} nt before 3′ end"
            elif g["p5"] >= 10:
                line2 = f"{g['p5']:.0f}% start > {a.off_bp} nt after 5′ end"
            else:
                line2 = "read ends match isoform"
        ax.text(mm_x(W - right + 2), (ya + yb) / 2, f"{line1}\n{line2}", transform=tf, fontsize=fs_s,
                color=INK if g["iso"] else INK3, ha="left", va="center", linespacing=1.25)

    # ---- x axis: breaks, true lengths, scale bar, window coordinates -----------------------------------------------
    ax = axes["axis"]
    ax.set_ylim(-1, 1)
    ax.plot([x0, x1], [0.75, 0.75], color=INK2, lw=0.5)
    per_mm = (x1 - x0) / (W - left - right)
    placed = []
    for d0, d1, length in sorted(gaps, key=lambda g: -g[2]):   # longest gaps get their length label first
        mid = (d0 + d1) / 2
        if all(abs(mid - m) >= 8.5 * per_mm for m in placed):
            placed.append(mid)
            ax.text(mid, 0.0, f"{length / 1000:.1f} kb", fontsize=5.5, color=INK3, ha="center", va="center")
    for d0, d1, length in gaps:
        ax.text((d0 + d1) / 2, 0.75, "//", fontsize=5.5, color=INK2, ha="center", va="center",
                bbox=dict(facecolor="white", edgecolor="none", pad=0.25))
    ax.text(x0, -0.8, f"{rs + 1:,}", fontsize=fs_s, color=INK2, ha="left", va="center")
    ax.text(x1, -0.8, f"{re_:,}", fontsize=fs_s, color=INK2, ha="right", va="center")
    nice = [50, 100, 200, 250, 500, 1000, 2000, 5000, 10000]
    bar = max([v for v in nice if v <= 12 * per_mm] or [nice[0]])   # scale bar <= 12 mm, in the right margin
    xb0 = x1 + 3 * per_mm
    ax.plot([xb0, xb0 + bar], [-0.8, -0.8], color=INK, lw=0.8, solid_capstyle="butt", clip_on=False)
    ax.text(xb0 + bar / 2, -0.55, f"{bar:,} bp" if bar < 1000 else f"{bar // 1000} kb", fontsize=fs_s,
            color=INK, ha="center", va="bottom", clip_on=False)
    ax.text(xb0, 0.75, f"// intron > {a.intron_display_bp} bp, shortened",
            fontsize=fs_s - 0.5, color=INK3, ha="left", va="center", clip_on=False)

    # ---- legend row: filter patterns and read marks ---------------------------------------------------------------
    if want_legend:
        rend = fig.canvas.get_renderer()
        yl = (dict(heights)["legend"] / 2) / H          # the legend row is the last row
        xpos = [mm_x(left)]

        def put_text(s, color=INK2):
            tt = fig.text(xpos[0], yl, s, fontsize=fs_s, color=color, ha="left", va="center")
            xpos[0] = tt.get_window_extent(rend).transformed(fig.transFigure.inverted()).x1 + 1.2 / W

        def put_box(w_mm, h_mm, **kw):
            fig.add_artist(Rectangle((xpos[0], yl - h_mm / 2 / H), w_mm / W, h_mm / H, transform=fig.transFigure, **kw))
            xpos[0] += (w_mm + 0.8) / W

        if reasons:
            put_text("† not in TED's filtered set" + (f" ({a.filter_note})" if a.filter_note else "") + ":")
            for name, hatch in [("5′ end", "//////"), ("3′ end", "\\\\\\"), ("both ends", "xxxxxx"),
                                ("filtered", "......")]:
                if any(v[0] == name for v in reasons.values()):
                    put_box(4.0, 1.7, facecolor=GREY, edgecolor="white", linewidth=0, hatch=hatch)
                    put_text(name)
            xpos[0] += 2.5 / W
        if a.show_mismatches:
            put_box(0.35, 1.7, facecolor=MARK, edgecolor="none")
            put_text("mismatch")
        if a.show_softclips:
            put_box(3.0, 0.7, facecolor=MARK, edgecolor="none")
            put_text("soft clip (true length)")

    out = Path(a.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out) + ".pdf", bbox_inches=None)
    fig.savefig(str(out) + ".png", dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"{a.gene}: {len(shown)} isoforms shown, {len(hidden)} hidden ({hidden_reads} reads), "
          f"{sum(len(g['reads']) for g in groups)} reads drawn, {len(gaps)} gaps compressed, {H:.0f} mm tall")
    for g in groups:
        print("  %-8s n=%-6d shown=%-3d %s" % (label.get(g["iso"], "unassigned"), g["n"], len(g["reads"]),
                                             "" if g["iso"] is None else "3'short %.0f%% 5'short %.0f%%  %s %s" % (
                                                 g["p3"], g["p5"], g["iso"], mark[g["iso"]])))


def draw_compact(a, v):
    """Streamlined layout; v = main()'s locals after the shared data preparation."""
    tx, x0, x1, rs, re_, chrom, st = v["tx"], v["x0"], v["x1"], v["rs"], v["re_"], v["chrom"], v["st"]
    isoforms, groups, shown, label, color = v["isoforms"], v["groups"], v["shown"], v["label"], v["color"]
    reasons, ann, ann_rows, alns = v["reasons"], v["ann"], v["ann_rows"], v["alns"]
    five, three, gaps, n_assigned = v["five"], v["three"], v["gaps"], v["n_assigned"]
    unassigned = v["unassigned"]

    W, left, right = 183.0, 14.0, 9.0
    fs, fs_s = 6.5, 5.5
    mrow, rrow, rh, bgap, trk, arow = 1.8, 0.42, 0.30, 1.3, 5.0, 1.8   # mm

    # vertical plan, mm from the top of the plotting body
    y, plan, blocks = 0.0, {}, []
    if a.cage_signal:
        plan["cage"] = (y, y + trk)
        y += trk + (1.0 if a.procap_signal else 1.8)
    if a.procap_signal:
        plan["procap"] = (y, y + trk)
        y += trk + 1.8
    if ann_rows:
        plan["gencode"] = (y, y + len(ann_rows) * arow)
        y += len(ann_rows) * arow + 1.8
    for g in groups[:-1]:
        yb0 = y
        y += mrow
        ry = [y + (j + 0.5) * rrow for j in range(len(g["reads"]))]
        y += len(g["reads"]) * rrow
        blocks.append((g, yb0, ry, y))
        y += bgap
    plan["other"] = (y, y + 2.2)
    y += 2.2 + 1.6
    if a.drna_signal:
        plan["drna"] = (y, y + trk)
        y += trk + 1.2
    plan["axis"] = (y, y + 3.4)
    body = y + 3.4
    head, leg = 5.0, (3.8 if reasons else 1.2)   # without the hatch key, keep a small bottom margin
    H = head + body + leg
    fig = plt.figure(figsize=(W * MM, H * MM))
    ax = fig.add_axes([left / W, leg / H, (W - left - right) / W, body / H])
    ax.set_xlim(x0, x1)
    ax.set_ylim(body, 0)
    ax.axis("off")
    gut = blend(fig.transFigure, ax.transData)
    wmin = (x1 - x0) * 0.0008

    def mmx(mm):
        return mm / W

    def gutter(yc, s, color=INK):
        ax.text(mmx(1.5), yc, s, transform=gut, fontsize=fs_s, color=color, ha="left", va="center")

    # header: gene, library, strand (coordinates go to the caption)
    yh = (H - 2.6) / H
    t = fig.text(mmx(1.5), yh, a.gene, fontsize=7.5, fontstyle="italic", fontweight="bold", color=INK, va="center")
    fig.canvas.draw()
    bb = t.get_window_extent(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bb.x1 + 1.5 / W, yh, f"{a.sample_label}   {'+' if st == '+' else chr(0x2212)} strand "
             f"{'→' if st == '+' else '←'}", fontsize=fs, color=INK2, va="center")

    def signal(path, key, name):
        ya, yb = plan[key]
        xs, ys, _ = IPT.compute_binned_signal_from_bedgraphs([Path(path)], chrom, rs, re_,
                                                             max(300, min(40000, (re_ - rs) // 10)))
        ymax = max(ys) if ys and max(ys) > 0 else 1.0
        ax.fill_between([tx(x) for x in xs], [yb - q / ymax * (yb - ya) for q in ys], yb, color=SIGNAL, linewidth=0)
        ax.plot([x0, x1], [yb, yb], color="#C8C8C8", lw=0.35)
        gutter(ya + (yb - ya) * 0.6, name)
        return ymax

    smax = {}
    if a.cage_signal:
        smax["CAGE"] = signal(a.cage_signal, "cage", "CAGE")
    if a.procap_signal:
        smax["PRO-cap"] = signal(a.procap_signal, "procap", "PRO-cap")
    if a.drna_signal:
        smax["dRNA"] = signal(a.drna_signal, "drna", "dRNA")

    if ann_rows:
        ya, yb = plan["gencode"]
        for i, t_ in enumerate(ann_rows):
            yc = ya + (i + 0.5) * arow
            b, (ts, te) = ann[t_]["blocks"], ann[t_]["thick"]
            ax.plot([tx(max(b[0][0], rs)), tx(min(b[-1][1], re_))], [yc, yc], color=ANNOT, lw=0.45)
            thin, thick = [], []
            for s, e in b:
                cs, ce = max(s, ts), min(e, te)
                if ce > cs:
                    thick.append((cs, ce))
                    thin += [(s, cs)] if cs > s else []
                    thin += [(ce, e)] if e > ce else []
                else:
                    thin.append((s, e))
            for parts, hh in ((thin, 0.5), (thick, 1.0)):
                if parts:
                    ax.broken_barh([(tx(s), max(tx(e) - tx(s), wmin)) for s, e in parts], (yc - hh / 2, hh),
                                   facecolors=ANNOT, edgecolor="none")
            if ann[t_]["mane"]:
                ax.text(mmx(W - right + 1.0), yc, "MANE", transform=gut, fontsize=fs_s, color=INK2, va="center")
        gutter((ya + yb) / 2, "GENCODE")

    for g, yb0, ry, yend in blocks:
        k = g["iso"]
        col, b = color[k], isoforms[k]["blocks"]
        ym = yb0 + mrow * 0.42
        ax.plot([tx(b[0][0]), tx(b[-1][1])], [ym, ym], color=col, lw=0.6)
        kw = (dict(facecolors=col, edgecolor="white", linewidth=0, hatch=reasons[k][1]) if k in reasons
              else dict(facecolors=col, edgecolor="none"))
        ax.broken_barh([(tx(s), tx(e) - tx(s)) for s, e in b], (ym - 0.55, 1.1), **kw)
        ax.add_patch(Rectangle((mmx(1.5), ym - 0.5), mmx(1.5), 1.0, transform=gut, facecolor=col, edgecolor="none",
                               clip_on=False))
        ax.text(mmx(3.6), ym, f"{label[k]}{'†' if k in reasons else ''} · {n_assigned[k]:,}", transform=gut,
                fontsize=fs_s, color=INK, ha="left", va="center")
        for r, yr in zip(g["reads"], ry):
            rb = r["blocks"]
            ax.plot([tx(rb[0][0]), tx(rb[-1][1])], [yr, yr], color=col, lw=0.2, alpha=0.5, solid_capstyle="butt")
            ax.broken_barh([(tx(s), max(tx(e) - tx(s), wmin)) for s, e in rb], (yr - rh / 2, rh), facecolors=col,
                           edgecolor="none")
            if a.show_mismatches or a.show_softclips:
                cl, cr, mm = aln_marks(alns[r["name"]], a.show_mismatches)
                if mm:
                    ax.broken_barh([(tx(q), max(tx(q + 1) - tx(q), (x1 - x0) * 0.0006)) for q in mm],
                                   (yr - rh / 2, rh), facecolors=MARK, edgecolor="none")
                if a.show_softclips:
                    xa, xb = tx(rb[0][0]), tx(rb[-1][1])
                    clips = (([(xa - max(cl, wmin), max(cl, wmin))] if cl else [])
                             + ([(xb, max(cr, wmin))] if cr else []))
                    if clips:
                        ax.broken_barh(clips, (yr - rh * 0.3, rh * 0.6), facecolors=MARK, edgecolor="none")

    # reads in no isoform of the read map: one coverage strip
    ya, yb = plan["other"]
    nb = 800
    bw = (x1 - x0) / nb
    cov = [0] * nb
    for name in unassigned:
        for s, e in alns[name].get_blocks():
            s, e = max(s, rs), min(e, re_)
            if e > s:
                for i in range(max(0, int((tx(s) - x0) / bw)), min(nb - 1, int((tx(e) - x0) / bw)) + 1):
                    cov[i] += 1
    cmax = max(cov) or 1   # square-root scale, so a local pile does not flatten the rest of the strip
    ax.fill_between([x0 + (i + 0.5) * bw for i in range(nb)], [yb - (c / cmax) ** 0.5 * (yb - ya) for c in cov], yb,
                    step="mid", color=GREY, linewidth=0)
    gutter((ya + yb) / 2, f"other · {len(unassigned):,}", INK2)

    # axis: break marks and a scale bar only
    ya, yb = plan["axis"]
    yl = ya + 0.7
    ax.plot([x0, x1], [yl, yl], color=INK2, lw=0.45)
    for d0, d1, _ in gaps:
        ax.text((d0 + d1) / 2, yl, "//", fontsize=5, color=INK2, ha="center", va="center",
                bbox=dict(facecolor="white", edgecolor="none", pad=0.2))
    per_mm = (x1 - x0) / (W - left - right)
    nice = [50, 100, 200, 250, 500, 1000, 2000, 5000, 10000]
    bar = max([q for q in nice if q <= 12 * per_mm] or [nice[0]])
    ys_ = yl + 1.9
    ax.plot([x1 - bar, x1], [ys_, ys_], color=INK, lw=0.7, solid_capstyle="butt")
    ax.text(x1 - bar - 1.0 * per_mm, ys_, f"{bar:,} bp" if bar < 1000 else f"{bar // 1000} kb", fontsize=fs_s,
            color=INK, ha="right", va="center")

    # legend: filter hatches only
    if reasons:
        rend = fig.canvas.get_renderer()
        yl_ = (leg / 2) / H
        xpos = [mmx(1.5)]

        def put_text(s):
            tt = fig.text(xpos[0], yl_, s, fontsize=fs_s, color=INK2, ha="left", va="center")
            xpos[0] = tt.get_window_extent(rend).transformed(fig.transFigure.inverted()).x1 + 1.2 / W

        put_text("† not in TED's filtered set" + (f" ({a.filter_note})" if a.filter_note else "") + ":")
        for name, hatch in [("5′ end", "//////"), ("3′ end", "\\\\\\"), ("both ends", "xxxxxx"),
                            ("filtered", "......")]:
            if any(q[0] == name for q in reasons.values()):
                fig.add_artist(Rectangle((xpos[0], yl_ - 0.8 / H), 4.0 / W, 1.6 / H, transform=fig.transFigure,
                                         facecolor=GREY, edgecolor="white", linewidth=0, hatch=hatch))
                xpos[0] += 4.8 / W
                put_text(name)

    out = Path(a.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out) + ".pdf", bbox_inches=None)
    fig.savefig(str(out) + ".png", dpi=600, bbox_inches=None)
    plt.close(fig)
    # caption material
    print(f"{a.gene} (compact): {chrom}:{rs + 1:,}-{re_:,} ({st}), {H:.0f} mm tall; "
          + "; ".join(f"{k_} 0-{q:.3g}" for k_, q in smax.items())
          + "; shortened gaps (kb): " + ", ".join(f"{L / 1000:.1f}" for _, _, L in gaps))
    for g in groups[:-1]:
        print("  %-4s n=%-6d shown=%-3d 3'short %.0f%% 5'short %.0f%%  %s %s" % (
            label[g["iso"]], g["n"], len(g["reads"]), g["p3"], g["p5"], g["iso"], v["mark"][g["iso"]]))
    print("  other n=%d" % len(unassigned))


if __name__ == "__main__":
    main()
