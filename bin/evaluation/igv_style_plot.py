#!/usr/bin/env python3
"""IGV-style locus renderer (matplotlib). Each isoform gets its own row — the thing
headless IGV's `expand` refused to do. Track types: coverage (BAM), isoforms (BED12, one
row per transcript), annotation (BED12 gene models), peaks (BED ticks), region (boxes),
sequence (base-zoom colored letters). Config-driven so the same code renders any locus."""
import argparse, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
import pysam

# ---------- parsers ----------
def parse_bed12(path, chrom, start, end):
    out = []
    for line in open(path):
        if line.startswith(("track", "#", "browser")):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 12 or c[0] != chrom:
            continue
        s, e = int(c[1]), int(c[2])
        if e < start or s > end:
            continue
        sizes = [int(x) for x in c[10].rstrip(",").split(",")]
        st = [int(x) for x in c[11].rstrip(",").split(",")]
        exons = [(s + st[i], s + st[i] + sizes[i]) for i in range(len(sizes))]
        out.append({"name": c[3], "start": s, "end": e, "strand": c[5], "exons": exons})
    return out

def feat_jc(f):
    ex = sorted(f["exons"])
    if len(ex) < 2:
        return None
    return format(hash(tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))) & 0xFFFF, "04x")


def load_counts(path, col):
    d = {}
    for L in open(path):
        c = L.rstrip("\n").split("\t")
        if not c or c[0].startswith("#") or c[0] in ("TXNAME", "feature_id"):
            continue
        try:
            d[c[0]] = float(c[col])
        except (IndexError, ValueError):
            pass
    return d


def parse_bed_peaks(path, chrom, start, end):
    out = []
    for line in open(path):
        if line.startswith(("track", "#")):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 3 or c[0] != chrom:
            continue
        s, e = int(c[1]), int(c[2])
        if e < start or s > end:
            continue
        out.append({"start": s, "end": e, "name": c[3] if len(c) > 3 else "",
                    "strand": c[5] if len(c) > 5 else "."})
    return out

# ---------- drawing ----------
def draw_tx(ax, y, h, tx, color, x0, x1, arrow_px, span_data):
    a, b = max(tx["start"], x0), min(tx["end"], x1)
    if b > a:
        ax.plot([a, b], [y, y], color=color, lw=0.9, zorder=2, solid_capstyle="butt")
    for es, ee in tx["exons"]:
        es2, ee2 = max(es, x0), min(ee, x1)
        if ee2 > es2:
            ax.add_patch(Rectangle((es2, y - h / 2), ee2 - es2, h, facecolor=color,
                                   edgecolor="none", zorder=4))
    # direction chevrons: white inside exons, transcript-color on introns
    if b > a:
        step = max(span_data * arrow_px, 1)
        ch = ">" if tx["strand"] == "+" else "<"
        n = int((b - a) / step)
        for i in range(1, n + 1):
            cx = a + step * i
            on_exon = any(es <= cx <= ee for es, ee in tx["exons"])
            ax.text(cx, y, ch, color=("white" if on_exon else color),
                    ha="center", va="center", fontsize=6.5, zorder=5,
                    fontweight="bold", clip_on=True)

def draw_ideogram(ax, chrom, clen, rs, re, cen=None):
    ax.set_xlim(-clen * 0.06, clen * 1.02)
    ax.set_ylim(0, 1)
    ax.axis("off")
    # chromosome body
    ax.add_patch(Rectangle((0, 0.30), clen, 0.40, facecolor="#e0e0e0",
                           edgecolor="#888888", lw=0.8, zorder=1))
    if cen:  # centromere notch
        ax.add_patch(Rectangle((cen - clen * 0.006, 0.30), clen * 0.012, 0.40,
                               facecolor="#b0392b", edgecolor="none", zorder=2))
    # red locus marker (widened for visibility)
    mw = max(re - rs, clen * 0.006)
    mc = (rs + re) / 2
    ax.add_patch(Rectangle((mc - mw / 2, 0.14), mw, 0.72, facecolor="none",
                           edgecolor="#d62728", lw=1.6, zorder=3))
    ax.text(-clen * 0.01, 0.5, chrom, ha="right", va="center", fontsize=8.5, fontweight="bold")


def draw_gene_context(ax, gc, chrom, x0, x1, feats):
    """Mini full-gene overview with a red box marking the current zoom window."""
    gs, ge = gc["region"]
    span = ge - gs
    ax.set_xlim(gs - span * 0.13, ge + span * 0.02)
    ax.axis("off")
    color = gc.get("color", "#555555")
    n = max(1, len(feats))
    ax.set_ylim(-0.7, n - 0.3)
    for i, f in enumerate(feats):
        y = n - 1 - i
        ax.plot([f["start"], f["end"]], [y, y], color=color, lw=0.7, zorder=2)
        for es, ee in f["exons"]:
            ax.add_patch(Rectangle((es, y - 0.28), max(ee - es, span * 0.0015), 0.56,
                                   facecolor=color, edgecolor="none", zorder=3))
    ax.add_patch(Rectangle((x0, -0.62), x1 - x0, n - 0.1, facecolor="none",
                           edgecolor="#d62728", lw=1.7, zorder=6))
    ax.text(gs - span * 0.02, (n - 1) / 2.0, gc.get("name", ""), ha="right", va="center",
            fontsize=8, fontweight="bold")


def plot_locus(cfg, out):
    chrom, x0, x1 = cfg["region"]
    span = x1 - x0
    ROW = 1.0            # one isoform row
    GAP = 0.6            # gap between tracks
    LBL = 0.0
    # ---- precompute track heights (in row-units) ----
    blocks = []
    for t in cfg["tracks"]:
        ty = t["type"]
        if ty == "coverage":
            cov_h = t.get("rows", 4)
            blocks.append((t, cov_h))
        elif ty == "isoforms":
            feats = parse_bed12(t["bed"], chrom, x0, x1)
            if t.get("strand"):
                feats = [f for f in feats if f["strand"] == t["strand"]]
            if t.get("jc_filter"):   # keep only isoforms on a specific junction chain
                feats = [f for f in feats if feat_jc(f) == t["jc_filter"]]
            if t.get("counts"):   # keep only isoforms with >= min_count reads assigned
                cnt = load_counts(t["counts"], t.get("count_col", 1))
                feats = [f for f in feats if cnt.get(f["name"], 0) >= t.get("min_count", 1)]
            sk = t.get("sort", "5p")
            if sk in ("5p", "3p"):
                def keyf(f, sk=sk):
                    if sk == "5p":
                        return f["end"] if f["strand"] == "-" else f["start"]
                    return f["start"] if f["strand"] == "-" else f["end"]
                feats.sort(key=keyf)
            t["_total"] = len(feats)
            mr = t.get("max_rows")
            if mr and len(feats) > mr:   # evenly sample to keep the end-range visible
                idx = sorted(set(round(i * (len(feats) - 1) / (mr - 1)) for i in range(mr)))
                feats = [feats[i] for i in idx]
            t["_feats"] = feats
            blocks.append((t, max(1, len(feats))))
        elif ty == "annotation":
            feats = parse_bed12(t["bed"], chrom, x0, x1)
            if t.get("strand"):
                feats = [f for f in feats if f["strand"] == t["strand"]]
            mr = t.get("max_rows")
            if mr and len(feats) > mr:
                feats = feats[:mr]
            t["_feats"] = feats
            blocks.append((t, max(1, len(feats))))
        elif ty in ("peaks", "region"):
            blocks.append((t, 1))
        elif ty == "sequence":
            blocks.append((t, 1.2))
    total = sum(h + GAP for _, h in blocks) + 1
    fig_h = max(3.0, total * cfg.get("row_scale", 0.26))   # row_scale up = taller/more vertical
    W = cfg.get("width", 13)
    ideo = cfg.get("ideogram")
    gc = cfg.get("gene_context")
    rs = cfg.get("row_scale", 0.26)
    if gc:
        gcf = parse_bed12(gc["bed"], chrom, gc["region"][0], gc["region"][1])
        if gc.get("strand"):
            gcf = [f for f in gcf if f["strand"] == gc["strand"]]
        gcf = gcf[:4]
        ntx = max(1, len(gcf))
        top_rows = ntx + 1.2
        fig = plt.figure(figsize=(W, (total + top_rows) * rs + 0.3))
        gsp = fig.add_gridspec(2, 1, height_ratios=[top_rows, total], hspace=0.07)
        ax = fig.add_subplot(gsp[1])
        draw_gene_context(fig.add_subplot(gsp[0]), gc, chrom, x0, x1, gcf)
    elif ideo:
        fig = plt.figure(figsize=(W, fig_h + 0.55))
        gs = fig.add_gridspec(2, 1, height_ratios=[0.5, fig_h], hspace=0.05)
        ax_ideo = fig.add_subplot(gs[0])
        ax = fig.add_subplot(gs[1])
        draw_ideogram(ax_ideo, chrom, ideo["chrom_length"], x0, x1, ideo.get("centromere"))
    else:
        fig, ax = plt.subplots(figsize=(W, fig_h))
    arrow_px = cfg.get("arrow_every", 0.045)

    y = total
    label_x = x0 - span * 0.005
    for t, h in blocks:
        ty = t["type"]
        top = y
        mid = y - h / 2
        color = t.get("color", "#333333")
        if ty == "coverage":
            bam = pysam.AlignmentFile(t["bam"])
            cov = np.sum(bam.count_coverage(chrom, x0, x1, quality_threshold=0), axis=0)
            xs = np.arange(x0, x1)
            cm = cov.max() if cov.max() > 0 else 1
            ax.fill_between(xs, y - h, y - h + (cov / cm) * h, color=color, lw=0, zorder=2)
            ax.text(x0 - span * 0.012, y - h / 2, t["name"], ha="right", va="center", fontsize=8)
            ax.text(x0 + span * 0.002, y - 0.25, f"[0–{int(cm)}]", fontsize=5.5, color="#666", va="top")
        elif ty in ("isoforms", "annotation"):
            feats = t["_feats"]
            n = max(1, len(feats))
            for i, tx in enumerate(feats):
                ry = top - 0.5 - i  # one row each, top-down
                draw_tx(ax, ry, 0.42, tx, color, x0, x1, arrow_px, span)
            ax.text(x0 - span * 0.012, top - n / 2, t["name"], ha="right", va="center",
                    fontsize=8.5, color=color if ty == "isoforms" else "#222",
                    fontweight="bold" if ty == "isoforms" else "normal")
            if ty == "isoforms" and cfg.get("show_n", True):
                tot = t.get("_total", n)
                lbl = f"n={tot}" + (f"\n({n} shown)" if n < tot else "")
                ax.text(x1 + span * 0.004, top - n / 2, lbl, ha="left", va="center",
                        fontsize=7, color=color)
        elif ty == "peaks":
            for pk in parse_bed_peaks(t["bed"], chrom, x0, x1):
                cx = (pk["start"] + pk["end"]) / 2
                ax.add_patch(Rectangle((pk["start"], mid - 0.32), max(pk["end"] - pk["start"], span * 0.001),
                                       0.64, facecolor=color, edgecolor="none", zorder=3))
                if pk["name"] and t.get("peak_labels"):
                    ax.text(cx, mid - 0.55, pk["name"], ha="center", va="top", fontsize=5.5, color=color)
            ax.text(x0 - span * 0.012, mid, t["name"], ha="right", va="center", fontsize=8, color=color)
        elif ty == "region":
            for (rs, re) in t["intervals"]:
                ax.add_patch(Rectangle((rs, mid - 0.3), re - rs, 0.6, facecolor=color,
                                       edgecolor="none", zorder=3, alpha=0.9))
                ax.text((rs + re) / 2, mid - 0.5, t.get("label", ""), ha="center", va="top",
                        fontsize=6, color="#8a6d00")
            ax.text(x0 - span * 0.012, mid, t["name"], ha="right", va="center", fontsize=8)
        y = top - h - GAP

    # vertical guide lines
    for gx, gc, gl in cfg.get("vlines", []):
        ax.axvline(gx, color=gc, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.set_xlim(x0, x1)
    ax.set_ylim(0, total + 0.5)
    ax.set_yticks([])
    for sp in ("top", "left", "right"):
        ax.spines[sp].set_visible(False)
    if cfg.get("hide_xaxis"):
        ax.set_xticks([])
        ax.spines["bottom"].set_visible(False)
    else:
        ax.tick_params(axis="x", labelsize=8)
        ax.ticklabel_format(axis="x", style="plain", useOffset=False)
        ax.set_xlabel(f"{chrom}:{x0:,}-{x1:,}", fontsize=9)
    if cfg.get("title"):
        ax.set_title(cfg["title"], fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print("wrote", out, f"({fig_h:.1f}in)")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    plot_locus(json.load(open(a.config)), a.out)
