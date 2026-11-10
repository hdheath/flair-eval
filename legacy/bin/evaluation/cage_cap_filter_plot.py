#!/usr/bin/env python3
"""cage_cap_filter_plot.py — visualize the CAGE cap-signature (uG) filter.

The A549 CAGE truth peaks (ReCappable-seq, 2-rep reproducible) carry a per-peak
untemplated-G fraction (uG) = the m7G-cap signature.  cap_filter.py keeps a peak
if it is annotation-trusted (class in --exempt, default near_TSS) OR its uG >=
--threshold (default 0.45).  This script shows that filter:

  Panel A: uG distribution per class, with the threshold line — near_TSS peaks
           are cap-rich; novel peaks split into capped (real TSS) vs uncapped
           (recapping / internal artifacts).
  Panel B: kept vs dropped per class (exempt classes kept regardless).

Input is the uG-scored BED (name = consensus_X|class|uGY.YY).
"""
import argparse
from collections import defaultdict, Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CLASS_COLOR = {"near_TSS": "#1d9e75", "novel_intragenic": "#d4537e",
               "splice_acceptor": "#b19a00"}


def parse_ug_bed(path):
    """Return list of (class, uG_or_None) per peak."""
    out = []
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            parts = f[3].split("|")
            cls = parts[1] if len(parts) > 2 else "NA"
            ug = None
            for tok in parts:
                if tok.startswith("uG"):
                    try:
                        ug = float(tok[2:])
                    except ValueError:
                        ug = None
            out.append((cls, ug))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ug-bed",
                    default="/private/home/hdheath/projects/ortho-peaks/cage/"
                            "a549_cage_results/a549_cage_stringent_uG.bed")
    ap.add_argument("--threshold", type=float, default=0.45)
    ap.add_argument("--exempt", default="near_TSS",
                    help="comma list of annotation-trusted classes kept regardless of uG")
    ap.add_argument("--outpng", required=True)
    ap.add_argument("--title", default="A549 CAGE cap-signature (uG) filter")
    a = ap.parse_args()
    exempt = set(c for c in a.exempt.split(",") if c)

    peaks = parse_ug_bed(a.ug_bed)
    by_class = defaultdict(list)
    for cls, ug in peaks:
        by_class[cls].append(ug)
    # order classes by size (largest first), known colors first
    classes = sorted(by_class, key=lambda c: -len(by_class[c]))

    kept = Counter()
    dropped = Counter()
    tsv_rows = []
    for cls in classes:
        ugs = by_class[cls]
        for ug in ugs:
            keep = (cls in exempt) or (ug is not None and ug >= a.threshold)
            (kept if keep else dropped)[cls] += 1
        tsv_rows.append((cls, len(ugs), kept[cls], dropped[cls],
                         "exempt" if cls in exempt else f"uG>={a.threshold}"))

    fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))

    # ── Panel A: uG distribution per class + threshold ──
    bins = np.linspace(0, 1, 41)
    for cls in classes:
        vals = [u for u in by_class[cls] if u is not None]
        if not vals:
            continue
        ax[0].hist(vals, bins=bins, histtype="step", linewidth=2,
                   color=CLASS_COLOR.get(cls, "#888888"),
                   label=f"{cls} (n={len(by_class[cls])})")
    ax[0].axvline(a.threshold, color="black", ls="--", lw=1.4)
    ax[0].text(a.threshold + 0.01, ax[0].get_ylim()[1] * 0.92,
               f"keep ≥ {a.threshold}", fontsize=9, va="top")
    ax[0].set_xlabel("untemplated-G fraction (uG = cap signature)", fontsize=10)
    ax[0].set_ylabel("# CAGE peaks", fontsize=10)
    ax[0].set_title("uG distribution by peak class", fontsize=11)
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.25)
    ax[0].spines["top"].set_visible(False); ax[0].spines["right"].set_visible(False)

    # ── Panel B: kept vs dropped per class ──
    x = np.arange(len(classes))
    k = [kept[c] for c in classes]
    d = [dropped[c] for c in classes]
    ax[1].bar(x, k, color="#1d9e75", label="kept (cap-confirmed / exempt)")
    ax[1].bar(x, d, bottom=k, color="#cccccc", label="dropped (no cap signal)")
    for i, c in enumerate(classes):
        tot = k[i] + d[i]
        note = "exempt" if c in exempt else f"−{d[i]} ({100*d[i]/tot:.0f}%)"
        ax[1].text(i, tot, f"n={tot}\n{note}", ha="center", va="bottom", fontsize=8)
    ax[1].set_xticks(x); ax[1].set_xticklabels(classes, rotation=15, ha="right", fontsize=9)
    ax[1].set_ylabel("# CAGE peaks", fontsize=10)
    ax[1].set_title(f"kept vs dropped (total kept {sum(k)} / {sum(k)+sum(d)})", fontsize=11)
    ax[1].legend(fontsize=8, loc="upper right")
    ax[1].grid(axis="y", alpha=0.25)
    ax[1].spines["top"].set_visible(False); ax[1].spines["right"].set_visible(False)

    fig.suptitle(a.title + f"  —  exempt={sorted(exempt)}, uG≥{a.threshold}", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(a.outpng, dpi=200, bbox_inches="tight")
    print("wrote", a.outpng)

    tsv = a.outpng.rsplit(".", 1)[0] + ".tsv"
    with open(tsv, "w") as fh:
        fh.write("class\tn_total\tn_kept\tn_dropped\trule\n")
        for r in tsv_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print("wrote", tsv)


if __name__ == "__main__":
    main()
