#!/usr/bin/env python3
"""
end_signal_metaplot.py — Meta-profile of orthogonal signal at called isoform ends.

For each method, averages CAGE/dRNA bedGraph signal in a ±flank window
centred on every called TSS (5′) and TTS (3′).  A good method shows a
sharp peak at position 0; a poor/permissive method shows a flat or
off-centre profile.

Outputs (all in --output dir):
    metaplot_5prime.png   — per-method CAGE profiles at TSS, all on one axes
    metaplot_3prime.png   — per-method dRNA profiles at TTS, all on one axes
    metaplot_combined.png — 2-panel (5′ left, 3′ right) summary figure
    metaplot_metrics.tsv  — enrichment & FWHM per method per end

Usage:
    python end_signal_metaplot.py \\
        --bed label1:bed1.bed label2:bed2.bed ... \\
        --cage-plus  cage_plus.bg  --cage-minus  cage_minus.bg \\
        --qs-plus    qs_plus.bg    --qs-minus    qs_minus.bg \\
        --output     output_dir/   \\
        [--flank 500] [--bin-size 10] [--smooth 11]
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.signal import savgol_filter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
    from signal_utils import parse_isoforms, load_signal_tracks, BedGraphTrack
except ImportError:
    from evaluation.pub_style import apply_rc, style_ax, savefig, W1, W2, ModeStyler, legend_outside
    from evaluation.signal_utils import parse_isoforms, load_signal_tracks, BedGraphTrack

apply_rc()


# ── Boundary extraction ───────────────────────────────────────────────────────

def _get_boundaries(isos: List[dict]) -> Tuple[
    List[Tuple[str, int, str]],   # tss: (chrom, pos, strand)
    List[Tuple[str, int, str]],   # tts: (chrom, pos, strand)
]:
    tss, tts = [], []
    for iso in isos:
        ch, strand = iso["chrom"], iso["strand"]
        if strand == "+":
            tss.append((ch, iso["start"], strand))
            tts.append((ch, iso["end"],   strand))
        else:
            tss.append((ch, iso["end"],   strand))
            tts.append((ch, iso["start"], strand))
    return tss, tts


# ── Meta-profile builder ─────────────────────────────────────────────────────

def _compute_profile(
    boundaries: List[Tuple[str, int, str]],
    track_plus: BedGraphTrack,
    track_minus: BedGraphTrack,
    flank: int,
    bin_size: int,
) -> Optional[np.ndarray]:
    """Average signal profile centred on boundary positions.

    Minus-strand profiles are reversed so upstream is always left.
    Returns a 1-D array of length (2*flank // bin_size) or None.
    """
    n_bins = (2 * flank) // bin_size
    acc = np.zeros(n_bins, dtype=np.float64)
    n_valid = 0

    for chrom, pos, strand in boundaries:
        track = track_plus if strand == "+" else track_minus
        raw = track.window_values(chrom, pos - flank, pos + flank)
        if raw.sum() == 0:
            continue
        if strand == "-":
            raw = raw[::-1]
        # Bin-average
        usable = (len(raw) // bin_size) * bin_size
        binned = raw[:usable].reshape(-1, bin_size).mean(axis=1)
        if len(binned) == n_bins:
            acc += binned
            n_valid += 1

    return (acc / n_valid).astype(np.float32) if n_valid > 0 else None


def _profile_metrics(profile: np.ndarray, bin_size: int) -> dict:
    n = len(profile)
    c = n // 2
    hw = max(2, n // 40)
    centre = float(profile[c - hw: c + hw + 1].mean())
    flank_n = max(2, n // 5)
    flank_sig = float((profile[:flank_n].mean() + profile[-flank_n:].mean()) / 2)
    enrichment = centre / (flank_sig + 1e-12)
    peak = profile.max()
    above = np.where(profile >= peak / 2)[0]
    fwhm_bp = int(above[-1] - above[0] + 1) * bin_size if len(above) else 0
    return {"enrichment": round(enrichment, 3), "fwhm_bp": fwhm_bp,
            "centre_signal": round(centre, 6)}


def _smooth(arr: np.ndarray, window: int = 11, order: int = 3) -> np.ndarray:
    if len(arr) < window:
        return arr
    return savgol_filter(arr, window_length=window, polyorder=order)


# ── Plotting ─────────────────────────────────────────────────────────────────

def _plot_profiles(
    profiles: Dict[str, np.ndarray],
    flank: int,
    bin_size: int,
    smooth: int,
    ax: plt.Axes,
    end_label: str,
    signal_label: str,
    styler: ModeStyler,
) -> List:
    """Draw all method profiles onto ax. Returns legend handles."""
    n_bins = (2 * flank) // bin_size
    x = np.linspace(-flank, flank, n_bins, endpoint=False) + bin_size / 2
    handles = []
    for method, profile in profiles.items():
        sm = _smooth(profile, smooth)
        line, = ax.plot(x, sm, color=styler.color(method), linewidth=1.2, alpha=0.9)
        line.set_dashes(styler.dash(method))
        handles.append(styler.legend_handle(method, label=method, markersize=5))

    ax.axvline(0, color="0.3", linewidth=0.8, linestyle="--", alpha=0.6)
    ax.set_xlim(-flank, flank)
    ax.set_ylim(bottom=0)
    style_ax(ax,
             xlabel=f"Distance from {end_label} (bp)",
             ylabel=f"Mean {signal_label} signal")
    ax.set_axisbelow(True)
    return handles


def _save_single(
    profiles: Dict[str, np.ndarray],
    flank: int, bin_size: int, smooth: int,
    end_label: str, signal_label: str,
    styler: ModeStyler,
    output_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(W2, W2 * 0.65))
    handles = _plot_profiles(profiles, flank, bin_size, smooth,
                             ax, end_label, signal_label, styler)
    legend_outside(fig, handles=handles, loc="upper left",
                   bbox_to_anchor=(1.02, 1.0), ncol=1, fontsize=6)
    fig.tight_layout(pad=0.3)
    savefig(fig, output_path)


def plot_combined(
    tss_profiles: Dict[str, np.ndarray],
    tts_profiles: Dict[str, np.ndarray],
    flank: int, bin_size: int, smooth: int,
    output_dir: Path,
) -> None:
    methods = sorted(set(list(tss_profiles.keys()) + list(tts_profiles.keys())))
    styler = ModeStyler(methods)

    # Individual plots
    if tss_profiles:
        _save_single(tss_profiles, flank, bin_size, smooth,
                     "TSS", "CAGE", styler,
                     output_dir / "metaplot_5prime.png")
    if tts_profiles:
        _save_single(tts_profiles, flank, bin_size, smooth,
                     "TTS", "dRNA", styler,
                     output_dir / "metaplot_3prime.png")

    # Combined 2-panel
    if tss_profiles and tts_profiles:
        fig, (ax5, ax3) = plt.subplots(1, 2, figsize=(W2 * 2, W2 * 0.65))
        handles = _plot_profiles(tss_profiles, flank, bin_size, smooth,
                                 ax5, "TSS", "CAGE", styler)
        _plot_profiles(tts_profiles, flank, bin_size, smooth,
                       ax3, "TTS", "dRNA", styler)
        ax5.set_title("5′ end (CAGE)", fontsize=7)
        ax3.set_title("3′ end (dRNA)", fontsize=7)
        legend_outside(fig, handles=handles, loc="upper left",
                       bbox_to_anchor=(1.02, 1.0), ncol=1, fontsize=6)
        fig.tight_layout(pad=0.3)
        savefig(fig, output_dir / "metaplot_combined.png")


def write_metrics_tsv(
    tss_profiles: Dict[str, np.ndarray],
    tts_profiles: Dict[str, np.ndarray],
    bin_size: int,
    output_path: Path,
) -> None:
    rows = []
    for method in sorted(set(list(tss_profiles.keys()) + list(tts_profiles.keys()))):
        row = {"method": method}
        if method in tss_profiles:
            m = _profile_metrics(tss_profiles[method], bin_size)
            row["tss_enrichment"]    = m["enrichment"]
            row["tss_fwhm_bp"]       = m["fwhm_bp"]
            row["tss_centre_signal"] = m["centre_signal"]
        if method in tts_profiles:
            m = _profile_metrics(tts_profiles[method], bin_size)
            row["tts_enrichment"]    = m["enrichment"]
            row["tts_fwhm_bp"]       = m["fwhm_bp"]
            row["tts_centre_signal"] = m["centre_signal"]
        rows.append(row)

    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with open(output_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"Saved metaplot metrics to {output_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed",        nargs="+", required=True,
                   help="label:path pairs for BED12/GTF isoform files")
    p.add_argument("--cage-plus",  required=True)
    p.add_argument("--cage-minus", required=True)
    p.add_argument("--qs-plus",    required=True)
    p.add_argument("--qs-minus",   required=True)
    p.add_argument("--output",     required=True, help="Output directory")
    p.add_argument("--flank",      type=int, default=500,
                   help="Flank in bp each side of called end (default 500)")
    p.add_argument("--bin-size",   type=int, default=10,
                   help="Bin size in bp for averaging (default 10)")
    p.add_argument("--smooth",     type=int, default=11,
                   help="Savitzky-Golay smoothing window in bins (default 11)")
    p.add_argument("--verbose",    action="store_true")
    args = p.parse_args()

    # Parse BED files
    beds_by_method: Dict[str, List[dict]] = {}
    for entry in args.bed:
        if ":" not in entry:
            print(f"WARNING: skipping malformed entry '{entry}'", file=sys.stderr)
            continue
        label, path = entry.split(":", 1)
        if not Path(path).exists():
            print(f"WARNING: file not found: {path}", file=sys.stderr)
            continue
        isos = parse_isoforms(path)
        if isos:
            beds_by_method[label] = isos
            if args.verbose:
                print(f"  {label}: {len(isos)} isoforms")

    if not beds_by_method:
        print("No isoform data loaded — skipping", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print("Loading signal tracks...")
    cage_p, cage_m, qs_p, qs_m = load_signal_tracks(
        args.cage_plus, args.cage_minus, args.qs_plus, args.qs_minus,
    )

    tss_profiles: Dict[str, np.ndarray] = {}
    tts_profiles: Dict[str, np.ndarray] = {}

    for method, isos in beds_by_method.items():
        if args.verbose:
            print(f"  Computing profiles for {method}...")
        tss_boundaries, tts_boundaries = _get_boundaries(isos)

        p5 = _compute_profile(tss_boundaries, cage_p, cage_m,
                              args.flank, args.bin_size)
        p3 = _compute_profile(tts_boundaries, qs_p, qs_m,
                              args.flank, args.bin_size)
        if p5 is not None:
            tss_profiles[method] = p5
        if p3 is not None:
            tts_profiles[method] = p3

    if not tss_profiles and not tts_profiles:
        print("No profiles computed (no signal overlap?) — skipping", file=sys.stderr)
        sys.exit(1)

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    plot_combined(tss_profiles, tts_profiles, args.flank, args.bin_size,
                  args.smooth, out)
    write_metrics_tsv(tss_profiles, tts_profiles, args.bin_size,
                      out / "metaplot_metrics.tsv")
    print(f"Saved metaplot outputs to {args.output}")


if __name__ == "__main__":
    main()
