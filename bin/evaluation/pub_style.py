"""
Nature Portfolio–compliant matplotlib style for the flair-eval pipeline.

All plotting scripts should import from this module rather than setting
their own rcParams.  See ``signal_utils.py`` for shared data helpers.

Enforces:
  • Arial / Helvetica (falls back to DejaVu Sans)
  • 8 pt axis labels, 7 pt tick labels (Nature minimum 5 pt)
  • 0.5 pt axes & tick-mark widths, 0.8 pt data lines
  • Inward-facing ticks, no top/right spines
  • Okabe-Ito colorblind-safe categorical palette
  • 300 DPI PNG + SVG dual output
  • No background fill, no grid, no legend frames

Guidelines: https://www.nature.com/documents/nature-final-artwork.pdf
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional, Sequence, Tuple, Union

import matplotlib as mpl
import matplotlib.pyplot as plt


# ── Dimensions (Nature Portfolio) ───────────────────────────────────────────
# Single-column 89 mm, double-column 183 mm, max height 247 mm.

_MM = 1.0 / 25.4                       # mm → inches conversion factor
W1   = 89  * _MM                        # single-column width  (≈3.50 in)
W2   = 183 * _MM                        # double-column width  (≈7.20 in)
HMAX = 247 * _MM                        # max figure height    (≈9.72 in)


# ── Colorblind-safe palette (Okabe-Ito + extended) ──────────────────────────

PALETTE = [
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#999999",  # grey
]

# Assembler colours (kept distinct & colourblind-safe)
ASSEMBLER_COLORS = {
    "flair":    "#0072B2",  # blue
    "bambu":    "#E69F00",  # orange
    "isoquant": "#009E73",  # bluish green
    "unknown":  "#999999",
}

# Sample/dataset palette (cycles through PALETTE)
SAMPLE_PALETTE = PALETTE

# Library-type marker shapes (unchanged, but provided here for consistency)
LIBRARY_SHAPES = {
    "pacbio_cDNA": "o",
    "ont_cDNA":    "^",
    "ont_dRNA":    "s",
    "unknown":     "D",
}

# Assembler initials for labels
ASSEMBLER_INITIALS = {
    "flair":    "FL",
    "bambu":    "BU",
    "isoquant": "IS",
    "unknown":  "??",
}

# Mode colours (used in concordance / signal-support scripts)
MODE_COLORS = {
    "default":                     "#0072B2",
    "density-asymmetric":          "#009E73",
    "density-asymmetric-softclip": "#56B4E9",
    "density-plain":               "#E69F00",
    "density-strict":              "#D55E00",
    "k-means":                     "#CC79A7",
    "more-ends":                   "#F0E442",
    "bambu_default":               "#E69F00",
    "isoquant_pacbio":             "#009E73",
    # Trust-ends modes
    "baseline":                    "#332288",  # indigo (distinct from default)
    "trust-ends":                  "#D55E00",
    "trust-ends-sensitive":        "#E69F00",
    "asymmetric-3p-only":          "#009E73",
    "refine-tts-only":             "#CC79A7",
    # TED modes — extended palette (Paul Tol vibrant + muted)
    "ted-default":                 "#0072B2",  # blue
    "ted-2d":                      "#E69F00",  # orange
    "ted-1d2d":                    "#999999",  # grey
    "ted-2d-cluster":              "#E69F00",  # orange
    "ted-2d-reassign-refext":      "#CC6677",  # rose
    "ted-allow-single":            "#44AA99",  # teal
    "ted-epsilon-cluster":         "#AA4499",  # purple
    "ted-high-annot-weight":       "#882255",  # wine
    "ted-high-depth-weight":       "#009E73",  # bluish green
    "ted-keep-noise":              "#D55E00",  # vermillion
    "ted-leaf-robust-strict":      "#56B4E9",  # sky blue
    "ted-leaf-selection":          "#DDCC77",  # sand
    "ted-lenient-threshold":       "#999933",  # olive
    "ted-minmax-norm":             "#117733",  # forest green
    "ted-reassign-noise":          "#88CCEE",  # cyan
    "ted-ref-extension":           "#CC79A7",  # reddish purple
    "ted-robust-norm":             "#661100",  # brown
    "ted-strict-threshold":        "#F0E442",  # yellow
}

# Reason colours (peak recovery; already fairly accessible, minor tweaks)
REASON_COLORS = {
    "recovered":           "#009E73",
    "no_reads":            "#999999",
    "single_exon_only":    "#56B4E9",
    "near_miss":           "#F0E442",
    "reads_unassigned":    "#D55E00",
    "reads_redirected":    "#CC79A7",
    "proximal_apa":        "#E69F00",
    "trailing_truncation": "#0072B2",
    "no_isoform_model":    "#882255",
    "alignment_filtered":  "#AA4499",
    "end_absorbed":        "#44AA99",
    "end_spread":          "#DDCC77",
    "low_signal":          "#88CCEE",
    "other_missed":        "#999999",
}


# ── Nature Portfolio rcParams (applied once at import time) ─────────────────

_RC_OVERRIDES = {
    # Font: Arial / Helvetica, Nature 8 pt labels / 7 pt ticks
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Arial", "Helvetica Neue", "Helvetica",
                           "DejaVu Sans"],
    "font.size":          8,
    "axes.titlesize":     8,
    "axes.labelsize":     8,
    "xtick.labelsize":    7,
    "ytick.labelsize":    7,
    "legend.fontsize":    7,
    # Line weights: 0.5 pt axes/ticks, 0.8 pt data
    "axes.linewidth":     0.5,
    "xtick.major.width":  0.5,
    "ytick.major.width":  0.5,
    "xtick.minor.width":  0.25,
    "ytick.minor.width":  0.25,
    "xtick.major.size":   3,
    "ytick.major.size":   3,
    "xtick.minor.size":   1.5,
    "ytick.minor.size":   1.5,
    "lines.linewidth":    0.8,
    "lines.markersize":   3,
    "patch.linewidth":    0.35,
    # Ticks: inward, no top/right spines
    "xtick.direction":    "in",
    "ytick.direction":    "in",
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    # Background: clean white, no grid, no legend frame
    "figure.facecolor":   "white",
    "axes.facecolor":     "white",
    "savefig.facecolor":  "white",
    "axes.grid":          False,
    "legend.frameon":     False,
    # Output quality
    "figure.dpi":         300,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
    "savefig.pad_inches": 0.02,
    # Embed fonts for vector formats
    "svg.fonttype":       "none",
    "pdf.fonttype":       42,
    "ps.fonttype":        42,
    "mathtext.default":   "regular",
}


# ── Plotly layout template (for scripts using Plotly instead of mpl) ────────

PLOTLY_LAYOUT = dict(
    font=dict(family="Arial, Helvetica, sans-serif", size=8),
    title_font_size=8,
    plot_bgcolor="white",
    paper_bgcolor="white",
    margin=dict(l=50, r=50, t=50, b=40),
)


class ModeStyler:
    """Assign unique (color, marker, linestyle) to an arbitrary set of modes.

    Modes are grouped by their family prefix (everything before the first
    hyphen or underscore).  Within a family, colors are spread across a
    gradient so family members share a hue but stay distinct.  Markers and
    linestyles cycle independently, ensuring every mode is unique even if
    colors overlap.

    Usage
    -----
    >>> styler = ModeStyler(list_of_mode_names)
    >>> color  = styler.color(mode)
    >>> marker = styler.marker(mode)
    >>> dash   = styler.dash(mode)   # (on, off, …) tuple for set_dashes()
    >>> handle = styler.legend_handle(mode, label=mode)
    """

    _MARKERS = ["o", "s", "^", "v", "D", "P", "X", "*", "h", "p", "<", ">", "H"]

    _DASHES = [
        (1, 0),          # solid
        (4, 2),          # dashed
        (1, 2),          # dotted
        (6, 2, 2, 2),    # dash-dot
        (3, 1, 1, 1),    # dash-dot-dot
        (8, 2),          # long dash
        (2, 2, 6, 2),    # dot-long-dash
    ]

    # One color family per common prefix.  More-specific prefixes must be
    # listed before shorter ones so "TED-asym" matches before "TED".
    _FAMILY_PALETTE = {
        "TED-asym":  ["#E69F00", "#D55E00", "#CC79A7", "#882255", "#F0A050"],
        "TED":       ["#88CCEE", "#0072B2", "#003A6B", "#5BA3CC", "#56B4E9"],
        "FLAIR":     ["#44AA99", "#009E73", "#005740"],
        "isoquant":  ["#AAAAAA", "#555555"],
        "bambu":     ["#FFCC80", "#E69F00", "#8B5E00"],
        "isoseq":    ["#DDB8E8", "#CC79A7", "#6B2D54"],
        "flames":    ["#F7E988", "#F0E442", "#A89900"],
        "stringtie": ["#FFAA88", "#D55E00", "#8B3A00"],
        "_other":    PALETTE,
    }

    def __init__(self, modes: list):
        seen = set()
        self._modes = [m for m in modes if not (m in seen or seen.add(m))]
        self._color  = {}
        self._marker = {}
        self._dash   = {}
        self._build()

    def _family(self, mode: str) -> str:
        mode_lower = mode.lower()
        for key in sorted(self._FAMILY_PALETTE.keys(), key=len, reverse=True):
            if key == "_other":
                continue
            if mode_lower.startswith(key.lower()):
                return key
        return "_other"

    def _build(self):
        from collections import defaultdict
        from matplotlib.colors import LinearSegmentedColormap
        families: dict = defaultdict(list)
        for mode in self._modes:
            families[self._family(mode)].append(mode)

        for fam, fam_modes in families.items():
            pal = self._FAMILY_PALETTE.get(fam, PALETTE)
            n = len(fam_modes)
            if n == 1:
                colors = [pal[len(pal) // 2]]
            else:
                cmap = LinearSegmentedColormap.from_list(f"_ms_{fam}", pal, N=max(n, 2))
                colors = [
                    "#{:02x}{:02x}{:02x}".format(int(r*255), int(g*255), int(b*255))
                    for r, g, b, *_ in (cmap(i / (n-1)) for i in range(n))
                ]
            for i, mode in enumerate(fam_modes):
                self._color[mode] = colors[i]

        for i, mode in enumerate(self._modes):
            self._marker[mode] = self._MARKERS[i % len(self._MARKERS)]
            self._dash[mode]   = self._DASHES[i % len(self._DASHES)]

    def color(self, mode: str) -> str:
        return self._color.get(mode, PALETTE[hash(mode) % len(PALETTE)])

    def marker(self, mode: str) -> str:
        return self._marker.get(mode, "o")

    def dash(self, mode: str) -> tuple:
        return self._dash.get(mode, (1, 0))

    def legend_handle(self, mode: str, label: str = None, markersize: float = 6):
        """Return a Line2D showing both line style and marker for use in legend."""
        from matplotlib.lines import Line2D
        dash = self.dash(mode)
        ls = "-" if dash == (1, 0) else (0, dash)
        return Line2D(
            [], [],
            color=self.color(mode),
            marker=self.marker(mode),
            markersize=markersize,
            markerfacecolor=self.color(mode),
            markeredgecolor="white",
            markeredgewidth=0.4,
            linestyle=ls,
            linewidth=1.2,
            label=label if label is not None else mode,
        )

    def all_modes(self):
        return list(self._modes)


def apply_rc() -> None:
    """Apply Nature-compliant rcParams globally (idempotent)."""
    mpl.rcParams.update(_RC_OVERRIDES)


# Auto-apply at import so all downstream code inherits the settings.
apply_rc()


# ── Per-Axes helper ─────────────────────────────────────────────────────────

def style_ax(
    ax: plt.Axes,
    ylabel: Optional[str] = None,
    xlabel: Optional[str] = None,
    title: Optional[str] = None,
    faint_y_grid: bool = False,
) -> None:
    """Apply Nature-quality styling to a single Axes.

    Parameters
    ----------
    ax : matplotlib Axes
    ylabel, xlabel, title : optional override text
    faint_y_grid : if True, add faint horizontal grid behind data
    """
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="in", width=0.5, length=3)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    if title:
        ax.set_title(title, fontweight="normal")
    if faint_y_grid:
        ax.yaxis.grid(True, alpha=0.15, linewidth=0.5, color="#888888")
        ax.set_axisbelow(True)


def legend_outside(
    fig_or_ax,
    handles=None,
    labels=None,
    loc: str = "upper left",
    bbox_to_anchor=(1.02, 1.0),
    ncol: int = 1,
    **kw,
):
    """Place a frameless legend outside the plot area."""
    kw.setdefault("frameon", False)
    kw.setdefault("fontsize", 7)
    kw.setdefault("borderaxespad", 0)
    return fig_or_ax.legend(
        handles=handles,
        labels=labels,
        loc=loc,
        bbox_to_anchor=bbox_to_anchor,
        ncol=ncol,
        **kw,
    )


def savefig(
    fig: plt.Figure,
    path: Union[str, Path],
    dpi: int = 300,
    formats: Sequence[str] = ("png", "svg"),
    close: bool = True,
    **kw,
) -> None:
    """Save figure as PNG + SVG (Nature requires vector), then close.

    Parameters
    ----------
    fig : matplotlib Figure
    path : output path (extension is replaced per format)
    dpi : raster resolution (default 300, Nature minimum)
    formats : iterable of format strings; default ``("png", "svg")``
    close : whether to close the figure after saving
    **kw : forwarded to ``fig.savefig``
    """
    kw.setdefault("bbox_inches", "tight")
    kw.setdefault("facecolor", "white")
    path = Path(path)
    stem = path.with_suffix("")          # strip original extension
    for fmt in formats:
        out = stem.with_suffix(f".{fmt}")
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, format=fmt, dpi=dpi, **kw)
    if close:
        plt.close(fig)
