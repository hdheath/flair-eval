"""Publication-oriented signal hexbin plots.

The signal axes use a zero-safe transform:

    log10((TPM + 0.1) / 0.1)

This keeps true zero on the plot while giving useful resolution below 1 TPM.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

from pub_style import W1, style_ax


SIGNAL_AXIS_OFFSET = 0.1
DENSITY_CMAP = plt.cm.viridis
X_HIST_COLOR = "#0072B2"  # TTS / dRNA marginal
Y_HIST_COLOR = "#D55E00"  # TSS / CAGE marginal
ZERO_HIST_COLOR = "#CC2936"  # dedicated bar for true-zero TPM in marginals


def as_signal_array(values: Iterable[float]) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr) & (arr >= 0)]


def signal_to_axis(values: Iterable[float] | float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return np.log10((values + SIGNAL_AXIS_OFFSET) / SIGNAL_AXIS_OFFSET)


def raw_signal_upper(values: Iterable[float], floor: float = 10.0) -> float:
    arr = as_signal_array(values)
    if arr.size == 0:
        return float(floor)
    max_v = float(arr.max())
    if max_v <= 0:
        return float(floor)
    return float(max(floor, 10.0 ** np.ceil(np.log10(max_v))))


def format_count_label(value: float) -> str:
    value = float(value)
    if value >= 1_000_000:
        label = f"{value / 1_000_000:.1f}M"
    elif value >= 1000:
        label = f"{value / 1000:.1f}k"
    else:
        label = f"{int(value)}"
    return label.replace(".0M", "M").replace(".0k", "k")


def _format_tpm_tick_label(value: float) -> str:
    value = float(value)
    if abs(value) < 1e-12:
        return "0"
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:g}"


def _raw_tpm_ticks(raw_upper: float) -> np.ndarray:
    if raw_upper <= 12 and abs(raw_upper - round(raw_upper)) < 1e-9:
        return np.arange(0, int(round(raw_upper)) + 1, 1, dtype=float)
    candidate_ticks = np.array(
        [0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0],
        dtype=float,
    )
    ticks = candidate_ticks[candidate_ticks <= raw_upper * 1.000001]
    if ticks.size == 0 or ticks[0] != 0.0:
        ticks = np.insert(ticks, 0, 0.0)
    if ticks[-1] < raw_upper:
        ticks = np.append(ticks, float(raw_upper))
    return ticks


def _apply_raw_tpm_ticks(ax, axis: str, raw_upper: float) -> None:
    ticks_raw = _raw_tpm_ticks(raw_upper)
    ticks_axis = signal_to_axis(ticks_raw)
    emphasized_ticks = {0.0, 0.1, 0.5, 1.0, 10.0, 100.0}
    labels: list[str] = []
    for tick in ticks_raw:
        is_integer = abs(tick - round(tick)) < 1e-9
        if raw_upper > 12:
            is_endpoint = abs(tick - raw_upper) < 1e-9
            if any(abs(tick - keep) < 1e-9 for keep in emphasized_ticks) or is_endpoint:
                labels.append(_format_tpm_tick_label(tick))
            else:
                labels.append("")
        elif is_integer and (int(round(tick)) == 0 or int(round(tick)) % 2 == 0):
            labels.append(str(int(round(tick))))
        else:
            labels.append("")
    if axis == "x":
        ax.set_xticks(ticks_axis)
        ax.set_xticklabels(labels)
    elif axis == "y":
        ax.set_yticks(ticks_axis)
        ax.set_yticklabels(labels)
    else:
        raise ValueError(f"unknown axis: {axis!r}")


def _step_from_hist(edges: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    step_edges = np.repeat(edges, 2)[1:-1]
    step_values = np.repeat(values, 2)
    return step_edges, step_values


def _apply_count_limit_ticks(ax, axis: str, count_max: float, color: str) -> None:
    upper = max(float(count_max), 1.0)
    positions = [0.0, float(np.log10(upper + 1.0))]
    labels = ["0", format_count_label(upper)]
    if axis == "x":
        ax.set_xticks(positions)
        ax.set_xticklabels(labels)
        ax.xaxis.tick_top()
        ax.tick_params(
            axis="x", top=True, labeltop=True, bottom=False, labelbottom=False,
            direction="out", labelsize=5, pad=1, length=2, width=0.4, colors=color,
        )
    elif axis == "y":
        ax.set_yticks(positions)
        ax.set_yticklabels(labels)
        ax.tick_params(
            axis="y", left=True, labelleft=True, direction="out", labelsize=5,
            pad=1, length=2, width=0.4, colors=color,
        )
    else:
        raise ValueError(f"unknown axis: {axis!r}")


def signal_hexbin_with_marginals(
    x,
    y,
    *,
    fig=None,
    gs=None,
    xlabel="TTS signal (dRNA TPM)",
    ylabel="TSS signal (CAGE TPM)",
    title: Optional[str] = None,
    auto_range: bool = False,
    fixed_range: Optional[tuple[float, float]] = (100.0, 100.0),
    cmap=None,
    show_colorbar: bool = True,
    show_marginal_limits: bool = True,
    hex_gridsize: int = 46,
    color_vmax_quantile: float = 0.995,
    fig_size=None,
    hist_bins: int = 48,
):
    """Draw a zero-safe log hexbin signal scatter with marginal histograms."""
    if cmap is None:
        cmap = DENSITY_CMAP

    x_raw = np.asarray(x, dtype=float)
    y_raw = np.asarray(y, dtype=float)
    if x_raw.shape != y_raw.shape:
        raise ValueError("x and y must have the same shape")

    valid = np.isfinite(x_raw) & np.isfinite(y_raw) & (x_raw >= 0) & (y_raw >= 0)
    dropped = int(valid.size - np.count_nonzero(valid))
    x_raw = x_raw[valid]
    y_raw = y_raw[valid]

    if fixed_range is not None:
        xmax_raw, ymax_raw = float(fixed_range[0]), float(fixed_range[1])
    elif auto_range:
        xmax_raw = raw_signal_upper(x_raw)
        ymax_raw = raw_signal_upper(y_raw)
    else:
        xmax_raw = ymax_raw = 10.0

    xmax_raw = max(float(xmax_raw), 1.0)
    ymax_raw = max(float(ymax_raw), 1.0)
    xmax_axis = float(signal_to_axis(xmax_raw))
    ymax_axis = float(signal_to_axis(ymax_raw))
    x_clipped_high = int(np.count_nonzero(x_raw > xmax_raw))
    y_clipped_high = int(np.count_nonzero(y_raw > ymax_raw))
    visible = (x_raw <= xmax_raw) & (y_raw <= ymax_raw)
    n_visible = int(np.count_nonzero(visible))
    positive_pair = visible & (x_raw > 0.0) & (y_raw > 0.0)
    x_zero_y_positive = visible & (x_raw == 0.0) & (y_raw > 0.0)
    y_zero_x_positive = visible & (y_raw == 0.0) & (x_raw > 0.0)
    both_zero = visible & (x_raw == 0.0) & (y_raw == 0.0)
    x_pos = x_raw[positive_pair]
    y_pos = y_raw[positive_pair]
    n_positive_pair = int(np.count_nonzero(positive_pair))

    standalone = gs is None
    if standalone:
        if fig_size is None:
            fig_size = (W1, W1)
        fig = plt.figure(figsize=fig_size)
        outer = gridspec.GridSpec(
            2, 3, figure=fig,
            width_ratios=[5.0, 0.95, 0.28],
            height_ratios=[0.95, 5.0],
            wspace=0.07, hspace=0.06,
            left=0.16, right=0.94, bottom=0.14, top=0.93,
        )
    else:
        if fig is None:
            raise ValueError("must pass fig when supplying gs")
        outer = gs.subgridspec(
            2, 3,
            width_ratios=[5.0, 0.95, 0.28],
            height_ratios=[0.95, 5.0],
            wspace=0.07, hspace=0.06,
        )

    ax_top = fig.add_subplot(outer[0, 0])
    ax_sc = fig.add_subplot(outer[1, 0])
    ax_right = fig.add_subplot(outer[1, 1])
    cax = fig.add_subplot(outer[1, 2])
    fig.add_subplot(outer[0, 1]).axis("off")
    fig.add_subplot(outer[0, 2]).axis("off")

    hex_artist = None
    hex_count_max = 0
    hex_color_vmax = 1.0
    if n_positive_pair:
        hex_artist = ax_sc.hexbin(
            signal_to_axis(x_pos),
            signal_to_axis(y_pos),
            gridsize=hex_gridsize,
            extent=(0.0, xmax_axis, 0.0, ymax_axis),
            mincnt=1,
            cmap=cmap,
            linewidths=0.0,
            edgecolors="none",
        )
        counts = np.asarray(hex_artist.get_array(), dtype=float)
        if counts.size:
            hex_count_max = int(counts.max())
            robust = float(np.quantile(counts, color_vmax_quantile))
            hex_color_vmax = max(2.0, robust)
            if hex_count_max > 1 and robust <= 1.0:
                hex_color_vmax = float(hex_count_max)
            hex_artist.set_norm(LogNorm(vmin=1.0, vmax=hex_color_vmax))
            hex_artist.set_clim(1.0, hex_color_vmax)
    elif not n_visible:
        ax_sc.text(
            0.5, 0.5, "No points within range", transform=ax_sc.transAxes,
            ha="center", va="center", fontsize=7, color="#555555",
        )

    bins_x = np.linspace(0.0, xmax_axis, hist_bins)
    bins_y = np.linspace(0.0, ymax_axis, hist_bins)
    bin_w_x = xmax_axis / max(hist_bins - 1, 1)
    bin_w_y = ymax_axis / max(hist_bins - 1, 1)
    zero_rail_w_x = max(0.055 * xmax_axis, 2.2 * bin_w_x)
    zero_rail_w_y = max(0.055 * ymax_axis, 2.2 * bin_w_y)
    zero_rail_gap_x = max(0.015 * xmax_axis, 0.35 * bin_w_x)
    zero_rail_gap_y = max(0.015 * ymax_axis, 0.35 * bin_w_y)
    x_min_axis = -(zero_rail_w_x + zero_rail_gap_x)
    y_min_axis = -(zero_rail_w_y + zero_rail_gap_y)

    def _scaled_rail(values: np.ndarray, max_width: float) -> np.ndarray:
        if values.size == 0 or values.max() <= 0:
            return np.zeros_like(values, dtype=float)
        denom = float(np.log10(values.max() + 1.0))
        return max_width * np.log10(values + 1.0) / denom

    y_axis_for_x_zero = signal_to_axis(y_raw[x_zero_y_positive])
    x_zero_counts_by_y, _ = np.histogram(y_axis_for_x_zero, bins=bins_y)
    x_axis_for_y_zero = signal_to_axis(x_raw[y_zero_x_positive])
    y_zero_counts_by_x, _ = np.histogram(x_axis_for_y_zero, bins=bins_x)
    x_zero_y_positive_count = int(np.count_nonzero(x_zero_y_positive))
    y_zero_x_positive_count = int(np.count_nonzero(y_zero_x_positive))
    both_zero_count = int(np.count_nonzero(both_zero))

    if x_zero_y_positive_count or both_zero_count:
        ax_sc.axvspan(
            x_min_axis, -zero_rail_gap_x,
            color=ZERO_HIST_COLOR, alpha=0.055, linewidth=0, zorder=0,
        )
    if y_zero_x_positive_count or both_zero_count:
        ax_sc.axhspan(
            y_min_axis, -zero_rail_gap_y,
            color=ZERO_HIST_COLOR, alpha=0.055, linewidth=0, zorder=0,
        )
    if x_zero_y_positive_count:
        rail_widths = _scaled_rail(x_zero_counts_by_y.astype(float), zero_rail_w_x * 0.92)
        centers_y = 0.5 * (bins_y[:-1] + bins_y[1:])
        heights_y = np.diff(bins_y) * 0.86
        keep = rail_widths > 0
        ax_sc.barh(
            centers_y[keep],
            rail_widths[keep],
            left=-zero_rail_gap_x - rail_widths[keep],
            height=heights_y[keep],
            color=ZERO_HIST_COLOR,
            alpha=0.70,
            linewidth=0,
            align="center",
            zorder=2.2,
        )
    if y_zero_x_positive_count:
        rail_heights = _scaled_rail(y_zero_counts_by_x.astype(float), zero_rail_w_y * 0.92)
        centers_x = 0.5 * (bins_x[:-1] + bins_x[1:])
        widths_x = np.diff(bins_x) * 0.86
        keep = rail_heights > 0
        ax_sc.bar(
            centers_x[keep],
            rail_heights[keep],
            bottom=-zero_rail_gap_y - rail_heights[keep],
            width=widths_x[keep],
            color=ZERO_HIST_COLOR,
            alpha=0.70,
            linewidth=0,
            align="center",
            zorder=2.2,
        )
    if both_zero_count:
        ax_sc.add_patch(Rectangle(
            (x_min_axis, y_min_axis),
            zero_rail_w_x,
            zero_rail_w_y,
            facecolor=ZERO_HIST_COLOR,
            edgecolor="none",
            alpha=0.80,
            zorder=2.4,
        ))

    ax_sc.axvline(0.0, color="#555555", linewidth=0.45, alpha=0.70, zorder=2.6)
    ax_sc.axhline(0.0, color="#555555", linewidth=0.45, alpha=0.70, zorder=2.6)
    ax_sc.set_xlim(x_min_axis, xmax_axis)
    ax_sc.set_ylim(y_min_axis, ymax_axis)
    _apply_raw_tpm_ticks(ax_sc, "x", xmax_raw)
    _apply_raw_tpm_ticks(ax_sc, "y", ymax_raw)
    style_ax(ax_sc, xlabel=xlabel, ylabel=ylabel)
    ax_sc.tick_params(axis="both", direction="out", labelsize=7, pad=2)

    if show_colorbar and hex_artist is not None:
        cbar = fig.colorbar(hex_artist, cax=cax)
        cbar.ax.tick_params(direction="out", labelsize=5, length=2, width=0.4, pad=1)
        cbar.set_label("Count", fontsize=6)
    else:
        cax.axis("off")

    # True-zero isoforms are rendered as a dedicated bar (see below).  Count
    # them separately and exclude exact zeros from the regular histogram so
    # they don't contaminate the leftmost positive bin.
    x_zero_count = int(np.sum(x_raw == 0.0))
    y_zero_count = int(np.sum(y_raw == 0.0))
    x_hist_axis = signal_to_axis(x_raw[(x_raw > 0.0) & (x_raw <= xmax_raw)])
    y_hist_axis = signal_to_axis(y_raw[(y_raw > 0.0) & (y_raw <= ymax_raw)])
    counts_x, edges_x = np.histogram(x_hist_axis, bins=bins_x)
    counts_y, edges_y = np.histogram(y_hist_axis, bins=bins_y)
    hist_x = np.log10(counts_x + 1.0)
    hist_y = np.log10(counts_y + 1.0)

    # Zero-bar geometry: a narrow dedicated bar to the left of axis 0 (top
    # marginal) and below axis 0 (right marginal), separated by a small gap
    # from the regular histogram.  Heights use the same log10(count + 1)
    # mapping as the regular bars so they're directly comparable.
    zero_bar_w_x = 0.40 * bin_w_x
    zero_bar_w_y = 0.40 * bin_w_y
    zero_bar_gap_x = 0.20 * bin_w_x
    zero_bar_gap_y = 0.20 * bin_w_y
    zero_bar_center_x = -(zero_bar_gap_x + 0.5 * zero_bar_w_x)
    zero_bar_center_y = -(zero_bar_gap_y + 0.5 * zero_bar_w_y)
    zero_h_x = float(np.log10(x_zero_count + 1.0))
    zero_h_y = float(np.log10(y_zero_count + 1.0))

    x_step, top_step = _step_from_hist(edges_x, hist_x)
    ax_top.fill_between(x_step, 0.0, top_step, color=X_HIST_COLOR, alpha=0.26)
    ax_top.plot(x_step, top_step, color=X_HIST_COLOR, alpha=0.80, linewidth=0.45)
    if x_zero_count > 0:
        ax_top.bar(
            zero_bar_center_x, zero_h_x, width=zero_bar_w_x,
            color=ZERO_HIST_COLOR, alpha=0.85, linewidth=0.0, align="center",
        )
    # Decouple top x-axis from scatter to make room for the zero bar on the left.
    ax_top.set_xlim(x_min_axis, xmax_axis)
    y_top_max = max(float(hist_x.max()) if hist_x.size else 0.0, zero_h_x)
    ax_top.set_ylim(0.0, max(1.0, y_top_max * 1.08))
    ax_top.tick_params(axis="x", bottom=False, labelbottom=False, direction="out")
    if show_marginal_limits:
        _apply_count_limit_ticks(
            ax_top, "y",
            float(max(counts_x.max() if counts_x.size else 0, x_zero_count)),
            X_HIST_COLOR,
        )
    else:
        ax_top.tick_params(axis="y", left=False, labelleft=False, direction="out")
    ax_top.spines[["top", "right", "bottom", "left"]].set_visible(False)
    if title:
        ax_top.set_title(title, fontsize=8, pad=2)

    y_step, right_step = _step_from_hist(edges_y, hist_y)
    ax_right.fill_betweenx(y_step, 0.0, right_step, color=Y_HIST_COLOR, alpha=0.26)
    ax_right.plot(right_step, y_step, color=Y_HIST_COLOR, alpha=0.80, linewidth=0.45)
    if y_zero_count > 0:
        ax_right.barh(
            zero_bar_center_y, zero_h_y, height=zero_bar_w_y,
            color=ZERO_HIST_COLOR, alpha=0.85, linewidth=0.0, align="center",
        )
    # Decouple right y-axis from scatter to make room for the zero bar at the bottom.
    ax_right.set_ylim(y_min_axis, ymax_axis)
    x_right_max = max(float(hist_y.max()) if hist_y.size else 0.0, zero_h_y)
    ax_right.set_xlim(0.0, max(1.0, x_right_max * 1.08))
    ax_right.tick_params(axis="y", left=False, labelleft=False, direction="out")
    if show_marginal_limits:
        _apply_count_limit_ticks(
            ax_right, "x",
            float(max(counts_y.max() if counts_y.size else 0, y_zero_count)),
            Y_HIST_COLOR,
        )
    else:
        ax_right.tick_params(axis="x", bottom=False, labelbottom=False, direction="out")
    ax_right.spines[["top", "right", "left", "bottom"]].set_visible(False)

    return fig, {
        "ax_sc": ax_sc,
        "ax_top": ax_top,
        "ax_right": ax_right,
        "cax": cax,
        "xmax_raw": xmax_raw,
        "ymax_raw": ymax_raw,
        "n": int(len(x_raw)),
        "dropped": dropped,
        "n_visible": n_visible,
        "n_positive_pair": n_positive_pair,
        "x_clipped_high": x_clipped_high,
        "y_clipped_high": y_clipped_high,
        "x_hist_max": int(counts_x.max()) if len(counts_x) else 0,
        "y_hist_max": int(counts_y.max()) if len(counts_y) else 0,
        "x_zero_count": x_zero_count,
        "y_zero_count": y_zero_count,
        "x_zero_y_positive_count": x_zero_y_positive_count,
        "y_zero_x_positive_count": y_zero_x_positive_count,
        "both_zero_count": both_zero_count,
        "hex_count_max": hex_count_max,
        "hex_color_vmax": hex_color_vmax,
        "transform": "log10((TPM + 0.1) / 0.1)",
        "density_norm": "positive-positive hexbin-log-count; exact-zero axes as red rails",
    }


def write_signal_summary(path: str | Path, rows: list[dict]) -> None:
    """Write per-panel hexbin metadata as a TSV."""
    if not rows:
        return
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(rows[0].keys())
    with open(path, "w") as out:
        out.write("\t".join(keys) + "\n")
        for row in rows:
            out.write("\t".join(str(row.get(k, "")) for k in keys) + "\n")
