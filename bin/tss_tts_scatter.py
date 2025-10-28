"""TSS/TTS scatter plotting utilities for TED QC."""
from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter


def _load_coords_bed(path: Path) -> tuple[list[str], np.ndarray, np.ndarray, list[str], list[str]]:
    """Load BED6 and compute TSS/TTS.

    Returns lists/arrays for chromosome, TTS (x), TSS (y), transcript ID and strand.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        comment="#",
        usecols=[0, 1, 2, 3, 5],
        names=["chrom", "start", "end", "transcript_id", "strand"],
        dtype={"chrom": str, "start": np.int64, "end": np.int64, "transcript_id": str, "strand": str},
    )
    tss = np.where(df["strand"] == "+", df["start"].to_numpy(), df["end"].to_numpy())
    tts = np.where(df["strand"] == "+", df["end"].to_numpy(), df["start"].to_numpy())
    return (
        df["chrom"].tolist(),
        tts.astype(int),
        tss.astype(int),
        df["transcript_id"].tolist(),
        df["strand"].tolist(),
    )


def _interval_overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """Return ``True`` if two intervals overlap."""
    return not (a_end < b_start or a_start > b_end)


def _parse_joint_window(jw_str: str) -> tuple[int, int, int, int]:
    parts = jw_str.split(":", 1)[1].split("|")
    tss_min, tss_max = map(int, parts[0].split("-", 1))
    tts_min, tts_max = map(int, parts[1].split("-", 1))
    return tss_min, tss_max, tts_min, tts_max


def _load_mapping(path: Optional[Path], chrom: str, region_start: Optional[int], region_end: Optional[int]):
    if path is None or not path.exists():
        return []
    df = pd.read_csv(path, sep="\t", dtype=str)
    required = {"most_freq_coord", "all_tx_ids", "joint_window", "strand"}
    if not required.issubset(df.columns):
        raise ValueError(
            f"{path} must have columns: most_freq_coord, all_tx_ids, joint_window, strand"
        )
    mf_split = df["most_freq_coord"].str.split(":", n=1, expand=True)
    df["map_chrom"] = mf_split[0]
    df = df[df["map_chrom"] == chrom].copy()
    if df.empty:
        return []

    mapping_list = []
    for idx, row in df.reset_index(drop=True).iterrows():
        try:
            tss_min, tss_max, tts_min, tts_max = _parse_joint_window(row["joint_window"])
        except Exception:
            continue
        if region_start is not None and region_end is not None:
            if not (
                _interval_overlaps(tss_min, tss_max, region_start, region_end)
                or _interval_overlaps(tts_min, tts_max, region_start, region_end)
            ):
                continue
        all_ids = row["all_tx_ids"].split(",")
        tx_set = {tx.strip() for tx in all_ids if tx.strip()}
        strand = row["strand"].strip()
        mapping_list.append(
            {
                "cluster_idx": idx,
                "tss_min": tss_min,
                "tss_max": tss_max,
                "tts_min": tts_min,
                "tts_max": tts_max,
                "strand": strand,
                "tx_set": tx_set,
            }
        )
    return mapping_list


def _compute_limits(x: np.ndarray, y: np.ndarray) -> tuple[int, int, int, int]:
    if len(x) == 0:
        return 0, 1, 0, 1
    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()
    xpad = max(int((xmax - xmin) * 0.02), 1)
    ypad = max(int((ymax - ymin) * 0.02), 1)
    return xmin - xpad, xmax + xpad, ymin - ypad, ymax + ypad


def _compute_hist(arr: np.ndarray, nbins: int, data_min: int, data_max: int):
    if len(arr) == 0:
        edges = np.linspace(data_min, data_max, nbins + 1)
        return np.zeros(nbins, dtype=int), edges
    counts, edges = np.histogram(arr, bins=nbins, range=(data_min, data_max))
    return counts, edges


def _coord_formatter(scale: int):
    def _formatter(x, pos):
        return f"{x/scale:.3g}"

    return FuncFormatter(_formatter)


def generate(
    coords: Path,
    chrom: str,
    output: Path,
    region: Optional[Tuple[int, int]] = None,
    mapping: Optional[Path] = None,
) -> None:
    """Create TTS vs TSS scatter plots per strand.

    Parameters
    ----------
    coords : Path
        BED6 file with collapsed isoform coordinates.
    chrom : str
        Chromosome to plot.
    output : Path
        Output PNG filename (base name used for each strand).
    region : tuple[int, int], optional
        Limit plotting to ``(start, end)`` within ``chrom``.
    mapping : Path, optional
        Optional mapping TSV for cluster highlighting.
    """

    chroms_all, x_all, y_all, tx_ids_all, strands_all = _load_coords_bed(coords)

    idxs_chr = [i for i, c in enumerate(chroms_all) if c == chrom]
    if not idxs_chr:
        raise ValueError(f"No entries found for chromosome '{chrom}'.")
    x_chr = x_all[idxs_chr]
    y_chr = y_all[idxs_chr]
    tx_ids_chr = [tx_ids_all[i] for i in idxs_chr]
    strs_chr = [strands_all[i] for i in idxs_chr]

    region_start = region_end = None
    if region is not None:
        region_start, region_end = region
        mask = (
            (x_chr >= region_start)
            & (x_chr <= region_end)
            & (y_chr >= region_start)
            & (y_chr <= region_end)
        )
        x_chr = x_chr[mask]
        y_chr = y_chr[mask]
        tx_ids_chr = [tx_ids_chr[i] for i in range(len(mask)) if mask[i]]
        strs_chr = [strs_chr[i] for i in range(len(mask)) if mask[i]]
        if len(x_chr) == 0:
            raise ValueError(f"No points remain in region {region_start}-{region_end}.")

    x = x_chr
    y = y_chr
    tx_ids = tx_ids_chr
    strs = strs_chr
    N = len(x)

    if region_start is not None:
        xmin, xmax = region_start, region_end
        ymin, ymax = region_start, region_end
    else:
        xmin, xmax, ymin, ymax = _compute_limits(x, y)

    mapping_list = _load_mapping(mapping, chrom, region_start, region_end)

    cluster_ids = np.full(N, -1, dtype=int)
    if mapping_list:
        for i in range(N):
            for entry in mapping_list:
                if (
                    y[i] >= entry["tss_min"]
                    and y[i] <= entry["tss_max"]
                    and x[i] >= entry["tts_min"]
                    and x[i] <= entry["tts_max"]
                    and strs[i].strip() == entry["strand"]
                    and tx_ids[i] in entry["tx_set"]
                ):
                    cluster_ids[i] = entry["cluster_idx"]
                    break

    overlap = np.zeros(N, dtype=int)
    if N > 1:
        fig_w, fig_h = 10.0, 5.0
        panel_w = (fig_w / 2) * 0.5
        panel_h = fig_h * 0.5
        dx_scale = (xmax - xmin) or 1
        dy_scale = (ymax - ymin) or 1
        r_pt = math.sqrt(9 / math.pi)
        r_in = r_pt / 75
        d = 2 * r_in * 0.74
        thresh2 = d * d
        for i in range(N - 1):
            dx = (x[i] - x[i + 1 :]) / dx_scale * panel_w
            dy = (y[i] - y[i + 1 :]) / dy_scale * panel_h
            m = dx * dx + dy * dy <= thresh2
            overlap[i] += m.sum()
            overlap[i + 1 :][m] += 1
    pct95 = np.percentile(overlap, 95) if N > 0 else 0
    color_cap = max(5, int(math.ceil(pct95 / 5.0) * 5))
    norm = Normalize(vmin=0, vmax=color_cap)

    max_coord = max(abs(xmin), abs(xmax), abs(ymin), abs(ymax))
    scale = 1_000_000 if max_coord >= 1_000_000 else 1_000
    unit_label = "Mb" if scale == 1_000_000 else "Kb"
    coord_formatter = _coord_formatter(scale)

    out_dir = output.parent
    out_base = output.name
    prefix = (
        f"{chrom}_{region_start}-{region_end}_" if region_start is not None else f"{chrom}_"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    strand_groups = defaultdict(list)
    for i, s in enumerate(strs):
        strand_groups[s.strip()].append(i)

    for strand_symbol, idxs in strand_groups.items():
        mapping_strand = [e for e in mapping_list if e["strand"] == strand_symbol]

        x_s = x[idxs]
        y_s = y[idxs]
        overlap_s = overlap[idxs]
        clust_s = cluster_ids[idxs]
        N_s = len(x_s)

        hx_s, xed_s = _compute_hist(x_s, nbins=N_s, data_min=xmin, data_max=xmax)
        hy_s, yed_s = _compute_hist(y_s, nbins=N_s, data_min=ymin, data_max=ymax)
        hx_s = np.clip(hx_s, None, 240)
        hy_s = np.clip(hy_s, None, 240)
        # Dynamic count cap: use actual histogram max rounded to 5, but no more than 250
        data_max_count = int(max(hx_s.max() if hx_s.size > 0 else 0, hy_s.max() if hy_s.size > 0 else 0))
        cap = max(5, int(math.ceil(data_max_count / 5.0) * 5))
        cap = min(250, cap)

        fig = plt.figure(figsize=(10, 5), dpi=600)
        orig_mx, orig_left_w, orig_bot_h, orig_bot_w, orig_top_h = (
            0.1273,
            0.083,
            0.50,
            0.50,
            0.083,
        )
        orig_my = (1 - (orig_bot_h + 0.025 + orig_top_h)) / 2
        orig_sc_x = orig_mx + orig_left_w + 0.025
        orig_top_y = orig_my + orig_bot_h + 0.025
        orig_cax_x = orig_mx + orig_left_w + orig_bot_w + 0.0565

        def left_coords(o):
            return [o[0] * 0.5, o[1], o[2] * 0.5, o[3]]

        def right_coords(o):
            return [0.5 + o[0] * 0.5, o[1], o[2] * 0.5, o[3]]

        ax_left_l = fig.add_axes(left_coords([orig_mx, orig_my, orig_left_w, orig_bot_h]))
        ax_sc_l = fig.add_axes(left_coords([orig_sc_x, orig_my, orig_bot_w, orig_bot_h]))
        ax_top_l = fig.add_axes(left_coords([orig_sc_x, orig_top_y, orig_bot_w, orig_top_h]))
        cax_l = fig.add_axes(left_coords([orig_cax_x, orig_my, 0.0333, orig_bot_h]))
        ax_left_r = fig.add_axes(right_coords([orig_mx, orig_my, orig_left_w, orig_bot_h]))
        ax_sc_r = fig.add_axes(right_coords([orig_sc_x, orig_my, orig_bot_w, orig_bot_h]))
        ax_top_r = fig.add_axes(right_coords([orig_sc_x, orig_top_y, orig_bot_w, orig_top_h]))

        if N_s > 0:
            order = np.argsort(overlap_s)
            ax_sc_l.scatter(
                x_s[order],
                y_s[order],
                s=9,
                c=overlap_s[order],
                cmap="plasma",
                norm=norm,
                marker="o",
                linewidths=0,
                alpha=0.1,
                rasterized=True,
            )
        ax_sc_l.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
        ax_sc_l.xaxis.set_major_formatter(coord_formatter)
        ax_sc_l.yaxis.set_major_formatter(coord_formatter)
        ax_sc_l.set_xticks(np.linspace(xmin, xmax, 5))
        ax_sc_l.set_yticks(np.linspace(ymin, ymax, 5))
        ax_sc_l.tick_params(axis="x", labelsize=9)
        ax_sc_l.tick_params(axis="y", labelleft=False)
        ax_sc_l.set_xlabel(f"TTS ({unit_label})", fontsize=9)
        sm = plt.cm.ScalarMappable(norm=norm, cmap="plasma")
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cax_l, orientation="vertical", ticks=[0, color_cap])
        cbar.ax.set_yticklabels(["0", str(color_cap)], fontsize=9)
        cbar.ax.tick_params(axis="y", length=4)
        bw = xed_s[1] - xed_s[0] if len(xed_s) > 1 else 1
        ax_top_l.bar(
            xed_s[:-1],
            hx_s,
            width=bw,
            align="edge",
            edgecolor="black",
            linewidth=0.3,
            color=(120 / 255, 172 / 255, 145 / 255),
        )
        ax_top_l.set(xlim=(xmin, xmax))
        ax_top_l.xaxis.set_major_formatter(coord_formatter)
        ax_top_l.set_xticks([])
        ax_top_l.set_yticks([0, cap])
        ax_top_l.set_yticklabels(["0", str(cap)], fontsize=9)
        ax_top_l.tick_params(axis="y", labelsize=9)
        ax_top_l.set_ylabel("Count", fontsize=9)
        bh = yed_s[1] - yed_s[0] if len(yed_s) > 1 else 1
        ax_left_l.barh(
            yed_s[:-1],
            -hy_s,
            height=bh,
            align="edge",
            color="grey",
            edgecolor="black",
            linewidth=0.3,
        )
        ax_left_l.set(ylim=(ymin, ymax))
        ax_left_l.yaxis.set_major_formatter(coord_formatter)
        ax_left_l.set_yticks(np.linspace(ymin, ymax, 5))
        ax_left_l.set_xticks([-cap, 0])
        ax_left_l.set_xticklabels([str(cap), "0"], fontsize=9)
        ax_left_l.tick_params(axis="both", labelsize=9)
        ax_left_l.set_ylabel(f"TSS ({unit_label})", fontsize=9)

        if mapping_strand:
            assigned = clust_s >= 0
            uniq = np.unique(clust_s[assigned]) if np.any(assigned) else np.array([])
            n_cl = len(uniq)
            cmap = (
                ListedColormap(plt.get_cmap("tab20").colors[:n_cl])
                if n_cl <= 20
                else ListedColormap(plt.get_cmap("hsv", n_cl).colors)
            )
            ax_top_r.set(xlim=(xmin, xmax), ylim=(0, cap))
            ax_top_r.xaxis.set_major_formatter(coord_formatter)
            ax_top_r.set_xticks([])
            ax_top_r.set_yticks([])
            ax_top_r.set_ylabel("Windows", fontsize=9)
            ax_top_r.yaxis.set_label_coords(-0.11, 0.5)
            for e in mapping_strand:
                cl = e["cluster_idx"]
                color = cmap(list(uniq).index(cl)) if cl in uniq else "black"
                ax_top_r.add_patch(
                    Rectangle(
                        (e["tts_min"], 0),
                        e["tts_max"] - e["tts_min"],
                        cap,
                        facecolor=color,
                        edgecolor=color,
                        alpha=0.5,
                        linewidth=0,
                    )
                )
            ax_left_r.set(ylim=(ymin, ymax), xlim=(-cap, 0))
            ax_left_r.yaxis.set_major_formatter(coord_formatter)
            ax_left_r.set_yticks(np.linspace(ymin, ymax, 5))
            ax_left_r.set_xticks([])
            ax_left_r.set_ylabel(f"TSS ({unit_label})", fontsize=9)
            for e in mapping_strand:
                cl = e["cluster_idx"]
                color = cmap(list(uniq).index(cl)) if cl in uniq else "black"
                ax_left_r.add_patch(
                    Rectangle(
                        (-cap, e["tss_min"]),
                        cap,
                        e["tss_max"] - e["tss_min"],
                        facecolor=color,
                        edgecolor=color,
                        alpha=0.5,
                        linewidth=0,
                    )
                )
        else:
            ax_top_r.set(xlim=(xmin, xmax), ylim=(0, cap))
            ax_top_r.xaxis.set_major_formatter(coord_formatter)
            ax_top_r.set_xticks([])
            ax_top_r.set_yticks([])
            ax_left_r.set(ylim=(ymin, ymax), xlim=(-cap, 0))
            ax_left_r.yaxis.set_major_formatter(coord_formatter)
            ax_left_r.set_yticks(np.linspace(ymin, ymax, 5))
            ax_left_r.set_xticks([])

        if N_s > 0:
            if mapping_strand:
                assigned = clust_s >= 0
                uniq = np.unique(clust_s[assigned]) if np.any(assigned) else np.array([])
                n_cl = len(uniq)
                cmap2 = (
                    ListedColormap(plt.get_cmap("tab20").colors[:n_cl])
                    if n_cl <= 20
                    else ListedColormap(plt.get_cmap("hsv", n_cl).colors)
                )
            for i in range(N_s):
                xi, yi = x_s[i], y_s[i]
                if mapping_strand and clust_s[i] >= 0:
                    color = cmap2(list(uniq).index(clust_s[i]))
                else:
                    color = "black"
                ax_sc_r.scatter(
                    [xi],
                    [yi],
                    s=20,
                    facecolors=color,
                    edgecolors="none",
                    linewidths=0.5,
                    alpha=0.6,
                    marker="o",
                    rasterized=True,
                )
        ax_sc_r.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
        ax_sc_r.xaxis.set_major_formatter(coord_formatter)
        ax_sc_r.yaxis.set_major_formatter(coord_formatter)
        ax_sc_r.set_xticks(np.linspace(xmin, xmax, 5))
        ax_sc_r.set_yticks(np.linspace(ymin, ymax, 5))
        ax_sc_r.tick_params(axis="x", labelsize=9)
        ax_sc_r.tick_params(axis="y", labelleft=False)
        ax_sc_r.set_xlabel(f"TTS ({unit_label})", fontsize=9)

        new_base = f"{prefix}{strand_symbol}_{out_base}"
        plt.savefig(out_dir / new_base, dpi=600)
        plt.close(fig)
