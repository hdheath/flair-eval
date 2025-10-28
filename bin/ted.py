from __future__ import annotations

import io
import json
import logging
import os
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from . import register, write_metrics  # QC plumbing
from .qc_utils import count_lines
try:  # pragma: no cover - optional dependency
    from ..plotting import tss_tts_scatter
except Exception:  # pragma: no cover - handle missing heavy deps
    tss_tts_scatter = None
try:  # pragma: no cover - optional dependency
    from ..plotting import transcriptome_browser
except Exception:  # pragma: no cover - handle missing heavy deps
    transcriptome_browser = None
from ..lib import PathBuilder


# Use module logger
logger = logging.getLogger(__name__)

# ────────────────────────── tiny utils ──────────────────────────
def _which(exe: str) -> bool:
    from shutil import which
    return which(exe) is not None


def _resolve(p: Optional[str], data_dir: Path) -> Optional[Path]:
    """Resolve ``p`` against ``data_dir`` with extra debug logging.

    This helper expands user/environment variables and, when ``p`` is a
    relative path, joins it to ``data_dir``.  Logging the before/after values
    makes it easier to diagnose misconfigured paths that lead to missing peak
    files and unset precision/recall metrics.
    """


    if not p:
        logger.debug("[TED] _resolve called with empty path; returning None")
        return None

    expanded = os.path.expanduser(os.path.expandvars(p))
    pp = Path(expanded)
    resolved = (data_dir / pp).resolve() if not pp.is_absolute() else pp
    logger.debug(
        f"[TED] Resolving path '{p}' against data_dir '{data_dir}' -> '{resolved}'"
    )
    return resolved


def _cfg_get(cfg_obj, path: List[str], default=None):
    """Safely walk a mixed object/``dict`` config structure."""

    cur = cfg_obj
    for k in path:
        if isinstance(cur, dict):
            cur = cur.get(k)
        else:
            cur = getattr(cur, k, None)
        if cur is None:
            return default
    return cur


def _run(cmd: List[str]) -> subprocess.CompletedProcess:
    logger.debug(f"[TED] RUN: {' '.join(cmd)}")
    try:
        return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"[TED] Command failed ({' '.join(cmd)}): {e.stderr.strip()}")
        raise


def _read_bed6(path: Path) -> pd.DataFrame:
    """
    Read peaks as BED6 if possible; tolerate BED3 by adding Strand='.' (strand-insensitive).
    """
    try:
        return pd.read_csv(
            path, sep="\t", header=None, usecols=[0, 1, 2, 5],
            names=["Chrom", "Start", "End", "Strand"],
            dtype={"Chrom": str, "Start": np.int32, "End": np.int32, "Strand": str},
            comment="#",
        )
    except Exception:
        # Fallback: BED3
        logger.warning(f"[TED] {path} does not look like BED6; treating as BED3 (strand='.') which may zero-out recall.")
        df = pd.read_csv(
            path, sep="\t", header=None, usecols=[0, 1, 2],
            names=["Chrom", "Start", "End"],
            dtype={"Chrom": str, "Start": np.int32, "End": np.int32},
            comment="#",
        )
        df["Strand"] = "."
        return df


def _prepare_bed6_sorted(path: Path, tmps: List[Path]) -> Path:
    trimmed = path.with_suffix(path.suffix + ".trimmed.tmp")
    sorted_p = path.with_suffix(path.suffix + ".sorted.tmp")

    # Trim to first six columns of BED and write to temporary file
    res = _run(["cut", "-f", "1-6", str(path)])
    trimmed.write_text(res.stdout)

    # Sort the trimmed BED and capture output
    res = _run(["bedtools", "sort", "-i", str(trimmed)])
    sorted_p.write_text(res.stdout)

    try:
        trimmed.unlink()
    except Exception:
        pass
    tmps.append(sorted_p)
    logger.debug(f"[TED] Prepared sorted BED6: {sorted_p}")
    return sorted_p


def _run_closest(a: Path, b: Path) -> pd.DataFrame:
    res = _run(["bedtools", "closest", "-a", str(a), "-b", str(b), "-s", "-d"])
    return pd.read_csv(io.StringIO(res.stdout), sep="\t", header=None, comment="#")


def _extract_distance_and_peak(df: pd.DataFrame, label: str, max_dist: int) -> pd.DataFrame:
    orig_cols = df.shape[1]
    tx_series = df.iloc[:, 3].astype(str)
    dist_series = pd.to_numeric(df.iloc[:, orig_cols - 1], errors="coerce").fillna(max_dist + 1).astype(int)
    # Peak tuple (Chrom, Start, End) for -b
    peak_series = (df.iloc[:, 6].astype(str) + "_" + df.iloc[:, 7].astype(str) + "_" + df.iloc[:, 8].astype(str))
    out = pd.DataFrame({"isoform_id": tx_series, f"{label}_dist": dist_series, f"{label}_peak": peak_series})
    return out.groupby("isoform_id", as_index=False).first()


def _vectorized_overlap_counts(bed_df: pd.DataFrame, peaks_df: pd.DataFrame, window: int) -> Tuple[int, int]:
    """Return (#peaks matched by ≥1 isoform within window, total peaks)."""
    consumed = np.zeros(len(peaks_df), dtype=bool)
    for (c, s), chunk in bed_df.groupby(["Chrom", "Strand"], sort=False):
        mask = (peaks_df["Chrom"] == c) & (
            (peaks_df["Strand"] == s) | (peaks_df["Strand"] == ".")
        )
        idxs = np.flatnonzero(mask)
        if idxs.size == 0:
            continue
        sub = peaks_df.iloc[idxs]
        starts = sub["Start"].to_numpy()
        ends = sub["End"].to_numpy()
        order = np.argsort(starts)
        starts = starts[order]
        ends = ends[order]
        global_i = idxs[order]
        for _, row in chunk.iterrows():
            lo, hi = int(row["Start"]) - window, int(row["End"]) + window
            r = np.searchsorted(starts, hi, side="right")
            if r == 0:
                continue
            hits = np.nonzero(ends[:r] >= lo)[0]
            if hits.size == 0:
                continue
            matched = global_i[hits]
            new_hits = matched[~consumed[matched]]
            consumed[new_hits] = True
    return int(consumed.sum()), int(len(peaks_df))


def _safe_f1(p: Optional[float], r: Optional[float]) -> Optional[float]:
    if p is None or r is None:
        return None
    s = p + r
    return (2 * p * r / s) if s > 0 else 0.0


def _extract_gene_id(name: str) -> Optional[str]:
    """
    Isoform name like: <read>_ENSG00000310526.1 → return ENSG00000310526 (drop version).
    If no recognizable gene suffix, return None.
    """
    if not name:
        return None
    tail = name.split("_")[-1]
    if tail.startswith(("ENS", "ENSG", "ENST", "ENSMUSG", "ENSMUST")):
        return tail.split(".")[0]
    m = re.search(r"(ENS[GMT][A-Z]*\d+)(?:\.\d+)?", name)
    return m.group(1) if m else None


def _isoform_counts(iso_bed: Path) -> Tuple[int, int]:
    df = pd.read_csv(iso_bed, sep="\t", header=None, comment="#", usecols=[3], names=["Name"], dtype=str)
    n_iso = len(df)
    genes = {g for g in (_extract_gene_id(nm) for nm in df["Name"].tolist()) if g}
    return n_iso, len(genes)


def _read_map_unique_reads(map_path: Path) -> int:
    """Count unique read IDs across all isoforms in *.isoform.read.map.txt."""
    if not map_path.exists() or map_path.stat().st_size == 0:
        return 0
    uniq = set()
    with open(map_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or "\t" not in line:
                continue
            _, rhs = line.split("\t", 1)
            for rid in rhs.split(","):
                rid = rid.strip()
                if rid:
                    uniq.add(rid)
    return len(uniq)


def _synthesize_mapping_from_map_and_bam(
    map_path: Path,
    bam_path: Optional[Path],
    iso_bed: Path,
    chrom: str,
    region_start: Optional[int],
    region_end: Optional[int],
    out_dir: Path,
) -> Optional[Path]:
    """Create a mapping TSV compatible with tss_tts_scatter from an
    isoform.read.map.txt and a region BAM.

    Returns path to generated TSV or None on failure.
    """
    if not map_path or not map_path.exists():
        logger.debug(f"[TED] No isoform.read.map.txt to synthesize mapping: {map_path}")
        return None
    if not bam_path or not bam_path.exists():
        logger.warning(f"[TED] Cannot synthesize mapping: BAM not found: {bam_path}")
        return None

    # Parse isoform->reads map
    iso_to_reads: Dict[str, set] = {}
    with open(map_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or "\t" not in line:
                continue
            left, rhs = line.split("\t", 1)
            iso = left.strip()
            reads = {r.strip() for r in rhs.split(",") if r.strip()}
            if reads:
                iso_to_reads.setdefault(iso, set()).update(reads)

    if not iso_to_reads:
        logger.warning(f"[TED] isoform.read.map.txt contained no assignments: {map_path}")
        return None

    # Load isoform strands from BED (fallback to read strand if missing)
    iso_strand: Dict[str, str] = {}
    try:
        df_iso = _read_bed6(iso_bed)
        if not df_iso.empty:
            for _, r in df_iso.iterrows():
                # Name column may be present at col 3 in BED; try to extract from file name
                # We only used _read_bed6 which returns Chrom,Start,End,Strand; need Name separately
                pass
    except Exception:
        # ignore — we'll infer strand from reads
        pass

    # Fetch reads from BAM in region and collect coords
    try:
        import pysam
    except Exception:
        logger.warning("[TED] pysam not available; cannot synthesize mapping from BAM")
        return None

    read_names_needed = set(r for s in iso_to_reads.values() for r in s)
    read_coords: Dict[str, Tuple[str, int, int, str]] = {}
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bf:
            # iterate over region if provided, else whole file (but limit to reads we need)
            iterator = bf.fetch(chrom, region_start, region_end) if (region_start is not None and region_end is not None) else bf.fetch()
            for rd in iterator:
                if rd.is_unmapped or rd.is_secondary or rd.is_supplementary:
                    continue
                qn = rd.query_name
                if qn not in read_names_needed:
                    continue
                start = rd.reference_start
                end = rd.reference_end - 1 if rd.reference_end is not None else rd.reference_end
                strand = '-' if rd.is_reverse else '+'
                tss = start if strand == '+' else end
                tts = end if strand == '+' else start
                read_coords[qn] = (bf.get_reference_name(rd.reference_id), int(tss), int(tts), strand)
                # early exit if we've seen all
                if len(read_coords) >= len(read_names_needed):
                    break
    except Exception as e:
        logger.warning(f"[TED] Error reading BAM while synthesizing mapping: {e}")
        return None

    # Build mapping rows per-isoform
    rows = []
    from collections import Counter
    for iso, reads in iso_to_reads.items():
        coords = [read_coords[r] for r in reads if r in read_coords]
        if not coords:
            continue
        # filter to same chrom
        coords = [c for c in coords if c[0] == chrom]
        if not coords:
            continue
        tss_list = [c[1] for c in coords]
        tts_list = [c[2] for c in coords]
        strands = [c[3] for c in coords]
        # most frequent tss coordinate
        most_tss = Counter(tss_list).most_common(1)[0][0]
        strand_mode = Counter(strands).most_common(1)[0][0]
        tss_min, tss_max = min(tss_list), max(tss_list)
        tts_min, tts_max = min(tts_list), max(tts_list)
        # small padding
        pad = max(1, int((tss_max - tss_min) * 0.05))
        tss_min_p, tss_max_p = tss_min - pad, tss_max + pad
        pad2 = max(1, int((tts_max - tts_min) * 0.05))
        tts_min_p, tts_max_p = tts_min - pad2, tts_max + pad2
        most_freq_coord = f"{chrom}:{most_tss}"
        all_tx_ids = iso
        joint_window = f"{tss_min_p}-{tss_max_p}|{tts_min_p}-{tts_max_p}"
        rows.append({
            "most_freq_coord": most_freq_coord,
            "all_tx_ids": all_tx_ids,
            "joint_window": joint_window,
            "strand": strand_mode,
        })

    if not rows:
        logger.warning("[TED] No mapping rows synthesized from reads")
        return None

    out_dir.mkdir(parents=True, exist_ok=True)
    out_p = out_dir / (map_path.stem + ".synth.mapping.tsv")
    pd.DataFrame(rows).to_csv(out_p, sep="\t", index=False)
    logger.info(f"[TED] Wrote synthesized mapping TSV: {out_p}")
    return out_p


def _build_peaks_cfg(cfg, stage_name: str) -> Tuple[Dict[str, str], int]:
    """Collect peak file configuration for a stage.

    TSS/TTS peak BED files are now sourced exclusively from the shared
    run-level inputs.  ``experiment_*`` and ``reference_*`` keys are required
    under ``[run]`` and no longer honored from ``[qc.<stage>.TED]`` blocks or
    per-stage flags.  ``window`` may still be overridden under the QC block.
    """

    qc_block = _cfg_get(cfg, ["qc", stage_name, "TED"], {}) or {}
    window = int(qc_block.get("window", 50))

    run_cfg = getattr(cfg, "run", None)
    peaks_cfg = {
        "prime5": getattr(run_cfg, "experiment_5_prime_regions_bed_file", None) if run_cfg else None,
        "prime3": getattr(run_cfg, "experiment_3_prime_regions_bed_file", None) if run_cfg else None,
        "ref_prime5": getattr(run_cfg, "reference_5_prime_regions_bed_file", None) if run_cfg else None,
        "ref_prime3": getattr(run_cfg, "reference_3_prime_regions_bed_file", None) if run_cfg else None,
    }

    missing = [k for k, v in peaks_cfg.items() if not v]
    if missing:
        raise RuntimeError(
            "TSS/TTS peak BED files must be provided under [run]; missing: "
            + ", ".join(missing)
        )

    return peaks_cfg, window


def _audit(
    rows: List[Dict],
    context: str,
    tag: str,
    label: str,
    path: Optional[Path],
    *,
    line_counter=count_lines,
    lines: Optional[int] = None,
    derived: Optional[Dict] = None,
    note: str = "",
) -> None:
    """Append a standardized audit record for a file."""
    exists = bool(path and path.exists())
    size = path.stat().st_size if exists else None
    if exists:
        if lines is not None:
            lines_est = lines
        elif line_counter:
            try:
                lines_est = line_counter(path)
            except Exception:
                lines_est = None
        else:
            lines_est = None
    else:
        lines_est = None
    rows.append(
        {
            "context": context,
            "tag": tag,
            "label": label,
            "path": str(path) if path else None,
            "exists": exists,
            "size_bytes": size,
            "lines_est": lines_est,
            "derived_count": derived,
            "note": note,
        }
    )


# ────────────────────────── discovery helpers ──────────────────────────
def _dirs_with_file(root: Path, stage: str, fname: str) -> List[Path]:
    base = root / stage
    if not base.exists():
        return []
    dirs = [p for p in base.iterdir() if p.is_dir() and (p / fname).exists()]
    dirs.sort(key=lambda p: (p / fname).stat().st_mtime, reverse=True)
    logger.debug(f"[TED] Searching {base} for {fname} -> {[str(d) for d in dirs]}")
    return dirs


def _build_region_metrics_index(run_root: Path, regionalize_dir: Path | None = None) -> Dict[str, Dict]:
    """
    Build index: region_tag -> {'gene_count', 'transcript_count', 'reg_dir'}.

    If ``regionalize_dir`` is provided, only that directory is scanned.
    Otherwise, scan all ``regionalize/*`` directories.
    """
    idx: Dict[str, Dict] = {}
    if regionalize_dir is not None:
        dirs = [regionalize_dir]
    else:
        regionalize_root = run_root / "regionalize"
        if not regionalize_root.exists():
            return idx
        dirs = [d for d in regionalize_root.iterdir() if d.is_dir()]
    for d in dirs:
        metrics = d / "region_metrics.tsv"
        if not metrics.exists():
            continue
        try:
            df = pd.read_csv(metrics, sep="\t")
        except Exception as e:
            logger.warning(f"[TED] Failed reading {metrics}: {e}")
            continue
        if "region_tag" not in df.columns:
            continue
        for _, row in df.iterrows():
            tag = str(row["region_tag"])
            rec = {
                "gene_count": int(row.get("gene_count", 0)) if pd.notna(row.get("gene_count", np.nan)) else 0,
                "transcript_count": int(row.get("transcript_count", 0)) if pd.notna(row.get("transcript_count", np.nan)) else 0,
                "reg_dir": d,
            }
            idx[tag] = rec
    logger.debug(f"[TED] Region index entries: {len(idx)}")
    return idx


def _find_correct_bed_for_tag(run_root: Path, tag: str) -> Optional[Path]:
    for d in _dirs_with_file(run_root, "correct", f"{tag}_all_corrected.bed"):
        return d / f"{tag}_all_corrected.bed"
    return None


def _find_correct_bed_single(run_root: Path, run_id: str) -> Optional[Path]:
    for d in _dirs_with_file(run_root, "correct", f"{run_id}_all_corrected.bed"):
        return d / f"{run_id}_all_corrected.bed"
    return None


def _find_align_bam(run_root: Path, run_id: str) -> Optional[Path]:
    for d in _dirs_with_file(run_root, "align", f"{run_id}_flair.bam"):
        return d / f"{run_id}_flair.bam"
    return None

def _count_primary_alignments_bam(bam: Path) -> int:
    if not _which("samtools"):
        raise RuntimeError("TED requires 'samtools' on PATH to count primary alignments from BAM.")
    if not bam.exists() or bam.stat().st_size == 0:
        return 0
    # exclude unmapped(0x4), secondary(0x100), supplementary(0x800)
    res = _run(["samtools", "view", "-c", "-F", "2308", str(bam)])
    try:
        return int(res.stdout.strip() or "0")
    except Exception:
        return 0


def _is_region_tag(stem: str) -> bool:
    """chr_start_end[...].isoforms.bed → True if start/end are integers."""
    m = re.match(r"^(.+?)_(\d+)_(\d+)$", stem)  # stem has no .isoforms suffix when we call this
    return m is not None


def _tss_tts_metrics_full(iso_bed: Path, peaks: Dict[str, Optional[Path]], window: int,
                          audit_rows: List[Dict], tag_ctx: str) -> Dict[str, Optional[float]]:
    """
    Compute precision/recall/F1 for experimental and reference TSS/TTS:
      - experimental: keys 'prime5' → 5prime_*, 'prime3' → 3prime_*
      - reference:    keys 'ref_prime5' → ref5prime_*, 'ref_prime3' → ref3prime_*
    If a file is missing, that side returns None for all three metrics.
    """
    have_any = any(peaks.values())
    metrics: Dict[str, Optional[float]] = {
        "5prime_precision": None, "5prime_recall": None, "5prime_f1": None,
        "3prime_precision": None, "3prime_recall": None, "3prime_f1": None,
        "ref5prime_precision": None, "ref5prime_recall": None, "ref5prime_f1": None,
        "ref3prime_precision": None, "ref3prime_recall": None, "ref3prime_f1": None,
    }
    if not have_any:
        logging.warning("[TED] No peak files found; precision/recall metrics will be None")
        for key_cfg, label in [
            ("prime5", "5"), ("prime3", "3"), ("ref_prime5", "ref5"), ("ref_prime3", "ref3"),
        ]:
            _audit(
                audit_rows,
                tag_ctx,
                label,
                f"{key_cfg}_peaks",
                peaks.get(key_cfg),
                note="missing",
            )
        return metrics

    if not _which("bedtools"):
        raise RuntimeError("TED requires 'bedtools' on PATH when TSS/TTS files are provided.")


    full_bed_df = pd.read_csv(
        iso_bed, sep="\t", header=None, comment="#",
        usecols=[0, 1, 2, 3, 5],
        names=["Chrom", "Start", "End", "Name", "Strand"],
        dtype={"Chrom": str, "Start": np.int32, "End": np.int32, "Name": str, "Strand": str},
    )
    n_tx = len(full_bed_df)
    if n_tx == 0:
        logger.warning(f"[TED] Isoforms BED '{iso_bed}' has 0 transcripts; "
                       f"precision undefined for {tag_ctx}.")
    bed_df = full_bed_df[["Chrom", "Start", "End", "Strand"]]

    tmp_local: List[Path] = []
    try:
        bed_sorted = _prepare_bed6_sorted(iso_bed, tmp_local)


        def _side(key_cfg: str, label: str) -> Tuple[Optional[float], Optional[float], Optional[float]]:
            pth = peaks.get(key_cfg)
            if pth is None or not pth.exists() or pth.stat().st_size == 0:
                note = "missing" if (pth is None or not pth.exists()) else "empty"
                logger.warning(
                    f"[TED] Peaks file for key '{key_cfg}' {note}: {pth}"
                )
                _audit(audit_rows, tag_ctx, label, f"{key_cfg}_peaks", pth, note=note)
                return None, None, None

            # Sort/trim peaks, too
            peaks_sorted = _prepare_bed6_sorted(pth, tmp_local)
            dfc = _run_closest(bed_sorted, peaks_sorted)
            ddf = _extract_distance_and_peak(dfc, label, window)
            m = int((ddf[f"{label}_dist"] <= window).sum())
            precision = (m / n_tx) if n_tx else None

            peaks_df = _read_bed6(pth)
            if peaks_df.empty:
                logger.warning(f"[TED] Peaks file '{pth}' parsed to 0 intervals; "
                               f"recall undefined for {tag_ctx}/{key_cfg}.")
            c, t = _vectorized_overlap_counts(bed_df, peaks_df, window)
            recall = (c / t) if t else None
            f1 = _safe_f1(precision, recall)

            logger.debug(f"[TED] {tag_ctx}/{key_cfg}: peaks rows={len(peaks_df)}, "
                         f"isoforms={n_tx}, window={window}")

            _audit(
                audit_rows,
                tag_ctx,
                label,
                f"{key_cfg}_peaks",
                pth,
                derived={"matches_le_window": m, "peaks_total": t},
                note=f"precision={precision}, recall={recall}, f1={f1}",
            )
            return precision, recall, f1

        p5, r5, f5 = _side("prime5", "5")
        p3, r3, f3 = _side("prime3", "3")
        rp5, rr5, rf5 = _side("ref_prime5", "ref5")
        rp3, rr3, rf3 = _side("ref_prime3", "ref3")

        metrics.update({
            "5prime_precision": p5, "5prime_recall": r5, "5prime_f1": f5,
            "3prime_precision": p3, "3prime_recall": r3, "3prime_f1": f3,
            "ref5prime_precision": rp5, "ref5prime_recall": rr5, "ref5prime_f1": rf5,
            "ref3prime_precision": rp3, "ref3prime_recall": rr3, "ref3prime_f1": rf3,
        })
        return metrics
    finally:
        for f in tmp_local:
            try:
                f.unlink()
            except Exception:
                pass


# ────────────────────────── main collector ──────────────────────────
@register("TED")
@register("ted")
@register("collapse")
@register("transcriptome")
def collect(
    stage_dir: Path,
    cfg,
    upstreams: Dict[str, "PathBuilder"] | None = None,
    *,
    out_dir: Optional[Path] = None,
) -> None:
    """
    Writes:
      - TED.tsv: per-region rows if regionalized else one row
      - TED.audit.tsv: everything we touched + counts
      - Transcriptome browser images under ``out_dir/transcriptome_browser``
      - TED sidecar with simple totals (no means)
    Implements:
      * region_genes_expected / region_isoforms_expected via an index over ALL regionalize dirs
      * assigned_pct denominator: 
          - collapse: corrected BED line count (per-region if regionalized, else single)
          - transcriptome: primary alignments in BAM (per-region regional BAM, else full align BAM)
    Parameters
    ----------
    stage_dir : Path
        Stage directory containing isoform BEDs/GTfs.
    cfg : object
        Parsed config object.
    upstreams : dict[str, PathBuilder] | None
        Optional upstream PathBuilders.
    out_dir : Path | None
        Directory where outputs should be written. Defaults to ``stage_dir``.
    """
    logging.debug("[TED] Starting TED QC collection.")
    stage_name = stage_dir.parent.name           # collapse | transcriptome
    run_root = stage_dir.parent.parent           # .../outputs/<run_id>
    run_id = run_root.name
    out_dir = out_dir or stage_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    logging.debug(f"[TED] Stage name: {stage_name}, Run ID: {run_id}")

    data_dir_cfg = Path(_cfg_get(cfg, ["run", "data_dir"], "."))
    if data_dir_cfg.is_absolute():
        data_dir = data_dir_cfg
    else:
        # ``stage_dir`` lives at ``<repo>/outputs/<run>/<stage>/<hash>``. The
        # project root is therefore two levels above ``run_root``.  Earlier the
        # code only ascended once, resolving relative ``data_dir`` entries
        # against ``outputs`` and causing non-regionalized runs to miss the
        # configured TSS/TTS BED files.  Walk up to the project root before
        # joining with ``data_dir``.
        repo_root = run_root.parents[1]
        data_dir = (repo_root / data_dir_cfg).resolve()
    logging.debug(f"[TED] Data directory set to: {data_dir}")



    peaks_cfg, window = _build_peaks_cfg(cfg, stage_name)
    logging.info(f"[TED] Stage name for peak config: {stage_name}")
    logging.info(f"[TED] Parsed peaks config: {peaks_cfg}")

    # Determine if this is regionalized by matching chr_start_end pattern
    iso_beds = sorted(stage_dir.glob("*.isoforms.bed"))
    regional_files = [
        p for p in iso_beds if _is_region_tag(p.stem.replace(".isoforms", ""))
    ]
    is_regionalized = len(regional_files) > 0
    logging.debug(f"[TED] Isoform beds: {[p.name for p in iso_beds]}")
    logging.debug(f"[TED] Regional isoform beds: {[p.name for p in regional_files]}")
    logging.debug(f"[TED] is_regionalized: {is_regionalized}")

    # If no peak config for this stage and we're regionalized, fallback to regionalize stage config
    if is_regionalized and not any(peaks_cfg.values()):
        logger.info(
            f"[TED] No peak config for stage '{stage_name}', falling back to 'regionalize' config for regionalized run."
        )
        peaks_cfg, window = _build_peaks_cfg(cfg, "regionalize")
    logging.debug(f"[TED] Window: {window}")
    logging.debug(f"[TED] Peaks configuration (final): {peaks_cfg}")

    regionalize_pb = (upstreams or {}).get("regionalize")
    reg_index = _build_region_metrics_index(
        run_root, regionalize_pb.stage_dir if regionalize_pb else None
    )

    rows: List[Dict] = []
    audit_rows: List[Dict] = []
    browser_paths: Dict[str, str] = {}

    if is_regionalized:
        # ── Regionalized: iterate over discovered isoform files
        # Prefer the explicit regionalize PathBuilder when provided
        if regionalize_pb is not None:
            reg_dir_base = regionalize_pb.stage_dir
        else:
            # Fallback: try to infer from outputs/<run_id>/regionalize/<hash>
            # ``stage_dir`` is <run>/<stage>/<hash>; the hash component is the
            # name of ``stage_dir`` itself, not its parent.  Previously the
            # parent directory name (e.g. ``collapse``) was used which produced
            # paths like ``regionalize/collapse`` and caused sliced peak files
            # to be missed.  Use the stage hash instead.
            hash_dir = stage_dir.name
            reg_dir_base = stage_dir.parent.parent / "regionalize" / hash_dir
        for iso_bed in regional_files:
            tag = iso_bed.stem.replace(".isoforms", "")
            logging.debug(f"[TED] Processing region tag: {tag}")
            try:
                chrom, start, end = tag.split("_", 3)
                start_i, end_i = int(start), int(end)
            except Exception:
                logging.warning(f"[TED] Could not parse region tag from {iso_bed.name}; skipping")
                continue

            if not iso_bed.exists():
                logging.warning(f"[TED] Isoforms BED for region {tag} not found: {iso_bed}")

            map_txt = stage_dir / f"{tag}.isoform.read.map.txt"
            if not map_txt.exists():
                logging.warning(f"[TED] Isoform read map for region {tag} not found: {map_txt}")

            n_iso, n_genes_obs = _isoform_counts(iso_bed)
            if n_iso == 0:
                logger.warning(f"[TED] {iso_bed.name} contains 0 isoforms; "
                               "precision will be None for this row.")
            assigned_reads = _read_map_unique_reads(map_txt)
            logging.debug(f"[TED] Region {tag} - n_iso: {n_iso}, n_genes_obs: {n_genes_obs}, assigned_reads: {assigned_reads}")

            # Use region-index information when available to locate sliced
            # peak files and expected counts.  Fall back to ``reg_dir_base``
            # (best effort) when the tag is absent from the index.
            rec = reg_index.get(tag)
            if rec:
                reg_dir_for_tag = rec.get('reg_dir')
                reg_gene_cnt = rec.get('gene_count')
                reg_tx_cnt = rec.get('transcript_count')
            else:
                reg_dir_for_tag = reg_dir_base
                reg_gene_cnt = None
                reg_tx_cnt = None
            reg_bam = (reg_dir_for_tag / f"{tag}.bam") if reg_dir_for_tag else None

            # peaks: prefer sliced if available, otherwise global run-level path
            peaks = {}
            for key, conf in peaks_cfg.items():
                base = Path(conf).name
                sliced = (reg_dir_for_tag / f"{tag}_{base}") if reg_dir_for_tag else None
                logger.debug(f"[TED] Looking for sliced peak file for key '{key}': {sliced}")
                if sliced and sliced.exists() and sliced.stat().st_size > 0:
                    logger.debug(f"[TED] Found sliced peak file for key '{key}': {sliced}")
                    peaks[key] = sliced
                else:
                    rp = _resolve(conf, data_dir)
                    logger.debug(f"[TED] Looking for global peak file for key '{key}': {rp}")
                    if not (rp and rp.exists() and rp.stat().st_size > 0):
                        logger.warning(f"[TED] Peaks file for key '{key}' and region {tag} not found: {conf}")
                        peaks[key] = None
                    else:
                        logger.debug(f"[TED] Found global peak file for key '{key}': {rp}")
                        peaks[key] = rp

            logging.debug(f"[TED] For region {tag}, resolved peaks: { {k: str(v) if v else None for k,v in peaks.items()} }")

            tss_tts = _tss_tts_metrics_full(iso_bed, peaks, window, audit_rows, f"region:{tag}")
            logging.debug(f"[TED] Region {tag} TSS/TTS metrics: {tss_tts}")

            # Denominator per rule
            if stage_name == "collapse":
                corr_bed = _find_correct_bed_for_tag(run_root, tag)
                if corr_bed is None:
                    logging.warning(f"[TED] Corrected BED for {tag} not found; assigned_pct will be None")
                    denom = None
                else:
                    denom = count_lines(corr_bed)
                    logging.debug(f"[TED] For region {tag}, corrected BED: {corr_bed.name}, lines: {denom}")
            else:
                if reg_bam and reg_bam.exists():
                    denom = _count_primary_alignments_bam(reg_bam)
                else:
                    logging.warning(f"[TED] Regional BAM for {tag} not found; assigned_pct will be None")
                    denom = None

            row = {
                "run_id": run_id,
                "stage": stage_name,
                "region_tag": tag,
                "chrom": chrom, "start": start_i, "end": end_i,
                "span_bp": int(end_i - start_i + 1),
                "region_genes_expected": reg_gene_cnt,
                "region_isoforms_expected": reg_tx_cnt,
                "isoforms_observed": n_iso,
                "genes_observed": n_genes_obs,
                "assigned_primary_reads": int(assigned_reads),
                "align_primary_total": int(denom) if (denom is not None) else None,
                "assigned_pct": (assigned_reads / denom) if (denom and denom > 0) else None,
                **tss_tts,
            }
            logging.debug(f"[TED] Region {tag} final row: {row}")
            if reg_dir_for_tag is None:
                logging.debug(f"[TED] No regionalize index entry for {tag}; expected metrics defaulted to 0")
            rows.append(row)

            # ── Optional scatter plot for collapsed isoforms ──
            try:
                span = int(end_i - start_i + 1)
                if span > 100000:
                    logging.info(
                        f"[TED] Skipping TSS/TTS scatter plot for {tag}: region span {span} > 100000bp"
                    )
                elif tss_tts_scatter:
                    outdir = out_dir
                    out_png = outdir / "scatter.png"
                    # Validate existing mapping TSV; synthesize if missing or malformed
                    map_arg = map_txt if map_txt.exists() else None
                    if map_arg is not None:
                        try:
                            tmpdf = pd.read_csv(map_arg, sep="\t", dtype=str, nrows=1)
                            required = {"most_freq_coord", "all_tx_ids", "joint_window", "strand"}
                            if not required.issubset(set(tmpdf.columns)):
                                logger.debug(f"[TED] Mapping file {map_arg} missing required columns; attempting to synthesize from BAM")
                                synth = _synthesize_mapping_from_map_and_bam(map_txt, reg_bam, iso_bed, chrom, start_i, end_i, outdir)
                                if synth:
                                    map_arg = synth
                                else:
                                    map_arg = None
                        except Exception:
                            logger.debug(f"[TED] Could not read mapping file {map_arg}; attempting to synthesize from BAM")
                            synth = _synthesize_mapping_from_map_and_bam(map_txt, reg_bam, iso_bed, chrom, start_i, end_i, outdir)
                            if synth:
                                map_arg = synth
                            else:
                                map_arg = None
                    else:
                        # Try to synthesize from isoform.read.map.txt and BAM
                        synth = _synthesize_mapping_from_map_and_bam(map_txt, reg_bam, iso_bed, chrom, start_i, end_i, outdir)
                        if synth:
                            map_arg = synth
                        else:
                            map_arg = None

                    tss_tts_scatter.generate(
                        coords=iso_bed,
                        chrom=chrom,
                        region=(int(start_i), int(end_i)),
                        mapping=map_arg,
                        output=out_png,
                    )
                else:
                    logging.warning(
                        f"[TED] Skipping scatter plot for {tag}: dependency missing"
                    )
                # Also attempt transcriptome browser plot for small regions
                try:
                    if span <= 20000:
                        if transcriptome_browser:
                            browser_dir = out_dir / "transcriptome_browser"
                            # Resolve GTF from config if present
                            gtf_cfg = _cfg_get(cfg, ["run", "gtf"], None)
                            gtf_path = _resolve(gtf_cfg, data_dir) if gtf_cfg else None
                            # Use regional BAM if available; transcriptome_browser will fall back to provided BAM
                            bam_for_browser = reg_bam if (reg_bam and reg_bam.exists()) else None
                            cfg_obj = transcriptome_browser.Config(
                                gtf=gtf_path if gtf_path else Path("."),
                                bam=bam_for_browser,
                                genome="",
                                outdir=browser_dir,
                                mapping=map_txt if map_txt.exists() else None,
                                collapsed_isoforms=iso_bed if iso_bed.exists() else None,
                            )
                            region_str = f"{chrom}:{int(start_i)}-{int(end_i)}"
                            png_path = transcriptome_browser.generate(cfg_obj, region=region_str)
                            if png_path:
                                browser_paths[tag] = str(png_path)
                        else:
                            logging.warning(f"[TED] Skipping transcriptome browser plot for {tag}: dependency missing")
                except Exception as e:
                    logging.warning(f"[TED] transcriptome browser plot failed for region {tag}: {e}")
            except Exception as e:
                logging.warning(
                    f"[TED] scatter plot failed for region {tag}: {e}"
                )
    else:
        # ── Single (non-regionalized): one row ──
        iso_bed = stage_dir / f"{run_id}.isoforms.bed"
        if not iso_bed.exists() or iso_bed.stat().st_size == 0:
            raise RuntimeError(f"[TED] Isoforms BED not found or empty: {iso_bed}")
        map_txt = stage_dir / f"{run_id}.isoform.read.map.txt"
        if not map_txt.exists():
            logging.warning(f"[TED] Isoform read map not found: {map_txt}")

        n_iso, n_genes_obs = _isoform_counts(iso_bed)
        if n_iso == 0:
            logger.warning(f"[TED] {iso_bed.name} contains 0 isoforms; "
                           "precision will be None for this row.")
        assigned_reads = _read_map_unique_reads(map_txt)
        logging.debug(f"[TED] Single mode - n_iso: {n_iso}, n_genes_obs: {n_genes_obs}, assigned_reads: {assigned_reads}")

        # peaks: global only
        peaks = {}
        missing_peaks = []
        logging.info(f"[TED] Looking for peak files in directory: {data_dir}")
        for key, conf in peaks_cfg.items():
            if not conf:
                logging.error(f"[TED] Peaks config missing for key '{key}' in non-regionalized run; cannot proceed.")
                missing_peaks.append(key)
                peaks[key] = None
                continue
            rp = _resolve(conf, data_dir)
            logging.info(f"[TED] Resolving peak file for key '{key}': config='{conf}', resolved='{rp}'")
            if not (rp and rp.exists() and rp.stat().st_size > 0):
                logging.error(f"[TED] Peaks file for key '{key}' not found at {rp}; required for non-regionalized run.")
                missing_peaks.append(key)
                peaks[key] = None
            else:
                peaks[key] = rp
        logging.debug("[TED] Single mode peaks resolved: %s", {k: str(v) if v else None for k, v in peaks.items()})
        if missing_peaks:
            raise RuntimeError(f"[TED] Required peak files missing for keys: {', '.join(missing_peaks)}. Aborting non-regionalized run.")
        tss_tts = _tss_tts_metrics_full(iso_bed, peaks, window, audit_rows, "single")

        if stage_name == "collapse":
            corr_bed = _find_correct_bed_single(run_root, run_id)
            if corr_bed:
                denom = count_lines(corr_bed)
                logging.debug(f"[TED] Found single corrected BED: {corr_bed.name} with {denom} lines")
            else:
                logging.warning("[TED] Single corrected BED not found; assigned_pct will be None")
                denom = None
        else:
            aln_bam = _find_align_bam(run_root, run_id)
            if aln_bam:
                denom = _count_primary_alignments_bam(aln_bam)
                logging.debug(f"[TED] Found single align BAM: {aln_bam.name}")
            else:
                logging.warning("[TED] Align BAM not found; assigned_pct will be None")
                denom = None

        row = {
            "run_id": run_id,
            "stage": stage_name,
            "region_tag": None,
            "chrom": None, "start": None, "end": None, "span_bp": None,
            "region_genes_expected": None,
            "region_isoforms_expected": None,
            "isoforms_observed": n_iso,
            "genes_observed": n_genes_obs,
            "assigned_primary_reads": int(assigned_reads),
            "align_primary_total": int(denom) if (denom is not None) else None,
            "assigned_pct": (assigned_reads / denom) if (denom and denom > 0) else None,
            **tss_tts,
        }
        logging.debug(f"[TED] Single mode final row: {row}")
        rows.append(row)

    # write TSV
    out_tsv = out_dir / "TED.tsv"
    df = pd.DataFrame(rows)
    if is_regionalized:
        df.to_csv(out_tsv, sep="\t", index=False)
    else:
        drop_cols = [
            "region_tag",
            "chrom",
            "start",
            "end",
            "span_bp",
            "region_genes_expected",
            "region_isoforms_expected",
        ]
        df.drop(columns=drop_cols, inplace=True)
        df.to_csv(out_tsv, sep="\t", index=False)
    logging.info(f"[TED] Wrote TED.tsv with {len(rows)} rows at {out_tsv}")

    # Write audit TSV
    audit_tsv = out_dir / "TED.audit.tsv"
    pd.DataFrame(
        audit_rows,
        columns=["context","tag","label","path","exists","size_bytes","lines_est","derived_count","note"]
    ).to_csv(audit_tsv, sep="\t", index=False)
    logger.info(f"[TED] Wrote audit log at {audit_tsv}")

    # sidecar summary (simple totals—no means)
    assigned_total = int(sum((r.get("assigned_primary_reads") or 0) for r in rows))
    metrics = {
        "window": int(window),
        "rows": len(rows),
        "assigned_primary_reads_total": assigned_total,
        "tsv": str(out_tsv),
    }
    write_metrics(out_dir, "TED", metrics)
    logging.info(f"[TED] Sidecar metrics written: {metrics}")

    if browser_paths:
        br_dir = out_dir / "transcriptome_browser"
        br_dir.mkdir(parents=True, exist_ok=True)
        (br_dir / "region_map.json").write_text(json.dumps(browser_paths, indent=2))
