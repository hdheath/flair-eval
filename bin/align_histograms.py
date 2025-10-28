from __future__ import annotations

from pathlib import Path
from typing import Sequence, Dict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _hist(vals: Sequence[float] | Sequence[int], out: Path, xlabel: str) -> str | None:
    if not vals:
        return None
    plt.figure()
    plt.hist(vals, bins=50)
    plt.xlabel(xlabel)
    plt.ylabel("count")
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()
    return out.name


def generate_histograms(
    mapq_vals: Sequence[int],
    identity_vals: Sequence[float],
    read_len_vals: Sequence[int],
    out_dir: Path,
) -> Dict[str, str | None]:
    """Generate align QC histograms and return mapping of type->filename."""
    out_dir.mkdir(parents=True, exist_ok=True)
    mapping: Dict[str, str | None] = {}
    mapping["mapq"] = _hist(mapq_vals, out_dir / "align_mapq_hist.png", "MAPQ")
    mapping["identity"] = _hist(identity_vals, out_dir / "align_identity_hist.png", "Read identity (%)")
    mapping["length"] = _hist(read_len_vals, out_dir / "align_length_hist.png", "Read length (bp)")
    return mapping
