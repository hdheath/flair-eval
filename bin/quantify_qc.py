from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from . import register, write_metrics
from .qc_utils import count_lines

__all__ = ["collect"]


@register("quantify")
def collect(
    primary: Path,
    out_dir: Path,
    n_input_reads: Optional[int] = None,
    runtime_sec: float | None = None,
) -> dict:
    """QC collector for the quantify stage."""
    qc_start = time.time()

    isoform_tpm = primary
    gene_counts = primary.with_name(primary.name.replace("isoform.tpm", "gene.counts"))

    n_isoforms = max(count_lines(isoform_tpm) - 1, 0)
    n_genes = max(count_lines(gene_counts) - 1, 0) if gene_counts.exists() else 0

    tpm_vals: list[float] = []
    with open(isoform_tpm) as fh:
        for line in fh:
            if not line.strip() or line.startswith("isoform"):
                continue
            parts = line.rstrip().split("\t")[1:]
            for p in parts:
                try:
                    tpm_vals.append(float(p))
                except ValueError:
                    pass

    png_name = None
    if tpm_vals:
        plt.figure()
        plt.hist(tpm_vals, bins=50)
        plt.xlabel("TPM")
        plt.ylabel("count")
        plt.tight_layout()
        png_path = Path(out_dir) / "quantify_tpm_hist.png"
        plt.savefig(png_path, dpi=150)
        plt.close()
        png_name = png_path.name

    metrics = {
        "n_isoforms": n_isoforms,
        "n_genes": n_genes,
        "quantify_runtime_sec": round(runtime_sec, 2) if runtime_sec else None,
        "qc_runtime_sec": round(time.time() - qc_start, 2),
    }
    write_metrics(out_dir, "quantify", metrics)
    if png_name:
        with open(Path(out_dir) / "quantify_plot_manifest.json", "w") as fh:
            json.dump({"tpm": png_name}, fh, indent=2)
    return metrics
