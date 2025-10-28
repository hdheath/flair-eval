#!/usr/bin/env python3
"""Placeholder QC generation for FLAIR collapse runs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def _summarise_file(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"exists": False, "size": 0, "lines": 0}
    size = path.stat().st_size
    line_breaks = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            line_breaks += chunk.count(b"\n")
    has_trailing_newline = False
    if size > 0:
        with path.open("rb") as handle:
            handle.seek(-1, 2)
            has_trailing_newline = handle.read(1) == b"\n"
    lines = line_breaks if has_trailing_newline else line_breaks + (1 if size > 0 else 0)
    return {"exists": True, "size": size, "lines": lines}


def _load_metadata(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return {}


def _write_rows(output: Path, rows: list[tuple[str, Any]]) -> None:
    output.write_text(
        "metric\tvalue\n"
        + "\n".join(f"{metric}\t{value}" for metric, value in rows),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate placeholder QC for FLAIR collapse runs.")
    parser.add_argument("--dataset", required=True, help="Dataset name for this collapse run.")
    parser.add_argument("--align-idx", required=True, type=int, help="Align command index (1-based).")
    parser.add_argument("--region-idx", required=True, type=int, help="Regionalize command index (1-based). Use 0 for non-regionalized runs.")
    parser.add_argument("--mode", required=True, help="Source mode for collapse (region or all).")
    parser.add_argument("--region-tag", required=False, default="", help="Region tag used for naming collapse outputs.")
    parser.add_argument("--stdout", required=True, type=Path, help="Collapse stdout file.")
    parser.add_argument("--stderr", required=True, type=Path, help="Collapse stderr file.")
    parser.add_argument("--metadata", required=True, type=Path, help="Collapse metadata JSON file.")
    parser.add_argument("--isoforms-bed", required=True, type=Path, help="Collapse isoforms BED output.")
    parser.add_argument("--isoforms-gtf", required=True, type=Path, help="Collapse isoforms GTF output.")
    parser.add_argument("--isoforms-fa", required=True, type=Path, help="Collapse isoforms FASTA output.")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV path.")
    args = parser.parse_args()

    stdout_stats = _summarise_file(args.stdout)
    stderr_stats = _summarise_file(args.stderr)
    isoforms_bed_stats = _summarise_file(args.isoforms_bed)
    isoforms_gtf_stats = _summarise_file(args.isoforms_gtf)
    isoforms_fa_stats = _summarise_file(args.isoforms_fa)
    metadata = _load_metadata(args.metadata)

    rows: list[tuple[str, Any]] = [
        ("dataset", args.dataset),
        ("align_index", args.align_idx),
        ("region_index", args.region_idx),
        ("mode", args.mode),
        ("region_tag", args.region_tag or "NA"),
        ("metadata_present", "yes" if metadata else "no"),
        ("stdout_lines", stdout_stats["lines"]),
        ("stdout_bytes", stdout_stats["size"]),
        ("stderr_lines", stderr_stats["lines"]),
        ("stderr_bytes", stderr_stats["size"]),
        ("isoforms_bed_exists", "yes" if isoforms_bed_stats["exists"] else "no"),
        ("isoforms_bed_lines", isoforms_bed_stats["lines"]),
        ("isoforms_gtf_exists", "yes" if isoforms_gtf_stats["exists"] else "no"),
        ("isoforms_gtf_lines", isoforms_gtf_stats["lines"]),
        ("isoforms_fa_exists", "yes" if isoforms_fa_stats["exists"] else "no"),
        ("isoforms_fa_bytes", isoforms_fa_stats["size"]),
    ]

    _write_rows(args.output, rows)


if __name__ == "__main__":
    main()
