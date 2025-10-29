#!/usr/bin/env python3
"""Placeholder QC generation for FLAIR regionalize runs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def _summarise_file(path: Path) -> dict[str, Any]:
    if path is None or not path.exists():
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
    parser = argparse.ArgumentParser(description="Generate placeholder QC for FLAIR regionalize runs.")
    parser.add_argument("--dataset", required=True, help="Dataset name for this regionalize run.")
    parser.add_argument("--align-idx", required=True, type=int, help="Align command index (1-based).")
    parser.add_argument("--region-idx", required=True, type=int, help="Regionalize command index (1-based).")
    parser.add_argument("--mode", required=True, help="Regionalize mode (e.g. region, all).")
    parser.add_argument(
        "--command-text",
        required=False,
        default="",
        help="Original regionalize command text for reference.",
    )
    parser.add_argument("--stdout", required=True, type=Path, help="Regionalize stdout file.")
    parser.add_argument("--stderr", required=True, type=Path, help="Regionalize stderr file.")
    parser.add_argument("--metadata", required=True, type=Path, help="Regionalize metadata JSON file.")
    parser.add_argument(
        "--region-details",
        required=False,
        type=Path,
        default=None,
        help="Optional region details TSV produced by regionalize.",
    )
    parser.add_argument("--output", required=True, type=Path, help="Output TSV path.")
    args = parser.parse_args()

    stdout_stats = _summarise_file(args.stdout)
    stderr_stats = _summarise_file(args.stderr)
    region_stats = _summarise_file(args.region_details)
    metadata = _load_metadata(args.metadata)

    rows: list[tuple[str, Any]] = [
        ("dataset", args.dataset),
        ("align_index", args.align_idx),
        ("region_index", args.region_idx),
        ("mode", args.mode),
        ("metadata_present", "yes" if metadata else "no"),
        ("stdout_lines", stdout_stats["lines"]),
        ("stdout_bytes", stdout_stats["size"]),
        ("stderr_lines", stderr_stats["lines"]),
        ("stderr_bytes", stderr_stats["size"]),
        ("region_details_present", "yes" if region_stats["exists"] else "no"),
        ("region_details_bytes", region_stats["size"]),
    ]
    command_text = args.command_text.strip()
    if command_text:
        rows.append(("command_text_snippet", command_text[:80]))

    _write_rows(args.output, rows)


if __name__ == "__main__":
    main()
