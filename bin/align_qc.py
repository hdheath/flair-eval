#!/usr/bin/env python3
"""Placeholder QC generation for FLAIR align runs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def _summarise_file(path: Path) -> dict[str, Any]:
    """Return simple file metrics used for placeholder QC rows."""
    if not path.exists():
        return {"exists": False, "size": 0, "lines": 0}
    size = path.stat().st_size
    line_breaks = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):  # 1 MiB chunks
            line_breaks += chunk.count(b"\n")
    # Estimate final line if file does not end with newline.
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
    parser = argparse.ArgumentParser(description="Generate placeholder QC for FLAIR align runs.")
    parser.add_argument("--dataset", required=True, help="Dataset name for this align run.")
    parser.add_argument(
        "--command-idx",
        required=True,
        type=int,
        help="Align command index (1-based).",
    )
    parser.add_argument("--stdout", required=True, type=Path, help="Path to command stdout file.")
    parser.add_argument("--stderr", required=True, type=Path, help="Path to command stderr file.")
    parser.add_argument("--metadata", required=True, type=Path, help="Path to metadata JSON file.")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV path.")
    args = parser.parse_args()

    stdout_stats = _summarise_file(args.stdout)
    stderr_stats = _summarise_file(args.stderr)
    metadata = _load_metadata(args.metadata)

    rows: list[tuple[str, Any]] = [
        ("dataset", args.dataset),
        ("command_index", args.command_idx),
        ("metadata_present", "yes" if metadata else "no"),
        ("stdout_lines", stdout_stats["lines"]),
        ("stdout_bytes", stdout_stats["size"]),
        ("stderr_lines", stderr_stats["lines"]),
        ("stderr_bytes", stderr_stats["size"]),
    ]
    if metadata:
        command_rendered = metadata.get("command_rendered") or ""
        rows.append(("command_rendered_snippet", command_rendered[:80]))

    _write_rows(args.output, rows)


if __name__ == "__main__":
    main()
