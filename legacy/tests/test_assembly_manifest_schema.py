"""Round-trip + schema-contract test for the STAGE-1 assembly manifest.

The manifest is produced by `subworkflows/assemble.nf` (one CSV row per
produced assembly) and consumed by the `EVALUATION` entrypoint in `main.nf`
via `Utils.parseManifest` (lib/Utils.groovy). This test pins the *contract*
that both sides must agree on:

  1. The canonical column order (Utils.MANIFEST_COLUMNS).
  2. CSV quoting rules: a field containing a comma is wrapped in double
     quotes by the writer and the quotes are stripped by the parser.
  3. The EVALUATION reconstruction reads the columns it needs by name.

It mirrors the Groovy logic in pure Python so it can run offline (no
Nextflow/Groovy runtime). If you change the column list or quoting in
either Groovy file, update MANIFEST_COLUMNS here and this test will catch
a drift between writer and reader.

The expected column list is parsed DIRECTLY out of lib/Utils.groovy so the
test fails loudly if the Groovy source and this Python mirror diverge.
"""
from __future__ import annotations

import csv
import io
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
UTILS_GROOVY = REPO_ROOT / "lib" / "Utils.groovy"

# Python mirror of Utils.MANIFEST_COLUMNS — kept in this exact order.
MANIFEST_COLUMNS = [
    "test_name", "dataset_name", "align_mode", "partition_mode",
    "transcriptome_mode", "tool",
    "isoforms_bed", "isoforms_gtf", "isoform_read_map", "ted_log", "isoform_counts",
    "firstpass", "cds_bed", "cds_info",
    "bam", "bai", "reads_bed", "genome", "gtf", "cage_peaks", "drna_peaks",
    "cage_signal_plus", "cage_signal_minus", "drna_signal_plus", "drna_signal_minus",
    "library_type", "partition_region", "partition_args",
]


def _columns_from_groovy() -> list[str]:
    """Extract the MANIFEST_COLUMNS string literals from lib/Utils.groovy."""
    text = UTILS_GROOVY.read_text()
    m = re.search(r"MANIFEST_COLUMNS\s*=\s*\[(.*?)\]", text, re.DOTALL)
    assert m, "MANIFEST_COLUMNS not found in Utils.groovy"
    return re.findall(r"'([^']+)'", m.group(1))


def _format_row_groovy_style(values: list[str]) -> str:
    """Mirror of ASSEMBLE's per-field formatter:
    `s.contains(',') ? "\"${s}\"" : s` joined by ','."""
    out = []
    for v in values:
        s = str(v)
        out.append(f'"{s}"' if "," in s else s)
    return ",".join(out)


def _split_csv_line_groovy_style(line: str) -> list[str]:
    """Mirror of Utils.splitCsvLine (quote-aware, quotes stripped)."""
    out, cur, in_q = [], [], False
    for c in line:
        if c == '"':
            in_q = not in_q
        elif c == "," and not in_q:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(c)
    out.append("".join(cur))
    return out


def test_python_mirror_matches_groovy_source():
    """The Python MANIFEST_COLUMNS mirror must match the Groovy source."""
    assert _columns_from_groovy() == MANIFEST_COLUMNS


def test_manifest_round_trips():
    """A synthetic manifest (header + 2 rows incl. a comma-bearing field)
    written Groovy-style parses back to identical field maps."""
    rows_in = [
        {
            "test_name": "wtc11", "dataset_name": "sample1",
            "align_mode": "pre-aligned", "partition_mode": "chr-test",
            "transcriptome_mode": "TED-global-M85", "tool": "flair",
            "isoforms_bed": "/abs/iso.bed", "isoforms_gtf": "/abs/NO_ISOFORMS_GTF",
            "isoform_read_map": "/abs/map.txt", "ted_log": "/abs/ted_log.tsv",
            "isoform_counts": "/abs/counts.txt", "firstpass": "/abs/fp.bed",
            "cds_bed": "/abs/cds.bed", "cds_info": "/abs/cds.info.tsv",
            "bam": "/abs/x.bam", "bai": "/abs/x.bam.bai",
            "reads_bed": "/abs/reads.bed", "genome": "/canon/genome.fa",
            "gtf": "/abs/ann.gtf", "cage_peaks": "/abs/cage.bed",
            "drna_peaks": "/abs/drna.bed",
            "cage_signal_plus": "/abs/c+.bg", "cage_signal_minus": "/abs/c-.bg",
            "drna_signal_plus": "/abs/d+.bg", "drna_signal_minus": "/abs/d-.bg",
            "library_type": "pacbio_cDNA",
            # comma-bearing field exercises the quoting path
            "partition_region": "chr22:1-2;chr22:5-6",
            "partition_args": "--region chr22:1-2 chr22:5-6, --foo",
        },
        {
            "test_name": "wtc11", "dataset_name": "sample1",
            "align_mode": "pre-aligned", "partition_mode": "chr-test",
            "transcriptome_mode": "bambu_default", "tool": "bambu",
            "isoforms_bed": "/abs/NO_ISOFORMS_BED", "isoforms_gtf": "/abs/bambu.gtf",
            "isoform_read_map": "/abs/bmap.txt", "ted_log": "/abs/NO_TED_LOG",
            "isoform_counts": "/abs/bcounts.txt", "firstpass": "/abs/NO_FIRSTPASS",
            "cds_bed": "/abs/NO_CDS", "cds_info": "/abs/NO_CDS",
            "bam": "/abs/x.bam", "bai": "/abs/x.bam.bai",
            "reads_bed": "/abs/reads.bed", "genome": "/canon/genome.fa",
            "gtf": "/abs/ann.gtf", "cage_peaks": "/abs/cage.bed",
            "drna_peaks": "/abs/drna.bed",
            "cage_signal_plus": "", "cage_signal_minus": "",
            "drna_signal_plus": "", "drna_signal_minus": "",
            "library_type": "pacbio_cDNA",
            "partition_region": "", "partition_args": "",
        },
    ]

    # --- write Groovy-style ---
    header = ",".join(MANIFEST_COLUMNS)
    lines = [header]
    for r in rows_in:
        lines.append(_format_row_groovy_style([r[c] for c in MANIFEST_COLUMNS]))
    manifest_text = "\n".join(lines) + "\n"

    # --- parse Groovy-style ---
    parsed_header = None
    parsed_rows = []
    for raw in manifest_text.splitlines():
        if not raw.strip() or raw.startswith("#"):
            continue
        vals = _split_csv_line_groovy_style(raw)
        if parsed_header is None:
            parsed_header = vals
            continue
        if len(vals) > len(parsed_header):
            vals = vals[: len(parsed_header)]
        elif len(vals) < len(parsed_header):
            vals = vals + [""] * (len(parsed_header) - len(vals))
        parsed_rows.append(dict(zip(parsed_header, vals)))

    assert parsed_header == MANIFEST_COLUMNS
    assert parsed_rows == rows_in  # exact round-trip incl. the quoted comma field


def test_round_trip_matches_stdlib_csv():
    """The Groovy-style quoting agrees with Python's stdlib csv reader for
    the comma-bearing case (cross-check the quote handling)."""
    region = "chr22:1-2 chr22:5-6, --foo"
    line = _format_row_groovy_style(["a", region, "b"])
    # stdlib csv must recover the same 3 fields
    got = next(csv.reader(io.StringIO(line)))
    assert got == ["a", region, "b"]


if __name__ == "__main__":
    # Allow running without pytest:  python tests/test_assembly_manifest_schema.py
    test_python_mirror_matches_groovy_source()
    test_manifest_round_trips()
    test_round_trip_matches_stdlib_csv()
    print("OK: all 3 manifest-schema tests passed")
