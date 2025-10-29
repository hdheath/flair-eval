#!/usr/bin/env python3
"""
Utility to regionalize FLAIR alignment outputs.

Given an aligned BAM/BED pair and dataset metadata, slice the inputs to a
specific genomic window or pass the original files through unchanged when
`--all` is requested.
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import shutil
import subprocess
import sys
import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

LOG_FORMAT = "[regionalize] %(message)s"
logger = logging.getLogger("regionalize")

Region = Tuple[str, int, int]


def run_cmd(cmd: List[str]) -> None:
    """Run a command and raise an informative error on failure."""
    logger.info("Running: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"Command failed with exit code {exc.returncode}: {' '.join(cmd)}") from exc


def ensure_bam_index(bam_path: Path) -> Optional[Path]:
    """Ensure the BAM has an index, returning the index path when present."""
    candidates = [
        bam_path.with_suffix(".bam.bai"),
        bam_path.with_suffix(".bai"),
        bam_path.parent / f"{bam_path.name}.bai",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    run_cmd(["samtools", "index", str(bam_path)])
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def parse_region(arg: str) -> Region:
    match = re.match(r"^([A-Za-z0-9._-]+):(\d+)-(\d+)$", arg.strip())
    if not match:
        raise argparse.ArgumentTypeError(
            "Region must be formatted as chr:start-end (e.g., chr3:100-200)."
        )
    chrom, start, end = match.group(1), int(match.group(2)), int(match.group(3))
    if start > end:
        start, end = end, start
    return chrom, start, end


def write_region_details(path: Path, regions: Iterable[Region]) -> None:
    lines = ["chrom\tstart\tend\tspan_bp"]
    for chrom, start, end in regions:
        lines.append(f"{chrom}\t{start}\t{end}\t{end - start + 1}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_with_name(src: Path, dest_dir: Path, name: Optional[str] = None) -> Optional[str]:
    if not src or not src.exists():
        logger.warning("Optional file missing, skipping copy: %s", src)
        return None
    destination = dest_dir / (name or src.name)
    if destination.exists():
        try:
            if os.path.samefile(src, destination):
                logger.info("Destination already populated with upstream symlink; recreating materialized copy: %s", destination)
                destination.unlink()
        except FileNotFoundError:
            pass
    logger.info("Copying %s -> %s", src, destination)
    shutil.copy2(src, destination)
    return destination.name


def slice_tabular(
    src: Path,
    dest: Path,
    chrom: str,
    start: int,
    end: int,
    start_idx: int,
    end_idx: int,
    one_based: bool,
    keep_comments: bool = False,
) -> str:
    written = False
    with src.open("r", encoding="utf-8") as fin, dest.open("w", encoding="utf-8") as fout:
        for raw in fin:
            line = raw.rstrip("\n")
            if not line:
                continue
            if keep_comments and line.startswith("#"):
                fout.write(line + "\n")
                continue
            parts = line.split("\t")
            if len(parts) <= max(start_idx, end_idx):
                continue
            if parts[0] != chrom:
                continue
            try:
                s_val = int(parts[start_idx])
                e_val = int(parts[end_idx])
            except ValueError:
                continue
            if one_based:
                if s_val >= start and e_val <= end:
                    fout.write(line + "\n")
                    written = True
            else:
                if s_val >= start and e_val <= end:
                    fout.write(line + "\n")
                    written = True
    if not written:
        dest.write_text("", encoding="utf-8")
    return dest.name


def slice_bam(
    bam_path: Path,
    chrom: str,
    start: int,
    end: int,
    output_dir: Path,
) -> Tuple[str, Optional[str]]:
    tag = f"{chrom}_{start}_{end}"
    output_bam = output_dir / f"{tag}.bam"
    tmp_prefix = output_dir / f"tmp_sort_{tag}"
    cmd = [
        "bash",
        "-lc",
        f"samtools view -b '{bam_path}' '{chrom}:{start}-{end}' | "
        f"samtools sort -o '{output_bam}' -T '{tmp_prefix}'",
    ]
    run_cmd(cmd)
    index_path = ensure_bam_index(output_bam)
    return output_bam.name, index_path.name if index_path else None


def run_fasta_slice(genome: Path, chrom: str, start: int, end: int, output_dir: Path) -> Optional[str]:
    destination = output_dir / f"{chrom}_{start}_{end}.fa"
    cmd = ["samtools", "faidx", str(genome), f"{chrom}:{start}-{end}"]
    try:
        with destination.open("w", encoding="utf-8") as fout:
            logger.info("Extracting FASTA slice -> %s", destination)
            subprocess.run(cmd, check=True, stdout=fout)
    except subprocess.CalledProcessError:
        logger.warning("samtools faidx failed; creating empty FASTA for %s", destination)
        destination.write_text("", encoding="utf-8")
    return destination.name


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Regionalize FLAIR align outputs.")
    parser.add_argument("--dataset", required=True, help="Dataset name.")
    parser.add_argument("--align-index", type=int, required=True, help="Align command index (1-based).")
    parser.add_argument("--command-index", type=int, required=True, help="Regionalize command index (1-based).")
    parser.add_argument("--command-text", required=True, help="Original regionalize command text.")
    parser.add_argument("--conda-env-label", required=True, help="Condensed conda environment label.")
    parser.add_argument("--output-dir", required=True, help="Directory to materialize outputs.")
    parser.add_argument("--align-bam", required=True, help="Aligned BAM from flair align.")
    parser.add_argument("--align-bed", required=True, help="Aligned BED from flair align.")
    parser.add_argument("--align-metadata", required=False, help="Align command metadata JSON.")
    parser.add_argument("--gtf", required=True, help="Reference GTF file.")
    parser.add_argument("--genome", required=False, help="Reference genome FASTA (optional).")
    parser.add_argument("--junctions", required=False, help="STAR SJ.out.tab or equivalent (optional).")
    parser.add_argument("--cage", required=False, help="CAGE peak BED (optional).")
    parser.add_argument("--quantseq", required=False, help="polyA peak BED (optional).")
    parser.add_argument("--target-regions", required=False, help="Target regions BED/TSV (optional).")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--all", action="store_true", help="Bypass regionalization; return original files.")
    group.add_argument("--region", metavar="CHR:START-END", help="Slice to the specified region.")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(stream=sys.stderr, level=logging.INFO, format=LOG_FORMAT)

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    align_bam = Path(args.align_bam).resolve()
    align_bed = Path(args.align_bed).resolve()
    if not align_bam.exists():
        parser.error(f"Align BAM not found: {align_bam}")
    if not align_bed.exists():
        parser.error(f"Align BED not found: {align_bed}")

    align_metadata = Path(args.align_metadata).resolve() if args.align_metadata else None
    align_index_path = ensure_bam_index(align_bam)

    gtf_path = Path(args.gtf).resolve()
    if not gtf_path.exists():
        parser.error(f"GTF not found: {gtf_path}")

    optional_inputs: Dict[str, Optional[Path]] = {
        "junctions": Path(args.junctions).resolve() if args.junctions else None,
        "cage": Path(args.cage).resolve() if args.cage else None,
        "quantseq": Path(args.quantseq).resolve() if args.quantseq else None,
        "target_regions": Path(args.target_regions).resolve() if args.target_regions else None,
    }
    genome_path = Path(args.genome).resolve() if args.genome else None

    outputs: List[str] = []
    if align_metadata and align_metadata.exists():
        copied = copy_with_name(align_metadata, output_dir, "align_command_metadata.json")
        if copied:
            outputs.append(copied)

    if args.all:
        tag = f"align{args.align_index}_region{args.command_index}_all"
        logger.info("Bypass requested (--all); copying original alignment artifacts.")
        bam_name = f"{tag}.bam"
        copied_bam = copy_with_name(align_bam, output_dir, bam_name)
        if copied_bam:
            outputs.append(copied_bam)
        if align_index_path and align_index_path.exists():
            bai_name = f"{bam_name}.bai"
            copied_index = copy_with_name(align_index_path, output_dir, bai_name)
            if copied_index:
                outputs.append(copied_index)
        copied_bed = copy_with_name(align_bed, output_dir, f"{tag}.bed")
        if copied_bed:
            outputs.append(copied_bed)
        # Always provide reference annotations for downstream expectations
        gtf_copied = copy_with_name(gtf_path, output_dir, f"{tag}.gtf")
        if gtf_copied:
            outputs.append(gtf_copied)
        if genome_path and genome_path.exists():
            genome_copied = copy_with_name(genome_path, output_dir, f"{tag}.fa")
            if genome_copied:
                outputs.append(genome_copied)

        for key, path in optional_inputs.items():
            if not path:
                continue
            opt_name = f"{tag}_{path.name}"
            copied = copy_with_name(path, output_dir, opt_name)
            if copied:
                outputs.append(copied)
        details = output_dir / "region_details.tsv"
        write_region_details(details, [("ALL", 0, 0)])
        outputs.append(details.name)
        region_info = None
        mode = "all"
    else:
        chrom, start, end = parse_region(args.region)
        region_info = {"chrom": chrom, "start": start, "end": end}
        mode = "region"
        details = output_dir / "region_details.tsv"
        write_region_details(details, [(chrom, start, end)])
        outputs.append(details.name)

        bam_name, bai_name = slice_bam(align_bam, chrom, start, end, output_dir)
        outputs.append(bam_name)
        if bai_name:
            outputs.append(bai_name)

        bed_name = slice_tabular(
            align_bed,
            output_dir / f"{chrom}_{start}_{end}.bed",
            chrom,
            start,
            end,
            start_idx=1,
            end_idx=2,
            one_based=False,
        )
        outputs.append(bed_name)

        gtf_name = slice_tabular(
            gtf_path,
            output_dir / f"{chrom}_{start}_{end}.gtf",
            chrom,
            start,
            end,
            start_idx=3,
            end_idx=4,
            one_based=True,
            keep_comments=True,
        )
        outputs.append(gtf_name)

        if genome_path and genome_path.exists():
            fasta_name = run_fasta_slice(genome_path, chrom, start, end, output_dir)
            if fasta_name:
                outputs.append(fasta_name)

        junctions_path = optional_inputs.get("junctions")
        if junctions_path and junctions_path.exists():
            sj_name = slice_tabular(
                junctions_path,
                output_dir / f"{chrom}_{start}_{end}.SJ.out.tab",
                chrom,
                start,
                end,
                start_idx=1,
                end_idx=2,
                one_based=True,
            )
            outputs.append(sj_name)
        elif junctions_path:
            logger.warning("Junctions file not found, skipping: %s", junctions_path)

        for key in ("cage", "quantseq", "target_regions"):
            path = optional_inputs.get(key)
            if not path:
                continue
            if not path.exists():
                logger.warning("Optional file not found for %s, skipping: %s", key, path)
                continue
            dest_name = output_dir / f"{chrom}_{start}_{end}_{path.name}"
            written_name = slice_tabular(
                path,
                dest_name,
                chrom,
                start,
                end,
                start_idx=1,
                end_idx=2,
                one_based=False,
            )
            outputs.append(written_name)

    metadata = {
        "dataset": args.dataset,
        "align_index": args.align_index,
        "regionalize_index": args.command_index,
        "regionalize_command": args.command_text,
        "mode": mode,
        "region": region_info,
        "inputs": {
            "align_bam": align_bam.name,
            "align_bed": align_bed.name,
            "align_metadata": "align_command_metadata.json" if align_metadata else None,
            "gtf": gtf_path.name,
            "genome": genome_path.name if genome_path else None,
            "junctions": optional_inputs["junctions"].name if optional_inputs["junctions"] else None,
            "cage": optional_inputs["cage"].name if optional_inputs["cage"] else None,
            "quantseq": optional_inputs["quantseq"].name if optional_inputs["quantseq"] else None,
            "target_regions": optional_inputs["target_regions"].name if optional_inputs["target_regions"] else None,
        },
        "outputs": sorted(outputs),
        "runtime": {
            "tool": "regionalize.py",
            "conda_env_label": args.conda_env_label,
        },
    }

    metadata_path = output_dir / "regionalize_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    logger.info("Wrote metadata -> %s", metadata_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
