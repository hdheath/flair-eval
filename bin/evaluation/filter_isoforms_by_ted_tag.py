#!/usr/bin/env python3
"""Subset a FLAIR/TED `.isoforms.bed` by the per-end TED tag (column 31,
`ted_end_tag`, written by `flair transcriptome --ted_tag_ends`).

The tag is a comma-joined string of `5_prime_TED` and/or `3_prime_TED` marking
which END(s) of the isoform TED approved. This lets the nextflow eval pipeline
evaluate each end independently or paired FROM ONE TAGGED RUN:

  --mode 5prime  keep isoforms whose 5' (TSS) is TED-approved (5_prime_TED)
  --mode 3prime  keep isoforms whose 3' (TTS) is TED-approved (3_prime_TED)
  --mode paired  keep isoforms with BOTH (== today's default-filter output)
  --mode any     keep isoforms with at least one approved end
  --mode none    no filtering (pass-through)

A 30-column (untagged) BED — i.e. a run WITHOUT --ted_tag_ends — passes through
UNFILTERED for every mode, so wiring this into the pipeline is byte-identical for
existing (untagged) outputs.
"""
import sys
import argparse

TSS_TAG = "5_prime_TED"
TTS_TAG = "3_prime_TED"


def _keep(cols, mode):
    if mode == "none":
        return True
    if len(cols) <= 30:          # untagged BED (no ted_end_tag column) -> pass through
        return True
    parts = set(p for p in cols[30].split(",") if p)
    if mode == "5prime":
        return TSS_TAG in parts
    if mode == "3prime":
        return TTS_TAG in parts
    if mode == "paired":
        return TSS_TAG in parts and TTS_TAG in parts
    if mode == "any":
        return bool(parts)
    raise ValueError(f"unknown mode {mode}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-bed", required=True, help="input .isoforms.bed")
    ap.add_argument("--out-bed", required=True, help="output filtered .isoforms.bed")
    ap.add_argument("--mode", required=True,
                    choices=["5prime", "3prime", "paired", "any", "none"],
                    help="which end(s) must be TED-approved to keep an isoform")
    a = ap.parse_args()
    n_in = n_out = 0
    with open(a.in_bed) as fh, open(a.out_bed, "w") as out:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                out.write(line)
                continue
            n_in += 1
            cols = line.rstrip("\n").split("\t")
            if _keep(cols, a.mode):
                out.write(line)
                n_out += 1
    sys.stderr.write(f"filter_isoforms_by_ted_tag: mode={a.mode} kept {n_out}/{n_in} isoforms\n")


if __name__ == "__main__":
    main()
