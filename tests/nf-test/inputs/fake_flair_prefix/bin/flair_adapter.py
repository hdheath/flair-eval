#!/usr/bin/env python3
"""Test double for bin/flair_adapter.py (nf-test only). `<step> --help` lists every optional
flag of CONTRACTS section 7 (A7-A9) so core's probing passes them; a run records its argv in
OUT/argv.json and writes the normalised outputs of CONTRACTS section 2."""
import json
import sys
from pathlib import Path

argv = sys.argv[1:]
if "--help" in argv:
    print("--reads --bam --genome --isoforms-bed --flair-bed --sample --platform --library")
    sys.exit(0)
step = argv[0]
opts, key = {}, None
for a in argv[1:]:
    if key is None and a.startswith("--") and "=" in a:
        k, v = a[2:].split("=", 1)
        opts[k] = v
    elif key is None and a.startswith("--"):
        key = a[2:]
    elif key is not None:
        opts[key] = a
        key = None
out = Path(opts["out"])
out.mkdir(parents=True, exist_ok=True)
(out / "argv.json").write_text(json.dumps({"argv": argv, "opts": opts}))
files = {
    "align": {"aligned.bam": "", "aligned.bam.bai": ""},
    "assemble": {"isoforms.bed": "chr1\t0\t20\tiso1\t1000\t+\t0\t20\t0\t1\t20,\t0,\n",
                 "isoforms.flair.bed": "chr1\t0\t20\tiso1\t1000\t+\t0\t20\t0\t1\t20,\t0,\n",
                 "isoforms.gtf": 'chr1\tFLAIR\texon\t1\t20\t.\t+\t.\tgene_id "g1"; transcript_id "iso1";\n',
                 "isoforms.fa": ">iso1\nACGTACGTACGTACGTACGT\n"},
    "quantify": {"counts.tsv": "isoform_id\tcount\tcpm\niso1\t1\t1000000\n"},
    "productivity": {"productivity.tsv": "isoform_id\tcds_start\tcds_end\tclass\niso1\t0\t20\tPRO\n"},
}[step]
for name, text in files.items():
    (out / name).write_text(text)
(out / "adapter_log.json").write_text(json.dumps({"step": step, "fake": True, "argv": argv}))
