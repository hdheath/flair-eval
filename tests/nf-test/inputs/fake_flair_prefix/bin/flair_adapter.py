#!/usr/bin/env python3
"""Test double for bin/flair_adapter.py (nf-test only). Records its argv in OUT/argv.json and
writes the normalised outputs of docs/CONTRACTS.md section 2 for each step."""
import json
import sys
from pathlib import Path

argv = sys.argv[1:]
step = argv[0]
opts, key = {}, None
for a in argv[1:]:
    if a.startswith("--") and key is None:
        key = a[2:]
    elif key is not None:
        opts[key] = a
        key = None
out = Path(opts["out"])
out.mkdir(parents=True, exist_ok=True)
(out / "argv.json").write_text(json.dumps(argv))
files = {
    "align": {"aligned.bam": "", "aligned.bam.bai": ""},
    "assemble": {"isoforms.bed": "chr1\t0\t20\tiso1\t1000\t+\t0\t20\t0\t1\t20,\t0,\n",
                 "isoforms.gtf": 'chr1\tFLAIR\texon\t1\t20\t.\t+\t.\tgene_id "g1"; transcript_id "iso1";\n',
                 "isoforms.fa": ">iso1\nACGTACGTACGTACGTACGT\n"},
    "quantify": {"counts.tsv": "isoform_id\tcount\tcpm\niso1\t1\t1000000\n"},
    "productivity": {"productivity.tsv": "isoform_id\tcds_start\tcds_end\tclass\niso1\t0\t20\tPRO\n"},
}[step]
for name, text in files.items():
    (out / name).write_text(text)
(out / "adapter_log.json").write_text(json.dumps({"step": step, "fake": True, "argv": argv}))
