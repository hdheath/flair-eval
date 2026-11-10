#!/usr/bin/env python3
"""nf-test double: `<step> --help` prints the REAL bin/flair_adapter.py help (so core's flag
probing sees exactly what the merged adapter offers); any other call is recorded by the fake
adapter next door (OUT/argv.json + contract outputs). No FLAIR is run."""
import os
import sys
from pathlib import Path

here = Path(__file__).resolve()
repo = here.parents[5]                      # tests/nf-test/inputs/real_help_prefix/bin -> repo
real = repo / "bin" / "flair_adapter.py"
fake = here.parents[2] / "fake_flair_prefix" / "bin" / "flair_adapter.py"
target = real if "--help" in sys.argv[1:] else fake
os.execv(sys.executable, [sys.executable, str(target)] + sys.argv[1:])
