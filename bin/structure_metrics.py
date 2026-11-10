#!/usr/bin/env python3
"""PLACEHOLDER owned by the 'report' stream (docs/CONTRACTS.md). Replace wholesale.

Contract usage:
  structure_metrics.py --isoforms-bed B --isoforms-gtf T --gtf A --sample S --run-json RUN.json --env-json flair_env.json --out BUNDLE_DIR [--regions BED]
"""
import sys
if __name__ == "__main__":
    if "--help" in sys.argv or "-h" in sys.argv:
        print(__doc__); sys.exit(0)
    sys.exit("structure_metrics.py: not implemented yet (placeholder, see docs/CONTRACTS.md)")
