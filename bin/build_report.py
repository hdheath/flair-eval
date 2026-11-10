#!/usr/bin/env python3
"""PLACEHOLDER owned by the 'report' stream (docs/CONTRACTS.md). Replace wholesale.

Contract usage:
  build_report.py --bundles DIR_OR_GLOB [...] --out report.html [--title T] [--baseline RUN_ID]
"""
import sys
if __name__ == "__main__":
    if "--help" in sys.argv or "-h" in sys.argv:
        print(__doc__); sys.exit(0)
    sys.exit("build_report.py: not implemented yet (placeholder, see docs/CONTRACTS.md)")
