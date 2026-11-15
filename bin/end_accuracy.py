#!/usr/bin/env python3
"""Score one FLAIR run's isoform ends against its sample's end sets (CONTRACTS section 12).

Runs flair4_ted's end judge, ted_end_precision.py (vendored into
modules/end_accuracy/resources/usr/bin/, legacy columns unchanged), once per planned call,
and writes:

  OUT/end_accuracy.json        sections.end_accuracy, merged into bundle.json by BUNDLE
  OUT/tables/end_accuracy.tsv  one row per end set, then one row per 5' x 3' pair
  OUT/judge/<call>/            the judge's own outputs and log, kept for the record

--sets is a JSON list of end sets: {id, kind: truth|reference, end: tss|tes, assay, path,
source, provisional}. bin/flair_eval_lib/end_accuracy.py has the call plan and the rules.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flair_eval_lib.annotation import InputError  # noqa: E402
from flair_eval_lib.end_accuracy import (  # noqa: E402
    JUDGE, WINDOW, build_section, plan_calls, read_sets, regions_args, run_calls, write_table)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--isoforms-bed", required=True, help="FLAIR isoforms, BED12")
    p.add_argument("--gtf", required=True, help="reference annotation GTF (the judge needs it)")
    p.add_argument("--sets", required=True, help="end sets, JSON list")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--regions", help="BED of regions; passed to the judge as --region")
    p.add_argument("--judge", help=f"path to {JUDGE} (default: found on PATH)")
    p.add_argument("--window", type=int, default=WINDOW, help="bp window (default %(default)s)")
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = parse_args(argv)
    try:
        judge = a.judge or shutil.which(JUDGE)
        if not judge or not os.path.isfile(judge):
            raise InputError(f"{JUDGE} not found (give --judge, or put it on PATH)")
        sets = read_sets(a.sets)
        if not sets:
            raise InputError(f"--sets {a.sets}: no end sets")
        regions = regions_args(a.regions)
        calls = plan_calls(sets)
        os.makedirs(a.out, exist_ok=True)
        run_calls(judge, a.isoforms_bed, a.gtf, calls, os.path.join(a.out, "judge"), regions,
                  a.window)
        section = build_section(sets, calls, judge, a.window, regions)
    except InputError as e:
        print(f"end_accuracy.py: error: {e}", file=sys.stderr)
        return 2
    write_table(os.path.join(a.out, "tables", "end_accuracy.tsv"), section)
    with open(os.path.join(a.out, "end_accuracy.json"), "w") as fh:
        json.dump(section, fh, indent=1)
        fh.write("\n")
    for r in section["ends"]:
        pr = "" if r["precision"] is None else f"P {r['precision']:.3f}"
        rc = "" if r["recall"] is None else f" R {r['recall']:.3f}"
        print(f"end_accuracy.py: {r['end']} {r['set']}: {pr}{rc}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
