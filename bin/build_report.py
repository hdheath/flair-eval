#!/usr/bin/env python3
"""Build one self-contained report.html from any set of flair-eval bundles.

Every bundle.json under the given paths is read (a path may be a bundle directory, a folder
of bundles, a bundle.json file or a glob). When two bundles share (sample, run_id), the one
with the newest pipeline.created wins. The HTML inlines its CSS/JS and embeds the bundle
data as JSON, so it opens offline. Per-isoform tables are not embedded.

The baseline is --baseline (core always passes the first run of runs.yaml). Without it, the
first run id in name order is used (CONTRACTS A5). Runs keep one colour throughout the report:
the baseline first, then the other runs in the order they are first seen.

Round-2 sections (CONTRACTS section 8) are computed here when their inputs sit beside each
bundle.json: run-vs-run concordance from isoforms.bed (+ tables/isoforms.tsv), and CPM
agreement between runs from counts.tsv. Per-bundle sections (quant, productivity) and the
compute block are passed through.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import glob
import html
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flair_eval_lib import __version__  # noqa: E402
from flair_eval_lib.concordance import (  # noqa: E402
    RunModels, concordance, quant_agreement, read_counts, read_gene_names)
from flair_eval_lib.schema import load_schema, validate  # noqa: E402

ASSET_DIR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                          "..", "assets", "report"))
# Large, per-input fields that the report never shows.
_DROP_INPUT_KEYS = ()


class ReportError(Exception):
    pass


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bundles", nargs="+", required=True, metavar="DIR_OR_GLOB",
                   help="bundle dirs, folders of bundles, bundle.json files or globs")
    p.add_argument("--out", required=True, help="output HTML path")
    p.add_argument("--title", default="FLAIR comparison", help="report heading")
    p.add_argument("--baseline", help="run id to compute differences against "
                                      "(default: first run)")
    return p.parse_args(argv)


def find_bundle_files(paths):
    """bundle.json paths in discovery order, without duplicates."""
    out, seen = [], set()
    for arg in paths:
        matches = sorted(glob.glob(arg, recursive=True)) if glob.has_magic(arg) else [arg]
        if not matches:
            raise ReportError(f"--bundles {arg!r}: no match")
        for m in matches:
            if os.path.isdir(m):
                found = []
                for root, dirs, files in os.walk(m):
                    dirs.sort()
                    if "bundle.json" in files:
                        found.append(os.path.join(root, "bundle.json"))
                if not found and not glob.has_magic(arg):
                    raise ReportError(f"--bundles {arg!r}: no bundle.json found under it")
            elif os.path.isfile(m):
                found = [m]
            else:
                raise ReportError(f"--bundles {arg!r}: no such file or directory")
            for f in found:
                real = os.path.realpath(f)
                if real not in seen:
                    seen.add(real)
                    out.append(f)
    if not out:
        raise ReportError("no bundle.json files found")
    return out


def _created_key(b):
    raw = str(b.get("pipeline", {}).get("created", ""))
    try:
        dt = _dt.datetime.fromisoformat(raw.replace("Z", "+00:00"))
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=_dt.timezone.utc)
        return dt.timestamp()
    except ValueError:
        return float("-inf")


def load_bundles(files):
    schema = load_schema()
    chosen = {}
    order = []
    for f in files:
        try:
            with open(f) as fh:
                b = json.load(fh)
        except (OSError, json.JSONDecodeError) as e:
            raise ReportError(f"{f}: cannot read bundle ({e})") from None
        errs = validate(b, schema)
        if errs:
            raise ReportError(f"{f}: not a valid bundle:\n  " + "\n  ".join(errs[:10]))
        key = (b["sample"], b["run_id"])
        if key not in chosen:
            order.append(key)
            chosen[key] = (b, f)
        elif _created_key(b) > _created_key(chosen[key][0]):
            chosen[key] = (b, f)
    return [chosen[k] for k in order]


def sample_sections(entries, baseline):
    """Section-8 data computed across the runs of each sample: {sample: {...}}."""
    out = {}
    samples = []
    for b, _ in entries:
        if b["sample"] not in samples:
            samples.append(b["sample"])
    for sample in samples:
        rows = [(b, os.path.dirname(f)) for b, f in entries if b["sample"] == sample]
        names = {}
        models = []
        cpm = {}
        for b, d in rows:
            names.update(read_gene_names(d))
            try:
                m = RunModels.from_bundle_dir(b["run_id"], d)
            except (ValueError, OSError) as e:
                print(f"build_report.py: warning: {d}: concordance skipped ({e})", file=sys.stderr)
                m = None
            models.append(m)
            counts = os.path.join(d, "counts.tsv")
            if m is not None and os.path.exists(counts):
                try:
                    cpm[b["run_id"]] = read_counts(counts)
                except (ValueError, OSError) as e:
                    print(f"build_report.py: warning: {counts}: skipped ({e})", file=sys.stderr)
        entry = {}
        conc = concordance(models, baseline, names)
        if conc:
            entry["concordance"] = conc
        quant = quant_agreement(models, cpm, baseline)
        if quant:
            entry["quant"] = quant
        if entry:
            out[sample] = entry
    return out


def slim(b):
    """The bundle as embedded: everything the report reads, nothing per-isoform."""
    out = dict(b)
    inputs = dict(out.get("inputs") or {})
    for k in _DROP_INPUT_KEYS:
        inputs.pop(k, None)
    out["inputs"] = inputs
    return out


def report_data(entries, title, baseline=None):
    bundles = [b for b, _ in entries]
    seen, samples = [], []
    for b in bundles:
        if b["run_id"] not in seen:
            seen.append(b["run_id"])
        if b["sample"] not in samples:
            samples.append(b["sample"])
    if baseline is not None and baseline not in seen:
        raise ReportError(f"--baseline {baseline!r} is not a run in these bundles "
                          f"(runs: {', '.join(seen)})")
    if baseline is None:
        baseline = sorted(seen)[0]          # A5: name order only when no baseline is given
    runs = [baseline] + [r for r in seen if r != baseline]
    ordered = sorted(bundles, key=lambda b: (samples.index(b["sample"]), runs.index(b["run_id"])))
    return {
        "title": title,
        "generated": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "tool_version": __version__,
        "baseline": baseline,
        "runs": runs,
        "samples": samples,
        "bundles": [slim(b) for b in ordered],
        "by_sample": sample_sections(entries, baseline),
    }


def render_html(data):
    with open(os.path.join(ASSET_DIR, "template.html")) as fh:
        tpl = fh.read()
    with open(os.path.join(ASSET_DIR, "report.css")) as fh:
        css = fh.read()
    with open(os.path.join(ASSET_DIR, "report.js")) as fh:
        js = fh.read()
    if "</script" in js.lower() or "</style" in css.lower():
        raise ReportError("report assets must not contain </script or </style")
    payload = json.dumps(data, separators=(",", ":"), ensure_ascii=False)
    # Keep the JSON inert inside <script>: no "</" or "<!--" sequences.
    payload = payload.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    out = tpl.replace("/*__CSS__*/", css).replace("/*__JS__*/", js)
    out = out.replace("__TITLE__", html.escape(data["title"]))
    return out.replace("__DATA__", payload)


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        files = find_bundle_files(args.bundles)
        entries = load_bundles(files)
        data = report_data(entries, args.title, args.baseline)
        page = render_html(data)
    except ReportError as e:
        print(f"build_report.py: error: {e}", file=sys.stderr)
        return 2
    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(page)
    print(f"build_report.py: {len(entries)} bundles ({len(data['runs'])} runs, "
          f"{len(data['samples'])} samples) -> {args.out} "
          f"({os.path.getsize(args.out) / 1024:.0f} KB)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
