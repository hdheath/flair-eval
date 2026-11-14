#!/usr/bin/env python3
"""Build one self-contained report.html from any set of flair-eval bundles.

Every bundle.json under the given paths is read (a path may be a bundle directory, a folder
of bundles, a bundle.json file or a glob). When two bundles share (sample, run_id), the one
with the newest pipeline.created wins. The HTML inlines its CSS/JS and embeds the bundle
data as JSON, so it opens offline. Per-isoform tables are not embedded.

The baseline is --baseline when given. Without it (CONTRACTS A22) it is the first run of the
earliest launch: bundles are grouped into launches by pipeline.launch (A24), or, for bundles
without it, where run_index starts over; the run with the lowest run_index in the launch holding
the oldest pipeline.created wins, ties by id. With no
run_index at all, the first run id in name order is used (A5). The page says which rule picked
the baseline and lets the reader switch it. Runs are ordered baseline first, then by
bundle.json run_index, then by id (A17), and keep that colour throughout the report.

Cross-run data (CONTRACTS section 8 and A33) is computed here per sample when its inputs sit
beside each bundle.json: run-vs-run concordance from isoforms.bed (+ tables/isoforms.tsv), CPM
agreement from counts.tsv, TSS/TTS end-cluster agreement from isoforms.bed, TUSCO incremental
recall from sections.tusco.per_gene, and lost genes by baseline expression decile from
counts.tsv. Per-bundle sections (quant, productivity, ends, tusco, priming, variation) and the
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
from flair_eval_lib.crossrun import end_agreement, lost_by_decile, tusco_incremental  # noqa: E402
from flair_eval_lib.schema import load_schema, validate  # noqa: E402
from flair_eval_lib.sections import quant_input_from_dir  # noqa: E402

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
                                      "(default: first run of the earliest launch)")
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


def sample_sections(entries):
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
        tusco = {}
        for b, d in rows:
            ts = (b.get("sections") or {}).get("tusco") or {}
            pg = ts.get("per_gene")
            if isinstance(pg, dict) and pg:
                tusco[b["run_id"]] = pg
            tnames = ts.get("gene_names")
            if isinstance(tnames, dict):
                names.update({g: n for g, n in tnames.items() if n and g not in names})
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
        conc = concordance(models, names)
        if conc:
            entry["concordance"] = conc
        quant = quant_agreement(models, cpm)
        if quant:
            entry["quant"] = quant
        ends = end_agreement(models)
        if ends:
            entry["end_agreement"] = ends
        inc = tusco_incremental(tusco, names)
        if inc:
            entry["tusco_incremental"] = inc
        lost = lost_by_decile(models, cpm)
        if lost:
            entry["lost_deciles"] = lost
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


def fill_quant_input(entries):
    """Older bundles lack sections.quant.input (A21); derive it from the adapter log beside them."""
    for b, f in entries:
        q = (b.get("sections") or {}).get("quant")
        if isinstance(q, dict) and not q.get("input"):
            qi = quant_input_from_dir(os.path.dirname(f))
            if qi:
                q["input"] = qi
                q["input_inferred"] = True


def run_order(bundles, baseline):
    """A17: baseline first, then run_index, then id; runs without run_index after, by id."""
    index = {}
    for b in bundles:
        ri = b.get("run_index")
        if isinstance(ri, int) and not isinstance(ri, bool):
            index[b["run_id"]] = min(ri, index.get(b["run_id"], ri))
    ids = sorted({b["run_id"] for b in bundles},
                 key=lambda r: (r != baseline, r not in index, index.get(r, 0), r))
    return ids


def _run_index(b):
    ri = b.get("run_index")
    return ri if isinstance(ri, int) and not isinstance(ri, bool) else None


def _launch_id(b):
    v = (b.get("pipeline") or {}).get("launch")
    return v if isinstance(v, str) and v else None


def _guess_launches(bundles):
    """A22 fallback for bundles without pipeline.launch. pipeline.created is stamped per bundle,
    so a launch is recovered from time order: walking the bundles oldest first, a new launch
    starts when a (sample, run_index) pair repeats, since run_index counts runs.yaml from 0 in
    every launch. Bundles without run_index join the launch they fall in. Limit: a later launch
    whose run indexes do not repeat any of the earlier launch's is merged into it."""
    out, seen = [], None
    for b in sorted(bundles, key=lambda b: (_created_key(b), b["run_id"], b["sample"])):
        ri = _run_index(b)
        key = (b["sample"], ri)
        if out and (ri is None or key not in seen):
            out[-1].append(b)
        else:
            out.append([b])
            seen = set()
        if ri is not None:
            seen.add(key)
    return out


def launches(bundles):
    """Bundles grouped by pipeline launch, earliest first (CONTRACTS A24).

    Bundles with pipeline.launch (Nextflow's session id) are grouped by it; the rest by the A22
    heuristic. The earliest launch is the one whose oldest pipeline.created is oldest."""
    by_id = {}
    for b in bundles:
        lid = _launch_id(b)
        if lid is not None:
            by_id.setdefault(lid, []).append(b)
    groups = list(by_id.values()) + _guess_launches([b for b in bundles if _launch_id(b) is None])
    for g in groups:
        g.sort(key=lambda b: (_created_key(b), b["run_id"], b["sample"]))
    return sorted(groups, key=lambda g: (_created_key(g[0]), g[0]["run_id"], g[0]["sample"]))


def default_baseline(bundles):
    """(run_id, reason) for the A22 default baseline. reason: earliest_launch | name_order."""
    if not any(_run_index(b) is not None for b in bundles):
        return sorted({b["run_id"] for b in bundles})[0], "name_order"
    first = next(g for g in ([b for b in l if _run_index(b) is not None] for l in launches(bundles)) if g)
    best = min(first, key=lambda b: (_run_index(b), b["run_id"]))
    return best["run_id"], "earliest_launch"


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
        baseline, rule = default_baseline(bundles)
    else:
        rule = "given"
    n_launches = len(launches(bundles))
    fill_quant_input(entries)
    runs = run_order(bundles, baseline)
    ordered = sorted(bundles, key=lambda b: (samples.index(b["sample"]), runs.index(b["run_id"])))
    return {
        "title": title,
        "generated": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "tool_version": __version__,
        "baseline": baseline,
        "baseline_rule": {"rule": rule, "launches": n_launches},
        "runs": runs,
        "samples": samples,
        "bundles": [slim(b) for b in ordered],
        "by_sample": sample_sections(entries),
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


def _n(k, word):
    """'1 sample', '2 samples'."""
    return f"{k} {word}{'' if k == 1 else 's'}"


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
    print(f"build_report.py: {_n(len(entries), 'bundle')} ({_n(len(data['runs']), 'run')}, "
          f"{_n(len(data['samples']), 'sample')}) -> {args.out} "
          f"({os.path.getsize(args.out) / 1024:.0f} KB)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
