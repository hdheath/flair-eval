"""Tests for the round-3 report (CONTRACTS section 11): the cross-run A33 computations, the
round-3 sections in the fixtures (A28-A32 shapes), and the viewer page itself, driven in node
through a small DOM stub over every tab, sample, baseline and layout."""

import copy
import json
import os
import re
import shutil
import subprocess
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
BIN = os.path.join(ROOT, "bin")
FIXTURES = os.path.join(ROOT, "tests", "fixtures", "bundles")
sys.path.insert(0, BIN)

from flair_eval_lib.annotation import Transcript  # noqa: E402
from flair_eval_lib.concordance import RunModels  # noqa: E402
from flair_eval_lib.crossrun import (  # noqa: E402
    cluster_jaccard, end_agreement, lost_by_decile, run_ends, tusco_incremental)

E = [-1000, -200, -50, -10, 0, 10, 50, 200, 1000]
WINDOWS = [0, 10, 25, 50, 100, 200, 300, 500]
ENDS_COLUMNS = ("isoform_id category subcategory strand tss tes d_reftss d_pas d_ref_tss "
                "d_ref_tes polya_motif down_a_frac").split()
NODE = shutil.which("node")

# Renders every tab for every sample x baseline x layout at desktop and phone width, plus one
# run alone and print mode; prints a JSON summary. Exits 1 on an exception, or on NaN /
# undefined written into an attribute or text.
DOM_STUB = r"""
'use strict';
var fs = require('fs');
var html = fs.readFileSync(process.argv[2], 'utf8');
var data = /<script type="application\/json" id="report-data">([\s\S]*?)<\/script>/.exec(html)[1];
var code = /<script>\n([\s\S]*?)<\/script>/.exec(html)[1];
var bad = [];
function Node(tag) {
  this.tagName = tag; this.childNodes = []; this.attrs = {}; this.listeners = {}; this.parentNode = null;
  this.style = {}; this.hidden = false; this._text = '';
}
Node.prototype.appendChild = function (c) { if (c.parentNode) c.parentNode.removeChild(c); c.parentNode = this; this.childNodes.push(c); return c; };
Node.prototype.removeChild = function (c) { var i = this.childNodes.indexOf(c); if (i >= 0) this.childNodes.splice(i, 1); c.parentNode = null; return c; };
Node.prototype.setAttribute = function (k, v) {
  v = String(v);
  if (/NaN|undefined|Infinity|\[object/.test(v)) bad.push(this.tagName + '[' + k + ']=' + v);
  this.attrs[k] = v;
};
Node.prototype.getAttribute = function (k) { return k in this.attrs ? this.attrs[k] : null; };
Node.prototype.removeAttribute = function (k) { delete this.attrs[k]; };
Node.prototype.addEventListener = function (k, f) { (this.listeners[k] = this.listeners[k] || []).push(f); };
Node.prototype.getBoundingClientRect = function () { return { width: global.__W, height: 100, left: 0, top: 0 }; };
Node.prototype.closest = function () { return null; };
Node.prototype.focus = function () {};
Object.defineProperty(Node.prototype, 'textContent', {
  get: function () { return this._text + this.childNodes.map(function (c) { return c.textContent; }).join(''); },
  set: function (v) {
    v = String(v);
    if (/NaN|undefined|\[object/.test(v)) bad.push(this.tagName + ' text=' + v);
    this.childNodes = []; this._text = v;
  }
});
Object.defineProperty(Node.prototype, 'offsetWidth', { get: function () { return 100; } });
Object.defineProperty(Node.prototype, 'offsetHeight', { get: function () { return 20; } });
var byId = {};
['report-data', 'controls', 'basenote', 'tabs', 'view', 'footer', 'subtitle', 'theme', 'print', 'tooltip'].forEach(function (id) { byId[id] = new Node('div'); });
byId['report-data']._text = data.replace(/\\u003c/g, '<');
var root = new Node('html');
global.document = {
  documentElement: root,
  getElementById: function (id) {
    if (byId[id]) return byId[id];
    var hit = null;
    (function walk(n) { if (hit) return; if (n.attrs && n.attrs.id === id) { hit = n; return; } (n.childNodes || []).forEach(walk); })(byId.tabs);
    return hit;
  },
  createElement: function (t) { return new Node(t); },
  createElementNS: function (ns, t) { return new Node(t); },
  createTextNode: function (t) { var n = new Node('#text'); n._text = String(t); if (/NaN|undefined/.test(n._text)) bad.push('text=' + t); return n; },
  addEventListener: function () {}
};
global.window = global;
global.__W = 1100;
window.innerWidth = 1100; window.innerHeight = 800;
window.addEventListener = function () {};
window.scrollTo = function () {};
window.print = function () {};
global.location = { hash: '' };
global.history = { replaceState: function (a, b, h) { location.hash = h; } };
global.localStorage = { getItem: function () { return null; }, setItem: function () {} };
global.setTimeout = setTimeout; global.clearTimeout = clearTimeout;
new Function(code)();
var R = window.flairEvalReport;
var D = JSON.parse(byId['report-data']._text);
var n = 0, errors = [];
function count(node, tag) { var c = node.tagName === tag ? 1 : 0; node.childNodes.forEach(function (k) { c += count(k, tag); }); return c; }
[1100, 360].forEach(function (W) {
  global.__W = W; window.innerWidth = W;
  D.samples.forEach(function (sample) {
    var runs = D.bundles.filter(function (b) { return b.sample === sample; }).map(function (b) { return b.run_id; });
    runs.forEach(function (base) {
      ['overlay', 'side'].forEach(function (layout) {
        R.views.forEach(function (tab) {
          R.state.sample = sample; R.state.base = base; R.state.layout = layout; R.state.tab = tab; R.state.hidden = {};
          try { R.render(); n++; } catch (e) { errors.push([W, sample, base, layout, tab, e.stack.split('\n').slice(0, 3).join(' | ')].join(' ')); }
        });
      });
    });
    // one run only shown, and print mode
    R.state.hidden = {}; runs.slice(1).forEach(function (r) { R.state.hidden[r] = 1; }); R.state.base = runs[0];
    R.views.forEach(function (tab) { R.state.tab = tab; try { R.render(); n++; } catch (e) { errors.push(['single', sample, tab, e.stack.split('\n').slice(0, 3).join(' | ')].join(' ')); } });
    R.state.hidden = {}; R.state.printing = true;
    try { R.render(); n++; } catch (e) { errors.push(['print', sample, e.stack.split('\n').slice(0, 3).join(' | ')].join(' ')); }
    R.state.printing = false;
  });
});
// what each tab shows, per sample, with the default baseline
var summary = {};
D.samples.forEach(function (sample) {
  R.state.sample = sample; R.state.base = D.baseline; R.state.layout = 'overlay'; R.state.hidden = {};
  summary[sample] = {};
  R.views.forEach(function (tab) {
    R.state.tab = tab; R.render();
    var cards = [], off = [];
    (function walk(nd) {
      if (nd.tagName === 'section') {
        var h = null;
        (function f(x) { if (h) return; if (x.tagName === 'h3') { h = x.textContent; return; } x.childNodes.forEach(f); })(nd);
        if ((nd.attrs['class'] || '').indexOf('off') >= 0) off.push([h, nd.childNodes[1] ? nd.childNodes[1].textContent : '']);
        else cards.push(h);
      }
      nd.childNodes.forEach(walk);
    })(byId.view);
    summary[sample][tab] = { cards: cards, hidden: off, svgs: count(byId.view, 'svg'), text: byId.view.textContent };
  });
});
console.log(JSON.stringify({ renders: n, errors: errors.slice(0, 20), n_errors: errors.length, bad_attrs: bad.slice(0, 20), n_bad: bad.length, summary: summary }, null, 1));
process.exitCode = errors.length || bad.length ? 1 : 0;
"""


def T(name, exons, strand="+", chrom="chr1"):
    return Transcript(name, "", chrom, strand, tuple(exons))


def run(rid, isoforms, info=None):
    return RunModels(rid, isoforms, info or {t.tid: ("FSM", "") for t in isoforms})


def build(*args):
    return subprocess.run([sys.executable, os.path.join(BIN, "build_report.py"), *map(str, args)],
                          capture_output=True, text=True)


def embedded(path):
    html = open(path, encoding="utf-8").read()
    return json.loads(re.search(r'id="report-data">(.*?)</script>', html, re.S).group(1))


def drive(tmp_path, html):
    """Run the page's script in node over every view; return the JSON summary."""
    stub = tmp_path / "dom_stub.js"
    stub.write_text(DOM_STUB)
    res = subprocess.run([NODE, str(stub), str(html)], capture_output=True, text=True, timeout=300)
    assert res.stdout, res.stderr
    out = json.loads(res.stdout)
    assert res.returncode == 0 and out["n_errors"] == 0 and out["n_bad"] == 0, \
        json.dumps({k: out[k] for k in ("errors", "bad_attrs")}, indent=1) + res.stderr
    return out


def fixture_dirs():
    for sample in sorted(os.listdir(FIXTURES)):
        sd = os.path.join(FIXTURES, sample)
        if os.path.isdir(sd):
            for rid in sorted(os.listdir(sd)):
                yield sample, rid, os.path.join(sd, rid)


# --------------------------------------------------------------------- A33: end clusters
def test_run_ends_are_strand_aware():
    e = run_ends(run("a", [T("p", [(100, 200), (300, 400)], "+"), T("m", [(100, 200)], "-"), T("u", [(5, 9)], ".")]))
    assert e["tss"] == {("chr1", "+"): [100], ("chr1", "-"): [199]}
    assert e["tes"] == {("chr1", "+"): [399], ("chr1", "-"): [100]}


def test_cluster_jaccard_single_linkage_within_window():
    a = {("chr1", "+"): [100, 1000]}
    b = {("chr1", "+"): [150, 2000]}           # 100~150 link (50 bp); 1000 and 2000 alone
    assert cluster_jaccard(a, b) == {"shared": 1, "only_a": 1, "only_b": 1, "jaccard": 0.3333}
    b2 = {("chr1", "+"): [151]}                # 51 bp: no link
    assert cluster_jaccard({("chr1", "+"): [100]}, b2)["shared"] == 0
    # chaining: a 100, b 140, a 180 is one cluster (single linkage), shared
    assert cluster_jaccard({("chr1", "+"): [100, 180]}, {("chr1", "+"): [140]})["shared"] == 1
    # strand matters
    assert cluster_jaccard({("chr1", "+"): [100]}, {("chr1", "-"): [100]})["jaccard"] == 0


def test_end_agreement_pairs():
    a = run("a", [T("a1", [(100, 200), (300, 400)])])
    b = run("b", [T("b1", [(120, 200), (300, 900)])])
    c = end_agreement([a, b, None])
    assert c["runs"] == ["a", "b"] and c["window"] == 50
    p = c["pairs"][0]
    assert p["tss"]["jaccard"] == 1.0 and p["tes"]["jaccard"] == 0.0
    assert end_agreement([a]) is None


# --------------------------------------------------------------------- A33: TUSCO incremental
def test_tusco_incremental_rules():
    per = {"base": {"g1": "FN", "g2": "TP", "g3": "PTP", "g4": "FP", "g5": "TP"},
           "run": {"g1": "PTP", "g2": "FN", "g3": "TP", "g4": "FN", "g5": "FP"}}
    t = tusco_incremental(per, {"g1": "GENE1"})
    p = next(x for x in t["pairs"] if x["base"] == "base" and x["run"] == "run")
    assert (p["gained"], p["both"], p["lost"]) == (1, 1, 1)       # FP -> FN is neither
    assert p["gained_genes"] == ["GENE1"] and p["lost_genes"] == ["g2"]
    q = next(x for x in t["pairs"] if x["base"] == "run")
    assert (q["gained"], q["lost"]) == (1, 1)
    assert tusco_incremental({"base": per["base"]}, {}) is None


# --------------------------------------------------------------------- A33: lost genes by decile
def test_lost_by_decile():
    iso, info, cpm = [], {}, {}
    for i in range(20):
        t = T(f"i{i}", [(i * 1000, i * 1000 + 100)])
        iso.append(t)
        info[t.tid] = ("FSM", f"G{i:02d}")
        cpm[t.tid] = float(i + 1)                 # G00 lowest ... G19 highest
    base = RunModels("base", iso, info)
    other = RunModels("other", iso[2:], {k: v for k, v in info.items() if k not in ("i0", "i1")})
    third = RunModels("third", iso[:-2], {k: v for k, v in info.items() if k not in ("i18", "i19")})
    d = lost_by_decile([base, other, third], {"base": cpm})
    b = d["by_base"]["base"]
    assert d["deciles"] == 10 and b["genes"] == 20 and b["n"] == [2] * 10
    assert b["lost"]["other"] == [2] + [0] * 9          # the two lowest genes
    assert b["lost"]["third"] == [0] * 9 + [2]          # the two highest
    assert "other" not in d["by_base"]                  # no counts.tsv for the other runs
    assert lost_by_decile([base], {"base": cpm}) is None


def test_lost_by_decile_splits_multi_gene_cpm():
    t = T("f", [(0, 100), (200, 300)])
    base = RunModels("base", [t], {"f": ("fusion", "A,B")})
    other = RunModels("other", [], {})
    b = lost_by_decile([base, other], {"base": {"f": 10.0}})["by_base"]["base"]
    assert b["genes"] == 2 and b["lost"]["other"][0] == b["lost"]["other"][5] == 1 and b["cpm_max"][5] == 5.0


# --------------------------------------------------------------------- fixtures: A28-A32 shapes
def test_fixtures_carry_round3_sections():
    seen_tusco = 0
    for sample, rid, d in fixture_dirs():
        b = json.load(open(os.path.join(d, "bundle.json")))
        n = b["summary"]["isoforms"]
        e = b["sections"]["ends"]
        assert e["windows"] == WINDOWS
        assert set(e["refs"]) == {"tss", "pas"} and all(set(r) == {"file", "md5", "sites"} for r in e["refs"].values())
        for half in ("tss", "tes"):
            h = e[half]
            assert h["n"] == n and len(h["within"]) == len(WINDOWS)
            assert h["within"] == sorted(h["within"]) and h["within"][-1] <= n
            assert h["dist_hist"]["edges"] == E and len(h["dist_hist"]["counts"]) == len(E) + 1
        for k in ("fsm_tss_hist", "fsm_tes_hist"):
            assert e[k]["edges"] == E and len(e[k]["counts"]) == len(E) + 1
            assert sum(e[k]["counts"]) <= b["categories"]["FSM"]
        assert sum(c["n"] for c in e["by_category"].values()) == n
        assert all(set(c) == {"n", "tss_within_50", "tes_within_50"} for c in e["by_category"].values())
        assert sum(c["tss_within_50"] for c in e["by_category"].values()) == e["tss"]["within"][WINDOWS.index(50)]
        assert set(e["polya_motif"]) == {"n_with_motif", "n"}
        p = b["sections"]["priming"]
        assert set(p) == {"n", "risk", "down_a_hist"} and sum(p["down_a_hist"]["counts"]) == p["n"]
        assert p["down_a_hist"]["edges"][0] == 0 and p["down_a_hist"]["edges"][-1] == 1.0
        v = b["sections"]["variation"]
        assert set(v) == {"genes_multi", "ends_only", "splice_only", "both"}
        assert v["ends_only"] + v["splice_only"] + v["both"] == v["genes_multi"]
        t = b["sections"].get("tusco")
        if t:
            seen_tusco += 1
            assert set(t) == {"genes_total", "genes_in_regions", "TP", "PTP", "FP", "FN", "sensitivity", "precision", "per_gene"}
            assert t["TP"] + t["PTP"] + t["FP"] + t["FN"] == len(t["per_gene"]) == t["genes_in_regions"]
            assert set(t["per_gene"].values()) <= {"TP", "PTP", "FP", "FN"}
            assert 0 <= t["sensitivity"] <= 1 and 0 <= t["precision"] <= 1
        ends = open(os.path.join(d, "tables", "ends.tsv")).read().splitlines()
        assert ends[0].split("\t") == ENDS_COLUMNS and len(ends) - 1 == n
    assert seen_tusco == 5           # wtc11_pacbio only: wtc11_ont tests the hidden TUSCO panels


# --------------------------------------------------------------------- report payload
def test_report_embeds_a33(tmp_path):
    out = tmp_path / "r.html"
    assert build("--bundles", FIXTURES, "--out", out).returncode == 0
    d = embedded(out)
    for sample in d["samples"]:
        x = d["by_sample"][sample]
        n = len(x["concordance"]["runs"])
        assert len(x["end_agreement"]["pairs"]) == n * (n - 1) // 2
        # lost genes by decile agree with the concordance gene lists, for every baseline
        g = x["concordance"]["genes_detected"]
        for base, bb in x["lost_deciles"]["by_base"].items():
            assert sum(bb["n"]) == bb["genes"]
            for r, row in bb["lost"].items():
                assert sum(row) == len(set(g["runs"][base]) - set(g["runs"][r]))
    assert "tusco_incremental" in d["by_sample"]["wtc11_pacbio"]
    assert "tusco_incremental" not in d["by_sample"]["wtc11_ont"]
    inc = d["by_sample"]["wtc11_pacbio"]["tusco_incremental"]
    assert len(inc["pairs"]) == 5 * 4                         # every ordered pair


# --------------------------------------------------------------------- the page in node
needs_node = pytest.mark.skipif(NODE is None, reason="node not installed")


@needs_node
def test_viewer_renders_every_view(tmp_path):
    out = tmp_path / "r.html"
    assert build("--bundles", FIXTURES, "--out", out).returncode == 0
    s = drive(tmp_path, out)["summary"]
    assert set(s) == {"wtc11_ont", "wtc11_pacbio"}
    for sample, tabs in s.items():
        assert list(tabs) == ["diff", "structure", "recovery", "diversity", "agreement", "ends", "reads", "compute"]
        assert {"All metrics", "Terms"} <= set(tabs["diff"]["cards"])
        assert {"Precision at 50 bp", "5′ precision as the window widens", "Distance to the nearest refTSS site",
                "FSM 3′ end vs the reference transcript", "PolyA signal before the 3′ end",
                "End precision by category", "Internal-priming risk"} <= set(tabs["ends"]["cards"])
        assert not any("recall" in c.lower() for c in tabs["ends"]["cards"])       # precision only (A29)
        assert "End precision by category" not in tabs["structure"]["cards"]          # moved to Ends
        assert {"How genes vary"} <= set(tabs["diversity"]["cards"])
        assert {"5′ and 3′ end agreement", "Expression on shared isoforms"} <= set(tabs["agreement"]["cards"])
        assert {"Who are the lost genes?", "Genes gained and lost vs baseline"} <= set(tabs["recovery"]["cards"])
        assert "Not like for like" in tabs["reads"]["text"] and "genome alignment (BAM)" in tabs["reads"]["text"]  # A21
        assert tabs["reads"]["text"].count("Not like for like") == 1                 # once per tab
        assert "A21" not in tabs["reads"]["text"] and "A21" not in tabs["agreement"]["text"]
        assert tabs["ends"]["svgs"] >= 6
    assert "TUSCO" in s["wtc11_pacbio"]["recovery"]["cards"]
    assert "TUSCO genes vs baseline" in s["wtc11_pacbio"]["agreement"]["cards"]
    hidden = dict(s["wtc11_ont"]["recovery"]["hidden"])
    assert hidden["TUSCO"] == "Hidden: needs the TUSCO gene list (--genome hg38, or --tusco)."  # one-line reason (A34)


@needs_node
def test_old_bundles_still_render(tmp_path):
    """Bundles from before round 3 (no ends/tusco/priming/variation, no isoforms.bed beside them)."""
    for sample, rid, d in fixture_dirs():
        b = json.load(open(os.path.join(d, "bundle.json")))
        for k in ("ends", "tusco", "priming", "variation"):
            b["sections"].pop(k, None)
        dst = tmp_path / "old" / sample / rid
        os.makedirs(dst)
        json.dump(b, open(dst / "bundle.json", "w"))
        if sample == "wtc11_pacbio":                    # one sample keeps the files beside it
            for f in ("isoforms.bed", "counts.tsv"):
                shutil.copy(os.path.join(d, f), dst / f)
            shutil.copytree(os.path.join(d, "tables"), dst / "tables")
    out = tmp_path / "r.html"
    assert build("--bundles", tmp_path / "old", "--out", out).returncode == 0
    s = drive(tmp_path, out)["summary"]
    ends = s["wtc11_pacbio"]["ends"]
    assert ends["cards"] == [] and ends["hidden"][0][1].startswith("Hidden: needs refTSS or PolyA_DB sites (--genome hg38, or --ref_tss")
    assert dict(ends["hidden"])["Internal-priming risk"].startswith("Hidden: not in these results")
    assert "Isoforms shared between runs (Jaccard)" in s["wtc11_pacbio"]["agreement"]["cards"]
    ont = dict(s["wtc11_ont"]["agreement"]["hidden"])                                # no isoforms.bed
    assert "Isoforms shared between runs (Jaccard)" in ont and "5′ and 3′ end agreement" in ont
    assert "Isoforms by structural category" in s["wtc11_ont"]["structure"]["cards"]


@needs_node
def test_single_run_sample_renders(tmp_path):
    sample, rid, d = next(fixture_dirs())
    dst = tmp_path / "one" / rid
    shutil.copytree(d, dst)
    out = tmp_path / "r.html"
    res = build("--bundles", tmp_path / "one", "--out", out)
    assert res.returncode == 0 and "1 bundle (1 run, 1 sample)" in res.stderr
    s = drive(tmp_path, out)["summary"][sample]
    assert "Only one run on this sample" in s["diff"]["text"]
    assert any(h[0] == "Isoforms shared between runs (Jaccard)" for h in s["agreement"]["hidden"])


@needs_node
def test_neither_genes_get_their_own_part(tmp_path):
    """A35: genes_multi includes `neither`; the stacked bar shows it rather than dropping it."""
    for sample, rid, d in fixture_dirs():
        if sample != "wtc11_pacbio":
            continue
        b = json.load(open(os.path.join(d, "bundle.json")))
        v = b["sections"]["variation"]
        v["neither"] = 3
        v["genes_multi"] += 3
        dst = tmp_path / "b" / sample / rid
        os.makedirs(dst)
        json.dump(b, open(dst / "bundle.json", "w"))
    out = tmp_path / "r.html"
    assert build("--bundles", tmp_path / "b", "--out", out).returncode == 0
    assert "Neither (same splicing, ends within 50 bp)" in drive(tmp_path, out)["summary"]["wtc11_pacbio"]["diversity"]["text"]


def test_colour_slots_do_not_follow_the_in_page_baseline():
    """A34: colour = DATA.runs position, fixed for the page; the in-page baseline only reorders."""
    js = open(os.path.join(ROOT, "assets", "report", "report.js")).read()
    assert "function slot(r) { var i = RUN_ORDER.indexOf(r);" in js
    assert "RUN_ORDER = DATA.runs.slice()" in js


def test_page_shell_has_viewer_controls(tmp_path):
    out = tmp_path / "r.html"
    assert build("--bundles", FIXTURES, "--out", out).returncode == 0
    html = open(out).read()
    for needle in ('id="controls"', 'id="tabs"', 'role="tablist"', 'id="view"', 'id="basenote"',
                   "repeat(auto-fill", "grid-template-columns: repeat(2, minmax(0, 1fr))",
                   'prefers-color-scheme: dark', ':root[data-theme="dark"]'):
        assert needle in html, needle
    # no concept-only badges or out-of-scope views
    for gone in ("Illustrative", "badge real", "Plot map", "drawCompact", "metaplot"):
        assert gone not in html, gone
