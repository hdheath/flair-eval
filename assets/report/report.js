/* flair-eval report: the transcriptome viewer (CONTRACTS section 11).
 * Tabs: What differs, Structure, Recovery, Diversity, Agreement, Ends, Reads, Compute.
 * One sample at a time; every comparison is against the baseline, which the reader can switch.
 * Plain JS + inline SVG; no libraries and no network. Data: <script id="report-data">.
 * Text from the data is inserted with textContent, never innerHTML.
 * Colour means run: a run's --series-N slot is its position in DATA.runs (A17 order), fixed for
 * the page, so colours never move when the baseline changes. Shared / other is neutral grey. */
(function () {
  'use strict';

  var DATA = JSON.parse(document.getElementById('report-data').textContent);
  var NS = 'http://www.w3.org/2000/svg';

  // ------------------------------------------------------------------ runs, samples, state
  function prettyId(id) {
    var m;
    if ((m = /^flair-(.+)$/.exec(id))) return 'FLAIR ' + m[1];
    if ((m = /^env-(.+)$/.exec(id))) return m[1] + ' (env)';
    if ((m = /^ctr-(.+)$/.exec(id))) return m[1] + ' (container)';
    return id;
  }
  var RUN_ORDER = DATA.runs.slice();
  var SAMPLES = DATA.samples.slice();
  var BY = {};
  var LABEL = {};
  DATA.bundles.forEach(function (b) {
    (BY[b.sample] = BY[b.sample] || {})[b.run_id] = b;
    if (!LABEL[b.run_id]) LABEL[b.run_id] = b.label && b.label !== b.run_id ? b.label : prettyId(b.run_id);
  });
  var CROSS = DATA.by_sample || {};
  function slot(r) { var i = RUN_ORDER.indexOf(r); return i >= 0 && i < 8 ? i + 1 : 0; }
  function col(r) { var k = slot(r); return k ? 'var(--series-' + k + ')' : 'var(--series-other)'; }
  function lab(r) { return LABEL[r] || r; }
  function trunc(t, n) { return t.length > n ? t.slice(0, Math.max(1, n - 1)) + '…' : t; }

  var VIEW_IDS = ['diff', 'structure', 'recovery', 'diversity', 'agreement', 'ends', 'reads', 'compute'];
  var state = { tab: 'diff', sample: SAMPLES[0], base: DATA.baseline, hidden: {}, layout: 'overlay', printing: false };

  function sampleRuns() { var m = BY[state.sample] || {}; return RUN_ORDER.filter(function (r) { return m[r]; }); }
  function baseId() { var rs = sampleRuns(); return rs.indexOf(state.base) >= 0 ? state.base : rs[0]; }
  // The baseline always comes first, so it is drawn first everywhere.
  function sel() { var b = baseId(); return [b].concat(sampleRuns().filter(function (r) { return r !== b && !state.hidden[r]; })); }
  function others() { return sel().slice(1); }
  function B(r) { return (BY[state.sample] || {})[r]; }
  function XS() { return CROSS[state.sample] || {}; }

  // ------------------------------------------------------------------ small helpers
  function isNum(v) { return typeof v === 'number' && isFinite(v); }
  function get(o, path) { return path.split('.').reduce(function (x, k) { return x === null || x === undefined ? undefined : x[k]; }, o); }
  function sum(a) { return a.reduce(function (x, y) { return x + (isNum(y) ? y : 0); }, 0); }
  function maxOf(a) { var m = -Infinity; a.forEach(function (v) { if (isNum(v) && v > m) m = v; }); return m; }
  function pc(a, b) { return isNum(a) && isNum(b) && b > 0 ? 100 * a / b : null; }
  var nf = new Intl.NumberFormat('en-US');
  function N(v) { return isNum(v) ? nf.format(Math.round(v)) : '–'; }
  function pct(v, d) { return isNum(v) ? v.toFixed(d === undefined ? 1 : d) + '%' : '–'; }
  function f2(v) { return isNum(v) ? v.toFixed(2) : '–'; }
  function f3(v) { return isNum(v) ? v.toFixed(3) : '–'; }
  function niceMax(v) {
    if (!(v > 0)) return 1;
    var p = Math.pow(10, Math.floor(Math.log(v) / Math.LN10)), m = v / p;
    return (m <= 1 ? 1 : m <= 2 ? 2 : m <= 2.5 ? 2.5 : m <= 5 ? 5 : 10) * p;
  }
  // Axis with 3-6 round ticks: {max, step}; cap = upper limit (100 for percentages).
  function niceAxis(v, cap) {
    if (!(v > 0)) return { max: 1, step: 0.25 };
    var raw = v / 4, p = Math.pow(10, Math.floor(Math.log(raw) / Math.LN10)), m = raw / p;
    var step = (m <= 1 ? 1 : m <= 2 ? 2 : m <= 5 ? 5 : 10) * p, max = Math.ceil(v / step - 1e-9) * step;
    if (cap && max > cap) { max = cap; step = cap / 4; }
    return { max: max, step: step };
  }
  function tickPct(v) { var r = Math.round(v * 10) / 10; return (r % 1 ? r.toFixed(1) : String(r)) + '%'; }
  function seqIdx(f) { return Math.max(0, Math.min(7, Math.round(f * 7))); }
  function fmtWall(v) { return !isNum(v) ? '–' : v >= 120 ? (v / 60).toFixed(1) + ' h' : v >= 10 ? v.toFixed(0) + ' min' : v.toFixed(1) + ' min'; }
  function fmtCpu(v) { return !isNum(v) ? '–' : v < 1 ? (v * 60).toFixed(v * 60 >= 10 ? 0 : 1) + ' CPU-min' : v.toFixed(2) + ' CPU-h'; }
  function fmtMem(v) { return !isNum(v) ? '–' : v < 1 ? Math.round(v * 1024) + ' MB' : v.toFixed(v >= 10 ? 0 : 1) + ' GB'; }
  // TUSCO sensitivity/precision: fractions per the contract; a value above 1 is read as a percent.
  function asPct(v) { return isNum(v) ? (v <= 1 ? 100 * v : v) : null; }

  // ------------------------------------------------------------------ DOM helpers
  function el(tag, attrs, kids) {
    var e = document.createElement(tag);
    Object.keys(attrs || {}).forEach(function (k) {
      var v = attrs[k];
      if (v === null || v === undefined) return;
      if (k === 'text') e.textContent = v;
      else if (k.slice(0, 2) === 'on') e.addEventListener(k.slice(2), v);
      else e.setAttribute(k, v);
    });
    (kids || []).forEach(function (c) { if (c !== null && c !== undefined) e.appendChild(typeof c === 'string' ? document.createTextNode(c) : c); });
    return e;
  }
  function s(tag, attrs, parent) {
    var e = document.createElementNS(NS, tag);
    Object.keys(attrs || {}).forEach(function (k) {
      if (attrs[k] === null || attrs[k] === undefined) return;
      if (k === 'text') e.textContent = attrs[k]; else e.setAttribute(k, attrs[k]);
    });
    if (parent) parent.appendChild(e);
    return e;
  }
  function svgBox(box, W, H, label) {
    return s('svg', { width: W, height: H, viewBox: '0 0 ' + W + ' ' + H, class: 'chart', role: 'img', 'aria-label': label || null }, box);
  }
  function fill(c) { return 'fill:' + c; }
  function width(box) { return Math.max(240, Math.floor(box.getBoundingClientRect().width)); }
  function swatch(r) { return el('span', { class: 'sw', style: 'background:' + col(r) }); }
  function runTag(r, short, noBase) {
    return el('span', { class: 'runtag' }, [swatch(r), short ? trunc(lab(r), short === true ? 18 : short) : lab(r), r === baseId() && !noBase ? el('span', { class: 'basetag', text: 'baseline' }) : null]);
  }
  function segment(label, opts, value, onPick) {
    var g = el('div', { class: 'seg', role: 'group', 'aria-label': label });
    opts.forEach(function (o) { g.appendChild(el('button', { type: 'button', 'aria-pressed': String(o[0] === value), onclick: function () { onPick(o[0]); }, text: o[1] })); });
    return g;
  }
  function picker(label, opts, value, onPick) {
    var sl = el('select', { 'aria-label': label });
    opts.forEach(function (o) { var op = el('option', { value: o[0], text: o[1] }); if (o[0] === value) op.selected = true; sl.appendChild(op); });
    sl.addEventListener('change', function () { onPick(sl.value); });
    return sl;
  }

  // ------------------------------------------------------------------ cards
  function grid(parent, cls) { var g = el('div', { class: 'grid' + (cls ? ' ' + cls : '') }); parent.appendChild(g); return g; }
  function section(root, title, dek) {
    root.appendChild(el('div', { class: 'sec-head' }, [el('h2', { text: title }), dek ? el('p', { class: 'dek', text: dek }) : null]));
  }
  function card(parent, o) {
    var c = el('section', { class: 'card' + (o.wide ? ' wide' : '') });
    var tools = el('div', { class: 'card-tools no-print' });
    c.appendChild(el('div', { class: 'card-head' }, [el('div', { class: 'card-titles' }, [el('h3', { text: o.title }), o.dek ? el('p', { class: 'dek', text: o.dek }) : null]), tools]));
    var body = el('div', { class: 'card-body' });
    c.appendChild(body);
    body._card = { node: c, tools: tools, table: false };
    parent.appendChild(c);
    return body;
  }
  // A panel whose section is missing is hidden, with a one-line reason (A34).
  function off(parent, title, reason) {
    var c = el('section', { class: 'card off' }, [el('h3', { text: title }), el('p', { class: 'small', text: 'Hidden: ' + reason + '.' })]);
    parent.appendChild(c);
    return null;
  }
  function cardOf(n) { while (n && !n._card) n = n.parentNode; return n ? n._card : null; }
  // A Table toggle with the chart's numbers. Tables are hidden until the reader asks for them.
  function addTable(box, head, rows) {
    var c = cardOf(box);
    var t = el('div', { class: 'tbl tablewrap' }, [el('table', { class: 'data' }, [
      el('thead', {}, [el('tr', {}, head.map(function (h, i) { return el('th', { class: i ? 'num' : null }, [typeof h === 'string' ? h : h]); }))]),
      el('tbody', {}, rows.map(function (r) { return el('tr', {}, r.map(function (v, i) { return el('td', { class: i ? 'num' : null, text: v }); })); }))])]);
    box.appendChild(t);
    if (c && !c.table) {
      c.table = true;
      var btn = el('button', { type: 'button', class: 'btn small-btn', 'aria-pressed': 'false', text: 'Table' });
      btn.addEventListener('click', function () {
        var on = c.node.getAttribute('data-view') === 'table';
        if (on) c.node.removeAttribute('data-view'); else c.node.setAttribute('data-view', 'table');
        btn.setAttribute('aria-pressed', String(!on));
      });
      c.tools.appendChild(btn);
    }
  }
  function note(parent, text, cls) { if (text) parent.appendChild(el('p', { class: cls || 'small', text: text })); }
  function missingNote(box, all, ok, reason) {
    var miss = all.filter(function (r) { return ok.indexOf(r) < 0; });
    if (miss.length) note(box, 'Not shown for ' + miss.map(lab).join(', ') + ': ' + reason + '.', 'small missing');
  }
  // Runs (of the selected ones) for which f(run) gives data.
  function withData(runs, f) { return runs.filter(function (r) { var v = f(r); return v !== null && v !== undefined; }); }
  // A card drawn for the runs that have its data; hidden with a reason when none do.
  function panel(parent, o, f, reason, draw) {
    var rs = o.runs || sel(), ok = withData(rs, f);
    if (!ok.length) return off(parent, o.title, reason);
    var body = card(parent, o);
    draw(body, ok);
    missingNote(body, rs, ok, reason);
    return body;
  }
  function tiles(parent, list) {
    var t = el('div', { class: 'tiles' });
    list.forEach(function (x) {
      t.appendChild(el('div', { class: 'tile' }, [el('div', { class: 'tile-k' }, [x[0]]), el('div', { class: 'tile-v num', text: x[1] }), x[2] ? el('div', { class: 'tile-s', text: x[2] }) : null]));
    });
    parent.appendChild(t);
    return t;
  }

  // overlay = one chart, every run; side = one panel per run, on the same scale (A34)
  function multi(box, runs, draw) {
    var viz = el('div', { class: 'viz' });
    box.appendChild(viz);
    if (state.layout === 'side' && runs.length > 1) {
      var g = el('div', { class: 'multiples' });
      viz.appendChild(g);
      runs.forEach(function (r) {
        var m = el('div', { class: 'mult' }, [el('div', { class: 'mult-title' }, [runTag(r)])]);
        g.appendChild(m);
        draw(m, [r]);
      });
    } else draw(viz, runs);
  }

  // ------------------------------------------------------------------ charts
  // Horizontal grouped bars. o: cats, vals {run: [..]}, fmt, tip(run, i), max, runs, labelW, valW
  function hbars(box, o) {
    var runs = (o.runs || sel()).filter(function (r) { return o.vals[r]; });
    var all = [];
    runs.forEach(function (r) { all = all.concat(o.vals[r]); });
    var max = o.max || niceMax(maxOf(all));
    multi(box, runs, function (b, rs) {
      var W = width(b), lw = Math.min(o.labelW || 150, W * 0.42), vw = o.valW || 58;
      var bh = rs.length > 1 ? 9 : 13, gap = 2, gg = 11;
      var gh = rs.length * bh + (rs.length - 1) * gap;
      var H = o.cats.length * (gh + gg);
      var svg = svgBox(b, W, H, o.label);
      var x = function (v) { return Math.max(0, (W - lw - vw) * Math.min(1, v / max)); };
      o.cats.forEach(function (c, i) {
        var y0 = i * (gh + gg) + 2;
        s('text', { x: lw - 8, y: y0 + gh / 2, 'text-anchor': 'end', 'dominant-baseline': 'central', class: 'lab', text: c }, svg);
        rs.forEach(function (r, j) {
          var v = o.vals[r][i], y = y0 + j * (bh + gap);
          if (!isNum(v)) return;
          s('rect', { x: lw, y: y, width: Math.max(1, x(v)), height: bh, rx: 2, style: fill(col(r)),
            'data-tip': lab(r) + ' · ' + c + ': ' + (o.tip ? o.tip(r, i) : o.fmt(v)) }, svg);
          s('text', { x: lw + x(v) + 5, y: y + bh / 2, 'dominant-baseline': 'central', class: 'val', text: o.fmt(v) }, svg);
        });
      });
    });
    addTable(box, [''].concat(runs.map(lab)), o.cats.map(function (c, i) {
      return [c].concat(runs.map(function (r) { var v = o.vals[r][i]; return isNum(v) ? (o.tip ? o.tip(r, i) : o.fmt(v)) : '–'; }));
    }));
  }
  // Vertical grouped bars. o: bins (tooltip names), ticks (axis labels), vals {run: [..]},
  // yfmt (axis), fmt (value), tip(run, i), xlab, mark (bin index whose left edge gets a dashed line)
  function vbars(box, o) {
    var runs = (o.runs || sel()).filter(function (r) { return o.vals[r]; });
    var all = [];
    runs.forEach(function (r) { all = all.concat(o.vals[r]); });
    var ax = niceAxis(o.max || maxOf(all), o.cap), max = ax.max;
    multi(box, runs, function (b, rs) {
      var W = width(b), H = o.h || 200, ml = 42, mr = 6, mt = 8, mb = o.xlab ? 46 : 30;
      var svg = svgBox(b, W, H, o.label);
      var n = o.bins.length, bw = (W - ml - mr) / n, inner = bw * 0.8, per = inner / rs.length;
      var y = function (v) { return mt + (H - mt - mb) * (1 - Math.min(1, v / max)); };
      for (var v = 0; v <= max + 1e-9; v += ax.step) {
        s('line', { x1: ml, x2: W - mr, y1: y(v), y2: y(v), class: 'gl' }, svg);
        s('text', { x: ml - 6, y: y(v), 'text-anchor': 'end', 'dominant-baseline': 'central', class: 'tick', text: o.yfmt(v, max) }, svg);
      }
      var longest = maxOf((o.ticks || o.bins).map(function (t) { return String(t).length; }));
      var step = Math.max(1, Math.ceil((longest * 6.3 + 6) / bw));
      o.bins.forEach(function (name, i) {
        var x0 = ml + i * bw + (bw - inner) / 2;
        rs.forEach(function (r, j) {
          var v = o.vals[r][i];
          if (!isNum(v)) return;
          s('rect', { x: x0 + j * per + 0.5, y: y(v), width: Math.max(1, per - 1), height: Math.max(0, y(0) - y(v)), rx: 1.5, style: fill(col(r)),
            'data-tip': lab(r) + ' · ' + name + ': ' + (o.tip ? o.tip(r, i) : o.fmt(v)) }, svg);
        });
        if (i % step === 0) s('text', { x: ml + i * bw + bw / 2, y: H - mb + 15, 'text-anchor': 'middle', class: 'tick', text: o.ticks ? o.ticks[i] : name }, svg);
      });
      if (isNum(o.mark)) s('line', { x1: ml + o.mark * bw, x2: ml + o.mark * bw, y1: mt, y2: H - mb, class: 'refl' }, svg);
      s('line', { x1: ml, x2: W - mr, y1: y(0), y2: y(0), class: 'ax' }, svg);
      if (o.xlab) s('text', { x: ml + (W - ml - mr) / 2, y: H - 6, 'text-anchor': 'middle', class: 'axlab', text: o.xlab }, svg);
    });
    addTable(box, [o.xname || 'Bin'].concat(runs.map(lab)), o.bins.map(function (name, i) {
      return [name].concat(runs.map(function (r) { var v = o.vals[r][i]; return isNum(v) ? (o.tip ? o.tip(r, i) : o.fmt(v)) : '–'; }));
    }));
  }
  // Histogram shown as the share of each run's own total, so runs of different sizes compare.
  function vhist(box, o) {
    var runs = (o.runs || sel()).filter(function (r) { return o.counts[r]; });
    var share = {};
    runs.forEach(function (r) { var t = sum(o.counts[r]) || 1; share[r] = o.counts[r].map(function (c) { return 100 * (c || 0) / t; }); });
    var cfg = {};
    Object.keys(o).forEach(function (k) { cfg[k] = o[k]; });
    cfg.runs = runs;
    cfg.vals = share;
    cfg.fmt = function (v) { return pct(v); };
    cfg.yfmt = tickPct;
    cfg.cap = 100;
    cfg.tip = function (r, i) { return N(o.counts[r][i]) + ' (' + pct(share[r][i]) + ')'; };
    vbars(box, cfg);
  }
  // Line chart over categorical x positions. o: xs (labels), series(run) -> [y|null], ymax, yfmt,
  // xlab, ylab, ref (x index for a dashed line), tip(run, i)
  function lines(box, o) {
    var runs = o.runs || sel();
    multi(box, runs, function (b, rs) {
      var W = width(b), H = o.h || 220, ml = 44, mr = 14, mt = 10, mb = 44;
      var svg = svgBox(b, W, H, o.label);
      var n = o.xs.length;
      var X = function (i) { return ml + (W - ml - mr) * (n > 1 ? i / (n - 1) : 0.5); };
      var ax = niceAxis(o.ymax, o.cap);
      var Y = function (v) { return mt + (H - mt - mb) * (1 - Math.min(1, v / ax.max)); };
      for (var v = 0; v <= ax.max + 1e-9; v += ax.step) {
        s('line', { x1: ml, x2: W - mr, y1: Y(v), y2: Y(v), class: 'gl' }, svg);
        s('text', { x: ml - 6, y: Y(v), 'text-anchor': 'end', 'dominant-baseline': 'central', class: 'tick', text: o.yfmt(v) }, svg);
      }
      var step = (W - ml - mr) / Math.max(1, n - 1) < 26 ? 2 : 1;
      o.xs.forEach(function (x, i) { if (i % step === 0 || i === n - 1) s('text', { x: X(i), y: H - mb + 16, 'text-anchor': 'middle', class: 'tick', text: x }, svg); });
      if (isNum(o.ref)) s('line', { x1: X(o.ref), x2: X(o.ref), y1: mt, y2: H - mb, class: 'refl' }, svg);
      rs.forEach(function (r) {
        var ys = o.series(r), d = '';
        ys.forEach(function (y, i) { if (isNum(y)) d += (d ? 'L' : 'M') + X(i).toFixed(1) + ',' + Y(y).toFixed(1); });
        if (d) s('path', { d: d, class: 'ln', style: 'stroke:' + col(r) }, svg);
      });
      // dots: the baseline first, larger and outlined (A34)
      rs.forEach(function (r) {
        var isBase = r === baseId();
        o.series(r).forEach(function (y, i) {
          if (!isNum(y)) return;
          s('circle', { cx: X(i), cy: Y(y), r: isBase ? 5 : 3.5, class: isBase ? 'dot base-dot' : 'dot', style: fill(col(r)),
            'data-tip': lab(r) + ' · ' + o.xs[i] + ': ' + (o.tip ? o.tip(r, i) : o.yfmt(y)) }, svg);
        });
      });
      s('line', { x1: ml, x2: W - mr, y1: H - mb, y2: H - mb, class: 'ax' }, svg);
      s('text', { x: ml + (W - ml - mr) / 2, y: H - 6, 'text-anchor': 'middle', class: 'axlab', text: o.xlab }, svg);
      if (rs.length > 1) b.appendChild(el('div', { class: 'legend' }, rs.map(function (r) { return runTag(r, 22); })));
    });
    addTable(box, [o.xname || ''].concat(runs.map(lab)), o.xs.map(function (x, i) {
      return [x].concat(runs.map(function (r) { var y = o.series(r)[i]; return isNum(y) ? (o.tip ? o.tip(r, i) : o.yfmt(y)) : '–'; }));
    }));
  }
  // Stacked 100% bars, one row per run. parts [{label, seq}] (grey ramp: parts of a whole), vals(run) -> [n..]
  function stack100(box, o) {
    var runs = o.runs || sel(), viz = el('div', { class: 'viz' });
    box.appendChild(viz);
    var W = width(viz), lw = Math.min(150, W * 0.32), bh = 18, gap = 10;
    var H = runs.length * (bh + gap) + 30;
    var svg = svgBox(viz, W, H, o.label);
    runs.forEach(function (r, i) {
      var v = o.vals(r), t = sum(v) || 1, x = lw, y = i * (bh + gap);
      s('rect', { x: 0, y: y + 4, width: 10, height: 10, rx: 2, style: fill(col(r)) }, svg);
      s('text', { x: 16, y: y + bh / 2, 'dominant-baseline': 'central', class: 'lab', text: trunc(lab(r), Math.floor((lw - 20) / 6.5)) }, svg);
      v.forEach(function (n, j) {
        var w = (W - lw) * (n || 0) / t;
        s('rect', { x: x, y: y, width: Math.max(0, w - 1), height: bh, style: fill('var(--seq-' + o.parts[j].seq + ')'),
          'data-tip': lab(r) + ' · ' + o.parts[j].label + ': ' + N(n) + ' (' + pct(100 * n / t) + ')' }, svg);
        if (w > 34) s('text', { x: x + w / 2, y: y + bh / 2, 'text-anchor': 'middle', 'dominant-baseline': 'central', class: 'inbar',
          style: fill('var(--seq-on-' + o.parts[j].seq + ')'), text: Math.round(100 * n / t) + '%' }, svg);
        x += w;
      });
    });
    var lx = lw, ly = runs.length * (bh + gap) + 8;
    o.parts.forEach(function (p) {
      if (lx + 22 + p.label.length * 6.2 > W) { lx = lw; ly += 16; }
      s('rect', { x: lx, y: ly, width: 10, height: 10, rx: 2, style: fill('var(--seq-' + p.seq + ')') }, svg);
      s('text', { x: lx + 14, y: ly + 5, 'dominant-baseline': 'central', class: 'tick', text: p.label }, svg);
      lx += 22 + p.label.length * 6.2;
    });
    if (ly + 14 > H) { svg.setAttribute('height', ly + 14); svg.setAttribute('viewBox', '0 0 ' + W + ' ' + (ly + 14)); }
    addTable(box, [''].concat(o.parts.map(function (p) { return p.label; })), runs.map(function (r) {
      var v = o.vals(r), t = sum(v) || 1;
      return [lab(r)].concat(v.map(function (n) { return N(n) + ' (' + pct(100 * n / t) + ')'; }));
    }));
  }
  // Grouped bars where each group lists its own (colour, value) bars. groups: [{label, bars: [{run, v, tip}]}]
  function pairBars(box, groups, fmt, labelW) {
    var viz = el('div', { class: 'viz' });
    box.appendChild(viz);
    var W = width(viz), lw = Math.min(labelW || 170, W * 0.45), vw = 50, bh = 10, gap = 2, gg = 14;
    var max = niceMax(maxOf([].concat.apply([], groups.map(function (g) { return g.bars.map(function (x) { return x.v; }); }))));
    var H = sum(groups.map(function (g) { return g.bars.length * (bh + gap) - gap + gg; }));
    var svg = svgBox(viz, W, H), y0 = 2;
    groups.forEach(function (g) {
      var gh = g.bars.length * (bh + gap) - gap;
      s('text', { x: lw - 8, y: y0 + gh / 2, 'text-anchor': 'end', 'dominant-baseline': 'central', class: 'lab', text: g.label }, svg);
      g.bars.forEach(function (x, j) {
        var y = y0 + j * (bh + gap), w = Math.max(1, (W - lw - vw) * x.v / max);
        s('rect', { x: lw, y: y, width: w, height: bh, rx: 2, style: fill(col(x.run)), 'data-tip': x.tip }, svg);
        s('text', { x: lw + w + 5, y: y + bh / 2, 'dominant-baseline': 'central', class: 'val', text: fmt(x.v) }, svg);
      });
      y0 += gh + gg;
    });
    return viz;
  }
  // One bar per run: baseline only · shared · run only (shared is neutral grey).
  function sharedBars(box, rows, what) {
    var viz = el('div', { class: 'viz' });
    box.appendChild(viz);
    var W = width(viz), bh = 20, gap = 30, H = rows.length * (bh + gap) + 4, svg = svgBox(viz, W, H);
    var b = baseId();
    rows.forEach(function (row, i) {
      var r = row.run, p = row.p, t = p.only_a + p.shared + p.only_b || 1, y = i * (bh + gap) + 16, x = 0;
      s('text', { x: 0, y: y - 6, class: 'lab', text: trunc(lab(b), 22) + ' only · shared · ' + trunc(lab(r), 22) + ' only' }, svg);
      [[p.only_a, col(b), lab(b) + ' only'], [p.shared, 'var(--neutral-bar)', 'Shared'], [p.only_b, col(r), lab(r) + ' only']].forEach(function (seg) {
        var w = W * seg[0] / t;
        s('rect', { x: x, y: y, width: Math.max(0, w - 1.5), height: bh, rx: 2, style: fill(seg[1]), 'data-tip': seg[2] + ': ' + N(seg[0]) + ' ' + what + ' (' + pct(100 * seg[0] / t) + ')' }, svg);
        if (w > 48) s('text', { x: x + w / 2, y: y + bh / 2, 'text-anchor': 'middle', 'dominant-baseline': 'central', class: 'inbar onbar', text: N(seg[0]) }, svg);
        x += w;
      });
    });
    addTable(box, ['Run', lab(b) + ' only', 'Shared', 'Run only', 'Jaccard'], rows.map(function (row) {
      return [lab(row.run), N(row.p.only_a), N(row.p.shared), N(row.p.only_b), f2(row.p.jaccard)];
    }));
  }

  // ------------------------------------------------------------------ cross-run data, oriented to the baseline
  // Pairs are stored once per unordered pair; orient(p, a) returns it with a = the given run.
  function transpose(g) { return g[0].map(function (_, i) { return g.map(function (row) { return row[i]; }); }); }
  function orient(p, a) {
    if (!p || p.a === a) return p;
    var q = {};
    Object.keys(p).forEach(function (k) { q[k] = p[k]; });
    q.a = p.b; q.b = p.a;
    var swap = [['only_a', 'only_b'], ['cpm_shared_a', 'cpm_shared_b']];
    swap.forEach(function (kv) { if (kv[0] in p) { q[kv[0]] = p[kv[1]]; q[kv[1]] = p[kv[0]]; } });
    ['tss', 'tes'].forEach(function (k) { if (p[k]) q[k] = { shared: p[k].shared, only_a: p[k].only_b, only_b: p[k].only_a, jaccard: p[k].jaccard }; });
    if (p.grid) q.grid = transpose(p.grid);
    return q;
  }
  function pairIn(block, a, b) {
    if (!block || !block.pairs) return null;
    var p = block.pairs.find(function (q) { return (q.a === a && q.b === b) || (q.a === b && q.b === a); });
    return p ? orient(p, a) : null;
  }
  function conc() { return XS().concordance; }
  function concPair(r) { return r === baseId() ? null : pairIn(conc(), baseId(), r); }
  function quantPair(r) { return r === baseId() ? null : pairIn(XS().quant, baseId(), r); }
  function endPair(r) { return r === baseId() ? null : pairIn(XS().end_agreement, baseId(), r); }
  function tuscoPair(r) {
    var t = XS().tusco_incremental;
    if (!t || r === baseId()) return null;
    return t.pairs.find(function (p) { return p.base === baseId() && p.run === r; }) || null;
  }
  // Genes gained / lost by run r against the baseline, from concordance.genes_detected.
  function genesVs(r) {
    var c = conc(), g = c && c.genes_detected, b = baseId();
    if (!g || !g.runs[r] || !g.runs[b] || r === b) return null;
    var inR = {}, inB = {};
    g.runs[r].forEach(function (i) { inR[i] = 1; });
    g.runs[b].forEach(function (i) { inB[i] = 1; });
    var gained = g.runs[r].filter(function (i) { return !inB[i]; }).map(function (i) { return g.names[i]; });
    var lost = g.runs[b].filter(function (i) { return !inR[i]; }).map(function (i) { return g.names[i]; });
    return { gained: gained, lost: lost };
  }
  function names(list) { return list.length ? ': ' + list.slice(0, 12).join(', ') + (list.length > 12 ? ' and ' + (list.length - 12) + ' more' : '') : ''; }

  // ------------------------------------------------------------------ section accessors
  function sec(r, k) { var b = B(r); return b && b.sections ? b.sections[k] : undefined; }
  function endsHalf(r, k) { var e = sec(r, 'ends'); return e && e[k] && isNum(e[k].n) ? e[k] : null; }
  function windowsOf(r) { var e = sec(r, 'ends'); return e && e.windows ? e.windows : null; }
  function within(r, k, w) {
    var h = endsHalf(r, k), ws = windowsOf(r);
    if (!h || !ws) return null;
    var i = ws.indexOf(w);
    return i >= 0 && h.within && isNum(h.within[i]) ? pc(h.within[i], h.n) : null;
  }
  // A TUSCO section with no gene inside the evaluated regions has nothing to show.
  function tusco(r) { var t = sec(r, 'tusco'); return t && isNum(t.TP) && t.genes_in_regions !== 0 ? t : null; }
  function setSampleNeeds() {
    var empty = sampleRuns().some(function (r) { var t = sec(r, 'tusco'); return t && t.genes_in_regions === 0; });
    NEEDS.tusco = empty ? 'no TUSCO gene lies inside the evaluated regions (--regions)' : 'needs --tusco';
    NEEDS.tuscoInc = empty ? NEEDS.tusco : 'needs --tusco for the baseline and the run';
  }
  function prodCounts(r) { var p = sec(r, 'productivity'); return p && p.class_counts ? p.class_counts : null; }
  function prodShare(r, k) { var cc = prodCounts(r); if (!cc) return null; var t = sum(Object.keys(cc).map(function (x) { return cc[x]; })); return pc(cc[k] || 0, t); }
  function variation(r) { var v = sec(r, 'variation'); return v && isNum(v.genes_multi) ? v : null; }
  function priming(r) { var p = sec(r, 'priming'); return p && isNum(p.n) ? p : null; }
  function quant(r) { var q = sec(r, 'quant'); return q && isNum(q.reads_assigned) ? q : null; }
  function compute(r, k) { var c = (B(r) || {}).compute; return c && isNum(c[k]) ? c[k] : null; }

  var NEEDS = {
    ends: 'needs --ref-tss or --ref-pas (sections.ends)',
    tss: 'needs --ref-tss (refTSS)',
    tes: 'needs --ref-pas (PolyA_DB)',
    motif: 'needs --ref-tss or --ref-pas with the genome (sections.ends.polya_motif)',
    tusco: 'needs --tusco',
    priming: 'needs a round-3 bundle with the genome (sections.priming)',
    variation: 'needs a round-3 bundle (sections.variation)',
    quant: 'needs quantification (counts.tsv)',
    productivity: 'needs productivity (productivity.tsv)',
    conc: 'needs isoforms.bed beside each bundle and two or more runs',
    cpm: 'needs counts.tsv beside the baseline and the run',
    endsAgree: 'needs isoforms.bed beside each bundle and two or more runs',
    tuscoInc: 'needs --tusco for the baseline and the run',
    lostDec: 'needs isoforms.bed, and counts.tsv beside the baseline',
    compute: 'needs the Nextflow trace (compute block)',
    hist: 'the bundle has no such histogram',
    others: 'needs a second run on this sample'
  };

  // ------------------------------------------------------------------ bins
  function kfmt(v) { return Math.abs(v) >= 1000 ? (v / 1000) + 'k' : String(v); }
  function sgn(v) { return v < 0 ? '−' + kfmt(-v) : v > 0 ? '+' + kfmt(v) : '0'; }
  // Bin names for fixed-edge histograms (bundle section 3): the last bin is open-ended when
  // len(counts) == len(edges).
  function binNames(edges, n, discrete) {
    var out = [];
    for (var i = 0; i < n; i++) {
      var a = edges[i], b = edges[i + 1];
      if (b === undefined) out.push('≥' + kfmt(a));
      else if (discrete) out.push(b - a === 1 ? String(a) : a + '–' + (b - 1));
      else out.push(kfmt(a) + '–' + kfmt(b));
    }
    return out;
  }
  // Signed-distance bins (A29): open below the first edge and above the last, [e_i, e_i+1) between.
  function distBins(E) {
    var full = [], ticks = [];
    for (var i = 0; i <= E.length; i++) {
      if (i === 0) { full.push('below ' + sgn(E[0]) + ' bp'); ticks.push('<' + sgn(E[0])); }
      else if (i === E.length) { full.push(sgn(E[i - 1]) + ' bp or more'); ticks.push('≥' + sgn(E[i - 1])); }
      else { full.push(sgn(E[i - 1]) + ' to ' + sgn(E[i]) + ' bp'); ticks.push(sgn(E[i - 1])); }
    }
    return { full: full, ticks: ticks, zero: E.indexOf(0) >= 0 ? E.indexOf(0) + 1 : null };
  }
  function histOf(r, key) { var h = get(B(r), 'histograms.' + key); return h && h.counts && h.counts.length ? h : null; }

  // ------------------------------------------------------------------ metric registry (A34)
  // Each family registers its own rows. lowerIsWorse: true = higher is better, false = lower is
  // better, null = no direction. vsBase rows compare a run with the baseline (none for the baseline).
  var METRICS = [];
  var FAMS = [];
  function family(fam, tab, rows) { FAMS.push({ fam: fam, tab: tab }); rows.forEach(function (m) { m.fam = fam; m.tab = tab; METRICS.push(m); }); }
  family('Structure', 'structure', [
    { name: 'Isoforms', top: true, k: 'n', lowerIsWorse: null, get: function (r) { return B(r).summary.isoforms; } },
    { name: 'Genes (FLAIR gene ids)', k: 'n', lowerIsWorse: null, get: function (r) { return B(r).summary.genes; } },
    { name: 'Mono-exon share', k: 'p', lowerIsWorse: null, get: function (r) { var x = B(r).summary; return pc(x.mono_exon, x.isoforms); } },
    { name: 'FSM share', top: true, k: 'p', lowerIsWorse: null, get: function (r) { var x = B(r); return pc(x.categories.FSM, x.summary.isoforms); } },
    { name: 'NIC + NNC share', top: true, k: 'p', lowerIsWorse: null, get: function (r) { var x = B(r); return pc((x.categories.NIC || 0) + (x.categories.NNC || 0), x.summary.isoforms); } },
    { name: 'Known junctions', top: true, k: 'p', lowerIsWorse: true, get: function (r) { var j = B(r).junctions; return j ? pc(j.known, j.total) : null; } }
  ]);
  family('Recovery', 'recovery', [
    { name: 'Reference genes detected', top: true, k: 'n', lowerIsWorse: true, get: function (r) { return B(r).summary.ref_genes_detected; } },
    { name: 'Reference transcripts with an FSM', k: 'n', lowerIsWorse: true, get: function (r) { return B(r).summary.ref_transcripts_fsm; } },
    { name: 'Reference transcripts with a CIM', top: true, k: 'n', lowerIsWorse: true, get: function (r) { return B(r).summary.ref_transcripts_cim; } },
    { name: 'TUSCO sensitivity', top: true, k: 'p', lowerIsWorse: true, needs: 'tusco', get: function (r) { var t = tusco(r); return t ? asPct(t.sensitivity) : null; } },
    { name: 'TUSCO precision', top: true, k: 'p', lowerIsWorse: true, needs: 'tusco', get: function (r) { var t = tusco(r); return t ? asPct(t.precision) : null; } },
    { name: 'TUSCO false negatives', k: 'n', lowerIsWorse: false, needs: 'tusco', get: function (r) { var t = tusco(r); return t ? t.FN : null; } },
    { name: 'Genes lost vs baseline', k: 'n', lowerIsWorse: false, vsBase: true, needs: 'conc', get: function (r) { var g = genesVs(r); return g ? g.lost.length : null; } },
    { name: 'Genes gained vs baseline', k: 'n', lowerIsWorse: null, vsBase: true, needs: 'conc', get: function (r) { var g = genesVs(r); return g ? g.gained.length : null; } }
  ]);
  family('Diversity', 'diversity', [
    { name: 'Isoforms per gene (mean)', k: 'f2', lowerIsWorse: null, get: function (r) { var x = B(r).summary; return x.genes ? x.isoforms / x.genes : null; } },
    { name: 'Single-isoform genes', k: 'p', lowerIsWorse: null, needs: 'hist', get: function (r) { var h = histOf(r, 'isoforms_per_gene'); return h ? pc(h.counts[0], sum(h.counts)) : null; } },
    { name: 'Productive (PRO)', k: 'p', lowerIsWorse: null, needs: 'productivity', get: function (r) { return prodShare(r, 'PRO'); } },
    { name: 'Premature stop, NMD candidate (PTC)', k: 'p', lowerIsWorse: null, needs: 'productivity', get: function (r) { return prodShare(r, 'PTC'); } },
    { name: 'Multi-isoform genes varying by ends only', k: 'p', lowerIsWorse: null, needs: 'variation', get: function (r) { var v = variation(r); return v ? pc(v.ends_only, v.genes_multi) : null; } },
    { name: 'Multi-isoform genes varying by splicing only', k: 'p', lowerIsWorse: null, needs: 'variation', get: function (r) { var v = variation(r); return v ? pc(v.splice_only, v.genes_multi) : null; } }
  ]);
  family('Agreement', 'agreement', [
    { name: 'Intron-chain Jaccard vs baseline', k: 'f2', lowerIsWorse: null, vsBase: true, needs: 'conc', get: function (r) { var p = concPair(r); return p ? p.jaccard : null; } },
    { name: 'TSS-cluster Jaccard vs baseline', k: 'f2', lowerIsWorse: null, vsBase: true, needs: 'endsAgree', get: function (r) { var p = endPair(r); return p ? p.tss.jaccard : null; } },
    { name: 'TTS-cluster Jaccard vs baseline', k: 'f2', lowerIsWorse: null, vsBase: true, needs: 'endsAgree', get: function (r) { var p = endPair(r); return p ? p.tes.jaccard : null; } },
    { name: 'CPM agreement with baseline (Pearson r)', k: 'f3', lowerIsWorse: null, vsBase: true, needs: 'cpm', get: function (r) { var p = quantPair(r); return p ? p.pearson : null; } },
    { name: 'TUSCO genes found only here', k: 'n', lowerIsWorse: null, vsBase: true, needs: 'tuscoInc', get: function (r) { var p = tuscoPair(r); return p ? p.gained : null; } }
  ]);
  family('Ends', 'ends', [
    { name: '5′ end within 50 bp of refTSS', top: true, k: 'p', lowerIsWorse: true, needs: 'tss', get: function (r) { return within(r, 'tss', 50); } },
    { name: '3′ end within 50 bp of PolyA_DB', top: true, k: 'p', lowerIsWorse: true, needs: 'tes', get: function (r) { return within(r, 'tes', 50); } },
    { name: 'PolyA motif before the 3′ end', k: 'p', lowerIsWorse: true, needs: 'motif', get: function (r) { var m = get(B(r), 'sections.ends.polya_motif'); return m ? pc(m.n_with_motif, m.n) : null; } },
    { name: 'FSM with both ends ≤ 50 bp (CIM)', k: 'p', lowerIsWorse: true, get: function (r) { var x = B(r); return pc(x.fsm_subcategories.CIM, x.categories.FSM); } }
  ]);
  family('Reads', 'reads', [
    { name: 'Reads assigned', top: true, k: 'n', lowerIsWorse: null, needs: 'quant', get: function (r) { var q = quant(r); return q ? q.reads_assigned : null; } },
    { name: 'Isoforms with reads', k: 'n', lowerIsWorse: null, needs: 'quant', get: function (r) { var q = quant(r); return q ? q.isoforms_quantified : null; } },
    { name: 'Internal-priming risk (≥ 60% A after the 3′ end)', top: true, k: 'p', lowerIsWorse: false, needs: 'priming', get: function (r) { var p = priming(r); return p ? pc(p.risk, p.n) : null; } }
  ]);
  family('Compute', 'compute', [
    { name: 'Wall time', k: 'wall', lowerIsWorse: false, needs: 'compute', get: function (r) { return compute(r, 'wall_minutes'); } },
    { name: 'CPU time', top: true, k: 'cpu', lowerIsWorse: false, needs: 'compute', get: function (r) { return compute(r, 'cpu_hours'); } },
    { name: 'Peak memory', k: 'mem', lowerIsWorse: false, needs: 'compute', get: function (r) { return compute(r, 'peak_rss_gb'); } }
  ]);
  function mval(m, r) { try { var v = m.get(r); return isNum(v) ? v : null; } catch (e) { return null; } }
  function mfmt(m, v) {
    if (!isNum(v)) return '–';
    return m.k === 'p' ? pct(v) : m.k === 'f2' ? f2(v) : m.k === 'f3' ? f3(v) : m.k === 'wall' ? fmtWall(v) : m.k === 'cpu' ? fmtCpu(v) : m.k === 'mem' ? fmtMem(v) : N(v);
  }
  function mdelta(m, v, b) {
    if (!isNum(v) || !isNum(b) || m.vsBase) return '';
    var d = v - b, sign = d > 0 ? '+' : d < 0 ? '−' : '±';
    if (m.k === 'p') return sign + Math.abs(d).toFixed(1) + ' pp';
    if (m.k === 'f2') return sign + Math.abs(d).toFixed(2);
    if (!b) return v ? 'new' : '±0';
    return sign + Math.abs(100 * d / b).toFixed(0) + '%';
  }
  function mdir(m, v, b) {
    if (!isNum(v) || !isNum(b) || m.lowerIsWorse === null || m.lowerIsWorse === undefined || v === b || m.vsBase) return '';
    return (v > b) === m.lowerIsWorse ? 'better' : 'worse';
  }
  function mscore(m, v, b, withCompute) {
    if (!isNum(v) || !isNum(b) || m.vsBase || (m.fam === 'Compute' && !withCompute)) return 0;
    if (m.k === 'p') return Math.abs(v - b) / 15;
    return 2 * Math.abs(v - b) / Math.max(Math.abs(v), Math.abs(b), 1e-9);
  }

  // ------------------------------------------------------------------ notes
  var QUANT_INPUT = { reads: 'from reads', bam: 'from the genome BAM' };
  function quantNote(parent, runs) {
    var groups = {}, order = [];
    runs.forEach(function (r) {
      var q = sec(r, 'quant');
      if (!q) return;
      var k = q.input || 'unknown';
      if (!groups[k]) { groups[k] = []; order.push(k); }
      groups[k].push(lab(r));
    });
    if (!order.length) return;
    var parts = order.map(function (k) { return (QUANT_INPUT[k] ? 'Counted ' + QUANT_INPUT[k] : 'Input not recorded') + ': ' + groups[k].join(', '); });
    var mixed = order.filter(function (k) { return QUANT_INPUT[k]; }).length > 1;
    parent.appendChild(el('p', { class: 'callout', role: 'note' }, [el('strong', { text: mixed ? 'Not like for like. ' : 'How reads are counted. ' }),
      'FLAIR dev counts reads from the genome BAM; 2.x and 3.0.1 count from the reads (A21). ' + parts.join('; ') + '.']));
  }
  var RULE = DATA.baseline_rule || {};
  function baselineWhy(id) {
    if (id !== DATA.baseline) return 'chosen in this page';
    if (RULE.rule === 'given') return 'set with --baseline';
    if (RULE.rule === 'earliest_launch') return RULE.launches > 1 ? 'first run of the earliest of ' + RULE.launches + ' launches' : 'first run of the launch';
    if (RULE.rule === 'name_order') return 'first run id in name order (these bundles record no run order)';
    return 'default';
  }

  // ------------------------------------------------------------------ views
  function viewDiff(root) {
    var b = baseId();
    section(root, 'What differs', 'Every run against the baseline. Pick another baseline, or hide runs, in the bar above.');
    var gtf = (B(b).inputs || {}).gtf;
    root.appendChild(el('div', { class: 'callout' }, [el('strong', { text: 'Read with care. ' }),
      'Runs are scored against the reference annotation they were built with' + (gtf ? ' (' + gtf + ')' : '') +
      ', so FSM and gene recovery are recall checks, not accuracy. TUSCO and the refTSS / PolyA_DB end databases are independent evidence when present; end databases give precision only.']));
    if (!others().length) note(root, sampleRuns().length > 1 ? 'Every other run is hidden: show one in the bar above to compare.' : 'Only one run on this sample, so there is nothing to compare. Its numbers are in the table below.', 'note');
    if (others().length) whereDiffer(root);
    var g = grid(root);
    others().forEach(function (r) {
      var body = card(g, { title: lab(r) + ' vs ' + lab(b) });
      var p = concPair(r), q = quantPair(r), t = tusco(r), tb = tusco(b);
      var list = [];
      if (p) list.push(['Intron chains shared', N(p.shared), 'Jaccard ' + f2(p.jaccard)]);
      list.push(['Isoforms', N(B(r).summary.isoforms), mdelta(METRICS[0], B(r).summary.isoforms, B(b).summary.isoforms) + ' vs baseline']);
      if (t && tb) list.push(['TUSCO sensitivity', pct(asPct(t.sensitivity)), mdelta({ k: 'p' }, asPct(t.sensitivity), asPct(tb.sensitivity)) + ' vs baseline']);
      else list.push(['Reference genes detected', N(B(r).summary.ref_genes_detected), mdelta({ k: 'n' }, B(r).summary.ref_genes_detected, B(b).summary.ref_genes_detected) + ' vs baseline']);
      [['tss', '5′ within 50 bp of refTSS'], ['tes', '3′ within 50 bp of PolyA_DB']].forEach(function (k) {
        var v = within(r, k[0], 50), w = within(b, k[0], 50);
        if (isNum(v)) list.push([k[1], pct(v), isNum(w) ? mdelta({ k: 'p' }, v, w) + ' vs baseline' : '']);
      });
      if (q && isNum(q.pearson)) list.push(['CPM agreement', 'r = ' + f3(q.pearson), 'log CPM, ' + N(q.n) + ' shared structures']);
      if (list.length % 2) list = list.slice(0, list.length - (list.length > 4 ? 1 : 0));
      tiles(body, list);
      var ranked = METRICS.map(function (m) { var v = mval(m, r), w = mval(m, b); return { m: m, v: v, w: w, sc: mscore(m, v, w) }; })
        .filter(function (x) { return x.sc > 0.05; }).sort(function (x, y) { return y.sc - x.sc; }).slice(0, 5);
      body.appendChild(el('p', { class: 'eyebrow', style: 'margin:16px 0 6px', text: 'Biggest differences' }));
      if (!ranked.length) note(body, 'No metric differs by more than a few percent.');
      var ul = el('ul', { class: 'diffs' });
      ranked.forEach(function (x) {
        var dir = mdir(x.m, x.v, x.w);
        ul.appendChild(el('li', {}, [el('button', { type: 'button', class: 'linkish', onclick: function () { go(x.m.tab); }, text: x.m.name }),
          el('span', { class: 'num' }, [' ' + mfmt(x.m, x.v) + ' vs ' + mfmt(x.m, x.w) + ' ']), el('span', { class: 'delta num', text: mdelta(x.m, x.v, x.w) }),
          dir ? el('span', { class: 'dir', text: dir }) : null]));
      });
      body.appendChild(ul);
    });
    metricTable(root);
  }
  // The headline metrics (top: true) as dot rows. Each row runs from its own lowest to highest
  // value (min-max, padded 8%), so small gaps are visible; the range at the right gives the scale.
  function whereDiffer(root) {
    var runs = sel(), b = baseId();
    var rows = METRICS.filter(function (m) { return m.top && !m.vsBase; }).map(function (m) {
      return { m: m, vals: runs.map(function (r) { return { r: r, v: mval(m, r) }; }).filter(function (d) { return d.v !== null; }) };
    }).filter(function (x) { return x.vals.length; });
    if (!rows.length) return;
    var body = card(root, { title: 'Where the runs differ', wide: true,
      dek: 'One row per metric. Each row runs from its lowest to its highest value, so read the range at the right for how big a gap is. The ringed dot is the baseline.' });
    var viz = el('div', { class: 'viz' });
    body.appendChild(viz);
    var W = width(viz), narrow = W < 560;
    var lw = narrow ? 0 : Math.min(250, Math.floor(W * 0.36)), rw = narrow ? 0 : 128, pad = 10, lane = 13;
    var tx0 = lw + pad, tx1 = W - rw - pad;
    var layout = rows.map(function (x) {
      var vs = x.vals.map(function (d) { return d.v; });
      var lo = Math.min.apply(null, vs), hi = Math.max.apply(null, vs), span = hi - lo;
      var dlo = lo - span * 0.08, dhi = hi + span * 0.08;
      var X = function (v) { return span > 0 ? tx0 + (tx1 - tx0) * (v - dlo) / (dhi - dlo) : (tx0 + tx1) / 2; };
      var sorted = x.vals.slice().sort(function (p, q) { return p.v - q.v; }), ends = [];
      sorted.forEach(function (d) { var px = X(d.v), l = 0; while (ends[l] !== undefined && px - ends[l] < lane) l++; ends[l] = px; d.px = px; d.lane = l; });
      var nl = Math.max(1, ends.length);
      return { m: x.m, sorted: sorted, lo: lo, hi: hi, nl: nl, h: (narrow ? 18 : 0) + Math.max(30, 16 + nl * lane) };
    });
    var H = layout.reduce(function (t, x) { return t + x.h; }, 0);
    var svg = svgBox(viz, W, H, 'Where the runs differ: one dot per run on each metric’s own scale');
    var y = 0;
    layout.forEach(function (x, i) {
      var m = x.m, top = y + (narrow ? 18 : 0), mid = top + (x.h - (narrow ? 18 : 0)) / 2;
      var range = x.lo === x.hi ? mfmt(m, x.lo) : mfmt(m, x.lo) + ' – ' + mfmt(m, x.hi);
      var all = x.sorted.slice().reverse().map(function (d) { return lab(d.r) + (d.r === b ? ' (baseline)' : '') + ' ' + mfmt(m, d.v); }).join(' · ');
      if (i) s('line', { x1: 0, x2: W, y1: y, y2: y, class: 'gl rowrule' }, svg);
      s('rect', { x: 0, y: y, width: W, height: x.h, class: 'hit', 'data-tip': m.name + ': ' + all }, svg);
      if (narrow) {
        s('text', { x: 0, y: y + 12, class: 'lab', text: trunc(m.name, Math.floor(W / 7.5)) + ' (' + range + ')' }, svg);
      } else {
        s('text', { x: lw - 4, y: mid, 'text-anchor': 'end', 'dominant-baseline': 'central', class: 'lab', text: trunc(m.name, Math.floor(lw / 6.8)) }, svg);
        s('text', { x: W, y: mid, 'text-anchor': 'end', 'dominant-baseline': 'central', class: 'val', text: range }, svg);
      }
      s('line', { x1: tx0, x2: tx1, y1: mid, y2: mid, class: 'gl' }, svg);
      var cy0 = mid - (x.nl - 1) * lane / 2;
      x.sorted.forEach(function (d) {
        var isBase = d.r === b;
        s('circle', { cx: d.px, cy: cy0 + d.lane * lane, r: isBase ? 6.5 : 5.5, class: isBase ? 'base-dot' : 'dot', style: fill(col(d.r)),
          'data-tip': lab(d.r) + (isBase ? ' (baseline)' : '') + ' · ' + m.name + ': ' + mfmt(m, d.v) }, svg);
      });
      y += x.h;
    });
    addTable(body, ['Metric'].concat(runs.map(lab)), layout.map(function (x) {
      return [x.m.name].concat(runs.map(function (r) { return mfmt(x.m, mval(x.m, r)); }));
    }));
  }
  function metricTable(root) {
    var runs = sel(), b = baseId();
    var body = card(root, { title: 'Every metric, one row each', wide: true,
      dek: 'Dots place each run on the row’s own scale (zero to the largest value); the baseline is the larger, outlined dot. Click a family to open its plots.' });
    var head = el('tr', {}, [el('th', { text: 'Metric' })].concat(runs.map(function (r) { return el('th', { class: 'num runhead' }, [runTag(r, 16)]); }))
      .concat([el('th', { text: 'Spread' }), el('th', { class: 'num', text: 'Largest gap vs baseline' })]));
    var tbody = el('tbody');
    FAMS.forEach(function (F) {
      var rows = METRICS.filter(function (m) { return m.fam === F.fam; });
      var shown = rows.filter(function (m) { return runs.some(function (r) { return mval(m, r) !== null; }); });
      var hidden = rows.filter(function (m) { return shown.indexOf(m) < 0; });
      tbody.appendChild(el('tr', { class: 'famrow' }, [el('td', { colspan: runs.length + 3 }, [el('button', { type: 'button', class: 'linkish', onclick: function () { go(F.tab); }, text: F.fam + ' →' })])]));
      shown.forEach(function (m) {
        var vals = runs.map(function (r) { return mval(m, r); });
        var tr = el('tr', {}, [el('td', { text: m.name })].concat(vals.map(function (v, i) {
          var txt = m.vsBase && runs[i] === b ? 'baseline' : mfmt(m, v);
          var d = i ? mdelta(m, v, vals[0]) : '';
          return el('td', { class: 'num' + (m.vsBase && runs[i] === b ? ' muted' : '') }, [txt, d ? el('span', { class: 'cell-delta', text: d }) : null]);
        })));
        var strip = el('td', { class: 'strip' });
        var sv = s('svg', { width: 130, height: 16, viewBox: '0 0 130 16', class: 'chart' }, null);
        var mx = maxOf(vals);
        if (!(mx > 0)) mx = 1;
        s('line', { x1: 4, x2: 126, y1: 8, y2: 8, class: 'gl' }, sv);
        vals.forEach(function (v, i) {
          if (!isNum(v)) return;
          s('circle', { cx: 4 + 122 * Math.max(0, v) / mx, cy: 8, r: i === 0 ? 5 : 3.5, class: i === 0 ? 'base-dot' : 'dot', style: fill(col(runs[i])), 'data-tip': lab(runs[i]) + ': ' + mfmt(m, v) }, sv);
        });
        strip.appendChild(sv);
        tr.appendChild(strip);
        var best = '', bs = 0;
        if (m.vsBase) {
          var lo = null, who = null;
          vals.slice(1).forEach(function (v, i) { if (isNum(v) && (lo === null || v < lo)) { lo = v; who = runs[i + 1]; } });
          best = lo === null ? '' : 'lowest ' + mfmt(m, lo) + (runs.length > 2 ? ' (' + trunc(lab(who), 16) + ')' : '');
        } else {
          vals.slice(1).forEach(function (v, i) {
            var sc = mscore(m, v, vals[0], true);
            if (sc > bs) { bs = sc; best = mdelta(m, v, vals[0]) + (mdir(m, v, vals[0]) ? ' ' + mdir(m, v, vals[0]) : '') + (runs.length > 2 ? ' (' + trunc(lab(runs[i + 1]), 16) + ')' : ''); }
          });
        }
        tr.appendChild(el('td', { class: 'num delta', text: best }));
        tbody.appendChild(tr);
      });
      if (hidden.length) {
        var why = {};
        hidden.forEach(function (m) { var k = m.vsBase && runs.length < 2 ? NEEDS.others : (NEEDS[m.needs] || m.needs || 'no data'); (why[k] = why[k] || []).push(m.name); });
        tbody.appendChild(el('tr', { class: 'hiddenrow' }, [el('td', { colspan: runs.length + 3, text: 'Hidden: ' + Object.keys(why).map(function (k) { return why[k].join(', ') + ' (' + k + ')'; }).join('; ') + '.' })]));
      }
    });
    body.appendChild(el('div', { class: 'tablewrap' }, [el('table', { class: 'data metrics' }, [el('thead', {}, [head]), tbody])]));
  }

  var CATS = ['FSM', 'ISM', 'NIC', 'NNC', 'genic', 'genic_intron', 'antisense', 'intergenic', 'fusion'];
  var CAT_LABEL = { FSM: 'FSM', ISM: 'ISM', NIC: 'NIC', NNC: 'NNC', genic: 'Genic', genic_intron: 'Genic intron', antisense: 'Antisense', intergenic: 'Intergenic', fusion: 'Fusion' };
  function byRun(runs, f) { var o = {}; runs.forEach(function (r) { o[r] = f(r); }); return o; }

  function viewStructure(root) {
    section(root, 'Structure', 'What each transcriptome is made of. Bars are shares of each run’s own isoforms, so runs of different sizes compare directly.');
    var runs = sel(), g = grid(root);
    var cats = CATS.filter(function (c) { return runs.some(function (r) { return B(r).categories[c] > 0; }); });
    var b1 = card(g, { title: 'Isoforms by structural category', dek: 'SQANTI-style classes against the reference annotation.' + (cats.length < CATS.length ? ' Classes no run has are left out.' : '') });
    hbars(b1, { cats: cats.map(function (c) { return CAT_LABEL[c]; }), vals: byRun(runs, function (r) { var x = B(r); return cats.map(function (c) { return pc(x.categories[c] || 0, x.summary.isoforms); }); }),
      fmt: function (v) { return pct(v); }, labelW: 96,
      tip: function (r, i) { var x = B(r); return N(x.categories[cats[i]] || 0) + ' of ' + N(x.summary.isoforms) + ' (' + pct(pc(x.categories[cats[i]] || 0, x.summary.isoforms)) + ')'; } });
    var b2 = card(g, { title: 'How many isoforms', dek: 'Total, multi-exon and mono-exon isoforms.' });
    hbars(b2, { cats: ['Isoforms', 'Multi-exon', 'Mono-exon'], vals: byRun(runs, function (r) { var x = B(r).summary; return [x.isoforms, x.multi_exon, x.mono_exon]; }), fmt: N, labelW: 100 });
    [['transcript_length', 'Transcript length', 'Transcript length (nt)', false], ['exons_per_isoform', 'Exons per isoform', 'Exons', true]].forEach(function (h) {
      panel(g, { title: h[1] }, function (r) { return histOf(r, h[0]); }, NEEDS.hist, function (body, ok) {
        var e = histOf(ok[0], h[0]);
        var nm = binNames(e.edges, e.counts.length, h[3]);
        vhist(body, { runs: ok, bins: nm, ticks: h[3] ? nm : e.edges.slice(0, e.counts.length).map(function (x, i) { return i === e.edges.length - 1 ? '≥' + kfmt(x) : kfmt(x); }), counts: byRun(ok, function (r) { return histOf(r, h[0]).counts; }), xlab: h[2] + (h[3] ? '' : ', bin start'), xname: h[2] });
      });
    });
    panel(g, { title: 'End precision by category', wide: true, dek: 'Share of each category’s isoforms whose 5′ end is within 50 bp of a refTSS site and whose 3′ end is within 50 bp of a PolyA_DB site. Shows whether novel classes carry worse ends.' },
      function (r) { return get(B(r), 'sections.ends.by_category'); }, NEEDS.ends, function (body, ok) {
        var present = CATS.filter(function (c) { return ok.some(function (r) { var bc = sec(r, 'ends').by_category[c]; return bc && bc.n > 0; }); });
        var two = el('div', { class: 'grid sub' });
        body.appendChild(two);
        [['tss', 'tss_within_50', '5′ end vs refTSS'], ['tes', 'tes_within_50', '3′ end vs PolyA_DB']].forEach(function (k) {
          var rs = ok.filter(function (r) { return endsHalf(r, k[0]); });
          var h = el('div', {}, [el('p', { class: 'eyebrow', style: 'margin:0 0 6px', text: k[2] })]);
          two.appendChild(h);
          if (!rs.length) { note(h, 'Hidden: ' + NEEDS[k[0]] + '.'); return; }
          hbars(h, { runs: rs, cats: present.map(function (c) { return CAT_LABEL[c]; }), max: 100, labelW: 96, fmt: function (v) { return pct(v, 0); },
            vals: byRun(rs, function (r) { var bc = sec(r, 'ends').by_category; return present.map(function (c) { var x = bc[c]; return x && x.n ? pc(x[k[1]], x.n) : null; }); }),
            tip: function (r, i) { var x = sec(r, 'ends').by_category[present[i]]; return N(x[k[1]]) + ' of ' + N(x.n) + ' (' + pct(pc(x[k[1]], x.n)) + ')'; } });
        });
      });
  }

  function viewRecovery(root) {
    var runs = sel(), b = baseId();
    section(root, 'Recovery', 'First against TUSCO, a truth set the transcriptomes were not built from, then against the reference annotation.');
    root.appendChild(el('p', { class: 'eyebrow', style: 'margin:0 0 10px', text: 'Truth set' }));
    panel(root, { title: 'TUSCO', wide: true, dek: 'Genes expected to have exactly one isoform (SQANTI3 TUSCO rules). TP: the isoform is reported; PTP: a partial match; FP: only other isoforms; FN: the gene is missed.' },
      tusco, NEEDS.tusco, function (body, ok) {
        tiles(body, ok.map(function (r) {
          var t = tusco(r);
          return [runTag(r, true), pct(asPct(t.sensitivity)), 'sensitivity · precision ' + pct(asPct(t.precision))];
        }));
        var t0 = tusco(ok[0]);
        if (isNum(t0.genes_in_regions) && isNum(t0.genes_total)) note(body, N(t0.genes_in_regions) + ' of ' + N(t0.genes_total) + ' TUSCO genes are inside the evaluated regions.');
        var two = el('div', { class: 'grid sub', style: 'margin-top:16px' });
        body.appendChild(two);
        var h1 = el('div', {}, [el('p', { class: 'eyebrow', style: 'margin:0 0 6px', text: 'Sensitivity and precision' })]);
        var h2 = el('div', {}, [el('p', { class: 'eyebrow', style: 'margin:0 0 6px', text: 'TUSCO genes by call' })]);
        two.appendChild(h1); two.appendChild(h2);
        hbars(h1, { runs: ok, cats: ['Sensitivity', 'Precision'], vals: byRun(ok, function (r) { var t = tusco(r); return [asPct(t.sensitivity), asPct(t.precision)]; }), fmt: function (v) { return pct(v); }, max: 100, labelW: 84 });
        hbars(h2, { runs: ok, cats: ['True positive', 'Partial', 'False positive', 'False negative'], vals: byRun(ok, function (r) { var t = tusco(r); return [t.TP, t.PTP, t.FP, t.FN]; }), fmt: N, labelW: 104 });
      });
    root.appendChild(el('p', { class: 'eyebrow', style: 'margin:18px 0 10px', text: 'Against the reference annotation' }));
    var g = grid(root);
    var b3 = card(g, { title: 'Genes detected', dek: 'Reference genes with at least one FSM, ISM, NIC, NNC or fusion isoform.' });
    hbars(b3, { cats: ['Genes'], vals: byRun(runs, function (r) { return [B(r).summary.ref_genes_detected]; }), fmt: N, labelW: 60,
      tip: function (r) { var x = B(r).summary; return N(x.ref_genes_detected) + (isNum(x.ref_genes_total) ? ' of ' + N(x.ref_genes_total) + ' reference genes' : ''); } });
    if (!others().length) off(g, 'Genes gained and lost vs baseline', NEEDS.others);
    else panel(g, { title: 'Genes gained and lost vs baseline', dek: 'Reference genes one run detects and the baseline does not, and the reverse.', runs: others() }, genesVs, NEEDS.conc, function (body, ok) {
      pairBars(body, ok.map(function (r) {
        var gv = genesVs(r);
        return { label: trunc(lab(r), 24), bars: [
          { run: b, v: gv.lost.length, tip: 'Only in ' + lab(b) + ', not ' + lab(r) + ': ' + N(gv.lost.length) + ' genes' + names(gv.lost) },
          { run: r, v: gv.gained.length, tip: 'Only in ' + lab(r) + ': ' + N(gv.gained.length) + ' genes' + names(gv.gained) }] };
      }), N, 170);
      note(body, 'Per run: the first bar (baseline colour) counts genes only the baseline detects, the second (run colour) genes only the run detects.');
      addTable(body, ['Run', 'Only in baseline', 'Only in run'], ok.map(function (r) { var gv = genesVs(r); return [lab(r), N(gv.lost.length), N(gv.gained.length)]; }));
    });
    var b5 = card(g, { title: 'FSM isoforms with exact ends', dek: 'FSM isoforms, and those with both ends within 50 bp of the matched reference transcript (CIM). A recall check: the runs were built with this annotation.' });
    hbars(b5, { cats: ['FSM', 'CIM: both ends ≤ 50 bp'], vals: byRun(runs, function (r) { var x = B(r); return [x.categories.FSM, x.fsm_subcategories.CIM]; }), fmt: N, labelW: 150,
      tip: function (r, i) { var x = B(r); return i ? N(x.fsm_subcategories.CIM) + ' of ' + N(x.categories.FSM) + ' FSM (' + pct(pc(x.fsm_subcategories.CIM, x.categories.FSM)) + ')' : N(x.categories.FSM) + ' isoforms'; } });
    var ld = XS().lost_deciles, lb = ld && ld.by_base ? ld.by_base[b] : null;
    if (!others().length) off(g, 'Who are the lost genes?', NEEDS.others);
    else if (!lb) off(g, 'Who are the lost genes?', NEEDS.lostDec);
    else panel(g, { title: 'Who are the lost genes?', dek: 'Baseline genes a run misses, by the baseline’s expression decile (gene CPM, low to high). Losses piled in the low deciles are coverage dropouts; losses in high deciles are real.', runs: others() },
      function (r) { return lb.lost[r]; }, NEEDS.lostDec, function (body, ok) {
        var bins = [], ticks = [];
        for (var i = 0; i < lb.n.length; i++) { ticks.push(String(i + 1)); bins.push('Decile ' + (i + 1) + ' (' + N(lb.n[i]) + ' genes, CPM ≤ ' + (isNum(lb.cpm_max[i]) ? lb.cpm_max[i].toFixed(1) : '–') + ')'); }
        vbars(body, { runs: ok, cap: 100, bins: bins, ticks: ticks, xname: 'Baseline decile', xlab: 'Baseline expression decile (low → high)',
          vals: byRun(ok, function (r) { return lb.lost[r].map(function (n, i) { return lb.n[i] ? 100 * n / lb.n[i] : null; }); }),
          fmt: function (v) { return pct(v, 0); }, yfmt: tickPct,
          tip: function (r, i) { return N(lb.lost[r][i]) + ' of ' + N(lb.n[i]) + ' genes lost (' + pct(pc(lb.lost[r][i], lb.n[i])) + ')'; } });
        note(body, 'Share of each decile’s ' + N(lb.genes) + ' baseline genes that the run does not detect. CPM from the baseline’s counts.tsv.');
      });
  }

  var PROD_ORDER = ['PRO', 'PTC', 'NST', 'NGO'];
  var PROD_LABEL = { PRO: 'Productive (PRO)', PTC: 'Premature stop · NMD (PTC)', NST: 'No stop codon (NST)', NGO: 'No annotated start (NGO)', NA: 'Not classified (NA)' };
  function viewDiversity(root) {
    section(root, 'Diversity', 'How many isoforms each gene gets, what they code for, and what separates them.');
    var g = grid(root);
    panel(g, { title: 'Isoforms per gene', dek: 'FLAIR genes by how many isoforms they carry.' }, function (r) { return histOf(r, 'isoforms_per_gene'); }, NEEDS.hist, function (body, ok) {
      var e = histOf(ok[0], 'isoforms_per_gene');
      vhist(body, { runs: ok, bins: binNames(e.edges, e.counts.length, true), counts: byRun(ok, function (r) { return histOf(r, 'isoforms_per_gene').counts; }), xlab: 'Isoforms in the gene', xname: 'Isoforms' });
    });
    panel(g, { title: 'Coding and NMD', dek: 'FLAIR productivity classes, as shares of classified isoforms. PTC isoforms are NMD candidates.' }, prodCounts, NEEDS.productivity, function (body, ok) {
      var ks = PROD_ORDER.slice();
      ok.forEach(function (r) { Object.keys(prodCounts(r)).forEach(function (k) { if (ks.indexOf(k) < 0) ks.push(k); }); });
      ks = ks.filter(function (k) { return ok.some(function (r) { return prodCounts(r)[k] > 0; }); });
      hbars(body, { runs: ok, cats: ks.map(function (k) { return PROD_LABEL[k] || k; }), vals: byRun(ok, function (r) { return ks.map(function (k) { return prodShare(r, k); }); }),
        fmt: function (v) { return pct(v); }, labelW: 170,
        tip: function (r, i) { var cc = prodCounts(r); return N(cc[ks[i]] || 0) + ' isoforms (' + pct(prodShare(r, ks[i])) + ')'; } });
    });
    panel(g, { title: 'How genes vary', wide: true, dek: 'Each gene with two or more isoforms, classed by what separates its isoforms: ends only (one intron chain, ends more than 50 bp apart), splicing only (several chains, ends within 50 bp for a shared chain), or both.' },
      variation, NEEDS.variation, function (body, ok) {
        stack100(body, { runs: ok, parts: [{ label: 'Ends only', seq: 2 }, { label: 'Splicing only', seq: 4 }, { label: 'Both', seq: 6 }],
          vals: function (r) { var v = variation(r); return [v.ends_only, v.splice_only, v.both]; } });
        note(body, ok.map(function (r) { return lab(r) + ': ' + N(variation(r).genes_multi) + ' multi-isoform genes'; }).join(' · ') + '.');
      });
  }

  function viewAgreement(root) {
    var runs = sel(), b = baseId(), c = conc();
    section(root, 'Agreement', 'Isoforms matched by intron chain; mono-exon isoforms by at least 50% reciprocal overlap on the same strand.');
    var g = grid(root);
    var inConc = c ? runs.filter(function (r) { return c.runs.indexOf(r) >= 0; }) : [];
    if (inConc.length < 2) off(g, 'Intron-chain Jaccard', runs.length < 2 ? NEEDS.others : NEEDS.conc);
    else {
      var b1 = card(g, { title: 'Intron-chain Jaccard', dek: 'Shared structures over all structures, for every pair of runs.' });
      var viz = el('div', { class: 'viz' });
      b1.appendChild(viz);
      var W = width(viz), n = inConc.length, lw = Math.min(120, W * 0.3), cs = Math.min(92, (W - lw - 4) / n), H = 24 + cs * n;
      var svg = svgBox(viz, W, H, 'Jaccard matrix');
      var chars = Math.max(3, Math.floor(cs / 7.4));
      inConc.forEach(function (r, i) {
        s('text', { x: lw + cs * i + cs / 2, y: 12, 'text-anchor': 'middle', class: 'lab', text: trunc(lab(r), chars) }, svg);
        s('rect', { x: lw - 14, y: 24 + cs * i + cs / 2 - 5, width: 10, height: 10, rx: 2, style: fill(col(r)) }, svg);
        s('text', { x: lw - 20, y: 24 + cs * i + cs / 2, 'text-anchor': 'end', 'dominant-baseline': 'central', class: 'lab', text: trunc(lab(r), Math.floor((lw - 20) / 6.5)) }, svg);
        inConc.forEach(function (q, j) {
          var v = i === j ? 1 : (pairIn(c, r, q) || {}).jaccard;
          if (!isNum(v)) return;
          var k = seqIdx(v);
          s('rect', { x: lw + cs * j + 1, y: 24 + cs * i + 1, width: cs - 2, height: cs - 2, rx: 3, style: fill('var(--seq-' + k + ')'), 'data-tip': lab(r) + ' × ' + lab(q) + ': ' + f3(v) }, svg);
          if (cs > 26) s('text', { x: lw + cs * j + cs / 2, y: 24 + cs * i + cs / 2, 'text-anchor': 'middle', 'dominant-baseline': 'central', class: 'inbar', style: fill('var(--seq-on-' + k + ')'), text: f2(v) }, svg);
        });
      });
      addTable(b1, [''].concat(inConc.map(lab)), inConc.map(function (r) { return [lab(r)].concat(inConc.map(function (q) { return r === q ? '1.00' : f2((pairIn(c, r, q) || {}).jaccard); })); }));
    }
    if (!others().length) off(g, 'Shared and unique intron chains', NEEDS.others);
    else panel(g, { title: 'Shared and unique vs baseline', dek: 'Structures only in the baseline, shared, and only in the run.', runs: others() }, concPair, NEEDS.conc, function (body, ok) {
      sharedBars(body, ok.map(function (r) { return { run: r, p: concPair(r) }; }), 'structures');
      if (c && c.runs.length > 2) note(body, 'All ' + c.runs.length + ' runs share ' + N(c.shared_by_all) + ' of ' + N(c.union) + ' structures (' + pct(pc(c.shared_by_all, c.union)) + ').');
    });
    if (!others().length) off(root, 'Expression on shared structures', NEEDS.others);
    else panel(root, { title: 'Expression on shared structures', wide: true, runs: others(), dek: 'Density of log10(CPM + 1) for structures both runs report, baseline on x. Darker = more structures; dashed = equal CPM.' },
      quantPair, NEEDS.cpm, function (body, ok) {
        quantNote(body, [b].concat(ok));
        var lm = XS().quant.log_max || 5;
        var row = el('div', { class: 'multiples dens viz' });
        body.appendChild(row);
        ok.forEach(function (r) {
          var p = quantPair(r);
          var m = el('div', { class: 'mult' }, [el('div', { class: 'mult-title' }, [runTag(r), el('span', { class: 'small', text: '  r = ' + f3(p.pearson) + ' · ' + N(p.n) + ' structures' })])]);
          row.appendChild(m);
          var W = Math.min(300, width(m)), S = W - 40, nb = p.grid.length, cw = S / nb, svg = svgBox(m, W, S + 40, 'CPM density');
          var mx = 0;
          p.grid.forEach(function (rw) { rw.forEach(function (v) { if (v > mx) mx = v; }); });
          s('rect', { x: 34, y: 2, width: S, height: S, style: fill('var(--surface-2)') }, svg);
          p.grid.forEach(function (rw, yi) {
            rw.forEach(function (v, xi) {
              if (!v) return;
              var k = Math.max(1, seqIdx(Math.log(v + 1) / Math.log(mx + 1)));
              s('rect', { x: 34 + xi * cw, y: 2 + S - (yi + 1) * cw, width: cw + 0.3, height: cw + 0.3, style: fill('var(--seq-' + k + ')'), 'data-tip': N(v) + ' structures' }, svg);
            });
          });
          s('line', { x1: 34, y1: 2 + S, x2: 34 + S, y2: 2, class: 'refl' }, svg);
          var nt = Math.max(1, Math.floor(lm));
          for (var t = 0; t <= nt; t++) {
            s('text', { x: 34 + S * t / lm, y: S + 16, 'text-anchor': 'middle', class: 'tick', text: t }, svg);
            s('text', { x: 28, y: 2 + S - S * t / lm, 'text-anchor': 'end', 'dominant-baseline': 'central', class: 'tick', text: t }, svg);
          }
          s('text', { x: 34 + S / 2, y: S + 34, 'text-anchor': 'middle', class: 'axlab', text: trunc(lab(b), 24) + ' log10 CPM+1' }, svg);
        });
        addTable(body, ['Run', 'Shared structures', 'Pearson r', 'Spearman ρ', 'Baseline CPM on shared', 'Run CPM on shared'], ok.map(function (r) {
          var p = quantPair(r);
          return [lab(r), N(p.n), f3(p.pearson), f3(p.spearman), pct(100 * p.cpm_shared_a), pct(100 * p.cpm_shared_b)];
        }));
      });
    var g2 = grid(root);
    if (!others().length) off(g2, 'TSS and TTS agreement', NEEDS.others);
    else panel(g2, { title: 'TSS and TTS agreement', runs: others(), dek: 'Jaccard of 5′ and 3′ end clusters against the baseline (same strand, ends within ' + ((XS().end_agreement || {}).window || 50) + ' bp linked), next to the intron-chain Jaccard.' },
      endPair, NEEDS.endsAgree, function (body, ok) {
        hbars(body, { runs: ok, cats: ['Intron chains', 'TSS clusters', 'TTS clusters'], max: 1, labelW: 104, fmt: f2,
          vals: byRun(ok, function (r) { var p = endPair(r), q = concPair(r); return [q ? q.jaccard : null, p.tss.jaccard, p.tes.jaccard]; }),
          tip: function (r, i) {
            var p = i === 0 ? concPair(r) : endPair(r)[i === 1 ? 'tss' : 'tes'];
            return f3(p.jaccard) + ' · ' + N(p.shared) + ' shared, ' + N(p.only_a) + ' only baseline, ' + N(p.only_b) + ' only run';
          } });
      });
    if (!others().length) off(g2, 'TUSCO incremental recall', NEEDS.others);
    else panel(g2, { title: 'TUSCO incremental recall', runs: others(), dek: 'TUSCO genes a run finds (TP or PTP) that the baseline misses (FN), found by both, and found by the baseline but missed by the run.' },
      tuscoPair, NEEDS.tuscoInc, function (body, ok) {
        hbars(body, { runs: ok, cats: ['Found only here', 'Found by both', 'Missed by this run'], fmt: N, labelW: 126,
          vals: byRun(ok, function (r) { var p = tuscoPair(r); return [p.gained, p.both, p.lost]; }),
          tip: function (r, i) { var p = tuscoPair(r); return i === 0 ? N(p.gained) + ' genes' + names(p.gained_genes) : i === 2 ? N(p.lost) + ' genes' + names(p.lost_genes) : N(p.both) + ' genes'; } });
      });
  }

  function refLine(runs, k) {
    var seen = {}, out = [];
    runs.forEach(function (r) {
      var f = get(B(r), 'sections.ends.refs.' + k);
      if (!f) return;
      var key = (f.file || '') + '|' + (f.md5 || '');
      if (!seen[key]) { seen[key] = 1; out.push(f); }
    });
    return out;
  }
  function viewEnds(root) {
    var runs = sel(), b = baseId();
    section(root, 'Ends', 'How close each isoform’s 5′ and 3′ ends sit to known sites: refTSS for transcription start sites, PolyA_DB for cleavage sites. Distances are signed in transcript orientation: negative means the isoform end is upstream (5′) of the site.');
    root.appendChild(el('p', { class: 'callout' }, [el('strong', { text: 'Precision only. ' }),
      'A database site near an end supports it, so the share of ends near a site is a precision. These databases pool many tissues, so a missing site does not make an end wrong, and the report never scores recall against them.']));
    var tssRuns = withData(runs, function (r) { return endsHalf(r, 'tss'); });
    var tesRuns = withData(runs, function (r) { return endsHalf(r, 'tes'); });
    if (!tssRuns.length && !tesRuns.length) {
      off(root, 'Ends against refTSS and PolyA_DB', NEEDS.ends);
    } else {
      [['tss', 'refTSS'], ['pas', 'PolyA_DB']].forEach(function (k) {
        var refs = refLine(runs, k[0]);
        if (refs.length > 1) root.appendChild(el('p', { class: 'callout' }, [el('strong', { text: 'Not like for like. ' }), 'Runs were scored against different ' + k[1] + ' files: ' + refs.map(function (f) { return f.file + (f.md5 ? ' (md5 ' + f.md5.slice(0, 8) + ')' : ''); }).join(', ') + '.']));
      });
      var top = card(root, { title: 'Precision at 50 bp', wide: true, dek: 'Share of isoforms whose 5′ end is within 50 bp of a refTSS site (left) and whose 3′ end is within 50 bp of a PolyA_DB site (right).' });
      var list = [];
      runs.forEach(function (r) {
        [['tss', '5′ · refTSS'], ['tes', '3′ · PolyA_DB']].forEach(function (k) {
          var h = endsHalf(r, k[0]), v = within(r, k[0], 50), w = within(b, k[0], 50);
          var sub = h ? N(h.within[windowsOf(r).indexOf(50)]) + ' of ' + N(h.n) + (r !== b && isNum(v) && isNum(w) ? ' · ' + mdelta({ k: 'p' }, v, w) + ' vs baseline' : '') : NEEDS[k[0]];
          list.push([el('span', { class: 'runtag' }, [swatch(r), trunc(lab(r), 20) + ' · ' + k[1]]), h ? pct(v) : '–', sub]);
        });
      });
      tiles(top, list);
      var g = grid(root);
      [['tss', '5′ precision as the window widens', 'refTSS'], ['tes', '3′ precision as the window widens', 'PolyA_DB']].forEach(function (k) {
        panel(g, { title: k[1], dek: 'Share of isoforms with an end within the window of a ' + k[2] + ' site. Dashed: 50 bp.' }, function (r) { return endsHalf(r, k[0]); }, NEEDS[k[0]], function (body, ok) {
          var ws = [];
          ok.forEach(function (r) { windowsOf(r).forEach(function (w) { if (ws.indexOf(w) < 0) ws.push(w); }); });
          ws.sort(function (x, y) { return x - y; });
          var ser = function (r) { var h = endsHalf(r, k[0]), own = windowsOf(r); return ws.map(function (w) { var i = own.indexOf(w); return i >= 0 ? pc(h.within[i], h.n) : null; }); };
          var all = [];
          ok.forEach(function (r) { all = all.concat(ser(r)); });
          lines(body, { runs: ok, xs: ws.map(String), xname: 'Window (bp)', series: ser, ymax: maxOf(all), cap: 100, yfmt: tickPct, xlab: 'Window (± bp)', ref: ws.indexOf(50) >= 0 ? ws.indexOf(50) : null,
            tip: function (r, i) { var h = endsHalf(r, k[0]), j = windowsOf(r).indexOf(ws[i]); return pct(pc(h.within[j], h.n)) + ' (' + N(h.within[j]) + ' of ' + N(h.n) + ')'; } });
        });
      });
      [['tss', 'Distance to the nearest refTSS site', '5′ end to refTSS'], ['tes', 'Distance to the nearest PolyA_DB site', '3′ end to PolyA_DB']].forEach(function (k) {
        panel(g, { title: k[1], dek: 'Signed distance, isoforms with a same-strand site within 10 kb.' }, function (r) { var h = endsHalf(r, k[0]); return h && h.dist_hist ? h : null; }, NEEDS[k[0]], function (body, ok) {
          var E = endsHalf(ok[0], k[0]).dist_hist.edges, db = distBins(E);
          vhist(body, { runs: ok, bins: db.full, ticks: db.ticks, mark: db.zero, xname: 'Distance', counts: byRun(ok, function (r) { return endsHalf(r, k[0]).dist_hist.counts; }), xlab: k[2] + ' (bp, bin start)' });
        });
      });
      [['fsm_tss_hist', 'FSM 5′ end vs the reference transcript', 'Distance to reference TSS'], ['fsm_tes_hist', 'FSM 3′ end vs the reference transcript', 'Distance to reference TES']].forEach(function (k) {
        panel(g, { title: k[1], dek: 'Signed distance from each FSM isoform’s end to its matched reference transcript’s end.' },
          function (r) { var h = get(B(r), 'sections.ends.' + k[0]); return h && h.counts && sum(h.counts) > 0 ? h : null; }, NEEDS.ends, function (body, ok) {
            var E = get(B(ok[0]), 'sections.ends.' + k[0]).edges, db = distBins(E);
            vhist(body, { runs: ok, bins: db.full, ticks: db.ticks, mark: db.zero, xname: 'Distance', counts: byRun(ok, function (r) { return get(B(r), 'sections.ends.' + k[0]).counts; }), xlab: k[2] + ' (bp, bin start)' });
          });
      });
      panel(g, { title: 'PolyA signal before the 3′ end', dek: 'Isoforms with a polyA hexamer (AATAAA, ATTAAA or one of ten variants) in the 50 nt upstream of the 3′ end.' },
        function (r) { var m = get(B(r), 'sections.ends.polya_motif'); return m && isNum(m.n) ? m : null; }, NEEDS.motif, function (body, ok) {
          hbars(body, { runs: ok, cats: ['With a polyA motif'], max: 100, labelW: 120, fmt: function (v) { return pct(v); },
            vals: byRun(ok, function (r) { var m = get(B(r), 'sections.ends.polya_motif'); return [pc(m.n_with_motif, m.n)]; }),
            tip: function (r) { var m = get(B(r), 'sections.ends.polya_motif'); return N(m.n_with_motif) + ' of ' + N(m.n) + ' (' + pct(pc(m.n_with_motif, m.n)) + ')'; } });
        });
      var refs = [['tss', 'refTSS'], ['pas', 'PolyA_DB']].map(function (k) {
        return refLine(runs, k[0]).map(function (f) { return k[1] + ': ' + f.file + (isNum(f.sites) ? ' (' + N(f.sites) + ' sites)' : '') + (f.md5 ? ', md5 ' + f.md5.slice(0, 8) : ''); }).join('; ');
      }).filter(function (x) { return x; });
      if (refs.length) note(root, 'References: ' + refs.join(' · ') + '.');
    }
  }

  function viewReads(root) {
    var runs = sel();
    section(root, 'Reads', 'Where 3′ ends sit on genomic A-runs, and how many reads each run assigns to its isoforms.');
    var g = grid(root);
    panel(g, { title: 'Genomic A content after the 3′ end', dek: 'Share of A in the 20 nt of genome just downstream of each 3′ end. High values point to internal priming on a genomic A-run. Dashed: 60%.' },
      function (r) { var p = priming(r); return p && p.down_a_hist ? p : null; }, NEEDS.priming, function (body, ok) {
        var h = priming(ok[0]).down_a_hist, n = h.counts.length, e = h.edges, bins = [], ticks = [];
        for (var i = 0; i < n; i++) {
          var a = Math.round(100 * e[i]), z = e[i + 1];
          if (z === undefined) { bins.push(a + '% A'); ticks.push(String(a)); } else { bins.push(a + '–' + Math.round(100 * z) + '% A'); ticks.push(String(a)); }
        }
        var mark = null;
        e.forEach(function (x, i) { if (Math.abs(x - 0.6) < 1e-9) mark = i; });
        vhist(body, { runs: ok, bins: bins, ticks: ticks, mark: mark, xname: 'A downstream', counts: byRun(ok, function (r) { return priming(r).down_a_hist.counts; }), xlab: '% A downstream (bin start)' });
      });
    panel(g, { title: 'Internal-priming risk', dek: 'Isoforms with at least 60% A in the 20 nt after the 3′ end.' }, priming, NEEDS.priming, function (body, ok) {
      hbars(body, { runs: ok, cats: ['≥ 60% A downstream'], fmt: function (v) { return pct(v); }, labelW: 136,
        vals: byRun(ok, function (r) { var p = priming(r); return [pc(p.risk, p.n)]; }),
        tip: function (r) { var p = priming(r); return N(p.risk) + ' of ' + N(p.n) + ' isoforms (' + pct(pc(p.risk, p.n)) + ')'; } });
    });
    panel(g, { title: 'Reads assigned', dek: 'Reads FLAIR quantify assigns to the run’s isoforms, and isoforms with at least one read.' }, quant, NEEDS.quant, function (body, ok) {
      quantNote(body, ok);
      hbars(body, { runs: ok, cats: ['Reads assigned'], fmt: N, labelW: 130,
        vals: byRun(ok, function (r) { return [quant(r).reads_assigned]; }) });
      hbars(body, { runs: ok, cats: ['Isoforms with reads'], fmt: N, labelW: 130,
        vals: byRun(ok, function (r) { return [quant(r).isoforms_quantified]; }) });
    });
    panel(g, { title: 'Isoform expression', dek: 'Isoforms with reads, by CPM.' }, function (r) { var q = quant(r); return q && q.cpm_hist && q.cpm_hist.counts ? q : null; }, NEEDS.quant, function (body, ok) {
      var h = quant(ok[0]).cpm_hist;
      quantNote(body, ok);
      vhist(body, { runs: ok, bins: binNames(h.edges, h.counts.length, false).map(function (x) { return x + ' CPM'; }), ticks: h.edges.slice(0, h.counts.length).map(function (x, i) { return i === h.edges.length - 1 ? '≥' + kfmt(x) : kfmt(x); }), xname: 'CPM', counts: byRun(ok, function (r) { return quant(r).cpm_hist.counts; }), xlab: 'CPM (bin start)' });
    });
  }

  function viewCompute(root) {
    section(root, 'Compute', 'FLAIR tasks of each run (assemble, quantify, productivity), from the Nextflow trace. The shared align step and the evaluation are not counted.');
    var g = grid(root, 'three');
    [['wall_minutes', 'Wall time', fmtWall], ['cpu_hours', 'CPU time', fmtCpu], ['peak_rss_gb', 'Peak memory', fmtMem]].forEach(function (m) {
      panel(g, { title: m[1] }, function (r) { return compute(r, m[0]); }, NEEDS.compute, function (body, ok) {
        hbars(body, { runs: ok, cats: [''], vals: byRun(ok, function (r) { return [compute(r, m[0])]; }), fmt: m[2], labelW: 4, valW: 72 });
      });
    });
  }

  var VIEWS = [
    { id: 'diff', label: 'What differs', fn: viewDiff },
    { id: 'structure', label: 'Structure', fn: viewStructure },
    { id: 'recovery', label: 'Recovery', fn: viewRecovery },
    { id: 'diversity', label: 'Diversity', fn: viewDiversity },
    { id: 'agreement', label: 'Agreement', fn: viewAgreement },
    { id: 'ends', label: 'Ends', fn: viewEnds },
    { id: 'reads', label: 'Reads', fn: viewReads },
    { id: 'compute', label: 'Compute', fn: viewCompute }
  ];

  // ------------------------------------------------------------------ shell
  // Deep link: #tab or #tab/sample.
  function hashFor() { return '#' + state.tab + (state.sample !== SAMPLES[0] ? '/' + encodeURIComponent(state.sample) : ''); }
  function fromHash() {
    var h = (location.hash || '').slice(1).split('/'), changed = false;
    if (VIEW_IDS.indexOf(h[0]) >= 0 && h[0] !== state.tab) { state.tab = h[0]; changed = true; }
    var smp = h[1] ? decodeURIComponent(h[1]) : null;
    if (smp && SAMPLES.indexOf(smp) >= 0 && smp !== state.sample) { state.sample = smp; changed = true; }
    return changed;
  }
  function setHash() { try { history.replaceState(null, '', hashFor()); } catch (e) { /* file:// in some browsers */ } }
  function go(id) {
    state.tab = id;
    setHash();
    render();
    try { window.scrollTo(0, 0); } catch (e) { /* ignore */ }
  }
  function renderControls() {
    var c = document.getElementById('controls');
    c.textContent = '';
    var pick = function (v) { state.sample = v; setHash(); render(); };
    c.appendChild(el('div', { class: 'ctl' }, [el('span', { class: 'ctl-label', text: 'Sample' }),
      SAMPLES.length <= 4 ? segment('Sample', SAMPLES.map(function (x) { return [x, x]; }), state.sample, pick) : picker('Sample', SAMPLES.map(function (x) { return [x, x]; }), state.sample, pick)]));
    var here = sampleRuns(), b = baseId();
    if (here.length > 1) {
      var chips = el('div', { class: 'ctl' }, [el('span', { class: 'ctl-label', text: 'Compare' })]);
      here.forEach(function (r) {
        var on = r === b || !state.hidden[r];
        chips.appendChild(el('button', { type: 'button', class: 'chip', 'aria-pressed': String(on), title: r === b ? 'The baseline is always shown' : null,
          onclick: function () { if (r === b) return; if (on) state.hidden[r] = 1; else delete state.hidden[r]; render(); } }, [swatch(r), lab(r)]));
      });
      c.appendChild(chips);
      var bopts = here.map(function (r) { return [r, lab(r) + (r === DATA.baseline ? ' (default)' : '')]; });
      var setB = function (v) { state.base = v; delete state.hidden[v]; render(); };
      c.appendChild(el('div', { class: 'ctl' }, [el('span', { class: 'ctl-label', text: 'Baseline' }),
        here.length <= 3 ? segment('Baseline', here.map(function (r) { return [r, trunc(lab(r), 18)]; }), b, setB) : picker('Baseline', bopts, b, setB)]));
      if (sel().length > 1) c.appendChild(el('div', { class: 'ctl' }, [el('span', { class: 'ctl-label', text: 'Layout' }),
        segment('Layout', [['overlay', 'Overlay'], ['side', 'Side by side']], state.layout, function (v) { state.layout = v; render(); })]));
    }
  }
  function renderBaseNote() {
    var n = document.getElementById('basenote'), b = baseId();
    n.textContent = '';
    n.appendChild(el('span', { class: 'ctl-label', text: 'Baseline' }));
    n.appendChild(runTag(b, false, true));
    var why = state.base !== b ? state.sample + ' has no run ' + lab(state.base) + ', so its first run stands in' : baselineWhy(b);
    n.appendChild(document.createTextNode(' · ' + why + '. Every comparison is against it.'));
    if (b !== DATA.baseline && sampleRuns().indexOf(DATA.baseline) >= 0) {
      n.appendChild(el('button', { type: 'button', class: 'linkbtn no-print', text: 'Use default: ' + lab(DATA.baseline), onclick: function () { state.base = DATA.baseline; render(); } }));
    }
  }
  function renderTabs() {
    var t = document.getElementById('tabs');
    t.textContent = '';
    VIEWS.forEach(function (v, i) {
      var on = v.id === state.tab;
      t.appendChild(el('button', { type: 'button', role: 'tab', id: 'tab-' + v.id, 'aria-selected': String(on), tabindex: on ? '0' : '-1', 'aria-controls': 'view',
        onclick: function () { go(v.id); },
        onkeydown: function (e) {
          var j = e.key === 'ArrowRight' ? i + 1 : e.key === 'ArrowLeft' ? i - 1 : e.key === 'Home' ? 0 : e.key === 'End' ? VIEWS.length - 1 : null;
          if (j === null) return;
          e.preventDefault();
          go(VIEWS[(j + VIEWS.length) % VIEWS.length].id);
          var nb = document.getElementById('tab-' + state.tab);
          if (nb && nb.focus) nb.focus();
        }, text: v.label }));
    });
  }
  function renderHead() {
    var when = DATA.generated ? DATA.generated.slice(0, 10) : '';
    document.getElementById('subtitle').textContent = DATA.bundles.length + ' bundles · ' + RUN_ORDER.length + ' runs · ' + SAMPLES.length + ' sample' + (SAMPLES.length === 1 ? '' : 's') +
      (when ? ' · built ' + when : '') + (DATA.tool_version ? ' · flair-eval ' + DATA.tool_version : '');
    document.getElementById('footer').textContent = 'Colour means run: each run keeps its colour in every tab, from its place in the run order; grey is shared or other. One sample is shown at a time and nothing is pooled across samples.';
  }
  function render() {
    hideTip();
    setSampleNeeds();
    renderControls(); renderBaseNote(); renderTabs();
    var root = document.getElementById('view');
    root.textContent = '';
    if (state.printing) {
      VIEWS.forEach(function (v) { var sct = el('div', { class: 'print-section' }); root.appendChild(sct); v.fn(sct); });
      return;
    }
    var v = VIEWS.find(function (x) { return x.id === state.tab; }) || VIEWS[0];
    root.setAttribute('aria-labelledby', 'tab-' + v.id);
    v.fn(root);
  }

  // ------------------------------------------------------------------ theme, print, tooltip
  function setTheme(t) {
    var de = document.documentElement;
    if (t === 'auto') de.removeAttribute('data-theme'); else de.setAttribute('data-theme', t);
    try { localStorage.setItem('flair-eval-theme', t); } catch (e) { /* storage blocked */ }
    renderTheme(t);
  }
  function renderTheme(t) {
    var box = document.getElementById('theme');
    box.textContent = '';
    box.appendChild(segment('Colour theme', [['auto', 'Auto'], ['light', 'Light'], ['dark', 'Dark']], t, setTheme));
  }
  var theme0 = 'auto';
  try { theme0 = localStorage.getItem('flair-eval-theme') || 'auto'; } catch (e) { /* storage blocked */ }
  if (theme0 !== 'auto') document.documentElement.setAttribute('data-theme', theme0);
  renderTheme(theme0);
  document.getElementById('print').addEventListener('click', function () {
    state.printing = true; render();
    try { window.print(); } finally { state.printing = false; render(); }
  });

  var tip = document.getElementById('tooltip');
  function hideTip() { tip.hidden = true; }
  function showTip(e) {
    var t = e.target && e.target.closest ? e.target.closest('[data-tip]') : null;
    if (!t) { hideTip(); return; }
    tip.textContent = t.getAttribute('data-tip');
    tip.hidden = false;
    var x = Math.min(e.clientX + 12, window.innerWidth - tip.offsetWidth - 8), y = e.clientY + 14;
    if (y + tip.offsetHeight > window.innerHeight - 8) y = e.clientY - tip.offsetHeight - 10;
    tip.style.left = Math.max(8, x) + 'px';
    tip.style.top = y + 'px';
  }
  document.addEventListener('pointermove', showTip);
  document.addEventListener('pointerdown', showTip);
  document.addEventListener('scroll', hideTip, true);

  var lastW = window.innerWidth, timer = null;
  window.addEventListener('resize', function () {
    if (Math.abs(window.innerWidth - lastW) < 8) return;
    lastW = window.innerWidth; clearTimeout(timer); timer = setTimeout(render, 120);
  });
  window.addEventListener('hashchange', function () { if (fromHash()) render(); });
  fromHash();
  renderHead();
  render();
  // Test hook: lets tests drive every tab, sample and layout without a browser.
  window.flairEvalReport = { state: state, render: render, views: VIEW_IDS };
})();
