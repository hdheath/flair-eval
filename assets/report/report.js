/* flair-eval report: the transcriptome viewer (CONTRACTS section 11).
 * Tabs: What differs, Structure, Recovery, Diversity, Agreement, Ends, Quantification (id 'reads'), Compute.
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
  var state = { tab: 'diff', sample: SAMPLES[0], base: DATA.baseline, hidden: {}, layout: 'overlay', printing: false, qcCat: 'ISM' };

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
  // One value per run, one labelled bar each (no colour lookup needed). o: runs, val(run), fmt,
  // tip(run), max, name (table column), label
  function runBars(box, o) {
    var runs = o.runs || sel(), viz = el('div', { class: 'viz' });
    box.appendChild(viz);
    var W = width(viz), lw = Math.min(o.labelW || 160, W * 0.42), vw = o.valW || 72, bh = 14, gap = 8;
    var vals = runs.map(function (r) { return o.val(r); });
    var max = o.max || niceMax(maxOf(vals)), H = runs.length * (bh + gap);
    var svg = svgBox(viz, W, H, o.label);
    runs.forEach(function (r, i) {
      var v = vals[i], y = i * (bh + gap) + 2, w = isNum(v) ? Math.max(1, (W - lw - vw) * Math.min(1, v / max)) : 0;
      s('text', { x: lw - 8, y: y + bh / 2, 'text-anchor': 'end', 'dominant-baseline': 'central', class: 'lab', text: trunc(lab(r), Math.floor((lw - 8) / 6.5)) }, svg);
      if (!isNum(v)) return;
      s('rect', { x: lw, y: y, width: w, height: bh, rx: 2, style: fill(col(r)), 'data-tip': lab(r) + ': ' + (o.tip ? o.tip(r) : o.fmt(v)) }, svg);
      s('text', { x: lw + w + 5, y: y + bh / 2, 'dominant-baseline': 'central', class: 'val', text: o.fmt(v) }, svg);
    });
    addTable(box, ['Run', o.name || ''], runs.map(function (r, i) { return [lab(r), isNum(vals[i]) ? (o.tip ? o.tip(r) : o.fmt(vals[i])) : '–']; }));
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
    NEEDS.tusco = empty ? 'no TUSCO gene lies inside the regions given with --regions' : 'needs the TUSCO gene list (--genome hg38, or --tusco)';
    NEEDS.tuscoInc = empty ? NEEDS.tusco : 'needs the TUSCO gene list for both the baseline and this run';
  }
  function prodCounts(r) { var p = sec(r, 'productivity'); return p && p.class_counts ? p.class_counts : null; }
  function prodShare(r, k) { var cc = prodCounts(r); if (!cc) return null; var t = sum(Object.keys(cc).map(function (x) { return cc[x]; })); return pc(cc[k] || 0, t); }
  function variation(r) { var v = sec(r, 'variation'); return v && isNum(v.genes_multi) ? v : null; }
  function priming(r) { var p = sec(r, 'priming'); return p && isNum(p.n) ? p : null; }
  // Section 13 (A42-A44): one key of sections.structure_qc, or null.
  function sqc(r, k) { var q = sec(r, 'structure_qc'); return q && q[k] ? q[k] : null; }
  function jsup(r) { var j = sqc(r, 'junction_support'); return j && isNum(j.junctions) ? j : null; }
  // A45: support in the sample's sr_junctions file (short reads, scoring only).
  function srsup(r) { var j = sqc(r, 'sr_junction_support'); return j && isNum(j.junctions) ? j : null; }
  // False when FLAIR may have seen the short-read file: it is also the sample's junctions file
  // (same md5) or the run's own junctions file (same name).
  function srIndependent(r) {
    var s = srsup(r), j = jsup(r), inj = ((B(r) || {}).inputs || {}).junctions;
    if (!s) return false;
    if (j && j.md5 && s.md5 && j.md5 === s.md5) return false;
    return !(inj && s.file && String(inj).split('/').pop() === s.file);
  }
  function quant(r) { var q = sec(r, 'quant'); return q && isNum(q.reads_assigned) ? q : null; }
  function compute(r, k) { var c = (B(r) || {}).compute; return c && isNum(c[k]) ? c[k] : null; }
  // End accuracy (A40): sections.end_accuracy, every number from ted_end_precision.py. Fractions
  // (0-1) in the section; shown as percentages, F1 as 0-1.
  function endAcc(r) { var e = sec(r, 'end_accuracy'); return e && e.ends && e.ends.length ? e : null; }
  function eaRow(r, id) { var e = endAcc(r); return e ? e.ends.find(function (x) { return x.set === id; }) || null : null; }
  function eaPair(r, a, z) { var e = endAcc(r); return e && e.pairs ? e.pairs.find(function (p) { return p.tss_set === a && p.tes_set === z; }) || null : null; }
  function asP(v) { return isNum(v) ? 100 * v : null; }
  // End sets over the runs, by id: 5' before 3', databases before truth sets, then section order.
  function eaSets(runs) {
    var seen = {}, out = [], ek = { tss: 0, tes: 1 }, kk = { reference: 0, truth: 1 };
    runs.forEach(function (r) {
      var e = endAcc(r);
      if (e) (e.sets || []).forEach(function (x) { if (!seen[x.id]) { seen[x.id] = 1; out.push(x); } });
    });
    return out.map(function (x, i) { return [x, i]; }).sort(function (a, z) {
      return ek[a[0].end] - ek[z[0].end] || kk[a[0].kind] - kk[z[0].kind] || a[1] - z[1];
    }).map(function (x) { return x[0]; });
  }
  function eaPairs(runs) {
    var seen = {}, out = [];
    runs.forEach(function (r) {
      var e = endAcc(r);
      if (e) (e.pairs || []).forEach(function (p) { var k = p.tss_set + '|' + p.tes_set; if (!seen[k]) { seen[k] = 1; out.push(p); } });
    });
    return out.sort(function (a, z) { return (a.kind === 'reference' ? 0 : 1) - (z.kind === 'reference' ? 0 : 1); });
  }
  function endWord(e) { return e === 'tss' ? '5′' : '3′'; }

  // One-line reasons for hidden panels, in the reader's terms: pipeline params, never section keys.
  var OLDER = 'not in these results (made by an older flair-eval); rerun to see it';
  var NEEDS = {
    ends: 'needs refTSS or PolyA_DB sites (--genome hg38, or --ref_tss / --ref_pas)',
    truth: 'needs truth sets for this sample (--truth)',
    tss: 'needs refTSS sites (--genome hg38, or --ref_tss)',
    tes: 'needs PolyA_DB sites (--genome hg38, or --ref_pas)',
    motif: 'needs refTSS or PolyA_DB sites and the genome',
    tusco: 'needs the TUSCO gene list (--genome hg38, or --tusco)',
    priming: OLDER,
    variation: OLDER,
    quant: 'the run has no read counts (FLAIR quantify did not run or failed)',
    productivity: 'the run has no coding predictions (FLAIR productivity did not run or failed)',
    conc: 'needs the isoform files of at least two runs next to their results',
    cpm: 'needs read counts for both the baseline and this run',
    endsAgree: 'needs the isoform files of at least two runs next to their results',
    tuscoInc: 'needs the TUSCO gene list for both the baseline and this run',
    lostDec: 'needs the isoform files, and read counts for the baseline',
    compute: 'no run-time record (the Nextflow trace was not found)',
    hist: OLDER,
    structureQc: OLDER,
    srJunctionFile: 'needs short-read junctions for this sample (the sr_junctions column of samples.csv)',
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
  // Share of reference transcripts at least minLen long that have an FSM (structure_qc.ref_recovery).
  function longRecovery(r, minLen) {
    var h = get(sqc(r, 'ref_recovery'), 'length');
    if (!h) return null;
    var a = 0, f = 0;
    h.edges.forEach(function (e, i) { if (e >= minLen) { a += h.all[i]; f += h.fsm[i]; } });
    return pc(f, a);
  }
  var METRICS = [];
  var FAMS = [];
  function family(fam, tab, rows) { FAMS.push({ fam: fam, tab: tab }); rows.forEach(function (m) { m.fam = fam; m.tab = tab; METRICS.push(m); }); }
  family('Structure', 'structure', [
    { name: 'Isoforms', top: true, k: 'n', lowerIsWorse: null, get: function (r) { return B(r).summary.isoforms; } },
    { name: 'FLAIR genes', k: 'n', lowerIsWorse: null, get: function (r) { return B(r).summary.genes; } },
    { name: 'Mono-exon isoforms', k: 'p', lowerIsWorse: null, get: function (r) { var x = B(r).summary; return pc(x.mono_exon, x.isoforms); } },
    { name: 'FSM isoforms', top: true, k: 'p', lowerIsWorse: null, get: function (r) { var x = B(r); return pc(x.categories.FSM, x.summary.isoforms); } },
    { name: 'Novel-junction isoforms (NIC + NNC)', top: true, k: 'p', lowerIsWorse: null, get: function (r) { var x = B(r); return pc((x.categories.NIC || 0) + (x.categories.NNC || 0), x.summary.isoforms); } },
    { name: 'Known splice junctions', top: true, k: 'p', lowerIsWorse: true, get: function (r) { var j = B(r).junctions; return j ? pc(j.known, j.total) : null; } },
    { name: 'ISM missing 5′ exons (3′ fragment)', k: 'p', lowerIsWorse: false, needs: 'structureQc', get: function (r) { var q = sqc(r, 'subcategories'); return q ? pc(q.ISM['3prime_fragment'], B(r).categories.ISM) : null; } },
    { name: 'Multi-exon isoforms with only canonical junctions', k: 'p', lowerIsWorse: true, needs: 'structureQc', get: function (r) { var c = sqc(r, 'canonical'); return c ? pc(c.all_canonical, c.multi_exon) : null; } },
    { name: 'RT-switching junctions', k: 'p', lowerIsWorse: false, needs: 'structureQc', get: function (r) { var t = sqc(r, 'rts'); return t ? pc(t.rts, t.junctions) : null; } },
    { name: 'Junctions supported by short reads', top: true, k: 'p', lowerIsWorse: true, needs: 'srJunctionFile', get: function (r) { var j = srsup(r); return j ? pc(j.supported, j.junctions) : null; } },
    { name: 'Multi-exon isoforms with every junction supported by short reads', k: 'p', lowerIsWorse: true, needs: 'srJunctionFile', get: function (r) { var j = srsup(r); return j ? pc(j.all_supported, j.multi_exon) : null; } }
  ]);
  // noRank: listed, but kept out of "Biggest differences" (it repeats the exact-ends row on Ends).
  // sameQuant: ranked only when the run and the baseline counted reads the same way (A21).
  family('Recovery', 'recovery', [
    { name: 'Reference genes detected', top: true, k: 'n', lowerIsWorse: true, get: function (r) { return B(r).summary.ref_genes_detected; } },
    { name: 'Reference transcripts with an FSM', k: 'n', lowerIsWorse: true, get: function (r) { return B(r).summary.ref_transcripts_fsm; } },
    { name: 'Reference transcripts with exact ends', top: true, noRank: true, k: 'n', lowerIsWorse: true, get: function (r) { return B(r).summary.ref_transcripts_cim; } },
    { name: 'TUSCO sensitivity (genes)', top: true, k: 'p', lowerIsWorse: true, needs: 'tusco', get: function (r) { var t = tusco(r); return t ? asPct(t.sensitivity) : null; } },
    { name: 'TUSCO precision (isoforms)', top: true, k: 'p', lowerIsWorse: true, needs: 'tusco', get: function (r) { var t = tusco(r); return t ? asPct(t.precision) : null; } },
    { name: 'TUSCO genes missed', k: 'n', lowerIsWorse: false, needs: 'tusco', get: function (r) { var t = tusco(r); return t ? t.FN : null; } },
    { name: 'Genes lost vs baseline', k: 'n', lowerIsWorse: false, vsBase: true, needs: 'conc', get: function (r) { var g = genesVs(r); return g ? g.lost.length : null; } },
    { name: 'Genes gained vs baseline', k: 'n', lowerIsWorse: null, vsBase: true, needs: 'conc', get: function (r) { var g = genesVs(r); return g ? g.gained.length : null; } },
    { name: 'Reference transcripts of 3 kb or more with an FSM', k: 'p', lowerIsWorse: true, needs: 'structureQc', get: function (r) { return longRecovery(r, 3000); } }
  ]);
  family('Diversity', 'diversity', [
    { name: 'Isoforms per gene (mean)', k: 'f2', lowerIsWorse: null, get: function (r) { var x = B(r).summary; return x.genes ? x.isoforms / x.genes : null; } },
    { name: 'Single-isoform genes', k: 'p', lowerIsWorse: null, needs: 'hist', get: function (r) { var h = histOf(r, 'isoforms_per_gene'); return h ? pc(h.counts[0], sum(h.counts)) : null; } },
    { name: 'Productive (PRO)', k: 'p', lowerIsWorse: null, needs: 'productivity', get: function (r) { return prodShare(r, 'PRO'); } },
    { name: 'Premature stop, NMD candidate (PTC)', k: 'p', lowerIsWorse: null, needs: 'productivity', get: function (r) { return prodShare(r, 'PTC'); } },
    { name: 'Genes differing by ends only', k: 'p', lowerIsWorse: null, needs: 'variation', get: function (r) { var v = variation(r); return v ? pc(v.ends_only, v.genes_multi) : null; } },
    { name: 'Genes differing by splicing only', k: 'p', lowerIsWorse: null, needs: 'variation', get: function (r) { var v = variation(r); return v ? pc(v.splice_only, v.genes_multi) : null; } },
    { name: 'Reference transcripts with more than one FSM', k: 'p', lowerIsWorse: false, needs: 'structureQc', get: function (r) { var d = sqc(r, 'redundancy'); return d ? pc(d.fsm_refs_multi, d.fsm_refs) : null; } }
  ]);
  family('Agreement', 'agreement', [
    { name: 'Isoforms shared with baseline (Jaccard)', k: 'f2', lowerIsWorse: null, vsBase: true, needs: 'conc', get: function (r) { var p = concPair(r); return p ? p.jaccard : null; } },
    { name: '5′ end Jaccard vs baseline', k: 'f2', lowerIsWorse: null, vsBase: true, needs: 'endsAgree', get: function (r) { var p = endPair(r); return p ? p.tss.jaccard : null; } },
    { name: '3′ end Jaccard vs baseline', k: 'f2', lowerIsWorse: null, vsBase: true, needs: 'endsAgree', get: function (r) { var p = endPair(r); return p ? p.tes.jaccard : null; } },
    { name: 'Expression agreement with baseline (Pearson r)', k: 'f3', lowerIsWorse: null, vsBase: true, needs: 'cpm', get: function (r) { var p = quantPair(r); return p ? p.pearson : null; } },
    { name: 'TUSCO genes found that the baseline misses', k: 'n', lowerIsWorse: null, vsBase: true, needs: 'tuscoInc', get: function (r) { var p = tuscoPair(r); return p ? p.gained : null; } }
  ]);
  family('Ends', 'ends', [
    { name: '5′ end within 50 bp of refTSS', top: true, k: 'p', lowerIsWorse: true, needs: 'tss', get: function (r) { return within(r, 'tss', 50); } },
    { name: '3′ end within 50 bp of PolyA_DB', top: true, k: 'p', lowerIsWorse: true, needs: 'tes', get: function (r) { return within(r, 'tes', 50); } },
    { name: 'PolyA motif before the 3′ end', k: 'p', lowerIsWorse: true, needs: 'motif', get: function (r) { var m = get(B(r), 'sections.ends.polya_motif'); return m ? pc(m.n_with_motif, m.n) : null; } },
    { name: 'FSM isoforms with exact ends', k: 'p', lowerIsWorse: true, get: function (r) { var x = B(r); return pc(x.fsm_subcategories.CIM, x.categories.FSM); } },
    { name: 'Internal-priming risk', top: true, k: 'p', lowerIsWorse: false, needs: 'priming', get: function (r) { var p = priming(r); return p ? pc(p.risk, p.n) : null; } },
    { name: 'Internal-priming risk, multi-exon isoforms', k: 'p', lowerIsWorse: false, needs: 'structureQc', get: function (r) { var p = sqc(r, 'priming'); return p ? pc(p.by_exons.multi.risk, p.by_exons.multi.n) : null; } }
  ]);
  family('Quantification', 'reads', [
    { name: 'Reads assigned', top: true, sameQuant: true, k: 'n', lowerIsWorse: null, needs: 'quant', get: function (r) { var q = quant(r); return q ? q.reads_assigned : null; } },
    { name: 'Isoforms with reads', sameQuant: true, k: 'n', lowerIsWorse: null, needs: 'quant', get: function (r) { var q = quant(r); return q ? q.isoforms_quantified : null; } }
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
  // How each run counted reads (A21: dev counts from the genome BAM, 2.x and 3.0.1 from the reads).
  var QUANT_INPUT = { reads: 'from the reads', bam: 'from the genome alignment (BAM)' };
  function quantInput(r) { var q = sec(r, 'quant'); return q ? q.input || 'unknown' : null; }
  // Mixed inputs get a callout; one shared input gets a one-line note.
  function quantNote(parent, runs) {
    var groups = {}, order = [];
    runs.forEach(function (r) {
      var k = quantInput(r);
      if (!k) return;
      if (!groups[k]) { groups[k] = []; order.push(k); }
      groups[k].push(lab(r));
    });
    if (!order.length) return;
    if (order.length === 1) {
      if (QUANT_INPUT[order[0]]) note(parent, 'All runs counted reads ' + QUANT_INPUT[order[0]] + '.');
      return;
    }
    var parts = order.map(function (k) { return (QUANT_INPUT[k] || 'not recorded') + ' for ' + groups[k].join(', '); });
    parent.appendChild(el('p', { class: 'callout', role: 'note' }, [el('strong', { text: 'Not like for like. ' }),
      'Reads were counted differently: ' + parts.join('; ') + '.']));
  }
  var RULE = DATA.baseline_rule || {};
  function baselineWhy(id) {
    if (id !== DATA.baseline) return 'chosen in this page';
    if (RULE.rule === 'given') return 'set with --baseline';
    if (RULE.rule === 'earliest_launch') return RULE.launches > 1 ? 'first run of the earliest of ' + RULE.launches + ' launches' : 'first run of the launch';
    if (RULE.rule === 'name_order') return 'first run in alphabetical order (the results do not record the run order)';
    return 'default';
  }

  // ------------------------------------------------------------------ views
  function viewDiff(root) {
    var b = baseId();
    section(root, 'What differs', 'Every run against the baseline. Pick another baseline, or hide runs, in the bar above.');
    var gtf = (B(b).inputs || {}).gtf;
    root.appendChild(el('div', { class: 'callout' }, [el('strong', { text: 'Read with care. ' }),
      'Runs were built with the same annotation they are scored against' + (gtf ? ' (' + gtf + ')' : '') +
      ', so FSM counts and gene recovery reward matching it, not being correct. TUSCO, refTSS and PolyA_DB are independent checks; refTSS and PolyA_DB measure precision only.']));
    if (!others().length) note(root, sampleRuns().length > 1 ? 'Every other run is hidden: show one in the bar above to compare.' : 'Only one run on this sample, so there is nothing to compare. Its numbers are in the table below.', 'note');
    if (others().length) whereDiffer(root);
    var g = grid(root);
    others().forEach(function (r) {
      var body = card(g, { title: lab(r) + ' vs ' + lab(b) });
      var p = concPair(r), q = quantPair(r), t = tusco(r), tb = tusco(b);
      var list = [];
      if (p) list.push(['Isoforms shared', N(p.shared), 'Jaccard ' + f2(p.jaccard) + ' (shared ÷ all)']);
      list.push(['Isoforms', N(B(r).summary.isoforms), mdelta(METRICS[0], B(r).summary.isoforms, B(b).summary.isoforms) + ' vs baseline']);
      if (t && tb) list.push(['TUSCO sensitivity', pct(asPct(t.sensitivity)), mdelta({ k: 'p' }, asPct(t.sensitivity), asPct(tb.sensitivity)) + ' vs baseline']);
      else list.push(['Reference genes detected', N(B(r).summary.ref_genes_detected), mdelta({ k: 'n' }, B(r).summary.ref_genes_detected, B(b).summary.ref_genes_detected) + ' vs baseline']);
      [['tss', '5′ within 50 bp of refTSS'], ['tes', '3′ within 50 bp of PolyA_DB']].forEach(function (k) {
        var v = within(r, k[0], 50), w = within(b, k[0], 50);
        if (isNum(v)) list.push([k[1], pct(v), isNum(w) ? mdelta({ k: 'p' }, v, w) + ' vs baseline' : '']);
      });
      if (q && isNum(q.pearson)) list.push(['Expression agreement', 'r = ' + f3(q.pearson), 'Pearson r on log CPM, ' + N(q.n) + ' shared isoforms']);
      tiles(body, list);
      // At most one row per family, so one signal does not fill the list; read counts only when
      // both runs counted reads the same way.
      var famSeen = {};
      var ranked = METRICS.filter(function (m) { return !m.noRank && !(m.sameQuant && quantInput(r) !== quantInput(b)); })
        .map(function (m) { var v = mval(m, r), w = mval(m, b); return { m: m, v: v, w: w, sc: mscore(m, v, w) }; })
        .filter(function (x) { return x.sc > 0.05; }).sort(function (x, y) { return y.sc - x.sc; })
        .filter(function (x) { if (famSeen[x.m.fam]) return false; famSeen[x.m.fam] = 1; return true; }).slice(0, 5);
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
    terms(root);
  }
  var TERMS = [
    ['FSM, ISM', 'full or partial match to every splice junction of a reference transcript'],
    ['NIC, NNC', 'novel isoforms: known junctions in a new combination (NIC), or at least one new junction (NNC)'],
    ['Exact ends', 'an FSM isoform with both ends within 50 bp of the matched reference transcript (CIM in the data files)'],
    ['TUSCO', 'SQANTI3’s list of genes expected to have exactly one isoform; an independent truth set'],
    ['refTSS, PolyA_DB', 'databases of known transcription start sites and poly(A) cleavage sites'],
    ['Jaccard', 'shared ÷ all, from 0 (nothing shared) to 1 (identical)'],
    ['CPM', 'counts per million assigned reads'],
    ['PTC, NMD', 'premature stop codon; such isoforms are candidates for nonsense-mediated decay'],
    ['3′ fragment', 'an ISM isoform that matches the 3′ junctions of a reference transcript but lacks its 5′ exons (SQANTI3 3prime_fragment), the usual sign of a 5′-truncated read'],
    ['RT switching', 'a junction whose donor and acceptor sides share a short repeat, so the reverse transcriptase may have jumped across it and made a false intron (SQANTI3 definition)']
  ];
  function terms(root) {
    var body = card(root, { title: 'Terms', wide: true });
    body.appendChild(el('dl', { class: 'terms' }, [].concat.apply([], TERMS.map(function (t) { return [el('dt', { text: t[0] }), el('dd', { text: t[1] })]; }))));
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
      dek: 'One row per metric. Each row runs from its lowest to its highest value, so read the range at the right for how big a gap is. The outlined dot is the baseline.' });
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
    var body = card(root, { title: 'All metrics', wide: true,
      dek: 'Every metric, one row each, with each run’s difference from the baseline under its value. Click a family to open its plots.' });
    var head = el('tr', {}, [el('th', { text: 'Metric' })].concat(runs.map(function (r) { return el('th', { class: 'num runhead' }, [runTag(r, 16)]); }))
      .concat([el('th', { class: 'num', text: 'Largest difference vs baseline' })]));
    var tbody = el('tbody');
    FAMS.forEach(function (F) {
      var rows = METRICS.filter(function (m) { return m.fam === F.fam; });
      var shown = rows.filter(function (m) { return runs.some(function (r) { return mval(m, r) !== null; }); });
      var hidden = rows.filter(function (m) { return shown.indexOf(m) < 0; });
      tbody.appendChild(el('tr', { class: 'famrow' }, [el('td', { colspan: runs.length + 2 }, [el('button', { type: 'button', class: 'linkish', onclick: function () { go(F.tab); }, text: F.fam + ' →' })])]));
      shown.forEach(function (m) {
        var vals = runs.map(function (r) { return mval(m, r); });
        var tr = el('tr', {}, [el('td', { text: m.name })].concat(vals.map(function (v, i) {
          var txt = m.vsBase && runs[i] === b ? 'baseline' : mfmt(m, v);
          var d = i ? mdelta(m, v, vals[0]) : '';
          return el('td', { class: 'num' + (m.vsBase && runs[i] === b ? ' muted' : '') }, [txt, d ? el('span', { class: 'cell-delta', text: d }) : null]);
        })));
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
        tbody.appendChild(el('tr', { class: 'hiddenrow' }, [el('td', { colspan: runs.length + 2, text: 'Hidden: ' + Object.keys(why).map(function (k) { return why[k].join(', ') + ' (' + k + ')'; }).join('; ') + '.' })]));
      }
    });
    // Closed by default: the dot rows and per-run cards above carry the headline; print opens it.
    body.appendChild(el('details', { class: 'more', open: state.printing ? '' : null }, [el('summary', { text: 'Show the table' }),
      el('div', { class: 'tablewrap' }, [el('table', { class: 'data metrics' }, [el('thead', {}, [head]), tbody])])]));
  }

  var CATS = ['FSM', 'ISM', 'NIC', 'NNC', 'genic', 'genic_intron', 'antisense', 'intergenic', 'fusion'];
  var CAT_LABEL = { FSM: 'FSM', ISM: 'ISM', NIC: 'NIC', NNC: 'NNC', genic: 'Genic', genic_intron: 'Genic intron', antisense: 'Antisense', intergenic: 'Intergenic', fusion: 'Fusion' };
  function byRun(runs, f) { var o = {}; runs.forEach(function (r) { o[r] = f(r); }); return o; }

  function viewStructure(root) {
    section(root, 'Structure', 'What each transcriptome is made of. Shares are of each run’s own isoforms, so runs of different sizes compare directly.');
    var runs = sel(), g = grid(root);
    var cats = CATS.filter(function (c) { return runs.some(function (r) { return B(r).categories[c] > 0; }); });
    var b1 = card(g, { title: 'Isoforms by structural category', dek: 'Classes against the reference annotation: FSM matches every splice junction of a reference transcript; ISM matches part of one; NIC uses known junctions in a new combination; NNC has at least one new junction. Fusion isoforms join two genes; the rest (genic, antisense, intergenic) match no reference splice junction.' + (cats.length < CATS.length ? ' Classes no run has are left out.' : '') });
    hbars(b1, { cats: cats.map(function (c) { return CAT_LABEL[c]; }), vals: byRun(runs, function (r) { var x = B(r); return cats.map(function (c) { return pc(x.categories[c] || 0, x.summary.isoforms); }); }),
      fmt: function (v) { return pct(v); }, labelW: 96,
      tip: function (r, i) { var x = B(r); return N(x.categories[cats[i]] || 0) + ' of ' + N(x.summary.isoforms) + ' (' + pct(pc(x.categories[cats[i]] || 0, x.summary.isoforms)) + ')'; } });
    var b2 = card(g, { title: 'How many isoforms', dek: 'Total, multi-exon and mono-exon isoforms.' });
    hbars(b2, { cats: ['Isoforms', 'Multi-exon', 'Mono-exon'], vals: byRun(runs, function (r) { var x = B(r).summary; return [x.isoforms, x.multi_exon, x.mono_exon]; }), fmt: N, labelW: 100 });
    [['transcript_length', 'Isoform length', 'Isoform length (nt)', false], ['exons_per_isoform', 'Exons per isoform', 'Exons', true]].forEach(function (h) {
      panel(g, { title: h[1] }, function (r) { return histOf(r, h[0]); }, NEEDS.hist, function (body, ok) {
        var e = histOf(ok[0], h[0]);
        var nm = binNames(e.edges, e.counts.length, h[3]);
        vhist(body, { runs: ok, bins: nm, ticks: h[3] ? nm : e.edges.slice(0, e.counts.length).map(function (x, i) { return i === e.edges.length - 1 ? '≥' + kfmt(x) : kfmt(x); }), counts: byRun(ok, function (r) { return histOf(r, h[0]).counts; }), xlab: h[2], xname: h[2] });
      });
    });
    structureQc(g);
  }

  // Section 13 (A44): SQANTI3-style checks on the Structure tab.
  var SUB_LABEL = {
    '3prime_fragment': '3′ fragment (missing 5′ exons)', '5prime_fragment': '5′ fragment (missing 3′ exons)',
    internal_fragment: 'Internal fragment', intron_retention: 'Intron retention', 'mono-exon': 'Mono-exon',
    combination_of_known_junctions: 'Known junctions, new combination', combination_of_known_splicesites: 'Known splice sites, new junction',
    'mono-exon_by_intron_retention': 'Mono-exon, retains an intron', at_least_one_novel_splicesite: 'At least one new splice site'
  };
  function structureQc(g) {
    var subs = [];
    ['ISM', 'NIC', 'NNC'].forEach(function (c) {
      var ks = {};
      sel().forEach(function (r) { var q = sqc(r, 'subcategories'); if (q && q[c]) Object.keys(q[c]).forEach(function (k) { if (q[c][k] > 0) ks[k] = 1; }); });
      Object.keys(ks).forEach(function (k) { subs.push([c, k]); });
    });
    panel(g, { title: 'ISM, NIC and NNC subcategories', wide: true, dek: 'SQANTI3 subcategories, as shares of each class’s isoforms. Many ISM 3′ fragments point to 5′-truncated reads. Intron retention: an exon spans a reference intron.' },
      function (r) { return sqc(r, 'subcategories'); }, NEEDS.structureQc, function (body, ok) {
        if (!subs.length) { note(body, 'No run has ISM, NIC or NNC isoforms.'); return; }
        hbars(body, { runs: ok, cats: subs.map(function (x) { return x[0] + ' · ' + (SUB_LABEL[x[1]] || x[1]); }), labelW: 250, fmt: function (v) { return pct(v); },
          vals: byRun(ok, function (r) { var q = sqc(r, 'subcategories'), cs = B(r).categories; return subs.map(function (x) { return pc(q[x[0]][x[1]], cs[x[0]]); }); }),
          tip: function (r, i) { var x = subs[i], n = sqc(r, 'subcategories')[x[0]][x[1]], d = B(r).categories[x[0]]; return N(n) + ' of ' + N(d) + ' ' + x[0] + ' (' + pct(pc(n, d)) + ')'; } });
      });
    var rows = [
      ['Multi-exon isoforms, only canonical junctions', function (r) { var c = sqc(r, 'canonical'); return c ? [c.all_canonical, c.multi_exon, 'multi-exon isoforms'] : null; }],
      ['Junctions with an RT-switching repeat', function (r) { var t = sqc(r, 'rts'); return t ? [t.rts, t.junctions, 'junctions'] : null; }],
      ['Isoforms with an RT-switching junction', function (r) { var t = sqc(r, 'rts'); return t ? [t.isoforms_with_rts, t.isoforms, 'multi-exon isoforms'] : null; }]
    ];
    // A45: short-read support. Support in the sample's junctions file (A43) stays in the bundle
    // but is not shown: FLAIR is given that file (A41), so it is not independent.
    var srRows = [
      ['Junctions supported by short reads', function (r) { var j = srsup(r); return j ? [j.supported, j.junctions, 'junctions'] : null; }],
      ['Known junctions supported by short reads', function (r) { var j = srsup(r); return j ? [j.known.supported, j.known.n, 'known junctions'] : null; }],
      ['Novel junctions supported by short reads', function (r) { var j = srsup(r); return j ? [j.novel.supported, j.novel.n, 'novel junctions'] : null; }],
      ['Multi-exon isoforms, every junction in short reads', function (r) { var j = srsup(r); return j ? [j.all_supported, j.multi_exon, 'multi-exon isoforms'] : null; }]
    ];
    var anySr = sel().some(function (r) { return srsup(r); });
    if (anySr) rows = rows.concat(srRows);
    panel(g, { title: 'Junction quality', wide: true, dek: 'Canonical: GT-AG, GC-AG or AT-AC. RT switching: the two sides of the junction share an 8-nt repeat (one mismatch allowed), so it may be a reverse-transcription artefact. Short reads: the junction is in the sample’s sr_junctions file (samples.csv) with at least one read; that file is used only for scoring.' },
      function (r) { return sqc(r, 'canonical') || srsup(r); }, NEEDS.structureQc, function (body, ok) {
        hbars(body, { runs: ok, cats: rows.map(function (x) { return x[0]; }), max: 100, labelW: 270, fmt: function (v) { return pct(v); },
          vals: byRun(ok, function (r) { return rows.map(function (x) { var v = x[1](r); return v ? pc(v[0], v[1]) : null; }); }),
          tip: function (r, i) { var v = rows[i][1](r); return v ? N(v[0]) + ' of ' + N(v[1]) + ' ' + v[2] + ' (' + pct(pc(v[0], v[1])) + ')' : 'not measured'; } });
        if (anySr) {
          var sr0 = ok.filter(srsup)[0], sf = srsup(sr0);
          var dep = ok.filter(function (r) { return srsup(r) && !srIndependent(r); });
          note(body, 'Short-read junctions: ' + (sf.file || 'the sample’s sr_junctions file') + ', ' + N(sf.file_junctions) + ' junctions in the file (sr_junctions column of samples.csv). A junction counts as supported with at least one read (STAR SJ.out.tab: uniquely mapped reads). ' +
            (dep.length ? 'FLAIR was also given this file in ' + dep.map(function (r) { return lab(r); }).join(', ') + ', so support there is partly built in.' : 'FLAIR was not given this file, so this is independent support.'));
        } else note(body, 'Short-read support: hidden, ' + NEEDS.srJunctionFile + '.');
      });
    var cats = CATS.filter(function (c) { return sel().some(function (r) { return B(r).categories[c] > 0; }); });
    if (cats.indexOf(state.qcCat) < 0) state.qcCat = cats.indexOf('ISM') >= 0 ? 'ISM' : cats[0];
    panel(g, { title: 'Length and exons by category', wide: true, dek: 'Isoform length and exon count for one structural category at a time. Pick the category below.' },
      function (r) { return sqc(r, 'by_category'); }, NEEDS.structureQc, function (body, ok) {
        body.appendChild(el('div', { class: 'no-print', style: 'margin:0 0 8px' }, [segment('Category', cats.map(function (c) { return [c, CAT_LABEL[c]]; }), state.qcCat, function (v) { state.qcCat = v; render(); })]));
        var two = el('div', { class: 'grid sub' });
        body.appendChild(two);
        [['length', 'Isoform length (nt)', false], ['exons', 'Exons', true]].forEach(function (h) {
          var e = sqc(ok[0], 'by_category')[h[0]], n = e.edges.length, nm = binNames(e.edges, n, h[2]);
          var box = el('div', {}, [el('p', { class: 'eyebrow', style: 'margin:0 0 6px', text: CAT_LABEL[state.qcCat] + ' · ' + h[1] })]);
          two.appendChild(box);
          vhist(box, { runs: ok, bins: nm, ticks: h[2] ? nm : e.edges.map(function (x, i) { return i === n - 1 ? '≥' + kfmt(x) : kfmt(x); }), xlab: h[1], xname: h[1],
            counts: byRun(ok, function (r) { return sqc(r, 'by_category')[h[0]].counts[state.qcCat] || e.edges.map(function () { return 0; }); }) });
        });
      });
  }

  // Ends by structural category (Ends tab): do novel classes have less supported ends?
  function endsByCategory(g) {
    panel(g, { title: 'End precision by category', wide: true, dek: 'Share of each category’s isoforms with an end within 50 bp of a known site. Left: 5′ end vs refTSS. Right: 3′ end vs PolyA_DB. Shows whether novel classes have less supported ends.' },
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
    // Sensitivity counts genes and precision counts isoforms (A35), so they are not drawn on one axis.
    panel(root, { title: 'TUSCO', wide: true, dek: 'Genes expected to have exactly one isoform (SQANTI3 TUSCO list). Sensitivity: share of these genes whose isoform is found. Precision: share of this run’s isoforms at these genes that match that isoform. Found: exact match; Partial: a partial match; Wrong isoform: only other isoforms; Missed: no isoform.' },
      tusco, NEEDS.tusco, function (body, ok) {
        tiles(body, ok.map(function (r) {
          var t = tusco(r);
          return [runTag(r, true), pct(asPct(t.sensitivity)), 'sensitivity (genes) · precision ' + pct(asPct(t.precision)) + ' (isoforms)'];
        }));
        var t0 = tusco(ok[0]);
        if (isNum(t0.genes_in_regions) && isNum(t0.genes_total)) note(body, N(t0.genes_in_regions) + ' of ' + N(t0.genes_total) + ' TUSCO genes are inside the evaluated regions.');
        var h2 = el('div', { style: 'margin-top:16px' }, [el('p', { class: 'eyebrow', style: 'margin:0 0 6px', text: 'TUSCO genes by call' })]);
        body.appendChild(h2);
        hbars(h2, { runs: ok, cats: ['Found', 'Partial', 'Wrong isoform', 'Missed'], vals: byRun(ok, function (r) { var t = tusco(r); return [t.TP, t.PTP, t.FP, t.FN]; }), fmt: N, labelW: 104,
          tip: function (r, i) { var t = tusco(r), v = [t.TP, t.PTP, t.FP, t.FN][i]; return N(v) + ' of ' + N(t.genes_in_regions) + ' genes'; } });
      });
    root.appendChild(el('p', { class: 'eyebrow', style: 'margin:18px 0 10px', text: 'Against the reference annotation' }));
    var g = grid(root);
    var tot = B(b).summary.ref_genes_total;
    var b3 = card(g, { title: 'Genes detected', dek: 'Reference genes' + (isNum(tot) ? ' (of ' + N(tot) + ' in the annotation, or in --regions)' : '') + ' with at least one FSM, ISM, NIC, NNC or fusion isoform.' });
    runBars(b3, { runs: runs, val: function (r) { return B(r).summary.ref_genes_detected; }, fmt: N, name: 'Genes detected',
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
    var b5 = card(g, { title: 'FSM isoforms with exact ends', dek: 'FSM isoforms, and those with both ends within 50 bp of the matched reference transcript. A recall check: the runs were built with this annotation.' });
    hbars(b5, { cats: ['FSM', 'Both ends within 50 bp'], vals: byRun(runs, function (r) { var x = B(r); return [x.categories.FSM, x.fsm_subcategories.CIM]; }), fmt: N, labelW: 150,
      tip: function (r, i) { var x = B(r); return i ? N(x.fsm_subcategories.CIM) + ' of ' + N(x.categories.FSM) + ' FSM (' + pct(pc(x.fsm_subcategories.CIM, x.categories.FSM)) + ')' : N(x.categories.FSM) + ' isoforms'; } });
    panel(g, { title: 'Reference transcripts found, by length', dek: 'Share of the reference transcripts in each length bin that have an FSM isoform. Shows whether long transcripts are recovered as often as short ones.' },
      function (r) { return get(sqc(r, 'ref_recovery'), 'length'); }, NEEDS.structureQc, function (body, ok) {
        var h0 = sqc(ok[0], 'ref_recovery').length, e = h0.edges, n = e.length;
        var ticks = e.map(function (x, i) { return i === n - 1 ? '≥' + kfmt(x) : kfmt(x); });
        var bins = binNames(e, n, false).map(function (nm, i) { return nm + ' nt (' + N(h0.all[i]) + ' reference transcripts)'; });
        vbars(body, { runs: ok, cap: 100, bins: bins, ticks: ticks, xname: 'Length', xlab: 'Reference transcript length (nt)',
          vals: byRun(ok, function (r) { var h = sqc(r, 'ref_recovery').length; return h.fsm.map(function (f, i) { return pc(f, h.all[i]); }); }),
          fmt: function (v) { return pct(v, 0); }, yfmt: tickPct,
          tip: function (r, i) { var h = sqc(r, 'ref_recovery').length; return N(h.fsm[i]) + ' of ' + N(h.all[i]) + ' reference transcripts (' + pct(pc(h.fsm[i], h.all[i])) + ')'; } });
      });
    var ld = XS().lost_deciles, lb = ld && ld.by_base ? ld.by_base[b] : null;
    if (!others().length) off(g, 'Who are the lost genes?', NEEDS.others);
    else if (!lb) off(g, 'Who are the lost genes?', NEEDS.lostDec);
    else panel(g, { title: 'Who are the lost genes?', dek: 'Baseline genes a run misses, by the baseline’s expression decile (gene CPM, low to high). Losses in the lowest deciles are likely coverage dropouts; losses among well-expressed genes are more likely real differences.', runs: others() },
      function (r) { return lb.lost[r]; }, NEEDS.lostDec, function (body, ok) {
        var bins = [], ticks = [];
        for (var i = 0; i < lb.n.length; i++) { ticks.push(String(i + 1)); bins.push('Decile ' + (i + 1) + ' (' + N(lb.n[i]) + ' genes, CPM ≤ ' + (isNum(lb.cpm_max[i]) ? lb.cpm_max[i].toFixed(1) : '–') + ')'); }
        vbars(body, { runs: ok, cap: 100, bins: bins, ticks: ticks, xname: 'Baseline decile', xlab: 'Baseline expression decile (low → high)',
          vals: byRun(ok, function (r) { return lb.lost[r].map(function (n, i) { return lb.n[i] ? 100 * n / lb.n[i] : null; }); }),
          fmt: function (v) { return pct(v, 0); }, yfmt: tickPct,
          tip: function (r, i) { return N(lb.lost[r][i]) + ' of ' + N(lb.n[i]) + ' genes lost (' + pct(pc(lb.lost[r][i], lb.n[i])) + ')'; } });
        note(body, 'Share of each decile’s baseline genes that the run does not detect (' + N(lb.genes) + ' genes, about ' + N(lb.genes / lb.n.length) + ' per decile). CPM from the baseline’s read counts.');
      });
  }

  var PROD_ORDER = ['PRO', 'PTC', 'NST', 'NGO'];
  var PROD_LABEL = { PRO: 'Productive (PRO)', PTC: 'Premature stop, NMD candidate (PTC)', NST: 'No stop codon (NST)', NGO: 'No annotated start (NGO)', NA: 'Not classified (NA)' };
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
        fmt: function (v) { return pct(v); }, labelW: 200,
        tip: function (r, i) { var cc = prodCounts(r); return N(cc[ks[i]] || 0) + ' isoforms (' + pct(prodShare(r, ks[i])) + ')'; } });
    });
    panel(g, { title: 'FSM isoforms per reference transcript', dek: 'Reference transcripts with an FSM, by how many FSM isoforms match them. More than one means redundant isoforms that differ only in their ends.' },
      function (r) { return sqc(r, 'redundancy'); }, NEEDS.structureQc, function (body, ok) {
        var e = sqc(ok[0], 'redundancy').fsm_hist.edges;
        vhist(body, { runs: ok, bins: binNames(e, e.length, true), counts: byRun(ok, function (r) { return sqc(r, 'redundancy').fsm_hist.counts; }), xlab: 'FSM isoforms per reference transcript', xname: 'FSM isoforms' });
        note(body, ok.map(function (r) { var d = sqc(r, 'redundancy'); return lab(r) + ': ' + N(d.fsm_isoforms) + ' FSM on ' + N(d.fsm_refs) + ' transcripts'; }).join(' · ') + '.');
      });
    // genes_multi = ends_only + splice_only + both + neither (A35); older bundles have no `neither`.
    panel(g, { title: 'How genes vary', wide: true, dek: 'Genes with two or more isoforms, by what makes the isoforms differ: only their ends (same splicing, ends more than 50 bp apart), only their splicing, or both.' },
      variation, NEEDS.variation, function (body, ok) {
        var neither = ok.some(function (r) { return variation(r).neither > 0; });
        var parts = [{ label: 'Ends only', seq: 2 }, { label: 'Splicing only', seq: 4 }, { label: 'Both', seq: 6 }];
        if (neither) parts.push({ label: 'Neither (same splicing, ends within 50 bp)', seq: 1 });
        stack100(body, { runs: ok, parts: parts,
          vals: function (r) { var v = variation(r), a = [v.ends_only, v.splice_only, v.both]; return neither ? a.concat([v.neither || 0]) : a; } });
        note(body, ok.map(function (r) { return lab(r) + ': ' + N(variation(r).genes_multi) + ' multi-isoform genes'; }).join(' · ') + '.');
      });
  }

  function viewAgreement(root) {
    var runs = sel(), b = baseId(), c = conc();
    section(root, 'Agreement', 'How many isoforms the runs share. Isoforms match when their intron chains are identical; mono-exon isoforms when they overlap by at least 50% on the same strand.');
    var g = grid(root);
    var inConc = c ? runs.filter(function (r) { return c.runs.indexOf(r) >= 0; }) : [];
    if (inConc.length < 2) off(g, 'Isoforms shared between runs (Jaccard)', runs.length < 2 ? NEEDS.others : NEEDS.conc);
    else {
      var b1 = card(g, { title: 'Isoforms shared between runs (Jaccard)', dek: 'Shared isoforms ÷ all isoforms of the two runs, for every pair: 0 = nothing shared, 1 = identical.' });
      var viz = el('div', { class: 'viz' });
      b1.appendChild(viz);
      var W = width(viz), n = inConc.length, lw = Math.min(120, W * 0.3), cs = Math.min(92, (W - lw - 4) / n), H = 38 + cs * n;
      var svg = svgBox(viz, W, H, 'Jaccard matrix');
      var chars = Math.max(3, Math.floor(cs / 7.4));
      inConc.forEach(function (r, i) {
        // a colour swatch above each column too, so truncated run names stay tellable apart
        s('rect', { x: lw + cs * i + cs / 2 - 5, y: 0, width: 10, height: 10, rx: 2, style: fill(col(r)), 'data-tip': lab(r) }, svg);
        s('text', { x: lw + cs * i + cs / 2, y: 26, 'text-anchor': 'middle', class: 'lab', text: trunc(lab(r), chars), 'data-tip': lab(r) }, svg);
        s('rect', { x: lw - 14, y: 38 + cs * i + cs / 2 - 5, width: 10, height: 10, rx: 2, style: fill(col(r)) }, svg);
        s('text', { x: lw - 20, y: 38 + cs * i + cs / 2, 'text-anchor': 'end', 'dominant-baseline': 'central', class: 'lab', text: trunc(lab(r), Math.floor((lw - 20) / 6.5)) }, svg);
        inConc.forEach(function (q, j) {
          var v = i === j ? 1 : (pairIn(c, r, q) || {}).jaccard;
          if (!isNum(v)) return;
          var k = seqIdx(v);
          s('rect', { x: lw + cs * j + 1, y: 38 + cs * i + 1, width: cs - 2, height: cs - 2, rx: 3, style: fill('var(--seq-' + k + ')'), 'data-tip': lab(r) + ' × ' + lab(q) + ': ' + f3(v) }, svg);
          if (cs > 26) s('text', { x: lw + cs * j + cs / 2, y: 38 + cs * i + cs / 2, 'text-anchor': 'middle', 'dominant-baseline': 'central', class: 'inbar', style: fill('var(--seq-on-' + k + ')'), text: f2(v) }, svg);
        });
      });
      addTable(b1, [''].concat(inConc.map(lab)), inConc.map(function (r) { return [lab(r)].concat(inConc.map(function (q) { return r === q ? '1.00' : f2((pairIn(c, r, q) || {}).jaccard); })); }));
    }
    if (!others().length) off(g, 'Shared and unique isoforms vs baseline', NEEDS.others);
    else panel(g, { title: 'Shared and unique isoforms vs baseline', dek: 'Isoforms only in the baseline, shared, and only in the run.', runs: others() }, concPair, NEEDS.conc, function (body, ok) {
      sharedBars(body, ok.map(function (r) { return { run: r, p: concPair(r) }; }), 'isoforms');
      if (c && c.runs.length > 2) note(body, 'All ' + c.runs.length + ' runs share ' + N(c.shared_by_all) + ' of ' + N(c.union) + ' isoforms (' + pct(pc(c.shared_by_all, c.union)) + ').');
    });
    if (!others().length) off(root, 'Expression on shared isoforms', NEEDS.others);
    else panel(root, { title: 'Expression on shared isoforms', wide: true, runs: others(), dek: 'Expression of isoforms both runs report: baseline on x, this run on y (log10 CPM+1). Darker squares hold more isoforms; the dashed line is equal expression.' },
      quantPair, NEEDS.cpm, function (body, ok) {
        quantNote(body, [b].concat(ok));
        var lm = XS().quant.log_max || 5;
        var row = el('div', { class: 'multiples dens viz' });
        body.appendChild(row);
        ok.forEach(function (r) {
          var p = quantPair(r);
          var m = el('div', { class: 'mult' }, [el('div', { class: 'mult-title' }, [runTag(r), el('span', { class: 'small', text: '  r = ' + f3(p.pearson) + ' · ' + N(p.n) + ' isoforms' })])]);
          row.appendChild(m);
          var W = Math.min(300, width(m)), S = W - 52, nb = p.grid.length, cw = S / nb, svg = svgBox(m, W, S + 40, 'CPM density');
          var mx = 0;
          p.grid.forEach(function (rw) { rw.forEach(function (v) { if (v > mx) mx = v; }); });
          s('rect', { x: 34, y: 2, width: S, height: S, style: fill('var(--surface-2)') }, svg);
          p.grid.forEach(function (rw, yi) {
            rw.forEach(function (v, xi) {
              if (!v) return;
              var k = Math.max(1, seqIdx(Math.log(v + 1) / Math.log(mx + 1)));
              s('rect', { x: 34 + xi * cw, y: 2 + S - (yi + 1) * cw, width: cw + 0.3, height: cw + 0.3, style: fill('var(--seq-' + k + ')'), 'data-tip': N(v) + ' isoforms' }, svg);
            });
          });
          s('line', { x1: 34, y1: 2 + S, x2: 34 + S, y2: 2, class: 'refl' }, svg);
          var nt = Math.max(1, Math.floor(lm));
          for (var t = 0; t <= nt; t++) {
            s('text', { x: 34 + S * t / lm, y: S + 16, 'text-anchor': 'middle', class: 'tick', text: t }, svg);
            s('text', { x: 28, y: 2 + S - S * t / lm, 'text-anchor': 'end', 'dominant-baseline': 'central', class: 'tick', text: t }, svg);
          }
          s('text', { x: 34 + S / 2, y: S + 34, 'text-anchor': 'middle', class: 'axlab', text: trunc(lab(b), 24) + ' log10(CPM+1)' }, svg);
          s('text', { x: 34 + S + 10, y: 2 + S / 2, 'text-anchor': 'middle', class: 'axlab', transform: 'rotate(90 ' + (34 + S + 10) + ' ' + (2 + S / 2) + ')', text: trunc(lab(r), 24) + ' log10(CPM+1)' }, svg);
        });
        addTable(body, ['Run', 'Shared isoforms', 'Pearson r', 'Spearman ρ', 'Baseline reads on shared isoforms (%)', 'Run reads on shared isoforms (%)'], ok.map(function (r) {
          var p = quantPair(r);
          return [lab(r), N(p.n), f3(p.pearson), f3(p.spearman), pct(100 * p.cpm_shared_a), pct(100 * p.cpm_shared_b)];
        }));
      });
    var g2 = grid(root);
    if (!others().length) off(g2, '5′ and 3′ end agreement', NEEDS.others);
    else panel(g2, { title: '5′ and 3′ end agreement', runs: others(), dek: 'Jaccard of 5′ and 3′ end clusters against the baseline (same strand, ends within ' + ((XS().end_agreement || {}).window || 50) + ' bp grouped), next to the isoform Jaccard.' },
      endPair, NEEDS.endsAgree, function (body, ok) {
        hbars(body, { runs: ok, cats: ['Isoforms', '5′ ends', '3′ ends'], max: 1, labelW: 104, fmt: f2,
          vals: byRun(ok, function (r) { var p = endPair(r), q = concPair(r); return [q ? q.jaccard : null, p.tss.jaccard, p.tes.jaccard]; }),
          tip: function (r, i) {
            var p = i === 0 ? concPair(r) : endPair(r)[i === 1 ? 'tss' : 'tes'];
            return f3(p.jaccard) + ' · ' + N(p.shared) + ' shared, ' + N(p.only_a) + ' only baseline, ' + N(p.only_b) + ' only run';
          } });
      });
    if (!others().length) off(g2, 'TUSCO genes vs baseline', NEEDS.others);
    else panel(g2, { title: 'TUSCO genes vs baseline', runs: others(), dek: 'TUSCO genes this run finds (exact or partial) that the baseline misses, genes both find, and genes only the baseline finds.' },
      tuscoPair, NEEDS.tuscoInc, function (body, ok) {
        hbars(body, { runs: ok, cats: ['Found by this run only', 'Found by both', 'Found by the baseline only'], fmt: N, labelW: 170,
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
  // Truth sets and end databases scored by ted_end_precision.py (A40). One table (rows = runs,
  // one column per set), a precision / recall / F1 panel per truth set, and the both-ends panel.
  function endAccuracy(root, runs, ok) {
    var b = baseId(), sets = eaSets(ok), pairs = eaPairs(ok);
    sets.forEach(function (x) {
      var m = {};
      ok.forEach(function (r) { var f = (endAcc(r).sets || []).find(function (y) { return y.id === x.id; }); if (f) m[(f.file || '') + '|' + (f.md5 || '')] = f; });
      var fs = Object.keys(m).map(function (k) { return m[k]; });
      if (fs.length > 1) root.appendChild(el('p', { class: 'callout' }, [el('strong', { text: 'Not like for like. ' }), 'Runs were scored against different ' + x.assay + ' files: ' + fs.map(function (f) { return f.file + (f.md5 ? ' (md5 ' + f.md5.slice(0, 8) + ')' : ''); }).join(', ') + '.']));
    });
    var top = card(root, { title: 'Precision at 50 bp', wide: true, dek: 'Isoform ends on the same strand within 50 bp of a site, per end set. Databases (refTSS, PolyA_DB): precision. Truth sets: precision (P) and recall (R), with F1.' });
    var head = function (x) {
      return el('th', { class: 'num' }, [endWord(x.end) + ' end vs ' + x.assay, el('span', { class: 'cell-delta', text: x.kind === 'truth' ? (x.provisional ? 'truth, provisional' : 'truth') : 'database' })]);
    };
    var cell = function (r, x) {
      var e = eaRow(r, x.id), w = eaRow(b, x.id);
      if (!e || !isNum(e.precision)) return el('td', { class: 'num muted', text: '–' });
      var p = asP(e.precision), d = [];
      if (r !== b && w && isNum(w.precision)) d.push('P ' + mdelta({ k: 'p' }, p, asP(w.precision)));
      if (x.kind !== 'truth') {
        return el('td', { class: 'num' }, [pct(p) + ' ', el('span', { class: 'muted', text: '(of ' + N(e.isoform_ends) + ' ends)' }), d.length ? el('span', { class: 'cell-delta', text: d[0] + ' vs baseline' }) : null]);
      }
      var rc = asP(e.recall);
      if (r !== b && w && isNum(w.recall) && isNum(rc)) d.push('R ' + mdelta({ k: 'p' }, rc, asP(w.recall)));
      return el('td', { class: 'num' }, ['P ' + pct(p) + ' · R ' + pct(rc),
        el('span', { class: 'cell-delta', text: 'F1 ' + f2(e.f1) + ' · ' + N(e.sites_hit) + ' of ' + N(e.sites) + ' sites' }),
        d.length ? el('span', { class: 'cell-delta', text: d.join(', ') + ' vs baseline' }) : null]);
    };
    top.appendChild(el('div', { class: 'tablewrap' }, [el('table', { class: 'data' }, [
      el('thead', {}, [el('tr', {}, [el('th', { text: 'Run' })].concat(sets.map(head)))]),
      el('tbody', {}, runs.map(function (r) { return el('tr', {}, [el('td', {}, [runTag(r)])].concat(sets.map(function (x) { return cell(r, x); }))); }))])]));
    missingNote(top, runs, ok, NEEDS.truth);
    var g = grid(root);
    sets.filter(function (x) { return x.kind === 'truth'; }).forEach(function (x) {
      var w = endWord(x.end);
      panel(g, { title: w + ' end vs ' + x.assay + (x.provisional ? ' (provisional)' : ''),
        dek: 'Precision: isoform ' + w + ' ends within 50 bp of a ' + x.assay + ' site. Recall: ' + x.assay + ' sites with an isoform ' + w + ' end within 50 bp. F1: their harmonic mean.' },
        function (r) { var e = eaRow(r, x.id); return e && isNum(e.precision) ? e : null; }, NEEDS.truth, function (body, rs) {
          body.appendChild(el('p', { class: 'eyebrow', style: 'margin:0 0 6px', text: 'Precision and recall' }));
          hbars(body, { runs: rs, cats: ['Precision', 'Recall'], max: 100, labelW: 90, fmt: function (v) { return pct(v); },
            vals: byRun(rs, function (r) { var e = eaRow(r, x.id); return [asP(e.precision), asP(e.recall)]; }),
            tip: function (r, i) {
              var e = eaRow(r, x.id);
              return i === 0 ? pct(asP(e.precision)) + ' of ' + N(e.isoform_ends) + ' isoform ' + w + ' ends'
                : pct(asP(e.recall)) + ' (' + N(e.sites_hit) + ' of ' + N(e.sites) + ' ' + x.assay + ' sites)';
            } });
          body.appendChild(el('p', { class: 'eyebrow', style: 'margin:14px 0 6px', text: 'F1' }));
          runBars(body, { runs: rs, max: 1, fmt: f2, name: 'F1', labelW: 120, val: function (r) { return eaRow(r, x.id).f1; },
            tip: function (r) { return f3(eaRow(r, x.id).f1); } });
        });
    });
    if (pairs.length) {
      var pname = function (p) { return p.tss_set + ' × ' + p.tes_set; };
      panel(g, { title: 'Both ends on a site', dek: 'Isoforms whose 5′ and 3′ ends both lie within 50 bp of a site, out of all isoforms. Within a junction chain, isoforms that land on the same (5′ site, 3′ site) pair count once.' },
        function (r) { return pairs.some(function (p) { var q = eaPair(r, p.tss_set, p.tes_set); return q && isNum(q.precision); }) ? 1 : null; }, NEEDS.truth, function (body, rs) {
          hbars(body, { runs: rs, cats: pairs.map(function (p) { return pname(p) + (p.kind === 'truth' ? '' : ' (databases)'); }), labelW: 170, fmt: function (v) { return pct(v); },
            vals: byRun(rs, function (r) { return pairs.map(function (p) { var q = eaPair(r, p.tss_set, p.tes_set); return q ? asP(q.precision) : null; }); }),
            tip: function (r, i) {
              var q = eaPair(r, pairs[i].tss_set, pairs[i].tes_set);
              return pct(asP(q.precision)) + ' (' + N(q.both_hit) + ' of ' + N(q.isoforms) + ' isoforms; ' + N(q.both_hit_raw) + ' before merging same-pair isoforms)';
            } });
        });
    }
    var sc = null;
    ok.forEach(function (r) { if (!sc) sc = endAcc(r).scorer; });
    var files = sets.map(function (x) { return x.assay + ': ' + (x.file || '–') + (x.md5 ? ', md5 ' + x.md5.slice(0, 8) : ''); });
    note(root, 'Scored by ' + (sc && sc.script ? sc.script : 'ted_end_precision.py') + (sc && sc.md5 ? ' (md5 ' + sc.md5.slice(0, 8) + ')' : '') +
      ', window ' + (sc && isNum(sc.window) ? sc.window : 50) + ' bp. End sets: ' + files.join(' · ') + '.');
  }
  function viewEnds(root) {
    var runs = sel(), b = baseId();
    section(root, 'Ends', 'How close each isoform’s 5′ and 3′ ends sit to known sites (refTSS for transcription start sites, PolyA_DB for poly(A) cleavage sites), and, with --truth, to the sample’s own truth sets (CAGE, PRO-cap, direct RNA …), and whether 3′ ends sit on genomic A-runs. Distances are signed in transcript orientation: negative means the isoform end is upstream (5′) of the site.');
    var tssRuns = withData(runs, function (r) { return endsHalf(r, 'tss'); });
    var tesRuns = withData(runs, function (r) { return endsHalf(r, 'tes'); });
    var eaRuns = withData(runs, endAcc);
    if (!tssRuns.length && !tesRuns.length && !eaRuns.length) {
      off(root, 'Ends against refTSS and PolyA_DB', NEEDS.ends);
    } else {
      var truthSets = eaSets(eaRuns).filter(function (x) { return x.kind === 'truth'; });
      root.appendChild(el('p', { class: 'callout' }, truthSets.length ? [el('strong', { text: 'Databases: precision only. ' }),
        'A refTSS or PolyA_DB site near an end supports it, but the databases pool many tissues, so an end with no site nearby is not necessarily wrong and recall is not scored against them. ',
        el('strong', { text: 'Truth sets: precision and recall. ' }),
        'Truth sets (' + truthSets.map(function (x) { return x.assay; }).join(', ') + ') are this sample’s own end calls from an orthogonal assay, so recall is scored too: the share of truth sites with an isoform end within 50 bp. Every truth site counts, also those of genes the run did not sample.']
        : [el('strong', { text: 'Precision only. ' }),
          'A known site near an end supports it. The databases pool many tissues, so an end with no site nearby is not necessarily wrong, and recall is not scored against them.']));
      var prov = truthSets.filter(function (x) { return x.provisional; });
      if (prov.length) root.appendChild(el('p', { class: 'callout' }, [el('strong', { text: 'Provisional truth. ' }),
        prov.map(function (x) { return x.assay; }).join(', ') + (prov.length > 1 ? ' are' : ' is') + ' marked provisional in the truth sheet: the calls have not passed a final quality gate yet. Read recall and F1 against ' + (prov.length > 1 ? 'them' : 'it') + ' as indicative, and compare runs rather than quote the numbers.']));
      [['tss', 'refTSS'], ['pas', 'PolyA_DB']].forEach(function (k) {
        var refs = refLine(runs, k[0]);
        if (refs.length > 1) root.appendChild(el('p', { class: 'callout' }, [el('strong', { text: 'Not like for like. ' }), 'Runs were scored against different ' + k[1] + ' files: ' + refs.map(function (f) { return f.file + (f.md5 ? ' (md5 ' + f.md5.slice(0, 8) + ')' : ''); }).join(', ') + '.']));
      });
      if (eaRuns.length) endAccuracy(root, runs, eaRuns);
      else {
        // One small table instead of a tile per run and end: the same numbers as the dashed 50 bp
        // point of the two line charts below, with counts and the change against the baseline.
        var top = card(root, { title: 'Precision at 50 bp', wide: true, dek: 'Share of isoforms whose 5′ end is within 50 bp of a refTSS site, and whose 3′ end is within 50 bp of a PolyA_DB site.' });
        var cell = function (r, k) {
          var h = endsHalf(r, k), v = within(r, k, 50), w = within(b, k, 50);
          if (!h || !isNum(v)) return el('td', { class: 'num muted', text: '–' });
          var d = r !== b && isNum(w) ? mdelta({ k: 'p' }, v, w) + ' vs baseline' : '';
          return el('td', { class: 'num' }, [pct(v) + ' ', el('span', { class: 'muted', text: '(' + N(h.within[windowsOf(r).indexOf(50)]) + ' of ' + N(h.n) + ')' }), d ? el('span', { class: 'cell-delta', text: d }) : null]);
        };
        top.appendChild(el('div', { class: 'tablewrap' }, [el('table', { class: 'data' }, [
          el('thead', {}, [el('tr', {}, [el('th', { text: 'Run' }), el('th', { class: 'num', text: '5′ end vs refTSS' }), el('th', { class: 'num', text: '3′ end vs PolyA_DB' })])]),
          el('tbody', {}, runs.map(function (r) { return el('tr', {}, [el('td', {}, [runTag(r)]), cell(r, 'tss'), cell(r, 'tes')]); }))])]));
        if (tssRuns.length || tesRuns.length) off(root, 'Ends against truth sets', NEEDS.truth);
      }
      if (eaRuns.length && (tssRuns.length || tesRuns.length))
        note(root, 'The charts below measure each end’s distance to the nearest refTSS or PolyA_DB site position, so their 50 bp values can differ a little from the table above, which counts an end within 50 bp of any part of a site’s interval.');
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
      endsByCategory(g);
      [['tss', 'Distance to the nearest refTSS site', '5′ end to nearest refTSS site (bp)'], ['tes', 'Distance to the nearest PolyA_DB site', '3′ end to nearest PolyA_DB site (bp)']].forEach(function (k) {
        panel(g, { title: k[1], dek: 'Signed distance to the nearest same-strand site; negative = upstream. Isoforms with no site within 10 kb are left out.' }, function (r) { var h = endsHalf(r, k[0]); return h && h.dist_hist ? h : null; }, NEEDS[k[0]], function (body, ok) {
          var E = endsHalf(ok[0], k[0]).dist_hist.edges, db = distBins(E);
          vhist(body, { runs: ok, bins: db.full, ticks: db.ticks, mark: db.zero, xname: 'Distance', counts: byRun(ok, function (r) { return endsHalf(r, k[0]).dist_hist.counts; }), xlab: k[2] });
        });
      });
      [['fsm_tss_hist', 'FSM 5′ end vs the reference transcript', 'Distance to the reference 5′ end (bp)'], ['fsm_tes_hist', 'FSM 3′ end vs the reference transcript', 'Distance to the reference 3′ end (bp)']].forEach(function (k) {
        panel(g, { title: k[1], dek: 'Signed distance from each FSM isoform’s end to its matched reference transcript’s end.' },
          function (r) { var h = get(B(r), 'sections.ends.' + k[0]); return h && h.counts && sum(h.counts) > 0 ? h : null; }, NEEDS.ends, function (body, ok) {
            var E = get(B(ok[0]), 'sections.ends.' + k[0]).edges, db = distBins(E);
            vhist(body, { runs: ok, bins: db.full, ticks: db.ticks, mark: db.zero, xname: 'Distance', counts: byRun(ok, function (r) { return get(B(r), 'sections.ends.' + k[0]).counts; }), xlab: k[2] });
          });
      });
      panel(g, { title: 'PolyA signal before the 3′ end', dek: 'Isoforms with a polyA hexamer (AATAAA, ATTAAA or one of ten variants) in the 50 nt upstream of the 3′ end.' },
        function (r) { var m = get(B(r), 'sections.ends.polya_motif'); return m && isNum(m.n) ? m : null; }, NEEDS.motif, function (body, ok) {
          runBars(body, { runs: ok, max: 100, fmt: function (v) { return pct(v); }, name: 'With a polyA motif',
            val: function (r) { var m = get(B(r), 'sections.ends.polya_motif'); return pc(m.n_with_motif, m.n); },
            tip: function (r) { var m = get(B(r), 'sections.ends.polya_motif'); return N(m.n_with_motif) + ' of ' + N(m.n) + ' (' + pct(pc(m.n_with_motif, m.n)) + ')'; } });
        });
      var refs = [['tss', 'refTSS'], ['pas', 'PolyA_DB']].map(function (k) {
        return refLine(runs, k[0]).map(function (f) { return k[1] + ': ' + f.file + (isNum(f.sites) ? ' (' + N(f.sites) + ' sites)' : '') + (f.md5 ? ', md5 ' + f.md5.slice(0, 8) : ''); }).join('; ');
      }).filter(function (x) { return x; });
      if (refs.length) note(root, 'References: ' + refs.join(' · ') + '.');
    }
    // Internal priming needs only the genome, not the end databases, so it is drawn either way.
    root.appendChild(el('p', { class: 'eyebrow', style: 'margin:18px 0 10px', text: '3′ ends on genomic A-runs' }));
    var g2 = grid(root);
    panel(g2, { title: 'Genomic A content after the 3′ end', dek: 'Share of A in the 20 nt of genome just downstream of each 3′ end. High values point to internal priming on a genomic A-run. Dashed: 60%.' },
      function (r) { var p = priming(r); return p && p.down_a_hist ? p : null; }, NEEDS.priming, function (body, ok) {
        var h = priming(ok[0]).down_a_hist, n = h.counts.length, e = h.edges, bins = [], ticks = [];
        for (var i = 0; i < n; i++) {
          var a = Math.round(100 * e[i]), z = e[i + 1];
          if (z === undefined) { bins.push(a + '% A'); ticks.push(String(a)); } else { bins.push(a + '–' + Math.round(100 * z) + '% A'); ticks.push(String(a)); }
        }
        var mark = null;
        e.forEach(function (x, i) { if (Math.abs(x - 0.6) < 1e-9) mark = i; });
        vhist(body, { runs: ok, bins: bins, ticks: ticks, mark: mark, xname: 'A downstream', counts: byRun(ok, function (r) { return priming(r).down_a_hist.counts; }), xlab: 'A in the 20 nt after the 3′ end (%)' });
      });
    panel(g2, { title: 'Internal-priming risk', dek: 'Isoforms with at least 60% A in the 20 nt after the 3′ end.' }, priming, NEEDS.priming, function (body, ok) {
      runBars(body, { runs: ok, fmt: function (v) { return pct(v); }, name: '≥ 60% A downstream',
        val: function (r) { var p = priming(r); return pc(p.risk, p.n); },
        tip: function (r) { var p = priming(r); return N(p.risk) + ' of ' + N(p.n) + ' isoforms (' + pct(pc(p.risk, p.n)) + ')'; } });
    });
    panel(g2, { title: 'Internal-priming risk by group', wide: true, dek: 'The same risk (at least 60% A after the 3′ end) for mono- and multi-exon isoforms, each structural category, and coding vs non-coding isoforms when productivity ran (coding: PRO, PTC, NST; non-coding: NGO).' },
      function (r) { return sqc(r, 'priming'); }, NEEDS.structureQc, function (body, ok) {
        var groups = [['by_exons', 'mono', 'Mono-exon'], ['by_exons', 'multi', 'Multi-exon']];
        CATS.forEach(function (c) { if (ok.some(function (r) { var x = sqc(r, 'priming').by_category[c]; return x && x.n; })) groups.push(['by_category', c, CAT_LABEL[c]]); });
        if (ok.some(function (r) { return sqc(r, 'priming').by_coding; })) groups.push(['by_coding', 'coding', 'Coding'], ['by_coding', 'noncoding', 'Non-coding']);
        var at = function (r, i) { var x = sqc(r, 'priming')[groups[i][0]]; return x ? x[groups[i][1]] : null; };
        hbars(body, { runs: ok, cats: groups.map(function (x) { return x[2]; }), labelW: 110, fmt: function (v) { return pct(v); },
          vals: byRun(ok, function (r) { return groups.map(function (_x, i) { var v = at(r, i); return v && v.n ? pc(v.risk, v.n) : null; }); }),
          tip: function (r, i) { var v = at(r, i); return v ? N(v.risk) + ' of ' + N(v.n) + ' isoforms (' + pct(pc(v.risk, v.n)) + ')' : 'not measured'; } });
      });
  }

  function viewReads(root) {
    var runs = sel();
    section(root, 'Quantification', 'How many reads each run assigns to its isoforms, and how expression is spread.');
    quantNote(root, withData(runs, quant));
    var g = grid(root);
    panel(g, { title: 'Reads assigned', dek: 'Reads FLAIR quantify assigns to the run’s isoforms, and isoforms with at least one read.' }, quant, NEEDS.quant, function (body, ok) {
      body.appendChild(el('p', { class: 'eyebrow', style: 'margin:0 0 6px', text: 'Reads assigned' }));
      runBars(body, { runs: ok, fmt: N, name: 'Reads assigned', val: function (r) { return quant(r).reads_assigned; } });
      body.appendChild(el('p', { class: 'eyebrow', style: 'margin:14px 0 6px', text: 'Isoforms with reads' }));
      runBars(body, { runs: ok, fmt: N, name: 'Isoforms with reads', val: function (r) { return quant(r).isoforms_quantified; } });
    });
    panel(g, { title: 'Isoform expression', dek: 'Isoforms with reads, by CPM.' }, function (r) { var q = quant(r); return q && q.cpm_hist && q.cpm_hist.counts ? q : null; }, NEEDS.quant, function (body, ok) {
      var h = quant(ok[0]).cpm_hist;
      vhist(body, { runs: ok, bins: binNames(h.edges, h.counts.length, false).map(function (x) { return x + ' CPM'; }), ticks: h.edges.slice(0, h.counts.length).map(function (x, i) { return i === h.edges.length - 1 ? '≥' + kfmt(x) : kfmt(x); }), xname: 'CPM', counts: byRun(ok, function (r) { return quant(r).cpm_hist.counts; }), xlab: 'CPM' });
    });
  }

  function viewCompute(root) {
    section(root, 'Compute', 'FLAIR tasks of each run (assemble, quantify, productivity), from the Nextflow trace. The shared align step and the evaluation are not counted.');
    var g = grid(root, 'three');
    [['wall_minutes', 'Wall time', fmtWall], ['cpu_hours', 'CPU time', fmtCpu], ['peak_rss_gb', 'Peak memory', fmtMem]].forEach(function (m) {
      panel(g, { title: m[1] }, function (r) { return compute(r, m[0]); }, NEEDS.compute, function (body, ok) {
        runBars(body, { runs: ok, val: function (r) { return compute(r, m[0]); }, fmt: m[2], name: m[1], labelW: 120 });
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
    { id: 'reads', label: 'Quantification', fn: viewReads },
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
    document.getElementById('subtitle').textContent = RUN_ORDER.length + ' run' + (RUN_ORDER.length === 1 ? '' : 's') + ' · ' + SAMPLES.length + ' sample' + (SAMPLES.length === 1 ? '' : 's') +
      (when ? ' · built ' + when : '') + (DATA.tool_version ? ' · flair-eval ' + DATA.tool_version : '');
    document.getElementById('footer').textContent = 'Each run keeps one colour in every tab; grey means shared. One sample is shown at a time; samples are never pooled.';
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
