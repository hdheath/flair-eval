/* flair-eval report: overview table, single-run page, comparison view.
 * Plain JS + inline SVG; no libraries, no network. Data: <script id="report-data">.
 * All text from the data goes in with textContent, never innerHTML. */
(function () {
  'use strict';

  var DATA = JSON.parse(document.getElementById('report-data').textContent);
  var SVGNS = 'http://www.w3.org/2000/svg';

  var CATS = ['FSM', 'ISM', 'NIC', 'NNC', 'genic', 'genic_intron', 'antisense', 'intergenic', 'fusion'];
  var CAT_LABEL = {
    FSM: 'FSM', ISM: 'ISM', NIC: 'NIC', NNC: 'NNC', genic: 'Genic', genic_intron: 'Genic intron',
    antisense: 'Antisense', intergenic: 'Intergenic', fusion: 'Fusion'
  };
  var CAT_LONG = {
    FSM: 'Full splice match', ISM: 'Incomplete splice match', NIC: 'Novel in catalog',
    NNC: 'Novel not in catalog', genic: 'Genic', genic_intron: 'Genic intron',
    antisense: 'Antisense', intergenic: 'Intergenic', fusion: 'Fusion'
  };
  var SUBS = ['CIM', 'FSM5', 'FSM3', 'FSMTE'];
  var SUB_LABEL = { CIM: 'CIM (both ends match)', FSM5: "FSM5 (5' end new)", FSM3: "FSM3 (3' end new)", FSMTE: 'FSMTE (both ends new)' };
  var HISTS = [
    { key: 'transcript_length', title: 'Transcript length', unit: 'nt', of: 'isoforms', cont: true },
    { key: 'exon_length', title: 'Exon length', unit: 'nt', of: 'exons', cont: true },
    { key: 'exons_per_isoform', title: 'Exons per isoform', unit: 'exons', of: 'isoforms', cont: false },
    { key: 'isoforms_per_gene', title: 'Isoforms per gene', unit: 'isoforms', of: 'genes', cont: false }
  ];
  // Overview category mix: five colours + "Other" in gray (fixed order, never cycled).
  var MIX = [
    { label: 'FSM', cats: ['FSM'], color: 'var(--series-1)', on: 'var(--on-series)' },
    { label: 'ISM', cats: ['ISM'], color: 'var(--series-2)', on: 'var(--on-series)' },
    { label: 'NIC', cats: ['NIC'], color: 'var(--series-3)', on: 'var(--on-series)' },
    { label: 'NNC', cats: ['NNC'], color: 'var(--series-4)', on: 'var(--on-series-4)' },
    { label: 'Genic', cats: ['genic'], color: 'var(--series-5)', on: 'var(--on-series)' },
    { label: 'Other', cats: ['genic_intron', 'antisense', 'intergenic', 'fusion'], color: 'var(--series-other)', on: 'var(--on-series-other)' }
  ];
  var SLOT = function (i) { return 'var(--series-' + (i + 1) + ')'; };
  var ON_SLOT = function (i) { return i === 3 ? 'var(--on-series-4)' : 'var(--on-series)'; };

  // Panels for later report sections (bundle.sections.<key>). v1 has none; unknown keys are
  // ignored and a panel is drawn only when its section is present.
  var SECTION_PANELS = {};

  // ------------------------------------------------------------------ data
  var bundles = DATA.bundles.map(function (b) {
    var o = Object.assign({}, b);
    o.key = b.sample + '/' + b.run_id;
    o.label = b.label || b.run_id;
    return o;
  });
  var runOrder = DATA.runs.slice();
  var samples = DATA.samples.slice();
  var runInfo = {};
  bundles.forEach(function (b) { if (!runInfo[b.run_id]) runInfo[b.run_id] = b; });

  function runColor(runId) {
    var i = runOrder.indexOf(runId);
    return i >= 0 && i < 8 ? SLOT(i) : 'var(--series-other)';
  }

  var state = {
    view: 'overview',
    filters: { sample: '', version: '' },
    tagFilters: {},
    baseline: DATA.baseline,
    runKey: null,
    compareSample: null,
    unticked: {},
    mode: 'share'
  };

  // ------------------------------------------------------------------ helpers
  function el(tag, attrs) {
    var node = document.createElement(tag);
    if (attrs) Object.keys(attrs).forEach(function (k) {
      var v = attrs[k];
      if (v === null || v === undefined || v === false) return;
      if (k === 'text') node.textContent = v;
      else if (k === 'class') node.className = v;
      else if (k === 'style') node.setAttribute('style', v);
      else if (k.slice(0, 2) === 'on') node.addEventListener(k.slice(2), v);
      else node.setAttribute(k, v === true ? '' : v);
    });
    for (var i = 2; i < arguments.length; i++) append(node, arguments[i]);
    return node;
  }
  function append(node, c) {
    if (c === null || c === undefined || c === false) return;
    if (Array.isArray(c)) { c.forEach(function (x) { append(node, x); }); return; }
    node.appendChild(typeof c === 'string' || typeof c === 'number' ? document.createTextNode(String(c)) : c);
  }
  function s(tag, attrs) {
    var node = document.createElementNS(SVGNS, tag);
    if (attrs) Object.keys(attrs).forEach(function (k) {
      var v = attrs[k];
      if (v === null || v === undefined) return;
      if (k === 'text') node.textContent = v;
      else if (k === 'fill' || k === 'stroke') node.style[k] = v;
      else node.setAttribute(k, v);
    });
    for (var i = 2; i < arguments.length; i++) if (arguments[i]) node.appendChild(arguments[i]);
    return node;
  }
  var nf = new Intl.NumberFormat('en-US');
  function fmtInt(v) { return v === null || v === undefined || isNaN(v) ? '—' : nf.format(Math.round(v)); }
  function fmtPct(v, d) { return v === null || v === undefined || isNaN(v) ? '—' : v.toFixed(d === undefined ? 1 : d) + '%'; }
  function fmtNum(v, d) { return v === null || v === undefined || isNaN(v) ? '—' : v.toFixed(d === undefined ? 1 : d); }
  function compact(v) {
    if (Math.abs(v) >= 1e6) return (v / 1e6).toFixed(v % 1e6 ? 1 : 0) + 'M';
    if (Math.abs(v) >= 1e4) return (v / 1e3).toFixed(0) + 'k';
    if (Math.abs(v) >= 1e3 && v % 1000 !== 0) return (v / 1e3).toFixed(1) + 'k';
    if (Math.abs(v) >= 1e3) return (v / 1e3).toFixed(0) + 'k';
    return String(v);
  }
  function pct(n, d) { return d ? (100 * n / d) : null; }
  function sum(a) { return a.reduce(function (x, y) { return x + (y || 0); }, 0); }
  function get(o, path) {
    return path.split('.').reduce(function (x, k) { return x === null || x === undefined ? undefined : x[k]; }, o);
  }
  function niceMax(v) {
    if (!(v > 0)) return 1;
    var p = Math.pow(10, Math.floor(Math.log10(v))), f = v / p;
    return (f <= 1 ? 1 : f <= 2 ? 2 : f <= 2.5 ? 2.5 : f <= 5 ? 5 : 10) * p;
  }
  function ticks(max, n) {
    var step = niceMax(max / (n || 4)), out = [];
    for (var t = 0; t <= max + 1e-9; t += step) out.push(+t.toFixed(10));
    return out;
  }

  // ------------------------------------------------------------------ tooltip
  var tip = document.getElementById('tooltip');
  function showTip(evt, title, rows) {
    tip.textContent = '';
    if (title) tip.appendChild(el('div', { class: 'tt-title', text: title }));
    rows.forEach(function (r) {
      var row = el('div', { class: 'tt-row' });
      if (r.color) row.appendChild(el('span', { class: 'key line', style: 'background:' + r.color }));
      row.appendChild(el('strong', { text: r.value }));
      if (r.label) row.appendChild(el('span', { class: 'lbl', text: r.label }));
      tip.appendChild(row);
    });
    tip.hidden = false;
    var x, y;
    if (evt && evt.clientX !== undefined && evt.type !== 'focus') { x = evt.clientX; y = evt.clientY; }
    else { var r = evt.target.getBoundingClientRect(); x = r.left + r.width / 2; y = r.top; }
    var w = tip.offsetWidth, h = tip.offsetHeight, vw = window.innerWidth;
    var left = Math.min(Math.max(8, x + 14), vw - w - 8);
    var top = y - h - 12 < 8 ? y + 16 : y - h - 12;
    tip.style.left = left + 'px'; tip.style.top = top + 'px';
  }
  function hideTip() { tip.hidden = true; }
  function hover(node, fn) {
    node.setAttribute('tabindex', '0');
    node.addEventListener('pointermove', function (e) { var t = fn(); showTip(e, t[0], t[1]); });
    node.addEventListener('focus', function (e) { var t = fn(); showTip(e, t[0], t[1]); });
    node.addEventListener('pointerleave', hideTip);
    node.addEventListener('blur', hideTip);
  }

  // ------------------------------------------------------------------ building blocks
  function legend(items, kind) {
    return el('ul', { class: 'legend' }, items.map(function (it) {
      return el('li', null, el('span', { class: 'key' + (kind === 'line' ? ' line' : ''), style: 'background:' + it.color }), it.label);
    }));
  }

  // A card with an optional chart/table toggle. table = {head: [...], rows: [[...]]}
  function panel(title, note, body, table, tools) {
    var card = el('section', { class: 'card' });
    var head = el('div', { class: 'card-head' }, el('div', null, el('h3', { text: title }), note ? el('p', { class: 'note', text: note }) : null));
    var toolsBox = el('div', { class: 'panel-tools' }, tools || null);
    head.appendChild(toolsBox);
    card.appendChild(head);
    var bodyBox = el('div', null, body);
    card.appendChild(bodyBox);
    if (table) {
      var tbl = dataTable(table.head, table.rows);
      tbl.hidden = true;
      card.appendChild(tbl);
      var btn = el('button', { class: 'btn', type: 'button', 'aria-pressed': 'false', text: 'Table' });
      btn.addEventListener('click', function () {
        var on = tbl.hidden;
        tbl.hidden = !on; bodyBox.hidden = on;
        btn.setAttribute('aria-pressed', on ? 'true' : 'false');
      });
      toolsBox.appendChild(btn);
    }
    return card;
  }
  function dataTable(head, rows) {
    return el('div', { class: 'table-wrap' }, el('table', null,
      el('thead', null, el('tr', null, head.map(function (h, i) { return el('th', { class: i ? 'num' : null, text: h }); }))),
      el('tbody', null, rows.map(function (r) {
        return el('tr', null, r.map(function (c, i) { return el('td', { class: i ? 'num' : null, text: c }); }));
      }))));
  }

  // Tall segmented bar: one whole split into parts. parts = [{label, value, color, on}]
  function segbar(parts, opts) {
    opts = opts || {};
    var total = sum(parts.map(function (p) { return p.value; }));
    var bar = el('div', { class: 'segbar' + (opts.thin ? ' thin' : ''), role: 'img' });
    bar.setAttribute('aria-label', (opts.name ? opts.name + ': ' : '') + parts.map(function (p) {
      return p.label + ' ' + fmtInt(p.value) + ' (' + fmtPct(pct(p.value, total)) + ')';
    }).join(', '));
    parts.forEach(function (p) {
      if (!p.value) return;
      var share = pct(p.value, total);
      var part = el('div', { class: 'part', style: 'flex:' + p.value + ' 1 0;background:' + p.color + ';color:' + (p.on || 'var(--on-series)') });
      bar.appendChild(part);
      if (!opts.thin) {
        requestAnimationFrame(function () {
          var txt = share >= 9.5 ? Math.round(share) + '%' : share.toFixed(1) + '%';
          if (part.clientWidth > txt.length * 7.5 + 12) part.textContent = txt;
        });
      }
      hover(part, function () {
        return [opts.name || '', [{ color: p.color, value: fmtInt(p.value) + ' (' + fmtPct(share) + ')', label: p.label }]];
      });
    });
    return bar;
  }

  function roundedBar(x, y, w, h, dir) {
    // 4px rounded data end, square at the baseline. dir 'right' (horizontal) or 'up'.
    var r = Math.max(0, Math.min(4, dir === 'right' ? w / 2 : h / 2, dir === 'right' ? h / 2 : w / 2));
    if (dir === 'right') {
      return 'M' + x + ',' + y + 'H' + (x + w - r) + 'Q' + (x + w) + ',' + y + ' ' + (x + w) + ',' + (y + r) +
        'V' + (y + h - r) + 'Q' + (x + w) + ',' + (y + h) + ' ' + (x + w - r) + ',' + (y + h) + 'H' + x + 'Z';
    }
    return 'M' + x + ',' + (y + h) + 'V' + (y + r) + 'Q' + x + ',' + y + ' ' + (x + r) + ',' + y +
      'H' + (x + w - r) + 'Q' + (x + w) + ',' + y + ' ' + (x + w) + ',' + (y + r) + 'V' + (y + h) + 'Z';
  }

  // Grouped horizontal bars. series = [{label, color}], groups = [{label, values: [v per series], extra?}]
  function hbars(container, cfg) {
    var width = Math.max(280, container.clientWidth || 600);
    var ns = cfg.series.length;
    var barH = ns === 1 ? 16 : ns <= 3 ? 11 : ns <= 6 ? 8 : 6;
    var groupH = ns * barH + (ns - 1) * 2;
    var pad = ns === 1 ? 8 : 12;
    var labelW = Math.min(150, Math.max(80, width * 0.28));
    var valueW = cfg.labelValues ? 92 : 12;
    var top = 4, axisH = 22;
    var plotW = width - labelW - valueW;
    var max = cfg.max || niceMax(Math.max.apply(null, cfg.groups.map(function (g) {
      return Math.max.apply(null, g.values.map(function (v) { return v || 0; }));
    })));
    var height = top + cfg.groups.length * (groupH + pad) + axisH;
    var svg = s('svg', { viewBox: '0 0 ' + width + ' ' + height, role: 'img', 'aria-label': cfg.name || 'bar chart' });
    var x = function (v) { return labelW + plotW * (v / max); };
    ticks(max, width < 420 ? 3 : 4).forEach(function (t) {
      svg.appendChild(s('line', { class: t === 0 ? 'baseline' : 'gridline', x1: x(t), x2: x(t), y1: top, y2: height - axisH + 4 }));
      svg.appendChild(s('text', { class: 'tick', x: x(t), y: height - 6, 'text-anchor': 'middle', text: cfg.tickFmt ? cfg.tickFmt(t) : compact(t) }));
    });
    cfg.groups.forEach(function (g, gi) {
      var gy = top + gi * (groupH + pad) + pad / 2;
      svg.appendChild(s('text', { class: 'cat-label', x: labelW - 10, y: gy + groupH / 2 + 4, 'text-anchor': 'end', text: g.label }));
      g.values.forEach(function (v, si) {
        if (v === null || v === undefined) return;
        var y = gy + si * (barH + 2), w = Math.max(0, x(v) - labelW);
        var path = s('path', { class: 'mark', d: roundedBar(labelW, y, Math.max(w, v > 0 ? 1 : 0), barH, 'right'), fill: cfg.series[si].color });
        var hit = s('rect', { class: 'hit', x: labelW, y: y - 1, width: plotW, height: barH + 2 });
        svg.appendChild(path);
        svg.appendChild(hit);
        hover(hit, function () {
          return [g.label, [{ color: cfg.series[si].color, value: cfg.fmt(v, g, si), label: ns > 1 ? cfg.series[si].label : '' }]];
        });
        if (cfg.labelValues) {
          svg.appendChild(s('text', { x: x(v) + 6, y: y + barH / 2 + 4, text: cfg.fmt(v, g, si) }));
        }
      });
    });
    container.appendChild(svg);
  }

  // Columns (one series) or lines (several), over categorical bins.
  // cfg = {bins: [label], binTitles: [label], series: [{label,color}], values: [[...] per series], fmt, yFmt}
  function binChart(container, cfg) {
    var width = Math.max(260, container.clientWidth || 560);
    var height = width < 420 ? 190 : 220;
    var left = 44, right = 8, top = 10, bottom = 26;
    var plotW = width - left - right, plotH = height - top - bottom;
    var nb = cfg.bins.length, ns = cfg.series.length;
    var all = [].concat.apply([], cfg.values).filter(function (v) { return v !== null && v !== undefined; });
    var max = niceMax(Math.max.apply(null, all.concat([0])));
    var band = plotW / nb;
    var y = function (v) { return top + plotH * (1 - v / max); };
    var cx = function (i) { return left + band * (i + 0.5); };
    var svg = s('svg', { viewBox: '0 0 ' + width + ' ' + height, role: 'img', 'aria-label': cfg.name || 'histogram' });
    ticks(max, 4).forEach(function (t) {
      svg.appendChild(s('line', { class: t === 0 ? 'baseline' : 'gridline', x1: left, x2: width - right, y1: y(t), y2: y(t) }));
      svg.appendChild(s('text', { class: 'tick', x: left - 6, y: y(t) + 4, 'text-anchor': 'end', text: cfg.yFmt(t) }));
    });
    var every = band < 26 ? Math.ceil(26 / band) : 1;
    cfg.bins.forEach(function (b, i) {
      if (i % every === 0 || i === nb - 1) {
        var lx = cfg.edgeLabels && i < nb - 1 ? left + band * i : cx(i);
        svg.appendChild(s('text', { class: 'tick', x: lx, y: height - 8, 'text-anchor': 'middle', text: b }));
      }
    });
    if (ns === 1) {
      var bw = Math.min(24, band - 2);
      cfg.values[0].forEach(function (v, i) {
        var h = Math.max(0, y(0) - y(v || 0));
        if (v > 0) svg.appendChild(s('path', { class: 'mark', d: roundedBar(cx(i) - bw / 2, y(0) - Math.max(h, 1), bw, Math.max(h, 1), 'up'), fill: cfg.series[0].color }));
        var hit = s('rect', { class: 'hit', x: left + band * i, y: top, width: band, height: plotH });
        svg.appendChild(hit);
        hover(hit, function () { return [cfg.binTitles[i], [{ value: cfg.fmt(v, 0, i), label: '' }]]; });
      });
    } else {
      var cross = s('line', { class: 'crosshair', x1: 0, x2: 0, y1: top, y2: top + plotH, visibility: 'hidden' });
      svg.appendChild(cross);
      cfg.values.forEach(function (vals, si) {
        var d = '';
        vals.forEach(function (v, i) { if (v !== null && v !== undefined) d += (d ? 'L' : 'M') + cx(i) + ',' + y(v); });
        svg.appendChild(s('path', { class: 'series-line', d: d, stroke: cfg.series[si].color }));
      });
      cfg.values.forEach(function (vals, si) {
        vals.forEach(function (v, i) {
          if (v === null || v === undefined) return;
          svg.appendChild(s('circle', { class: 'dot', cx: cx(i), cy: y(v), r: 4, fill: cfg.series[si].color }));
        });
      });
      var overlay = s('rect', { class: 'hit', x: left, y: top, width: plotW, height: plotH });
      svg.appendChild(overlay);
      var at = function (i) {
        cross.setAttribute('x1', cx(i)); cross.setAttribute('x2', cx(i)); cross.setAttribute('visibility', 'visible');
        return [cfg.binTitles[i], cfg.series.map(function (sr, si) {
          return { color: sr.color, value: cfg.fmt(cfg.values[si][i], si, i), label: sr.label };
        })];
      };
      var idxFromEvent = function (e) {
        var r = svg.getBoundingClientRect();
        var px = (e.clientX - r.left) * (width / r.width);
        return Math.max(0, Math.min(nb - 1, Math.floor((px - left) / band)));
      };
      var focusIdx = 0;
      overlay.setAttribute('tabindex', '0');
      overlay.addEventListener('pointermove', function (e) { var t = at(idxFromEvent(e)); showTip(e, t[0], t[1]); });
      overlay.addEventListener('pointerleave', function () { cross.setAttribute('visibility', 'hidden'); hideTip(); });
      overlay.addEventListener('focus', function (e) { var t = at(focusIdx); showTip(e, t[0], t[1]); });
      overlay.addEventListener('blur', function () { cross.setAttribute('visibility', 'hidden'); hideTip(); });
      overlay.addEventListener('keydown', function (e) {
        if (e.key !== 'ArrowRight' && e.key !== 'ArrowLeft') return;
        focusIdx = Math.max(0, Math.min(nb - 1, focusIdx + (e.key === 'ArrowRight' ? 1 : -1)));
        var t = at(focusIdx); showTip({ target: overlay, type: 'focus' }, t[0], t[1]); e.preventDefault();
      });
    }
    container.appendChild(svg);
  }

  function histBins(h, def) {
    var e = h.edges, n = e.length, labels = [], titles = [];
    for (var i = 0; i < n; i++) {
      var last = i === n - 1;
      if (def.cont) {
        labels.push(last ? compact(e[i]) + '+' : compact(e[i]));
        titles.push(last ? '≥ ' + fmtInt(e[i]) + ' ' + def.unit : fmtInt(e[i]) + ' to < ' + fmtInt(e[i + 1]) + ' ' + def.unit);
      } else {
        var hi = last ? null : e[i + 1] - 1;
        labels.push(last ? e[i] + '+' : hi === e[i] ? String(e[i]) : e[i] + '–' + hi);
        titles.push((last ? e[i] + ' or more' : hi === e[i] ? String(e[i]) : e[i] + ' to ' + hi) + ' ' + def.unit);
      }
    }
    return { labels: labels, titles: titles };
  }

  // Draw charts after the card is in the DOM so they can measure their width.
  var pending = [];
  function chartBox(draw) {
    var box = el('div', { class: 'chart' });
    pending.push(function () { draw(box); });
    return box;
  }
  function flush() { var p = pending; pending = []; p.forEach(function (f) { f(); }); }

  // ------------------------------------------------------------------ metrics
  function catPct(b, cats) { return pct(sum(cats.map(function (c) { return b.categories[c]; })), b.summary.isoforms); }
  var METRICS = [
    { key: 'isoforms', label: 'Isoforms', kind: 'count', get: function (b) { return b.summary.isoforms; } },
    { key: 'genes', label: 'Genes (FLAIR)', short: 'Genes', kind: 'count', get: function (b) { return b.summary.genes; } },
    { key: 'loci', label: 'Loci', kind: 'count', get: function (b) { return b.summary.loci; } },
    { key: 'fsm', label: 'FSM', kind: 'pct', note: '% of isoforms', get: function (b) { return catPct(b, ['FSM']); } },
    { key: 'ism', label: 'ISM', kind: 'pct', note: '% of isoforms', get: function (b) { return catPct(b, ['ISM']); } },
    { key: 'novel', label: 'NIC + NNC', kind: 'pct', note: '% of isoforms', get: function (b) { return catPct(b, ['NIC', 'NNC']); } },
    { key: 'other', label: 'Outside known structures', kind: 'pct', note: 'genic, genic intron, antisense, intergenic, fusion', get: function (b) { return catPct(b, ['genic', 'genic_intron', 'antisense', 'intergenic', 'fusion']); } },
    { key: 'mono', label: 'Mono-exon', kind: 'pct', note: '% of isoforms', get: function (b) { return pct(b.summary.mono_exon, b.summary.isoforms); } },
    { key: 'ref_fsm', label: 'Ref. transcripts with FSM', short: 'Ref. tx FSM', kind: 'count', get: function (b) { return b.summary.ref_transcripts_fsm; }, of: function (b) { return b.summary.ref_transcripts_total; } },
    { key: 'ref_cim', label: 'Ref. transcripts with CIM', short: 'Ref. tx CIM', kind: 'count', get: function (b) { return b.summary.ref_transcripts_cim; }, of: function (b) { return b.summary.ref_transcripts_total; } },
    { key: 'ref_genes', label: 'Ref. genes detected', short: 'Ref. genes', kind: 'count', get: function (b) { return b.summary.ref_genes_detected; }, of: function (b) { return b.summary.ref_genes_total; } },
    { key: 'junc_total', label: 'Distinct junctions', kind: 'count', get: function (b) { return get(b, 'junctions.total'); } },
    { key: 'junc_known', label: 'Known junctions', short: 'Known junc.', kind: 'pct', note: '% of junctions', get: function (b) { return b.junctions ? pct(b.junctions.known, b.junctions.total) : null; } },
    { key: 'junc_nc', label: 'Novel non-canonical junctions', kind: 'count', get: function (b) { return b.junctions && b.junctions.motif_checked !== false ? b.junctions.novel_noncanonical : null; } },
    { key: 'wall', label: 'Wall time (min)', kind: 'num', get: function (b) { return get(b, 'compute.wall_minutes'); } },
    { key: 'cpu', label: 'CPU hours', kind: 'num', get: function (b) { return get(b, 'compute.cpu_hours'); } },
    { key: 'rss', label: 'Peak memory (GB)', kind: 'num', get: function (b) { return get(b, 'compute.peak_rss_gb'); } }
  ];
  var OVERVIEW_KEYS = ['isoforms', 'genes', 'fsm', 'novel', 'mono', 'ref_fsm', 'ref_cim', 'ref_genes', 'junc_known'];

  function metricValue(m, b) {
    var v = m.get(b);
    return v === undefined ? null : v;
  }
  function fmtMetric(m, v, b) {
    if (v === null) return '—';
    if (m.kind === 'pct') return fmtPct(v);
    if (m.kind === 'num') return fmtNum(v);
    var s0 = fmtInt(v);
    if (m.of && b) { var d = m.of(b); if (d) s0 += ' (' + fmtPct(pct(v, d), 0) + ')'; }
    return s0;
  }
  function fmtDelta(m, v, bv) {
    if (v === null || bv === null || bv === undefined) return '';
    var d = v - bv;
    if (m.kind === 'pct') return (d > 0 ? '+' : d < 0 ? '−' : '±') + Math.abs(d).toFixed(1) + ' pp';
    if (m.kind === 'num') return (d > 0 ? '+' : d < 0 ? '−' : '±') + Math.abs(d).toFixed(1);
    var rel = bv ? ' (' + (d >= 0 ? '+' : '−') + Math.abs(100 * d / bv).toFixed(0) + '%)' : '';
    return d === 0 ? '±0' : (d > 0 ? '+' : '−') + fmtInt(Math.abs(d)) + rel;
  }
  function metricsPresent(list) {
    return METRICS.filter(function (m) { return list.some(function (b) { return metricValue(m, b) !== null; }); });
  }

  // ------------------------------------------------------------------ filtering
  var versions = uniq(bundles.map(function (b) { return get(b, 'flair.version') || ''; }));
  var tagKeys = uniq([].concat.apply([], bundles.map(function (b) { return Object.keys(b.tags || {}); })));
  tagKeys = tagKeys.filter(function (k) {
    // a "version" tag that just repeats flair.version would be a duplicate filter
    return !(k === 'version' && bundles.every(function (b) { return !b.tags || b.tags.version === undefined || b.tags.version === get(b, 'flair.version'); }));
  });
  function uniq(a) { var seen = {}, out = []; a.forEach(function (x) { if (x !== '' && !seen[x]) { seen[x] = 1; out.push(x); } }); return out; }

  function filtered() {
    return bundles.filter(function (b) {
      if (state.filters.sample && b.sample !== state.filters.sample) return false;
      if (state.filters.version && get(b, 'flair.version') !== state.filters.version) return false;
      return Object.keys(state.tagFilters).every(function (k) {
        var v = state.tagFilters[k];
        return !v || (b.tags && b.tags[k] === v);
      });
    });
  }
  function baselineFor(sample) {
    return bundles.find(function (b) { return b.sample === sample && b.run_id === state.baseline; }) || null;
  }

  function selectField(label, value, options, onChange) {
    var sel = el('select', { 'aria-label': label });
    options.forEach(function (o) {
      var opt = el('option', { value: o[0], text: o[1] });
      if (o[0] === value) opt.selected = true;
      sel.appendChild(opt);
    });
    sel.addEventListener('change', function () { onChange(sel.value); });
    return el('div', { class: 'field' }, el('label', { text: label }), sel);
  }

  function renderFilters() {
    var box = document.getElementById('filters');
    box.textContent = '';
    if (state.view === 'help') return;
    if (samples.length > 1) box.appendChild(selectField('Sample', state.filters.sample,
      [['', 'All samples']].concat(samples.map(function (x) { return [x, x]; })),
      function (v) { state.filters.sample = v; render(); }));
    if (versions.length > 1) box.appendChild(selectField('FLAIR version', state.filters.version,
      [['', 'All versions']].concat(versions.map(function (x) { return [x, x]; })),
      function (v) { state.filters.version = v; render(); }));
    tagKeys.forEach(function (k) {
      var vals = uniq(bundles.map(function (b) { return (b.tags || {})[k] || ''; }));
      if (vals.length < 2) return;
      box.appendChild(selectField(k.charAt(0).toUpperCase() + k.slice(1).replace(/_/g, ' '), state.tagFilters[k] || '',
        [['', 'All']].concat(vals.map(function (x) { return [x, x]; })),
        function (v) { state.tagFilters[k] = v; render(); }));
    });
    if (state.view !== 'run' && runOrder.length > 1) box.appendChild(selectField('Baseline run', state.baseline,
      runOrder.map(function (r) { return [r, runInfo[r].label]; }),
      function (v) { state.baseline = v; render(); }));
  }

  // ------------------------------------------------------------------ views
  function runCell(b, withKey) {
    return el('div', { class: 'runcell' },
      withKey ? el('span', { class: 'key', style: 'background:' + runColor(b.run_id) }) : null,
      el('div', null, el('span', { class: 'name', text: b.label }),
        el('span', { class: 'meta', text: b.run_id + ' · FLAIR ' + (get(b, 'flair.version') || '?') })));
  }

  function emptyState(msg) { return el('div', { class: 'card empty', text: msg }); }

  function viewOverview(root) {
    var list = filtered();
    if (!list.length) { root.appendChild(emptyState('No runs match these filters.')); return; }
    var bySample = uniq(list.map(function (b) { return b.sample; }));
    var keys = OVERVIEW_KEYS.map(function (k) { return METRICS.find(function (m) { return m.key === k; }); });
    bySample.forEach(function (sample) {
      var rows = list.filter(function (b) { return b.sample === sample; });
      var base = baselineFor(sample);
      var ref0 = rows[0].summary;
      var region = get(rows[0], 'inputs.regions_intervals');
      var note = rows.length + (rows.length === 1 ? ' run' : ' runs') + ' · ' + fmtInt(ref0.ref_transcripts_total) +
        ' reference transcripts in scope' + (region && region.length ? ' (' + region.map(function (r) { return r[0] + ':' + fmtInt(r[1]) + '–' + fmtInt(r[2]); }).join(', ') + ')' : '');
      var card = el('section', { class: 'card' });
      var cmp = rows.length > 1 ? el('button', { class: 'btn', type: 'button', text: 'Compare these runs', onclick: function () { state.compareSample = sample; setView('compare'); } }) : null;
      card.appendChild(el('div', { class: 'card-head' }, el('div', null, el('h2', { text: sample }), el('p', { class: 'note', text: note })), cmp));
      card.appendChild(legend(MIX.map(function (m) { return { label: m.label, color: m.color }; })));
      var thead = el('tr', null, el('th', { text: 'Run' }),
        el('th', { text: 'Category mix', style: 'text-align:left;min-width:120px' }),
        keys.map(function (m) { return el('th', { class: 'num', text: m.short || m.label, title: m.label + (m.note ? ' (' + m.note + ')' : '') }); }));
      var tbody = el('tbody');
      rows.forEach(function (b) {
        var isBase = base && b.key === base.key;
        var tr = el('tr', { class: 'clickable' + (isBase ? ' is-baseline' : ''), tabindex: '0', title: 'Open this run' });
        var open = function () { state.runKey = b.key; setView('run'); };
        tr.addEventListener('click', open);
        tr.addEventListener('keydown', function (e) { if (e.key === 'Enter') open(); });
        tr.appendChild(el('td', null, runCell(b, false)));
        tr.appendChild(el('td', { style: 'text-align:left;min-width:120px' }, segbar(MIX.map(function (mx) {
          return { label: mx.label, value: sum(mx.cats.map(function (c) { return b.categories[c]; })), color: mx.color, on: mx.on };
        }), { thin: true, name: b.label + ' categories' })));
        keys.forEach(function (m) {
          var v = metricValue(m, b);
          var td = el('td', { class: 'num' }, fmtMetric(m, v, b));
          if (base && !isBase) td.appendChild(el('span', { class: 'delta', text: fmtDelta(m, v, metricValue(m, base)) }));
          tr.appendChild(td);
        });
        tbody.appendChild(tr);
      });
      card.appendChild(el('div', { class: 'table-wrap' }, el('table', null, el('thead', null, thead), tbody)));
      if (base) card.appendChild(el('p', { class: 'note', style: 'margin-top:8px', text: 'Small numbers under each value are the difference from the baseline run (' + base.label + '). pp = percentage points. Click a row to open that run.' }));
      root.appendChild(card);
    });
  }

  function provenance(b) {
    var f = b.flair || {};
    var items = [
      ['Sample', b.sample], ['Run id', b.run_id],
      ['FLAIR', (f.version || '?') + ' (' + f.source + ')'], ['Requested', f.requested || '—'],
      ['Commit', f.commit || '—'], ['Settings', get(b, 'inputs.transcriptome_args') || '—'],
      ['Reference', get(b, 'inputs.gtf') || '—'], ['Regions', get(b, 'inputs.regions') || 'whole genome'],
      ['Built', get(b, 'pipeline.created') || '—'], ['flair-eval', (get(b, 'pipeline.version') || '') + (get(b, 'pipeline.commit') ? ' @ ' + String(b.pipeline.commit).slice(0, 10) : '')]
    ];
    return el('dl', { class: 'kv' }, items.map(function (kv) { return [el('dt', { text: kv[0] }), el('dd', { class: kv[0] === 'Commit' ? 'mono' : null, text: kv[1] })]; }));
  }

  function tile(label, value, detail) {
    return el('div', { class: 'tile' }, el('div', { class: 'label', text: label }), el('div', { class: 'value', text: value }), detail ? el('div', { class: 'detail', text: detail }) : null);
  }

  function viewRun(root) {
    var list = filtered();
    if (!list.length) { root.appendChild(emptyState('No runs match these filters.')); return; }
    var b = list.find(function (x) { return x.key === state.runKey; }) ||
      list.find(function (x) { return x.run_id === state.baseline; }) || list[0];
    state.runKey = b.key;
    var picker = selectField('Run', b.key, list.map(function (x) { return [x.key, (samples.length > 1 ? x.sample + ' · ' : '') + x.label]; }),
      function (v) { state.runKey = v; render(); });
    var head = el('section', { class: 'card' },
      el('div', { class: 'card-head' }, el('div', null, el('h2', { text: b.label }),
        el('p', { class: 'note' }, Object.keys(b.tags || {}).map(function (k) { return el('span', { class: 'tag', text: k + ': ' + b.tags[k] }); }))), picker),
      el('details', null, el('summary', { text: 'Provenance' }), provenance(b)));
    root.appendChild(head);

    var S = b.summary;
    root.appendChild(el('div', { class: 'tiles' },
      tile('Isoforms', fmtInt(S.isoforms), fmtInt(S.mono_exon) + ' mono-exon'),
      tile('Genes', fmtInt(S.genes), fmtInt(S.loci) + ' loci'),
      tile('Ref. transcripts with FSM', fmtInt(S.ref_transcripts_fsm), 'of ' + fmtInt(S.ref_transcripts_total) + ' (' + fmtPct(pct(S.ref_transcripts_fsm, S.ref_transcripts_total)) + ')'),
      tile('Ref. transcripts with CIM', fmtInt(S.ref_transcripts_cim), 'both ends within 50 bp'),
      tile('Ref. genes detected', fmtInt(S.ref_genes_detected), S.ref_genes_total ? 'of ' + fmtInt(S.ref_genes_total) : null)));

    var grid = el('div', { class: 'grid' });
    root.appendChild(grid);
    var color = 'var(--series-1)';
    var catRows = CATS.map(function (c) { return [CAT_LONG[c] + ' (' + CAT_LABEL[c] + ')', fmtInt(b.categories[c]), fmtPct(pct(b.categories[c], S.isoforms))]; });
    grid.appendChild(panel('Structural categories', 'Isoforms by SQANTI category against the reference.', chartBox(function (box) {
      hbars(box, {
        name: 'Structural categories', series: [{ label: b.label, color: color }],
        groups: CATS.map(function (c) { return { label: CAT_LABEL[c], values: [b.categories[c]] }; }),
        labelValues: true, fmt: function (v) { return fmtInt(v) + ' · ' + fmtPct(pct(v, S.isoforms)); }
      });
    }), { head: ['Category', 'Isoforms', 'Share'], rows: catRows }));

    var right = el('section', { class: 'card' });
    var fsmTotal = b.categories.FSM;
    right.appendChild(el('h3', { text: 'FSM by end match' }));
    right.appendChild(el('p', { class: 'note', text: fmtInt(fsmTotal) + ' FSM isoforms, split by whether their ends fall within 50 bp of the matched reference.' }));
    var subParts = SUBS.map(function (k, i) { return { label: SUB_LABEL[k], value: b.fsm_subcategories[k], color: SLOT(i), on: ON_SLOT(i) }; });
    right.appendChild(legend(subParts));
    right.appendChild(segbar(subParts, { name: 'FSM subcategories' }));
    right.appendChild(el('h3', { text: 'Mono- vs multi-exon', style: 'margin-top:14px' }));
    var mm = [{ label: 'Multi-exon', value: S.multi_exon, color: SLOT(0) }, { label: 'Mono-exon', value: S.mono_exon, color: SLOT(1) }];
    right.appendChild(legend(mm));
    right.appendChild(segbar(mm, { name: 'Exon count' }));
    if (b.junctions) {
      var J = b.junctions;
      right.appendChild(el('h3', { text: 'Splice junctions', style: 'margin-top:14px' }));
      right.appendChild(el('p', { class: 'note', text: fmtInt(J.total) + ' distinct junctions across all isoforms.' }));
      var jp = junctionParts(J);
      right.appendChild(legend(jp));
      right.appendChild(segbar(jp, { name: 'Junctions' }));
    }
    right.appendChild(el('details', null, el('summary', { text: 'Numbers' }), dataTable(['', 'Count', 'Share'],
      subParts.map(function (p) { return [p.label, fmtInt(p.value), fmtPct(pct(p.value, fsmTotal))]; })
        .concat(mm.map(function (p) { return [p.label, fmtInt(p.value), fmtPct(pct(p.value, S.isoforms))]; }))
        .concat(b.junctions ? junctionParts(b.junctions).map(function (p) { return [p.label + ' junctions', fmtInt(p.value), fmtPct(pct(p.value, b.junctions.total))]; }) : []))));
    grid.appendChild(right);

    var hgrid = el('div', { class: 'grid' });
    var anyHist = false;
    HISTS.forEach(function (def) {
      var h = get(b, 'histograms.' + def.key);
      if (!h || !h.counts || !h.counts.length) return;
      anyHist = true;
      var bins = histBins(h, def), total = sum(h.counts);
      hgrid.appendChild(panel(def.title, fmtInt(total) + ' ' + def.of + '.', chartBox(function (box) {
        binChart(box, {
          name: def.title, bins: bins.labels, binTitles: bins.titles, edgeLabels: def.cont, series: [{ label: b.label, color: color }],
          values: [h.counts], yFmt: compact, fmt: function (v) { return fmtInt(v) + ' ' + def.of + ' (' + fmtPct(pct(v, total)) + ')'; }
        });
      }), { head: [def.unit, def.of, 'Share'], rows: h.counts.map(function (c, i) { return [bins.titles[i], fmtInt(c), fmtPct(pct(c, total))]; }) }));
    });
    if (anyHist) root.appendChild(hgrid);

    var comp = ['wall', 'cpu', 'rss'].map(function (k) { return METRICS.find(function (m) { return m.key === k; }); })
      .filter(function (m) { return metricValue(m, b) !== null; });
    if (comp.length) {
      root.appendChild(el('section', { class: 'card' }, el('h3', { text: 'Compute' }),
        el('div', { class: 'tiles', style: 'margin:8px 0 0' }, comp.map(function (m) { return tile(m.label, fmtMetric(m, metricValue(m, b))); }))));
    }
    sectionPanels([b], root, 'run');
  }

  function junctionParts(J) {
    if (J.motif_checked === false) {
      return [{ label: 'Known', value: J.known, color: SLOT(0) },
        { label: 'Novel (motif not checked)', value: J.novel !== undefined ? J.novel : J.total - J.known, color: SLOT(1) }];
    }
    return [{ label: 'Known', value: J.known, color: SLOT(0) },
      { label: 'Novel canonical', value: J.novel_canonical, color: SLOT(1) },
      { label: 'Novel non-canonical', value: J.novel_noncanonical, color: SLOT(2) }];
  }

  function sectionPanels(list, root, mode) {
    Object.keys(SECTION_PANELS).forEach(function (k) {
      var present = list.filter(function (b) { return b.sections && b.sections[k]; });
      if (!present.length) return;
      var p = SECTION_PANELS[k];
      if (mode === 'run' && p.renderRun) root.appendChild(p.renderRun(present[0]));
      if (mode === 'compare' && p.renderCompare) root.appendChild(p.renderCompare(present));
    });
  }

  function viewCompare(root) {
    var list = filtered();
    if (!list.length) { root.appendChild(emptyState('No runs match these filters.')); return; }
    var sampleList = uniq(list.map(function (b) { return b.sample; }));
    if (sampleList.indexOf(state.compareSample) < 0) state.compareSample = sampleList[0];
    var sample = state.compareSample;
    var pool = list.filter(function (b) { return b.sample === sample; });
    var ticked = pool.filter(function (b) { return !state.unticked[b.key]; });
    var base = baselineFor(sample);

    var ctl = el('section', { class: 'card' });
    var headLeft = el('div', null, el('h2', { text: 'Compare runs on ' + sample }),
      el('p', { class: 'note', text: 'Tick the runs to show. Runs are compared within one sample, so they see the same reads and reference.' }));
    var headRight = el('div', { class: 'panel-tools' });
    if (sampleList.length > 1) headRight.appendChild(selectField('Sample', sample, sampleList.map(function (x) { return [x, x]; }), function (v) { state.compareSample = v; render(); }));
    var modeBtns = [['share', 'Share'], ['count', 'Count']].map(function (m) {
      return el('button', { class: 'btn', type: 'button', 'aria-pressed': state.mode === m[0] ? 'true' : 'false', text: m[1], onclick: function () { state.mode = m[0]; render(); } });
    });
    headRight.appendChild(el('div', { class: 'field' }, el('label', { text: 'Show values as' }), el('div', { class: 'panel-tools' }, modeBtns)));
    ctl.appendChild(el('div', { class: 'card-head' }, headLeft, headRight));
    ctl.appendChild(el('div', { class: 'chips' }, pool.map(function (b) {
      var on = !state.unticked[b.key];
      var cb = el('input', { type: 'checkbox' });
      cb.checked = on;
      cb.addEventListener('change', function () { if (cb.checked) delete state.unticked[b.key]; else state.unticked[b.key] = 1; render(); });
      return el('label', { class: 'chip' + (on ? '' : ' off') }, cb, el('span', { class: 'key', style: 'background:' + runColor(b.run_id) + ';margin:0' }), b.label);
    })));
    if (runOrder.length > 8) ctl.appendChild(el('p', { class: 'note', style: 'margin-top:8px', text: 'Only the first 8 runs get their own colour; later runs are drawn in gray. Use the filters or tick fewer runs.' }));
    root.appendChild(ctl);
    if (!ticked.length) { root.appendChild(emptyState('Tick at least one run.')); return; }

    var series = ticked.map(function (b) { return { label: b.label, color: runColor(b.run_id) }; });
    var share = state.mode === 'share';

    // run x metric table with differences from the baseline
    var ms = metricsPresent(ticked);
    var head = ['Metric'].concat(ticked.map(function (b) { return b.label + (base && b.key === base.key ? ' (baseline)' : ''); }));
    var trs = ms.map(function (m) {
      return el('tr', null, el('td', { title: m.note || '' }, m.label, m.note ? el('span', { class: 'delta', text: m.note }) : null),
        ticked.map(function (b) {
          var v = metricValue(m, b), td = el('td', { class: 'num' }, fmtMetric(m, v, b));
          if (base && b.key !== base.key) td.appendChild(el('span', { class: 'delta', text: fmtDelta(m, v, metricValue(m, base)) || '' }));
          return td;
        }));
    });
    var tcard = el('section', { class: 'card' });
    tcard.appendChild(el('div', { class: 'card-head' }, el('div', null, el('h3', { text: 'Run × metric' }),
      el('p', { class: 'note', text: base ? 'Differences are against the baseline run ' + base.label + ' (change it in the filter row). pp = percentage points.' : 'The baseline run is not in this sample, so no differences are shown.' }))));
    tcard.appendChild(el('div', { class: 'table-wrap' }, el('table', null,
      el('thead', null, el('tr', null, head.map(function (h, i) {
        return el('th', { class: i ? 'num' : null }, i ? el('span', { class: 'key', style: 'background:' + runColor(ticked[i - 1].run_id) + ';margin:0 6px 0 0;vertical-align:-1px' }) : null, h);
      }))), el('tbody', null, trs))));
    root.appendChild(tcard);

    var leg = legend(series);
    var valFmt = function (count, total) { return share ? fmtPct(pct(count, total)) + ' (' + fmtInt(count) + ')' : fmtInt(count) + ' (' + fmtPct(pct(count, total)) + ')'; };
    var grid = el('div', { class: 'grid' });
    root.appendChild(grid);

    var shown = CATS.filter(function (c) { return ticked.some(function (b) { return b.categories[c] > 0; }); });
    var hidden = CATS.filter(function (c) { return shown.indexOf(c) < 0; });
    grid.appendChild(panel('Structural categories', (share ? 'Share of each run\'s isoforms.' : 'Isoform counts.') + (hidden.length ? ' Not shown (0 in every run): ' + hidden.map(function (c) { return CAT_LABEL[c]; }).join(', ') + '.' : ''),
      [leg.cloneNode(true), chartBox(function (box) {
        hbars(box, {
          name: 'Structural categories by run', series: series,
          groups: shown.map(function (c) { return { label: CAT_LABEL[c], cat: c, values: ticked.map(function (b) { return share ? pct(b.categories[c], b.summary.isoforms) : b.categories[c]; }) }; }),
          tickFmt: share ? function (t) { return t + '%'; } : compact,
          fmt: function (v, g, si) { var b = ticked[si]; return valFmt(b.categories[g.cat], b.summary.isoforms); }
        });
      })], { head: ['Category'].concat(ticked.map(function (b) { return b.label; })), rows: CATS.map(function (c) { return [CAT_LABEL[c]].concat(ticked.map(function (b) { return valFmt(b.categories[c], b.summary.isoforms); })); }) }));

    var rec = [
      { label: 'Ref. transcripts with FSM', get: function (b) { return b.summary.ref_transcripts_fsm; }, of: function (b) { return b.summary.ref_transcripts_total; } },
      { label: 'Ref. transcripts with CIM', get: function (b) { return b.summary.ref_transcripts_cim; }, of: function (b) { return b.summary.ref_transcripts_total; } },
      { label: 'Ref. genes detected', get: function (b) { return b.summary.ref_genes_detected; }, of: function (b) { return b.summary.ref_genes_total; } }
    ];
    var pctTick = function (t) { return t + '%'; };
    grid.appendChild(panel('Reference recovery', share ? 'Share of the reference transcripts / genes in scope that each run finds.' : 'Reference transcripts / genes found.',
      [leg.cloneNode(true), chartBox(function (box) {
        hbars(box, {
          name: 'Reference recovery', series: series,
          groups: rec.map(function (r) { return { label: r.label.replace('Ref. transcripts with ', 'Ref. tx ').replace('Ref. genes detected', 'Ref. genes'), r: r, values: ticked.map(function (b) { return share ? pct(r.get(b), r.of(b)) : r.get(b); }) }; }),
          tickFmt: share ? pctTick : compact,
          fmt: function (v, g, si) { var b = ticked[si]; return valFmt(g.r.get(b), g.r.of(b)); }
        });
      })], { head: [''].concat(ticked.map(function (b) { return b.label; })),
        rows: rec.map(function (r) { return [r.label].concat(ticked.map(function (b) { return valFmt(r.get(b), r.of(b)); })); }) }));
    grid.appendChild(panel('FSM by end match', (share ? 'Share of each run\'s FSM isoforms' : 'FSM isoforms') + ' whose ends fall within 50 bp of the reference (CIM) or not.',
      [leg.cloneNode(true), chartBox(function (box) {
        hbars(box, {
          name: 'FSM subcategories', series: series,
          groups: SUBS.map(function (k) { return { label: k, k: k, values: ticked.map(function (b) { return share ? pct(b.fsm_subcategories[k], b.categories.FSM) : b.fsm_subcategories[k]; }) }; }),
          tickFmt: share ? pctTick : compact,
          fmt: function (v, g, si) { var b = ticked[si]; return valFmt(b.fsm_subcategories[g.k], b.categories.FSM); }
        });
      })], { head: [''].concat(ticked.map(function (b) { return b.label; })),
        rows: SUBS.map(function (k) { return [SUB_LABEL[k]].concat(ticked.map(function (b) { return valFmt(b.fsm_subcategories[k], b.categories.FSM); })); }) }));

    // junctions + exon count: run-coloured grouped bars (in Compare, colour always means run)
    var jRows = [
      { label: 'Known', get: function (b) { return b.junctions && b.junctions.known; }, of: function (b) { return b.junctions && b.junctions.total; } },
      { label: 'Novel canonical', motif: true, get: function (b) { return b.junctions && b.junctions.motif_checked !== false ? b.junctions.novel_canonical : null; }, of: function (b) { return b.junctions && b.junctions.total; } },
      { label: 'Novel non-canonical', motif: true, get: function (b) { return b.junctions && b.junctions.motif_checked !== false ? b.junctions.novel_noncanonical : null; }, of: function (b) { return b.junctions && b.junctions.total; } },
      { label: 'Novel', nomotif: true, get: function (b) { return b.junctions && b.junctions.motif_checked === false ? (b.junctions.total - b.junctions.known) : null; }, of: function (b) { return b.junctions && b.junctions.total; } }
    ].filter(function (r) { return ticked.some(function (b) { var v = r.get(b); return v !== null && v !== undefined; }); });
    if (jRows.length) {
      grid.appendChild(panel('Splice junctions', (share ? 'Share of each run\'s distinct junctions' : 'Distinct junctions') + ', by whether the reference has them. Canonical = GT-AG, GC-AG or AT-AC.',
        [leg.cloneNode(true), chartBox(function (box) {
          hbars(box, {
            name: 'Splice junctions', series: series,
            groups: jRows.map(function (r) { return { label: r.label, r: r, values: ticked.map(function (b) { var v = r.get(b); return v === null || v === undefined ? null : share ? pct(v, r.of(b)) : v; }) }; }),
            tickFmt: share ? pctTick : compact,
            fmt: function (v, g, si) { var b = ticked[si]; return valFmt(g.r.get(b), g.r.of(b)); }
          });
        })], { head: ['Junctions'].concat(ticked.map(function (b) { return b.label; })),
          rows: jRows.map(function (r) { return [r.label].concat(ticked.map(function (b) { var v = r.get(b); return v === null || v === undefined ? '—' : valFmt(v, r.of(b)); })); }) }));
    }
    var exRows = [
      { label: 'Multi-exon', get: function (b) { return b.summary.multi_exon; } },
      { label: 'Mono-exon', get: function (b) { return b.summary.mono_exon; } }
    ];
    grid.appendChild(panel('Mono- vs multi-exon', share ? 'Share of each run\'s isoforms.' : 'Isoform counts.',
      [leg.cloneNode(true), chartBox(function (box) {
        hbars(box, {
          name: 'Exon count', series: series,
          groups: exRows.map(function (r) { return { label: r.label, r: r, values: ticked.map(function (b) { return share ? pct(r.get(b), b.summary.isoforms) : r.get(b); }) }; }),
          tickFmt: share ? pctTick : compact,
          fmt: function (v, g, si) { var b = ticked[si]; return valFmt(g.r.get(b), b.summary.isoforms); }
        });
      })], { head: [''].concat(ticked.map(function (b) { return b.label; })),
        rows: exRows.map(function (r) { return [r.label].concat(ticked.map(function (b) { return valFmt(r.get(b), b.summary.isoforms); })); }) }));

    // histograms as overlaid lines
    var hgrid = el('div', { class: 'grid' });
    var anyH = false;
    HISTS.forEach(function (def) {
      var withH = ticked.filter(function (b) { var h = get(b, 'histograms.' + def.key); return h && h.counts && h.counts.length; });
      if (!withH.length) return;
      anyH = true;
      var h0 = withH[0].histograms[def.key], bins = histBins(h0, def);
      var ser = withH.map(function (b) { return { label: b.label, color: runColor(b.run_id) }; });
      var totals = withH.map(function (b) { return sum(b.histograms[def.key].counts); });
      var vals = withH.map(function (b, si) { return b.histograms[def.key].counts.map(function (c) { return share ? pct(c, totals[si]) : c; }); });
      hgrid.appendChild(panel(def.title, share ? 'Share of ' + def.of + ' in each bin, so runs of different size compare.' : 'Number of ' + def.of + ' in each bin.',
        [legend(ser, 'line'), chartBox(function (box) {
          binChart(box, {
            name: def.title, bins: bins.labels, binTitles: bins.titles, edgeLabels: def.cont, series: ser, values: vals,
            yFmt: share ? function (t) { return t + '%'; } : compact,
            fmt: function (v, si, i) { var c = withH[si].histograms[def.key].counts[i]; return valFmt(c, totals[si]); }
          });
        })], { head: [def.unit].concat(withH.map(function (b) { return b.label; })), rows: bins.titles.map(function (t, i) { return [t].concat(withH.map(function (b, si) { return valFmt(b.histograms[def.key].counts[i], totals[si]); })); }) }));
    });
    if (anyH) root.appendChild(hgrid);

    var comp = ['wall', 'cpu', 'rss'].map(function (k) { return METRICS.find(function (m) { return m.key === k; }); })
      .filter(function (m) { return ticked.some(function (b) { return metricValue(m, b) !== null; }); });
    if (comp.length) {
      var cg = el('div', { class: 'grid' });
      comp.forEach(function (m) {
        cg.appendChild(panel(m.label, 'From the Nextflow trace of each run\'s jobs.', [leg.cloneNode(true), chartBox(function (box) {
          hbars(box, { name: m.label, series: series, groups: [{ label: m.label.replace(/ \(.*\)/, ''), values: ticked.map(function (b) { return metricValue(m, b); }) }], fmt: function (v) { return fmtNum(v); }, tickFmt: function (t) { return String(t); } });
        })], { head: ['Run', m.label], rows: ticked.map(function (b) { return [b.label, fmtNum(metricValue(m, b))]; }) }));
      });
      root.appendChild(cg);
    }
    sectionPanels(ticked, root, 'compare');
  }

  function viewHelp(root) {
    var card = el('section', { class: 'card help' });
    var p = function (t) { return el('p', { text: t }); };
    var dl = function (pairs) { return el('dl', null, pairs.map(function (x) { return [el('dt', { text: x[0] }), el('dd', { text: x[1] })]; })); };
    card.appendChild(el('h2', { text: 'How to read this report' }));
    card.appendChild(p('Each row in the report is one FLAIR run on one sample. A run is one FLAIR version, branch or setting. Every run on a sample sees the same reads, the same alignment and the same reference annotation, so differences come from FLAIR itself.'));
    card.appendChild(el('h3', { text: 'Views' }));
    card.appendChild(dl([
      ['Overview', 'One table per sample with the headline numbers for every run. The small line under each value is the difference from the baseline run. Click a row to open that run.'],
      ['Single run', 'All panels for one run: categories, end matches, junctions and length distributions, plus where the run came from (Provenance).'],
      ['Compare', 'Tick runs of one sample to draw them together. "Share" divides by each run\'s own total, so a run that reports more isoforms does not look better just for being bigger.'],
      ['Filters', 'The row above each view filters by sample, FLAIR version and any tags given in runs.yaml. The baseline run sets what the differences are measured against.']
    ]));
    card.appendChild(el('h3', { text: 'Structural categories (SQANTI3 names)' }));
    card.appendChild(p('Each isoform is compared with the reference annotation (GTF) by its intron chain: the ordered list of its splice junctions.'));
    card.appendChild(dl([
      ['FSM', 'Full splice match. Same intron chain as a reference transcript. For mono-exon isoforms: overlaps a mono-exon reference transcript on the same strand.'],
      ['ISM', 'Incomplete splice match. A consecutive part of a reference intron chain (for example a 5\' truncated copy), or a mono-exon isoform inside one reference exon.'],
      ['NIC', 'Novel in catalog. Every splice site is annotated but the combination is new, or an intron is retained.'],
      ['NNC', 'Novel not in catalog. Uses at least one splice site the reference does not have, inside a known spliced gene.'],
      ['Genic', 'Overlaps exons of a known gene but cannot be matched by splice sites (for example over a mono-exon gene).'],
      ['Genic intron', 'Lies entirely inside an intron of a known gene.'],
      ['Antisense', 'Overlaps exons of a known gene only on the opposite strand.'],
      ['Intergenic', 'Does not overlap any reference transcript.'],
      ['Fusion', 'Shares splice sites with two or more separate genes.']
    ]));
    card.appendChild(el('h3', { text: 'FSM by end match' }));
    card.appendChild(dl([
      ['CIM', 'Complete isoform match: both transcript ends within 50 bp of the matched reference transcript (SQANTI3 reference_match).'],
      ['FSM5', 'The 5\' end is more than 50 bp away; the 3\' end matches.'],
      ['FSM3', 'The 3\' end is more than 50 bp away; the 5\' end matches.'],
      ['FSMTE', 'Both ends are more than 50 bp away.']
    ]));
    card.appendChild(el('h3', { text: 'Other numbers' }));
    card.appendChild(dl([
      ['Genes, loci', 'Genes are FLAIR\'s own gene ids. Loci are clusters of isoforms that overlap on the same strand.'],
      ['Ref. transcripts with FSM / CIM', 'How many reference transcripts have at least one FSM (or CIM) isoform. "In scope" means reference transcripts in the evaluated regions.'],
      ['Ref. genes detected', 'Reference genes with at least one FSM, ISM, NIC, NNC or fusion isoform.'],
      ['Junctions', 'Distinct splice junctions across all isoforms. Known = in the reference. Novel junctions are canonical when the intron is GT-AG, GC-AG or AT-AC.'],
      ['Histograms', 'Fixed bins, so runs compare bin by bin. The last bin is open-ended (for example 10k+).'],
      ['Compute', 'CPU hours, wall time and peak memory from the Nextflow trace, when recorded.']
    ]));
    card.appendChild(el('h3', { text: 'Good to know' }));
    card.appendChild(p('More isoforms is not automatically better. Compare the share of FSM and CIM, reference recovery, and how many isoforms fall outside known structures. Only main chromosomes and spike-in contigs are evaluated (no chrM or scaffolds).'));
    root.appendChild(card);
  }

  // ------------------------------------------------------------------ shell
  var VIEWS = { overview: viewOverview, run: viewRun, compare: viewCompare, help: viewHelp };
  function setView(v) {
    state.view = VIEWS[v] ? v : 'overview';
    try { history.replaceState(null, '', '#' + state.view); } catch (e) { /* file:// in some browsers */ }
    render();
    window.scrollTo(0, 0);
  }
  var lastWidth = 0;
  function render() {
    hideTip();
    document.querySelectorAll('.tabs [data-view]').forEach(function (t) {
      t.setAttribute('aria-selected', t.getAttribute('data-view') === state.view ? 'true' : 'false');
    });
    renderFilters();
    var root = document.getElementById('view');
    root.textContent = '';
    pending = [];
    VIEWS[state.view](root);
    flush();
    lastWidth = root.clientWidth;
  }

  function initTheme() {
    var btn = document.getElementById('theme');
    var modes = ['auto', 'light', 'dark'], cur = 'auto';
    try { cur = localStorage.getItem('flair-eval-theme') || 'auto'; } catch (e) { cur = 'auto'; }
    var apply = function () {
      if (cur === 'auto') document.documentElement.removeAttribute('data-theme');
      else document.documentElement.setAttribute('data-theme', cur);
      btn.textContent = 'Theme: ' + cur;
    };
    btn.addEventListener('click', function () {
      cur = modes[(modes.indexOf(cur) + 1) % modes.length];
      try { localStorage.setItem('flair-eval-theme', cur); } catch (e) { /* storage blocked */ }
      apply();
    });
    apply();
  }

  function init() {
    var nRuns = runOrder.length, nS = samples.length;
    document.getElementById('subtitle').textContent = nRuns + (nRuns === 1 ? ' run' : ' runs') + ' × ' + nS + (nS === 1 ? ' sample' : ' samples') +
      ' · ' + bundles.length + (bundles.length === 1 ? ' bundle' : ' bundles') + ' · built ' + (DATA.generated || '').replace('T', ' ').replace(/\+00:00$/, ' UTC');
    document.getElementById('footer').textContent = 'Made by flair-eval ' + (DATA.tool_version || '') + '. This file is self-contained and works offline.';
    document.querySelectorAll('.tabs [data-view]').forEach(function (t) {
      t.addEventListener('click', function () { setView(t.getAttribute('data-view')); });
    });
    initTheme();
    var h = (location.hash || '').replace('#', '');
    state.view = VIEWS[h] ? h : (bundles.length === 1 ? 'run' : 'overview');
    if (!runInfo[state.baseline]) state.baseline = runOrder[0];
    render();
    var t = null;
    window.addEventListener('resize', function () {
      clearTimeout(t);
      t = setTimeout(function () {
        var w = document.getElementById('view').clientWidth;
        if (Math.abs(w - lastWidth) > 8) render();
      }, 150);
    });
  }
  init();
})();
