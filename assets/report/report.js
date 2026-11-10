/* flair-eval report.
 * Views: Overview, Structure, Agreement, Expression, Runtime, Run detail, Guide.
 * Plain JS + inline SVG; no libraries and no network. Data: <script id="report-data">.
 * Text from the data is inserted with textContent, never innerHTML.
 * Each run keeps one colour everywhere: categorical slot = its index in DATA.runs (baseline first). */
(function () {
  'use strict';

  var DATA = JSON.parse(document.getElementById('report-data').textContent);
  var SVGNS = 'http://www.w3.org/2000/svg';

  var CATS = ['FSM', 'ISM', 'NIC', 'NNC', 'genic', 'genic_intron', 'antisense', 'intergenic', 'fusion'];
  var CAT_LABEL = { FSM: 'FSM', ISM: 'ISM', NIC: 'NIC', NNC: 'NNC', genic: 'Genic', genic_intron: 'Genic intron', antisense: 'Antisense', intergenic: 'Intergenic', fusion: 'Fusion' };
  var CAT_LONG = { FSM: 'Full splice match', ISM: 'Incomplete splice match', NIC: 'Novel in catalog', NNC: 'Novel not in catalog', genic: 'Genic', genic_intron: 'Genic intron', antisense: 'Antisense', intergenic: 'Intergenic', fusion: 'Fusion' };
  var SUBS = ['CIM', 'FSM5', 'FSM3', 'FSMTE'];
  var SUB_LABEL = { CIM: 'CIM · both ends match', FSM5: "FSM5 · 5' end new", FSM3: "FSM3 · 3' end new", FSMTE: 'FSMTE · both ends new' };
  var PROD_ORDER = ['PRO', 'PTC', 'NST', 'NGO'];
  var PROD_LABEL = { PRO: 'PRO · productive', PTC: 'PTC · premature stop', NST: 'NST · no stop codon', NGO: 'NGO · no annotated start', NA: 'NA · not classified' };
  var HISTS = [
    { key: 'transcript_length', title: 'Transcript length', unit: 'nt', of: 'isoforms', cont: true },
    { key: 'exons_per_isoform', title: 'Exons per isoform', unit: 'exons', of: 'isoforms', cont: false },
    { key: 'exon_length', title: 'Exon length', unit: 'nt', of: 'exons', cont: true },
    { key: 'isoforms_per_gene', title: 'Isoforms per gene', unit: 'isoforms', of: 'genes', cont: false }
  ];
  function SLOT(i) { return i < 0 || i > 7 ? 'var(--series-other)' : 'var(--series-' + (i + 1) + ')'; }
  // Parts of one whole (end match, junction classes, coding classes): neutral ramp, dark to light.
  var PART_STEPS = [7, 6, 4, 3];
  function PART(i) { return 'var(--seq-' + PART_STEPS[Math.min(i, 3)] + ')'; }
  function ON_PART(i) { return 'var(--seq-on-' + PART_STEPS[Math.min(i, 3)] + ')'; }

  // Later report sections can register here: key -> {renderRun(bundle) -> node}. Unknown keys
  // in bundle.sections are ignored; a panel is drawn only when its data is present.
  var SECTION_PANELS = {};

  // ------------------------------------------------------------------ data
  // Display name: the run's label, or a readable form of its id (CONTRACTS A1 id patterns).
  function prettyId(id) {
    var m;
    if ((m = /^flair-(.+)$/.exec(id))) return 'FLAIR ' + m[1];
    if ((m = /^env-(.+)$/.exec(id))) return m[1] + ' (env)';
    if ((m = /^ctr-(.+)$/.exec(id))) return m[1] + ' (container)';
    return id;
  }
  var bundles = DATA.bundles.map(function (b) {
    var o = Object.assign({}, b);
    o.key = b.sample + '/' + b.run_id;
    o.label = b.label && b.label !== b.run_id ? b.label : prettyId(b.run_id);
    return o;
  });
  var runOrder = DATA.runs.slice();
  var samples = DATA.samples.slice();
  var bySample = DATA.by_sample || {};
  var runInfo = {};
  bundles.forEach(function (b) { if (!runInfo[b.run_id]) runInfo[b.run_id] = b; });
  function runIndex(id) { return runOrder.indexOf(id); }
  function runColor(id) { return SLOT(runIndex(id)); }
  function runLabel(id) { return runInfo[id] ? runInfo[id].label : id; }
  // "3.0.1", or "3.0.0+dev @ 389748e" for git/env/container builds; long strings are cut.
  function versionText(b, full) {
    var f = b.flair || {}, v = f.version || '?';
    if (!full && v.length > 18) v = v.slice(0, 17) + '…';
    if (f.commit && f.source !== 'release') v += ' @ ' + String(f.commit).slice(0, 7);
    return v;
  }
  // A21: how each run's quantify counted reads.
  var QUANT_INPUT = { reads: 'from reads', bam: 'from the genome BAM' };
  function quantInput(b) { return get(b, 'sections.quant.input') || null; }
  function quantNote(list) {
    var withQ = list.filter(function (b) { return get(b, 'sections.quant'); });
    var kinds = uniq(withQ.map(function (b) { return quantInput(b) || 'unknown'; }));
    var known = kinds.filter(function (k) { return k !== 'unknown'; });
    if (known.length < 2 && !(known.length && kinds.length > 1)) return null;
    var parts = kinds.map(function (k) {
      var names = withQ.filter(function (b) { return (quantInput(b) || 'unknown') === k; }).map(function (b) { return b.label; });
      return names.join(', ') + (k === 'unknown' ? ' (input not recorded)' : (names.length === 1 ? ' counts ' : ' count ') + QUANT_INPUT[k]);
    });
    return parts.join('; ') + '. Reads assigned and CPM depend on what FLAIR quantify counts, so compare these runs with care.';
  }
  function noteBox(text) { return text ? el('p', { class: 'callout', role: 'note' }, el('strong', { text: 'Not like for like. ' }), text) : null; }

  var state = {
    view: 'overview', sample: samples[0], baseline: DATA.baseline, hidden: {},
    filters: { version: '' }, tagFilters: {}, mode: 'share', runKey: null, printing: false
  };

  // ------------------------------------------------------------------ helpers
  function el(tag, attrs) {
    var node = document.createElement(tag);
    if (attrs) Object.keys(attrs).forEach(function (k) {
      var v = attrs[k];
      if (v === null || v === undefined || v === false) return;
      if (k === 'text') node.textContent = v;
      else if (k === 'class') node.className = v;
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
    return node;
  }
  // SVG text that fits `px` (about 6.8 px per character at 12 px); full text stays in <title>.
  function label(attrs, text, px) {
    var max = Math.max(4, Math.floor(px / 6.8));
    var t = s('text', Object.assign({}, attrs, { text: text.length > max ? text.slice(0, max - 1) + '…' : text }));
    if (text.length > max) { var ti = s('title', { text: text }); t.appendChild(ti); }
    return t;
  }
  function rowLabel(svg, x, y, text, px, dotColor) {
    if (dotColor) {
      svg.appendChild(s('circle', { cx: x - 4, cy: y - 4, r: 4, fill: dotColor }));
      x -= 14; px -= 14;
    }
    svg.appendChild(label({ class: 'cat-label', x: x, y: y, 'text-anchor': 'end' }, text, px));
  }
  function swatch(color, kind) { return el('span', { class: 'swatch' + (kind ? ' ' + kind : ''), style: 'background:' + color }); }
  var nf = new Intl.NumberFormat('en-US');
  function isNum(v) { return typeof v === 'number' && !isNaN(v); }
  function fmtInt(v) { return isNum(v) ? nf.format(Math.round(v)) : '—'; }
  function fmtPct(v, d) { return isNum(v) ? v.toFixed(d === undefined ? 1 : d) + '%' : '—'; }
  function fmtNum(v, d) { return isNum(v) ? v.toFixed(d === undefined ? 2 : d) : '—'; }
  // CPU time in hours from the bundle; shown in CPU-minutes when every value is under an hour.
  function cpuUnit(vals) { return vals.filter(isNum).every(function (v) { return v < 1; }) ? 'min' : 'h'; }
  function fmtCpu(v, unit) { return !isNum(v) ? '—' : unit === 'min' ? (v * 60).toFixed(v * 60 >= 10 ? 0 : 1) + ' CPU-min' : v.toFixed(2) + ' CPU-h'; }
  function fmtWall(v) { return !isNum(v) ? '—' : v >= 120 ? (v / 60).toFixed(1) + ' h' : v >= 10 ? v.toFixed(0) + ' min' : v.toFixed(1) + ' min'; }
  function fmtMem(v) { return !isNum(v) ? '—' : v < 1 ? Math.round(v * 1024) + ' MB' : v.toFixed(v >= 10 ? 0 : 1) + ' GB'; }
  function fmtSmart(v) { return !isNum(v) ? '—' : Math.abs(v) >= 100 ? fmtInt(v) : Math.abs(v) >= 10 ? v.toFixed(1) : Math.abs(v) >= 1 ? v.toFixed(2) : v.toFixed(3); }
  function compact(v) {
    if (!isNum(v)) return '';
    var a = Math.abs(v);
    if (a >= 1e6) return +(v / 1e6).toFixed(1) + 'M';
    if (a >= 1e3) return +(v / 1e3).toFixed(a >= 1e4 ? 0 : 1) + 'k';
    return String(+v.toFixed(2));
  }
  function pct(n, d) { return d ? 100 * n / d : null; }
  function sum(a) { return a.reduce(function (x, y) { return x + (y || 0); }, 0); }
  function get(o, path) { return path.split('.').reduce(function (x, k) { return x === null || x === undefined ? undefined : x[k]; }, o); }
  function uniq(a) { var seen = {}, out = []; a.forEach(function (x) { if (x !== '' && x !== undefined && !seen[x]) { seen[x] = 1; out.push(x); } }); return out; }
  function niceStep(span, n) {
    var raw = span / (n || 4), p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p;
    return (f <= 1 ? 1 : f <= 2 ? 2 : f <= 2.5 ? 2.5 : f <= 5 ? 5 : 10) * p;
  }
  function niceMax(v) { if (!(v > 0)) return 1; var st = niceStep(v, 4); return Math.ceil(v / st - 1e-9) * st; }
  function ticks(max, n) { var st = niceStep(max, n || 4), out = []; for (var t = 0; t <= max + 1e-9; t += st) out.push(+t.toFixed(10)); return out; }

  // ------------------------------------------------------------------ tooltip
  var tip = document.getElementById('tooltip');
  function showTip(evt, title, rows) {
    tip.textContent = '';
    if (title) tip.appendChild(el('div', { class: 'tt-title', text: title }));
    rows.forEach(function (r) {
      tip.appendChild(el('div', { class: 'tt-row' }, r.color ? swatch(r.color, 'line') : null,
        el('strong', { text: r.value }), r.label ? el('span', { class: 'lbl', text: r.label }) : null));
    });
    tip.hidden = false;
    var x, y;
    if (evt && evt.clientX !== undefined && evt.type !== 'focus') { x = evt.clientX; y = evt.clientY; }
    else { var rc = evt.target.getBoundingClientRect(); x = rc.left + rc.width / 2; y = rc.top; }
    var w = tip.offsetWidth, h = tip.offsetHeight;
    tip.style.left = Math.min(Math.max(8, x + 14), window.innerWidth - w - 8) + 'px';
    tip.style.top = (y - h - 12 < 8 ? y + 16 : y - h - 12) + 'px';
  }
  function hideTip() { tip.hidden = true; }
  function hover(node, fn) {
    node.setAttribute('tabindex', '0');
    node.addEventListener('pointermove', function (e) { var t = fn(); showTip(e, t[0], t[1]); });
    node.addEventListener('focus', function (e) { var t = fn(); showTip(e, t[0], t[1]); });
    node.addEventListener('pointerleave', hideTip);
    node.addEventListener('blur', hideTip);
  }

  // ------------------------------------------------------------------ layout blocks
  var pending = [];
  function chartBox(draw) { var box = el('div', { class: 'chart' }); pending.push(function () { draw(box); }); return box; }
  function flush() { var p = pending; pending = []; p.forEach(function (f) { f(); }); }

  function legend(items, kind) {
    return el('ul', { class: 'legend' }, items.map(function (it) { return el('li', null, swatch(it.color, kind), it.label); }));
  }
  function runLegend(ids, kind) { return legend(ids.map(function (id) { return { label: runLabel(id), color: runColor(id) }; }), kind); }

  // Card with title, one-line dek, optional chart/table toggle (keyboard: it is a button).
  function card(title, dek, body, opts) {
    opts = opts || {};
    var c = el('section', { class: 'card' + (opts.span ? ' span-2' : '') });
    var tools = el('div', { class: 'card-tools' });
    c.appendChild(el('div', { class: 'card-head' }, el('div', null, el('h3', { text: title }), dek ? el('p', { class: 'dek', text: dek }) : null), tools));
    var bodyBox = el('div', null, body);
    c.appendChild(bodyBox);
    if (opts.table) {
      var tbl = dataTable(opts.table.head, opts.table.rows);
      tbl.hidden = true;
      c.appendChild(tbl);
      var btn = el('button', { class: 'btn', type: 'button', 'aria-pressed': 'false', text: 'Table' });
      btn.addEventListener('click', function () {
        var on = tbl.hidden; tbl.hidden = !on; bodyBox.hidden = on;
        btn.setAttribute('aria-pressed', on ? 'true' : 'false');
        btn.textContent = on ? 'Chart' : 'Table';
      });
      tools.appendChild(btn);
    }
    return c;
  }
  function dataTable(head, rows) {
    return el('div', { class: 'table-wrap' }, el('table', null,
      el('thead', null, el('tr', null, head.map(function (h) { return el('th', { text: h }); }))),
      el('tbody', null, rows.map(function (r) { return el('tr', null, r.map(function (c) { return el('td', { text: c }); })); }))));
  }
  function sectionHead(title, dek, tools) {
    return el('div', { class: 'section-head' }, el('div', null, el('h2', { text: title }), dek ? el('p', { class: 'dek', text: dek }) : null), tools || null);
  }
  function grid(children, cls) { return el('div', { class: 'grid' + (cls ? ' ' + cls : '') }, children); }
  function emptyCard(msg) { return el('section', { class: 'card empty', text: msg }); }
  function seg(label, options, value, onPick) {
    var g = el('div', { class: 'seg', role: 'group', 'aria-label': label });
    options.forEach(function (o) {
      g.appendChild(el('button', { type: 'button', 'aria-pressed': o[0] === value ? 'true' : 'false', text: o[1], onclick: function () { onPick(o[0]); } }));
    });
    return g;
  }
  function selectField(label, value, options, onChange) {
    var sel = el('select', { 'aria-label': label });
    options.forEach(function (o) { var op = el('option', { value: o[0], text: o[1] }); if (o[0] === value) op.selected = true; sel.appendChild(op); });
    sel.addEventListener('change', function () { onChange(sel.value); });
    return el('label', { class: 'field' }, label, sel);
  }
  function runName(id) { return el('span', { class: 'runname' }, swatch(runColor(id)), el('strong', { text: runLabel(id) })); }

  // ------------------------------------------------------------------ charts
  function roundedBar(x, y, w, h, dir) {
    var r = Math.max(0, Math.min(4, dir === 'up' ? h : w, dir === 'up' ? w / 2 : h / 2));
    if (dir === 'right') return 'M' + x + ',' + y + 'H' + (x + w - r) + 'Q' + (x + w) + ',' + y + ' ' + (x + w) + ',' + (y + r) + 'V' + (y + h - r) + 'Q' + (x + w) + ',' + (y + h) + ' ' + (x + w - r) + ',' + (y + h) + 'H' + x + 'Z';
    if (dir === 'left') return 'M' + (x + w) + ',' + y + 'H' + (x + r) + 'Q' + x + ',' + y + ' ' + x + ',' + (y + r) + 'V' + (y + h - r) + 'Q' + x + ',' + (y + h) + ' ' + (x + r) + ',' + (y + h) + 'H' + (x + w) + 'Z';
    return 'M' + x + ',' + (y + h) + 'V' + (y + r) + 'Q' + x + ',' + y + ' ' + (x + r) + ',' + y + 'H' + (x + w - r) + 'Q' + (x + w) + ',' + y + ' ' + (x + w) + ',' + (y + r) + 'V' + (y + h) + 'Z';
  }

  // Grouped horizontal bars. cfg: {series:[{label,color}], groups:[{label, values:[...]}], fmt(v,g,si), tickFmt, colorOf?, labelValues?}
  function hbars(box, cfg) {
    var width = Math.max(280, box.clientWidth || 560);
    var ns = cfg.series.length;
    var barH = ns === 1 ? 16 : ns <= 3 ? 11 : ns <= 6 ? 8 : 6;
    var groupH = ns * barH + (ns - 1) * 2, pad = ns === 1 ? 10 : 14;
    var labelW = Math.min(cfg.labelW || 150, Math.max(84, width * (cfg.labelFrac || 0.3)));
    var valueW = cfg.labelValues ? 96 : 14, top = 2, axisH = 24;
    var plotW = width - labelW - valueW;
    var all = [];
    cfg.groups.forEach(function (g) { g.values.forEach(function (v) { if (isNum(v)) all.push(v); }); });
    var max = cfg.max || niceMax(Math.max.apply(null, all.concat([0])));
    var height = top + cfg.groups.length * (groupH + pad) + axisH;
    var svg = s('svg', { viewBox: '0 0 ' + width + ' ' + height, role: 'img', 'aria-label': cfg.name || 'bar chart' });
    var x = function (v) { return labelW + plotW * Math.max(0, v) / max; };
    ticks(max, width < 420 ? 3 : 4).forEach(function (t) {
      svg.appendChild(s('line', { class: t === 0 ? 'baseline' : 'gridline', x1: x(t), x2: x(t), y1: top, y2: height - axisH + 4 }));
      svg.appendChild(s('text', { class: 'tick', x: x(t), y: height - 6, 'text-anchor': 'middle', text: (cfg.tickFmt || compact)(t) }));
    });
    cfg.groups.forEach(function (g, gi) {
      var gy = top + gi * (groupH + pad) + pad / 2;
      rowLabel(svg, labelW - 10, gy + groupH / 2 + 4, g.label, labelW - 14, cfg.dotOf ? cfg.dotOf(gi) : null);
      g.values.forEach(function (v, si) {
        if (!isNum(v)) return;
        var y = gy + si * (barH + 2), w = x(v) - labelW;
        var color = cfg.colorOf ? cfg.colorOf(gi, si) : cfg.series[si].color;
        if (v > 0) svg.appendChild(s('path', { d: roundedBar(labelW, y, Math.max(w, 1.5), barH, 'right'), fill: color }));
        var hit = s('rect', { class: 'hit', x: labelW, y: y - 1, width: plotW, height: barH + 2 });
        svg.appendChild(hit);
        hover(hit, function () { return [g.label, [{ color: color, value: cfg.fmt(v, g, si), label: ns > 1 ? cfg.series[si].label : '' }]]; });
        if (cfg.labelValues) svg.appendChild(s('text', { class: 'value-label', x: x(v) + 6, y: y + barH / 2 + 4, text: cfg.short ? cfg.short(v, g, si) : cfg.fmt(v, g, si) }));
      });
    });
    box.appendChild(svg);
  }

  // Columns (one series) or lines (several) over categorical bins.
  function binChart(box, cfg) {
    var width = Math.max(260, box.clientWidth || 520), height = width < 420 ? 180 : 210;
    var left = 46, right = 10, top = 10, bottom = 26, plotW = width - left - right, plotH = height - top - bottom;
    var nb = cfg.bins.length, ns = cfg.series.length;
    var all = [].concat.apply([], cfg.values).filter(isNum);
    var max = niceMax(Math.max.apply(null, all.concat([0])));
    var band = plotW / nb;
    var y = function (v) { return top + plotH * (1 - v / max); };
    var cx = function (i) { return left + band * (i + 0.5); };
    var svg = s('svg', { viewBox: '0 0 ' + width + ' ' + height, role: 'img', 'aria-label': cfg.name || 'histogram' });
    ticks(max, 4).forEach(function (t) {
      svg.appendChild(s('line', { class: t === 0 ? 'baseline' : 'gridline', x1: left, x2: width - right, y1: y(t), y2: y(t) }));
      svg.appendChild(s('text', { class: 'tick', x: left - 8, y: y(t) + 4, 'text-anchor': 'end', text: cfg.yFmt(t) }));
    });
    var every = band < 28 ? Math.ceil(28 / band) : 1;
    cfg.bins.forEach(function (b, i) {
      if (i % every === 0 || i === nb - 1) {
        var lx = cfg.edgeLabels && i < nb - 1 ? left + band * i : cx(i);
        svg.appendChild(s('text', { class: 'tick', x: lx, y: height - 7, 'text-anchor': 'middle', text: b }));
      }
    });
    if (ns === 1) {
      var bw = Math.min(24, band - 2);
      cfg.values[0].forEach(function (v, i) {
        var h = Math.max(0, y(0) - y(v || 0));
        if (v > 0) svg.appendChild(s('path', { d: roundedBar(cx(i) - bw / 2, y(0) - Math.max(h, 1.5), bw, Math.max(h, 1.5), 'up'), fill: cfg.series[0].color }));
        var hit = s('rect', { class: 'hit', x: left + band * i, y: top, width: band, height: plotH });
        svg.appendChild(hit);
        hover(hit, function () { return [cfg.binTitles[i], [{ value: cfg.fmt(v, 0, i) }]]; });
      });
      box.appendChild(svg);
      return;
    }
    var cross = s('line', { class: 'crosshair', x1: 0, x2: 0, y1: top, y2: top + plotH, visibility: 'hidden' });
    svg.appendChild(cross);
    cfg.values.forEach(function (vals, si) {
      var d = '';
      vals.forEach(function (v, i) { if (isNum(v)) d += (d ? 'L' : 'M') + cx(i) + ',' + y(v); });
      svg.appendChild(s('path', { class: 'series-line', d: d, stroke: cfg.series[si].color }));
    });
    cfg.values.forEach(function (vals, si) {
      vals.forEach(function (v, i) { if (isNum(v)) svg.appendChild(s('circle', { class: 'dot', cx: cx(i), cy: y(v), r: 4, fill: cfg.series[si].color })); });
    });
    var overlay = s('rect', { class: 'hit', x: left, y: top, width: plotW, height: plotH });
    svg.appendChild(overlay);
    var at = function (i) {
      cross.setAttribute('x1', cx(i)); cross.setAttribute('x2', cx(i)); cross.setAttribute('visibility', 'visible');
      return [cfg.binTitles[i], cfg.series.map(function (sr, si) { return { color: sr.color, value: cfg.fmt(cfg.values[si][i], si, i), label: sr.label }; })];
    };
    var idx = function (e) { var r = svg.getBoundingClientRect(); return Math.max(0, Math.min(nb - 1, Math.floor(((e.clientX - r.left) * width / r.width - left) / band))); };
    var fi = 0;
    overlay.setAttribute('tabindex', '0');
    overlay.setAttribute('aria-label', (cfg.name || 'chart') + ': use the arrow keys to read each bin');
    overlay.addEventListener('pointermove', function (e) { var t = at(idx(e)); showTip(e, t[0], t[1]); });
    overlay.addEventListener('pointerleave', function () { cross.setAttribute('visibility', 'hidden'); hideTip(); });
    overlay.addEventListener('focus', function (e) { var t = at(fi); showTip(e, t[0], t[1]); });
    overlay.addEventListener('blur', function () { cross.setAttribute('visibility', 'hidden'); hideTip(); });
    overlay.addEventListener('keydown', function (e) {
      if (e.key !== 'ArrowRight' && e.key !== 'ArrowLeft') return;
      fi = Math.max(0, Math.min(nb - 1, fi + (e.key === 'ArrowRight' ? 1 : -1)));
      var t = at(fi); showTip({ target: overlay, type: 'focus' }, t[0], t[1]); e.preventDefault();
    });
    box.appendChild(svg);
  }

  // Diverging bars around a centre line: left = neutral (only in the baseline), right = run colour.
  function diverging(box, cfg) {
    var width = Math.max(280, box.clientWidth || 560);
    var labelW = Math.min(170, Math.max(90, width * 0.3)), rowH = 30, top = 22, bottom = 6;
    var plot = width - labelW - 8, half = plot / 2, mid = labelW + half;
    var max = niceMax(Math.max.apply(null, cfg.rows.map(function (r) { return Math.max(r.left, r.right); }).concat([1])));
    var height = top + cfg.rows.length * rowH + bottom;
    var svg = s('svg', { viewBox: '0 0 ' + width + ' ' + height, role: 'img', 'aria-label': cfg.name });
    svg.appendChild(s('text', { class: 'tick', x: mid - 8, y: 12, 'text-anchor': 'end', text: '← ' + cfg.leftLabel }));
    svg.appendChild(s('text', { class: 'tick', x: mid + 8, y: 12, 'text-anchor': 'start', text: cfg.rightLabel + ' →' }));
    svg.appendChild(s('line', { class: 'baseline', x1: mid, x2: mid, y1: top - 4, y2: height - bottom }));
    var sc = function (v) { return (half - 44) * v / max; };
    cfg.rows.forEach(function (r, i) {
      var y = top + i * rowH, bh = 14, by = y + (rowH - bh) / 2;
      rowLabel(svg, labelW - 12, by + bh / 2 + 4, r.label, labelW - 16, r.color);
      if (r.left > 0) svg.appendChild(s('path', { d: roundedBar(mid - sc(r.left), by, Math.max(1.5, sc(r.left)), bh, 'left'), fill: 'var(--neutral-bar)' }));
      if (r.right > 0) svg.appendChild(s('path', { d: roundedBar(mid, by, Math.max(1.5, sc(r.right)), bh, 'right'), fill: r.color }));
      svg.appendChild(s('text', { class: 'value-label', x: mid - sc(r.left) - 6, y: by + bh / 2 + 4, 'text-anchor': 'end', text: fmtInt(r.left) }));
      svg.appendChild(s('text', { class: 'value-label', x: mid + sc(r.right) + 6, y: by + bh / 2 + 4, text: fmtInt(r.right) }));
      var hit = s('rect', { class: 'hit', x: labelW, y: y, width: plot, height: rowH });
      svg.appendChild(hit);
      hover(hit, function () { return [r.label, r.tip]; });
    });
    box.appendChild(svg);
  }

  // Sequential colour step 0..7 for t in [0,1]; text colour that stays readable on it.
  function seqStep(t) { return Math.max(0, Math.min(7, Math.round(t * 7))); }
  function seqFill(t) { return 'var(--seq-' + seqStep(t) + ')'; }
  function seqInk(t) { return 'var(--seq-on-' + seqStep(t) + ')'; }

  // n x n matrix heatmap with values in the cells.
  function matrix(box, cfg) {
    var n = cfg.ids.length, width = Math.max(280, box.clientWidth || 520);
    var labelW = Math.min(170, Math.max(90, width * 0.32));
    var cell = Math.min(64, Math.floor((width - labelW - 4) / n));
    var topH = 4, height = topH + n * cell + 64;
    var svg = s('svg', { viewBox: '0 0 ' + width + ' ' + height, role: 'img', 'aria-label': cfg.name });
    var lo = cfg.min, hi = cfg.max;
    cfg.ids.forEach(function (a, i) {
      rowLabel(svg, labelW - 10, topH + i * cell + cell / 2 + 4, runLabel(a), labelW - 14, runColor(a));
      cfg.ids.forEach(function (b, j) {
        var x = labelW + j * cell, y = topH + i * cell;
        var v = cfg.value(a, b);
        if (a === b) {
          svg.appendChild(s('rect', { x: x + 1, y: y + 1, width: cell - 2, height: cell - 2, rx: 4, fill: 'var(--surface-2)' }));
          svg.appendChild(s('circle', { cx: x + cell / 2, cy: y + cell / 2, r: 4, fill: runColor(a) }));
          return;
        }
        var t = isNum(v) ? (v - lo) / ((hi - lo) || 1) : 0;
        svg.appendChild(s('rect', { x: x + 1, y: y + 1, width: cell - 2, height: cell - 2, rx: 4, fill: isNum(v) ? seqFill(t) : 'var(--surface-2)' }));
        if (cell >= 34) svg.appendChild(s('text', { class: 'cell-label', x: x + cell / 2, y: y + cell / 2 + 4, 'text-anchor': 'middle', fill: isNum(v) ? seqInk(t) : 'var(--muted)', text: isNum(v) ? cfg.fmt(v) : '—' }));
        var hit = s('rect', { class: 'hit', x: x, y: y, width: cell, height: cell });
        svg.appendChild(hit);
        hover(hit, function () { return [runLabel(a) + ' vs ' + runLabel(b), cfg.tip(a, b)]; });
      });
      svg.appendChild(s('circle', { cx: labelW + i * cell + cell / 2, cy: topH + n * cell + 12, r: 4, fill: runColor(a) }));
    });
    var sx = labelW, sy = topH + n * cell + 32, sw = Math.min(200, n * cell);
    for (var k = 0; k < 8; k++) svg.appendChild(s('rect', { x: sx + k * sw / 8, y: sy, width: sw / 8 - 1, height: 8, fill: 'var(--seq-' + k + ')' }));
    svg.appendChild(s('text', { class: 'tick', x: sx, y: sy + 22, text: cfg.fmt(lo) }));
    svg.appendChild(s('text', { class: 'tick', x: sx + sw, y: sy + 22, 'text-anchor': 'end', text: cfg.fmt(hi) }));
    box.appendChild(svg);
  }

  // Density of shared isoforms: x = baseline log10(CPM+1), y = run.
  function density(box, cfg) {
    var width = Math.max(200, box.clientWidth || 300), pad = 38, size = Math.min(width - pad - 8, 300);
    var n = cfg.grid.length, cell = size / n, W = size + pad + 8, height = size + 30;
    var svg = s('svg', { viewBox: '0 0 ' + W + ' ' + height, role: 'img', 'aria-label': cfg.name, style: 'max-width:' + W + 'px' });
    var mx = 0;
    cfg.grid.forEach(function (row) { row.forEach(function (v) { if (v > mx) mx = v; }); });
    svg.appendChild(s('rect', { x: pad, y: 4, width: size, height: size, fill: 'var(--seq-0)', rx: 4 }));
    var lo = function (k) { return Math.pow(10, k * cfg.logMax / n) - 1; };
    cfg.grid.forEach(function (row, j) {
      row.forEach(function (v, i) {
        if (!v) return;
        var t = 0.15 + 0.85 * Math.log(1 + v) / Math.log(1 + mx);
        var x = pad + i * cell, y = 4 + size - (j + 1) * cell;
        svg.appendChild(s('rect', { x: x, y: y, width: cell, height: cell, fill: seqFill(t) }));
        var hit = s('rect', { class: 'hit', x: x, y: y, width: cell, height: cell });
        svg.appendChild(hit);
        hover(hit, function () { return ['Baseline CPM ' + fmtSmart(lo(i)) + '–' + fmtSmart(lo(i + 1)) + ', run CPM ' + fmtSmart(lo(j)) + '–' + fmtSmart(lo(j + 1)), [{ value: fmtInt(v), label: 'shared isoforms' }]]; });
      });
    });
    svg.appendChild(s('line', { class: 'diag', x1: pad, y1: 4 + size, x2: pad + size, y2: 4 }));
    var step = cfg.logMax > 3 ? 1 : 0.5;
    for (var t = 0; t <= cfg.logMax + 1e-9; t += step) {
      var p = t / cfg.logMax * size, lab = t === 0 ? '0' : compact(Math.round(Math.pow(10, t)));
      svg.appendChild(s('text', { class: 'tick', x: pad + p, y: size + 20, 'text-anchor': 'middle', text: lab }));
      svg.appendChild(s('text', { class: 'tick', x: pad - 6, y: 4 + size - p + 4, 'text-anchor': 'end', text: lab }));
    }
    box.appendChild(svg);
  }

  // Tall segmented bar (one whole split into parts).
  function segbar(parts, opts) {
    opts = opts || {};
    var total = sum(parts.map(function (p) { return p.value; }));
    var bar = el('div', { class: 'segbar' + (opts.thin ? ' thin' : ''), role: 'img' });
    bar.setAttribute('aria-label', (opts.name ? opts.name + ': ' : '') + parts.map(function (p) { return p.label + ' ' + fmtInt(p.value) + ' (' + fmtPct(pct(p.value, total)) + ')'; }).join(', '));
    parts.forEach(function (p) {
      if (!p.value) return;
      var share = pct(p.value, total);
      var part = el('div', { class: 'part', style: 'flex:' + p.value + ' 1 0;background:' + p.color + ';color:' + (p.on || 'var(--on-series)') });
      bar.appendChild(part);
      if (!opts.thin) pending.push(function () {
        var txt = share >= 9.5 ? Math.round(share) + '%' : share.toFixed(1) + '%';
        if (part.clientWidth > txt.length * 7.5 + 16) part.textContent = txt;
      });
      hover(part, function () { return [opts.name || '', [{ color: p.color, value: fmtInt(p.value) + ' (' + fmtPct(share) + ')', label: p.label }]]; });
    });
    return bar;
  }

  // Dot strip: one row per metric, each run a dot on the row's own axis.
  function dotStrip(rows, ids, baseId) {
    var wrap = el('div', { class: 'strip' });
    rows.forEach(function (row) {
      var vals = ids.map(function (id) { return { id: id, v: row.get(id) }; }).filter(function (d) { return isNum(d.v); });
      if (!vals.length) return;
      var lo = Math.min.apply(null, vals.map(function (d) { return d.v; }));
      var hi = Math.max.apply(null, vals.map(function (d) { return d.v; }));
      wrap.appendChild(el('div', { class: 'row-label' }, row.label, row.sub ? el('span', { class: 'small', text: row.sub }) : null));
      var track = el('div', { class: 'row-track' });
      wrap.appendChild(track);
      wrap.appendChild(el('div', { class: 'row-range', text: lo === hi ? row.fmt(lo) : rangeText(row.fmt(lo), row.fmt(hi)) }));
      pending.push(function () {
        var width = Math.max(120, track.clientWidth), padX = 12;
        var span = hi - lo || Math.abs(hi) || 1;
        var dlo = lo - span * 0.08, dhi = hi + span * 0.08;
        if (lo === hi) { dlo = lo - span; dhi = hi + span; }
        var x = function (v) { return padX + (width - 2 * padX) * (v - dlo) / (dhi - dlo); };
        var sorted = vals.slice().sort(function (a, b) { return a.v - b.v; });
        var lanes = [];
        sorted.forEach(function (d) {
          var px = x(d.v), l = 0;
          while (lanes[l] !== undefined && px - lanes[l] < 13) l++;
          lanes[l] = px; d.lane = l; d.px = px;
        });
        var nl = Math.max(1, lanes.length), height = Math.max(54, 26 + nl * 13);
        var svg = s('svg', { viewBox: '0 0 ' + width + ' ' + height, role: 'img', 'aria-label': row.label + ': ' + sorted.map(function (d) { return runLabel(d.id) + ' ' + row.fmt(d.v); }).join(', ') });
        var cy0 = height / 2 - (nl - 1) * 13 / 2;
        svg.appendChild(s('line', { class: 'gridline', x1: padX, x2: width - padX, y1: height / 2, y2: height / 2 }));
        sorted.forEach(function (d) {
          svg.appendChild(s('circle', { class: 'dot' + (d.id === baseId ? ' base' : ''), cx: d.px, cy: cy0 + d.lane * 13, r: 6, fill: runColor(d.id) }));
        });
        var hit = s('rect', { class: 'hit', x: 0, y: 0, width: width, height: height });
        svg.appendChild(hit);
        hover(hit, function () {
          return [row.label, sorted.slice().reverse().map(function (d) { return { color: runColor(d.id), value: row.fmt(d.v), label: runLabel(d.id) + (d.id === baseId ? ' (baseline)' : '') }; })];
        });
        track.appendChild(el('div', { class: 'chart' }, svg));
      });
    });
    return wrap;
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

  // ------------------------------------------------------------------ metrics
  function catPct(b, cats) { return pct(sum(cats.map(function (c) { return b.categories[c]; })), b.summary.isoforms); }
  function concFor(sample) { return (bySample[sample] || {}).concordance || null; }
  function quantFor(sample) { return (bySample[sample] || {}).quant || null; }
  function pairOf(block, a, b) {
    if (!block) return null;
    for (var i = 0; i < block.pairs.length; i++) {
      var p = block.pairs[i];
      if ((p.a === a && p.b === b) || (p.a === b && p.b === a)) return p;
    }
    return null;
  }
  // Genes one run detects that another does not, from the per-run lists of varying genes.
  var MAX_GENE_NAMES = 40;
  function genesVsBase(conc, id, base) {
    var gd = conc && conc.genes_detected;
    if (!gd || !gd.runs || !gd.runs[id] || !gd.runs[base]) return null;
    function minus(a, b) {
      var has = {};
      b.forEach(function (i) { has[i] = 1; });
      return a.filter(function (i) { return !has[i]; });
    }
    var gained = minus(gd.runs[id], gd.runs[base]), lost = minus(gd.runs[base], gd.runs[id]);
    var names = function (idx) { return idx.slice(0, MAX_GENE_NAMES).map(function (i) { return gd.names[i]; }); };
    return { gained: gained.length, lost: lost.length, gained_names: names(gained), lost_names: names(lost) };
  }
  // Density grid of a pair with the baseline across (x) and the other run up (y).
  function gridVsBase(p, base) {
    if (!p || !p.grid) return null;
    if (p.a === base) return p.grid;
    return p.grid[0].map(function (_, i) { return p.grid.map(function (row) { return row[i]; }); });
  }

  // Why the baseline is the baseline (CONTRACTS A22).
  var BASE_RULE = DATA.baseline_rule || {};
  function baselineWhy(id) {
    if (id !== DATA.baseline) return 'chosen in this page';
    if (BASE_RULE.rule === 'given') return 'set when the report was built';
    if (BASE_RULE.rule === 'earliest_launch') return BASE_RULE.launches > 1 ? 'first run of the earliest of ' + BASE_RULE.launches + ' launches' : 'first run of the launch';
    if (BASE_RULE.rule === 'name_order') return 'first run id in name order';
    return 'default';
  }
  function baselineNote() {
    var switched = state.baseline !== DATA.baseline;
    return el('p', { class: 'basenote' }, el('span', { class: 'label', text: 'Baseline' }), runName(state.baseline),
      el('span', { text: ' · ' + baselineWhy(state.baseline) + (switched ? ' (default: ' + runLabel(DATA.baseline) + ', ' + baselineWhy(DATA.baseline) + ')' : '') +
        '. Every difference is measured against it' + (state.printing ? '.' : '; change it with Baseline above.') }),
      switched ? el('button', { class: 'linkbtn no-print', type: 'button', text: 'Use default', onclick: function () { state.baseline = DATA.baseline; render(); } }) : null);
  }
  function setSubtitle() {
    var nR = runOrder.length, nS = samples.length;
    var when = (DATA.generated || '').slice(0, 10);
    document.getElementById('subtitle').textContent = nR + ' FLAIR ' + (nR === 1 ? 'run' : 'runs') + ' on ' + nS + ' ' + (nS === 1 ? 'sample' : 'samples') +
      ' · baseline ' + runLabel(state.baseline) + ', ' + baselineWhy(state.baseline) + (when ? ' · built ' + when : '');
  }

  // "3.1 CPU-min – 12 CPU-min" -> "3.1–12 CPU-min" when both ends share a unit.
  function rangeText(a, b) {
    var ma = /^(\S+) (.+)$/.exec(a), mb = /^(\S+) (.+)$/.exec(b);
    return ma && mb && ma[2] === mb[2] ? ma[1] + '–' + mb[1] + ' ' + ma[2] : a + ' – ' + b;
  }

  function jaccardVsBase(b) {
    if (b.run_id === state.baseline) return null;
    var p = pairOf(concFor(b.sample), b.run_id, state.baseline);
    return p ? p.jaccard : null;
  }
  function rVsBase(b) {
    if (b.run_id === state.baseline) return null;
    var p = pairOf(quantFor(b.sample), b.run_id, state.baseline);
    return p ? p.pearson : null;
  }
  function prodPct(b, cls) {
    var cc = get(b, 'sections.productivity.class_counts');
    if (!cc) return null;
    return pct(cc[cls] || 0, sum(Object.keys(cc).map(function (k) { return cc[k]; })));
  }
  var METRICS = [
    { key: 'isoforms', label: 'Isoforms', kind: 'count', get: function (b) { return b.summary.isoforms; } },
    { key: 'genes', label: 'Genes', sub: 'FLAIR gene ids', kind: 'count', get: function (b) { return b.summary.genes; } },
    { key: 'fsm', label: 'FSM', sub: '% of isoforms', kind: 'pct', get: function (b) { return catPct(b, ['FSM']); } },
    { key: 'novel', label: 'NIC + NNC', sub: '% of isoforms', kind: 'pct', get: function (b) { return catPct(b, ['NIC', 'NNC']); } },
    { key: 'other', label: 'Outside known structures', sub: '% of isoforms', kind: 'pct', get: function (b) { return catPct(b, ['genic', 'genic_intron', 'antisense', 'intergenic', 'fusion']); } },
    { key: 'mono', label: 'Mono-exon', sub: '% of isoforms', kind: 'pct', get: function (b) { return pct(b.summary.mono_exon, b.summary.isoforms); } },
    { key: 'ref_fsm', label: 'Ref. transcripts with FSM', short: 'Ref. tx FSM', kind: 'count', get: function (b) { return b.summary.ref_transcripts_fsm; }, of: function (b) { return b.summary.ref_transcripts_total; } },
    { key: 'ref_cim', label: 'Ref. transcripts with CIM', short: 'Ref. tx CIM', sub: 'both ends within 50 bp', kind: 'count', get: function (b) { return b.summary.ref_transcripts_cim; }, of: function (b) { return b.summary.ref_transcripts_total; } },
    { key: 'ref_genes', label: 'Ref. genes detected', short: 'Ref. genes', kind: 'count', get: function (b) { return b.summary.ref_genes_detected; }, of: function (b) { return b.summary.ref_genes_total; } },
    { key: 'junc_known', label: 'Known junctions', short: 'Known junc.', sub: '% of distinct junctions', kind: 'pct', get: function (b) { return b.junctions ? pct(b.junctions.known, b.junctions.total) : null; } },
    { key: 'jaccard', label: 'Structures shared with baseline', sub: 'Jaccard index', kind: 'ratio', get: jaccardVsBase },
    { key: 'cpm_r', label: 'CPM agreement with baseline', sub: 'Pearson r, log CPM', kind: 'ratio', get: rVsBase },
    { key: 'pro', label: 'Productive (PRO)', sub: '% of isoforms', kind: 'pct', get: function (b) { return prodPct(b, 'PRO'); } },
    { key: 'cpu', label: 'CPU time', sub: 'FLAIR tasks', kind: 'cpu', get: function (b) { return get(b, 'compute.cpu_hours'); } },
    { key: 'wall', label: 'Wall time', kind: 'wall', get: function (b) { return get(b, 'compute.wall_minutes'); } },
    { key: 'rss', label: 'Peak memory', kind: 'mem', get: function (b) { return get(b, 'compute.peak_rss_gb'); } }
  ];
  function M(key) { return METRICS.find(function (m) { return m.key === key; }); }
  function mval(m, b) { var v = m.get(b); return isNum(v) ? v : null; }
  function fmtMetric(m, v, b, withOf) {
    if (!isNum(v)) return '—';
    if (m.kind === 'pct') return fmtPct(v);
    if (m.kind === 'ratio') return v.toFixed(2);
    if (m.kind === 'num') return fmtSmart(v);
    if (m.kind === 'cpu') return fmtCpu(v, v < 1 ? 'min' : 'h');
    if (m.kind === 'wall') return fmtWall(v);
    if (m.kind === 'mem') return fmtMem(v);
    var out = fmtInt(v);
    if (withOf !== false && m.of && b) { var d = m.of(b); if (d) out += ' (' + fmtPct(pct(v, d), 0) + ')'; }
    return out;
  }
  function fmtDelta(m, v, bv) {
    if (!isNum(v) || !isNum(bv)) return '';
    var d = v - bv, sign = d > 0 ? '+' : d < 0 ? '−' : '±';
    if (m.kind === 'pct') return sign + Math.abs(d).toFixed(1) + ' pp';
    if (m.kind === 'ratio') return sign + Math.abs(d).toFixed(2);
    if (m.kind === 'num') return sign + fmtSmart(Math.abs(d));
    if (m.kind === 'cpu' || m.kind === 'wall' || m.kind === 'mem') return sign + fmtMetric(m, Math.abs(d));
    if (d === 0) return '±0';
    return sign + fmtInt(Math.abs(d)) + (bv ? ' (' + sign + Math.abs(100 * d / bv).toFixed(0) + '%)' : '');
  }
  function present(m, list) { return list.some(function (b) { return mval(m, b) !== null; }); }

  // ------------------------------------------------------------------ filtering
  var versions = uniq(bundles.map(function (b) { return get(b, 'flair.version') || ''; }));
  // A tag becomes a filter only when it groups runs: two or more values, at least one shared by
  // several runs, and not just a copy of the FLAIR version or source.
  var tagKeys = uniq([].concat.apply([], bundles.map(function (b) { return Object.keys(b.tags || {}); }))).filter(function (k) {
    var byRun = {};
    bundles.forEach(function (b) { if (b.tags && b.tags[k] !== undefined) byRun[b.run_id] = b.tags[k]; });
    var vals = Object.keys(byRun).map(function (r) { return byRun[r]; });
    var counts = {};
    vals.forEach(function (v) { counts[v] = (counts[v] || 0) + 1; });
    var groups = Object.keys(counts).length >= 2 && Object.keys(counts).some(function (v) { return counts[v] >= 2; });
    var dup = (k === 'version' || k === 'source') && bundles.every(function (b) { return !b.tags || b.tags[k] === undefined || b.tags[k] === get(b, 'flair.' + k); });
    return groups && !dup;
  });
  function passes(b) {
    if (state.filters.version && get(b, 'flair.version') !== state.filters.version) return false;
    return Object.keys(state.tagFilters).every(function (k) { var v = state.tagFilters[k]; return !v || (b.tags && b.tags[k] === v); });
  }
  function sampleRuns(sample) {
    return bundles.filter(function (b) { return b.sample === sample && passes(b); })
      .sort(function (a, b) { return runIndex(a.run_id) - runIndex(b.run_id); });
  }
  function visible(sample) { return sampleRuns(sample).filter(function (b) { return !state.hidden[b.run_id]; }); }
  function baseFor(sample) { return bundles.find(function (b) { return b.sample === sample && b.run_id === state.baseline; }) || null; }

  // ------------------------------------------------------------------ run bar
  function renderRunbar() {
    var bar = document.getElementById('runbar');
    bar.textContent = '';
    if (samples.length > 1) {
      var pick = function (v) { state.sample = v; state.runKey = null; render(); };
      bar.appendChild(el('div', { class: 'field' }, el('span', { class: 'label', text: 'Sample' }),
        samples.length <= 4 ? seg('Sample', samples.map(function (x) { return [x, x]; }), state.sample, pick)
          : selectField('', state.sample, samples.map(function (x) { return [x, x]; }), pick)));
    }
    var here = sampleRuns(state.sample).map(function (b) { return b.run_id; });
    var chips = el('div', { class: 'chips', role: 'group', 'aria-label': 'Runs shown (select to hide or show)' });
    runOrder.forEach(function (id) {
      if (here.indexOf(id) < 0) return;
      var on = !state.hidden[id];
      chips.appendChild(el('button', {
        class: 'chip', type: 'button', 'aria-pressed': on ? 'true' : 'false', title: (on ? 'Hide ' : 'Show ') + runLabel(id),
        onclick: function () { if (on) state.hidden[id] = 1; else delete state.hidden[id]; render(); }
      }, swatch(runColor(id)), runLabel(id), id === state.baseline ? el('span', { class: 'base', text: 'baseline' }) : null));
    });
    bar.appendChild(el('div', { class: 'field' }, el('span', { class: 'label', text: 'Runs' }), chips));
    var more = el('div', { class: 'more' });
    if (runOrder.length > 1) {
      var bsel = selectField('Baseline', state.baseline, runOrder.map(function (r) { return [r, runLabel(r) + (r === DATA.baseline ? ' (default)' : '')]; }), function (v) { state.baseline = v; render(); });
      bsel.title = 'Run every difference is measured against. Default: ' + runLabel(DATA.baseline) + ', ' + baselineWhy(DATA.baseline) + '.';
      more.appendChild(bsel);
    }
    // A "version" run tag (dev, ted, 3.0.1, ...) is finer than the version FLAIR reports, so it replaces that filter.
    if (versions.length > 1 && versions.length < bundles.length && tagKeys.indexOf('version') < 0) more.appendChild(selectField('FLAIR version', state.filters.version, [['', 'All']].concat(versions.map(function (x) { return [x, x.length > 18 ? x.slice(0, 17) + '…' : x]; })), function (v) { state.filters.version = v; render(); }));
    tagKeys.forEach(function (k) {
      var vals = uniq(bundles.map(function (b) { return (b.tags || {})[k] || ''; }));
      if (vals.length < 2) return;
      more.appendChild(selectField(k.charAt(0).toUpperCase() + k.slice(1).replace(/_/g, ' '), state.tagFilters[k] || '',
        [['', 'All']].concat(vals.map(function (x) { return [x, x]; })), function (v) { state.tagFilters[k] = v; render(); }));
    });
    bar.appendChild(more);
  }

  // ------------------------------------------------------------------ Overview
  function highlights(list) {
    var items = [];
    if (list.length < 2) return items;
    function ext(fn) {
      var vals = list.map(function (b) { return { b: b, v: fn(b) }; }).filter(function (d) { return isNum(d.v); });
      if (vals.length < 2) return null;
      vals.sort(function (a, b) { return b.v - a.v; });
      return vals[0].v === vals[vals.length - 1].v ? null : { hi: vals[0], lo: vals[vals.length - 1] };
    }
    var e = ext(function (b) { return b.summary.isoforms; });
    if (e) items.push(['Size', [runName(e.hi.b.run_id), ' reports the most isoforms (' + fmtInt(e.hi.v) + '); ', runName(e.lo.b.run_id), ' the fewest (' + fmtInt(e.lo.v) + ').']]);
    e = ext(function (b) { return b.summary.ref_transcripts_cim; });
    if (e) items.push(['Reference', [runName(e.hi.b.run_id), ' matches the most reference transcripts end to end (CIM): ' + fmtInt(e.hi.v) + ' of ' + fmtInt(e.hi.b.summary.ref_transcripts_total) + ', against ' + fmtInt(e.lo.v) + ' for ', runName(e.lo.b.run_id), '.']]);
    e = ext(function (b) { return catPct(b, ['NIC', 'NNC']); });
    if (e) items.push(['Novelty', [runName(e.hi.b.run_id), ' has the largest share of novel isoforms (NIC + NNC ' + fmtPct(e.hi.v) + '); ', runName(e.lo.b.run_id), ' the smallest (' + fmtPct(e.lo.v) + ').']]);
    var conc = concFor(list[0].sample);
    if (conc) {
      var ids = list.map(function (b) { return b.run_id; });
      var ps = conc.pairs.filter(function (p) { return ids.indexOf(p.a) >= 0 && ids.indexOf(p.b) >= 0 && isNum(p.jaccard); });
      if (ps.length) {
        ps.sort(function (a, b) { return b.jaccard - a.jaccard; });
        var hi = ps[0], lo = ps[ps.length - 1];
        var parts = ['Most alike: ', runName(hi.a), ' and ', runName(hi.b), ' (Jaccard ' + hi.jaccard.toFixed(2) + ')'];
        parts = ps.length > 1 ? parts.concat(['. Least alike: ', runName(lo.a), ' and ', runName(lo.b), ' (' + lo.jaccard.toFixed(2) + ').']) : parts.concat(['.']);
        items.push(['Agreement', parts]);
      }
    }
    var rs = list.map(rVsBase).filter(isNum);
    if (rs.length) items.push(['Expression', ['On isoforms shared with the baseline, CPM correlates at r = ' + Math.min.apply(null, rs).toFixed(2) + (rs.length > 1 ? '–' + Math.max.apply(null, rs).toFixed(2) : '') + ' (log scale).' + (quantNote(list) ? ' Quantification inputs differ between runs; see Expression.' : '')]]);
    e = ext(function (b) { return get(b, 'compute.cpu_hours'); });
    if (e) { var cu = cpuUnit([e.lo.v, e.hi.v]); items.push(['Cost', [runName(e.lo.b.run_id), ' used the least CPU time (' + fmtCpu(e.lo.v, cu) + '); ', runName(e.hi.b.run_id), ' the most (' + fmtCpu(e.hi.v, cu) + ').']]); }
    return items;
  }

  function viewOverview(root) {
    var all = sampleRuns(state.sample), list = visible(state.sample);
    if (!list.length) { root.appendChild(emptyCard('No runs selected. Use the run chips above to show runs.')); return; }
    var base = baseFor(state.sample);
    var b0 = list[0];
    var region = get(b0, 'inputs.regions_intervals');
    var scope = fmtInt(b0.summary.ref_transcripts_total) + ' reference transcripts in scope' +
      (region && region.length ? ' (' + region.slice(0, 3).map(function (r) { return r[0] + ':' + fmtInt(r[1] + 1) + '–' + fmtInt(r[2]); }).join(', ') + (region.length > 3 ? ', …' : '') + ')' : '');
    root.appendChild(sectionHead('How do these runs differ?', list.length + ' of ' + all.length + ' runs on ' + state.sample + ' · ' + scope + '.'));
    if (runOrder.length > 1) root.appendChild(baselineNote());

    var hl = highlights(list);
    if (hl.length) {
      root.appendChild(el('section', { class: 'card' }, el('ul', { class: 'highlights' }, hl.map(function (h) {
        return el('li', null, el('span', { class: 'tag', text: h[0] }), el('span', null, h[1]));
      }))));
    }

    var byId = {};
    list.forEach(function (b) { byId[b.run_id] = b; });
    var rows = ['isoforms', 'fsm', 'novel', 'other', 'ref_cim', 'ref_genes', 'junc_known', 'jaccard', 'cpm_r', 'pro', 'cpu'].map(M)
      .filter(function (m) { return present(m, list); }).map(function (m) {
        var asShare = m.key === 'ref_cim';
        var cu = cpuUnit(list.map(function (b) { return mval(m, b); }));
        return {
          label: m.label, sub: asShare ? '% of reference transcripts in scope' : m.sub,
          fmt: function (v) { return asShare ? fmtPct(v) : m.kind === 'cpu' ? fmtCpu(v, cu) : fmtMetric(m, v, null, false); },
          get: function (id) { var b = byId[id]; if (!b) return null; var v = mval(m, b); return asShare && isNum(v) ? pct(v, b.summary.ref_transcripts_total) : v; }
        };
      });
    root.appendChild(card('Where the runs differ', 'Each row has its own scale. Dots are runs and the ringed dot is the baseline. Hover or focus a row for every value.',
      dotStrip(rows, list.map(function (b) { return b.run_id; }), state.baseline)));

    var keys = ['isoforms', 'genes', 'fsm', 'novel', 'mono', 'ref_fsm', 'ref_cim', 'ref_genes', 'junc_known'].map(M);
    var thead = el('tr', null, el('th', { text: 'Run' }),
      keys.map(function (m) { return el('th', { text: m.short || m.label, title: m.label + (m.sub ? ' (' + m.sub + ')' : '') }); }));
    var tbody = el('tbody');
    list.forEach(function (b) {
      var isBase = base && b.key === base.key;
      var open = function () { state.runKey = b.key; setView('run'); };
      var tr = el('tr', { class: 'clickable', tabindex: '0', title: 'Open ' + b.label, onclick: open, onkeydown: function (ev) { if (ev.key === 'Enter') open(); } });
      tr.appendChild(el('td', null, el('div', { class: 'runcell' }, swatch(runColor(b.run_id)),
        el('div', null, el('span', { class: 'name', text: b.label }), isBase ? el('span', { class: 'pill', text: 'baseline' }) : null,
          el('span', { class: 'meta', text: 'FLAIR ' + versionText(b) + ' · ' + b.run_id, title: 'FLAIR ' + versionText(b, true) })))));
      keys.forEach(function (m) {
        var v = mval(m, b);
        tr.appendChild(el('td', null, fmtMetric(m, v, b), base && !isBase ? el('span', { class: 'delta', text: fmtDelta(m, v, mval(m, base)) }) : null));
      });
      tbody.appendChild(tr);
    });
    root.appendChild(card('All runs', 'The small line under each value is the change from the baseline (pp = percentage points). Select a row to open that run.',
      el('div', { class: 'table-wrap' }, el('table', null, el('thead', null, thead), tbody))));
  }

  // ------------------------------------------------------------------ Structure
  function modeTools() {
    return seg('Show values as', [['share', 'Share'], ['count', 'Count']], state.mode, function (v) { state.mode = v; render(); });
  }
  function viewStructure(root) {
    var list = visible(state.sample);
    root.appendChild(sectionHead('Structure', 'How each run’s isoforms compare with the reference annotation, in SQANTI3 categories.', state.printing ? null : modeTools()));
    if (!list.length) { root.appendChild(emptyCard('No runs selected.')); return; }
    var share = state.mode === 'share';
    var ser = list.map(function (b) { return { label: b.label, color: runColor(b.run_id) }; });
    var valFmt = function (c, t) { return share ? fmtPct(pct(c, t)) + ' (' + fmtInt(c) + ')' : fmtInt(c) + ' (' + fmtPct(pct(c, t)) + ')'; };
    var tickFmt = share ? function (t) { return t + '%'; } : compact;
    var heads = [''].concat(list.map(function (b) { return b.label; }));
    var lg = function () { return runLegend(list.map(function (b) { return b.run_id; })); };
    function groupedCard(title, dek, groups, getC, getT, span) {
      return card(title, dek, [lg(), chartBox(function (box) {
        hbars(box, {
          name: title, series: ser, tickFmt: tickFmt,
          groups: groups.map(function (g) { return { label: g.label, g: g, values: list.map(function (b) { var c = getC(b, g), t = getT(b, g); return isNum(c) ? (share ? pct(c, t) : c) : null; }) }; }),
          fmt: function (v, gg, si) { var b = list[si]; return valFmt(getC(b, gg.g), getT(b, gg.g)); }
        });
      })], { span: span, table: { head: heads, rows: groups.map(function (g) { return [g.long || g.label].concat(list.map(function (b) { var c = getC(b, g); return isNum(c) ? valFmt(c, getT(b, g)) : '—'; })); }) } });
    }
    var shown = CATS.filter(function (c) { return list.some(function (b) { return b.categories[c] > 0; }); });
    var hidden = CATS.filter(function (c) { return shown.indexOf(c) < 0; });
    var cards = [];
    cards.push(groupedCard('Structural categories', (share ? 'Share of each run’s isoforms.' : 'Isoforms per category.') + (hidden.length ? ' Not shown (0 in every run): ' + hidden.map(function (c) { return CAT_LABEL[c]; }).join(', ') + '.' : ''),
      shown.map(function (c) { return { label: CAT_LABEL[c], long: CAT_LONG[c] + ' (' + CAT_LABEL[c] + ')', c: c }; }),
      function (b, g) { return b.categories[g.c]; }, function (b) { return b.summary.isoforms; }, true));
    cards.push(groupedCard('Reference recovery', share ? 'Share of the reference transcripts and genes in scope that each run finds.' : 'Reference transcripts and genes found.',
      [{ label: 'Ref. tx FSM', long: 'Reference transcripts with an FSM', k: 'ref_transcripts_fsm', t: 'ref_transcripts_total' },
        { label: 'Ref. tx CIM', long: 'Reference transcripts with a CIM', k: 'ref_transcripts_cim', t: 'ref_transcripts_total' },
        { label: 'Ref. genes', long: 'Reference genes detected', k: 'ref_genes_detected', t: 'ref_genes_total' }],
      function (b, g) { return b.summary[g.k]; }, function (b, g) { return b.summary[g.t]; }));
    cards.push(groupedCard('FSM end match', (share ? 'Share of each run’s FSM isoforms' : 'FSM isoforms') + ' by whether both ends fall within 50 bp of the matched reference.',
      SUBS.map(function (k) { return { label: k, long: SUB_LABEL[k], k: k }; }),
      function (b, g) { return b.fsm_subcategories[g.k]; }, function (b) { return b.categories.FSM; }));
    var withJ = list.filter(function (b) { return b.junctions; });
    if (withJ.length) {
      var motif = withJ.some(function (b) { return b.junctions.motif_checked !== false; });
      var jgroups = motif ? [{ label: 'Known', k: 'known' }, { label: 'Novel canonical', long: 'Novel canonical (GT-AG, GC-AG, AT-AC)', k: 'novel_canonical' }, { label: 'Novel non-canon.', long: 'Novel non-canonical', k: 'novel_noncanonical' }]
        : [{ label: 'Known', k: 'known' }, { label: 'Novel', k: 'novel' }];
      cards.push(groupedCard('Splice junctions', (share ? 'Share of each run’s distinct junctions' : 'Distinct junctions') + ', by whether the reference has them.', jgroups,
        function (b, g) {
          var J = b.junctions; if (!J) return null;
          if (g.k === 'novel') return isNum(J.novel) ? J.novel : J.total - J.known;
          return J.motif_checked === false && g.k !== 'known' ? null : J[g.k];
        },
        function (b) { return b.junctions ? b.junctions.total : null; }));
    }
    var withP = list.filter(function (b) { return get(b, 'sections.productivity.class_counts'); });
    if (withP.length) {
      var classes = uniq(PROD_ORDER.concat([].concat.apply([], withP.map(function (b) { return Object.keys(b.sections.productivity.class_counts); }))))
        .filter(function (k) { return withP.some(function (b) { return b.sections.productivity.class_counts[k] > 0; }); });
      cards.push(groupedCard('Coding potential', 'FLAIR productivity class of each isoform (' + (share ? 'share of the isoforms classified' : 'counts') + ').',
        classes.map(function (k) { return { label: k, long: PROD_LABEL[k] || k, k: k }; }),
        function (b, g) { var cc = get(b, 'sections.productivity.class_counts'); return cc ? (cc[g.k] || 0) : null; },
        function (b) { var cc = get(b, 'sections.productivity.class_counts'); return cc ? sum(Object.keys(cc).map(function (k) { return cc[k]; })) : null; }));
    }
    root.appendChild(grid(cards));

    var hcards = [];
    HISTS.forEach(function (def) {
      var withH = list.filter(function (b) { var h = get(b, 'histograms.' + def.key); return h && h.counts && h.counts.length; });
      if (!withH.length) return;
      var bins = histBins(withH[0].histograms[def.key], def);
      var totals = withH.map(function (b) { return sum(b.histograms[def.key].counts); });
      var hs = withH.map(function (b) { return { label: b.label, color: runColor(b.run_id) }; });
      hcards.push(card(def.title, null, [legend(hs, 'line'), chartBox(function (box) {
        binChart(box, {
          name: def.title, bins: bins.labels, binTitles: bins.titles, edgeLabels: def.cont, series: hs,
          values: withH.map(function (b, i) { return b.histograms[def.key].counts.map(function (c) { return share ? pct(c, totals[i]) : c; }); }),
          yFmt: share ? function (t) { return t + '%'; } : compact,
          fmt: function (v, si, i) { return valFmt(withH[si].histograms[def.key].counts[i], totals[si]); }
        });
      })], { table: { head: [def.unit].concat(withH.map(function (b) { return b.label; })), rows: bins.titles.map(function (t, i) { return [t].concat(withH.map(function (b, si) { return valFmt(b.histograms[def.key].counts[i], totals[si]); })); }) } }));
    });
    if (hcards.length) {
      root.appendChild(el('div', { style: 'height:28px' }));
      root.appendChild(sectionHead('Distributions', share ? 'Share of each run’s isoforms (or exons, or genes) per bin, so runs of different size compare.' : 'Counts per bin.'));
      root.appendChild(grid(hcards));
    }
  }

  // ------------------------------------------------------------------ Agreement (concordance)
  function copyText(text, btn) {
    var done = function () { var t = btn.textContent; btn.textContent = 'Copied'; setTimeout(function () { btn.textContent = t; }, 1200); };
    var fallback = function () {
      var ta = el('textarea', { style: 'position:fixed;opacity:0' }); ta.value = text; document.body.appendChild(ta); ta.select();
      try { document.execCommand('copy'); done(); } catch (e) { /* ignore */ }
      document.body.removeChild(ta);
    };
    if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(text).then(done, fallback); else fallback();
  }
  function viewAgreement(root) {
    var conc = concFor(state.sample);
    var list = visible(state.sample).filter(function (b) { return conc && conc.runs.indexOf(b.run_id) >= 0; });
    root.appendChild(sectionHead('Agreement', conc ? conc.union + ' distinct isoform structures across ' + conc.runs.length + ' runs on ' + state.sample + '; ' + conc.shared_by_all + ' are found by every run. Multi-exon isoforms match by intron chain, mono-exon isoforms by at least 50% reciprocal overlap on the same strand.' : null));
    if (!conc || list.length < 2) { root.appendChild(emptyCard('Show at least two runs that have isoforms.bed to compare their structures.')); return; }
    var ids = list.map(function (b) { return b.run_id; });
    var others = ids.filter(function (id) { return id !== state.baseline; });
    var cards = [];
    cards.push(card('Structure overlap', 'Jaccard index for each pair of runs: shared structures ÷ structures in either run.', chartBox(function (box) {
      matrix(box, {
        name: 'Jaccard matrix', ids: ids, min: 0, max: 1, fmt: function (v) { return v.toFixed(2); },
        value: function (a, b) { var p = pairOf(conc, a, b); return p ? p.jaccard : null; },
        tip: function (a, b) {
          var p = pairOf(conc, a, b); if (!p) return [];
          var aa = p.a === a, oa = aa ? p.only_a : p.only_b, ob = aa ? p.only_b : p.only_a;
          return [{ value: p.jaccard.toFixed(3), label: 'Jaccard' }, { value: fmtInt(p.shared), label: 'shared' }, { color: runColor(a), value: fmtInt(oa), label: 'only in ' + runLabel(a) }, { color: runColor(b), value: fmtInt(ob), label: 'only in ' + runLabel(b) }];
        }
      });
    }), { table: { head: ['Pair', 'Shared', 'Only first', 'Only second', 'Jaccard'], rows: conc.pairs.filter(function (p) { return ids.indexOf(p.a) >= 0 && ids.indexOf(p.b) >= 0; }).map(function (p) { return [runLabel(p.a) + ' vs ' + runLabel(p.b), fmtInt(p.shared), fmtInt(p.only_a), fmtInt(p.only_b), fmtNum(p.jaccard, 3)]; }) } }));
    if (ids.indexOf(state.baseline) >= 0 && others.length) {
      var rows = others.map(function (id) {
        var p = pairOf(conc, id, state.baseline), aa = p.a === id;
        var mine = aa ? p.only_a : p.only_b, theirs = aa ? p.only_b : p.only_a;
        return { label: runLabel(id), color: runColor(id), left: theirs, right: mine, p: p,
          tip: [{ value: fmtInt(p.shared), label: 'shared with the baseline' }, { color: runColor(id), value: fmtInt(mine), label: 'only in ' + runLabel(id) }, { color: 'var(--neutral-bar)', value: fmtInt(theirs), label: 'only in ' + runLabel(state.baseline) }] };
      });
      cards.push(card('Against the baseline', 'Structures each run adds (right) or misses (left) compared with ' + runLabel(state.baseline) + '.', chartBox(function (box) {
        diverging(box, { name: 'Structures against the baseline', rows: rows, leftLabel: 'only in baseline', rightLabel: 'only in run' });
      }), { table: { head: ['Run', 'Shared', 'Only in run', 'Only in baseline', 'Jaccard'], rows: rows.map(function (r) { return [r.label, fmtInt(r.p.shared), fmtInt(r.right), fmtInt(r.left), fmtNum(r.p.jaccard, 3)]; }) } }));
    }
    root.appendChild(grid(cards));

    var gv = {};
    others.forEach(function (id) { var g = genesVsBase(conc, id, state.baseline); if (g) gv[id] = g; });
    if (ids.indexOf(state.baseline) >= 0) {
      var gids = others.filter(function (id) { return gv[id]; });
      if (gids.length) {
        var grows = gids.map(function (id) {
          return { label: runLabel(id), color: runColor(id), left: gv[id].lost, right: gv[id].gained,
            tip: [{ color: runColor(id), value: fmtInt(gv[id].gained), label: 'gained' }, { color: 'var(--neutral-bar)', value: fmtInt(gv[id].lost), label: 'lost' }] };
        });
        var lists = el('div', { style: 'margin-top:8px' }, gids.map(function (id) {
          var g = gv[id];
          var names = function (arr, n) { return arr.length ? arr.join(', ') + (n > arr.length ? ' … and ' + (n - arr.length) + ' more' : '') : 'none'; };
          return el('details', null, el('summary', null, runLabel(id) + ': ' + g.gained + ' gained, ' + g.lost + ' lost'),
            el('p', { class: 'names' }, el('strong', { text: 'Gained: ' }), names(g.gained_names, g.gained)),
            el('p', { class: 'names' }, el('strong', { text: 'Lost: ' }), names(g.lost_names, g.lost)));
        }));
        root.appendChild(card('Reference genes gained and lost', 'Reference genes a run detects (FSM, ISM, NIC, NNC or fusion) that the baseline does not (right), and the reverse (left). Open a run below for the gene names.',
          [chartBox(function (box) { diverging(box, { name: 'Genes against the baseline', rows: grows, leftLabel: 'lost', rightLabel: 'gained' }); }), lists]));
      }
    }

    var loci = conc.top_loci || [];
    if (loci.length) {
      var cidx = ids.map(function (id) { return conc.runs.indexOf(id); });
      var head = el('tr', null, el('th', { text: 'Gene' }), el('th', { class: 'left', text: 'Locus (1-based)' }), el('th', { text: 'Structures', title: 'Distinct structures in the locus (found by every run)' }),
        ids.map(function (id) { return el('th', { title: runLabel(id) }, swatch(runColor(id)), ' ', runLabel(id)); }));
      var body = el('tbody', null, loci.map(function (l) {
        var coord = l.chrom + ':' + (l.start + 1) + '-' + l.end;
        var counts = cidx.map(function (i) { return l.counts[i]; });
        var mx = Math.max.apply(null, counts);
        var btn = el('button', { class: 'copy', type: 'button', text: 'Copy', 'aria-label': 'Copy ' + coord });
        btn.addEventListener('click', function () { copyText(coord, btn); });
        var gl = l.genes && l.genes.length ? l.genes.slice(0, 2).join(', ') + (l.genes.length > 2 ? ' +' + (l.genes.length - 2) : '') : 'novel locus';
        return el('tr', null, el('td', { class: 'strong', text: gl, title: (l.genes || []).join(', ') }),
          el('td', { class: 'left' }, el('span', { class: 'mono', text: coord + ' ' + l.strand }), btn),
          el('td', { text: l.structures + ' (' + l.shared_by_all + ')' }),
          counts.map(function (c) { return el('td', { class: c === mx && mx > 0 ? 'strong' : null, text: fmtInt(c) }); }));
      }));
      root.appendChild(el('div', { style: 'height:20px' }));
      root.appendChild(card('Loci where the runs differ most', 'Top ' + loci.length + ' loci, ranked by structures not found in every run, then by the spread of isoform counts. Run columns give isoforms per run; the largest is bold. Structures: distinct (shared by every run). Copy puts the coordinates on the clipboard for IGV or a genome browser.',
        el('div', { class: 'table-wrap' }, el('table', null, el('thead', null, head), body))));
    }
  }

  // ------------------------------------------------------------------ Expression
  function viewExpression(root) {
    var list = visible(state.sample);
    var quant = quantFor(state.sample);
    var withQ = list.filter(function (b) { return get(b, 'sections.quant'); });
    root.appendChild(sectionHead('Expression', 'FLAIR quantify results. Runs are compared on the isoforms they share (same structure), with CPM summed per structure.'));
    if (!withQ.length && !quant) { root.appendChild(emptyCard('No quantification for these runs.')); return; }
    var qn = quantNote(list);
    if (qn) root.appendChild(noteBox(qn));
    var noQ = list.filter(function (b) { return !get(b, 'sections.quant'); });
    if (noQ.length) {
      var baseNoQ = noQ.some(function (b) { return b.run_id === state.baseline; });
      root.appendChild(el('p', { class: 'basenote' }, 'No quantification for ' + noQ.map(function (b) { return b.label; }).join(', ') +
        ' (FLAIR quantify gave no counts; see that run’s provenance).' + (baseNoQ ? ' The baseline has none, so the per-run plots against it are left out; pick another baseline to see them.' : '')));
    }
    var cards = [];
    if (withQ.length) {
      cards.push(card('Reads assigned', 'Reads that FLAIR quantify assigned to each run’s isoforms.', chartBox(function (box) {
        hbars(box, {
          name: 'Reads assigned', series: [{ label: 'Reads', color: 'var(--series-other)' }], labelValues: true, labelW: 200, labelFrac: 0.4,
          groups: withQ.map(function (b) { return { label: b.label, b: b, values: [b.sections.quant.reads_assigned] }; }),
          colorOf: function (gi) { return runColor(withQ[gi].run_id); }, dotOf: function (gi) { return runColor(withQ[gi].run_id); },
          fmt: function (v, g) { return fmtInt(v) + ' reads · ' + fmtInt(g.b.sections.quant.isoforms_quantified) + ' isoforms with reads' + (quantInput(g.b) ? ' · counted ' + QUANT_INPUT[quantInput(g.b)] : ''); },
          short: function (v) { return compact(v); }
        });
      }), { table: { head: ['Run', 'Reads assigned', 'Isoforms with reads', 'Counted'], rows: withQ.map(function (b) { return [b.label, fmtInt(b.sections.quant.reads_assigned), fmtInt(b.sections.quant.isoforms_quantified), quantInput(b) ? QUANT_INPUT[quantInput(b)] : 'not recorded']; }) } }));
    }
    var qids = quant ? list.map(function (b) { return b.run_id; }).filter(function (id) { return quant.runs.indexOf(id) >= 0; }) : [];
    if (quant && qids.length >= 2) {
      var vals = quant.pairs.map(function (p) { return p.pearson; }).filter(isNum);
      var lo = Math.min(0.8, Math.floor(Math.min.apply(null, vals.concat([1])) * 10) / 10);
      cards.push(card('CPM correlation', 'Pearson r of log10(CPM + 1) over shared structures, for each pair of runs. The colour scale starts at ' + lo.toFixed(1) + '.', chartBox(function (box) {
        matrix(box, {
          name: 'CPM correlation', ids: qids, min: lo, max: 1, fmt: function (v) { return v.toFixed(2); },
          value: function (a, b) { var p = pairOf(quant, a, b); return p ? p.pearson : null; },
          tip: function (a, b) { var p = pairOf(quant, a, b); return p ? [{ value: fmtNum(p.pearson, 3), label: 'Pearson r (log CPM)' }, { value: fmtNum(p.spearman, 3), label: 'Spearman ρ' }, { value: fmtInt(p.n), label: 'shared structures' }] : []; }
        });
      }), { table: { head: ['Pair', 'Shared', 'Pearson r', 'Spearman ρ', 'CPM in shared (first)', 'CPM in shared (second)'], rows: quant.pairs.filter(function (p) { return qids.indexOf(p.a) >= 0 && qids.indexOf(p.b) >= 0; }).map(function (p) { return [runLabel(p.a) + ' vs ' + runLabel(p.b), fmtInt(p.n), fmtNum(p.pearson, 3), fmtNum(p.spearman, 3), fmtPct(100 * p.cpm_shared_a), fmtPct(100 * p.cpm_shared_b)]; }) } }));
    }
    if (cards.length) root.appendChild(grid(cards));

    if (quant && qids.indexOf(state.baseline) >= 0) {
      var dcards = qids.filter(function (id) { return id !== state.baseline; }).map(function (id) {
        var p = pairOf(quant, id, state.baseline);
        var g = gridVsBase(p, state.baseline);
        if (!g) return null;
        var mine = p.a === id ? p.cpm_shared_a : p.cpm_shared_b;
        return el('section', { class: 'card' },
          el('div', { class: 'card-head' }, el('div', null, el('h3', null, swatch(runColor(id)), ' ', runLabel(id)),
            el('p', { class: 'dek', text: 'r = ' + fmtNum(p.pearson) + ' over ' + fmtInt(p.n) + ' shared isoforms, which hold ' + fmtPct(100 * mine, 0) + ' of its CPM.' }))),
          chartBox(function (box) { density(box, { name: runLabel(id) + ' against baseline CPM', grid: g, logMax: quant.log_max }); }),
          el('p', { class: 'small', text: 'Up: ' + runLabel(id) + ' CPM. Across: ' + runLabel(state.baseline) + ' CPM.' }));
      }).filter(Boolean);
      if (dcards.length) {
        root.appendChild(el('div', { style: 'height:28px' }));
        root.appendChild(sectionHead('Against the baseline', 'Shared isoforms, with run CPM (up) against baseline CPM (across). Both axes are log10(CPM + 1). Darker cells hold more isoforms; the diagonal is perfect agreement.'));
        root.appendChild(grid(dcards));
      }
    }

    var withH = withQ.filter(function (b) { return get(b, 'sections.quant.cpm_hist.counts'); });
    if (withH.length) {
      var bins = histBins(withH[0].sections.quant.cpm_hist, { cont: true, unit: 'CPM' });
      var totals = withH.map(function (b) { return sum(b.sections.quant.cpm_hist.counts); });
      var hs = withH.map(function (b) { return { label: b.label, color: runColor(b.run_id) }; });
      root.appendChild(el('div', { style: 'height:20px' }));
      root.appendChild(card('CPM distribution', 'Share of each run’s isoforms with reads, per CPM bin.', [legend(hs, 'line'), chartBox(function (box) {
        binChart(box, {
          name: 'CPM distribution', bins: bins.labels, binTitles: bins.titles, edgeLabels: true, series: hs,
          values: withH.map(function (b, i) { return b.sections.quant.cpm_hist.counts.map(function (c) { return pct(c, totals[i]); }); }),
          yFmt: function (t) { return t + '%'; },
          fmt: function (v, si, i) { return fmtPct(v) + ' (' + fmtInt(withH[si].sections.quant.cpm_hist.counts[i]) + ')'; }
        });
      })], { table: { head: ['CPM'].concat(withH.map(function (b) { return b.label; })), rows: bins.titles.map(function (t, i) { return [t].concat(withH.map(function (b) { return fmtInt(b.sections.quant.cpm_hist.counts[i]); })); }) } }));
    }
  }

  // ------------------------------------------------------------------ Runtime
  function viewRuntime(root) {
    var list = visible(state.sample).filter(function (b) { return ['cpu', 'wall', 'rss'].some(function (k) { return mval(M(k), b) !== null; }); });
    root.appendChild(sectionHead('Runtime', 'Compute recorded by the pipeline for each job’s FLAIR tasks (assemble, quantify, productivity).'));
    if (!list.length) { root.appendChild(emptyCard('No compute records for these runs.')); return; }
    var cards = ['cpu', 'wall', 'rss'].map(M).filter(function (m) { return present(m, list); }).map(function (m) {
      var lw = list.filter(function (b) { return mval(m, b) !== null; });
      var raw = lw.map(function (b) { return mval(m, b); });
      // one display unit per chart
      var u = m.kind === 'cpu' ? (cpuUnit(raw) === 'min' ? { f: 60, name: 'CPU-minutes' } : { f: 1, name: 'CPU-hours' })
        : m.kind === 'wall' ? (Math.max.apply(null, raw) >= 120 ? { f: 1 / 60, name: 'hours' } : { f: 1, name: 'minutes' })
          : (Math.max.apply(null, raw) < 1 ? { f: 1024, name: 'MB' } : { f: 1, name: 'GB' });
      var show = function (v) { var x = v * u.f; return x >= 100 ? fmtInt(x) : x >= 10 ? x.toFixed(1) : x.toFixed(2); };
      return card(m.label + ' (' + u.name + ')', null, chartBox(function (box) {
        hbars(box, {
          name: m.label, series: [{ label: m.label, color: 'var(--series-other)' }], labelValues: true, labelW: 190, labelFrac: 0.46,
          groups: lw.map(function (b) { return { label: b.label, values: [mval(m, b) * u.f] }; }),
          colorOf: function (gi) { return runColor(lw[gi].run_id); }, dotOf: function (gi) { return runColor(lw[gi].run_id); },
          fmt: function (v) { return show(v / u.f) + ' ' + u.name; }, short: function (v) { return show(v / u.f); }, tickFmt: function (t) { return String(+t.toFixed(3)); }
        });
      }), { table: { head: ['Run', m.label + ' (' + u.name + ')'], rows: lw.map(function (b) { return [b.label, show(mval(m, b))]; }) } });
    });
    root.appendChild(grid(cards, 'three'));
  }

  // ------------------------------------------------------------------ Run detail
  // Tags worth showing: those that do not just repeat the FLAIR version or source.
  function extraTags(b) {
    return Object.keys(b.tags || {}).filter(function (k) {
      var v = b.tags[k];
      return !((k === 'version' && v === get(b, 'flair.version')) || (k === 'source' && v === get(b, 'flair.source')));
    });
  }
  function tile(label, value, detail) {
    return el('div', { class: 'tile' }, el('div', { class: 'label', text: label }), el('div', { class: 'value', text: value }), detail ? el('div', { class: 'detail', text: detail }) : null);
  }
  function provenance(b) {
    var f = b.flair || {};
    var items = [['Sample', b.sample], ['Run id', b.run_id], ['FLAIR', (f.version || '?') + ' (' + f.source + ')'], ['Requested', f.requested || '—'],
      ['Quantify input', get(b, 'sections.quant.input') ? QUANT_INPUT[b.sections.quant.input] + (get(b, 'sections.quant.input_inferred') ? ' (from the adapter log)' : '') : '—'],
      ['Commit', f.commit || '—'], ['Settings', get(b, 'inputs.transcriptome_args') || '—'], ['Reference', get(b, 'inputs.gtf') || '—'],
      ['Regions', get(b, 'inputs.regions') || 'whole genome'], ['Built', get(b, 'pipeline.created') || '—'],
      ['flair-eval', (get(b, 'pipeline.version') || '') + (get(b, 'pipeline.commit') ? ' @ ' + String(b.pipeline.commit).slice(0, 10) : '')]];
    return el('dl', { class: 'kv' }, items.map(function (kv) { return [el('dt', { text: kv[0] }), el('dd', { class: kv[0] === 'Commit' ? 'mono' : null, text: kv[1] })]; }));
  }
  function viewRun(root) {
    var list = sampleRuns(state.sample);
    if (!list.length) { root.appendChild(emptyCard('No runs match these filters.')); return; }
    var b = list.find(function (x) { return x.key === state.runKey; }) || list.find(function (x) { return x.run_id === state.baseline; }) || list[0];
    state.runKey = b.key;
    var S = b.summary;
    var picker = state.printing ? null : selectField('Run', b.key, list.map(function (x) { return [x.key, x.label]; }), function (v) { state.runKey = v; render(); });
    root.appendChild(el('div', { class: 'section-head' },
      el('div', null, el('h2', null, swatch(runColor(b.run_id)), ' ', b.label),
        el('p', { class: 'dek', text: 'FLAIR ' + versionText(b, true) + ' · ' + b.run_id + ' on ' + b.sample + (extraTags(b).length ? ' · ' + extraTags(b).map(function (k) { return k + ': ' + b.tags[k]; }).join(' · ') : '') })),
      picker));
    var tiles = [tile('Isoforms', fmtInt(S.isoforms), fmtInt(S.mono_exon) + ' mono-exon'), tile('Genes', fmtInt(S.genes), fmtInt(S.loci) + ' loci'),
      tile('Ref. transcripts, FSM', fmtInt(S.ref_transcripts_fsm), 'of ' + fmtInt(S.ref_transcripts_total) + ' (' + fmtPct(pct(S.ref_transcripts_fsm, S.ref_transcripts_total)) + ')'),
      tile('Ref. transcripts, CIM', fmtInt(S.ref_transcripts_cim), 'both ends within 50 bp'),
      tile('Ref. genes detected', fmtInt(S.ref_genes_detected), S.ref_genes_total ? 'of ' + fmtInt(S.ref_genes_total) : null)];
    var j = jaccardVsBase(b);
    if (isNum(j)) tiles.push(tile('Shared with baseline', j.toFixed(2), 'Jaccard of structures'));
    if (get(b, 'sections.quant')) tiles.push(tile('Reads assigned', fmtInt(b.sections.quant.reads_assigned), fmtInt(b.sections.quant.isoforms_quantified) + ' isoforms with reads' + (quantInput(b) ? ', counted ' + QUANT_INPUT[quantInput(b)] : '')));
    if (isNum(get(b, 'compute.cpu_hours'))) tiles.push(tile('CPU time', fmtCpu(b.compute.cpu_hours, cpuUnit([b.compute.cpu_hours])), [isNum(get(b, 'compute.wall_minutes')) ? fmtWall(b.compute.wall_minutes) + ' wall' : null, isNum(get(b, 'compute.peak_rss_gb')) ? fmtMem(b.compute.peak_rss_gb) + ' peak memory' : null].filter(Boolean).join(' · ') || null));
    root.appendChild(el('div', { class: 'tiles' }, tiles));
    root.appendChild(el('div', { style: 'height:20px' }));

    var color = runColor(b.run_id);
    var cards = [];
    cards.push(card('Structural categories', 'Isoforms per SQANTI3 category.', chartBox(function (box) {
      hbars(box, { name: 'Structural categories', series: [{ label: b.label, color: color }], labelValues: true,
        groups: CATS.map(function (c) { return { label: CAT_LABEL[c], values: [b.categories[c]] }; }),
        fmt: function (v) { return fmtInt(v) + ' · ' + fmtPct(pct(v, S.isoforms)); } });
    }), { table: { head: ['Category', 'Isoforms', 'Share'], rows: CATS.map(function (c) { return [CAT_LONG[c] + ' (' + CAT_LABEL[c] + ')', fmtInt(b.categories[c]), fmtPct(pct(b.categories[c], S.isoforms))]; }) } }));
    var parts = el('section', { class: 'card' });
    function block(title, dek, ps, name) {
      parts.appendChild(el('h3', { text: title, style: parts.childNodes.length ? 'margin-top:18px' : null }));
      if (dek) parts.appendChild(el('p', { class: 'dek', text: dek, style: 'font-size:13px;margin-bottom:8px' }));
      parts.appendChild(legend(ps));
      parts.appendChild(segbar(ps, { name: name }));
    }
    block('FSM end match', fmtInt(b.categories.FSM) + ' FSM isoforms by end match (50 bp).', SUBS.map(function (k, i) { return { label: SUB_LABEL[k], value: b.fsm_subcategories[k], color: PART(i), on: ON_PART(i) }; }), 'FSM end match');
    if (b.junctions) {
      var J = b.junctions;
      var jp = (J.motif_checked === false ? [{ label: 'Known', value: J.known }, { label: 'Novel', value: isNum(J.novel) ? J.novel : J.total - J.known }]
        : [{ label: 'Known', value: J.known }, { label: 'Novel canonical', value: J.novel_canonical }, { label: 'Novel non-canonical', value: J.novel_noncanonical }])
        .map(function (p, i) { p.color = PART(i); p.on = ON_PART(i); return p; });
      block('Splice junctions', fmtInt(J.total) + ' distinct junctions.', jp, 'Junctions');
    }
    var cc = get(b, 'sections.productivity.class_counts');
    if (cc) {
      var ks = uniq(PROD_ORDER.concat(Object.keys(cc))).filter(function (k) { return cc[k] > 0; });
      block('Coding potential', 'FLAIR productivity classes.', ks.map(function (k, i) { return { label: PROD_LABEL[k] || k, value: cc[k], color: PART(i), on: ON_PART(i) }; }), 'Productivity');
    }
    cards.push(parts);
    root.appendChild(grid(cards));

    var hcards = [];
    HISTS.forEach(function (def) {
      var h = get(b, 'histograms.' + def.key);
      if (!h || !h.counts || !h.counts.length) return;
      var bins = histBins(h, def), total = sum(h.counts);
      hcards.push(card(def.title, fmtInt(total) + ' ' + def.of + '.', chartBox(function (box) {
        binChart(box, { name: def.title, bins: bins.labels, binTitles: bins.titles, edgeLabels: def.cont, series: [{ label: b.label, color: color }], values: [h.counts], yFmt: compact,
          fmt: function (v) { return fmtInt(v) + ' ' + def.of + ' (' + fmtPct(pct(v, total)) + ')'; } });
      }), { table: { head: [def.unit, def.of, 'Share'], rows: h.counts.map(function (c, i) { return [bins.titles[i], fmtInt(c), fmtPct(pct(c, total))]; }) } }));
    });
    if (hcards.length) root.appendChild(grid(hcards));
    Object.keys(SECTION_PANELS).forEach(function (k) { if (b.sections && b.sections[k] && SECTION_PANELS[k].renderRun) root.appendChild(SECTION_PANELS[k].renderRun(b)); });
    root.appendChild(el('section', { class: 'card', style: 'margin-top:20px' }, el('details', state.printing ? { open: true } : null, el('summary', { text: 'Provenance' }), provenance(b))));
  }

  // ------------------------------------------------------------------ Guide
  function viewGuide(root) {
    var c = el('section', { class: 'card guide' });
    var p = function (t) { return el('p', { text: t }); };
    var dl = function (pairs) { return el('dl', null, pairs.map(function (x) { return [el('dt', { text: x[0] }), el('dd', { text: x[1] })]; })); };
    c.appendChild(el('h2', { text: 'How to read this report' }));
    c.appendChild(p('Each run is one FLAIR version, branch or setting. All runs on a sample use the same reads, alignment and reference annotation, so the differences come from FLAIR. Each run keeps one colour throughout, and the baseline is the run that every difference is measured against.'));
    c.appendChild(el('h3', { text: 'Views' }));
    c.appendChild(dl([
      ['Overview', 'The headline differences, one dot row per metric, and a run table with changes from the baseline.'],
      ['Structure', 'SQANTI3 categories, reference recovery, FSM end match, junctions, coding potential and length distributions. "Share" divides by each run’s own total, so a bigger run does not look better just for being bigger.'],
      ['Agreement', 'Which isoform structures the runs share, which reference genes they gain or lose against the baseline, and the loci where they disagree most.'],
      ['Expression', 'FLAIR quantify: reads assigned, CPM agreement on shared isoforms, and CPM distributions.'],
      ['Runtime', 'CPU time, wall time and peak memory for each job’s FLAIR tasks.'],
      ['Run detail', 'Every panel for a single run, with its provenance.']
    ]));
    c.appendChild(el('h3', { text: 'Structural categories (SQANTI3 names)' }));
    c.appendChild(dl([
      ['FSM', 'Full splice match: the same intron chain as a reference transcript. For mono-exon isoforms: overlaps a mono-exon reference transcript on the same strand.'],
      ['ISM', 'Incomplete splice match: a consecutive part of a reference intron chain, or a mono-exon isoform inside one reference exon.'],
      ['NIC', 'Novel in catalog: every splice site is annotated but the combination is new, or an intron is retained.'],
      ['NNC', 'Novel not in catalog: at least one splice site the reference lacks, inside a known spliced gene.'],
      ['Genic', 'Overlaps exons of a known gene without matching its splice sites (for example a mono-exon gene).'],
      ['Genic intron', 'Lies inside an intron of a known gene.'],
      ['Antisense', 'Overlaps exons of a known gene only on the opposite strand.'],
      ['Intergenic', 'Overlaps no reference transcript.'],
      ['Fusion', 'Shares splice sites with two or more separate genes.']
    ]));
    c.appendChild(el('h3', { text: 'FSM end match' }));
    c.appendChild(dl([['CIM', 'Both ends within 50 bp of the matched reference transcript (SQANTI3 reference_match).'], ['FSM5', "5' end more than 50 bp away; 3' end matches."], ['FSM3', "3' end more than 50 bp away; 5' end matches."], ['FSMTE', 'Both ends more than 50 bp away.']]));
    c.appendChild(el('h3', { text: 'Other measures' }));
    c.appendChild(dl([
      ['Reference recovery', 'Reference transcripts with at least one FSM (or CIM) isoform, and reference genes with at least one FSM, ISM, NIC, NNC or fusion isoform. "In scope" means inside the evaluated regions.'],
      ['Junctions', 'Distinct splice junctions. Known means the reference has it. Novel junctions are canonical when the intron is GT-AG, GC-AG or AT-AC.'],
      ['Coding potential', 'FLAIR predictProductivity: PRO productive, PTC premature stop codon (likely NMD), NST start but no stop, NGO no annotated start.'],
      ['Structures', 'For agreement, multi-exon isoforms are the same structure when their intron chains match, and mono-exon isoforms when they overlap by at least half of both lengths on the same strand. Isoforms that differ only in their ends count once.'],
      ['Jaccard', 'Shared structures ÷ structures in either run. 1 means identical sets.'],
      ['CPM correlation', 'Pearson r of log10(CPM + 1), summed per structure, over structures both runs have.'],
      ['Histograms', 'Fixed bins, so runs compare bin by bin. The last bin is open-ended.']
    ]));
    c.appendChild(el('h3', { text: 'Using the report' }));
    c.appendChild(p('The run chips under the title show or hide runs in every view. Every chart has a Table button with its numbers. Hover a mark, or reach it with the Tab key, to read its value. Print / PDF prints every view in turn.'));
    root.appendChild(c);
  }

  // ------------------------------------------------------------------ shell
  var VIEWS = [
    { id: 'overview', label: 'Overview', fn: viewOverview },
    { id: 'structure', label: 'Structure', fn: viewStructure },
    { id: 'agreement', label: 'Agreement', fn: viewAgreement, when: function () { return Object.keys(bySample).some(function (k) { return bySample[k].concordance; }); } },
    { id: 'expression', label: 'Expression', fn: viewExpression, when: function () { return bundles.some(function (b) { return get(b, 'sections.quant'); }) || Object.keys(bySample).some(function (k) { return bySample[k].quant; }); } },
    { id: 'runtime', label: 'Runtime', fn: viewRuntime, when: function () { return bundles.some(function (b) { return ['cpu_hours', 'wall_minutes', 'peak_rss_gb'].some(function (k) { return isNum(get(b, 'compute.' + k)); }); }); } },
    { id: 'run', label: 'Run detail', fn: viewRun },
    { id: 'guide', label: 'Guide', fn: viewGuide }
  ].filter(function (v) { return !v.when || v.when(); });
  function viewById(id) { return VIEWS.find(function (v) { return v.id === id; }); }

  function renderTabs() {
    var nav = document.getElementById('tabs');
    nav.textContent = '';
    VIEWS.forEach(function (v, i) {
      var on = v.id === state.view;
      var b = el('button', { role: 'tab', type: 'button', id: 'tab-' + v.id, 'aria-selected': on ? 'true' : 'false', tabindex: on ? '0' : '-1', 'aria-controls': 'view', text: v.label });
      b.addEventListener('click', function () { setView(v.id); });
      b.addEventListener('keydown', function (e) {
        var k = e.key === 'ArrowRight' ? 1 : e.key === 'ArrowLeft' ? -1 : e.key === 'Home' ? -i : e.key === 'End' ? VIEWS.length - 1 - i : 0;
        if (!k) return;
        e.preventDefault();
        var nv = VIEWS[(i + k + VIEWS.length) % VIEWS.length];
        setView(nv.id);
        document.getElementById('tab-' + nv.id).focus();
      });
      nav.appendChild(b);
    });
  }
  function setView(id) {
    state.view = viewById(id) ? id : 'overview';
    try { history.replaceState(null, '', '#' + state.view); } catch (e) { /* file:// */ }
    render();
  }
  var lastWidth = 0;
  function render() {
    hideTip();
    setSubtitle();
    renderRunbar();
    renderTabs();
    var root = document.getElementById('view');
    root.setAttribute('role', 'tabpanel');
    root.setAttribute('aria-labelledby', 'tab-' + state.view);
    root.textContent = '';
    pending = [];
    if (state.printing) {
      VIEWS.forEach(function (v) {
        if (v.id === 'guide') return;
        var sec = el('div', { class: 'print-section' });
        root.appendChild(sec);
        v.fn(sec);
      });
    } else {
      viewById(state.view).fn(root);
    }
    flush();
    lastWidth = root.clientWidth;
  }

  function initTheme() {
    var box = document.getElementById('theme'), cur = 'auto';
    try { cur = localStorage.getItem('flair-eval-theme') || 'auto'; } catch (e) { cur = 'auto'; }
    function apply() {
      if (cur === 'auto') document.documentElement.removeAttribute('data-theme'); else document.documentElement.setAttribute('data-theme', cur);
      box.textContent = '';
      [['auto', 'Auto'], ['light', 'Light'], ['dark', 'Dark']].forEach(function (o) {
        box.appendChild(el('button', { type: 'button', 'aria-pressed': o[0] === cur ? 'true' : 'false', text: o[1], onclick: function () {
          cur = o[0]; try { localStorage.setItem('flair-eval-theme', cur); } catch (e) { /* storage blocked */ } apply();
        } }));
      });
    }
    apply();
  }

  function init() {
    document.getElementById('footer').textContent = 'Made with flair-eval ' + (DATA.tool_version || '') + ' from ' + bundles.length + (bundles.length === 1 ? ' bundle' : ' bundles') + '. This file is self-contained and works offline.';
    initTheme();
    document.getElementById('print').addEventListener('click', function () {
      state.printing = true; render();
      setTimeout(function () { window.print(); }, 150);
    });
    window.addEventListener('afterprint', function () { if (state.printing) { state.printing = false; render(); } });
    var h = (location.hash || '').replace('#', '');
    state.view = viewById(h) ? h : (bundles.length === 1 ? 'run' : 'overview');
    if (!runInfo[state.baseline]) state.baseline = runOrder[0];
    render();
    var t = null;
    window.addEventListener('resize', function () {
      clearTimeout(t);
      t = setTimeout(function () { if (Math.abs(document.getElementById('view').clientWidth - lastWidth) > 8) render(); }, 150);
    });
  }
  init();
})();
