/* flair-eval config builder: page wiring. All rules live in logic.js. */
(function () {
  "use strict";
  var FE = window.FlairEval;
  var STORE_KEY = "flair-eval-builder-v2"; // v2: phoenix defaults, genome presets, required outdir
  var RELEASES = ["3.0.1", "2.2.0"];
  var SOURCES = [["release", "Release"], ["git", "Git branch"], ["env", "Conda env"], ["container", "Container"]];
  var SAMPLE_LABELS = { sample: "sample name", reads: "reads or BAM file", bam: "BAM file", genome: "genome FASTA",
    gtf: "annotation GTF", platform: "platform", library: "library", cage: "CAGE", pas: "PAS", junctions: "junctions",
    spikein_gtf: "spike-in GTF", outdir: "output folder" };
  var GENOME_KINDS = Object.keys(FE.GENOME_PRESETS).map(function (k) { return [k, k + " (lab preset)"]; })
    .concat([["custom", "Custom FASTA"]]);
  var uidN = 0;
  function uid() { uidN += 1; return "u" + uidN; }

  /* ------------------------------------------------------------ state */
  function blankSample(prev) {
    prev = prev || {};
    return { uid: uid(), sample: "", nameAuto: true, kind: "reads", path: "",
      genomeKind: prev.genomeKind || "hg38", genome: prev.genome || "", gtf: prev.gtf || "",
      platform: prev.platform || "", library: prev.library || "",
      cage: "", pas: "", junctions: "", spikein_gtf: "" };
  }
  function blankRun(source) {
    return { uid: uid(), id: "", idAuto: true, source: source || "release", release: source === "release" || !source ? "3.0.1" : "",
      git_url: FE.DEFAULT_GIT_URL, ref: source === "git" ? "dev" : "", env: "", container: "",
      preset: "default", args: "", label: "", tags: "", quantify: true, productivity: true, own_align: false };
  }
  function defaultOpts() {
    return { outdir: "", executor: "phoenix", software: "conda", profile: "", ref: "", report_title: "", baseline: "", align_flair: "", regions: "",
      bundles: "", env_cache: "", resume: false, prefer_files: false };
  }
  function exampleState() {
    var s = blankSample();
    s.path = "/data/HG002_cdna.fq.gz"; s.genomeKind = "hg38";
    s.platform = "pacbio"; s.library = "cDNA";
    var a = blankRun("release");
    var b = blankRun("git");
    var opts = defaultOpts();
    opts.outdir = "/private/groups/brookslab/your-user/flair-eval/10_11_26_hg002";
    return { samples: [s], runs: [a, b], opts: opts, example: true, tab: "samples", touched: {} };
  }
  function blankState() {
    return { samples: [blankSample()], runs: [blankRun("release")], opts: defaultOpts(), example: false, tab: "samples", touched: {} };
  }

  var state = load() || exampleState();

  function load() {
    try {
      var raw = window.localStorage.getItem(STORE_KEY);
      if (!raw) return null;
      var s = JSON.parse(raw);
      if (!s || !Array.isArray(s.samples) || !Array.isArray(s.runs)) return null;
      s.samples.forEach(function (x) { x.uid = uid(); });
      s.runs.forEach(function (x) { x.uid = uid(); });
      s.opts = Object.assign(defaultOpts(), s.opts || {});
      s.touched = {};
      return s;
    } catch (e) { return null; }
  }
  function save() {
    try { window.localStorage.setItem(STORE_KEY, JSON.stringify(state)); } catch (e) { /* storage blocked */ }
  }

  /* ------------------------------------------------------------ derived data */
  function sampleName(s) {
    if (!s.nameAuto) return s.sample;
    var n = s.path ? FE.defaultSampleName(s.path) : "";
    return FE.ID_RE.test(n) ? n : n.replace(/[^A-Za-z0-9._-]/g, "-").replace(/^[^A-Za-z0-9]+/, "");
  }
  function toSample(s) {
    var row = { sample: sampleName(s), genome: s.genomeKind === "custom" ? s.genome : s.genomeKind, gtf: s.gtf, platform: s.platform, library: s.library,
      cage: s.cage, pas: s.pas, junctions: s.junctions, spikein_gtf: s.spikein_gtf };
    row[s.kind] = s.path;
    var clean = FE.cleanSample(row), out = {};
    FE.ALL_SAMPLE_COLUMNS.forEach(function (c) { if (c in clean) out[c] = clean[c]; });
    return out;
  }
  function flairOf(r) {
    if (r.source === "release") return { release: r.release.trim() };
    if (r.source === "git") return { git: r.git_url.trim() || FE.DEFAULT_GIT_URL, ref: r.ref.trim() };
    if (r.source === "env") return { env: r.env.trim() };
    return { container: r.container.trim() };
  }
  function parseTags(text) {
    var tags = {}, bad = [];
    String(text || "").split(",").forEach(function (part) {
      var p = part.trim();
      if (!p) return;
      var i = p.indexOf("=");
      if (i <= 0) { bad.push(p); return; }
      tags[p.slice(0, i).trim()] = p.slice(i + 1).trim();
    });
    return { tags: tags, bad: bad };
  }
  function autoId(r, taken) {
    var tok = FE.tokenForFlair(flairOf(r)), base = null;
    if (tok) { try { base = FE.parseToken(tok).id; } catch (e) { base = null; } }
    if (!base) base = r.source === "git" ? "flair-branch" : r.source === "env" ? "env-run" : r.source === "container" ? "ctr-run" : "flair-run";
    if (taken.indexOf(base) < 0) return base;
    var key = r.preset !== "default" ? r.preset.replace(/_/g, "-") : "2";
    var id = base + "-" + key, n = 2;
    while (taken.indexOf(id) >= 0) { n += 1; id = base + "-" + key + "-" + n; }
    return id;
  }
  function runIds() {
    var taken = [];
    return state.runs.map(function (r) {
      var id = r.idAuto ? autoId(r, taken) : r.id.trim();
      taken.push(id);
      return id;
    });
  }
  function toRuns() {
    var ids = runIds();
    return state.runs.map(function (r, i) {
      var run = { id: ids[i], flair: flairOf(r), transcriptome_args: r.args.trim(),
        quantify: !!r.quantify, productivity: !!r.productivity, own_align: !!r.own_align };
      if (!run.id) delete run.id;
      if (r.label.trim()) run.label = r.label.trim();
      var t = parseTags(r.tags).tags;
      if (Object.keys(t).length) run.tags = t;
      return run;
    });
  }
  function cmdOptions(runs) {
    var o = state.opts;
    var ids = runs.map(function (r) { return r.id; });
    var baseline = o.baseline && ids.indexOf(o.baseline) > 0 ? o.baseline : null; // first run is the default
    // phoenix is the pipeline's default profile, so it needs no -profile at all
    var profile = o.profile.trim() || (o.executor === "local" ? ["local", o.software].filter(Boolean).join(",") : "");
    return { ref: o.ref.trim(), profile: profile, outdir: o.outdir.trim(), resume: o.resume,
      prefer_files: o.prefer_files,
      params: { baseline: baseline, report_title: o.report_title.trim(), align_flair: o.align_flair.trim(),
        regions: o.regions.trim(), bundles: o.bundles.trim(), env_cache: o.env_cache.trim() } };
  }

  /* ------------------------------------------------------------ validation */
  function friendly(issue, kindOfRun) {
    var p = issue.path, field = p[p.length - 1];
    if (issue.code === "required") {
      if (p[0] === "runs" && field === "id") return "Give this run an id.";
      return "Fill in the " + (SAMPLE_LABELS[field] || field) + ".";
    }
    if (issue.code === "either_reads_bam") return "Give a reads file or a BAM.";
    if (issue.code === "pattern") return "Use letters, digits, '.', '_' or '-', starting with a letter or digit.";
    if (issue.code === "duplicate") return issue.message.replace(/^\w+ name |^run id /, "") + ".";
    if (issue.code === "flair_source" && kindOfRun) {
      return { release: "Enter a FLAIR version, e.g. 3.0.1.", git: "Enter a branch, tag or commit SHA.",
        env: "Enter the path of a conda env with FLAIR installed.", container: "Enter a container image (docker://… or a .sif path)." }[kindOfRun];
    }
    return issue.message.charAt(0).toUpperCase() + issue.message.slice(1) + ".";
  }
  function fieldKey(issue) {
    var p = issue.path;
    if (p[0] === "samples" && typeof p[1] === "number") {
      var f = p[2] === "bam" || p[2] === "reads" ? "path" : p[2];
      return "s:" + state.samples[p[1]].uid + ":" + f;
    }
    if (p[0] === "opts") return "o:" + p[1];
    if (p[0] === "runs" && typeof p[1] === "number") {
      var r = state.runs[p[1]];
      if (!r) return null;
      if (p[2] === "flair") return "r:" + r.uid + ":" + (r.source === "git" ? "ref" : r.source);
      if (p[2] === "tags") return "r:" + r.uid + ":tags";
      return "r:" + r.uid + ":" + p[2];
    }
    return null;
  }
  function where(issue) {
    var p = issue.path;
    if (p[0] === "samples" && typeof p[1] === "number") {
      var n = sampleName(state.samples[p[1]]);
      return "Sample " + (p[1] + 1) + (n ? " (" + n + ")" : "");
    }
    if (p[0] === "runs" && typeof p[1] === "number") {
      var id = runIds()[p[1]];
      return "Run " + (p[1] + 1) + (id ? " (" + id + ")" : "");
    }
    if (p[0] === "opts") return "Launch options";
    return p[0] === "samples" ? "Samples" : "Runs";
  }
  function collectIssues(samples, runs) {
    var list = FE.validateSamples(samples).filter(function (x) {
      // "fill in the sample name" adds nothing while the name still comes from an empty path
      var s = typeof x.path[1] === "number" ? state.samples[x.path[1]] : null;
      return !(x.code === "required" && x.path[2] === "sample" && s && s.nameAuto && !s.path.trim());
    }).concat(FE.validateRuns({ runs: runs }));
    state.runs.forEach(function (r, i) {
      var bad = parseTags(r.tags).bad;
      if (bad.length) list.push({ path: ["runs", i, "tags"], code: "tag_format", level: "error",
        message: "write tags as key=value, separated by commas (not '" + bad[0] + "')" });
      if (/--ted\b/.test(r.args) && r.source === "release") {
        list.push({ path: ["runs", i, "args"], code: "ted_release", level: "warning",
          message: "TED flags need a FLAIR build with TED (a git branch or conda env); releases do not have them" });
      }
    });
    state.samples.forEach(function (s, i) {
      var isBam = /\.bam$/i.test(s.path.trim());
      if (s.path.trim() && s.kind === "reads" && isBam) {
        list.push({ path: ["samples", i, "reads"], code: "looks_bam", level: "warning",
          message: "this looks like a BAM file; choose Aligned BAM as the input" });
      } else if (s.path.trim() && s.kind === "bam" && !isBam) {
        list.push({ path: ["samples", i, "bam"], code: "not_bam", level: "warning",
          message: "this does not end in .bam; choose Reads if it is a FASTQ/FASTA file" });
      }
    });
    var out = state.opts.outdir.trim();
    if (!out) {
      list.push({ path: ["opts", "outdir"], code: "required", level: "error", message: "missing output folder" });
    } else if (/[$~]/.test(out)) {
      list.push({ path: ["opts", "outdir"], code: "outdir_var", level: "warning",
        message: "the command quotes this path, so $USER or ~ is not expanded; type the full path" });
    } else if (/^\/(private\/)?home\//.test(out)) {
      list.push({ path: ["opts", "outdir"], code: "outdir_home", level: "warning",
        message: "this is in your home folder, which has a small quota; use groups or warm-archive storage" });
    }
    if (state.opts.executor === "local" && !state.opts.profile.trim()) {
      state.samples.forEach(function (s, i) {
        if (s.genomeKind !== "custom") {
          list.push({ path: ["samples", i, "genomeKind"], code: "preset_local", level: "warning",
            message: "the " + s.genomeKind + " preset uses lab files on phoenix; on another computer give a genome FASTA and GTF" });
        }
      });
    }
    var seen = {};
    runs.forEach(function (r, i) {
      var key = JSON.stringify([r.flair, r.transcriptome_args, r.quantify, r.productivity, r.own_align]);
      if (seen[key] !== undefined) {
        list.push({ path: ["runs", i, "args"], code: "same_run", level: "warning",
          message: "same FLAIR and settings as run " + (seen[key] + 1) + "; change the source or the settings" });
      } else seen[key] = i;
    });
    return list.map(function (x) {
      var r = x.path[0] === "runs" && typeof x.path[1] === "number" ? state.runs[x.path[1]] : null;
      x.level = x.level || "error";
      x.key = x.path[2] === "args" && r ? "r:" + r.uid + ":args" : fieldKey(x);
      x.text = friendly(x, r && x.code === "flair_source" ? r.source : null);
      return x;
    });
  }

  /* ------------------------------------------------------------ DOM helpers */
  function h(tag, attrs, kids) {
    var el = document.createElement(tag);
    Object.keys(attrs || {}).forEach(function (k) {
      var v = attrs[k];
      if (v === null || v === undefined || v === false) return;
      if (k === "class") el.className = v;
      else if (k === "text") el.textContent = v;
      else if (k.slice(0, 2) === "on") el.addEventListener(k.slice(2), v);
      else if (k === "value") el.value = v;
      else if (k === "checked") el.checked = !!v;
      else el.setAttribute(k, v === true ? "" : v);
    });
    (kids || []).forEach(function (c) { if (c !== null && c !== undefined) el.appendChild(typeof c === "string" ? document.createTextNode(c) : c); });
    return el;
  }
  function $(id) { return document.getElementById(id); }

  function textField(o) {
    var id = "f-" + o.key.replace(/:/g, "-");
    var input = h("input", { type: "text", id: id, class: o.mono ? "mono" : null, value: o.value, placeholder: o.placeholder || "",
      list: o.list || null, autocomplete: "off", spellcheck: "false", "data-key": o.key,
      oninput: function (e) { o.set(e.target.value); onEdit(); },
      onblur: function () { state.touched[o.key] = true; refreshOutputs(); } });
    return h("div", { class: "field" + (o.wide ? " wide" : "") }, [
      h("label", { for: id, text: o.label }), input,
      o.hint ? h("span", { class: "hint", text: o.hint }) : null,
      h("span", { class: "err", "data-err": o.key, hidden: true })
    ]);
  }
  function selectField(o) {
    var id = "f-" + o.key.replace(/:/g, "-");
    var sel = h("select", { id: id, "data-key": o.key, onchange: function (e) { o.set(e.target.value); onEdit(o.rerender); } },
      o.options.map(function (op) { return h("option", { value: op[0], text: op[1], selected: op[0] === o.value ? "selected" : null }); }));
    return h("div", { class: "field" }, [h("label", { for: id, text: o.label }), sel,
      o.hint ? h("span", { class: "hint", "data-hint": o.key, text: o.hint }) : null,
      h("span", { class: "err", "data-err": o.key, hidden: true })]);
  }
  function seg(label, options, value, set) {
    return h("div", { class: "seg", role: "group", "aria-label": label }, options.map(function (op) {
      return h("button", { type: "button", "aria-pressed": op[0] === value ? "true" : "false", text: op[1],
        onclick: function () { set(op[0]); onEdit(true); } });
    }));
  }
  function check(label, value, set) {
    return h("label", { class: "check" }, [h("input", { type: "checkbox", checked: value,
      onchange: function (e) { set(e.target.checked); onEdit(); } }), label]);
  }

  /* ------------------------------------------------------------ render: samples */
  function renderSamples() {
    var box = $("samples");
    box.textContent = "";
    state.samples.forEach(function (s, i) {
      var k = "s:" + s.uid + ":";
      var nameInput = textField({ key: k + "sample", label: "Sample name", mono: true,
        value: s.nameAuto ? "" : s.sample, placeholder: sampleName(s) || "e.g. HG002",
        hint: s.nameAuto ? "Taken from the file name until you type one." : null,
        set: function (v) { s.sample = v; s.nameAuto = v.trim() === ""; } });
      var card = h("div", { class: "card" }, [
        h("div", { class: "card-head" }, [
          h("span", { class: "title", text: "Sample " + (i + 1) }),
          h("span", { class: "chip", "data-sample-name": s.uid, text: sampleName(s) || "unnamed" }),
          h("span", { class: "spacer" }),
          state.samples.length > 1 ? h("button", { class: "btn ghost", type: "button", "aria-label": "Remove sample " + (i + 1), text: "Remove",
            onclick: function () { state.samples.splice(i, 1); onEdit(true); } }) : null
        ]),
        h("div", { class: "grid two" }, [
          h("div", { class: "field wide" }, [
            h("span", { class: "label", text: "Input" }),
            seg("Input type", [["reads", "Reads (FASTQ/FASTA)"], ["bam", "Aligned BAM"]], s.kind, function (v) { s.kind = v; })
          ]),
          textField({ key: k + "path", label: s.kind === "bam" ? "BAM file (sorted, with .bai next to it)" : "Reads file (.fq, .fa, .gz is fine)",
            mono: true, wide: true, value: s.path, placeholder: s.kind === "bam" ? "/data/HG002.bam" : "/data/HG002.fq.gz",
            set: function (v) { s.path = v; } })
        ]),
        h("div", { class: "field wide" }, [
          h("span", { class: "label", text: "Genome" }),
          seg("Genome", GENOME_KINDS, s.genomeKind, function (v) { s.genomeKind = v; }),
          s.genomeKind !== "custom" ? h("span", { class: "hint", text: FE.GENOME_PRESETS[s.genomeKind].note +
            " (lab reference on phoenix). The GTF below is optional." }) : null,
          h("span", { class: "err", "data-err": k + "genomeKind", hidden: true })
        ]),
        h("div", { class: "grid pair" }, [
          s.genomeKind === "custom" ? textField({ key: k + "genome", label: "Genome FASTA", mono: true, value: s.genome,
            placeholder: "/ref/GRCh38.fa", set: function (v) { s.genome = v; } }) : null,
          textField({ key: k + "gtf", label: s.genomeKind === "custom" ? "Reference annotation GTF" : "Annotation GTF (optional)",
            mono: true, value: s.gtf,
            placeholder: s.genomeKind === "custom" ? "/ref/gencode.gtf" : "from the preset: GENCODE v49 basic",
            set: function (v) { s.gtf = v; } })
        ]),
        h("div", { class: "grid" }, [
          nameInput,
          selectField({ key: k + "platform", label: "Platform", value: s.platform,
            options: [["", "Not set"], ["pacbio", "PacBio"], ["ont", "Oxford Nanopore"]], set: function (v) { s.platform = v; } }),
          selectField({ key: k + "library", label: "Library", value: s.library,
            options: [["", "Not set"], ["cDNA", "cDNA"], ["dRNA", "Direct RNA"]], set: function (v) { s.library = v; } })
        ]),
        h("details", { class: "more", open: (s.cage || s.pas || s.junctions || s.spikein_gtf) ? "open" : null }, [
          h("summary", { text: "Extra tracks (optional, for later report sections)" }),
          h("div", { class: "grid two" }, [
            textField({ key: k + "cage", label: "CAGE peaks (TSS)", mono: true, value: s.cage, placeholder: "optional", set: function (v) { s.cage = v; } }),
            textField({ key: k + "pas", label: "PAS / 3′-seq peaks (TES)", mono: true, value: s.pas, placeholder: "optional", set: function (v) { s.pas = v; } }),
            textField({ key: k + "junctions", label: "Short-read junctions (STAR SJ.out.tab)", mono: true, value: s.junctions, placeholder: "optional", set: function (v) { s.junctions = v; } }),
            textField({ key: k + "spikein_gtf", label: "Spike-in GTF (e.g. SIRVs)", mono: true, value: s.spikein_gtf, placeholder: "optional", set: function (v) { s.spikein_gtf = v; } })
          ])
        ])
      ]);
      box.appendChild(card);
    });
    $("sample-count").textContent = state.samples.length + (state.samples.length === 1 ? " sample" : " samples");
  }

  /* ------------------------------------------------------------ render: runs */
  function presetOf(args) {
    var a = args.trim();
    var p = FE.PRESETS.filter(function (x) { return x.args === a; })[0];
    return p ? p.key : "custom";
  }
  function renderRuns() {
    var box = $("runs");
    box.textContent = "";
    var ids = runIds();
    state.runs.forEach(function (r, i) {
      var k = "r:" + r.uid + ":";
      var src;
      if (r.source === "release") {
        src = [textField({ key: k + "release", label: "FLAIR version", mono: true, value: r.release, list: "releases",
          placeholder: "3.0.1", hint: "Installed from bioconda and cached.", set: function (v) { r.release = v; } })];
      } else if (r.source === "git") {
        src = [textField({ key: k + "ref", label: "Branch, tag or commit", mono: true, value: r.ref, list: "refs", placeholder: "dev",
                 hint: "Pinned to the exact commit when the run starts.", set: function (v) { r.ref = v; } }),
               textField({ key: k + "git_url", label: "Repository", mono: true, value: r.git_url, placeholder: FE.DEFAULT_GIT_URL,
                 set: function (v) { r.git_url = v; } })];
      } else if (r.source === "env") {
        src = [textField({ key: k + "env", label: "Conda env path", mono: true, wide: true, value: r.env,
          placeholder: "/path/to/conda/envs/flair-dev", hint: "An existing env with FLAIR installed. Used as is.",
          set: function (v) { r.env = v; } })];
      } else {
        src = [textField({ key: k + "container", label: "Container image", mono: true, wide: true, value: r.container,
          placeholder: "docker://quay.io/biocontainers/flair:3.0.1--pyhdfd78af_0 or /images/flair.sif",
          set: function (v) { r.container = v; } })];
      }
      var preset = FE.PRESETS.filter(function (p) { return p.key === r.preset; })[0];
      var card = h("div", { class: "card" }, [
        h("div", { class: "card-head" }, [
          h("span", { class: "title", text: "Run " + (i + 1) }),
          h("span", { class: "chip accent", "data-run-id": r.uid, text: ids[i] || "no id" }),
          h("span", { class: "chip", "data-run-token": r.uid }),
          h("span", { class: "spacer" }),
          h("button", { class: "btn ghost", type: "button", text: "Duplicate", "aria-label": "Duplicate run " + (i + 1),
            onclick: function () {
              var c = JSON.parse(JSON.stringify(r)); c.uid = uid(); c.idAuto = true; c.id = ""; c.label = "";
              state.runs.splice(i + 1, 0, c); onEdit(true);
            } }),
          state.runs.length > 1 ? h("button", { class: "btn ghost", type: "button", text: "Remove", "aria-label": "Remove run " + (i + 1),
            onclick: function () { state.runs.splice(i, 1); onEdit(true); } }) : null
        ]),
        h("div", { class: "field" }, [
          h("span", { class: "label", text: "Which FLAIR" }),
          seg("FLAIR source", SOURCES, r.source, function (v) {
            r.source = v;
            if (v === "release" && !r.release) r.release = "3.0.1";
            if (v === "git" && !r.ref) r.ref = "dev";
          })
        ]),
        h("div", { class: "grid two" }, src),
        h("div", { class: "grid two" }, [
          selectField({ key: k + "preset", label: "Settings", value: presetOf(r.args) === "custom" ? "custom" : r.preset,
            options: FE.PRESETS.map(function (p) { return [p.key, p.label]; }).concat([["custom", "Custom arguments"]]),
            hint: preset && presetOf(r.args) !== "custom" ? preset.note : "Your own flair transcriptome arguments.",
            rerender: true,
            set: function (v) {
              r.preset = v;
              var p = FE.PRESETS.filter(function (x) { return x.key === v; })[0];
              if (p) r.args = p.args;
            } }),
          textField({ key: k + "args", label: "flair transcriptome arguments", mono: true, value: r.args, placeholder: "none",
            set: function (v) { r.args = v; r.preset = presetOf(v); } })
        ]),
        h("details", { class: "more", open: (r.label || r.tags || !r.idAuto || !r.quantify || !r.productivity || r.own_align) ? "open" : null }, [
          h("summary", { text: "Name, tags and steps" }),
          h("div", { class: "col", style: "gap:12px" }, [
            h("div", { class: "grid" }, [
              textField({ key: k + "id", label: "Run id", mono: true, value: r.idAuto ? "" : r.id, placeholder: ids[i],
                hint: r.idAuto ? "Made from the FLAIR source until you type one." : null,
                set: function (v) { r.id = v; r.idAuto = v.trim() === ""; } }),
              textField({ key: k + "label", label: "Label in the report", value: r.label, placeholder: ids[i],
                set: function (v) { r.label = v; } }),
              textField({ key: k + "tags", label: "Tags", mono: true, value: r.tags, placeholder: "version=3.0.1, setting=stringent",
                hint: "Filters in the report.", set: function (v) { r.tags = v; } })
            ]),
            h("div", { class: "checks" }, [
              check("Quantify (flair quantify)", r.quantify, function (v) { r.quantify = v; }),
              check("Predict productivity", r.productivity, function (v) { r.productivity = v; }),
              check("Align its own BAM", r.own_align, function (v) { r.own_align = v; })
            ])
          ])
        ])
      ]);
      box.appendChild(card);
    });
    $("run-count").textContent = state.runs.length + (state.runs.length === 1 ? " run" : " runs");
  }

  function renderOptions() {
    Array.prototype.forEach.call(document.querySelectorAll("[data-opt]"), function (el) {
      var k = el.getAttribute("data-opt");
      if (el.type === "checkbox") el.checked = !!state.opts[k];
      else if (el.id !== "o-baseline") el.value = state.opts[k] || "";
    });
    var local = state.opts.executor === "local";
    $("software-field").hidden = !local;
    $("o-envcache").placeholder = local ? "./flair-env-cache" : "/private/groups/brookslab/flair-eval/env_cache (shared)";
  }
  function renderBaseline(ids) {
    var sel = $("o-baseline");
    var cur = state.opts.baseline;
    sel.textContent = "";
    ids.forEach(function (id, i) {
      sel.appendChild(h("option", { value: i === 0 ? "" : id, text: id + (i === 0 ? " (first run)" : ""),
        selected: (i === 0 && !cur) || id === cur ? "selected" : null }));
    });
  }

  /* ------------------------------------------------------------ outputs */
  var lastFiles = { samples: "", runs: "" };
  function refreshOutputs() {
    var samples = state.samples.map(toSample);
    var runs = toRuns();
    var ids = runs.map(function (r) { return r.id || ""; });
    var issues = collectIssues(samples, runs);
    var errors = issues.filter(function (x) { return x.level === "error"; });
    var warns = issues.filter(function (x) { return x.level === "warning"; });

    // live chips
    state.samples.forEach(function (s) {
      var chip = document.querySelector('[data-sample-name="' + s.uid + '"]');
      if (chip) chip.textContent = sampleName(s) || "unnamed";
      var inp = $("f-s-" + s.uid + "-sample");
      if (inp && s.nameAuto) inp.placeholder = sampleName(s) || "e.g. HG002";
    });
    state.runs.forEach(function (r, i) {
      var idChip = document.querySelector('[data-run-id="' + r.uid + '"]');
      if (idChip) idChip.textContent = ids[i] || "no id";
      var tokChip = document.querySelector('[data-run-token="' + r.uid + '"]');
      if (tokChip) {
        var tok = FE.tokenForRun(runs[i]);
        tokChip.textContent = tok ? "--flair " + tok : "runs.yaml only";
        tokChip.title = tok ? "This run can be written as a --flair token" : "Its id or source needs runs.yaml";
      }
      ["id", "label"].forEach(function (f) {
        var inp = $("f-r-" + r.uid + "-" + f);
        if (inp && (f === "label" || r.idAuto)) inp.placeholder = ids[i];
      });
    });
    renderBaseline(ids);

    // field errors
    Array.prototype.forEach.call(document.querySelectorAll("[data-err]"), function (el) { el.hidden = true; el.textContent = ""; });
    Array.prototype.forEach.call(document.querySelectorAll("[aria-invalid]"), function (el) { el.removeAttribute("aria-invalid"); });
    issues.forEach(function (x) {
      if (!x.key) return;
      var show = state.touched[x.key] || state.example || x.code === "duplicate" || x.level === "warning" || x.code === "pattern";
      if (!show) return;
      var e = document.querySelector('[data-err="' + x.key + '"]');
      var inp = document.querySelector('[data-key="' + x.key + '"]');
      if (e) { e.hidden = false; e.textContent = x.text; e.className = x.level === "warning" ? "warn-text" : "err"; }
      if (inp && x.level === "error") inp.setAttribute("aria-invalid", "true");
    });

    // status
    var pill = $("pill");
    if (errors.length) { pill.className = "pill bad"; pill.textContent = errors.length + (errors.length === 1 ? " thing to fix" : " things to fix"); }
    else if (warns.length) { pill.className = "pill warn"; pill.textContent = "Ready, " + warns.length + (warns.length === 1 ? " note" : " notes"); }
    else { pill.className = "pill ok"; pill.textContent = "Ready to run"; }
    var nj = samples.length * runs.length;
    $("jobs").textContent = samples.length + " sample" + (samples.length === 1 ? "" : "s") + " × " + runs.length + " run" +
      (runs.length === 1 ? "" : "s") + " = " + nj + " job" + (nj === 1 ? "" : "s");
    var ul = $("issues");
    ul.textContent = "";
    issues.forEach(function (x) {
      ul.appendChild(h("li", { class: x.level }, [h("button", { type: "button", onclick: function () { focusIssue(x); } },
        [h("span", { class: "where", text: where(x) + ": " }), x.text])]));
    });

    // command (missing required fields show as <placeholders> until filled in)
    var shown = samples.map(function (x, i) {
      var y = Object.assign({}, x);
      if (!y.reads && !y.bam) y[state.samples[i].kind] = state.samples[i].kind === "bam" ? "<aligned.bam>" : "<reads.fq.gz>";
      if (!y.sample) y.sample = "<sample-name>";
      if (!y.genome) y.genome = "<genome.fa>";
      if (!y.gtf && !FE.isGenomePreset(y.genome)) y.gtf = "<annotation.gtf>";
      return y;
    });
    var c = FE.buildCommand(shown, runs, cmdOptions(runs));
    $("cmd").textContent = FE.formatCommand(c.argv);
    var why = $("why");
    why.textContent = "";
    $("launch-hint").hidden = state.opts.executor === "local";
    if (c.mode === "quick") {
      $("mode").textContent = "One command, no files needed";
      $("mode-sub").textContent = "Run it where the paths exist.";
      why.hidden = true;
    } else {
      $("mode").textContent = "Uses samples.csv and runs.yaml";
      $("mode-sub").textContent = "Download both files next to where you launch.";
      var reasons = c.reasons.filter(function (x) { return x !== "files requested"; });
      why.hidden = reasons.length === 0;
      reasons.slice(0, 5).forEach(function (x) { why.appendChild(h("li", { text: "Needs files: " + x })); });
    }
    lastCommand = FE.oneLine(c.argv);

    // files
    lastFiles.samples = FE.samplesCsv(samples);
    lastFiles.runs = FE.runsYaml(runs);
    lastMode = c.mode;
    renderFile();
    save();
  }
  var lastCommand = "", lastMode = "quick";

  function renderFile() {
    var mode = lastMode;
    var tab = state.tab === "runs" ? "runs" : "samples";
    $("tab-samples").setAttribute("aria-selected", tab === "samples" ? "true" : "false");
    $("tab-runs").setAttribute("aria-selected", tab === "runs" ? "true" : "false");
    $("file-view").textContent = lastFiles[tab];
    $("dl-file").textContent = "Download " + (tab === "runs" ? "runs.yaml" : "samples.csv");
    $("file-hint").textContent = mode === "quick" ? "Optional for this command." : "Needed for this command.";
  }

  function focusIssue(x) {
    if (!x.key) return;
    state.touched[x.key] = true;
    var inp = document.querySelector('[data-key="' + x.key + '"]');
    if (inp) {
      var d = inp.closest("details");
      if (d) d.open = true;
      inp.focus();
      inp.scrollIntoView({ block: "center", behavior: "smooth" });
    }
    refreshOutputs();
  }

  function onEdit(rerender) {
    state.example = false;
    $("example-banner").hidden = true;
    if (rerender) renderAll();
    else refreshOutputs();
  }
  function renderAll() {
    renderSamples();
    renderRuns();
    renderOptions();
    $("example-banner").hidden = !state.example;
    refreshOutputs();
  }

  /* ------------------------------------------------------------ copy / download */
  function toast(msg) {
    var t = $("toast");
    t.textContent = msg;
    clearTimeout(toast.timer);
    toast.timer = setTimeout(function () { t.textContent = ""; }, 2500);
  }
  function selectText(el) {
    var range = document.createRange();
    range.selectNodeContents(el);
    var sel = window.getSelection();
    sel.removeAllRanges();
    sel.addRange(range);
  }
  function copy(text, el, label) {
    var done = function () { toast(label + " copied"); };
    var fail = function () { selectText(el); toast("Press Ctrl+C / ⌘C to copy"); };
    try {
      if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(text).then(done, fail);
      else fail();
    } catch (e) { fail(); }
  }
  function copySnip(text, el, btn) {
    var flash = function (msg) { btn.textContent = msg; setTimeout(function () { btn.textContent = "Copy"; }, 1500); };
    try {
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(text).then(function () { flash("Copied"); }, function () { selectText(el); flash("Ctrl+C"); });
      } else { selectText(el); flash("Ctrl+C"); }
    } catch (e) { selectText(el); flash("Ctrl+C"); }
  }
  function download(name, text) {
    try {
      var blob = new Blob([text], { type: name.slice(-4) === ".csv" ? "text/csv" : "text/yaml" });
      var url = URL.createObjectURL(blob);
      var a = h("a", { href: url, download: name });
      document.body.appendChild(a);
      a.click();
      setTimeout(function () { URL.revokeObjectURL(url); a.remove(); }, 1000);
      toast(name + " downloaded");
    } catch (e) { toast("Download blocked here; use Copy"); }
  }

  /* ------------------------------------------------------------ wiring */
  document.body.appendChild(h("datalist", { id: "releases" }, RELEASES.map(function (v) { return h("option", { value: v }); })));
  document.body.appendChild(h("datalist", { id: "refs" }, ["dev", "main"].map(function (v) { return h("option", { value: v }); })));

  $("add-sample").addEventListener("click", function () {
    state.samples.push(blankSample(state.samples[state.samples.length - 1]));
    onEdit(true);
  });
  Array.prototype.forEach.call(document.querySelectorAll("[data-add-run]"), function (b) {
    b.addEventListener("click", function () { state.runs.push(blankRun(b.getAttribute("data-add-run"))); onEdit(true); });
  });
  Array.prototype.forEach.call(document.querySelectorAll("[data-opt]"), function (el) {
    var k = el.getAttribute("data-opt");
    el.addEventListener(el.type === "checkbox" || el.tagName === "SELECT" ? "change" : "input", function () {
      state.opts[k] = el.type === "checkbox" ? el.checked : el.value;
      if (k === "executor") renderOptions();
      onEdit();
    });
    if (el.getAttribute("data-key")) {
      el.addEventListener("blur", function () { state.touched[el.getAttribute("data-key")] = true; refreshOutputs(); });
    }
  });
  function reset(toExample) {
    state = toExample ? exampleState() : blankState();
    renderAll();
  }
  $("btn-example").addEventListener("click", function () { reset(true); });
  $("btn-reset").addEventListener("click", function () { reset(false); });
  $("btn-reset-2").addEventListener("click", function () { reset(false); });
  $("copy-cmd").addEventListener("click", function () { copy(lastCommand, $("cmd"), "Command"); });
  $("copy-file").addEventListener("click", function () {
    var tab = state.tab === "runs" ? "runs" : "samples";
    copy(lastFiles[tab], $("file-view"), tab === "runs" ? "runs.yaml" : "samples.csv");
  });
  $("dl-file").addEventListener("click", function () {
    var tab = state.tab === "runs" ? "runs" : "samples";
    download(tab === "runs" ? "runs.yaml" : "samples.csv", lastFiles[tab]);
  });
  $("dl-both").addEventListener("click", function () {
    download("samples.csv", lastFiles.samples);
    setTimeout(function () { download("runs.yaml", lastFiles.runs); }, 400);
  });
  Array.prototype.forEach.call(document.querySelectorAll("[data-copy]"), function (b) {
    b.addEventListener("click", function () {
      var el = $(b.getAttribute("data-copy"));
      copySnip(el.textContent, el, b);
    });
  });
  Array.prototype.forEach.call(document.querySelectorAll("[data-tab]"), function (b) {
    b.addEventListener("click", function () { state.tab = b.getAttribute("data-tab"); renderFile(); save(); });
  });

  // On GitHub Pages the builder sits at /builder/ next to the landing page and demo report.
  if (/^https?:$/.test(location.protocol) && /\/builder\/(index\.html)?$/.test(location.pathname)) {
    $("home-link").hidden = false;
    $("demo-link").hidden = false;
  }

  renderAll();
})();
