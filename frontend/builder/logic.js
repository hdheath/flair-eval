/* flair-eval config builder: pure logic, no DOM.
 *
 * This file mirrors the Python package in cli/src/flair_eval (tokens.py, config.py,
 * validate.py, command.py, presets.py). cli/tests/test_builder_parity.py runs both on the
 * same inputs with Node and fails if they disagree, so change them together.
 * Plain ES2017 so it runs in any current browser and in old Node without a build step.
 */
(function (root, factory) {
  var api = factory();
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.FlairEval = api;
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  var PIPELINE = "hdheath/flair-eval";
  var DEFAULT_REF = "v2";
  var DEFAULT_GIT_URL = "https://github.com/BrooksLabUCSC/flair";
  var GIT_SHORTCUTS = ["dev", "main"];
  var ID_RE = /^[A-Za-z0-9][A-Za-z0-9._-]*$/;
  var NAME_HINT = "start with a letter or digit, then use only letters, digits, '.', '_' or '-'";
  var FLAIR_KEYS = ["release", "git", "env", "container"];

  var SAMPLE_COLUMNS = ["sample", "reads", "bam", "genome", "gtf", "platform", "library"];
  var TRACK_COLUMNS = ["cage", "pas", "junctions", "spikein_gtf"];
  var ALL_SAMPLE_COLUMNS = SAMPLE_COLUMNS.concat(TRACK_COLUMNS);
  var RUN_DEFAULTS = { transcriptome_args: "", quantify: true, productivity: true, own_align: false };

  var PRESETS = [
    { key: "default", label: "FLAIR defaults", args: "", note: "No extra arguments." },
    { key: "stringent", label: "Stringent", args: "--stringent",
      note: "Supporting reads must be full-length (span 25 bp of the first and last exons)." },
    { key: "check_splice", label: "Check splice sites", args: "--check_splice",
      note: "Reads must cover 4 of 6 bp around each splice site, with no insertion over 3 bp." },
    { key: "stringent_splice", label: "Stringent + check splice sites",
      args: "--stringent --check_splice", note: "Both of the above." },
    { key: "support5", label: "More read support (5)", args: "--support 5",
      note: "At least 5 supporting reads per isoform (default 3)." },
    { key: "noaligntoannot", label: "Skip annotation alignment", args: "--noaligntoannot",
      note: "Detect transcripts from the genome alignment only. Faster, less accurate." },
    { key: "ted_pacbio", label: "TED ends, PacBio cDNA", args: "--ted --ted_chemistry pacbio_cdna",
      note: "Transcript End Determination. Only for FLAIR builds with TED." },
    { key: "ted_ont", label: "TED ends, ONT cDNA", args: "--ted --ted_chemistry ont_cdna",
      note: "Transcript End Determination. Only for FLAIR builds with TED." }
  ];

  var PARAM_DEFAULTS = {
    outdir: "./flair-eval-results", report_title: "FLAIR comparison", baseline: null,
    align_flair: "3.0.1", regions: null, bundles: null, env_cache: null
  };
  var PARAM_ORDER = ["bundles", "baseline", "report_title", "align_flair", "regions", "env_cache"];

  function has(o, k) { return o !== null && typeof o === "object" && Object.prototype.hasOwnProperty.call(o, k); }
  function isObj(v) { return v !== null && typeof v === "object" && !Array.isArray(v); }
  function repr(s) { return "'" + String(s) + "'"; }

  /* ------------------------------------------------------------ tokens */
  function TokenError(message) { this.name = "TokenError"; this.message = message; }
  TokenError.prototype = Object.create(Error.prototype);

  function sanitize(s) { return s.replace(/[^A-Za-z0-9._-]/g, "-"); }
  function basename(p) { return p.replace(/\/+$/, "").split("/").pop(); }
  function looksLikeUrl(s) {
    return s.indexOf("://") >= 0 || s.indexOf("git@") === 0 || s.indexOf("/") === 0 || s.indexOf("file:") === 0;
  }
  function containerName(uri) {
    var i = uri.indexOf("://");
    var base = basename(i >= 0 ? uri.slice(i + 3) : uri);
    base = base.split("@")[0];
    var c = base.lastIndexOf(":");
    if (c >= 0) base = base.slice(0, c);
    return base;
  }
  function makeId(prefix, rest, token) {
    var rid = prefix + sanitize(rest);
    if (!ID_RE.test(rid) || !rest) {
      throw new TokenError("--flair token " + repr(token) + " gives an unusable run id " + repr(rid) +
        "; write this run in runs.yaml instead");
    }
    return rid;
  }
  function tr(token, id, flair) { return { token: token, id: id, flair: flair }; }

  function parseToken(token) {
    var t = String(token).trim();
    if (!t) throw new TokenError("empty --flair token (check for a stray comma)");
    var v;
    if (t.indexOf("release:") === 0) {
      v = t.slice(8).trim();
      if (!v) throw new TokenError("--flair token " + repr(t) + ": give a version, e.g. release:3.0.1");
      return tr(t, makeId("flair-", v, t), { release: v });
    }
    if (t.indexOf("git:") === 0) {
      var rest = t.slice(4).trim();
      if (!rest) throw new TokenError("--flair token " + repr(t) + ": give a branch, tag or SHA, e.g. git:dev");
      var url = DEFAULT_GIT_URL, ref = rest;
      if (looksLikeUrl(rest)) {
        var at = rest.lastIndexOf("@");
        if (at > 0) { url = rest.slice(0, at); ref = rest.slice(at + 1); } else { url = ""; ref = ""; }
        if (!looksLikeUrl(url) || !ref || (ref.indexOf("/") >= 0 && ref.indexOf(":") >= 0)) {
          throw new TokenError("--flair token " + repr(t) + ": a git URL needs @<ref>, " +
            "e.g. git:https://github.com/me/flair@my-branch");
        }
      }
      return tr(t, makeId("flair-", ref, t), { git: url, ref: ref });
    }
    if (t.indexOf("env:") === 0) {
      v = t.slice(4).trim();
      if (!v) throw new TokenError("--flair token " + repr(t) + ": give a conda env path, e.g. env:/path/to/env");
      return tr(t, makeId("env-", basename(v), t), { env: v });
    }
    if (t.indexOf("container:") === 0) {
      v = t.slice(10).trim();
      if (!v) throw new TokenError("--flair token " + repr(t) + ": give an image, e.g. " +
        "container:docker://org/flair:3.0.1 or container:/path/flair.sif");
      return tr(t, makeId("ctr-", containerName(v), t), { container: v });
    }
    if (GIT_SHORTCUTS.indexOf(t) >= 0) return tr(t, "flair-" + t, { git: DEFAULT_GIT_URL, ref: t });
    if (/^[0-9]/.test(t)) return tr(t, makeId("flair-", t, t), { release: t });
    throw new TokenError("unknown --flair token " + repr(t) + ". Use a release (3.0.1), dev, main, " +
      "git:" + t + " for another branch, env:/path or container:<uri>");
  }

  function parseTokens(value) {
    if (value === null || value === undefined || !String(value).trim()) {
      throw new TokenError("--flair is empty; give at least one token, e.g. --flair 3.0.1");
    }
    var runs = String(value).split(",").map(parseToken);
    var seen = {};
    runs.forEach(function (r) {
      if (has(seen, r.id)) {
        throw new TokenError("--flair tokens " + repr(seen[r.id]) + " and " + repr(r.token) +
          " both give run id " + repr(r.id) + "; write these runs in runs.yaml with distinct ids");
      }
      seen[r.id] = r.token;
    });
    return runs;
  }

  function sameFlat(a, b) {
    var ka = Object.keys(a).sort(), kb = Object.keys(b).sort();
    if (ka.join("\u0000") !== kb.join("\u0000")) return false;
    return ka.every(function (k) { return a[k] === b[k]; });
  }
  function keysAre(o, keys) { return Object.keys(o).sort().join(",") === keys.slice().sort().join(","); }

  function tokenForFlair(flair) {
    if (!isObj(flair) || Object.keys(flair).length === 0) return null;
    var tok;
    if (keysAre(flair, ["release"])) {
      var v = String(flair.release);
      tok = /^[0-9]/.test(v) ? v : "release:" + v;
    } else if (keysAre(flair, ["git", "ref"])) {
      var url = String(flair.git), ref = String(flair.ref);
      if (url === DEFAULT_GIT_URL) tok = GIT_SHORTCUTS.indexOf(ref) >= 0 ? ref : "git:" + ref;
      else tok = "git:" + url + "@" + ref;
    } else if (keysAre(flair, ["env"])) {
      tok = "env:" + String(flair.env);
    } else if (keysAre(flair, ["container"])) {
      tok = "container:" + String(flair.container);
    } else return null;
    if (tok.indexOf(",") >= 0 || tok !== tok.trim()) return null;
    var back;
    try { back = parseToken(tok); } catch (e) { return null; }
    return sameFlat(back.flair, flair) ? tok : null;
  }

  function tokenForRun(run) {
    var tok = tokenForFlair(run.flair);
    if (tok === null) return null;
    return parseToken(tok).id === run.id ? tok : null;
  }

  /* ------------------------------------------------------------ config files */
  var READ_EXT = /(\.(fastq|fq|fasta|fa|fna|bam|sam|cram))?(\.(gz|bgz|bz2|xz|zst))?$/i;
  function defaultSampleName(path) {
    var base = String(path).replace(/\/+$/, "").split("/").pop();
    var name = base.replace(READ_EXT, "");
    return name || base;
  }

  function cleanSample(row) {
    var out = {};
    Object.keys(row).forEach(function (k) {
      var v = row[k];
      if (v === null || v === undefined) return;
      v = String(v).trim();
      if (v !== "") out[String(k).trim()] = v;
    });
    return out;
  }

  function csvField(v) {
    var s = v === null || v === undefined ? "" : String(v);
    return /[",\r\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
  }
  function samplesCsv(samples) {
    var cols = SAMPLE_COLUMNS.concat(TRACK_COLUMNS.filter(function (c) {
      return samples.some(function (s) { return !!s[c]; });
    }));
    var lines = [cols.join(",")];
    samples.forEach(function (s) {
      lines.push(cols.map(function (c) { return csvField(s[c] || ""); }).join(","));
    });
    return lines.join("\n") + "\n";
  }

  function q(s) { return JSON.stringify(String(s)); }
  function runsYaml(runs) {
    var L = ["# flair-eval runs (schema flair-eval/runs/1). One entry per FLAIR configuration.", "runs:"];
    runs.forEach(function (r) {
      L.push("  - id: " + q(r.id));
      if (r.label) L.push("    label: " + q(r.label));
      L.push("    flair:");
      ["release", "git", "ref", "env", "container"].forEach(function (k) {
        if (has(r.flair || {}, k)) L.push("      " + k + ": " + q(r.flair[k]));
      });
      L.push("    transcriptome_args: " + q(has(r, "transcriptome_args") ? r.transcriptome_args : ""));
      ["quantify", "productivity", "own_align"].forEach(function (k) {
        var v = has(r, k) ? r[k] : RUN_DEFAULTS[k];
        L.push("    " + k + ": " + (v ? "true" : "false"));
      });
      var tags = r.tags || {};
      if (Object.keys(tags).length) {
        L.push("    tags:");
        Object.keys(tags).forEach(function (k) { L.push("      " + q(k) + ": " + q(tags[k])); });
      }
    });
    return L.join("\n") + "\n";
  }
  function runIsDefault(run, key) { return (has(run, key) ? run[key] : RUN_DEFAULTS[key]) === RUN_DEFAULTS[key]; }

  /* ------------------------------------------------------------ JSON Schema subset */
  function typeOk(v, t) {
    if (Array.isArray(t)) return t.some(function (x) { return typeOk(v, x); });
    switch (t) {
      case "string": return typeof v === "string";
      case "boolean": return typeof v === "boolean";
      case "object": return isObj(v);
      case "array": return Array.isArray(v);
      case "integer": return typeof v === "number" && Number.isInteger(v);
      case "number": return typeof v === "number";
      case "null": return v === null;
    }
    return true;
  }
  function deepEqual(a, b) { return JSON.stringify(a) === JSON.stringify(b); }

  /* Raw errors {kw, path, inst, sch, value}, one per failing keyword (like jsonschema's
   * iter_errors; errors inside anyOf/oneOf are folded into the anyOf/oneOf error). */
  function walk(inst, sch, path, out) {
    if (has(sch, "type") && !typeOk(inst, sch.type)) out.push({ kw: "type", path: path, inst: inst, value: sch.type });
    if (has(sch, "const") && !deepEqual(inst, sch.const)) out.push({ kw: "const", path: path, inst: inst, value: sch.const });
    if (has(sch, "enum") && !sch.enum.some(function (e) { return deepEqual(e, inst); })) {
      out.push({ kw: "enum", path: path, inst: inst, value: sch.enum });
    }
    if (typeof inst === "string") {
      if (has(sch, "minLength") && Array.from(inst).length < sch.minLength) out.push({ kw: "minLength", path: path, inst: inst });
      if (has(sch, "pattern") && !new RegExp(sch.pattern, "u").test(inst)) out.push({ kw: "pattern", path: path, inst: inst });
    }
    if (typeof inst === "number" && has(sch, "minimum") && inst < sch.minimum) out.push({ kw: "minimum", path: path, inst: inst });
    if (Array.isArray(inst)) {
      if (has(sch, "minItems") && inst.length < sch.minItems) out.push({ kw: "minItems", path: path, inst: inst });
      if (has(sch, "items")) inst.forEach(function (x, i) { walk(x, sch.items, path.concat([i]), out); });
    }
    if (isObj(inst)) {
      (sch.required || []).forEach(function (k) {
        if (!has(inst, k)) out.push({ kw: "required", path: path, inst: inst, field: k });
      });
      var props = sch.properties || {};
      Object.keys(props).forEach(function (k) { if (has(inst, k)) walk(inst[k], props[k], path.concat([k]), out); });
      if (has(sch, "additionalProperties")) {
        var extras = Object.keys(inst).filter(function (k) { return !has(props, k); });
        if (sch.additionalProperties === false && extras.length) {
          out.push({ kw: "additionalProperties", path: path, inst: inst, sch: sch });
        } else if (isObj(sch.additionalProperties)) {
          extras.forEach(function (k) { walk(inst[k], sch.additionalProperties, path.concat([k]), out); });
        }
      }
    }
    function passes(s) { var o = []; walk(inst, s, path, o); return o.length === 0; }
    if (has(sch, "anyOf") && !sch.anyOf.some(passes)) out.push({ kw: "anyOf", path: path, inst: inst });
    if (has(sch, "oneOf") && sch.oneOf.filter(passes).length !== 1) out.push({ kw: "oneOf", path: path, inst: inst });
  }

  function typeName(v) {
    if (typeof v === "boolean") return "true/false";
    if (typeof v === "number") return "a number";
    if (v === null || v === undefined) return "empty";
    if (Array.isArray(v)) return "a list";
    if (isObj(v)) return "a mapping";
    return "text";
  }
  var NICE_TYPE = { string: "text", boolean: "true or false", object: "a mapping", array: "a list",
    integer: "a whole number", number: "a number" };

  function issue(path, code, message, level) { return { path: path, code: code, message: message, level: level || "error" }; }

  function flairIssue(path, inst) {
    var msg = "give exactly one FLAIR source: release, git (with ref), env or container";
    if (isObj(inst)) {
      var given = FLAIR_KEYS.filter(function (k) { return has(inst, k); });
      var extra = Object.keys(inst).filter(function (k) { return FLAIR_KEYS.indexOf(k) < 0 && k !== "ref"; }).sort();
      var nonEmpty = function (v) { return typeof v === "string" && v.trim() !== ""; };
      if (given.length > 1) msg = "give only one FLAIR source, not " + given.join(" and ");
      else if (given.length === 1 && given[0] === "git" && !has(inst, "ref")) msg = "git needs a ref: a branch, tag or commit SHA (e.g. ref: dev)";
      else if (has(inst, "ref") && !(given.length === 1 && given[0] === "git")) msg = "ref only goes with git";
      else if (extra.length) msg = "unknown key '" + extra[0] + "' in flair; " + msg;
      else if (given.length && !nonEmpty(inst[given[0]])) {
        var v = inst[given[0]];
        msg = given[0] + " should be non-empty text";
        if (typeof v === "number") msg += ' (put it in quotes, e.g. ' + given[0] + ': "' + v + '")';
      } else if (given.length === 1 && given[0] === "git" && !nonEmpty(inst.ref)) {
        msg = "ref should be non-empty text";
        if (typeof inst.ref === "number") msg += ' (put it in quotes, e.g. ref: "' + inst.ref + '")';
      }
    }
    return issue(path, "flair_source", msg);
  }

  function issuesFromSchema(instance, sch, prefix, kind) {
    var raw = [];
    walk(instance, sch, prefix, raw);
    var out = [];
    raw.forEach(function (e) {
      var path = e.path, inst = e.inst;
      var last = path[path.length - 1];
      if (e.kw === "required") {
        out.push(issue(path.concat([e.field]), "required", "missing required " + kind + " '" + e.field + "'"));
      } else if (e.kw === "anyOf" && isObj(inst) && kind === "column") {
        out.push(issue(path.concat(["reads"]), "either_reads_bam", "give reads (FASTQ/FASTA) or bam (sorted, indexed BAM)"));
      } else if (e.kw === "oneOf" && last === "flair") {
        out.push(flairIssue(path, inst));
      } else if (e.kw === "additionalProperties") {
        var allowed = e.sch.properties || {};
        Object.keys(inst).filter(function (k) { return !has(allowed, k); }).sort().forEach(function (k) {
          var msg = "unknown " + kind + " '" + k + "'";
          if (kind === "key" && path.length === 2 && path[0] === "runs") {
            msg += "; allowed keys are id, flair, transcriptome_args, quantify, productivity, own_align, label, tags";
          } else if (kind === "column") msg += "; allowed columns are " + ALL_SAMPLE_COLUMNS.join(", ");
          out.push(issue(path.concat([k]), "unknown_key", msg));
        });
      } else if (e.kw === "pattern") {
        out.push(issue(path, "pattern", repr(inst) + " is not a valid name: " + NAME_HINT));
      } else if (e.kw === "enum") {
        out.push(issue(path, "enum", repr(inst) + " is not allowed; use one of " + e.value.join(", ")));
      } else if (e.kw === "type") {
        var want = Array.isArray(e.value) ? e.value.join("/") : e.value;
        var msg = "should be " + (NICE_TYPE[want] || want) + ", not " + typeName(inst);
        if (want === "string" && typeof inst === "number") msg += ' (put it in quotes, e.g. "' + inst + '")';
        out.push(issue(path, "type", msg));
      } else if (e.kw === "minItems") {
        out.push(issue(path, "min_items", "needs at least one run"));
      } else if (e.kw === "minLength") {
        out.push(issue(path, "min_length", "must not be empty"));
      } else if (e.kw === "const") {
        out.push(issue(path, "const", "should be " + repr(e.value)));
      } else {
        out.push(issue(path, e.kw, e.kw + " failed"));
      }
    });
    return out;
  }

  function duplicates(items, key, prefix, what) {
    var out = [], first = {};
    items.forEach(function (it, i) {
      var v = isObj(it) ? it[key] : undefined;
      if (typeof v !== "string") return;
      if (has(first, v)) {
        var one = prefix.slice(0, -1), noun = what.split(" ").pop();
        out.push(issue([prefix, i, key], "duplicate", what + " " + repr(v) + " is used by both " + one + " " +
          (first[v] + 1) + " and " + one + " " + (i + 1) + "; give each " + one + " its own " + noun));
      } else first[v] = i;
    });
    return out;
  }

  var SCHEMAS = null;
  function schemas() {
    if (SCHEMAS) return SCHEMAS;
    var g = typeof self !== "undefined" ? self : (typeof global !== "undefined" ? global : {});
    if (!g.FLAIR_EVAL_SCHEMAS) throw new Error("schemas.js is not loaded");
    return g.FLAIR_EVAL_SCHEMAS;
  }
  function setSchemas(s) { SCHEMAS = s; }

  function validateSamples(rows, header) {
    var out = [];
    if (header) {
      header.forEach(function (c) {
        if (c && ALL_SAMPLE_COLUMNS.indexOf(c) < 0) {
          out.push(issue(["samples", "header", c], "unknown_key", "unknown column '" + c + "'; allowed columns are " +
            ALL_SAMPLE_COLUMNS.join(", ")));
        }
      });
    }
    if (!rows.length) out.push(issue(["samples"], "min_items", "needs at least one sample row"));
    var sch = schemas().samples;
    rows.forEach(function (r, i) {
      var iss = issuesFromSchema(r, sch, ["samples", i], "column");
      if (header) iss = iss.filter(function (x) { return x.code !== "unknown_key"; });
      out = out.concat(iss);
    });
    return out.concat(duplicates(rows, "sample", "samples", "sample name"));
  }

  function validateRuns(doc) {
    if (Array.isArray(doc)) return [issue(["runs"], "type", "runs.yaml should start with a 'runs:' key that holds the list of runs")];
    if (!isObj(doc)) return [issue(["runs"], "type", "runs.yaml should be a mapping with a 'runs:' list")];
    var out = issuesFromSchema(doc, schemas().runs, [], "key");
    if (Array.isArray(doc.runs)) out = out.concat(duplicates(doc.runs, "id", "runs", "run id"));
    return out;
  }

  /* ------------------------------------------------------------ command */
  var SAFE = /^[A-Za-z0-9_@%+=:,./-]+$/;
  function sq(s) { return s && SAFE.test(s) ? s : "'" + s.split("'").join("'\"'\"'") + "'"; }
  function quote(arg) {
    var s = String(arg);
    var m = /^(--[A-Za-z0-9_-]+=)([\s\S]*)$/.exec(s);
    return m ? m[1] + sq(m[2]) : sq(s);
  }
  function opt(name, value) {
    var v = String(value);
    return v.charAt(0) === "-" ? [name + "=" + v] : [name, v];
  }

  function quickForm(samples, runs) {
    var reasons = [];
    if (samples.length !== 1) reasons.push(samples.length + " samples (one-command mode takes exactly one)");
    else {
      var s0 = samples[0];
      var used = TRACK_COLUMNS.filter(function (c) { return !!s0[c]; });
      if (used.length) reasons.push("extra tracks (" + used.join(", ") + ") need samples.csv");
      if (s0.reads && s0.bam) reasons.push("the sample has both reads and bam");
    }
    var tokens = [];
    runs.forEach(function (r) {
      var rid = has(r, "id") ? r.id : "?";
      var tok = tokenForRun(r);
      if (tok === null) reasons.push("run " + rid + ": its id or FLAIR source has no --flair token");
      tokens.push(tok);
      if (r.label && r.label !== r.id) reasons.push("run " + rid + ": labels need runs.yaml");
      if (r.tags && Object.keys(r.tags).length) reasons.push("run " + rid + ": tags need runs.yaml");
      ["quantify", "productivity", "own_align"].forEach(function (k) {
        if (!runIsDefault(r, k)) reasons.push("run " + rid + ": " + k + " is not the default");
      });
    });
    var argSet = [];
    runs.forEach(function (r) {
      var a = (has(r, "transcriptome_args") ? r.transcriptome_args : "") || "";
      if (argSet.indexOf(a) < 0) argSet.push(a);
    });
    if (argSet.length > 1) reasons.push("runs use different transcriptome_args");
    if (!runs.length) reasons.push("no runs");
    if (reasons.length) return { args: null, reasons: reasons };

    var s = samples[0];
    var args = s.reads ? ["--reads", s.reads] : ["--bam", s.bam];
    args = args.concat(["--genome", s.genome, "--gtf", s.gtf, "--sample", s.sample]);
    ["platform", "library"].forEach(function (k) { if (s[k]) args.push("--" + k, s[k]); });
    args.push("--flair", tokens.join(","));
    if (argSet[0]) args = args.concat(opt("--flair_args", argSet[0]));
    return { args: args, reasons: [] };
  }

  function options(o) {
    o = o || {};
    return {
      ref: o.ref || DEFAULT_REF, profile: o.profile || null, outdir: o.outdir || PARAM_DEFAULTS.outdir,
      params: o.params || {}, resume: !!o.resume, work_dir: o.work_dir || null,
      samples_path: o.samples_path || "samples.csv", runs_path: o.runs_path || "runs.yaml",
      prefer_files: !!o.prefer_files, extra: o.extra || []
    };
  }
  function head(o) {
    var argv = ["nextflow", "run", PIPELINE, "-r", o.ref || DEFAULT_REF];
    if (o.profile) argv.push("-profile", o.profile);
    return argv;
  }
  function tail(o) {
    var args = [];
    PARAM_ORDER.forEach(function (k) {
      var v = o.params[k];
      if (v === undefined || v === null || v === "" || v === PARAM_DEFAULTS[k]) return;
      args = args.concat(opt("--" + k, v === true ? "true" : v));
    });
    args.push("--outdir", o.outdir || PARAM_DEFAULTS.outdir);
    if (o.work_dir) args.push("-work-dir", o.work_dir);
    if (o.resume) args.push("-resume");
    return args.concat(o.extra);
  }

  function buildCommand(samples, runs, opts) {
    var o = options(opts);
    var qf = quickForm(samples, runs);
    if (qf.args !== null && !o.prefer_files) {
      return { argv: head(o).concat(qf.args, tail(o)), mode: "quick", reasons: [] };
    }
    var reasons = o.prefer_files && qf.args !== null ? ["files requested"] : qf.reasons;
    return { argv: head(o).concat(["--input", o.samples_path, "--runs", o.runs_path], tail(o)), mode: "files", reasons: reasons };
  }

  var NO_VALUE = ["-resume"];
  var VALUE_FLAGS = ["-r", "-profile", "-work-dir"];
  function formatCommand(argv) {
    var qd = argv.map(quote);
    var lines = [qd.slice(0, 5).join(" ")];
    var i = 5;
    while (i < argv.length) {
      var a = argv[i];
      var takes = NO_VALUE.indexOf(a) < 0 && a.indexOf("=") < 0 && (a.indexOf("--") === 0 || VALUE_FLAGS.indexOf(a) >= 0);
      if (takes && i + 1 < argv.length) { lines.push("  " + qd[i] + " " + qd[i + 1]); i += 2; }
      else { lines.push("  " + qd[i]); i += 1; }
    }
    return lines.join(" \\\n");
  }
  function oneLine(argv) { return argv.map(quote).join(" "); }

  return {
    PIPELINE: PIPELINE, DEFAULT_REF: DEFAULT_REF, DEFAULT_GIT_URL: DEFAULT_GIT_URL, ID_RE: ID_RE,
    NAME_HINT: NAME_HINT, PRESETS: PRESETS, SAMPLE_COLUMNS: SAMPLE_COLUMNS, TRACK_COLUMNS: TRACK_COLUMNS,
    ALL_SAMPLE_COLUMNS: ALL_SAMPLE_COLUMNS, PARAM_DEFAULTS: PARAM_DEFAULTS, RUN_DEFAULTS: RUN_DEFAULTS,
    TokenError: TokenError, parseToken: parseToken, parseTokens: parseTokens,
    tokenForFlair: tokenForFlair, tokenForRun: tokenForRun,
    defaultSampleName: defaultSampleName, cleanSample: cleanSample, samplesCsv: samplesCsv, runsYaml: runsYaml,
    validateSamples: validateSamples, validateRuns: validateRuns, setSchemas: setSchemas,
    quote: quote, opt: opt, quickForm: quickForm, buildCommand: buildCommand,
    formatCommand: formatCommand, oneLine: oneLine
  };
});
