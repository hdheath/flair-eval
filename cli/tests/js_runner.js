// Runs frontend/builder/logic.js on a batch of cases for test_builder_parity.py.
// usage: node js_runner.js <logic.js> <schemas.js>  (cases JSON on stdin, results on stdout)
"use strict";
var FE = require(process.argv[2]);
FE.setSchemas(require(process.argv[3]));
var input = "";
process.stdin.on("data", function (d) { input += d; });
process.stdin.on("end", function () {
  var c = JSON.parse(input);
  var out = {};
  out.tokens = c.tokens.map(function (t) {
    try { var r = FE.parseToken(t); return { id: r.id, flair: r.flair }; }
    catch (e) { if (e instanceof FE.TokenError) return { error: true }; throw e; }
  });
  out.token_lists = c.token_lists.map(function (t) {
    try { return FE.parseTokens(t).map(function (r) { return r.id; }); }
    catch (e) { if (e instanceof FE.TokenError) return { error: true }; throw e; }
  });
  out.token_for_flair = c.flairs.map(FE.tokenForFlair);
  out.sample_names = c.paths.map(FE.defaultSampleName);
  out.validate_samples = c.sample_sets.map(function (s) {
    return FE.validateSamples(s.rows, s.header || null).map(function (i) { return [i.path, i.code, i.message]; });
  });
  out.validate_runs = c.run_docs.map(function (d) {
    return FE.validateRuns(d).map(function (i) { return [i.path, i.code, i.message]; });
  });
  out.commands = c.commands.map(function (k) {
    var r = FE.buildCommand(k.samples, k.runs, k.opts);
    return { argv: r.argv, mode: r.mode, reasons: r.reasons, pretty: FE.formatCommand(r.argv),
      line: FE.oneLine(r.argv), cli: FE.formatCommand(FE.asCliRun(r.argv)), csv: FE.samplesCsv(k.samples), yaml: FE.runsYaml(k.runs) };
  });
  out.presets = FE.PRESETS;
  out.meta = { not_exposed: FE.NOT_EXPOSED, launch_options: FE.LAUNCH_OPTIONS,
    param_defaults: FE.PARAM_DEFAULTS, param_order: FE.PARAM_ORDER,
    genome_presets: FE.GENOME_PRESETS, outdir_placeholder: FE.OUTDIR_PLACEHOLDER,
    outdir_example: FE.OUTDIR_EXAMPLE };
  process.stdout.write(JSON.stringify(out));
});
