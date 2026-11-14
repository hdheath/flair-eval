# Report assets

`build_report.py` inlines these three files into one `report.html` and embeds the data as JSON in
`<script id="report-data">`:

- `template.html`: page shell with `__TITLE__`, `/*__CSS__*/`, `/*__JS__*/` and `__DATA__` slots.
- `report.css`: type scale, layout, and light/dark colour tokens (dark follows the OS, or the
  Auto/Light/Dark switch).
- `report.js`: the transcriptome viewer (CONTRACTS section 11). What each tab means for a
  reader is in [docs/output.md](../../docs/output.md#reading-the-report); this file is for
  developers.

Charts are inline SVG drawn by `report.js`. No third-party JavaScript or CSS is vendored, so
there is no licence file to ship, and the report makes no network requests.

## Tabs
One sample at a time: the Sample picker drives every tab, and nothing is pooled across samples.
`#<tab>` or `#<tab>/<sample>` in the URL opens a tab directly.

| Tab | Shows | Reads |
|---|---|---|
| What differs | "Where the runs differ" (headline metrics as dot rows, each on its own min–max scale); per run: tiles against the baseline and its biggest differences (at most one per family); "All metrics" (the full table, closed by default); "Terms" | the metric registry (`top: true` rows feed the dot rows) |
| Structure | categories (with one-line class definitions), isoform counts, isoform length and exons histograms | bundle summary, `histograms` |
| Recovery | TUSCO tiles (sensitivity in genes, precision in isoforms) and Found/Partial/Wrong isoform/Missed counts; genes detected; genes gained and lost; FSM isoforms with exact ends (CIM); lost genes by baseline expression decile | `sections.tusco`, `by_sample.concordance`, `by_sample.lost_deciles` |
| Diversity | isoforms per gene; coding and NMD classes; how genes vary (with "neither" when a run has any) | `histograms.isoforms_per_gene`, `sections.productivity`, `sections.variation` |
| Agreement | isoforms shared between runs (Jaccard matrix); shared and unique isoforms vs baseline; expression density + Pearson r; 5′ and 3′ end cluster Jaccard; TUSCO genes vs baseline | `by_sample.concordance`, `.quant`, `.end_agreement`, `.tusco_incremental` |
| Ends | precision at 50 bp (one table); precision as the window widens; end precision by category; signed distance to refTSS and PolyA_DB; FSM end distance to the reference; polyA motif; then "3′ ends on genomic A-runs": downstream-A histogram and internal-priming risk (drawn even without end databases) | `sections.ends`, `sections.priming` |
| Quantification (id `reads`) | read-count note once at the top (a callout only when runs counted reads differently, A21); reads assigned and isoforms with reads; isoform CPM | `sections.quant` |
| Compute | wall time, CPU time and peak memory of the FLAIR tasks | `compute` |

End databases give precision only: the report never shows recall against refTSS or PolyA_DB.

## Rules
- **Colour means run.** A run's `--series-N` slot is its position in `DATA.runs` (A17 order),
  fixed for the page, so colours do not move when the reader switches the baseline. Runs past
  the eighth, shared structures and "other" are grey. Magnitude (Jaccard matrix, CPM density) and
  parts of a whole (how genes vary) use the neutral grey ramp.
- **Baseline.** Drawn first; in line and dot charts it is the larger, outlined dot. The control
  bar switches it; the line under the bar says which rule picked the default (A22) and offers
  "Use default" after a switch. Every comparison is "vs baseline".
- **Overlay or side by side.** With two or more runs, charts overlay by default; the Layout
  toggle draws one panel per run on a shared scale.
- **Hidden panels.** A panel whose data is missing is replaced by a dashed card with a one-line
  reason in the reader's terms: the pipeline param and the preset that fills it ("needs refTSS
  sites (--genome hg38, or --ref_tss)"), never section keys or contract numbers. When only some
  runs have the data, the panel is drawn for those and a line names the others. Bundles from
  before round 3 render with those panels hidden ("not in these results (made by an older
  flair-eval)").
- **One value per run** is drawn with `runBars`: one labelled bar per run, so the reader does not
  need the colour key (Compute, genes detected, reads assigned, internal-priming risk, polyA
  signal).
- Card grids use `repeat(auto-fill, …)`; stat tiles are a two-column grid; the page works at
  phone width. Every chart has a Table toggle with its numbers.
- Print / PDF renders every tab in turn.

## Metric registry
`family(name, tab, rows)` in `report.js` registers a family's rows for "All metrics" and the
"Biggest differences" lists. A row is `{name, k, get(run), lowerIsWorse, needs, vsBase}`:
- `k` is the format (`n`, `p` percent, `f2`, `f3`, `wall`, `cpu`, `mem`);
- `lowerIsWorse` is `true` (higher is better), `false` (lower is better) or `null` (no
  direction); it labels a difference "better" or "worse";
- `top: true` also puts the row in "Where the runs differ" (never for `vsBase` rows, which the
  baseline lacks);
- `needs` is the one-line reason shown when no run has the row's data;
- `vsBase` rows compare a run with the baseline and are empty for the baseline itself;
- `noRank: true` keeps a row out of "Biggest differences" (it repeats another row's signal);
- `sameQuant: true` ranks a row only when the run and the baseline counted reads the same way.

"Biggest differences" takes at most one row per family, so one signal cannot fill the list.

A new section adds its rows to the registry and its panels to a view; section keys the page does
not know are ignored.

## Data `build_report.py` adds
Per sample, in `by_sample.<sample>` (none of it depends on the baseline, so the page can switch it):
- `concordance`, `quant` (CONTRACTS section 8), from `isoforms.bed`, `tables/isoforms.tsv`,
  `tables/genes.tsv` and `counts.tsv` beside each `bundle.json`
  (`bin/flair_eval_lib/concordance.py`). Every `quant` pair carries a density grid (`grid[j][i]`,
  i = run a); the page transposes it when the baseline is b.
- `end_agreement`, `tusco_incremental`, `lost_deciles` (A33), from `bin/flair_eval_lib/crossrun.py`:
  TSS/TTS end-cluster Jaccard for every run pair; TUSCO genes gained / both / lost for every
  ordered (baseline, run) pair; and, for every baseline with `counts.tsv`, its detected genes in
  ten CPM deciles with each other run's missing genes per decile.

`baseline` and `baseline_rule` (`{"rule": "given" | "earliest_launch" | "name_order",
"launches": n}`) record which run is the default baseline and why.

`window.flairEvalReport` exposes `{state, render, views}` so tests can drive every view in node
(`tests/python/test_report_viewer.py`).
