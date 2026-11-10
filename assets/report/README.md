# Report assets

`build_report.py` inlines these three files into one `report.html` and embeds the data as JSON in
`<script id="report-data">`:

- `template.html`: page shell with `__TITLE__`, `/*__CSS__*/`, `/*__JS__*/` and `__DATA__` slots.
- `report.css`: type scale, layout, and light/dark colour tokens (dark follows the OS, or the
  Auto/Light/Dark switch).
- `report.js`: the views and charts.

Charts are inline SVG drawn by `report.js`. No third-party JavaScript or CSS is vendored, so
there is no licence file to ship, and the report makes no network requests.

## Views
Overview, Structure, Agreement, Expression, Runtime, Run detail, Guide. A view whose data is
absent (no `isoforms.bed`, no quantification, no `compute`) is not shown.

## Design rules
- **Colour means run.** Each run takes a categorical slot from its position in `DATA.runs`
  (baseline first) and keeps it in every view. Runs past the eighth are gray.
- Magnitude (heatmaps) and parts of one whole (end match, junction and coding classes) use a
  neutral gray ramp, so they are never mistaken for a run.
- Every chart has a Table toggle with its numbers. Marks are keyboard-focusable and show the
  same tooltip as on hover. The tabs follow the ARIA tabs pattern (arrow keys, Home and End).
- Print / PDF renders every view in turn with print styles.

## Data `build_report.py` adds
`by_sample.<sample>.concordance` and `.quant` (CONTRACTS section 8), computed from
`isoforms.bed`, `tables/isoforms.tsv`, `tables/genes.tsv` (optional, for gene names) and
`counts.tsv` beside each `bundle.json`. See `bin/flair_eval_lib/concordance.py`.
None of it depends on the baseline, so the reader can switch the baseline in the page (A22):
`concordance.genes_detected` lists, per run, the reference genes it detects among those not every
run detects (gained/lost against any run follow), and every `quant` pair carries a density grid
(`grid[j][i]`, i = run a, j = run b; the page transposes it when the baseline is b).

`baseline` and `baseline_rule` (`{"rule": "given" | "earliest_launch" | "name_order",
"launches": n}`) record which run is the default baseline and why. The page shows it under the
Overview heading and in the subtitle, and offers "Use default" after a switch.

## Later sections
Register a panel in `SECTION_PANELS` in `report.js`, keyed by the `bundle.sections.<key>` it reads.
Unknown section keys are ignored.
