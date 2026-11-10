# Report assets

`build_report.py` inlines these three files into one `report.html` and embeds the bundle data as
JSON in `<script id="report-data">`:

- `template.html`: page shell with `__TITLE__`, `/*__CSS__*/`, `/*__JS__*/` and `__DATA__` slots.
- `report.css`: light and dark colour tokens (dark follows the OS, or the Theme button).
- `report.js`: views (overview, single run, compare, help) and the charts.

Charts are drawn as inline SVG by `report.js` itself. No third-party JavaScript or CSS is
vendored, so there is no licence file to ship. The report makes no network requests.

Later report sections: register a panel in `SECTION_PANELS` in `report.js`, keyed by the
`bundle.sections.<key>` it reads. Panels are drawn only when their section is present, and
unknown section keys are ignored.
