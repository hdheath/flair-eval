/*
 * REPORT: one self-contained report.html over this launch's bundles plus --bundles
 * (bin/build_report.py, owned by the report stream; CONTRACTS section 4).
 */
process REPORT {
    tag "${title}"
    label 'process_single'

    input:
    path(bundles, stageAs: 'bundles/b?/*')
    path(prior,   stageAs: 'prior/p?/*')
    val(title)
    val(baseline)

    output:
    path('report.html'), emit: html

    script:
    def safeTitle = title.toString().replace("'", "'\\''")
    """
    build_report.py \\
        --bundles ${bundles} ${prior} \\
        --out report.html \\
        --title '${safeTitle}' \\
        ${baseline ? "--baseline ${baseline}" : ''}
    """

    stub:
    def safeTitle = title.toString().replace("'", "'\\''")
    """
    build_report.py --help > /dev/null
    {
        echo '<!doctype html><html><head><meta charset="utf-8"><title>stub report</title></head><body>'
        echo '<pre id="command">build_report.py --bundles ${bundles} ${prior} --out report.html --title ${safeTitle} ${baseline ? "--baseline ${baseline}" : ''}</pre>'
        echo '<h1>${safeTitle}</h1><ul>'
        find -L ${bundles} ${prior} -name bundle.json | sort | sed 's|.*|<li>&</li>|'
        echo '</ul></body></html>'
    } > report.html
    """
}
