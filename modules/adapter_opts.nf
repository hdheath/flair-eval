/*
 * Optional flair_adapter.py flags (CONTRACTS section 7: A8 --flair-bed, A9 --platform/--library).
 *
 * adapterOptFlags(step, opts) returns bash that fills the array OPT with each flag the installed
 * adapter lists in `flair_adapter.py <step> --help`, so core works with adapters before and after
 * they add a flag. opts: list of [flag, value] or [flag, value, 'file'] (only when that file exists
 * and is non-empty). Null or empty values are skipped. Use it as:  ... ${OPT[@]+"${OPT[@]}"}
 */
def adapterOptFlags(String step, List opts) {
    def lines = ['OPT=()', "ADAPTER_HELP=\$(flair_adapter.py ${step} --help 2>/dev/null || true)".toString()]
    opts.findAll { o -> o[1] != null && o[1].toString() != '' }.each { o ->
        def flag = o[0]
        def value = o[1].toString().replace("'", "'\\''")
        def cond = "grep -q -e '${flag}' <<< \"\$ADAPTER_HELP\""
        if (o.size() > 2 && o[2] == 'file') {
            cond = "[ -s '${value}' ] && " + cond
        }
        lines << "if ${cond}; then OPT+=(${flag} '${value}'); fi".toString()
    }
    return lines.join('\n')
}
