// Module: Preflight Argument Validation
// Validates all assembler arg strings before any real work starts.
// Catches unrecognized flags early so the pipeline fails fast with a clear
// error message rather than discovering bad args after hours of compute.

// ---------------------------------------------------------------------------
// FLAIR transcriptome
// Runs: flair transcriptome -b /dev/null -g /dev/null <args>
// A non-zero exit (code 2) or the text "unrecognized arguments" in stderr
// means the arg string is invalid.
// ---------------------------------------------------------------------------
process PreflightValidateFlair {
    tag "preflight:flair"
    conda params.conda_flair
    executor 'local'
    // Never cache — always re-validate when the pipeline is invoked.
    cache false

    input:
    val(entries)   // list of [condition_name, transcriptome_args] pairs

    output:
    val(true)

    script:
    // Build a bash heredoc table of condition_name<TAB>args to iterate.
    def lines = entries.collect { e -> "${e[0]}\t${e[1]}" }.join('\n')
    """
    set -euo pipefail

    FAIL=0
    while IFS=\$'\t' read -r COND ARGS; do
        # Strip flags that require file paths so we don't need real files
        CLEANED_ARGS=\$(echo "\$ARGS" \\
            | sed 's|--junction_tab[[:space:]]*[^[:space:]]*||g' \\
            | sed 's|--tss_model[[:space:]]*[^[:space:]]*||g' \\
            | sed 's|--tts_model[[:space:]]*[^[:space:]]*||g')

        # Run flair with fake required args; capture stderr
        STDERR=\$(flair transcriptome -b /dev/null -g /dev/null \${CLEANED_ARGS} 2>&1 || true)

        if echo "\$STDERR" | grep -q "unrecognized arguments"; then
            UNRECOGNIZED=\$(echo "\$STDERR" | grep "unrecognized arguments")
            echo "PREFLIGHT ERROR [flair / \$COND]: \${UNRECOGNIZED}" >&2
            echo "  Full args were: \$ARGS" >&2
            FAIL=1
        else
            echo "PREFLIGHT OK [flair / \$COND]"
        fi
    done <<'ENTRIES'
${lines}
ENTRIES

    if [ "\$FAIL" -ne 0 ]; then
        exit 2
    fi
    """
}

// ---------------------------------------------------------------------------
// IsoQuant
// Runs: isoquant.py --reference /dev/null --genedb /dev/null
//                   --bam /dev/null --output /tmp/iq_preflight --data_type nanopore
//                   <extra_args>
// IsoQuant uses argparse; unrecognized args produce exit code 2.
// ---------------------------------------------------------------------------
process PreflightValidateIsoquant {
    tag "preflight:isoquant"
    conda params.conda_isoquant
    executor 'local'
    cache false

    input:
    val(entries)   // list of [condition_name, isoquant_args] pairs

    output:
    val(true)

    script:
    def lines = entries.collect { e -> "${e[0]}\t${e[1]}" }.join('\n')
    """
    set -euo pipefail

    FAIL=0
    while IFS=\$'\t' read -r COND ARGS; do
        OUTDIR=\$(mktemp -d)

        STDERR=\$(isoquant.py \\
            --reference /dev/null \\
            --genedb /dev/null \\
            --bam /dev/null \\
            --output "\${OUTDIR}" \\
            --data_type nanopore \\
            \${ARGS} 2>&1 || true)

        rm -rf "\${OUTDIR}"

        if echo "\$STDERR" | grep -qiE "unrecognized arguments|error: unrecognized"; then
            UNRECOGNIZED=\$(echo "\$STDERR" | grep -iE "unrecognized arguments|error: unrecognized" | head -1)
            echo "PREFLIGHT ERROR [isoquant / \$COND]: \${UNRECOGNIZED}" >&2
            echo "  Full args were: \$ARGS" >&2
            FAIL=1
        else
            echo "PREFLIGHT OK [isoquant / \$COND]"
        fi
    done <<'ENTRIES'
${lines}
ENTRIES

    if [ "\$FAIL" -ne 0 ]; then
        exit 2
    fi
    """
}
