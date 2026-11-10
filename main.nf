#!/usr/bin/env nextflow
/*
 * flair-eval v2: compare FLAIR runs (versions, branches or settings) on the same data and
 * write one self-contained HTML report.
 *
 *   nextflow run hdheath/flair-eval --input samples.csv --runs runs.yaml --outdir results
 *   nextflow run hdheath/flair-eval --reads r.fq.gz --genome g.fa --gtf a.gtf --flair 3.0.1,dev --outdir results
 *   nextflow run hdheath/flair-eval --report_only --bundles 'results/bundles' --outdir results
 *
 * Interfaces: docs/CONTRACTS.md. Design: docs/DESIGN.md. Outputs: docs/output.md.
 */

include { validateParameters; paramsSummaryLog; paramsHelp } from 'plugin/nf-schema'
include { loadInputs; parseBundlesParam        } from './subworkflows/inputs/main'
include { FLAIR_EVAL                           } from './subworkflows/flair_eval/main'
include { REPORT                               } from './modules/report/main'

workflow {
    if (params.help && params.help.toString() != 'false') {
        // --help lists all params; --help <param> explains one.
        def helpOpts = [command: 'nextflow run hdheath/flair-eval --reads reads.fq.gz --genome genome.fa --gtf annotation.gtf --flair 3.0.1,dev --outdir results']
        def topic = params.help.toString()
        log.info(topic == 'true' ? paramsHelp(helpOpts) : paramsHelp(helpOpts, topic))
    }
    else {
        checkFlairArgs(params.flair_args)
        validateParameters()
        log.info(paramsSummaryLog(workflow))
        writeParamsJson()

        def prior = parseBundlesParam(params.bundles)
        def ch_bundles
        def baseline = params.baseline ? params.baseline.toString() : ''

        if (params.report_only) {
            if (!prior) {
                error("--report_only needs --bundles: one or more comma-separated bundle directories or globs from earlier runs, e.g. --bundles 'results/bundles'.")
            }
            def ignored = ['input', 'runs', 'reads', 'bam', 'flair'].findAll { k -> params[k] }
            if (ignored) {
                log.warn("--report_only: ignoring ${ignored.collect { k -> '--' + k }.join(', ')}; no FLAIR jobs run.")
            }
            ch_bundles = channel.value([])
        }
        else {
            def inputs = loadInputs(params)
            def runIds = inputs.runs.collect { r -> r.id }
            log.info(
                "flair-eval: ${inputs.samples.size()} sample(s) x ${inputs.runs.size()} run(s) = ${inputs.samples.size() * inputs.runs.size()} job(s)\n" +
                "  samples: ${inputs.samples.collect { s -> s.sample }.join(', ')}\n" +
                "  runs:    ${runIds.join(', ')}"
            )
            if (params.baseline && !(params.baseline in runIds) && !prior) {
                log.warn("--baseline '${params.baseline}' is not one of this launch's runs (${runIds.join(', ')}) and no --bundles were given; the report will fall back to its default baseline.")
            }
            // A5: without --bundles the default baseline is this launch's first run.
            // A22: with --bundles and no --baseline, pass none; build_report.py picks the first
            // run of the earliest launch.
            baseline = baseline ?: (prior ? '' : runIds[0])
            // build_report.py orders runs as first seen, so hand it bundles in run order.
            def runIndex = runIds.withIndex().collectEntries { rid, i -> [(rid): i] }
            FLAIR_EVAL(inputs.samples, inputs.runs)
            ch_bundles = FLAIR_EVAL.out.bundles
                .toSortedList { a, b -> (runIndex[a[0].run_id] <=> runIndex[b[0].run_id]) ?: (a[0].sample <=> b[0].sample) }
                .map { items -> items.collect { item -> item[1] } }
        }

        REPORT(ch_bundles, prior, params.report_title, baseline)
    }
}

// A4: `--flair_args "--x"` makes Nextflow set flair_args=true and a stray --x param.
def checkFlairArgs(Object v) {
    if (v instanceof Boolean || (v != null && v.toString() == 'true')) {
        error("--flair_args has no value: Nextflow read the words after it as separate options. Quote it with '=', e.g. --flair_args='--stringent' (or --flair_args='--stringent --check_splice').")
    }
}

def writeParamsJson() {
    def dir = file("${params.outdir}/pipeline_info")
    dir.mkdirs()
    def p = params.collectEntries { k, v -> [(k): v instanceof CharSequence ? v.toString() : v] }
    dir.resolve('params.json').text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(p)) + '\n'
}
