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
        validateParameters()
        log.info(paramsSummaryLog(workflow))
        writeParamsJson()

        def prior = parseBundlesParam(params.bundles)
        def ch_bundles

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
            FLAIR_EVAL(inputs.samples, inputs.runs)
            ch_bundles = FLAIR_EVAL.out.bundles
                .toSortedList { a, b -> a[0].id <=> b[0].id }
                .map { items -> items.collect { item -> item[1] } }
        }

        REPORT(ch_bundles, prior, params.report_title, params.baseline ?: '')
    }
}

def writeParamsJson() {
    def dir = file("${params.outdir}/pipeline_info")
    dir.mkdirs()
    def p = params.collectEntries { k, v -> [(k): v instanceof CharSequence ? v.toString() : v] }
    dir.resolve('params.json').text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(p)) + '\n'
}
