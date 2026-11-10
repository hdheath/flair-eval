/*
 * FLAIR_EVAL: samples x runs -> one bundle per job (docs/DESIGN.md "Pipeline steps").
 *
 *   FLAIR_ENV           once per distinct FLAIR source, plus params.align_flair for ALIGN_SHARED
 *   ALIGN_SHARED        once per sample without a BAM (pinned FLAIR); skipped when a BAM is given
 *   ALIGN_OWN           per job, for runs with own_align: true (that run's FLAIR)
 *   ASSEMBLE            per job
 *   QUANTIFY            per job, when run.quantify and the sample has reads
 *   PRODUCTIVITY        per job, when run.productivity
 *   EVALUATE_STRUCTURE  per job
 *   BUNDLE              per job -> <outdir>/bundles/<sample>/<run_id>/
 *
 * Channel items carry a small job map ("meta"): id = "<sample>/<run_id>", sample, run_id, run.
 * Sample files are looked up by name, so no wide tuples are threaded through the graph.
 */

include { FLAIR_ENV                } from '../../modules/flair_env/main'
include { ALIGN as ALIGN_SHARED    } from '../../modules/align/main'
include { ALIGN as ALIGN_OWN       } from '../../modules/align/main'
include { ASSEMBLE                 } from '../../modules/assemble/main'
include { QUANTIFY                 } from '../../modules/quantify/main'
include { PRODUCTIVITY             } from '../../modules/productivity/main'
include { EVALUATE_STRUCTURE       } from '../../modules/evaluate_structure/main'
include { BUNDLE                   } from '../../modules/bundle/main'
include { flairSourceKey; safeId   } from '../inputs/keys'

// Steps whose task cost is counted in bundle.json "compute": the FLAIR run itself.
def COUNTED_STEPS() { ['align', 'assemble', 'quantify', 'productivity'] }

workflow FLAIR_EVAL {
    take:
    samples   // List of sample maps from loadInputs()
    runs      // List of run maps from loadInputs()

    main:
    def S = samples.collectEntries { s -> [(s.sample): s] }

    // ---- FLAIR environments: one FLAIR_ENV task per distinct source --------------------------
    def envRuns = [:]
    runs.each { r ->
        def k = flairSourceKey(r.flair)
        if (!envRuns.containsKey(k)) {
            envRuns[k] = [id: r.id, flair: r.flair]
        }
    }
    def sharedSamples = runs.any { r -> !r.own_align } ? samples.findAll { s -> !s.bam } : []
    def alignFlair = [release: params.align_flair.toString()]
    def alignKey = flairSourceKey(alignFlair)
    if (sharedSamples && !envRuns.containsKey(alignKey)) {
        envRuns[alignKey] = [id: "align-flair-${safeId(params.align_flair.toString())}".toString(), flair: alignFlair]
    }
    def keyById = envRuns.collectEntries { k, r -> [(r.id): k] }

    FLAIR_ENV(channel.fromList(envRuns.values() as List))

    // [source_key, flair_env.json, slim env map]. Only the fields tasks need go into the map, so
    // the task hashes downstream stay stable on -resume.
    ch_env = FLAIR_ENV.out.map { id, json ->
        def e = new groovy.json.JsonSlurper().parseText(json.text)
        [keyById[id], json, [
            conda_prefix   : e.conda_prefix ?: null,
            container      : e.container ?: null,
            flair_version  : e.flair_version ?: null,
            resolved_commit: e.resolved_commit ?: null,
        ]]
    }

    // ---- jobs = samples x runs --------------------------------------------------------------
    def jobs = []
    samples.each { s ->
        runs.each { r ->
            jobs << [id: "${s.sample}/${r.id}".toString(), sample: s.sample, run_id: r.id, run: r]
        }
    }
    def J = jobs.collectEntries { j -> [(j.id): j] }

    // [job, flair_env.json, env map]
    ch_jobs = channel.fromList(jobs)
        .map { j -> [flairSourceKey(j.run.flair), j] }
        .combine(ch_env, by: 0)
        .map { _k, j, json, env -> [j, json, env] }
    ch_jobenv = ch_jobs.map { j, json, env -> [j.id, json, env] }

    // ---- alignment ---------------------------------------------------------------------------
    ALIGN_SHARED(
        channel.fromList(sharedSamples)
            .map { s -> [alignKey, s] }
            .combine(ch_env, by: 0)
            .map { _k, s, _json, env -> [[id: s.sample, sample: s.sample], s.reads, s.genome, s.gtf, env] }
    )
    ch_sample_bam = ALIGN_SHARED.out.bam
        .map { m, bam, bai -> [m.sample, bam, bai] }
        .mix(channel.fromList(samples.findAll { s -> s.bam }).map { s -> [s.sample, s.bam, s.bai] })

    ALIGN_OWN(
        ch_jobs
            .filter { j, _json, _env -> j.run.own_align }
            .map { j, _json, env -> [j, S[j.sample].reads, S[j.sample].genome, S[j.sample].gtf, env] }
    )

    ch_job_bam = ch_jobs
        .filter { j, _json, _env -> !j.run.own_align }
        .map { j, _json, _env -> [j.sample, j.id] }
        .combine(ch_sample_bam, by: 0)
        .map { _sname, id, bam, bai -> [id, bam, bai] }
        .mix(ALIGN_OWN.out.bam.map { j, bam, bai -> [j.id, bam, bai] })

    // ---- assemble, quantify, productivity -----------------------------------------------------
    def regions = params.regions ? file(params.regions, checkIfExists: true) : []

    ASSEMBLE(
        ch_job_bam
            .join(ch_jobenv)
            .map { id, bam, bai, _json, env ->
                def j = J[id]
                def s = S[j.sample]
                [j, bam, bai, s.reads ?: [], s.genome, s.gtf, env]
            },
        regions
    )
    ch_asm = ASSEMBLE.out.dir
        .map { j, d -> [j.id, d] }
        .join(ch_jobenv)                       // [id, assemble_dir, flair_env.json, env]

    QUANTIFY(
        ch_asm
            .filter { id, _d, _json, _env -> J[id].run.quantify && S[J[id].sample].reads }
            .map { id, d, _json, env -> [J[id], d, S[J[id].sample].reads, env] }
    )
    PRODUCTIVITY(
        ch_asm
            .filter { id, _d, _json, _env -> J[id].run.productivity }
            .map { id, d, _json, env -> [J[id], d, S[J[id].sample].genome, S[J[id].sample].gtf, env] }
    )

    // ---- evaluate and bundle ------------------------------------------------------------------
    EVALUATE_STRUCTURE(
        ch_asm.map { id, d, json, _env -> [J[id], d, S[J[id].sample].gtf, json] },
        regions
    )

    ch_bundle_in = ASSEMBLE.out.dir.map { j, d -> [j.id, 'assemble', d] }
        .mix(
            QUANTIFY.out.dir.map { j, d -> [j.id, 'quantify', d] },
            PRODUCTIVITY.out.dir.map { j, d -> [j.id, 'productivity', d] },
            ALIGN_OWN.out.dir.map { j, d -> [j.id, 'align', d] },
            EVALUATE_STRUCTURE.out.dir.map { j, d -> [j.id, 'eval', d] },
        )
        .map { id, step, d -> [groupKey(id, nParts(J[id], S[J[id].sample])), [step: step, dir: d]] }
        .groupTuple()
        .map { gk, parts -> [gk.getGroupTarget().toString(), parts] }
        .join(ch_jobenv)
        .map { id, parts, json, _env ->
            def j = J[id]
            def sorted = parts.sort { p -> p.step }
            [j, json, sorted.collect { p -> p.dir }, bundleInfo(j, S[j.sample], sorted)]
        }

    BUNDLE(ch_bundle_in)

    emit:
    bundles = BUNDLE.out.bundle   // [job, bundle_dir]
}

/*
 * info.json for BUNDLE: what only the pipeline knows (inputs, pipeline version, compute).
 */
def bundleInfo(Map j, Map s, List parts) {
    def rows = parts.collect { p ->
        def t = readTaskTrace(p.dir)
        [
            step           : p.step,
            counted        : p.step in COUNTED_STEPS(),
            realtime_ms    : t.realtime,
            pct_cpu        : t.pct_cpu,
            peak_rss_bytes : t.peak_rss,
            peak_vmem_bytes: t.peak_vmem,
            rchar          : t.rchar,
            wchar          : t.wchar,
        ]
    }
    return [
        sample  : j.sample,
        run_id  : j.run_id,
        label   : j.run.label,
        tags    : j.run.tags,
        inputs  : [
            genome  : s.genome?.toString(),
            gtf     : s.gtf?.toString(),
            reads   : s.reads?.toString(),
            bam     : s.bam?.toString(),
            regions : params.regions ? file(params.regions).toString() : null,
            platform: s.platform,
            library : s.library,
        ],
        pipeline: [version: workflow.manifest.version ?: null, commit: workflow.commitId ?: null],
        compute : computeFromTrace(rows.findAll { r -> r.counted }),
        trace   : rows,
    ]
}

/*
 * Read <task workdir>/.command.trace (written by Nextflow when tracing is on; nextflow.config
 * turns it on). Returns a map of numbers, empty when the file is absent.
 * Units: realtime ms, pct_cpu percent, peak_rss/peak_vmem bytes.
 */
def readTaskTrace(Object stepDir) {
    def f = stepDir.parent.resolve('.command.trace')
    if (!f.exists()) {
        return [:]
    }
    def raw = [:]
    f.readLines().each { l ->
        def i = l.indexOf('=')
        if (i > 0) {
            raw[l.substring(0, i)] = l.substring(i + 1)
        }
    }
    def cpu  = traceNum(raw, '%cpu')
    def rss  = traceNum(raw, 'peak_rss')
    def vmem = traceNum(raw, 'peak_vmem')
    return [
        realtime : traceNum(raw, 'realtime'),
        pct_cpu  : cpu != null ? cpu / 10 : null,       // stored as tenths of a percent
        peak_rss : rss != null ? rss * 1024 : null,     // stored in KiB
        peak_vmem: vmem != null ? vmem * 1024 : null,
        rchar    : traceNum(raw, 'rchar'),
        wchar    : traceNum(raw, 'wchar'),
    ]
}

def traceNum(Map raw, String k) {
    def v = raw[k]
    return v != null && v.isNumber() ? v.toBigDecimal() : null
}

// Parts BUNDLE waits for: assemble + eval, plus the optional steps this job runs.
def nParts(Map j, Map s) {
    return 2 + (j.run.quantify && s.reads ? 1 : 0) + (j.run.productivity ? 1 : 0) + (j.run.own_align ? 1 : 0)
}

def computeFromTrace(List rows) {
    def timed = rows.findAll { r -> r.realtime_ms != null }
    if (!timed) {
        return [cpu_hours: null, peak_rss_gb: null, wall_minutes: null]
    }
    def cpuH = timed.findAll { r -> r.pct_cpu != null }.sum(0) { r -> r.realtime_ms / 3600000 * r.pct_cpu / 100 }
    def rss = timed.findAll { r -> r.peak_rss_bytes != null }.collect { r -> r.peak_rss_bytes }
    def wall = timed.sum(0) { r -> r.realtime_ms } / 60000
    return [
        cpu_hours   : (cpuH as BigDecimal).setScale(4, BigDecimal.ROUND_HALF_UP).doubleValue(),
        peak_rss_gb : rss ? (rss.max() / (1024 * 1024 * 1024) as BigDecimal).setScale(3, BigDecimal.ROUND_HALF_UP).doubleValue() : null,
        wall_minutes: (wall as BigDecimal).setScale(2, BigDecimal.ROUND_HALF_UP).doubleValue(),
    ]
}
