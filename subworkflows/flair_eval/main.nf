/*
 * FLAIR_EVAL: samples x runs -> one bundle per job (docs/DESIGN.md "Pipeline steps").
 *
 *   FLAIR_ENV           once per distinct FLAIR source, plus params.align_flair for ALIGN_SHARED
 *   PREP_REFERENCE      with --regions: genome + GTF cut to the regions' chromosomes
 *   ALIGN_SHARED        once per sample without a BAM (pinned FLAIR); skipped when a BAM is given
 *   ALIGN_OWN           per job, for runs with own_align: true (that run's FLAIR)
 *   ASSEMBLE            per job
 *   QUANTIFY            per job, when run.quantify and the sample has reads
 *   PRODUCTIVITY        per job, when run.productivity
 *   EVALUATE_STRUCTURE  per job (+ end references and TUSCO list when the sample has them, A26)
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
include { PREP_REFERENCE           } from '../../modules/prep_reference/main'
include { flairSourceKey; safeId   } from '../inputs/keys'

// Steps whose task cost is counted in bundle.json "compute": the FLAIR run itself.
// CONTRACTS section 7, A12: assemble, quantify, productivity (not align, not evaluation).
def COUNTED_STEPS() { ['assemble', 'quantify', 'productivity'] }

// Trace row columns of provenance/trace.tsv (A15).
def STEP_PROCESS() { [align: 'FLAIR_EVAL:ALIGN_OWN', assemble: 'FLAIR_EVAL:ASSEMBLE', quantify: 'FLAIR_EVAL:QUANTIFY', productivity: 'FLAIR_EVAL:PRODUCTIVITY', eval: 'FLAIR_EVAL:EVALUATE_STRUCTURE'] }

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
        runs.eachWithIndex { r, ri ->
            jobs << [id: "${s.sample}/${r.id}".toString(), sample: s.sample, run_id: r.id, run_index: ri, run: r, platform: s.platform, library: s.library]
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
            .map { _k, s, _json, env -> [[id: s.sample, sample: s.sample, platform: s.platform, library: s.library], s.reads, s.genome, s.gtf, env] }
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

    // ---- reference used by the FLAIR steps ---------------------------------------------------
    // With --regions, PREP_REFERENCE cuts genome + GTF to the chromosomes in the BED (same
    // coordinates). Alignment above keeps the full genome.
    def regions = params.regions ? file(params.regions, checkIfExists: true) : []
    def refs = samples.collect { s -> [refKey(s), s.genome, s.gtf] }.unique { r -> r[0] }
    // The genome's .fai is staged next to it when it exists: FLAIR dev quantify otherwise builds it
    // in several pool workers at once and fails with "sequence 'chrN' not present".
    if (params.regions) {
        PREP_REFERENCE(channel.fromList(refs), regions)
        ch_ref = PREP_REFERENCE.out.ref           // [key, genome, fai, gtf]
    }
    else {
        ch_ref = channel.fromList(refs).map { k, g, a -> [k, g, genomeIndex(g), a] }
    }
    ch_job_ref = channel.fromList(jobs)
        .map { j -> [refKey(S[j.sample]), j.id] }
        .combine(ch_ref, by: 0)
        .map { _k, id, g, fai, a -> [id, g, fai, a] }    // [id, genome, fai, gtf]

    // ---- assemble, quantify, productivity -----------------------------------------------------
    ASSEMBLE(
        ch_job_bam
            .join(ch_jobenv)
            .join(ch_job_ref)
            .map { id, bam, bai, _json, env, g, fai, a ->
                def j = J[id]
                [j, bam, bai, S[j.sample].reads ?: [], g, fai, a, env]
            },
        regions
    )
    ch_asm = ASSEMBLE.out.dir
        .map { j, d -> [j.id, d] }
        .join(ch_jobenv)                       // [id, assemble_dir, flair_env.json, env]
        .join(ch_job_bam)                      // + [bam, bai]
        .join(ch_job_ref)                      // + [genome, fai, gtf]

    // A6/A7: QUANTIFY runs for every quantify: true run. It gets reads when the sample has them,
    // plus the job's BAM, genome and isoforms; the adapter writes skipped.json if the FLAIR
    // version needs something that is missing.
    QUANTIFY(
        ch_asm
            .filter { id, _d, _json, _env, _bam, _bai, _g, _fai, _a -> J[id].run.quantify }
            .map { id, d, _json, env, bam, bai, g, fai, _a -> [J[id], d, S[J[id].sample].reads ?: [], bam, bai, g, fai, env] }
    )
    PRODUCTIVITY(
        ch_asm
            .filter { id, _d, _json, _env, _bam, _bai, _g, _fai, _a -> J[id].run.productivity }
            .map { id, d, _json, env, _bam, _bai, g, fai, a -> [J[id], d, g, fai, a, env] }
    )

    // ---- evaluate and bundle ------------------------------------------------------------------
    // A13: always the genome (junction motifs); reads/BAM paths are recorded in bundle inputs.
    // A26: the sample's end references and TUSCO list (explicit params or the genome preset);
    // a missing one is passed as [] and only turns its section off.
    EVALUATE_STRUCTURE(
        ch_asm.map { id, d, json, _env, _bam, _bai, g, fai, a ->
            def s = S[J[id].sample]
            [J[id], d, a, g, fai, json, s.ref_tss ?: [], s.ref_pas ?: [], s.tusco ?: [],
             [reads: s.reads?.toString(), bam: s.bam?.toString()]]
        },
        regions
    )

    ch_bundle_in = ASSEMBLE.out.dir.map { j, d -> [j.id, 'assemble', d] }
        .mix(
            QUANTIFY.out.dir.map { j, d -> [j.id, 'quantify', d] },
            PRODUCTIVITY.out.dir.map { j, d -> [j.id, 'productivity', d] },
            ALIGN_OWN.out.dir.map { j, d -> [j.id, 'align', d] },
            EVALUATE_STRUCTURE.out.dir.map { j, d -> [j.id, 'eval', d] },
        )
        .map { id, step, d -> [groupKey(id, nParts(J[id])), [step: step, dir: d]] }
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
            step       : p.step,
            counted    : p.step in COUNTED_STEPS(),
            task_id    : t.task_id,
            name       : "${STEP_PROCESS()[p.step]} (${j.id})".toString(),
            status     : t.exit == null ? 'UNKNOWN' : (t.exit == 0 ? 'COMPLETED' : 'FAILED'),
            exit       : t.exit,
            realtime   : t.realtime,
            '%cpu'     : t.pct_cpu,
            peak_rss   : t.peak_rss,
            rss        : t.rss,
            duration   : t.duration,
        ]
    }
    return [
        sample  : j.sample,
        run_id  : j.run_id,
        run_index: j.run_index,
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
            // A26: original paths of the end references / TUSCO list this job used (null = off)
            ref_tss : s.ref_tss?.toString(),
            ref_pas : s.ref_pas?.toString(),
            tusco   : s.tusco?.toString(),
            // with --regions FLAIR saw PREP_REFERENCE's copy cut to these chromosomes
            reference_cut_to_regions: params.regions ? true : false,
        ],
        // A24: launch = Nextflow session id, shared by every bundle of one launch (and kept on -resume)
        pipeline: [version: workflow.manifest.version ?: null, commit: workflow.commitId ?: null, launch: workflow.sessionId?.toString()],
        compute : computeFromTrace(rows.findAll { r -> r.counted }),
        trace   : rows,
    ]
}

/*
 * One task's numbers, read from its work dir (the step dir's parent):
 *   .command.trace  realtime (ms), %cpu (tenths of a percent), peak_rss/rss (KiB); written by
 *                   Nextflow because nextflow.config turns tracing on
 *   .exitcode       exit status
 *   .command.run / .exitcode mtimes  ~ submit and completion, for duration (ms)
 * Returned units: ms, percent, bytes. Missing files give nulls.
 */
def readTaskTrace(Object stepDir) {
    def wd = stepDir.parent
    def raw = [:]
    def f = wd.resolve('.command.trace')
    if (f.exists()) {
        f.readLines().each { l ->
            def i = l.indexOf('=')
            if (i > 0) {
                raw[l.substring(0, i)] = l.substring(i + 1)
            }
        }
    }
    def ec = wd.resolve('.exitcode')
    def run = wd.resolve('.command.run')
    def exit = ec.exists() && ec.text.trim().isInteger() ? ec.text.trim().toInteger() : null
    def duration = ec.exists() && run.exists() ? ec.lastModified() - run.lastModified() : null
    def cpu  = traceNum(raw, '%cpu')
    def prss = traceNum(raw, 'peak_rss')
    def rss  = traceNum(raw, 'rss')
    return [
        task_id : "${wd.parent.name}/${wd.name.take(6)}".toString(),   // work dir hash, as in the Nextflow log
        exit    : exit,
        realtime: traceNum(raw, 'realtime'),
        pct_cpu : cpu != null ? cpu / 10 : null,
        peak_rss: prss != null ? prss * 1024 : null,
        rss     : rss != null ? rss * 1024 : null,
        duration: duration != null && duration >= 0 ? duration : null,
    ]
}

def traceNum(Map raw, String k) {
    def v = raw[k]
    return v != null && v.isNumber() ? v.toBigDecimal() : null
}

// <genome>.fai next to the genome, or [] (nothing to stage).
def genomeIndex(Object genome) {
    def f = file(genome.toString() + '.fai')
    return f.exists() ? f : []
}

// Samples with the same genome + GTF share one PREP_REFERENCE task.
def refKey(Map s) {
    return "${s.genome}|${s.gtf}".toString()
}

// Parts BUNDLE waits for: assemble + eval, plus the optional steps this job runs.
def nParts(Map j) {
    return 2 + (j.run.quantify ? 1 : 0) + (j.run.productivity ? 1 : 0) + (j.run.own_align ? 1 : 0)
}

def computeFromTrace(List rows) {
    def timed = rows.findAll { r -> r.realtime != null }
    if (!timed) {
        return [cpu_hours: null, peak_rss_gb: null, wall_minutes: null]
    }
    def cpuH = timed.findAll { r -> r['%cpu'] != null }.sum(0) { r -> r.realtime / 3600000 * r['%cpu'] / 100 }
    def rss = timed.findAll { r -> r.peak_rss != null }.collect { r -> r.peak_rss }
    def wall = timed.sum(0) { r -> r.realtime } / 60000
    return [
        cpu_hours   : (cpuH as BigDecimal).setScale(4, BigDecimal.ROUND_HALF_UP).doubleValue(),
        peak_rss_gb : rss ? (rss.max() / (1024 * 1024 * 1024) as BigDecimal).setScale(3, BigDecimal.ROUND_HALF_UP).doubleValue() : null,
        wall_minutes: (wall as BigDecimal).setScale(2, BigDecimal.ROUND_HALF_UP).doubleValue(),
    ]
}
