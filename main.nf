#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import java.util.concurrent.atomic.AtomicInteger
import nextflow.extension.CH
import nextflow.Channel

// Utility to normalise optional strings to empty strings
String resolveString(def value) {
    return value ? value.toString() : ''
}

// Derive a fallback publish directory name when none is provided in the config
String inferPublishName(step, idx) {
    def raw = step.cmd?.tokenize()?.first() ?: "step${idx}"
    def base = raw.replaceAll(/[^A-Za-z0-9._-]/, '_')
    return base ? "${base}_out" : "step${idx}_out"
}

// Replace human-friendly placeholders with actual paths prior to execution
String renderCommand(def cmdTemplate, String prevDir, String reuseDir) {
    def result = cmdTemplate?.toString() ?: ''
    result = result.replace('{{prev}}', prevDir ?: '')
    result = result.replace('{{reuse}}', reuseDir ?: '')
    result = result.replaceAll(/\$\{params\.([A-Za-z0-9_]+)\}/) { full, key ->
        def value = params.containsKey(key) ? params.get(key) : null
        return value != null ? value.toString() : ''
    }
    return result
}

String escapeForSingleQuotes(String value) {
    if( !value ) { return '' }
    return value.replace("'", "'\\''")
}

Map buildStepContext(Map plan, int idx, String prevDir) {
    def stepInfo = plan.steps[idx]
    return [
        caseName   : plan.caseName,
        publishName: stepInfo.publishName,
        step       : stepInfo.step,
        prevDir    : resolveString(prevDir),
        reuseDir   : resolveString(plan.reusePath ?: ''),
        stepIndex  : idx,
        resultsDir : plan.resultsDirStr,
        targetDir  : stepInfo.targetDir,
        targetParent: stepInfo.targetParent ?: plan.resultsDirStr
    ]
}

workflow {
    def infraProfiles = ['standard','slurm','singularity','conda','docker','local'] as Set
    def profileTokens = (workflow.profile ?: '').tokenize(',').collect { it.trim() }.findAll { it }
    def activeProfile = profileTokens.find { !(it in infraProfiles) }
    if( !activeProfile && profileTokens ) {
        activeProfile = profileTokens[0]
    }

    // Derive dataset namespace from active profile when user didn't override with --dataset
    if( (!params.dataset || params.dataset == 'default') && activeProfile ) {
        params.dataset = activeProfile
    }

    def resultsNamespace = activeProfile ?: params.dataset ?: 'default'

    def cases = params.test_cases ?: []
    if( !cases ) {
        log.warn "No test cases configured. Nothing to do."
        return
    }

    def casePlans = [:]
    def caseOutputs = [:]
    def completionSignals = []
    def totalCasesWithSteps = new AtomicInteger(0)

    cases.each { caseConfig ->
        if( !caseConfig?.name ) {
            throw new IllegalArgumentException("Each test case requires a unique 'name' (missing entry: ${caseConfig}).")
        }

        def caseName = caseConfig.name.toString()
        def rawSteps = (caseConfig.steps ?: []) as List
        def resultsDirPath = file("${params.results_dir ?: 'results'}/${resultsNamespace}/${caseName}").toAbsolutePath()
        def planSteps = rawSteps.withIndex().collect { step, idx ->
            if( !step.cmd ) {
                throw new IllegalArgumentException("Step ${idx} in case '${caseName}' is missing a 'cmd' entry.")
            }
            def publishName = step.publish ?: inferPublishName(step, idx)
            def resolvedStep = step instanceof Map ? new LinkedHashMap(step) : [:]

            def caseConda = step.containsKey('conda') && step.conda ? step.conda : caseConfig.conda
            def caseContainer = step.containsKey('container') && step.container ? step.container : caseConfig.container

            def runtimeTag = (step.runtime ?: step.tool ?: '').toString().toLowerCase()
            def wantsSqanti = runtimeTag == 'sqanti' || runtimeTag == 'sqanti3'

            if( wantsSqanti ) {
                if( params.sqanti_conda ) {
                    caseConda = params.sqanti_conda
                }
                if( params.sqanti_container ) {
                    caseContainer = params.sqanti_container
                }
            }

            if( caseConda ) {
                resolvedStep.conda = caseConda
            }
            else {
                resolvedStep.remove('conda')
            }

            if( caseContainer ) {
                resolvedStep.container = caseContainer
            }
            else {
                resolvedStep.remove('container')
            }

            def targetDirPath = resultsDirPath.resolve(publishName)
            [
                publishName: publishName,
                step       : resolvedStep,
                targetDir  : targetDirPath.toString(),
                targetParent: targetDirPath.parent.toString()
            ]
        }

        if( planSteps ) {
            totalCasesWithSteps.incrementAndGet()
        }

        def plan = [
            caseName     : caseName,
            reuseName    : caseConfig.reuse_from ?: null,
            steps        : planSteps,
            resultsDir   : resultsDirPath,
            resultsDirStr: resultsDirPath.toString()
        ]
        casePlans[caseName] = plan

        def resultChannel = CH.queue()
        caseOutputs[caseName] = resultChannel
        completionSignals << resultChannel.map { path -> [ name: caseName, path: path ] }
    }

    def stepQueue = CH.queue()
    def remainingCasesWithSteps = new AtomicInteger(totalCasesWithSteps.get())
    def stepResults = runStep(stepQueue).step_output

    stepResults.subscribe { ctx ->
        def plan = casePlans[ctx.caseName]
        def nextIndex = ctx.stepIndex + 1
        def prevDirForNext = ctx.targetDir ?: file("${plan.resultsDirStr}/${ctx.publishName}").toString()

        if( nextIndex < plan.steps.size() ) {
            stepQueue << buildStepContext(plan, nextIndex, prevDirForNext)
        }
        else {
            def resultChannel = caseOutputs[ctx.caseName]
            resultChannel << plan.resultsDirStr
            resultChannel << Channel.STOP

            if( remainingCasesWithSteps.decrementAndGet() == 0 ) {
                stepQueue << Channel.STOP
            }
        }
    }

    casePlans.each { caseName, plan ->
        def resultChannel = caseOutputs[caseName]
        def reuseName = plan.reuseName
        def reuseChannel = reuseName ? caseOutputs[reuseName] : null

        if( reuseName && !reuseChannel ) {
            throw new IllegalArgumentException("Test case '${caseName}' declares reuse_from='${reuseName}' but no case with that name exists.")
        }

        reuseChannel = reuseChannel ?: Channel.value('')

        reuseChannel.subscribe { reuseValue ->
            plan.reusePath = resolveString(reuseValue)
            if( plan.steps ) {
                stepQueue << buildStepContext(plan, 0, '')
            }
            else {
                resultChannel << plan.resultsDirStr
                resultChannel << Channel.STOP
            }
        }

    }

    if( remainingCasesWithSteps.get() == 0 ) {
        stepQueue << Channel.STOP
    }

    completionSignals.each { ch ->
        ch.view { data -> "[nextflow] completed case ${data.name} → ${data.path}" }
    }
}

process runStep {
    tag { "${ctx.caseName}:${ctx.publishName}" }

    input:
        val ctx

    cpus { ctx.step.cpus ?: 1 }
    memory { ctx.step.mem ?: '2 GB' }
    time { ctx.step.time ?: '1h' }
    conda { ctx.step.conda ?: null }
    container { ctx.step.container ?: null }

    output:
        val(ctx), emit: step_output

    script:
        def prevExport = ctx.prevDir ? "export PREV_OUT='${escapeForSingleQuotes(ctx.prevDir)}'" : 'unset PREV_OUT'
        def reuseExport = ctx.reuseDir ? "export REUSE_OUT='${escapeForSingleQuotes(ctx.reuseDir)}'" : 'unset REUSE_OUT'
        def cmd = renderCommand(ctx.step.cmd, ctx.prevDir, ctx.reuseDir)
        """
        set -euo pipefail

        ${prevExport}
        ${reuseExport}
        export STEP_NAME='${ctx.publishName}'
        export CASE_NAME='${ctx.caseName}'

        mkdir -p ${ctx.publishName}
        pushd ${ctx.publishName} >/dev/null

        ${cmd}

        popd >/dev/null

        link_src="\$(realpath \"${ctx.publishName}\")"
        mkdir -p "${ctx.resultsDir}"
        mkdir -p "${ctx.targetParent}"
        if [ -e "${ctx.targetDir}" ] || [ -L "${ctx.targetDir}" ]; then
            rm -rf "${ctx.targetDir}"
        fi
        ln -s "\${link_src}" "${ctx.targetDir}"

        if [ -d "${ctx.resultsDir}" ]; then
            for item in "\${link_src}"/*; do
                [ -f "\${item}" ] || continue
                base="\$(basename "\${item}")"
                dest="${ctx.resultsDir}/\${base}"
                if [ -e "\${dest}" ] || [ -L "\${dest}" ]; then
                    rm -f "\${dest}"
                fi
                ln -s "\${item}" "\${dest}"
            done
        fi
        """
}
