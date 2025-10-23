#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import groovy.json.JsonOutput
import groovy.util.ConfigSlurper
import java.io.File
import java.nio.file.Path
import java.util.regex.Matcher
import nextflow.Channel

List<String> normalizeList(def value) {
    if (!value) {
        return []
    }
    if (value instanceof Collection) {
        return value.collect { it?.toString()?.trim() }
            .findAll { it }
    }
    if (value instanceof CharSequence) {
        def text = value.toString().trim()
        if (!text) {
            return []
        }
        return text.split(',')
            .collect { it.trim() }
            .findAll { it }
    }
    return [value.toString()]
}

List<String> normalizeCommands(def value) {
    if (!value) {
        return []
    }
    if (value instanceof Collection) {
        return value.collect { it?.toString()?.trim() }
            .findAll { it }
    }
    def text = value.toString()
    return text
        .split(/[\n;]/)
        .collect { it.trim() }
        .findAll { it }
}

List<Path> flattenPathList(def value) {
    if (!value) {
        return []
    }
    if (value instanceof Path) {
        return [value as Path]
    }
    if (value instanceof File) {
        return [((File) value).toPath()]
    }
    if (value instanceof Collection) {
        List<Path> result = []
        value.each { result.addAll(flattenPathList(it)) }
        return result
    }
    return []
}

String escapeForSingleQuotes(String value) {
    if (!value) {
        return ''
    }
    return value.replace("'", "'\\''")
}

Map<String, Object> parseRegionalizeCommand(String text) {
    def raw = text?.toString()?.trim()
    if (!raw) {
        throw new IllegalArgumentException("Regionalize command cannot be empty")
    }
    if (raw.equalsIgnoreCase('--all') || raw.equalsIgnoreCase('all')) {
        return [
            mode         : 'all',
            command_text : raw,
            region_string: null
        ]
    }
    def matcher = (raw =~ /(?i)^--region\s+([A-Za-z0-9._-]+):(\d+)-(\d+)$/)
    if (matcher.matches()) {
        def chrom = matcher.group(1)
        long start = matcher.group(2).toLong()
        long end = matcher.group(3).toLong()
        if (start > end) {
            long tmp = start
            start = end
            end = tmp
        }
        String regionString = "${chrom}:${start}-${end}"
        return [
            mode         : 'region',
            chrom        : chrom,
            start        : start,
            end          : end,
            region_string: regionString,
            command_text : raw
        ]
    }
    throw new IllegalArgumentException("Unsupported regionalize command '${raw}'. Expected '--all' or '--region chr:start-end'.")
}

String renderFlairCommand(String template, Map<String, Object> spec, String datasetName) {
    if (!template) {
        throw new IllegalArgumentException("Flair command template cannot be empty")
    }

    Map<String, String> substitutions = new LinkedHashMap<>()
    spec.each { key, value ->
        substitutions[key.toString()] = value?.toString() ?: ''
    }
    substitutions['dataset'] = datasetName

    Set<String> unknownKeys = new LinkedHashSet<>()
    def matcher = (template =~ /\{([A-Za-z0-9_]+)\}/)
    matcher.each { full, key ->
        if (!substitutions.containsKey(key)) {
            unknownKeys << key
        }
    }
    if (unknownKeys) {
        throw new IllegalArgumentException("Command '${template}' references unknown parameter(s): ${unknownKeys.join(', ')}")
    }

    String rendered = template
    substitutions.each { key, value ->
        rendered = rendered.replaceAll("\\{${key}\\}", Matcher.quoteReplacement(value ?: ''))
    }
    return rendered
}

Path locateDatasetsConfig() {
    if (params.datasets_config) {
        def provided = file(params.datasets_config)
        if (!provided.exists()) {
            throw new IllegalArgumentException("Provided --datasets_config '${params.datasets_config}' does not exist")
        }
        return provided
    }

    def underConfigs = file("${projectDir}/configs/datasets.config")
    if (underConfigs.exists()) {
        return underConfigs
    }

    def atRoot = file("${projectDir}/datasets.config")
    if (atRoot.exists()) {
        return atRoot
    }

    throw new IllegalStateException("Unable to locate datasets.config; provide explicit --datasets_config path")
}

Map<String, Map<String, Object>> loadDatasetCatalog(Path cfgPath) {
    def parsed = new ConfigSlurper().parse(cfgPath.toUri().toURL())
    def profileSection = (parsed?.profiles instanceof Map) ? parsed.profiles : [:]

    Map<String, Map<String, Object>> catalog = new LinkedHashMap<>()
    profileSection.each { name, body ->
        Map<String, Object> paramsMap = [:]
        def paramsSection = body?.params
        if (paramsSection instanceof Map) {
            paramsSection.each { key, value ->
                paramsMap[key.toString()] = value
            }
        }
        catalog[name.toString()] = paramsMap
    }
    return catalog
}

List<String> resolveRequestedDatasets(Map<String, Map<String, Object>> catalog) {
    def requested = new LinkedHashSet<String>()

    normalizeList(params.datasets).each { requested << it }

    def singleDataset = params.dataset?.toString()?.trim()
    if (singleDataset && singleDataset != 'default') {
        requested << singleDataset
    }

    def profileTokens = (workflow.profile ?: '')
        .tokenize(',')
        .collect { it.trim() }
        .findAll { it }

    profileTokens.each { token ->
        if (catalog.containsKey(token)) {
            requested << token
        }
    }

    if (requested.isEmpty()) {
        return catalog.keySet().toList()
    }

    return requested.toList()
}

String selectCondaEnvironment() {
    def cliSupplied = params.conda_env?.toString()?.trim()
    if (cliSupplied) {
        return cliSupplied
    }

    def configDefault = params.default_conda_env?.toString()?.trim()
    if (configDefault) {
        return configDefault
    }

    return null
}

final String RESOLVED_CONDA_ENV = selectCondaEnvironment()
if (RESOLVED_CONDA_ENV && !params.conda_env) {
    params.conda_env = RESOLVED_CONDA_ENV
}

String deriveCondaEnvLabel(String env) {
    def trimmed = env?.toString()?.trim()
    if (!trimmed) {
        return 'conda_env'
    }
    def tokens = trimmed.tokenize('/')
    def last = tokens ? tokens.last() : trimmed
    def sanitized = (last ?: trimmed).replaceAll(/[^A-Za-z0-9._-]/, '_')
    sanitized ? sanitized : 'conda_env'
}

final String RESOLVED_CONDA_ENV_LABEL = deriveCondaEnvLabel(RESOLVED_CONDA_ENV)

process MaterializeDatasetSpec {
    tag { datasetName }
    conda RESOLVED_CONDA_ENV
    publishDir { "results/datasets/${datasetName}" }, mode: 'copy'
    storeDir { "cache/datasets/${datasetName}" }

    input:
    tuple val(datasetName), val(datasetSpec)

    output:
    tuple val(datasetName), val(datasetSpec), path('dataset_manifest.json')

    script:
    def manifestJson = JsonOutput.prettyPrint(JsonOutput.toJson([
        dataset: datasetName,
        parameters: datasetSpec
    ]))
    """
    |cat <<'EOF' > dataset_manifest.json
    |${manifestJson}
    |EOF
    |
    |python - <<'PY'
    |import json
    |import os
    |import sys
    |from pathlib import Path
    |
    |manifest_path = Path("dataset_manifest.json")
    |data = json.loads(manifest_path.read_text())
    |data["runtime"] = {
    |    "python_executable": sys.executable,
    |    "working_directory": os.getcwd(),
    |    "conda_env": ${RESOLVED_CONDA_ENV.inspect()},
    |    "conda_env_label": ${RESOLVED_CONDA_ENV_LABEL.inspect()}
    |}
    |manifest_path.write_text(json.dumps(data, indent=2))
    |print(f"[flair-eval] Prepared manifest for {data['dataset']}")
    |PY
    """.stripMargin().trim()
}

process RunFlairAlign {
    tag { "${datasetName}::${commandIdx}" }
    conda RESOLVED_CONDA_ENV
    publishDir { "results/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${commandIdx}" }, mode: 'copy'
    storeDir { "cache/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${commandIdx}" }

    input:
    tuple val(datasetName), val(datasetSpec), val(commandIdx), val(commandTemplate), val(renderedCommand)

    output:
    tuple val(datasetName), val(datasetSpec), val(commandIdx), val(commandTemplate), val(renderedCommand),
        path('command_stdout.txt'), path('command_stderr.txt'), path('command_metadata.json'),
        path('*.bam', optional: true), path('*.bed', optional: true), path('*.bai', optional: true)

    script:
    def metadataJson = JsonOutput.prettyPrint(JsonOutput.toJson([
        dataset: datasetName,
        command_index: commandIdx,
        command_template: commandTemplate,
        command_rendered: renderedCommand,
        parameters: datasetSpec,
        runtime: [
            conda_env: RESOLVED_CONDA_ENV,
            conda_env_label: RESOLVED_CONDA_ENV_LABEL
        ]
    ]))
    """
    |set -euo pipefail
    |
    |cat <<'CMD' > command_to_run.sh
    |#!/usr/bin/env bash
    |set -euo pipefail
    |${renderedCommand}
    |CMD
    |
    |chmod +x command_to_run.sh
    |
    |./command_to_run.sh > command_stdout.txt 2> command_stderr.txt
    |
    |cat <<'EOF' > command_metadata.json
    |${metadataJson}
    |EOF
    |
    |python - <<'PY'
    |import glob
    |import json
    |from pathlib import Path
    |
    |meta_path = Path("command_metadata.json")
    |data = json.loads(meta_path.read_text())
    |bam_files = sorted(Path(p).name for p in glob.glob("*.bam"))
    |bed_files = sorted(Path(p).name for p in glob.glob("*.bed"))
    |data["generated_files"] = {
    |    "bam": bam_files,
    |    "bed": bed_files,
    |}
    |meta_path.write_text(json.dumps(data, indent=2))
    |PY
    """.stripMargin().trim()
}

process RunAlignQC {
    tag { "${datasetName}::${commandIdx}" }
    conda RESOLVED_CONDA_ENV
    publishDir { "results/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${commandIdx}" }, mode: 'copy'
    storeDir { "cache/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${commandIdx}/qc" }

    input:
    tuple val(datasetName), val(datasetSpec), val(commandIdx),
        path(stdoutFile), path(stderrFile), path(metadataFile)

    output:
    tuple val(datasetName), val(datasetSpec), val(commandIdx), path('align_qc.tsv')

    script:
    def datasetEsc = escapeForSingleQuotes(datasetName)
    def stdoutName = stdoutFile.getFileName().toString()
    def stderrName = stderrFile.getFileName().toString()
    def metadataName = metadataFile.getFileName().toString()

    """
    |set -euo pipefail
    |
    |python ${projectDir}/bin/align_qc.py \\
    |    --dataset '${datasetEsc}' \\
    |    --command-idx ${commandIdx} \\
    |    --stdout '${stdoutName}' \\
    |    --stderr '${stderrName}' \\
    |    --metadata '${metadataName}' \\
    |    --output 'align_qc.tsv'
    """.stripMargin().trim()
}

process RunFlairRegionalize {
    tag { "${datasetName}::align${alignIdx}::region${regionIdx}" }
    conda RESOLVED_CONDA_ENV
    publishDir { "results/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}/flair_regionalize${regionIdx}" }, mode: 'copy'
    storeDir { "cache/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}/flair_regionalize${regionIdx}" }

    input:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx), val(regionSpec),
        path(alignBam), path(alignBed), path(alignMetadata)

    output:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx), val(regionSpec),
        path('regionalize_stdout.txt'), path('regionalize_stderr.txt'), path('regionalize_metadata.json'),
        path('region_details.tsv', optional: true),
        path('*.bam', optional: true), path('*.bam.bai', optional: true),
        path('*.bed', optional: true), path('*.gtf', optional: true),
        path('*.fa', optional: true), path('*.SJ.out.tab', optional: true)

    script:
    def gtfPath = datasetSpec?.gtf?.toString()
    if (!gtfPath) {
        throw new IllegalArgumentException("Dataset '${datasetName}' is missing 'gtf' in its specification, required for regionalize.")
    }

    def regionMode = regionSpec.mode?.toString()
    def commandText = regionSpec.command_text?.toString() ?: ''
    def datasetEsc = escapeForSingleQuotes(datasetName)
    def commandTextEsc = escapeForSingleQuotes(commandText)
    def gtfEsc = escapeForSingleQuotes(gtfPath)
    def alignBamName = alignBam.getFileName().toString()
    def alignBedName = alignBed.getFileName().toString()
    def alignMetadataName = alignMetadata.getFileName().toString()

    def regionArg = (regionMode == 'region')
        ? "--region '${escapeForSingleQuotes(regionSpec.region_string?.toString())}'"
        : "--all"

    Map<String, String> optionalMapping = [
        genome         : '--genome',
        junctions      : '--junctions',
        cage           : '--cage',
        quantseq       : '--quantseq',
        target_regions : '--target-regions'
    ]

    List<String> optionalArgs = []
    optionalMapping.each { key, flag ->
        def value = datasetSpec[key]
        if (value) {
            optionalArgs << "${flag} '${escapeForSingleQuotes(value.toString())}'"
        }
    }

    def commandParts = [
        "python ${projectDir}/bin/regionalize.py",
        "--dataset '${datasetEsc}'",
        "--align-index ${alignIdx}",
        "--command-index ${regionIdx}",
        "--command-text='${commandTextEsc}'",
        "--conda-env-label '${escapeForSingleQuotes(RESOLVED_CONDA_ENV_LABEL)}'",
        "--output-dir '.'",
        "--align-bam '${alignBamName}'",
        "--align-bed '${alignBedName}'",
        "--align-metadata '${alignMetadataName}'",
        "--gtf '${gtfEsc}'",
        regionArg
    ] + optionalArgs

    def commandString = commandParts.join(" \\\n    ")

    """
    |set -euo pipefail
    |
    |${commandString} \\
    |    > regionalize_stdout.txt 2> regionalize_stderr.txt
    """.stripMargin().trim()
}

process RunRegionalizeQC {
    tag { "${datasetName}::align${alignIdx}::region${regionIdx}" }
    conda RESOLVED_CONDA_ENV
    publishDir { "results/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}/flair_regionalize${regionIdx}" }, mode: 'copy'
    storeDir { "cache/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}/flair_regionalize${regionIdx}/qc" }

    input:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx), val(regionSpec),
        path(stdoutFile), path(stderrFile), path(metadataFile), val(regionDetailsRef)

    output:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx), val(regionSpec),
        path('regionalize_qc.tsv')

    script:
    def datasetEsc = escapeForSingleQuotes(datasetName)
    def modeValue = regionSpec?.mode?.toString() ?: 'unknown'
    def modeEsc = escapeForSingleQuotes(modeValue)
    def commandText = regionSpec?.command_text?.toString() ?: ''
    def commandEsc = escapeForSingleQuotes(commandText)
    def stdoutName = stdoutFile.getFileName().toString()
    def stderrName = stderrFile.getFileName().toString()
    def metadataName = metadataFile.getFileName().toString()
    def regionDetailsPath = null
    if (regionDetailsRef) {
        if (regionDetailsRef instanceof Collection && !regionDetailsRef.isEmpty()) {
            def first = regionDetailsRef[0]
            regionDetailsPath = first ? first.toString() : null
        } else if (!(regionDetailsRef instanceof Collection)) {
            regionDetailsPath = regionDetailsRef.toString()
        }
    }

    List<String> argsList = [
        "python ${projectDir}/bin/regionalize_qc.py",
        "--dataset '${datasetEsc}'",
        "--align-idx ${alignIdx}",
        "--region-idx ${regionIdx}",
        "--mode '${modeEsc}'",
        "--command-text '${commandEsc}'",
        "--stdout '${stdoutName}'",
        "--stderr '${stderrName}'",
        "--metadata '${metadataName}'"
    ]
    if (regionDetailsPath) {
        argsList << "--region-details '${escapeForSingleQuotes(regionDetailsPath)}'"
    }
    argsList << "--output 'regionalize_qc.tsv'"
    def cmd = argsList.join(" \\\n    ")

    """
    |set -euo pipefail
    |
    |${cmd}
    """.stripMargin().trim()
}

workflow flair_eval {
    take:
        // none

    main:
        def condaEnv = RESOLVED_CONDA_ENV
        if (!condaEnv) {
            throw new IllegalArgumentException("Specify the conda environment path with --conda_env or set params.default_conda_env in a config file")
        }

        def cfgPath = locateDatasetsConfig()
        def catalog = loadDatasetCatalog(cfgPath)
        if (catalog.isEmpty()) {
            throw new IllegalStateException("No dataset profiles discovered in ${cfgPath}")
        }

        def requested = resolveRequestedDatasets(catalog)
        def missing = requested.findAll { !catalog.containsKey(it) }
        if (missing) {
            throw new IllegalArgumentException("Unknown dataset profile(s): ${missing.join(', ')}")
        }

        def datasetTuples = requested.collect { name ->
            def spec = new LinkedHashMap<String, Object>(catalog[name])
            def summary = spec.collect { k, v -> "${k}=${v}" }.join(', ')
            log.info("[flair-eval] Dataset ${name}: ${summary}")
            tuple(name, spec)
        }

        def dataset_for_manifest = Channel.from(datasetTuples)
        def dataset_for_commands = Channel.from(datasetTuples)

        dataset_manifests = MaterializeDatasetSpec(dataset_for_manifest)

        def flairCommands = normalizeCommands(params.flair_commands)
        if (!flairCommands) {
            log.info("No --flair_commands supplied; skipping flair align runs.")
        }

        def commandEntries = []
        flairCommands.eachWithIndex { cmd, idx ->
            commandEntries << [idx + 1, cmd]
        }

        def align_payloads = Channel.empty()
        def align_qc_results_channel = Channel.empty()
        if (commandEntries) {
            def align_inputs = dataset_for_commands.flatMap { datasetName, spec ->
                commandEntries.collect { info ->
                    def commandIdx = info[0] as int
                    def template = info[1] as String
                    def rendered = renderFlairCommand(template, spec, datasetName)
                    tuple(datasetName, spec, commandIdx, template, rendered)
                }
            }
            def raw_align_results = RunFlairAlign(align_inputs)
            def align_results_for_payload = raw_align_results.map { it }
            def align_results_for_qc = raw_align_results.map { it }

            align_payloads = align_results_for_payload.map { values ->
                def datasetName = values[0]
                def datasetSpec = values[1]
                def alignIdx = values[2]
                def commandTemplate = values[3]
                def renderedCommand = values[4]
                def stdoutFile = values[5]
                def stderrFile = values[6]
                def metadataFile = values[7]
                def bamPaths = values.size() > 8 ? values[8] : []
                def bedPaths = values.size() > 9 ? values[9] : []
                def bamFiles = flattenPathList(bamPaths)
                    .findAll { it.getFileName().toString().toLowerCase().endsWith('.bam') && !it.getFileName().toString().toLowerCase().endsWith('.bam.bai') }
                def bedFiles = flattenPathList(bedPaths)
                    .findAll { it.getFileName().toString().toLowerCase().endsWith('.bed') }
                if (bamFiles.isEmpty()) {
                    throw new IllegalStateException("Align command ${alignIdx} for dataset ${datasetName} did not produce a BAM file.")
                }
                if (bedFiles.isEmpty()) {
                    throw new IllegalStateException("Align command ${alignIdx} for dataset ${datasetName} did not produce a BED file.")
                }
                [
                    datasetName    : datasetName,
                    datasetSpec    : datasetSpec,
                    alignIndex     : alignIdx,
                    commandTemplate: commandTemplate,
                    renderedCommand: renderedCommand,
                    stdoutPath     : stdoutFile,
                    stderrPath     : stderrFile,
                    metadataPath   : metadataFile,
                    bamFiles       : bamFiles,
                    bedFiles       : bedFiles
                ]
            }

            def align_qc_inputs = align_results_for_qc.map { values ->
                def datasetName = values[0]
                def datasetSpec = values[1]
                def alignIdx = values[2]
                def stdoutFile = values[5]
                def stderrFile = values[6]
                def metadataFile = values[7]
                tuple(datasetName, datasetSpec, alignIdx, stdoutFile, stderrFile, metadataFile)
            }
            align_qc_results_channel = RunAlignQC(align_qc_inputs)
        }

        align_results = align_payloads
        align_qc_results = align_qc_results_channel

        def regionalizeCommandsRaw = normalizeCommands(params.regionalize_commands)
        if (!regionalizeCommandsRaw) {
            regionalizeCommandsRaw = ['--all']
        }

        List<Map<String, Object>> regionalizeSpecs = []
        regionalizeCommandsRaw.eachWithIndex { cmd, idx ->
            def spec = parseRegionalizeCommand(cmd)
            spec['command_index'] = idx + 1
            regionalizeSpecs << spec
        }

        def regionalize_results_channel = Channel.empty()
        def regionalize_qc_results_channel = Channel.empty()
        if (commandEntries && !regionalizeSpecs.isEmpty()) {
            def requests = align_payloads.flatMap { payload ->
                regionalizeSpecs.collect { spec ->
                    def copy = new LinkedHashMap<String, Object>(spec)
                    def mode = (copy.mode ?: 'region').toString().toLowerCase()
                    copy['_align_payload'] = payload
                    [
                        datasetName : payload.datasetName,
                        datasetSpec : payload.datasetSpec,
                        alignIndex  : payload.alignIndex,
                        regionIndex : copy.command_index,
                        regionSpec  : copy,
                        mode        : mode,
                        alignPayload: payload
                    ]
                }
            }

            def run_requests = requests.filter { it.mode != 'all' }
            def bypass_requests = requests.filter { it.mode == 'all' }

            def run_inputs = run_requests.map { req ->
                def alignPayload = req.alignPayload
                def bam = alignPayload.bamFiles[0]
                def bed = alignPayload.bedFiles[0]
                tuple(
                    req.datasetName,
                    req.datasetSpec,
                    req.alignIndex,
                    req.regionIndex,
                    req.regionSpec,
                    bam,
                    bed,
                    alignPayload.metadataPath
                )
            }
            def raw_run_results = RunFlairRegionalize(run_inputs)
            def run_results_for_map = raw_run_results.map { it }
            def run_results_for_qc = raw_run_results.map { it }

            def run_results = run_results_for_map.map { values ->
                def regionSpec = new LinkedHashMap<String, Object>(values[4])
                def alignPayload = regionSpec.remove('_align_payload')
                [
                    datasetName  : values[0],
                    datasetSpec  : values[1],
                    alignIndex   : values[2],
                    regionIndex  : values[3],
                    regionSpec   : regionSpec,
                    mode         : 'region',
                    stdoutPath   : values[5],
                    stderrPath   : values[6],
                    metadataPath : values[7],
                    regionDetails: values.size() > 8 ? values[8] : null,
                    bamFiles     : flattenPathList(values.size() > 9 ? values[9] : []),
                    baiFiles     : flattenPathList(values.size() > 10 ? values[10] : []),
                    bedFiles     : flattenPathList(values.size() > 11 ? values[11] : []),
                    gtfFiles     : flattenPathList(values.size() > 12 ? values[12] : []),
                    fastaFiles   : flattenPathList(values.size() > 13 ? values[13] : []),
                    sjFiles      : flattenPathList(values.size() > 14 ? values[14] : []),
                    alignPayload : alignPayload
                ]
            }

            def regionalize_qc_inputs = run_results_for_qc.map { values ->
                def datasetName = values[0]
                def datasetSpec = values[1]
                def alignIdx = values[2]
                def regionIdx = values[3]
                def regionSpec = new LinkedHashMap<String, Object>(values[4])
                regionSpec.remove('_align_payload')
                def stdoutFile = values[5]
                def stderrFile = values[6]
                def metadataFile = values[7]
                def regionDetails = values.size() > 8 ? values[8] : null
                tuple(datasetName, datasetSpec, alignIdx, regionIdx, regionSpec, stdoutFile, stderrFile, metadataFile, regionDetails)
            }
            regionalize_qc_results_channel = RunRegionalizeQC(regionalize_qc_inputs)

            def bypass_results = bypass_requests.map { req ->
                def regionSpec = new LinkedHashMap<String, Object>(req.regionSpec)
                def alignPayload = regionSpec.remove('_align_payload') ?: req.alignPayload
                [
                    datasetName  : req.datasetName,
                    datasetSpec  : req.datasetSpec,
                    alignIndex   : req.alignIndex,
                    regionIndex  : req.regionIndex,
                    regionSpec   : regionSpec,
                    mode         : 'all',
                    stdoutPath   : alignPayload.stdoutPath,
                    stderrPath   : alignPayload.stderrPath,
                    metadataPath : alignPayload.metadataPath,
                    regionDetails: null,
                    bamFiles     : alignPayload.bamFiles,
                    baiFiles     : [],
                    bedFiles     : alignPayload.bedFiles,
                    gtfFiles     : [],
                    fastaFiles   : [],
                    sjFiles      : [],
                    alignPayload : alignPayload
                ]
            }

            regionalize_results_channel = run_results.mix(bypass_results)
        }
        regionalize_results = regionalize_results_channel
        regionalize_qc_results = regionalize_qc_results_channel

    emit:
        dataset_manifests
        align_results
        align_qc_results
        regionalize_results
        regionalize_qc_results
}

workflow {
    flair_eval()
}
