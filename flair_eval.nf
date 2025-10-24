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
    publishDir "results/datasets/${datasetName}", mode: 'copy'
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
    publishDir "results/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${commandIdx}", mode: 'copy'
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
    publishDir "results/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${commandIdx}", mode: 'copy'
    storeDir { "cache/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${commandIdx}/qc" }

    input:
    tuple val(datasetName), val(datasetSpec), val(commandIdx),
        path(stdoutFile), path(stderrFile), path(metadataFile), path(alignQcScript)

    output:
    tuple val(datasetName), val(datasetSpec), val(commandIdx), path('align_qc.tsv')

    script:
    def datasetEsc = escapeForSingleQuotes(datasetName)
    def stdoutName = stdoutFile.getFileName().toString()
    def stderrName = stderrFile.getFileName().toString()
    def metadataName = metadataFile.getFileName().toString()
    def alignQcPath = alignQcScript.toString()

    """
    |set -euo pipefail
    |
    |python ${alignQcPath} \\
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
    publishDir "results/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}/flair_regionalize${regionIdx}", mode: 'copy'
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
    publishDir "results/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}/flair_regionalize${regionIdx}", mode: 'copy'
    storeDir { "cache/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}/flair_regionalize${regionIdx}/qc" }

    input:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx), val(regionSpec),
        path(stdoutFile), path(stderrFile), path(metadataFile), val(regionDetailsRef), path(regionalizeQcScript)

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
    def regionalizeQcPath = regionalizeQcScript.toString()
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
    def options = argsList.join(" \\\n    ")

    """
    |set -euo pipefail
    |
    |python ${regionalizeQcPath} \\
    |    ${options}
    """.stripMargin().trim()
}

process RunFlairCorrect {
    tag {
        def parts = ["${datasetName}", "align${alignIdx}"]
        parts << ((mode == 'region') ? "region${regionIdx}" : mode)
        parts << "correct${commandIdx}"
        parts.join("::")
    }
    conda RESOLVED_CONDA_ENV
    publishDir "results/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}${(mode == 'region') ? "/flair_regionalize${regionIdx}" : ""}/correct/flair_correct${commandIdx}", mode: 'copy'
    storeDir "cache/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}${(mode == 'region') ? "/flair_regionalize${regionIdx}" : ""}/correct/flair_correct${commandIdx}"

    input:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx), val(regionSpec),
        val(mode), val(regionTag), val(commandIdx), val(commandTemplate), val(renderedCommand),
        path(correctInputBed),
        val(correctInputBam),
        val(alignMetadata),
        val(alignPayload)

    output:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx), val(regionSpec),
        val(mode), val(regionTag), val(commandIdx), val(commandTemplate), val(renderedCommand),
        path('correct_stdout.txt'), path('correct_stderr.txt'), path('correct_metadata.json'),
        path('*_all_corrected.bed', optional: true),
        path('*_all_inconsistent.bed', optional: true),
        path('*_cannot_verify.bed', optional: true),
        val(alignPayload)

    script:
    def bedPath = correctInputBed as Path
    def bamPath = correctInputBam ? (correctInputBam as Path) : null
    def alignMetaPath = alignMetadata ? (alignMetadata as Path) : null

    def bedName = bedPath.getFileName().toString()
    def bedNameEsc = escapeForSingleQuotes(bedName)
    def bamName = bamPath ? bamPath.getFileName().toString() : null
    def bamSourceEsc = bamPath ? escapeForSingleQuotes(bamPath.toAbsolutePath().toString()) : null
    def bamNameEsc = bamName ? escapeForSingleQuotes(bamName) : null
    def alignMetaName = alignMetaPath ? alignMetaPath.getFileName().toString() : null
    def alignMetaSourceEsc = alignMetaPath ? escapeForSingleQuotes(alignMetaPath.toAbsolutePath().toString()) : null
    def alignMetaNameEsc = alignMetaName ? escapeForSingleQuotes(alignMetaName) : null
    def bedStem = bedName.toLowerCase().endsWith('.bed') ? bedName.substring(0, bedName.length() - 4) : bedName
    def cannotVerifyName = "${bedStem}_cannot_verify.bed"
    def cannotVerifyNameEsc = escapeForSingleQuotes(cannotVerifyName)

    def metadataMap = [
        dataset          : datasetName,
        align_index      : alignIdx,
        region_index     : regionIdx,
        region_mode      : mode,
        region_tag       : regionTag,
        command_index    : commandIdx,
        command_template : commandTemplate,
        command_rendered : renderedCommand,
        correct_input_bed: bedName,
        correct_input_bam: bamName,
        align_metadata   : alignMetaName,
        region_spec      : regionSpec,
        runtime          : [
            conda_env      : RESOLVED_CONDA_ENV,
            conda_env_label: RESOLVED_CONDA_ENV_LABEL
        ]
    ]
    def metadataJson = JsonOutput.prettyPrint(JsonOutput.toJson(metadataMap))

    """
    |set -euo pipefail
    |
    |if [ ${bamSourceEsc ? 1 : 0} -eq 1 ]; then
    |    ln -sf '${bamSourceEsc}' '${bamNameEsc}'
    |fi
    |if [ ${alignMetaSourceEsc ? 1 : 0} -eq 1 ]; then
    |    ln -sf '${alignMetaSourceEsc}' '${alignMetaNameEsc}'
    |fi
    |
    |cat <<'CMD' > command_to_run.sh
    |#!/usr/bin/env bash
    |set -euo pipefail
    |${renderedCommand}
    |CMD
    |
    |chmod +x command_to_run.sh
    |
    |./command_to_run.sh > correct_stdout.txt 2> correct_stderr.txt
    |
    |if [ ! -f '${cannotVerifyNameEsc}' ]; then
    |    : > '${cannotVerifyNameEsc}'
    |fi
    |
    |cat <<'EOF' > correct_metadata.json
    |${metadataJson}
    |EOF
    |
    |python - <<'PY'
    |from pathlib import Path
    |import json
    |
    |meta = Path('correct_metadata.json')
    |data = json.loads(meta.read_text())
    |stage = Path('.')
    |expected = {
    |    "corrected": "_all_corrected.bed",
    |    "inconsistent": "_all_inconsistent.bed",
    |    "cannot_verify": "_cannot_verify.bed",
    |}
    |observed = {}
    |for key, suffix in expected.items():
    |    files = sorted(p.name for p in stage.glob(f"*{suffix}"))
    |    observed[key] = files
    |    if key in ("corrected", "inconsistent") and not files:
    |        raise SystemExit(f"Expected output matching '*{suffix}' not produced.")
    |data["observed_outputs"] = observed
    |meta.write_text(json.dumps(data, indent=2))
    |PY
    """.stripMargin().trim()
}

process RunCorrectQC {
    tag {
        def parts = ["${datasetName}", "align${alignIdx}"]
        parts << ((mode == 'region') ? "region${regionIdx}" : mode)
        parts << "correct${commandIdx}"
        parts.join("::")
    }
    conda RESOLVED_CONDA_ENV
    publishDir "results/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}${(mode == 'region') ? "/flair_regionalize${regionIdx}" : ""}/correct/flair_correct${commandIdx}", mode: 'copy'
    storeDir "cache/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}${(mode == 'region') ? "/flair_regionalize${regionIdx}" : ""}/correct/flair_correct${commandIdx}/qc"

    input:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx), val(regionSpec),
        val(mode), val(regionTag), val(commandIdx),
        path(stdoutFile), path(stderrFile), path(metadataFile),
        val(correctedBed), val(inconsistentBed), val(cannotVerifyBed),
        val(correctQcScript)

    output:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx),
        val(mode), val(regionTag), val(commandIdx), path('correct_qc.tsv')

    script:
    def stdoutPath = stdoutFile as Path
    def stderrPath = stderrFile as Path
    def metadataPath = metadataFile as Path
    def correctedPath = correctedBed ? (correctedBed as Path) : null
    def inconsistentPath = inconsistentBed ? (inconsistentBed as Path) : null
    def cannotVerifyPath = cannotVerifyBed ? (cannotVerifyBed as Path) : null

    def datasetEsc = escapeForSingleQuotes(datasetName)
    def modeEsc = escapeForSingleQuotes(mode?.toString() ?: '')
    def regionTagEsc = escapeForSingleQuotes(regionTag?.toString() ?: '')
    def stdoutNameEsc = escapeForSingleQuotes(stdoutPath.getFileName().toString())
    def stderrNameEsc = escapeForSingleQuotes(stderrPath.getFileName().toString())
    def metadataNameEsc = escapeForSingleQuotes(metadataPath.getFileName().toString())
    def correctedEsc = correctedPath ? escapeForSingleQuotes(correctedPath.toAbsolutePath().toString()) : null
    def inconsistentEsc = inconsistentPath ? escapeForSingleQuotes(inconsistentPath.toAbsolutePath().toString()) : null
    def cannotVerifyEsc = cannotVerifyPath ? escapeForSingleQuotes(cannotVerifyPath.toAbsolutePath().toString()) : null
    def correctQcPathEsc = escapeForSingleQuotes((correctQcScript as Path).toAbsolutePath().toString())

    List<String> argsList = [
        "--dataset '${datasetEsc}'",
        "--align-idx ${alignIdx}",
        "--region-idx ${regionIdx != null ? regionIdx : 0}",
        "--mode '${modeEsc}'",
        "--region-tag '${regionTagEsc}'",
        "--stdout '${stdoutNameEsc}'",
        "--stderr '${stderrNameEsc}'",
        "--metadata '${metadataNameEsc}'",
        "--output 'correct_qc.tsv'"
    ]
    if (correctedEsc) {
        argsList << "--corrected '${correctedEsc}'"
    }
    if (inconsistentEsc) {
        argsList << "--inconsistent '${inconsistentEsc}'"
    }
    if (cannotVerifyEsc) {
        argsList << "--cannot-verify '${cannotVerifyEsc}'"
    }
    def options = argsList.join(" \\\n    ")

    """
    |set -euo pipefail
    |
    |python ${correctQcPathEsc} \\
    |    ${options}
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

        def flairCommands = normalizeCommands(params.flair_align_commands)
        if (!flairCommands) {
            log.info("No --flair_align_commands supplied; skipping flair align runs.")
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
                def alignQcScript = file("${projectDir}/bin/align_qc.py")
                tuple(datasetName, datasetSpec, alignIdx, stdoutFile, stderrFile, metadataFile, alignQcScript)
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
                def regionalizeQcScript = file("${projectDir}/bin/regionalize_qc.py")
                tuple(datasetName, datasetSpec, alignIdx, regionIdx, regionSpec, stdoutFile, stderrFile, metadataFile, regionDetails, regionalizeQcScript)
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

        def correctCommandsRaw = normalizeCommands(params.flair_correct_commands)
        if (!correctCommandsRaw) {
            log.info("No --flair_correct_commands supplied; skipping flair correct runs.")
        }

        def correctCommandEntries = []
        correctCommandsRaw.eachWithIndex { cmd, idx ->
            correctCommandEntries << [idx + 1, cmd]
        }

        def correct_results_channel = Channel.empty()
        def correct_qc_results_channel = Channel.empty()
        if (!correctCommandEntries.isEmpty()) {
            def correct_inputs = regionalize_results.flatMap { payload ->
                correctCommandEntries.collect { info ->
                    def commandIdx = info[0] as int
                    def template = info[1] as String

                    def regionSpec = new LinkedHashMap<String, Object>(payload.regionSpec ?: [:])
                    def alignPayload = payload.alignPayload

                    List<Path> bedCandidates = flattenPathList(payload.bedFiles ?: [])
                    if (!bedCandidates && alignPayload) {
                        bedCandidates = flattenPathList(alignPayload.bedFiles ?: [])
                    }
                    if (!bedCandidates) {
                        throw new IllegalStateException("Correct stage requires a BED file for dataset ${payload.datasetName}, align ${payload.alignIndex}, region ${payload.regionIndex}.")
                    }
                    Path bedFile = null
                    if (payload.mode == 'region') {
                        bedFile = bedCandidates.find { it.getFileName().toString() ==~ /(?i)[A-Za-z0-9._-]+_\d+_\d+\.bed$/ }
                    }
                    bedFile = bedFile ?: bedCandidates.find { it.getFileName().toString().toLowerCase().endsWith('.bed') }
                    if (!bedFile) {
                        throw new IllegalStateException("Correct stage could not select a BED file for dataset ${payload.datasetName}, align ${payload.alignIndex}, region ${payload.regionIndex}.")
                    }
                    def regionTagRaw = bedFile.getFileName().toString()
                    def regionTag = regionTagRaw.toLowerCase().endsWith('.bed')
                        ? regionTagRaw.substring(0, regionTagRaw.length() - 4)
                        : regionTagRaw

                    List<Path> bamCandidates = flattenPathList(payload.bamFiles ?: [])
                    if (!bamCandidates && alignPayload) {
                        bamCandidates = flattenPathList(alignPayload.bamFiles ?: [])
                    }
                    Path bamFile = null
                    if (payload.mode == 'region') {
                        bamFile = bamCandidates.find { it.getFileName().toString() ==~ /(?i)[A-Za-z0-9._-]+_\d+_\d+\.bam$/ }
                    }
                    bamFile = bamFile ?: (bamCandidates ? bamCandidates[0] : null)

                    def substitution = new LinkedHashMap<String, Object>(payload.datasetSpec ?: [:])
                    substitution['correct_bed'] = bedFile.getFileName().toString()
                    substitution['correct_bam'] = bamFile ? bamFile.getFileName().toString() : null
                    substitution['correct_output_prefix'] = regionTag
                    substitution['correct_mode'] = payload.mode
                    substitution['correct_region'] = regionTag
                    substitution['align_index'] = payload.alignIndex
                    substitution['region_index'] = payload.regionIndex
                    def rendered = renderFlairCommand(template, substitution, payload.datasetName)

                    def alignMetadata = alignPayload?.metadataPath
                    tuple(
                        payload.datasetName,
                        payload.datasetSpec,
                        payload.alignIndex,
                        payload.regionIndex,
                        regionSpec,
                        payload.mode,
                        regionTag,
                        commandIdx,
                        template,
                        rendered,
                        bedFile,
                        bamFile,
                        alignMetadata,
                        alignPayload
                    )
                }
            }
            def raw_correct_results = RunFlairCorrect(correct_inputs)
            def correct_results_for_map = raw_correct_results.map { it }
            def correct_results_for_qc = raw_correct_results.map { it }

            correct_results_channel = correct_results_for_map.map { values ->
                def datasetName = values[0]
                def datasetSpec = values[1]
                def alignIdx = values[2]
                def regionIdx = values[3]
                def regionSpec = values[4]
                def mode = values[5]
                def regionTag = values[6]
                def commandIdx = values[7]
                def commandTemplate = values[8]
                def renderedCommand = values[9]
                def stdoutFile = values[10]
                def stderrFile = values[11]
                def metadataFile = values[12]
                def correctedPaths = flattenPathList(values.size() > 13 ? values[13] : [])
                def inconsistentPaths = flattenPathList(values.size() > 14 ? values[14] : [])
                def cannotVerifyPaths = flattenPathList(values.size() > 15 ? values[15] : [])
                def alignPayload = values[16]
                [
                    datasetName      : datasetName,
                    datasetSpec      : datasetSpec,
                    alignIndex       : alignIdx,
                    regionIndex      : regionIdx,
                    regionSpec       : regionSpec,
                    mode             : mode,
                    regionTag        : regionTag,
                    commandIndex     : commandIdx,
                    commandTemplate  : commandTemplate,
                    renderedCommand  : renderedCommand,
                    stdoutPath       : stdoutFile,
                    stderrPath       : stderrFile,
                    metadataPath     : metadataFile,
                    correctedBeds    : correctedPaths,
                    inconsistentBeds : inconsistentPaths,
                    cannotVerifyBeds : cannotVerifyPaths,
                    alignPayload     : alignPayload
                ]
            }

            def correct_qc_inputs = correct_results_for_qc.map { values ->
                def datasetName = values[0]
                def datasetSpec = values[1]
                def alignIdx = values[2]
                def regionIdx = values[3]
                def regionSpec = values[4]
                def mode = values[5]
                def regionTag = values[6]
                def commandIdx = values[7]
                def stdoutFile = values[10]
                def stderrFile = values[11]
                def metadataFile = values[12]
                def extractSinglePath = { obj ->
                    if (!obj) {
                        return null
                    }
                    if (obj instanceof Collection) {
                        def first = obj.find { it }
                        return first ?: null
                    }
                    obj
                }
                def corrected = extractSinglePath(values.size() > 13 ? values[13] : null)
                def inconsistent = extractSinglePath(values.size() > 14 ? values[14] : null)
                def cannotVerify = extractSinglePath(values.size() > 15 ? values[15] : null)
                def correctQcScript = file("${projectDir}/bin/correct_qc.py")
                tuple(datasetName, datasetSpec, alignIdx, regionIdx, regionSpec, mode, regionTag, commandIdx, stdoutFile, stderrFile, metadataFile, corrected, inconsistent, cannotVerify, correctQcScript)
            }
            correct_qc_results_channel = RunCorrectQC(correct_qc_inputs)
        }
        correct_results = correct_results_channel
        correct_qc_results = correct_qc_results_channel

    emit:
        dataset_manifests
        align_results
        align_qc_results
        regionalize_results
        regionalize_qc_results
        correct_results
        correct_qc_results
}

workflow {
    flair_eval()
}
