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
        val(regionalizePayload),
        val(alignPayload)

    output:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx), val(regionSpec),
        val(mode), val(regionTag), val(commandIdx), val(commandTemplate), val(renderedCommand),
        path('correct_stdout.txt'), path('correct_stderr.txt'), path('correct_metadata.json'),
        path('*_all_corrected.bed', optional: true),
        path('*_all_inconsistent.bed', optional: true),
        path('*_cannot_verify.bed', optional: true),
        val(regionalizePayload),
        val(alignPayload)

    script:
    def bedPath = correctInputBed as Path
    def bamPath = correctInputBam ? (correctInputBam as Path) : null
    def alignMetaPath = alignMetadata ? (alignMetadata as Path) : null

    def bedName = bedPath.getFileName().toString()
    def bedNameEsc = escapeForSingleQuotes(bedName)
    def bedEmpty = !bedPath.toFile().exists() || bedPath.toFile().length() == 0L
    def bedEmptyFlag = bedEmpty ? 1 : 0
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

process RunFlairTranscriptome {
    tag {
        def parts = ["${datasetName}", "align${alignIdx}"]
        parts << ((mode == 'region') ? "region${regionIdx}" : mode)
        parts << "transcriptome${commandIdx}"
        parts.join("::")
    }
    conda RESOLVED_CONDA_ENV
    publishDir "results/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}${(mode == 'region') ? "/flair_regionalize${regionIdx}" : ""}/transcriptome/flair_transcriptome${commandIdx}", mode: 'copy'
    storeDir "cache/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}${(mode == 'region') ? "/flair_regionalize${regionIdx}" : ""}/transcriptome/flair_transcriptome${commandIdx}"

    input:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx), val(regionSpec),
        val(mode), val(regionTag), val(commandIdx), val(commandTemplate), val(renderedCommand),
        path(transcriptomeInputBam),
        val(regionalizeMetadata),
        val(alignMetadata),
        val(regionalizePayload),
        val(alignPayload)

    output:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx), val(regionSpec),
        val(mode), val(regionTag), val(commandIdx), val(commandTemplate), val(renderedCommand),
        path('transcriptome_stdout.txt'), path('transcriptome_stderr.txt'), path('transcriptome_metadata.json'),
        path('flair.isoforms.bed'), path('flair.isoforms.gtf'), path('flair.isoforms.fa'), path('flair.read.map.txt'),
        path('transcriptome_exit_code.txt'), path('transcriptome_error.log', optional: true),
        val(regionalizePayload), val(alignPayload)

    script:
    def bamPath = transcriptomeInputBam as Path
    def bamName = bamPath.getFileName().toString()
    def bamNameEsc = escapeForSingleQuotes(bamName)
    def bamSourceEsc = escapeForSingleQuotes(bamPath.toAbsolutePath().toString())

    def regionMetaPath = regionalizeMetadata ? (regionalizeMetadata as Path) : null
    def regionMetaName = regionMetaPath ? regionMetaPath.getFileName().toString() : null
    def regionMetaSourceEsc = regionMetaPath ? escapeForSingleQuotes(regionMetaPath.toAbsolutePath().toString()) : null
    def regionMetaNameEsc = regionMetaName ? escapeForSingleQuotes(regionMetaName) : null

    def alignMetaPath = alignMetadata ? (alignMetadata as Path) : null
    def alignMetaName = alignMetaPath ? alignMetaPath.getFileName().toString() : null
    def alignMetaSourceEsc = alignMetaPath ? escapeForSingleQuotes(alignMetaPath.toAbsolutePath().toString()) : null
    def alignMetaNameEsc = alignMetaName ? escapeForSingleQuotes(alignMetaName) : null

    def genomePath = datasetSpec?.genome?.toString()
    if (!genomePath) {
        throw new IllegalArgumentException("Dataset '${datasetName}' is missing 'genome' required for transcriptome.")
    }

    def exitCodeName = 'transcriptome_exit_code.txt'
    def exitCodeNameEsc = escapeForSingleQuotes(exitCodeName)
    def errorLogName = 'transcriptome_error.log'
    def errorLogNameEsc = escapeForSingleQuotes(errorLogName)

    def metadataMap = [
        dataset              : datasetName,
        align_index          : alignIdx,
        region_index         : regionIdx,
        region_mode          : mode,
        region_tag           : regionTag,
        command_index        : commandIdx,
        command_template     : commandTemplate,
        command_rendered     : renderedCommand,
        transcriptome_input_bam: bamName,
        genome               : genomePath,
        regionalize_metadata : regionMetaName,
        align_metadata       : alignMetaName,
        region_spec          : regionSpec,
        runtime              : [
            conda_env      : RESOLVED_CONDA_ENV,
            conda_env_label: RESOLVED_CONDA_ENV_LABEL
        ]
    ]
    def metadataJson = JsonOutput.prettyPrint(JsonOutput.toJson(metadataMap))

    """
    |set -euo pipefail
    |
    |ln -sf '${bamSourceEsc}' '${bamNameEsc}'
    |if [ ${regionMetaSourceEsc ? 1 : 0} -eq 1 ]; then
    |    ln -sf '${regionMetaSourceEsc}' '${regionMetaNameEsc}'
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
    |set +e
    |./command_to_run.sh > transcriptome_stdout.txt 2> transcriptome_stderr.txt
    |status=\$?
    |set -euo pipefail
    |
    |echo "\${status}" > ${exitCodeNameEsc}
    |if [ "\${status}" -ne 0 ]; then
    |    printf 'Transcriptome command failed with exit code %s\n' "\${status}" > ${errorLogNameEsc}
    |    touch flair.isoforms.bed flair.isoforms.gtf flair.isoforms.fa flair.read.map.txt
    |fi
    |
    |cat <<'EOF' > transcriptome_metadata.json
    |${metadataJson}
    |EOF
    |
    |python - <<'PY'
    |from pathlib import Path
    |import json
    |
    |meta = Path('transcriptome_metadata.json')
    |data = json.loads(meta.read_text())
    |stage = Path('.')
    |exit_code_path = Path('${exitCodeName}')
    |exit_code = 0
    |if exit_code_path.exists():
    |    try:
    |        exit_code = int(exit_code_path.read_text().strip())
    |    except ValueError:
    |        exit_code = 0
    |data["exit_code"] = exit_code
    |expected = {
    |    "isoforms_bed": "flair.isoforms.bed",
    |    "isoforms_gtf": "flair.isoforms.gtf",
    |    "isoforms_fa": "flair.isoforms.fa",
    |    "read_map": "flair.read.map.txt",
    |}
    |observed = {}
    |missing = []
    |for key, name in expected.items():
    |    path = stage / name
    |    if not path.exists():
    |        observed[key] = None
    |        if exit_code == 0:
    |            missing.append(name)
    |    else:
    |        observed[key] = name
    |if missing:
    |    raise SystemExit(f"Expected output(s) {', '.join(missing)} not produced.")
    |data["observed_outputs"] = observed
    |meta.write_text(json.dumps(data, indent=2))
    |PY
    """.stripMargin().trim()
}

process RunTranscriptomeQC {
    tag {
        def parts = ["${datasetName}", "align${alignIdx}"]
        parts << ((mode == 'region') ? "region${regionIdx}" : mode)
        parts << "transcriptome${commandIdx}"
        parts.join("::")
    }
    conda RESOLVED_CONDA_ENV
    publishDir "results/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}${(mode == 'region') ? "/flair_regionalize${regionIdx}" : ""}/transcriptome/flair_transcriptome${commandIdx}", mode: 'copy'
    storeDir "cache/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}${(mode == 'region') ? "/flair_regionalize${regionIdx}" : ""}/transcriptome/flair_transcriptome${commandIdx}/qc"

    input:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx), val(regionSpec),
        val(mode), val(regionTag), val(commandIdx),
        path(stdoutFile), path(stderrFile), path(metadataFile),
        path(isoformsBed), path(isoformsGtf), path(isoformsFa), path(readMap),
        path(qcScript)

    output:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx),
        val(mode), val(regionTag), val(commandIdx), path('transcriptome_qc.tsv')

    script:
    def datasetEsc = escapeForSingleQuotes(datasetName)
    def modeEsc = escapeForSingleQuotes(mode?.toString() ?: '')
    def regionTagEsc = escapeForSingleQuotes(regionTag?.toString() ?: '')
    def stdoutNameEsc = escapeForSingleQuotes(stdoutFile.getFileName().toString())
    def stderrNameEsc = escapeForSingleQuotes(stderrFile.getFileName().toString())
    def metadataNameEsc = escapeForSingleQuotes(metadataFile.getFileName().toString())
    def isoformsBedEsc = escapeForSingleQuotes(isoformsBed.getFileName().toString())
    def isoformsGtfEsc = escapeForSingleQuotes(isoformsGtf.getFileName().toString())
    def isoformsFaEsc = escapeForSingleQuotes(isoformsFa.getFileName().toString())
    def readMapEsc = escapeForSingleQuotes(readMap.getFileName().toString())
    def qcScriptEsc = escapeForSingleQuotes(qcScript.toString())

    """
    |set -euo pipefail
    |
    |python ${qcScriptEsc} \\
    |    --dataset '${datasetEsc}' \\
    |    --align-idx ${alignIdx} \\
    |    --region-idx ${regionIdx != null ? regionIdx : 0} \\
    |    --mode '${modeEsc}' \\
    |    --region-tag '${regionTagEsc}' \\
    |    --stdout '${stdoutNameEsc}' \\
    |    --stderr '${stderrNameEsc}' \\
    |    --metadata '${metadataNameEsc}' \\
    |    --isoforms-bed '${isoformsBedEsc}' \\
    |    --isoforms-gtf '${isoformsGtfEsc}' \\
    |    --isoforms-fa '${isoformsFaEsc}' \\
    |    --read-map '${readMapEsc}' \\
    |    --output 'transcriptome_qc.tsv'
    """.stripMargin().trim()
}

process RunFlairCollapse {
    tag {
        def parts = ["${datasetName}", "align${alignIdx}"]
        parts << ((mode == 'region') ? "region${regionIdx}" : mode)
        parts << "collapse${commandIdx}"
        parts.join("::")
    }
    conda RESOLVED_CONDA_ENV
    publishDir "results/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}${(mode == 'region') ? "/flair_regionalize${regionIdx}" : ""}/collapse/flair_collapse${commandIdx}", mode: 'copy'
    storeDir "cache/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}${(mode == 'region') ? "/flair_regionalize${regionIdx}" : ""}/collapse/flair_collapse${commandIdx}"

    input:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx), val(regionSpec),
        val(mode), val(regionTag), val(commandIdx), val(commandTemplate), val(renderedCommand),
        path(collapseInputBed),
        val(readsPath),
        val(regionalizePayload),
        val(alignPayload)

    output:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx), val(regionSpec),
        val(mode), val(regionTag), val(commandIdx), val(commandTemplate), val(renderedCommand),
        path('collapse_stdout.txt'), path('collapse_stderr.txt'), path('collapse_metadata.json'),
        path('collapse_exit_code.txt'), path('collapse_error.log', optional: true),
        path('isoforms.bed'), path('isoforms.gtf'), path('isoforms.fa'),
        val(regionalizePayload), val(alignPayload)

    script:
    def bedPath = collapseInputBed as Path
    def bedName = bedPath.getFileName().toString()
    def bedNameEsc = escapeForSingleQuotes(bedName)

    def bedFile = bedPath.toFile()
    def bedEmpty = !bedFile.exists() || bedFile.length() == 0L
    def bedEmptyFlag = bedEmpty ? 1 : 0

    def readsPathValue = readsPath ? readsPath.toString().trim() : null
    if (!readsPathValue) {
        throw new IllegalArgumentException("Collapse stage requires a reads input for dataset '${datasetName}', align ${alignIdx}, region ${regionIdx}.")
    }
    def readsFile = new File(readsPathValue)
    if (!readsFile.exists()) {
        throw new IllegalArgumentException("Collapse reads input '${readsPathValue}' does not exist for dataset '${datasetName}'.")
    }
    def readsName = readsFile.name ?: "collapse_reads_input"
    def readsNameEsc = escapeForSingleQuotes(readsName)
    def readsSourceEsc = escapeForSingleQuotes(readsFile.toPath().toAbsolutePath().toString())

    def metadataMap = [
        dataset            : datasetName,
        align_index        : alignIdx,
        region_index       : regionIdx,
        region_mode        : mode,
        region_tag         : regionTag,
        command_index      : commandIdx,
        command_template   : commandTemplate,
        command_rendered   : renderedCommand,
        collapse_input_bed : bedName,
        collapse_reads     : readsName,
        collapse_reads_path: readsPathValue,
        collapse_input_empty: bedEmpty,
        region_spec        : regionSpec,
        runtime            : [
            conda_env      : RESOLVED_CONDA_ENV,
            conda_env_label: RESOLVED_CONDA_ENV_LABEL
        ]
    ]
    def metadataJson = JsonOutput.prettyPrint(JsonOutput.toJson(metadataMap))

    """
    |set -euo pipefail
    |
    |status=0
    |if [ ${bedEmptyFlag} -eq 1 ]; then
    |    printf 'Collapse input bed is empty; skipping flair collapse.\n' > collapse_stdout.txt
    |    : > collapse_stderr.txt
    |else
    |    ln -sf '${readsSourceEsc}' '${readsNameEsc}'
    |
    |    cat <<'CMD' > command_to_run.sh
    |#!/usr/bin/env bash
    |set -euo pipefail
    |${renderedCommand}
    |CMD
    |
    |    chmod +x command_to_run.sh
    |
    |    set +e
    |    ./command_to_run.sh > collapse_stdout.txt 2> collapse_stderr.txt
    |    status=\$?
    |    set -euo pipefail
    |fi
    |
    |echo "\${status}" > collapse_exit_code.txt
    |if [ "\${status}" -ne 0 ]; then
    |    printf 'Collapse command failed with exit code %s\n' "\${status}" > collapse_error.log
    |elif [ ${bedEmptyFlag} -eq 1 ]; then
    |    printf 'Collapse skipped: input bed was empty.\n' > collapse_error.log
    |else
    |    printf 'Collapse completed successfully.\n' > collapse_error.log
    |fi
    |
    |cat <<'EOF' > collapse_metadata.json
    |${metadataJson}
    |EOF
    |
    |python - <<'PY'
    |from pathlib import Path
    |import json
    |
    |meta = Path('collapse_metadata.json')
    |data = json.loads(meta.read_text())
    |stage = Path('.')
    |exit_code_path = stage / 'collapse_exit_code.txt'
    |exit_code = 0
    |if exit_code_path.exists():
    |    try:
    |        exit_code = int(exit_code_path.read_text().strip())
    |    except ValueError:
    |        exit_code = 0
    |data['exit_code'] = exit_code
    |error_log_path = stage / 'collapse_error.log'
    |data['error_log'] = error_log_path.name if error_log_path.exists() else None
    |bed_empty = bool(data.get('collapse_input_empty'))
    |expected = {
    |    'isoforms_bed': 'isoforms.bed',
    |    'isoforms_gtf': 'isoforms.gtf',
    |    'isoforms_fa': 'isoforms.fa',
    |}
    |observed = {}
    |missing = []
    |placeholder = []
    |for key, name in expected.items():
    |    path = stage / name
    |    if not path.exists():
    |        if exit_code != 0 or bed_empty:
    |            path.touch()
    |            observed[key] = name
    |            placeholder.append(name)
    |        else:
    |            missing.append(name)
    |            observed[key] = None
    |    else:
    |        observed[key] = name
    |if missing:
    |    raise SystemExit(f"Expected output(s) {', '.join(missing)} not produced.")
    |if placeholder:
    |    data['placeholder_outputs'] = placeholder
    |data['observed_outputs'] = observed
    |meta.write_text(json.dumps(data, indent=2))
    |PY
    """.stripMargin().trim()
}

process RunCollapseQC {
    tag {
        def parts = ["${datasetName}", "align${alignIdx}"]
        parts << ((mode == 'region') ? "region${regionIdx}" : mode)
        parts << "collapse${commandIdx}"
        parts.join("::")
    }
    conda RESOLVED_CONDA_ENV
    publishDir "results/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}${(mode == 'region') ? "/flair_regionalize${regionIdx}" : ""}/collapse/flair_collapse${commandIdx}", mode: 'copy'
    storeDir "cache/${RESOLVED_CONDA_ENV_LABEL}/${datasetName}/flair_align${alignIdx}${(mode == 'region') ? "/flair_regionalize${regionIdx}" : ""}/collapse/flair_collapse${commandIdx}/qc"

    input:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx), val(regionSpec),
        val(mode), val(regionTag), val(commandIdx),
        path(stdoutFile), path(stderrFile), path(metadataFile),
        path(isoformsBed), path(isoformsGtf), path(isoformsFa),
        path(qcScript)

    output:
    tuple val(datasetName), val(datasetSpec), val(alignIdx), val(regionIdx),
        val(mode), val(regionTag), val(commandIdx), path('collapse_qc.tsv')

    script:
    def datasetEsc = escapeForSingleQuotes(datasetName)
    def modeEsc = escapeForSingleQuotes(mode?.toString() ?: '')
    def regionTagEsc = escapeForSingleQuotes(regionTag?.toString() ?: '')
    def stdoutNameEsc = escapeForSingleQuotes(stdoutFile.getFileName().toString())
    def stderrNameEsc = escapeForSingleQuotes(stderrFile.getFileName().toString())
    def metadataNameEsc = escapeForSingleQuotes(metadataFile.getFileName().toString())
    def isoformsBedEsc = escapeForSingleQuotes(isoformsBed.getFileName().toString())
    def isoformsGtfEsc = escapeForSingleQuotes(isoformsGtf.getFileName().toString())
    def isoformsFaEsc = escapeForSingleQuotes(isoformsFa.getFileName().toString())
    def qcScriptEsc = escapeForSingleQuotes(qcScript.toString())

    """
    |set -euo pipefail
    |
    |python ${qcScriptEsc} \\
    |    --dataset '${datasetEsc}' \\
    |    --align-idx ${alignIdx} \\
    |    --region-idx ${regionIdx != null ? regionIdx : 0} \\
    |    --mode '${modeEsc}' \\
    |    --region-tag '${regionTagEsc}' \\
    |    --stdout '${stdoutNameEsc}' \\
    |    --stderr '${stderrNameEsc}' \\
    |    --metadata '${metadataNameEsc}' \\
    |    --isoforms-bed '${isoformsBedEsc}' \\
    |    --isoforms-gtf '${isoformsGtfEsc}' \\
    |    --isoforms-fa '${isoformsFaEsc}' \\
    |    --output 'collapse_qc.tsv'
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
        def regionalize_results_for_transcriptome = regionalize_results_channel.map { it }
        def regionalize_results_for_correct = regionalize_results_channel.map { it }
        regionalize_results = regionalize_results_channel.map { it }
        regionalize_qc_results = regionalize_qc_results_channel

        def transcriptomeCommandsRaw = normalizeCommands(params.flair_transcriptome_commands)
        if (!transcriptomeCommandsRaw) {
            log.info("No --flair_transcriptome_commands supplied; skipping flair transcriptome runs.")
        }

        def transcriptomeCommandEntries = []
        transcriptomeCommandsRaw.eachWithIndex { cmd, idx ->
            transcriptomeCommandEntries << [idx + 1, cmd]
        }

        def transcriptome_results_channel = Channel.empty()
        def transcriptome_qc_results_channel = Channel.empty()
        if (!transcriptomeCommandEntries.isEmpty()) {
            def transcriptome_inputs = regionalize_results_for_transcriptome.flatMap { payload ->
                transcriptomeCommandEntries.collect { info ->
                    def commandIdx = info[0] as int
                    def template = info[1] as String

                    def regionSpec = new LinkedHashMap<String, Object>(payload.regionSpec ?: [:])
                    def alignPayload = payload.alignPayload

                    List<Path> bamCandidates = flattenPathList(payload.bamFiles ?: [])
                    if (!bamCandidates && alignPayload) {
                        bamCandidates = flattenPathList(alignPayload.bamFiles ?: [])
                    }
                    if (!bamCandidates) {
                        throw new IllegalStateException("Transcriptome stage requires a BAM file for dataset ${payload.datasetName}, align ${payload.alignIndex}, region ${payload.regionIndex}.")
                    }
                    Path bamFile = null
                    if (payload.mode == 'region') {
                        bamFile = bamCandidates.find { it.getFileName().toString() ==~ /(?i)[A-Za-z0-9._-]+_\d+_\d+\.bam$/ }
                    }
                    bamFile = bamFile ?: bamCandidates[0]

                    def regionTag = payload.regionTag
                    if (!regionTag) {
                        if (payload.mode == 'region' && payload.regionIndex != null) {
                            regionTag = "align${payload.alignIndex}_region${payload.regionIndex}"
                        } else {
                            regionTag = "align${payload.alignIndex}_all"
                        }
                    }

                    def genomePath = payload.datasetSpec?.genome?.toString()
                    if (!genomePath) {
                        throw new IllegalStateException("Transcriptome stage requires 'genome' in dataset specification for dataset ${payload.datasetName}.")
                    }

                    def substitution = new LinkedHashMap<String, Object>(payload.datasetSpec ?: [:])
                    substitution['transcriptome_bam'] = bamFile.getFileName().toString()
                    substitution['transcriptome_bam_path'] = bamFile.toAbsolutePath().toString()
                    substitution['transcriptome_output_prefix'] = regionTag
                    substitution['transcriptome_mode'] = payload.mode
                    substitution['transcriptome_region'] = regionTag
                    substitution['align_index'] = payload.alignIndex
                    substitution['region_index'] = payload.regionIndex
                    def rendered = renderFlairCommand(template, substitution, payload.datasetName)

                    def regionalizeMetadata = payload.metadataPath
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
                        bamFile,
                        regionalizeMetadata,
                        alignMetadata,
                        payload,
                        alignPayload
                    )
                }
            }

            def raw_transcriptome_results = RunFlairTranscriptome(transcriptome_inputs)
            def transcriptome_results_for_map = raw_transcriptome_results.map { it }
            def transcriptome_results_for_qc = raw_transcriptome_results.map { it }

            transcriptome_results_channel = transcriptome_results_for_map.map { values ->
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
                def isoformsBed = values[13]
                def isoformsGtf = values[14]
                def isoformsFa = values[15]
                def readMap = values[16]
                def exitCodePath = values[17]
                def errorLogPath = (values.size() > 18) ? values[18] : null
                def regionalizePayload = values.size() > 19 ? values[19] : null
                def alignPayload = values.size() > 20 ? values[20] : null
                [
                    datasetName       : datasetName,
                    datasetSpec       : datasetSpec,
                    alignIndex        : alignIdx,
                    regionIndex       : regionIdx,
                    regionSpec        : regionSpec,
                    mode              : mode,
                    regionTag         : regionTag,
                    commandIndex      : commandIdx,
                    commandTemplate   : commandTemplate,
                    renderedCommand   : renderedCommand,
                    stdoutPath        : stdoutFile,
                    stderrPath        : stderrFile,
                    metadataPath      : metadataFile,
                    isoformsBedPath   : isoformsBed,
                    isoformsGtfPath   : isoformsGtf,
                    isoformsFaPath    : isoformsFa,
                    readMapPath       : readMap,
                    exitCodePath      : exitCodePath,
                    errorLogPath      : errorLogPath,
                    regionalizePayload: regionalizePayload,
                    alignPayload      : alignPayload
                ]
            }

            def transcriptome_qc_inputs = transcriptome_results_for_qc.map { values ->
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
                def isoformsBed = values[13]
                def isoformsGtf = values[14]
                def isoformsFa = values[15]
                def readMap = values[16]
                def qcScript = file("${projectDir}/bin/transcriptome_qc.py")
                tuple(datasetName, datasetSpec, alignIdx, regionIdx, regionSpec, mode, regionTag, commandIdx, stdoutFile, stderrFile, metadataFile, isoformsBed, isoformsGtf, isoformsFa, readMap, qcScript)
            }
            transcriptome_qc_results_channel = RunTranscriptomeQC(transcriptome_qc_inputs)
        }
        transcriptome_results = transcriptome_results_channel
        transcriptome_qc_results = transcriptome_qc_results_channel

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
        def correct_results_for_collapse = Channel.empty()
        if (!correctCommandEntries.isEmpty()) {
            def correct_inputs = regionalize_results_for_correct.flatMap { payload ->
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
                        payload,
                        alignPayload
                    )
                }
            }
            def raw_correct_results = RunFlairCorrect(correct_inputs)
            def correct_results_for_map = raw_correct_results.map { it }
            def correct_results_for_qc = raw_correct_results.map { it }
            correct_results_for_collapse = raw_correct_results.map { it }

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
                def regionalizePayload = values.size() > 16 ? values[16] : null
                def alignPayload = values.size() > 17 ? values[17] : null
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
                    regionalizePayload: regionalizePayload,
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

        def collapseCommandsRaw = normalizeCommands(params.flair_collapse_commands)
        if (!collapseCommandsRaw) {
            log.info("No --flair_collapse_commands supplied; skipping flair collapse runs.")
        }

        def collapseCommandEntries = []
        collapseCommandsRaw.eachWithIndex { cmd, idx ->
            collapseCommandEntries << [idx + 1, cmd]
        }

        def collapse_results_channel = Channel.empty()
        def collapse_qc_results_channel = Channel.empty()
        if (!collapseCommandEntries.isEmpty()) {
            def collapse_inputs = correct_results_for_collapse.flatMap { values ->
                collapseCommandEntries.collect { info ->
                    def commandIdx = info[0] as int
                    def template = info[1] as String

                    def datasetName = values[0]
                    def datasetSpec = values[1]
                    def alignIdx = values[2]
                    def regionIdx = values[3]
                    def regionSpec = values[4]
                    def mode = values[5]
                    def regionTag = values[6]
                    def correctedBeds = flattenPathList(values.size() > 13 ? values[13] : [])
                    if (!correctedBeds) {
                        throw new IllegalStateException("Collapse stage requires a corrected BED file for dataset ${datasetName}, align ${alignIdx}, region ${regionIdx}.")
                    }
                    Path correctedBed = correctedBeds[0]

                    def regionalizePayload = values.size() > 16 ? values[16] : null
                    def alignPayload = values.size() > 17 ? values[17] : null

                    def readsPath = datasetSpec?.reads?.toString()
                    if (!readsPath && alignPayload) {
                        readsPath = alignPayload?.datasetSpec?.reads?.toString()
                    }
                    if (!readsPath && regionalizePayload) {
                        readsPath = regionalizePayload?.datasetSpec?.reads?.toString()
                    }
                    if (!readsPath) {
                        throw new IllegalStateException("Collapse stage requires 'reads' in dataset specification for dataset ${datasetName}.")
                    }

                    def substitution = new LinkedHashMap<String, Object>(datasetSpec ?: [:])
                    substitution['collapse_bed'] = correctedBed.getFileName().toString()
                    substitution['collapse_reads'] = new File(readsPath).name
                    substitution['collapse_output_prefix'] = regionTag ?: "align${alignIdx}${(mode == 'region' && regionIdx != null) ? "_region${regionIdx}" : "_all"}"
                    substitution['collapse_mode'] = mode
                    substitution['collapse_region'] = regionTag ?: ''
                    substitution['align_index'] = alignIdx
                    substitution['region_index'] = regionIdx
                    def rendered = renderFlairCommand(template, substitution, datasetName)

                    tuple(
                        datasetName,
                        datasetSpec,
                        alignIdx,
                        regionIdx,
                        regionSpec,
                        mode,
                        regionTag,
                        commandIdx,
                        template,
                        rendered,
                        correctedBed,
                        readsPath,
                        regionalizePayload,
                        alignPayload
                    )
                }
            }

            def raw_collapse_results = RunFlairCollapse(collapse_inputs)
            def collapse_results_for_map = raw_collapse_results.map { it }
            def collapse_results_for_qc = raw_collapse_results.map { it }

            collapse_results_channel = collapse_results_for_map.map { values ->
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
                def exitCodeFile = values[13]
                def errorLogFile = values[14]
                def isoformsBed = values[15]
                def isoformsGtf = values[16]
                def isoformsFa = values[17]
                def regionalizePayload = values.size() > 18 ? values[18] : null
                def alignPayload = values.size() > 19 ? values[19] : null
                [
                    datasetName       : datasetName,
                    datasetSpec       : datasetSpec,
                    alignIndex        : alignIdx,
                    regionIndex       : regionIdx,
                    regionSpec        : regionSpec,
                    mode              : mode,
                    regionTag         : regionTag,
                    commandIndex      : commandIdx,
                    commandTemplate   : commandTemplate,
                    renderedCommand   : renderedCommand,
                    stdoutPath        : stdoutFile,
                    stderrPath        : stderrFile,
                    metadataPath      : metadataFile,
                    exitCodePath      : exitCodeFile,
                    errorLogPath      : errorLogFile,
                    isoformsBedPath   : isoformsBed,
                    isoformsGtfPath   : isoformsGtf,
                    isoformsFaPath    : isoformsFa,
                    regionalizePayload: regionalizePayload,
                    alignPayload      : alignPayload
                ]
            }

            def collapse_qc_inputs = collapse_results_for_qc.map { values ->
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
                def isoformsBed = values[15]
                def isoformsGtf = values[16]
                def isoformsFa = values[17]
                def qcScript = file("${projectDir}/bin/collapse_qc.py")
                tuple(datasetName, datasetSpec, alignIdx, regionIdx, regionSpec, mode, regionTag, commandIdx, stdoutFile, stderrFile, metadataFile, isoformsBed, isoformsGtf, isoformsFa, qcScript)
            }
            collapse_qc_results_channel = RunCollapseQC(collapse_qc_inputs)
        }
        collapse_results = collapse_results_channel
        collapse_qc_results = collapse_qc_results_channel

    emit:
        dataset_manifests
        align_results
        align_qc_results
        regionalize_results
        regionalize_qc_results
        transcriptome_results
        transcriptome_qc_results
        collapse_results
        collapse_qc_results
        correct_results
        correct_qc_results
}

workflow {
    flair_eval()
}
