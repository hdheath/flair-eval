/*
 * Input parsing for flair-eval v2 (docs/CONTRACTS.md section 1).
 *
 *   loadInputs(params) -> [samples: [...], runs: [...]]
 *
 * Samples come from --input samples.csv or from the one-command flags (--reads/--bam --genome
 * --gtf ...). Runs come from --runs runs.yaml or from --flair tokens. Every problem found is
 * collected and reported at once, in plain English. The JSON schemas in assets/schema/ are then
 * applied as a final check, so this file and the schemas can never silently disagree.
 */

include { validate                } from 'plugin/nf-schema'
include { safeId                  } from './keys'

def DEFAULT_FLAIR_GIT() { 'https://github.com/BrooksLabUCSC/flair' }
def ID_PATTERN()        { '^[A-Za-z0-9][A-Za-z0-9._-]*$' }
def SAMPLE_COLUMNS()    { ['sample', 'reads', 'bam', 'genome', 'gtf', 'platform', 'library', 'cage', 'pas', 'junctions', 'spikein_gtf'] }
def SAMPLE_FILE_COLUMNS() { ['reads', 'bam', 'genome', 'gtf', 'cage', 'pas', 'junctions', 'spikein_gtf'] }
def QUICK_SAMPLE_PARAMS() { ['reads', 'bam', 'genome', 'gtf', 'sample', 'platform', 'library'] }
def RUN_KEYS()          { ['id', 'flair', 'transcriptome_args', 'quantify', 'productivity', 'own_align', 'label', 'tags'] }
def FLAIR_SOURCES()     { ['release', 'git', 'env', 'container'] }

/*
 * ---------------------------------------------------------------------------------------------
 * Entry point
 * ---------------------------------------------------------------------------------------------
 */
def loadInputs(Map p) {
    def errors = []

    // Mutual exclusions (CONTRACTS section 1).
    def quickGiven = QUICK_SAMPLE_PARAMS().findAll { k -> hasParam(p, k) }
    if (hasParam(p, 'input') && quickGiven) {
        errors << "--input cannot be combined with the one-command sample flags (${quickGiven.collect { k -> '--' + k }.join(', ')}). Put the sample in samples.csv, or drop --input."
    }
    if (hasParam(p, 'runs') && hasParam(p, 'flair')) {
        errors << "--runs and --flair cannot be used together. List every run in runs.yaml, or drop --runs and use --flair tokens."
    }
    if (hasParam(p, 'runs') && hasParam(p, 'flair_args')) {
        errors << "--flair_args only applies to runs built from --flair. With --runs, set transcriptome_args on each run in runs.yaml."
    }
    if (!hasParam(p, 'input') && !quickGiven) {
        errors << "No samples given. Use --input samples.csv, or --reads (or --bam) with --genome and --gtf."
    }
    if (!hasParam(p, 'runs') && !hasParam(p, 'flair')) {
        errors << "No FLAIR runs given. Use --runs runs.yaml, or --flair with comma-separated tokens, e.g. --flair 3.0.1,dev."
    }
    if (errors) {
        failInputs(errors)
    }

    def samples = hasParam(p, 'input') ? loadSamplesCsv(p.input.toString(), errors) : quickSample(p, errors)
    def runs    = hasParam(p, 'runs')  ? loadRunsYaml(p.runs.toString(), errors)    : quickRuns(p, errors)

    // Cross-checks between samples and runs.
    def ownAlign = runs.findAll { r -> r.own_align }
    def bamOnly  = samples.findAll { s -> s.bam && !s.reads }
    if (ownAlign && bamOnly) {
        errors << "Run(s) ${ownAlign.collect { r -> r.id }.join(', ')} set own_align: true, but sample(s) ${bamOnly.collect { s -> s.sample }.join(', ')} only have a BAM, so there are no reads to align. Add reads for those samples, or set own_align: false."
    }
    if (errors) {
        failInputs(errors)
    }
    return [samples: samples, runs: runs]
}

def hasParam(Map p, String k) {
    return p[k] != null && p[k] != false && p[k].toString() != ''
}

def failInputs(List errors) {
    error("flair-eval found ${errors.size()} problem${errors.size() == 1 ? '' : 's'} with the inputs:\n  - " + errors.join('\n  - ') + '\n')
}

/*
 * ---------------------------------------------------------------------------------------------
 * samples.csv
 * ---------------------------------------------------------------------------------------------
 */
def loadSamplesCsv(String csvPath, List errors) {
    def csv = file(csvPath)
    if (!csv.exists()) {
        errors << "samples file '${csvPath}' does not exist."
        failInputs(errors)
    }
    def lines = csv.readLines().findAll { l -> l.trim() && !l.trim().startsWith('#') }
    if (!lines) {
        errors << "samples file '${csvPath}' is empty. It needs a header line and one row per sample."
        failInputs(errors)
    }
    def header = splitCsvLine(lines[0]).collect { h -> h.trim() }
    def unknown = header.findAll { h -> h && !(h in SAMPLE_COLUMNS()) }
    if (unknown) {
        errors << "samples file '${csvPath}' has unknown column(s) ${unknown.collect { u -> "'${u}'" }.join(', ')}. Allowed columns: ${SAMPLE_COLUMNS().join(', ')}."
    }
    ['sample', 'genome', 'gtf'].each { c ->
        if (!(c in header)) {
            errors << "samples file '${csvPath}' has no '${c}' column (required)."
        }
    }
    if (!('reads' in header) && !('bam' in header)) {
        errors << "samples file '${csvPath}' needs a 'reads' column, a 'bam' column, or both."
    }
    if (errors) {
        failInputs(errors)
    }
    if (lines.size() < 2) {
        errors << "samples file '${csvPath}' has a header but no sample rows."
        failInputs(errors)
    }

    def baseDir = csv.toAbsolutePath().parent
    def seen = [:]
    def samples = []
    lines.drop(1).eachWithIndex { line, i ->
        def rowNo = i + 2   // 1-based line number, header is line 1
        def cells = splitCsvLine(line)
        if (cells.size() > header.size()) {
            errors << "samples file line ${rowNo} has ${cells.size()} values but the header has ${header.size()} columns."
            return
        }
        def row = [:]
        header.eachWithIndex { h, j ->
            def v = j < cells.size() ? cells[j].trim() : ''
            if (h && v) {
                row[h] = v
            }
        }
        def where = row.sample ? "samples file line ${rowNo} (sample '${row.sample}')" : "samples file line ${rowNo}"
        def nBefore = errors.size()
        def s = checkSample(row, where, baseDir, errors)
        if (s && seen.containsKey(s.sample)) {
            errors << "${where}: sample name '${s.sample}' is already used on line ${seen[s.sample]}. Sample names must be unique."
        }
        if (s) {
            seen[s.sample] = rowNo
        }
        if (errors.size() == nBefore) {
            schemaCheck(row, 'samples.schema.json', where, errors)
            samples << s
        }
    }
    return samples
}

// Minimal RFC-4180 split: commas, double-quoted fields, "" escapes.
def splitCsvLine(String line) {
    def st = line.replaceAll('\r$', '').toList().inject([out: [], cur: '', q: false, closed: false]) { acc, ch ->
        if (acc.q) {
            if (ch == '"') {
                acc.q = false
                acc.closed = true
            }
            else {
                acc.cur += ch
            }
        }
        else if (ch == '"') {
            if (acc.closed) {
                acc.cur += '"'      // "" inside a quoted field
            }
            acc.q = true
            acc.closed = false
        }
        else if (ch == ',') {
            acc.out << acc.cur
            acc.cur = ''
            acc.closed = false
        }
        else {
            acc.cur += ch
            acc.closed = false
        }
        acc
    }
    return st.out + [st.cur]
}

/*
 * Check one sample (from a CSV row or the one-command flags). Returns the normalised sample map
 * with resolved file objects, or null when a required value is missing.
 */
/*
 * Genome presets (params.genomes, lab references on phoenix): a genome value that is a preset key
 * becomes that preset's FASTA, and an empty gtf becomes the preset's GTF. Returns false on error.
 */
def genomePresetKeys() {
    def g = params.genomes instanceof Map ? params.genomes : [:]
    return g ? g.keySet().join(', ') : 'none'
}

def expandGenomePreset(Map row, String where, List errors) {
    def presets = params.genomes instanceof Map ? params.genomes : [:]
    def key = row.genome ? row.genome.toString() : ''
    if (!presets.containsKey(key)) {
        return true
    }
    def preset = presets[key]
    ['fasta', 'gtf'].each { k ->
        if (!file(preset[k].toString()).exists()) {
            errors << "${where}: genome preset '${key}' needs ${preset[k]}, which is not on this machine. The presets are lab references on phoenix; elsewhere pass --genome and --gtf files."
        }
    }
    row.genome = preset.fasta.toString()
    if (!row.gtf) {
        row.gtf = preset.gtf.toString()
    }
    return !errors.any { e -> e.startsWith("${where}: genome preset") }
}

def checkSample(Map row, String where, Object baseDir, List errors) {
    def ok = true
    if (!row.sample) {
        errors << "${where}: 'sample' is empty. Give each sample a name."
        ok = false
    }
    else if (!(row.sample ==~ ID_PATTERN())) {
        errors << "${where}: sample name '${row.sample}' may only use letters, digits, '.', '_' and '-', and must start with a letter or digit."
        ok = false
    }
    if (!row.reads && !row.bam) {
        errors << "${where}: give 'reads' (FASTQ/FASTA) or 'bam' (sorted and indexed), or both."
        ok = false
    }
    if (!expandGenomePreset(row, where, errors)) {
        ok = false
    }
    ['genome', 'gtf'].each { c ->
        if (!row[c]) {
            errors << "${where}: '${c}' is required" + (c == 'genome' ? " (a FASTA path, or a preset: ${genomePresetKeys()})." : '.')
            ok = false
        }
    }
    if (row.platform && !(row.platform in ['pacbio', 'ont'])) {
        errors << "${where}: platform '${row.platform}' is not allowed. Use 'pacbio' or 'ont' (lower case), or leave it empty."
        ok = false
    }
    if (row.library && !(row.library in ['cDNA', 'dRNA'])) {
        errors << "${where}: library '${row.library}' is not allowed. Use 'cDNA' or 'dRNA', or leave it empty."
        ok = false
    }
    def files = [:]
    SAMPLE_FILE_COLUMNS().each { c ->
        if (row[c]) {
            def f = resolveInputPath(row[c].toString(), baseDir)
            if (!f.exists()) {
                errors << "${where}: ${c} file '${row[c]}' does not exist${isRelative(row[c].toString()) && baseDir ? " (looked next to the samples file and in the launch directory)" : ''}${c == 'genome' && !row[c].toString().contains('/') ? "; genome presets: ${genomePresetKeys()}" : ''}."
                ok = false
            }
            files[c] = f
        }
    }
    def bai = null
    if (files.bam && files.bam.exists()) {
        bai = findBamIndex(files.bam)
        if (bai == null) {
            errors << "${where}: BAM '${row.bam}' has no index. Run 'samtools index ${row.bam}' so '${files.bam.name}.bai' sits next to it."
            ok = false
        }
    }
    if (!ok) {
        return null
    }
    return [
        sample     : row.sample.toString(),
        reads      : files.reads,
        bam        : files.bam,
        bai        : bai,
        genome     : files.genome,
        gtf        : files.gtf,
        platform   : row.platform ?: null,
        library    : row.library ?: null,
        cage       : files.cage,
        pas        : files.pas,
        junctions  : files.junctions,
        spikein_gtf: files.spikein_gtf,
    ]
}

def isRelative(String v) {
    return !(v.contains('://') || v.startsWith('/'))
}

// Relative paths in samples.csv resolve next to the CSV first, then in the launch directory.
def resolveInputPath(String v, Object baseDir) {
    if (!isRelative(v) || baseDir == null) {
        return file(v)
    }
    def nextToCsv = file(baseDir.resolve(v).toString())
    return nextToCsv.exists() ? nextToCsv : file(v)
}

def findBamIndex(Object bam) {
    def candidates = [bam.toString() + '.bai', bam.toString().replaceAll(/\.bam$/, '') + '.bai', bam.toString() + '.csi']
    return candidates.collect { c -> file(c) }.find { f -> f.exists() }
}

/*
 * ---------------------------------------------------------------------------------------------
 * One-command sample
 * ---------------------------------------------------------------------------------------------
 */
def quickSample(Map p, List errors) {
    def row = [:]
    ['reads', 'bam', 'genome', 'gtf', 'platform', 'library'].each { k ->
        if (p[k] != null && p[k].toString() != '') {
            row[k] = p[k].toString()
        }
    }
    def nameFrom = row.reads ?: row.bam
    row.sample = p.sample ? p.sample.toString() : (nameFrom ? stripReadExtensions(file(nameFrom).name) : null)
    def where = 'one-command sample'
    if (!p.sample && row.sample && !(row.sample ==~ ID_PATTERN())) {
        errors << "${where}: the sample name taken from the file name ('${row.sample}') is not valid. Pass --sample NAME using letters, digits, '.', '_' or '-'."
        return []
    }
    def s = checkSample(row, where, null, errors)
    if (s) {
        schemaCheck(row, 'samples.schema.json', where, errors)
    }
    return s ? [s] : []
}

def READ_EXTENSIONS() { ['.gz', '.bgz', '.bz2', '.xz', '.zst', '.fastq', '.fq', '.fasta', '.fa', '.fna', '.bam', '.cram', '.sam'] }

def stripReadExtensions(String name) {
    def ext = READ_EXTENSIONS().find { e -> name.toLowerCase().endsWith(e) && name.length() > e.length() }
    return ext ? stripReadExtensions(name.substring(0, name.length() - ext.length())) : name
}

/*
 * ---------------------------------------------------------------------------------------------
 * runs.yaml
 * ---------------------------------------------------------------------------------------------
 */
def loadRunsYaml(String yamlPath, List errors) {
    def f = file(yamlPath)
    if (!f.exists()) {
        errors << "runs file '${yamlPath}' does not exist."
        failInputs(errors)
    }
    def doc = null
    try {
        doc = new org.yaml.snakeyaml.Yaml().load(f.text)
    }
    catch (Exception e) {
        errors << "runs file '${yamlPath}' is not valid YAML: ${e.message?.readLines()?.take(3)?.join(' ')}"
        failInputs(errors)
    }
    if (!(doc instanceof Map) || !doc.containsKey('runs')) {
        errors << "runs file '${yamlPath}' must start with 'runs:' followed by a list of runs (see docs/CONTRACTS.md)."
        failInputs(errors)
    }
    def extraTop = doc.keySet().findAll { k -> k != 'runs' }
    if (extraTop) {
        errors << "runs file '${yamlPath}' has unknown top-level key(s) ${extraTop.join(', ')}. Only 'runs' is allowed."
    }
    if (!(doc.runs instanceof List) || !doc.runs) {
        errors << "runs file '${yamlPath}': 'runs' must be a non-empty list, each entry starting with '- id: ...'."
        failInputs(errors)
    }
    def runs = []
    def seen = [:]
    doc.runs.eachWithIndex { r, i ->
        def where = (r instanceof Map && r.id) ? "runs file entry ${i + 1} (id '${r.id}')" : "runs file entry ${i + 1}"
        def nBefore = errors.size()
        def run = checkRun(r, where, errors)
        if (run && seen.containsKey(run.id)) {
            errors << "runs file entries ${seen[run.id]} and ${i + 1} both use run id '${run.id}'. Run ids must be unique: give each run in runs.yaml its own explicit id."
        }
        if (run) {
            seen[run.id] = i + 1
        }
        if (errors.size() == nBefore && run) {
            runs << run
        }
    }
    if (errors.size() == 0) {
        schemaCheck(doc, 'runs.schema.json', "runs file '${yamlPath}'", errors)
    }
    return runs
}

def checkRun(Object r, String where, List errors) {
    if (!(r instanceof Map)) {
        errors << "${where}: each run must be a mapping with at least 'id' and 'flair'."
        return null
    }
    def ok = true
    def unknown = r.keySet().findAll { k -> !(k in RUN_KEYS()) }
    if (unknown) {
        errors << "${where}: unknown key(s) ${unknown.join(', ')}. Allowed: ${RUN_KEYS().join(', ')}."
        ok = false
    }
    if (r.id == null || r.id.toString() == '') {
        errors << "${where}: 'id' is required."
        ok = false
    }
    else if (!(r.id instanceof String)) {
        errors << "${where}: id ${r.id} must be text. Put it in quotes, e.g. id: \"${r.id}\"."
        ok = false
    }
    else if (!(r.id ==~ ID_PATTERN())) {
        errors << "${where}: id '${r.id}' may only use letters, digits, '.', '_' and '-', and must start with a letter or digit."
        ok = false
    }
    if (!checkFlairSource(r.flair, where, errors)) {
        ok = false
    }
    ['transcriptome_args', 'label'].each { k ->
        if (r.containsKey(k) && !(r[k] instanceof String)) {
            errors << "${where}: '${k}' must be text. Put it in quotes."
            ok = false
        }
    }
    ['quantify', 'productivity', 'own_align'].each { k ->
        if (r.containsKey(k) && !(r[k] instanceof Boolean)) {
            errors << "${where}: '${k}' must be true or false, not '${r[k]}'."
            ok = false
        }
    }
    if (r.containsKey('tags')) {
        if (!(r.tags instanceof Map)) {
            errors << "${where}: 'tags' must be a mapping of names to text, e.g. tags: { version: \"3.0.1\" }."
            ok = false
        }
        else {
            def bad = r.tags.findAll { _k, v -> !(v instanceof String) }
            if (bad) {
                errors << "${where}: tag value(s) for ${bad.keySet().join(', ')} must be text. Put them in quotes, e.g. ${bad.keySet()[0]}: \"${bad.values()[0]}\"."
                ok = false
            }
        }
    }
    if (!ok) {
        return null
    }
    return normaliseRun(r)
}

def checkFlairSource(Object fl, String where, List errors) {
    def hint = "Use exactly one of: release: \"3.0.1\" | git: <url> with ref: <branch/tag/sha> | env: /abs/conda/env | container: <uri or .sif>."
    if (!(fl instanceof Map)) {
        errors << "${where}: 'flair' is required and must be a mapping. ${hint}"
        return false
    }
    def given = FLAIR_SOURCES().findAll { k -> fl.containsKey(k) }
    def extra = fl.keySet().findAll { k -> !(k in FLAIR_SOURCES() + ['ref']) }
    if (extra) {
        errors << "${where}: flair has unknown key(s) ${extra.join(', ')}. ${hint}"
        return false
    }
    if (given.size() != 1) {
        errors << "${where}: flair gives ${given ? given.join(' and ') : 'no source'}. ${hint}"
        return false
    }
    def src = given[0]
    def v = fl[src]
    if (!(v instanceof String) || !v) {
        errors << "${where}: flair ${src} must be non-empty text. Put it in quotes, e.g. ${src}: \"${v ?: '...'}\"" + (v instanceof Number ? " (YAML read ${v} as a number)." : '.')
        return false
    }
    if (src == 'git') {
        if (!fl.containsKey('ref') || !(fl.ref instanceof String) || !fl.ref) {
            errors << "${where}: flair git needs 'ref' (a branch, tag or commit SHA), as text in quotes."
            return false
        }
    }
    else if (fl.containsKey('ref')) {
        errors << "${where}: 'ref' only goes with 'git'."
        return false
    }
    if (src == 'env') {
        if (!v.startsWith('/')) {
            errors << "${where}: flair env '${v}' must be an absolute path to a conda env."
            return false
        }
        if (!file(v).exists()) {
            errors << "${where}: flair env '${v}' does not exist."
            return false
        }
    }
    if (src == 'container' && !v.contains('://') && v.startsWith('/') && !file(v).exists()) {
        errors << "${where}: container image '${v}' does not exist."
        return false
    }
    return true
}

def normaliseRun(Map r) {
    def fl = r.flair.collectEntries { k, v -> [(k): v.toString()] }
    return [
        id                : r.id.toString(),
        flair             : fl,
        transcriptome_args: (r.transcriptome_args ?: '').toString(),
        quantify          : r.containsKey('quantify') ? r.quantify : true,
        productivity      : r.containsKey('productivity') ? r.productivity : true,
        own_align         : r.containsKey('own_align') ? r.own_align : false,
        label             : (r.label ?: r.id).toString(),
        tags              : (r.tags ?: [:]).collectEntries { k, v -> [(k.toString()): v.toString()] },
    ]
}

/*
 * ---------------------------------------------------------------------------------------------
 * --flair tokens (CONTRACTS section 1 table)
 * ---------------------------------------------------------------------------------------------
 */
def quickRuns(Map p, List errors) {
    def value = p.flair == null ? '' : p.flair.toString()
    if (!value.trim()) {
        errors << "--flair is empty; give at least one token, e.g. --flair 3.0.1"
        return []
    }
    def runs = []
    def seen = [:]
    value.split(',', -1).each { tok ->
        def t = parseToken(tok)
        if (t.error) {
            errors << t.error
            return
        }
        if (seen.containsKey(t.id)) {
            errors << "--flair tokens '${seen[t.id]}' and '${t.token}' both give run id '${t.id}'; write these runs in runs.yaml with distinct ids"
            return
        }
        seen[t.id] = t.token
        def nBefore = errors.size()
        if (!checkFlairSource(t.flair, "--flair token '${t.token}'", errors) || errors.size() != nBefore) {
            return
        }
        runs << normaliseRun([id: t.id, flair: t.flair, transcriptome_args: (p.flair_args ?: '').toString()])
    }
    if (errors.size() == 0) {
        schemaCheck([runs: runs], 'runs.schema.json', '--flair runs', errors)
    }
    return runs
}

/*
 * One --flair token -> [token, id, flair] or [token, error]. Mirrors cli/src/flair_eval/tokens.py
 * (parse_token) exactly, including its messages (CONTRACTS section 7, A1).
 */
def parseToken(String token) {
    def t = token.trim()
    if (!t) {
        return [token: t, error: 'empty --flair token (check for a stray comma)']
    }
    if (t.startsWith('release:')) {
        def v = t.substring('release:'.length()).trim()
        if (!v) {
            return [token: t, error: "--flair token '${t}': give a version, e.g. release:3.0.1".toString()]
        }
        return tokenRun(t, 'flair-', v, [release: v])
    }
    if (t.startsWith('git:')) {
        def rest = t.substring('git:'.length()).trim()
        if (!rest) {
            return [token: t, error: "--flair token '${t}': give a branch, tag or SHA, e.g. git:dev".toString()]
        }
        def url = DEFAULT_FLAIR_GIT()
        def ref = rest
        if (looksLikeUrl(rest)) {
            def at = rest.lastIndexOf('@')
            url = at > 0 ? rest.substring(0, at) : ''
            ref = at > 0 ? rest.substring(at + 1) : ''
            if (!looksLikeUrl(url) || !ref || (ref.contains('/') && ref.contains(':'))) {
                return [token: t, error: "--flair token '${t}': a git URL needs @<ref>, e.g. git:https://github.com/me/flair@my-branch".toString()]
            }
        }
        return tokenRun(t, 'flair-', ref, [git: url, ref: ref])
    }
    if (t.startsWith('env:')) {
        def path = t.substring('env:'.length()).trim()
        if (!path) {
            return [token: t, error: "--flair token '${t}': give a conda env path, e.g. env:/path/to/env".toString()]
        }
        return tokenRun(t, 'env-', tokenBasename(path), [env: path])
    }
    if (t.startsWith('container:')) {
        def uri = t.substring('container:'.length()).trim()
        if (!uri) {
            return [token: t, error: "--flair token '${t}': give an image, e.g. container:docker://org/flair:3.0.1 or container:/path/flair.sif".toString()]
        }
        return tokenRun(t, 'ctr-', containerName(uri), [container: uri])
    }
    if (t in ['dev', 'main']) {
        return [token: t, id: 'flair-' + t, flair: [git: DEFAULT_FLAIR_GIT(), ref: t]]
    }
    if (t ==~ /^[0-9].*/) {
        return tokenRun(t, 'flair-', t, [release: t])
    }
    return [token: t, error: "unknown --flair token '${t}'. Use a release (3.0.1), dev, main, git:${t} for another branch, env:/path or container:<uri>".toString()]
}

// Test helper: parseToken over a list (one nf-test call for a whole table).
def quickRunsCheck(List values, String flairArgs) {
    return values.collect { v ->
        def errors = []
        def runs = quickRuns([flair: v, flair_args: flairArgs], errors)
        [ids: runs.collect { r -> r.id }, args: runs.collect { r -> r.transcriptome_args }, errors: errors]
    }
}

def parseTokenEach(List tokens) {
    return tokens.collect { tok -> parseToken(tok.toString()) }
}

def tokenRun(String token, String prefix, String rest, Map flair) {
    def rid = prefix + safeId(rest)
    if (!rest || !(rid ==~ ID_PATTERN())) {
        return [token: token, error: "--flair token '${token}' gives an unusable run id '${rid}'; write this run in runs.yaml instead".toString()]
    }
    return [token: token, id: rid, flair: flair]
}

def looksLikeUrl(String s) {
    return s.contains('://') || s.startsWith('git@') || s.startsWith('/') || s.startsWith('file:')
}

def tokenBasename(String path) {
    def p = path.replaceAll('/+$', '')
    def i = p.lastIndexOf('/')
    return i >= 0 ? p.substring(i + 1) : p
}

def containerName(String uri) {
    def i = uri.indexOf('://')
    def base = tokenBasename(i >= 0 ? uri.substring(i + 3) : uri)
    def at = base.indexOf('@')
    base = at >= 0 ? base.substring(0, at) : base
    def colon = base.lastIndexOf(':')
    return colon >= 0 ? base.substring(0, colon) : base
}

/*
 * ---------------------------------------------------------------------------------------------
 * Shared helpers
 * ---------------------------------------------------------------------------------------------
 */

// Final check against assets/schema/*.json. Messages from nf-schema are kept, minus noise.
def schemaCheck(Object obj, String schemaName, String where, List errors) {
    def schema = "${projectDir}/assets/schema/${schemaName}"
    def found = validate(obj, schema, exitOnError: false)
    def useful = found.findAll { m -> !m.contains('False schema always fails') }.unique()
    if (found && !useful) {
        useful = found.take(1)
    }
    useful.each { m -> errors << "${where}: does not match ${schemaName}: ${m}" }
}

// --bundles: comma-separated directories or globs of earlier bundles.
def parseBundlesParam(Object value) {
    if (value == null || value == false || value.toString().trim() == '') {
        return []
    }
    def errors = []
    def out = []
    value.toString().split(',').collect { s -> s.trim() }.findAll { s -> s }.each { item ->
        def isGlob = item.contains('*') || item.contains('?') || item.contains('[') || item.contains('{')
        if (isGlob) {
            def hits = file(item, type: 'any')
            hits = hits instanceof List ? hits : [hits]
            hits = hits.findAll { h -> h.exists() }
            if (!hits) {
                errors << "--bundles pattern '${item}' matches nothing."
            }
            out.addAll(hits)
        }
        else {
            def f = file(item)
            if (!f.exists()) {
                errors << "--bundles path '${item}' does not exist."
            }
            else {
                out << f
            }
        }
    }
    if (errors) {
        failInputs(errors)
    }
    return out.unique { f -> f.toString() }
}
