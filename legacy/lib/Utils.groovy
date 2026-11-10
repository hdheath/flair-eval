import groovy.json.JsonSlurper

/**
 * Utils - Shared utility functions for the FLAIR evaluation pipeline.
 *
 * Provides mode name sanitization and samplesheet parsing.
 * Nextflow auto-loads all .groovy files from lib/, making these
 * available in the main workflow and all modules.
 */
class Utils {

    /**
     * Canonical column order for the STAGE-1 assembly manifest CSV.
     * The ASSEMBLE subworkflow writes rows in EXACTLY this order; the
     * EVALUATION entrypoint parses by these header names. Keep in sync.
     */
    static final List<String> MANIFEST_COLUMNS = [
        'test_name', 'dataset_name', 'align_mode', 'partition_mode',
        'transcriptome_mode', 'tool',
        'isoforms_bed', 'isoforms_gtf', 'isoform_read_map', 'ted_log', 'isoform_counts',
        'firstpass', 'cds_bed', 'cds_info',
        'bam', 'bai', 'reads_bed', 'genome', 'gtf', 'cage_peaks', 'drna_peaks',
        'cage_signal_plus', 'cage_signal_minus', 'drna_signal_plus', 'drna_signal_minus',
        'library_type', 'partition_region', 'partition_args'
    ]

    /** Header line (no trailing newline) for the manifest CSV. */
    static String manifestHeader() {
        return MANIFEST_COLUMNS.join(',')
    }

    /**
     * Parse the STAGE-1 assembly manifest CSV into a list of row maps.
     *
     * Mirrors the samplesheet-reading pattern: reads the header, splits each
     * data line, and returns [colName: value] maps. Handles quoted fields
     * (values that contained a comma were wrapped in double quotes by the
     * ASSEMBLE writer). Skips blank lines and a leading '#'-comment line.
     *
     * @param manifestPath  Path/String to the manifest CSV
     * @return              List<Map> one map per assembly row
     */
    static List<Map> parseManifest(manifestPath) {
        def rows = []
        def header = null
        new File(manifestPath.toString()).withReader { reader ->
            reader.eachLine { rawLine ->
                def line = rawLine
                if (line == null) return
                if (!line.trim()) return
                if (line.startsWith('#')) return
                def values = splitCsvLine(line)
                if (header == null) {
                    header = values
                    return
                }
                // Pad/trim to header width
                if (values.size() > header.size()) {
                    values = values[0..<header.size()]
                } else if (values.size() < header.size()) {
                    values = values + ([''] * (header.size() - values.size()))
                }
                rows.add([header, values].transpose().collectEntries())
            }
        }
        return rows
    }

    /**
     * Split a single CSV line honoring double-quoted fields (no embedded
     * newlines). Quotes are stripped from the returned values.
     */
    static List<String> splitCsvLine(String line) {
        def out = []
        def cur = new StringBuilder()
        def inQuotes = false
        for (int i = 0; i < line.length(); i++) {
            def c = line.charAt(i)
            if (c == ('"' as char)) {
                inQuotes = !inQuotes
            } else if (c == (',' as char) && !inQuotes) {
                out.add(cur.toString())
                cur = new StringBuilder()
            } else {
                cur.append(c)
            }
        }
        out.add(cur.toString())
        return out
    }

    /**
     * Sanitize mode names to use hyphens instead of underscores.
     * This ensures underscores can be used as delimiters in filenames.
     */
    static String sanitizeModeName(String name) {
        return name.replaceAll('_', '-').replaceAll(' ', '-')
    }

    /**
     * Sanitize all keys in a mode map.
     */
    static Map sanitizeModeMap(Map modes) {
        if (modes == null || modes.isEmpty()) {
            return modes
        }
        def result = [:]
        modes.each { k, v ->
            result[sanitizeModeName(k)] = v
        }
        return result
    }

    /**
     * Parse a samplesheet CSV and pipeline config JSON into a list of TestSet objects.
     *
     * @param csvPath       Path to samplesheet CSV
     * @param paramsFile    Path to JSON params file
     * @param testName      Base test name (prefixed to each dataset name)
     * @return              List of TestSet objects
     */
    static List<TestSet> parseSamplesheet(csvPath, paramsFile, String testName) {
        def jsonSlurper = new JsonSlurper()
        def pipeline_config = jsonSlurper.parse(new File(paramsFile.toString()))

        // Extract mode maps from JSON config (with defaults if not specified)
        def align_modes = pipeline_config.align ?: [default: '']
        def partition_modes = pipeline_config.partition ?: [all: '--all']
        def transcriptome_modes = pipeline_config.transcriptome ?: [default: '']
        def bambu_modes = pipeline_config.bambu ?: [:]
        def isoquant_modes = pipeline_config.isoquant ?: [:]
        def isoseq_modes = pipeline_config.isoseq ?: [:]
        def flames_modes = pipeline_config.flames ?: [:]
        def stringtie2_modes = pipeline_config.stringtie2 ?: [:]
        def lraa_modes = pipeline_config.lraa ?: [:]

        def test_sets_list = []
        new File(csvPath.toString()).withReader { reader ->
            def header = reader.readLine().split(',')
            reader.eachLine { line ->
                // Skip empty lines
                if (!line.trim()) return

                // Split and trim values, handling trailing commas
                def values = line.split(',', -1)  // -1 keeps trailing empty strings

                // Ensure values array matches header length
                if (values.size() > header.size()) {
                    values = values[0..<header.size()]
                } else if (values.size() < header.size()) {
                    values = values + ([''] * (header.size() - values.size()))
                }

                def row = [header, values].transpose().collectEntries()

                // Create Dataset object from CSV row
                def dataset = new Dataset(row.sample_id, ([
                    genome: row.genome,
                    gtf: row.gtf,
                    library_type: row.library_type && row.library_type != '' ? row.library_type : null,
                    bam: row.bam && row.bam != '' ? row.bam : null,
                    reads: row.reads && row.reads != '' ? row.reads : null,
                    cage: row.cage && row.cage != '' ? row.cage : null,
                    drna: row.drna && row.drna != '' ? row.drna : null,
                    junction_tab: (
                        row.junction_tab && row.junction_tab != '' ? row.junction_tab :
                        (row.junctions && row.junctions != '' ? row.junctions : null)
                    ),
                    cage_signal_plus: row.cage_signal_plus && row.cage_signal_plus != '' ? row.cage_signal_plus : null,
                    cage_signal_minus: row.cage_signal_minus && row.cage_signal_minus != '' ? row.cage_signal_minus : null,
                    drna_signal_plus: row.drna_signal_plus && row.drna_signal_plus != '' ? row.drna_signal_plus : null,
                    drna_signal_minus: row.drna_signal_minus && row.drna_signal_minus != '' ? row.drna_signal_minus : null
                ]))

                test_sets_list.add(new TestSet("${testName}_${dataset.name}", dataset, ([
                    align: sanitizeModeMap(align_modes),
                    partition: sanitizeModeMap(partition_modes),
                    transcriptome: sanitizeModeMap(transcriptome_modes),
                    bambu: sanitizeModeMap(bambu_modes),
                    isoquant: sanitizeModeMap(isoquant_modes),
                    isoseq: sanitizeModeMap(isoseq_modes),
                    flames: sanitizeModeMap(flames_modes),
                    stringtie2: sanitizeModeMap(stringtie2_modes),
                    lraa: sanitizeModeMap(lraa_modes)
                ])))
            }
        }
        return test_sets_list
    }
}
