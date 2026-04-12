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
                    stringtie2: sanitizeModeMap(stringtie2_modes)
                ])))
            }
        }
        return test_sets_list
    }
}
