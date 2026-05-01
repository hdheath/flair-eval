// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test workflow
include { SUMMARY_AND_VIZ } from '/private/home/hdheath/projects/flair4_ted/tests/subworkflows/../../subworkflows/summary_and_viz.nf'

workflow {

    // define custom rules for JSON that will be generated.
    def jsonOutput = createJsonOutput()
    def jsonWorkflowOutput = createJsonWorkflowOutput()

    def input = []

    // run dependencies
    

    // workflow mapping
    input = []
    
                // -- evaluation_results: tuple [test_name, dataset_name, align_mode, partition_mode, transcriptome_mode, eval_tsv]
                input[0] = Channel.of(
                    tuple('nf_test', 'sample1', 'pre-aligned', 'chr-test', 'default',
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_evaluation.tsv")),
                    tuple('nf_test', 'sample1', 'pre-aligned', 'chr-test', 'bambu_default',
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_evaluation_dup.tsv"))
                )

                // -- isoform_categories: tuple [test_name, dataset_name, align_mode, partition_mode, transcriptome_mode, categories.tsv]
                input[1] = Channel.of(
                    tuple('nf_test', 'sample1', 'pre-aligned', 'chr-test', 'default',
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_isoform_categories_default.tsv")),
                    tuple('nf_test', 'sample1', 'pre-aligned', 'chr-test', 'bambu_default',
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_isoform_categories_bambu.tsv"))
                )

                // -- cage_peak_reason_tsvs: tuple [test_name, path]
                input[2] = Channel.of(
                    tuple('nf_test', file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_cage_peak_reasons.tsv"))
                )

                // -- drna_peak_reason_tsvs: tuple [test_name, path]
                input[3] = Channel.of(
                    tuple('nf_test', file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_drna_peak_reasons.tsv"))
                )

                // -- all_eval_inputs: 23-field channel
                def no_bed   = file("NO_ISOFORMS_BED")
                no_bed.text  = ''
                def no_gtf   = file("NO_ISOFORMS_GTF")
                no_gtf.text  = ''
                def no_ted_log = file("NO_TED_LOG")
                no_ted_log.text = ''

                input[4] = Channel.of(
                    tuple('nf_test', 'sample1', 'pre-aligned', 'chr-test', 'default',
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_isoforms.bed"),
                          no_gtf,
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_read_map.txt"),
                          no_ted_log,
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_reads.bam"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_reads.bam.bai"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_reads.bed"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_genome.fa"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_annotation.gtf"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_cage_peaks.bed"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_drna_peaks.bed"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_ref_tss.bed"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_ref_tts.bed"),
                          'pacbio_cDNA',
                          '', '', '', ''),
                    tuple('nf_test', 'sample1', 'pre-aligned', 'chr-test', 'bambu_default',
                          no_bed,
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_isoforms.gtf"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_read_map_tool2.txt"),
                          no_ted_log,
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_reads.bam"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_reads.bam.bai"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_reads.bed"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_genome.fa"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_annotation.gtf"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_cage_peaks.bed"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_drna_peaks.bed"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_ref_tss.bed"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_ref_tts.bed"),
                          'pacbio_cDNA',
                          '', '', '', '')
                )

                // -- dataset_signal_ch
                input[5] = Channel.of(
                    tuple('nf_test', 'pacbio_cDNA', '', '', '', '')
                )

                // -- ted_precision_metrics
                input[6] = Channel.of(
                    tuple('nf_test', 'sample1', 'default',
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_precision_recall_default.tsv"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_per_junction_chain.tsv")),
                    tuple('nf_test', 'sample1', 'bambu_default',
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_precision_recall_bambu.tsv"),
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_per_junction_chain.tsv"))
                )

                // -- ted_gtf_precision
                input[7] = Channel.of(
                    tuple('nf_test', 'sample1', 'default',
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_gtf_precision_recall_default.tsv")),
                    tuple('nf_test', 'sample1', 'bambu_default',
                          file("/private/home/hdheath/projects/flair4_ted/tests/data/tiny_gtf_precision_recall_bambu.tsv"))
                )

                // -- flair_firstpass
                input[8] = Channel.empty()
                
    //----

    //run workflow
    SUMMARY_AND_VIZ.run(input.toArray())
    
    if (SUMMARY_AND_VIZ.output){

        // consumes all named output channels and stores items in a json file
        SUMMARY_AND_VIZ.out.getNames().each { name ->
            serializeChannel(name, SUMMARY_AND_VIZ.out.getProperty(name), jsonOutput, params.nf_test_output)
        }	  
    
        // consumes all unnamed output channels and stores items in a json file
        def array = SUMMARY_AND_VIZ.out as List<Object>
        def i = 0
        array.each { output ->
            serializeChannel(i, output, jsonOutput, params.nf_test_output)
            i += 1
        }    	

    }

    workflow.onComplete = {

        def result = [
            success: workflow.success,
            exitStatus: workflow.exitStatus,
            errorMessage: workflow.errorMessage,
            errorReport: workflow.errorReport
        ]
        new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
        
    }
}


def serializeChannel(name, channel, jsonOutput, outputDir) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: { entry ->
            list.add(entry)
        },
        onComplete: {
            def map = new HashMap()
            map[_name] = list
            def filename = "${outputDir}/output_${_name}.json"
            new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}

def createJsonOutput(_input = null) {
    // _input is needed because a closure is provided to all functions called in the process
    return [
        toJson: { obj ->
            def converted = convertPathsToStrings(obj)
            return groovy.json.JsonOutput.toJson(converted)
        }
    ]
}

def convertPathsToStrings(obj) {
    if (obj instanceof java.nio.file.Path) {
        return obj.toAbsolutePath().toString()
    } else if (obj instanceof Map) {
        return obj.collectEntries { k, v -> [k, convertPathsToStrings(v)] }
    } else if (obj instanceof Collection) {
        return obj.collect { it -> convertPathsToStrings(it) }
    } else {
        return obj
    }
}

def createJsonWorkflowOutput(_input = null) {
    // _input is needed because a closure is provided to all functions called in the workflow
    return [
        toJson: { obj ->
            def filtered = removeNullValues(obj)
            return groovy.json.JsonOutput.toJson(filtered)
        }
    ]
}

def removeNullValues(obj) {
    if (obj instanceof Map) {
        return obj.findAll { _k, v -> v != null }.collectEntries { k, v -> [k, removeNullValues(v)] }
    } else if (obj instanceof Collection) {
        return obj.findAll { it -> it != null }.collect { it -> removeNullValues(it) }
    } else {
        return obj
    }
}