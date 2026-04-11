// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test process
include { PlotIsoforms } from '/private/groups/brookslab/hdheath/projects/flair-eval/tests/modules/visualization/plot_isoforms/../../../../modules/visualization/plot_isoforms/main.nf'

workflow {

    // define custom rules for JSON that will be generated.
    def jsonOutput = createJsonOutput()
    def jsonWorkflowOutput = createJsonWorkflowOutput()

    def input = []

    // run dependencies
    

    // process mapping
    input = []
    
                def no_quantseq = file("NO_QUANTSEQ")
                no_quantseq.text = ''

                input[0] = Channel.of(
                    tuple(
                        'test_set',
                        'sample1',
                        'pre-aligned',
                        'chr-test',
                        'default',
                        file("/private/groups/brookslab/hdheath/projects/flair-eval/tests/data/tiny_reads.bam"),
                        file("/private/groups/brookslab/hdheath/projects/flair-eval/tests/data/tiny_reads.bam.bai"),
                        'chr_test:1000-40000',
                        file("/private/groups/brookslab/hdheath/projects/flair-eval/tests/data/tiny_cage_peaks.bed"),
                        no_quantseq,
                        '',
                        '',
                        '',
                        '',
                        file("/private/groups/brookslab/hdheath/projects/flair-eval/tests/data/tiny_isoforms.bed"),
                        file("/private/groups/brookslab/hdheath/projects/flair-eval/tests/data/tiny_read_map.txt")
                    )
                )
                
    //----

    //run process
    PlotIsoforms.run(input.toArray())

    if (PlotIsoforms.output){

        // consumes all named output channels and stores items in a json file
        PlotIsoforms.out.getNames().each { name ->
            serializeChannel(name, PlotIsoforms.out.getProperty(name), jsonOutput, params.nf_test_output)
        }	  

        // consumes all unnamed output channels and stores items in a json file
        def array = PlotIsoforms.out as List<Object>
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