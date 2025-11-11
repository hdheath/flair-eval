#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =============================================================================
// DATASET AND TEST SET CONFIGURATION MODULE
// =============================================================================
// This module contains reusable class definitions for creating datasets and 
// test configurations for FLAIR
// 
// =============================================================================

// =============================================================================
// CLASS DEFINITIONS
// =============================================================================

class Dataset {
    String name
    def reads  // Can be String (single file) or List (multiple files)
    String genome
    String gtf
    
    // Optional files for FLAIR processes
    String junction_bed
    String junction_tab
    
    // Constructor
    Dataset(String name, Map config) {
        this.name = name
        this.reads = config.reads
        this.genome = config.genome
        this.gtf = config.gtf
        this.junction_bed = config.junction_bed
        this.junction_tab = config.junction_tab
    }
    
    // Get reads as a list (handles both single file and multiple files)
    List<String> getReadsList() {
        return reads instanceof List ? reads : [reads]
    }
    
    // Validation methods
    boolean hasJunctionBed() { return junction_bed != null }
    boolean hasJunctionTab() { return junction_tab != null }
    boolean hasAnyJunctions() { return hasJunctionBed() || hasJunctionTab() }
    
    // Get file objects for Nextflow processes
    def getJunctionBedFile() {
        return junction_bed ? file(junction_bed) : file('NO_JUNCTION_BED')
    }
    
    def getJunctionTabFile() {
        return junction_tab ? file(junction_tab) : file('NO_JUNCTION_TAB')
    }
    
    String description() {
        def junctions = []
        if (hasJunctionBed()) junctions << "bed"
        if (hasJunctionTab()) junctions << "tab"
        def junction_info = junctions.isEmpty() ? "no junctions" : "junctions: ${junctions.join(', ')}"
        
        return "Dataset[${name}, ${junction_info}]"
    }
    
    // Validation method to check if all required files exist
    Map<String, Boolean> validateFiles() {
        def validation = [:]
        validation.reads = reads ? (reads instanceof List ? reads.every { file(it).exists() } : file(reads).exists()) : false
        validation.genome = genome ? file(genome).exists() : false  
        validation.gtf = gtf ? file(gtf).exists() : false
        validation.junction_bed = junction_bed ? file(junction_bed).exists() : true  // true if not specified
        validation.junction_tab = junction_tab ? file(junction_tab).exists() : true  // true if not specified
        
        return validation
    }
    
    // Check if dataset is valid (all specified files exist)
    boolean isValid() {
        def validation = validateFiles()
        return validation.values().every { it == true }
    }
}

class TestSet {
    String name
    Dataset dataset
    List<List> alignOptions
    List<List> transcriptomeOptions
    
    // Constructor
    TestSet(String name, Dataset dataset, Map options) {
        this.name = name
        this.dataset = dataset
        this.alignOptions = options.align
        this.transcriptomeOptions = options.transcriptome ?: []
    }
    
    // Calculate total number of jobs this test set will produce
    int totalJobs() {
        return alignOptions.size() * transcriptomeOptions.size()
    }
    
    // Get detailed job breakdown
    Map getJobBreakdown() {
        return [
            align_modes: alignOptions.size(),
            transcriptome_modes: transcriptomeOptions.size(), 
            total_jobs: totalJobs(),
            dataset: dataset.name
        ]
    }
    
    String description() {
        return "TestSet[${name}, ${totalJobs()} jobs on ${dataset.name}]"
    }
    
    // Validation method
    boolean isValid() {
        if (!dataset.isValid()) return false
        if (alignOptions.isEmpty()) return false
        if (transcriptomeOptions.isEmpty()) return false
        return true
    }
    
    // Get summary of what will be tested
    String getSummary() {
        def align_modes = alignOptions.collect { it[0] }.join(', ')
        def transcriptome_modes = transcriptomeOptions.collect { it[0] }.join(', ')
        
        return """
TestSet: ${name}
Dataset: ${dataset.description()}
Align modes: ${align_modes}
Transcriptome modes: ${transcriptome_modes}
Total jobs: ${totalJobs()}
""".trim()
    }
}

// =============================================================================
// HELPER FUNCTIONS FOR CREATING STANDARD CONFIGURATIONS
// =============================================================================

// Standard mode lists for easy reference
def getStandardAlignModes() {
    return ['default', 'nvrna', 'stringent', 'sensitive']
}

def getStandardTranscriptomeModes() {
    return ['with_gtf', 'denovo', 'stringent', 'permissive']
}

// Helper function to create align options from mode list
def createAlignOptions(List<String> modes) {
    return modes.collect { mode -> 
        [mode, getAlignParametersForMode(mode)] 
    }
}

// Helper function to create transcriptome options from mode list  
def createTranscriptomeOptions(List<String> modes) {
    return modes.collect { mode -> 
        [mode, getTranscriptomeParametersForMode(mode)] 
    }
}

// Convenience functions for single modes (avoid redundancy)
def createAlignOption(String mode) {
    return [mode, getAlignParametersForMode(mode)]
}

def createTranscriptomeOption(String mode) {
    return [mode, getTranscriptomeParametersForMode(mode)]
}

// =============================================================================
// PARAMETER DEFINITIONS FOR EACH STAGE
// =============================================================================

// FLAIR ALIGN parameter presets
def getAlignParameters() {
    return [
        'default': [
            '--output_dir': 'align_output',
            '--threads': '4'
        ],
        'nvrna': [
            '--output_dir': 'align_output', 
            '--threads': '4',
            '--nvrna': ''
        ],
        'stringent': [
            '--output_dir': 'align_output',
            '--threads': '4', 
            '--quality': '13',
            '--end_window': '20'
        ],
        'sensitive': [
            '--output_dir': 'align_output',
            '--threads': '4',
            '--quality': '1', 
            '--end_window': '5'
        ]
    ]
}

// FLAIR TRANSCRIPTOME parameter presets 
def getTranscriptomeParameters() {
    return [
        'with_gtf': [
            '--output_dir': 'transcriptome_output',
            '--threads': '4'
        ],
        'denovo': [
            '--output_dir': 'transcriptome_output',
            '--threads': '4',
            '--no_gtf': ''
        ],
        'stringent': [
            '--output_dir': 'transcriptome_output', 
            '--threads': '4',
            '--stringent': ''
        ],
        'permissive': [
            '--output_dir': 'transcriptome_output',
            '--threads': '4',
            '--no_redundant': 'none',
            '--end_window': '0'
        ]
    ]
}

// Helper function to get parameters for a specific mode
def getAlignParametersForMode(String mode) {
    def params = getAlignParameters()
    if (!params.containsKey(mode)) {
        throw new Exception("Unknown align mode: ${mode}. Available modes: ${params.keySet().join(', ')}")
    }
    return params[mode]
}

def getTranscriptomeParametersForMode(String mode) {
    def params = getTranscriptomeParameters()
    if (!params.containsKey(mode)) {
        throw new Exception("Unknown transcriptome mode: ${mode}. Available modes: ${params.keySet().join(', ')}")
    }
    return params[mode]
}

// Helper function to convert parameter map to command line arguments
def parametersToArgs(Map params) {
    def args = []
    params.each { key, value ->
        if (value == '') {
            // Flag parameter (no value)
            args << key
        } else {
            // Parameter with value
            args << key << value
        }
    }
    return args
}

// =============================================================================
// EXAMPLE DATASET DEFINITIONS
// =============================================================================

// function to create commonly used datasets
def createExampleDatasets() {
    def datasets = [:]
    
    // Human dataset example
    datasets.human = new Dataset('human', [
        reads: '/private/groups/brookslab/hdheath/projects/flair-eval/cache/datasets/human/WTC11.100reads.fasta',
        genome: '/private/groups/brookslab/reference_genomes/gencode_human/GRCh38.primary_assembly.genome.fa',
        gtf: '/private/groups/brookslab/reference_genomes/gencode_human/gencode.v47.annotation.gtf',
        junction_bed: '/private/groups/brookslab/hdheath/projects/flair-eval/cache/datasets/human/junctions.bed',
        junction_tab: null
    ])
    
    // A549 dataset example  
    datasets.a549 = new Dataset('A549_cDNA', [
        reads: '/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/WTC11.100reads.fasta',
        genome: '/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/GRCh38.primary_assembly.genome.fa',
        gtf: '/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/gencode.v48.annotation.gtf',
        junction_bed: null,
        junction_tab: '/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/WTC11_all.SJ.out.tab'
    ])
    
    return datasets
}

// =============================================================================
// EXAMPLE TEST SET CONFIGURATIONS  
// =============================================================================

def createStandardTestSets(Map datasets) {
    def test_sets = [:]
    
    // Quick validation test 
    test_sets.quick = new TestSet('quick_test', datasets.a549, [
        align: [createAlignOption('default')],
        transcriptome: [createTranscriptomeOption('with_gtf')]
    ])
    
    // Align method comparison 
    test_sets.align_comparison = new TestSet('align_comparison', datasets.human, [
        align: createAlignOptions(getStandardAlignModes()),
        transcriptome: [createTranscriptomeOption('with_gtf')]
    ])
    
    // Transcriptome method comparison 
    test_sets.transcriptome_comparison = new TestSet('transcriptome_comparison', datasets.a549, [
        align: [createAlignOption('default')],
        transcriptome: createTranscriptomeOptions(getStandardTranscriptomeModes())
    ])
    
    // Full matrix test - all combinations with parameters
    test_sets.full_matrix = new TestSet('full_matrix', datasets.human, [
        align: createAlignOptions(getStandardAlignModes()),
        transcriptome: createTranscriptomeOptions(getStandardTranscriptomeModes())
    ])
    
    return test_sets
}

// =============================================================================
// UTILITY FUNCTIONS
// =============================================================================

// Print summary of all datasets
def printDatasetSummary(Map datasets) {
    println("\n=== AVAILABLE DATASETS ===")
    datasets.each { name, dataset ->
        println(dataset.description())
        def validation = dataset.validateFiles()
        def missing = validation.findAll { k, v -> !v }.keySet()
        if (!missing.empty) {
            println("  ⚠️  Missing files: ${missing.join(', ')}")
        } else {
            println("  ✅ All files validated")
        }
    }
    println("==========================\n")
}

// Print summary of all test sets
def printTestSetSummary(Map test_sets) {
    println("\n=== AVAILABLE TEST SETS ===")
    test_sets.each { name, test_set ->
        println(test_set.getSummary())
        println("---")
    }
    println("============================\n")
}