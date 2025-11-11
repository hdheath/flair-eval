#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =============================================================================
// MINIMAL FLAIR WORKFLOW: HARDCODED EXAMPLE RUN
// =============================================================================
// This is the simplest possible version with hardcoded paths and commands
// Just runs: align → partition → transcriptome on a single dataset
// =============================================================================

process FlairAlign {
    conda 'flair-dev'
    publishDir "results/align", mode: 'symlink'
    
    output:
    tuple path("example.bam"), path("example.bed"), emit: alignments
    
    script:
    """
    flair align \
        -r /private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/WTC11.100reads.fasta \
        -g /private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/GRCh38.primary_assembly.genome.fa \
        -o example \
        --threads 4
    """
}

process FlairPartition {
    conda 'flair-dev'
    publishDir "results/partition", mode: 'symlink'
    
    input:
    tuple path(bam), path(bed)
    
    output:
    tuple path("example_chr1.bam"), path("example_chr1.bed"), emit: partitioned
    
    script:
    """
    python ${projectDir}/bin/simple_partition.py \
        --bam ${bam} \
        --bed ${bed} \
        --region chr1:1-100000 \
        --output-prefix example_chr1 \
        --gtf /private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/gencode.v48.annotation.gtf \
        --genome /private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/GRCh38.primary_assembly.genome.fa
    """
}

process FlairTranscriptome {
    conda 'flair-dev'
    publishDir "results/transcriptome", mode: 'symlink'
    
    input:
    tuple path(bam), path(bed)
    
    output:
    path("example.transcriptome.fa"), emit: transcriptome
    
    script:
    """
    flair transcriptome \
        -b ${bam} \
        -o example.transcriptome \
        --threads 4
    """
}

workflow {
    // Run the three processes in sequence
    FlairAlign()
    FlairPartition(FlairAlign.out.alignments)
    FlairTranscriptome(FlairAlign.out.alignments)
    
    // Print completion message
    workflow.onComplete {
        println """
        
        ====================================================
        FLAIR WORKFLOW COMPLETE!
        ====================================================
        
        What happened:
        1. ✅ Aligned reads to genome (chromosome 1, first 100k bases)
        2. ✅ Partitioned alignments by region  
        3. ✅ Generated transcriptome assembly
        
        Results:
        📁 results/align/example.{bam,bed}
        📁 results/partition/example_chr1.{bam,bed}
        📁 results/transcriptome/example.transcriptome.fa
        
        ====================================================
        """
    }
}

/*
====================================================
USAGE:
====================================================

Just run:
    nextflow run minimal_flair_example.nf

Or with resume:
    nextflow run minimal_flair_example.nf -resume

This runs a single example with:
- Input: WTC11.100reads.fasta (100 reads)
- Region: chr1:1-100000 (first 100k bases of chromosome 1)
- Output: Alignment, partition, and transcriptome files

Perfect for testing FLAIR installation or workflow setup!
====================================================
*/