#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =============================================================================
// MINIMAL FLAIR WORKFLOW: HARDCODED EXAMPLE RUN
// =============================================================================
// This is the simplest possible version with hardcoded paths and commands
// Just runs: align → partition → transcriptome on a single dataset
// =============================================================================

process FlairAlign {
    conda '/private/home/hdheath/miniforge3/envs/flair-dev'
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
    conda '/private/home/hdheath/miniforge3/envs/flair-dev'
    publishDir "results/partition", mode: 'symlink'
    
    input:
    tuple path(bam), path(bed)
    
    output:
    tuple path("example_chr1.bam"), path("example_chr1.bam.bai"), path("example_chr1.bed"), path("example_chr1_annotation.gtf"), emit: partitioned
    
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
    conda '/private/home/hdheath/miniforge3/envs/flair-dev'
    publishDir "results/transcriptome", mode: 'symlink'
    
    input:
    tuple path(bam), path(bai), path(bed), path(gtf)
    
    output:
    path("example_chr1.transcriptome.isoforms.fa"), emit: transcriptome
    
    script:
    """
    flair transcriptome \
        -b ${bam} \
        -g /private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/GRCh38.primary_assembly.genome.fa \
        -f ${gtf} \
        -o example_chr1.transcriptome \
        --threads 4
    """
}

workflow {
    // Run the three processes in sequence
    FlairAlign()
    FlairPartition(FlairAlign.out.alignments)
    FlairTranscriptome(FlairPartition.out.partitioned)
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

====================================================
*/