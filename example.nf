#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process FlairAlign {

    executor 'slurm'
    clusterOptions '--partition=short'

    time '1h'
    memory '32 GB'

    conda '/private/home/hdheath/miniforge3/envs/flair-dev'

    publishDir "results/align", mode: 'link'

    input:
    tuple val(sample), val(opt_label), val(opt_args), path(reads), path(genome), val(genome_ref)

    output:
    tuple val(sample),
          val(opt_label),
          val(opt_args),
          path("${sample}_${opt_label}.bam"),
          path("${sample}_${opt_label}.bam.bai"),
          path("${sample}_${opt_label}.bed"),
          val(genome_ref)

    script:
    """
    flair align \\
        --reads ${reads} \\
        --genome ${genome} \\
        ${opt_args} \\
        --output ${sample}_${opt_label}
    """
}

process FlairTranscriptome {

    executor 'slurm'
    clusterOptions '--partition=short'

    time '1h'
    memory '32 GB'

    conda '/private/home/hdheath/miniforge3/envs/flair-dev'

    publishDir "results/transcriptome", mode: 'link'

    input:
    tuple val(sample), val(align_label), val(trans_label), val(trans_args),
          path(align_bam), path(align_bai), path(genome)

    output:
    path "${sample}_${align_label}_${trans_label}.bed", optional: true
    path "${sample}_${align_label}_${trans_label}.gtf", optional: true
    path "${sample}_${align_label}_${trans_label}.fa", optional: true
    path "${sample}_${align_label}_${trans_label}.txt", optional: true

    script:
    """

    flair transcriptome \\
        -b ${align_bam} \\
        --genome ${genome} \\
        ${trans_args} \\
        --output ${sample}_${align_label}_${trans_label}
    """
}

workflow {
    def align_files_ch = Channel.of(
        tuple(
            'Test1',
            file('/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/WTC11.100reads.fasta'),
            file('/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/GRCh38.primary_assembly.genome.fa')
        ),
        tuple(
            'Test2',
            file('/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/WTC11.10reads.fasta'),
            file('/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/GRCh38.primary_assembly.genome.fa')
        )
    )

    
    def align_options_ch = Channel.of(
        // describing name, command line args
        tuple('default', ''),
        tuple('with_nvrna', '--nvrna')
    )

    def align_inputs = align_files_ch
        .combine(align_options_ch)
        .map { datasetName, readsPath, genomePath, optLabel, optArgs ->
            tuple(datasetName, optLabel, optArgs, readsPath, genomePath, genomePath)
        }

    def aligned_outputs = FlairAlign(align_inputs)
        .map { emit_tuple ->
            def (sample, alignLabel, alignArgs, bamPath, baiPath, bedPath, genomeRef) = emit_tuple
            def genomeFile = file(genomeRef).toAbsolutePath()
            println "[DEBUG] Align emit -> sample=${sample}, align_label=${alignLabel}, bam=${bamPath}, genome=${genomeFile}"
            tuple(sample, alignLabel, alignArgs, bamPath, baiPath, bedPath, genomeFile)
        }

    def transcriptome_options_ch = Channel.of(
        // describing name, command line args
        tuple('default', '--gtf /private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/gencode.v48.annotation.gtf')
    )

    def transcriptome_inputs = aligned_outputs
        .combine(transcriptome_options_ch)
        .map { sample, alignLabel, alignArgs, bamPath, baiPath, bedPath, genomePath, transLabel, transArgs ->
            println "[DEBUG] Transcriptome tuple -> sample=${sample}, align_label=${alignLabel}, trans_label=${transLabel}"
            println "[DEBUG]   BAM: ${bamPath}"
            println "[DEBUG]   BAI: ${baiPath}"
            println "[DEBUG]   Genome: ${genomePath}"
            println "[DEBUG]   Align args: ${alignArgs}"
            println "[DEBUG]   Trans args: ${transArgs}"
            tuple(sample, alignLabel, transLabel, transArgs, bamPath, baiPath, genomePath)
        }


    FlairTranscriptome(transcriptome_inputs)
}
