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
        --reads ${reads instanceof List ? reads.join(' ') : reads} \\
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

process FlairCorrect {

    executor 'slurm'
    clusterOptions '--partition=short'
    time '1h'
    memory '32 GB'
    conda '/private/home/hdheath/miniforge3/envs/flair-dev'
    publishDir "results/correct", mode: 'link'

    input:
    tuple val(sample), val(align_label),
          val(cor_label), val(cor_args),
          path(align_bed), path(genome), val(genome_ref)

    output:
    tuple val(sample),
          val(align_label),
          val(cor_label),
          val(cor_args),
          path("${sample}_${align_label}_${cor_label}_all_corrected.bed"),
          path("${sample}_${align_label}_${cor_label}_all_inconsistent.bed", optional: true),
          path("${sample}_${align_label}_${cor_label}_cannot_verify.bed", optional: true),
          val(genome_ref)

    script:
    """
    flair correct \\
        -q ${align_bed} \\
        ${cor_args} \\
        --output ${sample}_${align_label}_${cor_label}
    
    # Create empty files for optional outputs if they don't exist
    touch ${sample}_${align_label}_${cor_label}_all_inconsistent.bed
    touch ${sample}_${align_label}_${cor_label}_cannot_verify.bed
    """
}

process FlairCollapse {

    executor 'slurm'
    clusterOptions '--partition=short'
    time '1h'
    memory '32 GB'
    conda '/private/home/hdheath/miniforge3/envs/flair-dev'
    publishDir "results/collapse", mode: 'link'

    input:
    tuple val(sample), val(align_label),
          val(cor_label), val(cor_args),
          path(cor_bed), path(genome), val(genome_ref),
          path(reads), val(col_label), val(col_args)

    output:
    tuple val(sample),
          val(align_label),
          val(cor_label),
          val(col_label),
          path("${sample}_${align_label}_${cor_label}_${col_label}.bed"),
          path("${sample}_${align_label}_${cor_label}_${col_label}.gtf", optional: true),
          path("${sample}_${align_label}_${cor_label}_${col_label}.fa", optional: true),
          val(genome_ref)

    script:
    """
    flair collapse \\
        -g ${genome} \\
        -q ${cor_bed} \\
        -r ${reads instanceof List ? reads.join(' ') : reads} \\
        ${col_args} \\
        --output ${sample}_${align_label}_${cor_label}_${col_label}
    
    # Create empty files for outputs if they don't exist
    touch ${sample}_${align_label}_${cor_label}_${col_label}.bed
    touch ${sample}_${align_label}_${cor_label}_${col_label}.gtf
    touch ${sample}_${align_label}_${cor_label}_${col_label}.fa
    """
}

workflow {

    ///// Step 1: Create all input combinations /////
    
    def base_samples = Channel.of(
        tuple(
            'Test1',
            [
                file('/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/WTC11.100reads.fasta'),
                file('/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/WTC11.10reads.fasta')
            ],
            file('/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/GRCh38.primary_assembly.genome.fa')
        ),
        tuple(
            'Test2',
            file('/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/WTC11.10reads.fasta'),
            file('/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/GRCh38.primary_assembly.genome.fa')
        )
    )

    def align_options = Channel.of(
        tuple('default', ''),
        tuple('with_nvrna', '--nvrna')
    )

    ///// Step 2: Align /////
    
    def align_inputs = base_samples
        .combine(align_options)
        .map { sample, readsList, genome, label, args ->
            def readsFinal = readsList instanceof List ? readsList : [readsList]
            tuple(sample, label, args, readsFinal, genome, genome)
        }

    def aligned_outputs = FlairAlign(align_inputs)

    ///// Step 3: Transcriptome (using multiMap to split aligned_outputs) /////
    
    // Use multiMap to split the aligned outputs for different purposes
    def aligned_split = aligned_outputs.multiMap { sample, align_label, align_args, bam, bai, bed, genome_ref ->
        correct_path: tuple(sample, align_label, align_args, bam, bai, bed, genome_ref)
        transcriptome_path: tuple(sample, align_label, align_args, bam, bai, bed, genome_ref)
    }
    
    def transcriptome_options = Channel.of(
        tuple('default', '--gtf /private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/gencode.v48.annotation.gtf')
    )

    def transcriptome_inputs = aligned_split.transcriptome_path
        .combine(transcriptome_options)
        .map { sample, alignLabel, alignArgs, bamPath, baiPath, bedPath, genomePath, transLabel, transArgs ->
            tuple(sample, alignLabel, transLabel, transArgs, bamPath, baiPath, genomePath)
        }

    def transcriptome_outputs = FlairTranscriptome(transcriptome_inputs)

    ///// Step 4: Correct /////
    
    def correct_options = Channel.of(
        tuple('with_gtf', '--gtf /private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/gencode.v48.annotation.gtf'),
        tuple('with_gtf_and_nvrna', '--gtf /private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/gencode.v48.annotation.gtf --nvrna')
    )

    def correct_inputs = aligned_split.correct_path
        .combine(correct_options)
        .map { sample, align_label, align_args, bam, bai, bed, genome, cor_label, cor_args ->
            tuple(sample, align_label, cor_label, cor_args, bed, genome, genome)
        }

    def correct_outputs = FlairCorrect(correct_inputs)

    ///// Step 5: Collapse with reads mapping /////
    
    def collapse_options = Channel.of(
        tuple('with_gtf', '--gtf /private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/gencode.v48.annotation.gtf'),
        tuple('default', '')
    )

    def collapse_inputs = correct_outputs
        .combine(collapse_options)
        .map { sample, align_label, cor_label, cor_args, cor_bed, inconsistent_bed, cannot_verify_bed, genome_ref, col_label, col_args ->
            // Map reads directly based on sample name - this is the key solution
            def reads_files
            if (sample == 'Test1') {
                reads_files = [
                    file('/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/WTC11.100reads.fasta'),
                    file('/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/WTC11.10reads.fasta')
                ]
            } else { // Test2
                reads_files = [file('/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/WTC11.10reads.fasta')]
            }
            
            tuple(sample, align_label, cor_label, cor_args,
                  cor_bed, genome_ref, genome_ref, reads_files,
                  col_label, col_args)
        }

    def collapse_outputs = FlairCollapse(collapse_inputs)
    
    // Print final success
    collapse_outputs.view { "PIPELINE_SUCCESS: $it" }
}