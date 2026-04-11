// =============================================================================
// Module: Badread Simulation
// =============================================================================
// Generates synthetic long reads from a genome + GTF annotation using Badread.
//
// Three processes:
//   1. GenerateTranscriptomeFasta — extracts spliced cDNA sequences from genome+GTF
//   2. BadreadSimulate — simulates reads with controlled error profiles
//   3. AlignSimulatedReads — splice-aligns simulated reads back to genome
//
// The resulting BAM feeds into the standard FlairPartition → FlairTranscriptome
// → Evaluation pipeline, identical to real data.
// =============================================================================


/*
 * GenerateTranscriptomeFasta
 * --------------------------
 * Extract spliced cDNA sequences from genome+GTF for use as Badread reference.
 * Each transcript gets a depth header controlling relative abundance.
 */
process GenerateTranscriptomeFasta {
    tag "gen_txome_${dataset_name}"
    publishDir "${params.outdir}/simulation/${test_name}", mode: 'symlink'

    input:
    tuple val(test_name), val(dataset_name),
          path(genome_fa), path(genome_fai), path(annotation_gtf),
          val(depth), val(region)

    output:
    tuple val(test_name), val(dataset_name),
          path("transcriptome.fa"),
          path(genome_fa), path(genome_fai), path(annotation_gtf),
          emit: transcriptome_ref

    script:
    def region_arg = region ? "--region ${region}" : ""
    """
    python3 ${projectDir}/bin/transcriptome_testing/generate_transcriptome_fa.py \\
        --genome ${genome_fa} \\
        --gtf ${annotation_gtf} \\
        --output transcriptome.fa \\
        --depth ${depth} \\
        ${region_arg}
    """
}


/*
 * BadreadSimulate
 * ---------------
 * Generate synthetic long reads using Badread.  Each invocation uses a
 * specific "scenario" — a named combination of error model, identity,
 * chimera rate, glitch params, and fragment length distribution.
 */
process BadreadSimulate {
    tag "${dataset_name}_${scenario_name}"
    publishDir "${params.outdir}/simulation/${test_name}", mode: 'symlink'

    input:
    tuple val(test_name), val(dataset_name), val(scenario_name),
          path(transcriptome_fa),
          val(quantity),
          val(error_model),
          val(identity),
          val(length_params),
          val(chimeras),
          val(junk_reads),
          val(random_reads),
          val(glitches),
          val(seed)

    output:
    tuple val(test_name), val(dataset_name), val(scenario_name),
          path("${scenario_name}_reads.fastq.gz"),
          emit: simulated_reads

    script:
    """
    badread simulate \\
        --reference ${transcriptome_fa} \\
        --quantity ${quantity} \\
        --error_model ${error_model} \\
        --identity ${identity} \\
        --length ${length_params} \\
        --chimeras ${chimeras} \\
        --junk_reads ${junk_reads} \\
        --random_reads ${random_reads} \\
        --glitches ${glitches} \\
        --seed ${seed} \\
        | gzip > ${scenario_name}_reads.fastq.gz
    """
}


/*
 * SubsetGenomeForSimAlignment
 * ---------------------------
 * Extract the target region(s) from the full genome FASTA so that simulated
 * reads can be aligned to just a tiny reference (~200 KB vs ~3 GB).
 * This makes minimap2 indexing + alignment orders of magnitude faster.
 *
 * When badread_region is not set, this process is skipped and the full
 * genome is used (see main.nf wiring).
 */
process SubsetGenomeForSimAlignment {
    tag "subset_genome_${dataset_name}"
    publishDir "${params.outdir}/simulation/${test_name}", mode: 'symlink'

    input:
    tuple val(test_name), val(dataset_name),
          path(genome_fa), path(genome_fai),
          val(region)

    output:
    tuple val(test_name), val(dataset_name),
          path("sim_subset_genome.fa"),
          path("sim_subset_genome.fa.fai"),
          emit: subset_genome

    script:
    // region can be "chr12:6434516-6638373" or just "chr12"
    """
    samtools faidx ${genome_fa} ${region} > sim_subset_genome.fa
    samtools faidx sim_subset_genome.fa
    """
}


/*
 * AlignSimulatedReads
 * -------------------
 * Splice-align Badread output to the GENOME (not transcriptome) using
 * minimap2.  This mirrors real FLAIR usage where reads are genome-aligned.
 *
 * When a subset genome is available (via SubsetGenomeForSimAlignment),
 * alignment is dramatically faster since minimap2 only indexes the
 * target region.
 */
process AlignSimulatedReads {
    tag "${dataset_name}_${scenario_name}"
    publishDir "${params.outdir}/simulation/${test_name}", mode: 'symlink'

    input:
    tuple val(test_name), val(dataset_name), val(scenario_name),
          path(reads_fastq),
          path(genome_fa), path(genome_fai)

    output:
    tuple val(test_name), val(dataset_name), val(scenario_name),
          path("${scenario_name}.bam"),
          path("${scenario_name}.bam.bai"),
          emit: aligned

    script:
    """
    minimap2 -ax splice --secondary=no -t ${task.cpus} \\
        ${genome_fa} ${reads_fastq} \\
        | samtools sort -@ ${task.cpus} -o ${scenario_name}.bam

    samtools index ${scenario_name}.bam
    """
}
