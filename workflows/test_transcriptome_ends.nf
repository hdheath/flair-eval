#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * test_transcriptome_ends.nf
 * ==========================
 * Standalone Nextflow workflow for testing FLAIR's transcriptome end-detection
 * scoring logic.  Runnable independently from the main flair-eval pipeline.
 *
 * Two tracks:
 *   1. Simulated reads   — user-provided FASTQs with injected error profiles
 *   2. Real data sub-region — BAM subset from a targeted genomic locus
 *
 * Usage:
 *   # Track 1: simulated reads
 *   nextflow run workflows/test_transcriptome_ends.nf \
 *       --track simulated \
 *       --simulated_fastq reads.fq \
 *       --simulated_truth truth.bed \
 *       --error_profile "3prime_a_rich" \
 *       --genome GRCh38.fa \
 *       --gtf gencode.v48.gtf
 *
 *   # Track 2: real data sub-region
 *   nextflow run workflows/test_transcriptome_ends.nf \
 *       --track real_subregion \
 *       --bam sample.bam \
 *       --region "chr22:20000000-21000000" \
 *       --genome GRCh38.fa \
 *       --gtf gencode.v48.gtf \
 *       --library_type ont_cDNA
 *
 *   # Both tracks (requires all inputs)
 *   nextflow run workflows/test_transcriptome_ends.nf \
 *       --track both ...
 */

// -------------------------------------------------------------------------
// Module includes
// -------------------------------------------------------------------------
include {
    SubsetReference;
    IngestSimulatedReads;
    SubsetRealBam;
    FlairTranscriptomeTest;
    AssertBoundaries;
    AssertEndScoring;
} from '../modules/transcriptome_testing/main'


// -------------------------------------------------------------------------
// Default parameters
// -------------------------------------------------------------------------
params.track          = 'both'           // 'simulated', 'real_subregion', or 'both'
params.genome         = null
params.gtf            = null

// Simulated reads track
params.simulated_fastq = null
params.simulated_truth = null
params.error_profile   = 'default'

// Real data sub-region track
params.bam            = null
params.bai            = null
params.region         = 'chr22:20000000-21000000'
params.sample_id      = 'real_subregion'
params.library_type   = 'default'

// FLAIR transcriptome extra args
params.flair_extra_args = ''

// Output
params.outdir         = "${launchDir}/results/transcriptome_testing"


// -------------------------------------------------------------------------
// Validation
// -------------------------------------------------------------------------
def validateParams() {
    if (!params.genome) error "ERROR: --genome is required"
    if (!params.gtf)    error "ERROR: --gtf is required"

    if (params.track in ['simulated', 'both']) {
        if (!params.simulated_fastq) error "ERROR: --simulated_fastq required for simulated track"
        if (!params.simulated_truth) error "ERROR: --simulated_truth required for simulated track"
    }
    if (params.track in ['real_subregion', 'both']) {
        if (!params.bam) error "ERROR: --bam required for real_subregion track"
    }
}


// =========================================================================
// WORKFLOW
// =========================================================================
workflow {
    validateParams()

    genome_fa  = file(params.genome)
    genome_fai = file("${params.genome}.fai")
    gtf        = file(params.gtf)

    // =====================================================================
    // TRACK 1: Simulated reads
    // =====================================================================
    if (params.track in ['simulated', 'both']) {

        sim_fastq = file(params.simulated_fastq)
        sim_truth = file(params.simulated_truth)

        // Align simulated reads
        IngestSimulatedReads(
            Channel.of([
                params.error_profile,
                sim_fastq,
                sim_truth,
                genome_fa,
                genome_fai
            ])
        )

        // Run FLAIR transcriptome
        sim_flair_input = IngestSimulatedReads.out.aligned_sim
            .map { error_profile, bam, bai, truth_bed ->
                def label = "sim_${error_profile}"
                [label, bam, bai, genome_fa, genome_fai, gtf, params.flair_extra_args]
            }

        FlairTranscriptomeTest(sim_flair_input)

        // Assert boundary recovery
        sim_assert_input = FlairTranscriptomeTest.out.transcriptome_output
            .combine(
                IngestSimulatedReads.out.aligned_sim
                    .map { ep, bam, bai, truth -> ["sim_${ep}", truth] },
                by: 0
            )
            .map { label, isoforms_bed, isoforms_gtf, read_map, truth_bed ->
                [label, isoforms_bed, isoforms_gtf, read_map, truth_bed]
            }

        AssertBoundaries(sim_assert_input)

        // Assert end scoring on firstpass intermediate
        FlairTranscriptomeTest.out.firstpass_debug
            .map { label, firstpass_bed ->
                [label, firstpass_bed, genome_fa, genome_fai, params.library_type]
            }
            .set { sim_scoring_input }

        AssertEndScoring(sim_scoring_input)

        // Publish results
        AssertBoundaries.out.report
            .map { label, report -> report }
            .collectFile(name: 'sim_boundary_reports.json', storeDir: params.outdir)

        AssertEndScoring.out.report
            .map { label, report -> report }
            .collectFile(name: 'sim_scoring_reports.json', storeDir: params.outdir)
    }


    // =====================================================================
    // TRACK 2: Real data sub-region
    // =====================================================================
    if (params.track in ['real_subregion', 'both']) {

        bam_file = file(params.bam)
        bai_file = params.bai ? file(params.bai) : file("${params.bam}.bai")

        // Subset BAM + reference to the target region
        SubsetRealBam(
            Channel.of([
                params.sample_id,
                params.region,
                bam_file,
                bai_file,
                genome_fa,
                genome_fai,
                gtf
            ])
        )

        // Run FLAIR transcriptome on the subset
        real_flair_input = SubsetRealBam.out.subset_data
            .map { sample_id, region, bam, bai, sub_genome, sub_fai, sub_gtf, ref_bed ->
                def label = "real_${sample_id}_${region.replaceAll('[:-]', '_')}"
                [label, bam, bai, sub_genome, sub_fai, sub_gtf, params.flair_extra_args]
            }

        FlairTranscriptomeTest(real_flair_input)

        // Assert boundary recovery against annotation-derived ground truth
        real_assert_input = FlairTranscriptomeTest.out.transcriptome_output
            .combine(
                SubsetRealBam.out.subset_data
                    .map { sample_id, region, bam, bai, sub_genome, sub_fai, sub_gtf, ref_bed ->
                        ["real_${sample_id}_${region.replaceAll('[:-]', '_')}", ref_bed]
                    },
                by: 0
            )
            .map { label, isoforms_bed, isoforms_gtf, read_map, ref_bed ->
                [label, isoforms_bed, isoforms_gtf, read_map, ref_bed]
            }

        AssertBoundaries(real_assert_input)

        // Assert end scoring
        FlairTranscriptomeTest.out.firstpass_debug
            .map { label, firstpass_bed ->
                [label, firstpass_bed, genome_fa, genome_fai, params.library_type]
            }
            .set { real_scoring_input }

        AssertEndScoring(real_scoring_input)

        // Publish
        AssertBoundaries.out.report
            .map { label, report -> report }
            .collectFile(name: 'real_boundary_reports.json', storeDir: params.outdir)

        AssertEndScoring.out.report
            .map { label, report -> report }
            .collectFile(name: 'real_scoring_reports.json', storeDir: params.outdir)
    }
}
