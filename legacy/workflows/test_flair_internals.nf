#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * test_flair_internals.nf
 * =======================
 * Unified Nextflow workflow for testing FLAIR transcriptome internals.
 * Leverages the batching / workflow nature of Nextflow to run multiple
 * test tracks in parallel on SLURM.
 *
 * Three tracks:
 *   1. badread       — Badread-simulated reads with controlled error profiles
 *   2. real_subregion — BAM subset from a targeted genomic locus
 *   3. unit_tests     — pytest on flair-fusion's internal test suite
 *   4. end_trust      — Evaluate end-scoring trust across profiles × alphas
 *
 * Usage:
 *   # All tracks (tiny fixtures, fast CI-style)
 *   nextflow run workflows/test_flair_internals.nf \
 *       --track all \
 *       --genome tests/data/tiny_genome.fa \
 *       --gtf tests/data/tiny_annotation.gtf
 *
 *   # Badread simulation only
 *   nextflow run workflows/test_flair_internals.nf \
 *       --track badread \
 *       --genome tests/data/tiny_genome.fa \
 *       --gtf tests/data/tiny_annotation.gtf \
 *       --badread_scenarios clean_ont,clean_pacbio,high_chimera
 *
 *   # Pytest unit tests only
 *   nextflow run workflows/test_flair_internals.nf --track unit_tests
 *
 *   # Real data sub-region only
 *   nextflow run workflows/test_flair_internals.nf \
 *       --track real_subregion \
 *       --bam sample.bam \
 *       --region "chr22:20000000-21000000" \
 *       --genome GRCh38.fa \
 *       --gtf gencode.v48.gtf \
 *       --library_type ont_cDNA
 */

// -------------------------------------------------------------------------
// Module includes
// -------------------------------------------------------------------------
include {
    GenerateTranscriptomeFasta;
    BadreadSimulate;
    AlignSimulatedReads;
    FlairTranscriptomeTest as FlairSimTest;
    FlairTranscriptomeTest as FlairRealTest;
    FlairTranscriptomeTest as FlairEndTrustRun;
    AssertSimulatedResults;
    AssertBoundaries;
    AssertEndScoring;
    EvaluateEndTrust;
    PlotEndTrust;
    FlairWeightSweepRun;
    WeightSweepEval;
    CollectSweepResults;
    PlotWeightSweep;
    SubsetRealBam;
    RunPytest;
} from '../modules/transcriptome_testing/main'


// -------------------------------------------------------------------------
// Parameters
// -------------------------------------------------------------------------
params.track              = 'all'         // 'badread', 'real_subregion', 'unit_tests', 'end_trust', 'weight_sweep', or 'all'
params.genome             = null
params.gtf                = null

// Badread track
params.badread_depth      = 20            // reads per transcript
params.badread_quantity    = '5x'         // Badread --quantity (overrides depth-based calc)
params.badread_seed       = 42
params.badread_scenarios  = 'clean_ont,clean_pacbio,high_chimera,low_quality,truncated_5prime'

// Real data sub-region track
params.bam                = null
params.bai                = null
params.region             = 'chr22:20000000-21000000'
params.sample_id          = 'real_subregion'
params.library_type       = 'default'

// FLAIR arguments
params.flair_extra_args   = ''

// Weight sweep
params.sweep_profiles     = null          // list of profiles; null = all 5
params.sweep_alphas       = null          // list of floats; null = [0.0, 0.25, 0.5, 0.75, 1.0]
params.peaks_5prime       = null          // BED6 experimental 5' peaks (e.g. CAGE)
params.peaks_3prime       = null          // BED6 experimental 3' peaks (e.g. dRNA, dRNA)
params.peaks_3prime_alt   = null          // Optional second 3' peak set for comparison
params.peaks_3prime_alt_label = null      // Label for alt 3' peaks (e.g. 'dRNA')
params.tss_model          = null          // Path to XGBoost TSS boundary model (.pkl)
params.tts_model          = null          // Path to XGBoost TTS boundary model (.pkl)

// Internal test suite
params.flair_repo         = '/private/groups/brookslab/hdheath/tools/flair-fusion'

// Output
params.outdir             = "${launchDir}/results/flair_internals"


// -------------------------------------------------------------------------
// Scenario definitions
// -------------------------------------------------------------------------
// Each scenario is a map with all Badread parameters.
// Adding a new test case = adding one entry here.

def SCENARIO_CONFIGS = [
    clean_ont: [
        error_model: 'nanopore2023',
        identity:    '95,99,2.5',
        length:      '1500,1000',
        chimeras:    '0',
        junk_reads:  '0',
        random_reads:'0',
        glitches:    '0,0,0',
    ],
    clean_pacbio: [
        error_model: 'pacbio2021',
        identity:    '30,3',
        length:      '1500,1000',
        chimeras:    '0',
        junk_reads:  '0',
        random_reads:'0',
        glitches:    '0,0,0',
    ],
    high_chimera: [
        error_model: 'nanopore2023',
        identity:    '95,99,2.5',
        length:      '1500,1000',
        chimeras:    '15',
        junk_reads:  '1',
        random_reads:'1',
        glitches:    '10000,25,25',
    ],
    low_quality: [
        error_model: 'nanopore2023',
        identity:    '80,90,6',
        length:      '1500,1000',
        chimeras:    '1',
        junk_reads:  '5',
        random_reads:'5',
        glitches:    '1000,100,100',
    ],
    truncated_5prime: [
        error_model: 'nanopore2023',
        identity:    '95,99,2.5',
        length:      '500,250',
        chimeras:    '0',
        junk_reads:  '0',
        random_reads:'0',
        glitches:    '0,0,0',
    ],
    noisy_3prime: [
        error_model: 'nanopore2023',
        identity:    '95,99,2.5',
        length:      '1500,1000',
        chimeras:    '0',
        junk_reads:  '0',
        random_reads:'0',
        glitches:    '5000,50,50',
    ],
]


// -------------------------------------------------------------------------
// Validation
// -------------------------------------------------------------------------
def validateParams() {
    def tracks = params.track.tokenize(',')
    def valid = ['all', 'badread', 'real_subregion', 'unit_tests', 'end_trust', 'weight_sweep'] as Set
    tracks.each { t ->
        if (!(t in valid)) error "ERROR: Unknown track '${t}'. Valid: ${valid.join(', ')}"
    }

    def need_genome = tracks.intersect(['all', 'badread', 'real_subregion', 'end_trust'])
    if (need_genome) {
        if (!params.genome) error "ERROR: --genome is required for tracks: ${need_genome.join(', ')}"
        if (!params.gtf)    error "ERROR: --gtf is required for tracks: ${need_genome.join(', ')}"
    }

    if ('real_subregion' in tracks || 'all' in tracks) {
        // real_subregion is optional in 'all' mode — only runs if --bam provided
        if ('real_subregion' in tracks && !params.bam) {
            error "ERROR: --bam required for real_subregion track"
        }
    }

    if ('end_trust' in tracks) {
        if (!params.bam) error "ERROR: --bam required for end_trust track (needs firstpass BED from FLAIR)"
    }
}


// =========================================================================
// WORKFLOW
// =========================================================================
workflow {
    validateParams()

    def active_tracks = params.track.tokenize(',') as Set
    def run_badread     = 'all' in active_tracks || 'badread' in active_tracks
    def run_real        = ('all' in active_tracks && params.bam) || 'real_subregion' in active_tracks
    def run_pytest      = 'all' in active_tracks || 'unit_tests' in active_tracks
    def run_end_trust   = ('all' in active_tracks && params.bam) || 'end_trust' in active_tracks
    def run_weight_sweep = ('all' in active_tracks && params.bam) || 'weight_sweep' in active_tracks

    // Collect all assertion reports for final summary
    all_reports = Channel.empty()


    // =====================================================================
    // TRACK 1: Badread simulation
    // =====================================================================
    if (run_badread) {
        genome_fa  = file(params.genome)
        genome_fai = file("${params.genome}.fai")
        gtf        = file(params.gtf)

        // 1a. Generate transcriptome FASTA from genome + GTF
        GenerateTranscriptomeFasta(
            Channel.of([genome_fa, genome_fai, gtf, params.badread_depth])
        )

        // 1b. Build scenario channel from requested scenarios
        def requested = params.badread_scenarios.tokenize(',') as Set
        def scenario_list = SCENARIO_CONFIGS.findAll { name, cfg ->
            requested.contains(name)
        }.collect { name, cfg ->
            [name, cfg]
        }

        if (scenario_list.isEmpty()) {
            log.warn "No matching Badread scenarios found for: ${params.badread_scenarios}"
        }

        scenarios_ch = Channel.from(scenario_list)
            .combine(GenerateTranscriptomeFasta.out.transcriptome_ref)
            .map { name, cfg, txome_fa, gfa, gfai, an_gtf ->
                // [scenario_name, transcriptome_fa, quantity, error_model, identity,
                //  length, chimeras, junk, random, glitches, seed]
                [name, txome_fa,
                 params.badread_quantity,
                 cfg.error_model, cfg.identity, cfg.length,
                 cfg.chimeras, cfg.junk_reads, cfg.random_reads,
                 cfg.glitches, params.badread_seed]
            }

        // 1c. Simulate reads (each scenario runs as separate SLURM job)
        BadreadSimulate(scenarios_ch)

        // 1d. Align simulated reads to GENOME
        align_input = BadreadSimulate.out.simulated_reads
            .combine(GenerateTranscriptomeFasta.out.transcriptome_ref)
            .map { scenario_name, reads_fq, txome_fa, gfa, gfai, an_gtf ->
                [scenario_name, reads_fq, gfa, gfai]
            }

        AlignSimulatedReads(align_input)

        // 1e. Run FLAIR transcriptome on aligned reads
        flair_input = AlignSimulatedReads.out.aligned
            .combine(GenerateTranscriptomeFasta.out.transcriptome_ref)
            .map { scenario_name, bam, bai, txome_fa, gfa, gfai, an_gtf ->
                [scenario_name, bam, bai, gfa, gfai, an_gtf, params.flair_extra_args]
            }

        FlairSimTest(flair_input)

        // 1f. Assert results (scenario-aware)
        // The annotation GTF serves as ground truth since we simulated from it
        assert_input = FlairSimTest.out.transcriptome_output
            .combine(
                GenerateTranscriptomeFasta.out.transcriptome_ref
                    .map { txome_fa, gfa, gfai, an_gtf -> an_gtf }
            )
            .map { scenario_name, isoforms_bed, isoforms_gtf, read_map,
                   annotation_gtf ->
                [scenario_name, isoforms_bed, isoforms_gtf, read_map,
                 annotation_gtf]
            }

        AssertSimulatedResults(assert_input)

        all_reports = all_reports.mix(
            AssertSimulatedResults.out.report
                .map { name, report -> report }
        )

        // 1h. End scoring assertion on firstpass intermediates (if available)
        FlairSimTest.out.firstpass_debug
            .combine(
                GenerateTranscriptomeFasta.out.transcriptome_ref
                    .map { txome_fa, gfa, gfai, an_gtf -> [gfa, gfai] }
            )
            .map { label, firstpass_bed, gfa, gfai ->
                [label, firstpass_bed, gfa, gfai, params.library_type]
            }
            .set { scoring_input }

        AssertEndScoring(scoring_input)
    }


    // =====================================================================
    // TRACK 2: Real data sub-region
    // =====================================================================
    if (run_real) {
        genome_fa  = file(params.genome)
        genome_fai = file("${params.genome}.fai")
        gtf        = file(params.gtf)
        bam_file   = file(params.bam)
        bai_file   = params.bai ? file(params.bai) : file("${params.bam}.bai")

        SubsetRealBam(
            Channel.of([
                params.sample_id, params.region,
                bam_file, bai_file,
                genome_fa, genome_fai, gtf
            ])
        )

        real_flair_input = SubsetRealBam.out.subset_data
            .map { sample_id, region, bam, bai, sub_genome, sub_fai, sub_gtf, ref_bed ->
                def label = "real_${sample_id}_${region.replaceAll('[:-]', '_')}"
                [label, bam, bai, sub_genome, sub_fai, sub_gtf, params.flair_extra_args]
            }

        FlairRealTest(real_flair_input)

        real_assert_input = FlairRealTest.out.transcriptome_output
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

        all_reports = all_reports.mix(
            AssertBoundaries.out.report
                .map { name, report -> report }
        )
    }


    // =====================================================================
    // TRACK 3: Pytest unit tests on flair-fusion internals
    // =====================================================================
    if (run_pytest) {
        RunPytest(Channel.of(file(params.flair_repo)))
    }


    // =====================================================================
    // TRACK 4: End Trust Evaluation
    // =====================================================================
    // Runs the scoring sweep (all profiles × alphas) on FLAIR firstpass
    // output, then generates pub-quality plots.
    //
    // Re-uses real data track's FLAIR run when both are active; otherwise
    // runs its own FLAIR on the provided BAM.
    if (run_end_trust) {
        genome_fa  = file(params.genome)
        genome_fai = file("${params.genome}.fai")
        gtf        = file(params.gtf)
        bam_file   = file(params.bam)
        bai_file   = params.bai ? file(params.bai) : file("${params.bam}.bai")

        if (run_real) {
            // Reuse firstpass from the real-subregion FLAIR run
            end_trust_input = FlairRealTest.out.firstpass_debug
                .combine(
                    SubsetRealBam.out.subset_data
                        .map { sample_id, region, bam, bai, sub_genome, sub_fai, sub_gtf, ref_bed ->
                            ["real_${sample_id}_${region.replaceAll('[:-]', '_')}", sub_genome, sub_fai, sub_gtf]
                        },
                    by: 0
                )
                .map { label, firstpass_bed, sub_genome, sub_fai, sub_gtf ->
                    ["end_trust_${label}", firstpass_bed, sub_genome, sub_fai, sub_gtf]
                }
        } else {
            // Standalone: run FLAIR on the provided BAM to get firstpass
            def et_label = "end_trust_${params.sample_id}"

            et_flair_ch = Channel.of([
                et_label,
                bam_file, bai_file,
                genome_fa, genome_fai, gtf,
                "--keep_intermediate ${params.flair_extra_args}".trim()
            ])

            FlairEndTrustRun(et_flair_ch)

            end_trust_input = FlairEndTrustRun.out.firstpass_debug
                .map { label, firstpass_bed ->
                    [label, firstpass_bed, genome_fa, genome_fai, gtf]
                }
        }

        EvaluateEndTrust(end_trust_input)

        PlotEndTrust(EvaluateEndTrust.out.end_trust_results)
    }


    // =====================================================================
    // TRACK 5: Weight Sweep
    // =====================================================================
    // Runs FLAIR with every (profile, alpha) combination, evaluates each
    // via per-SQANTI-category precision/recall + end redundancy, then
    // collects results into a combined summary TSV and generates heatmap
    // / Pareto plots.
    if (run_weight_sweep) {
        genome_fa  = file(params.genome)
        genome_fai = file("${params.genome}.fai")
        gtf        = file(params.gtf)
        bam_file   = file(params.bam)
        bai_file   = params.bai ? file(params.bai) : file("${params.bam}.bai")

        // Resolve optional XGBoost model paths
        def tss_model = params.tss_model ? file(params.tss_model) : file('NO_TSS_MODEL')
        def tts_model = params.tts_model ? file(params.tts_model) : file('NO_TTS_MODEL')

        // Build the (profile, alpha) grid
        def profiles = params.sweep_profiles ?: [
            'default', 'ont_cDNA', 'ont_dRNA', 'pacbio_isoseq', 'pacbio_masseq'
        ]
        def alphas = params.sweep_alphas ?: [0.0, 0.25, 0.5, 0.75, 1.0]

        // Create channel of all (profile, alpha) pairs × shared inputs
        sweep_grid = Channel.fromList(
            [profiles, alphas].combinations().collect { combo ->
                [combo[0], combo[1],
                 bam_file, bai_file, genome_fa, genome_fai, gtf,
                 params.flair_extra_args ?: '']
            }
        )

        FlairWeightSweepRun(sweep_grid, tss_model, tts_model)

        // Resolve optional peaks files
        def peaks5  = params.peaks_5prime  ? file(params.peaks_5prime)  : file('NO_PEAKS_5')
        def peaks3  = params.peaks_3prime  ? file(params.peaks_3prime)  : file('NO_PEAKS_3')

        // Evaluate each (profile, alpha) output
        sweep_eval_input = FlairWeightSweepRun.out.sweep_isoforms
            .map { profile, alpha, isoforms_bed, isoforms_gtf ->
                [profile, alpha, isoforms_bed, gtf, peaks5, peaks3]
            }

        WeightSweepEval(sweep_eval_input)

        // Collect all per-run TSVs — rename to prevent filename collisions
        pr_tsvs = WeightSweepEval.out.sweep_metrics
            .map { profile, alpha, pr_tsv, redund_tsv ->
                // Stage with unique names so CollectSweepResults can distinguish them
                pr_tsv
            }
            .collect()

        redund_tsvs = WeightSweepEval.out.sweep_metrics
            .map { profile, alpha, pr_tsv, redund_tsv ->
                redund_tsv
            }
            .collect()

        CollectSweepResults(pr_tsvs, redund_tsvs)

        PlotWeightSweep(CollectSweepResults.out.sweep_summary)
    }


    // =====================================================================
    // Collect all reports
    // =====================================================================
    all_reports
        .collectFile(name: 'all_assertion_reports.txt', storeDir: params.outdir,
                     newLine: true)
}
