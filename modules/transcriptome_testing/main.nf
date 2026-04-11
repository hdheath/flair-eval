/*
 * Nextflow processes for testing FLAIR transcriptome internals.
 *
 * Three testing tracks:
 *   1. BadreadSimulation — Badread-generated reads with controlled error profiles
 *   2. RealDataSubregion — targeted genomic sub-regions from real BAMs
 *   3. UnitTests         — pytest on flair-fusion's internal test suite
 *
 * Track 1 flow:
 *   GenerateTranscriptomeFasta → BadreadSimulate → AlignSimulatedReads
 *     → FlairTranscriptomeTest → AssertSimulatedResults
 *
 * Track 2 flow:
 *   SubsetRealBam → FlairTranscriptomeTest → AssertBoundaries
 *
 * Track 3 flow:
 *   RunPytest → report
 */

nextflow.enable.dsl = 2


// =========================================================================
// TRACK 1: Badread simulation pipeline
// =========================================================================

/*
 * GenerateTranscriptomeFasta
 * --------------------------
 * Extract spliced cDNA sequences from genome + GTF.  The output FASTA has
 * per-transcript depth= headers that Badread reads natively.
 */
process GenerateTranscriptomeFasta {
    tag "gen_txome"
    label 'process_low'

    conda "${params.conda_flair}"

    input:
    tuple path(genome_fa), path(genome_fai), path(annotation_gtf), val(depth)

    output:
    tuple path("transcriptome.fa"), path(genome_fa), path(genome_fai), path(annotation_gtf),
          emit: transcriptome_ref

    script:
    """
    python3 ${projectDir}/../bin/transcriptome_testing/generate_transcriptome_fa.py \\
        --genome ${genome_fa} \\
        --gtf ${annotation_gtf} \\
        --output transcriptome.fa \\
        --depth ${depth}
    """
}


/*
 * BadreadSimulate
 * ---------------
 * Generate synthetic long reads using Badread.  Each invocation uses a
 * specific "scenario" — a named combination of error model, identity,
 * chimera rate, glitch params, and fragment length distribution.
 *
 * Scenarios are defined as channel tuples; adding a new test case is
 * one line in the workflow.
 */
process BadreadSimulate {
    tag "${scenario_name}"
    label 'process_low'

    conda "${params.conda_analysis}"

    input:
    tuple val(scenario_name),
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
    tuple val(scenario_name),
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
 * AlignSimulatedReads
 * -------------------
 * Splice-align Badread output to the GENOME (not transcriptome) using
 * minimap2.  This mirrors real FLAIR usage where reads are genome-aligned.
 */
process AlignSimulatedReads {
    tag "${scenario_name}"
    label 'process_medium'

    conda "${params.conda_flair}"

    input:
    tuple val(scenario_name),
          path(reads_fastq),
          path(genome_fa),
          path(genome_fai)

    output:
    tuple val(scenario_name),
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


// =========================================================================
// TRACK 3: Internal pytest suite
// =========================================================================

/*
 * RunPytest
 * ---------
 * Run pytest on flair-fusion's internal unit test suite.
 * Emits a JSON report and JUnit XML so results are collected alongside
 * the simulation tests.
 */
process RunPytest {
    tag "unit_tests"
    label 'process_low'

    conda "${params.conda_flair}"

    input:
    path(flair_repo)

    output:
    tuple path("pytest_report.json"),
          path("pytest_junit.xml"),
          emit: report

    script:
    """
    cd ${flair_repo}
    python -m pytest tests/unit/ -v \\
        --tb=short \\
        --json-report --json-report-file=\${OLDPWD}/pytest_report.json \\
        --junitxml=\${OLDPWD}/pytest_junit.xml \\
        || true
    cd \${OLDPWD}

    # Ensure files exist even if pytest-json-report not installed
    if [ ! -f pytest_report.json ]; then
        python -c "
import json, subprocess, sys
result = subprocess.run(
    [sys.executable, '-m', 'pytest', '${flair_repo}/tests/unit/', '-v', '--tb=short'],
    capture_output=True, text=True
)
report = {
    'exitcode': result.returncode,
    'stdout': result.stdout[-5000:] if len(result.stdout) > 5000 else result.stdout,
    'stderr': result.stderr[-2000:] if len(result.stderr) > 2000 else result.stderr,
    'passed': 'passed' in result.stdout,
}
with open('pytest_report.json', 'w') as f:
    json.dump(report, f, indent=2)
"
    fi
    if [ ! -f pytest_junit.xml ]; then
        echo '<testsuites><testsuite name=\"fallback\" tests=\"0\"/></testsuites>' > pytest_junit.xml
    fi
    """
}


// =========================================================================
// Assertion for Badread simulation results
// =========================================================================

/*
 * AssertSimulatedResults
 * ----------------------
 * Scenario-aware assertion: compares FLAIR output against the known
 * transcriptome used to generate Badread reads.
 *
 * Checks vary by scenario:
 *   - clean_*:       high junction recall, low TSS/TTS displacement
 *   - high_chimera:  chimeric isoforms should be filtered
 *   - low_quality:   graceful degradation, no crashes
 *   - truncated_*:   TSS displacement expected to be larger
 */
process AssertSimulatedResults {
    tag "${scenario_name}"
    label 'process_low'

    conda "${params.conda_flair}"

    input:
    tuple val(scenario_name),
          path(flair_isoforms_bed),
          path(flair_isoforms_gtf),
          path(flair_read_map),
          path(annotation_gtf)

    output:
    tuple val(scenario_name),
          path("${scenario_name}_assertion_report.json"),
          emit: report

    script:
    """
    python3 ${projectDir}/../bin/transcriptome_testing/assert_simulated_results.py \\
        --flair-bed ${flair_isoforms_bed} \\
        --annotation-gtf ${annotation_gtf} \\
        --read-map ${flair_read_map} \\
        --scenario ${scenario_name} \\
        --output ${scenario_name}_assertion_report.json
    """
}


// =========================================================================
// Existing processes (SubsetReference, IngestSimulatedReads, etc.)
// =========================================================================

/*
 * SubsetReference
 * ---------------
 * Extract a small genomic region from the full reference to serve as the
 * mock reference for simulation and assembly.
 */
process SubsetReference {
    tag "${region}"
    label 'process_low'

    conda "${params.conda_flair ?: 'flair-dev'}"

    input:
    tuple val(region), path(genome_fa), path(genome_fai), path(annotation_gtf)

    output:
    tuple val(region),
          path("subset_genome.fa"),
          path("subset_genome.fa.fai"),
          path("subset_annotation.gtf"),
          emit: subset_ref

    script:
    """
    # Extract region from genome
    samtools faidx ${genome_fa} ${region} > subset_genome.fa
    samtools faidx subset_genome.fa

    # Subset GTF to region
    chrom=\$(echo "${region}" | cut -d: -f1)
    start=\$(echo "${region}" | cut -d: -f2 | cut -d- -f1)
    end=\$(echo "${region}" | cut -d: -f2 | cut -d- -f2)

    awk -v chr="\$chrom" -v s="\$start" -v e="\$end" \\
        '\$1 == chr && \$4 >= s && \$5 <= e' \\
        ${annotation_gtf} > subset_annotation.gtf || true

    # If no annotation lines matched, create a minimal file
    if [ ! -s subset_annotation.gtf ]; then
        echo "# No annotation features in region ${region}" > subset_annotation.gtf
    fi
    """
}


/*
 * IngestSimulatedReads
 * --------------------
 * Wrapper process that accepts pre-generated simulated FASTQs with
 * explicitly injected error profiles.  Does NOT hardcode a simulation tool.
 *
 * The user provides:
 *   - simulated FASTQ (from badread, NanoSim, pbsim3, or a custom generator)
 *   - a ground-truth BED12 of the expected isoform boundaries
 *   - a label describing the error profile (e.g. "3prime_a_rich", "5prime_softclip")
 *
 * This process aligns the reads to the reference using minimap2.
 */
process IngestSimulatedReads {
    tag "${error_profile}"
    label 'process_medium'

    conda "${params.conda_flair ?: 'flair-dev'}"

    input:
    tuple val(error_profile),
          path(simulated_fastq),
          path(ground_truth_bed),
          path(genome_fa),
          path(genome_fai)

    output:
    tuple val(error_profile),
          path("simulated_aligned.bam"),
          path("simulated_aligned.bam.bai"),
          path(ground_truth_bed),
          emit: aligned_sim

    script:
    """
    minimap2 -ax splice --secondary=no -t ${task.cpus} \\
        ${genome_fa} ${simulated_fastq} \\
        | samtools sort -@ ${task.cpus} -o simulated_aligned.bam

    samtools index simulated_aligned.bam
    """
}


// =========================================================================
// TRACK 2: Real data sub-region extraction
// =========================================================================

/*
 * SubsetRealBam
 * -------------
 * Extract reads from a targeted genomic locus of a real BAM file.
 * This provides biological noise that simulated reads lack.
 *
 * Typical use: take a ~1 MB region from a WTC11 or A549 long-read alignment
 * that contains known complex overlapping genes.
 */
process SubsetRealBam {
    tag "${sample_id}:${region}"
    label 'process_medium'

    conda "${params.conda_flair ?: 'flair-dev'}"

    input:
    tuple val(sample_id),
          val(region),
          path(bam),
          path(bai),
          path(genome_fa),
          path(genome_fai),
          path(annotation_gtf)

    output:
    tuple val(sample_id),
          val(region),
          path("subset.bam"),
          path("subset.bam.bai"),
          path("subset_genome.fa"),
          path("subset_genome.fa.fai"),
          path("subset_annotation.gtf"),
          path("subset_reference_isoforms.bed"),
          emit: subset_data

    script:
    """
    # 1. Extract reads overlapping the region
    samtools view -b -h ${bam} ${region} > subset.bam
    samtools index subset.bam

    # 2. Extract reference genome for the region
    samtools faidx ${genome_fa} ${region} > subset_genome.fa
    samtools faidx subset_genome.fa

    # 3. Subset the GTF annotation
    chrom=\$(echo "${region}" | cut -d: -f1)
    start=\$(echo "${region}" | cut -d: -f2 | cut -d- -f1)
    end=\$(echo "${region}" | cut -d: -f2 | cut -d- -f2)

    awk -v chr="\$chrom" -v s="\$start" -v e="\$end" \\
        '\$1 == chr && \$4 >= s && \$5 <= e' \\
        ${annotation_gtf} > subset_annotation.gtf || true

    if [ ! -s subset_annotation.gtf ]; then
        echo "# No annotation in ${region}" > subset_annotation.gtf
    fi

    # 4. Convert GTF annotation to BED12 for ground-truth comparison
    python3 -c "
import sys
transcripts = {}
for line in open('subset_annotation.gtf'):
    if line.startswith('#'): continue
    f = line.strip().split('\\t')
    if len(f) < 9: continue
    chrom, _, feat, start, end, _, strand = f[0], f[1], f[2], int(f[3])-1, int(f[4]), f[5], f[6]
    tid = ''
    for token in f[8].split(';'):
        token = token.strip()
        if token.startswith('transcript_id'):
            tid = token.split('\"')[1]
            break
    if not tid: continue
    if feat == 'transcript':
        transcripts[tid] = {'chrom': chrom, 'start': start, 'end': end,
                            'strand': strand, 'exons': []}
    elif feat == 'exon' and tid in transcripts:
        transcripts[tid]['exons'].append((start, end))

with open('subset_reference_isoforms.bed', 'w') as out:
    for tid, info in transcripts.items():
        if not info['exons']: continue
        exons = sorted(info['exons'])
        t_start = exons[0][0]
        t_end = exons[-1][1]
        sizes = ','.join(str(e[1]-e[0]) for e in exons)
        starts = ','.join(str(e[0]-t_start) for e in exons)
        out.write(f\"{info['chrom']}\\t{t_start}\\t{t_end}\\t{tid}\\t0\\t{info['strand']}\\t{t_start}\\t{t_end}\\t0\\t{len(exons)}\\t{sizes},\\t{starts},\\n\")
"
    """
}


// =========================================================================
// Shared: Run FLAIR transcriptome on test inputs
// =========================================================================

/*
 * FlairTranscriptomeTest
 * ----------------------
 * Run `flair transcriptome` on either simulated or real sub-region data.
 * This is the code-under-test: it exercises the end-detection and
 * collapse_end_groups() logic that we are scoring.
 */
process FlairTranscriptomeTest {
    tag "${test_label}"
    label 'process_medium'

    conda "${params.conda_flair ?: 'flair-dev'}"

    input:
    tuple val(test_label),
          path(bam),
          path(bai),
          path(genome_fa),
          path(genome_fai),
          path(annotation_gtf),
          val(extra_args)

    output:
    tuple val(test_label),
          path("${test_label}.isoforms.bed"),
          path("${test_label}.isoforms.gtf"),
          path("${test_label}.isoform.read.map.txt"),
          emit: transcriptome_output
    tuple val(test_label),
          path("${test_label}.firstpass.unfiltered.bed"),
          optional: true,
          emit: firstpass_debug

    script:
    def keep_flag = "--keep_intermediate"
    """
    flair transcriptome \\
        -b ${bam} \\
        -g ${genome_fa} \\
        -f ${annotation_gtf} \\
        -o ${test_label} \\
        -t ${task.cpus} \\
        ${keep_flag} \\
        ${extra_args}
    """
}


// =========================================================================
// Assertion / validation processes
// =========================================================================

/*
 * AssertBoundaries
 * ----------------
 * Compare FLAIR output isoforms against ground-truth boundaries.
 *
 * For simulated reads: asserts that mock GTF boundaries are perfectly recovered
 *   and that synthetic artifacts (injected A-richness, soft-clips) are filtered.
 *
 * For real data: asserts that the pipeline doesn't crash and retains highly
 *   supported canonical transcripts (no false negatives on strong signal).
 */
process AssertBoundaries {
    tag "${test_label}"
    label 'process_low'

    conda "${params.conda_flair ?: 'flair-dev'}"

    input:
    tuple val(test_label),
          path(flair_isoforms_bed),
          path(flair_isoforms_gtf),
          path(flair_read_map),
          path(ground_truth_bed)

    output:
    tuple val(test_label),
          path("assertion_report.json"),
          emit: report

    script:
    """
    python3 ${projectDir}/../bin/transcriptome_testing/assert_boundaries.py \\
        --flair-bed ${flair_isoforms_bed} \\
        --truth-bed ${ground_truth_bed} \\
        --read-map ${flair_read_map} \\
        --output assertion_report.json \\
        --test-label ${test_label} \\
        --tss-tolerance 100 \\
        --tts-tolerance 100
    """
}


/*
 * AssertEndScoring
 * ----------------
 * Run the end_scoring module on candidate boundaries from the firstpass
 * output and validate that the confidence scores behave as expected.
 *
 * This is the unit-integration bridge: it runs the Python scoring logic
 * within the Nextflow framework on actual FLAIR intermediate files.
 */
process AssertEndScoring {
    tag "${test_label}"
    label 'process_low'

    conda "${params.conda_flair ?: 'flair-dev'}"

    input:
    tuple val(test_label),
          path(firstpass_bed),
          path(genome_fa),
          path(genome_fai),
          val(library_type)

    output:
    tuple val(test_label),
          path("scoring_report.json"),
          emit: report

    script:
    """
    python3 ${projectDir}/../bin/transcriptome_testing/assert_end_scoring.py \\
        --firstpass-bed ${firstpass_bed} \\
        --genome ${genome_fa} \\
        --library-type ${library_type} \\
        --output scoring_report.json \\
        --test-label ${test_label}
    """
}


// =========================================================================
// End Trust Evaluation — empirical alpha / profile sweep
// =========================================================================

/*
 * EvaluateEndTrust
 * ----------------
 * Score every candidate boundary from the firstpass BED under multiple
 * library profiles × alpha values.  Measures displacement from annotation.
 *
 * Outputs:
 *   - per-boundary TSV  (one row per boundary × profile × alpha)
 *   - aggregate summary TSV
 */
process EvaluateEndTrust {
    tag "${test_label}"
    label 'process_low'

    conda "${params.conda_flair ?: 'flair-dev'}"

    publishDir "${params.outdir}/${test_label}/end_trust", mode: 'copy'

    input:
    tuple val(test_label),
          path(firstpass_bed),
          path(genome_fa),
          path(genome_fai),
          path(annotation_gtf)

    output:
    tuple val(test_label),
          path("${test_label}_boundaries.tsv"),
          path("${test_label}_summary.tsv"),
          emit: end_trust_results

    script:
    """
    python3 ${projectDir}/bin/evaluation/evaluate_end_trust.py \\
        --bed ${firstpass_bed} \\
        --genome ${genome_fa} \\
        --gtf ${annotation_gtf} \\
        --output ${test_label}_boundaries.tsv \\
        --summary ${test_label}_summary.tsv
    """
}


/*
 * PlotEndTrust
 * ------------
 * Generate publication-quality plots from EvaluateEndTrust output.
 *
 * Produces 8 PNG files:
 *   1. end_trust_by_library.png
 *   2. displacement_distributions.png
 *   3. confidence_vs_displacement.png
 *   4. alpha_sweep_displacement.png
 *   5. trust_category_distribution.png
 *   6. weight_decomposition.png
 *   7. rescue_effectiveness.png
 *   8. end_trust_dashboard.png  (combined 2×2 panel)
 */
process PlotEndTrust {
    tag "${test_label}"
    label 'process_low'

    conda "${params.conda_analysis ?: 'nextflow_env'}"

    publishDir "${params.outdir}/${test_label}/end_trust/plots", mode: 'copy'

    input:
    tuple val(test_label),
          path(boundaries_tsv),
          path(summary_tsv)

    output:
    path "*.png", emit: plots

    script:
    """
    python3 ${projectDir}/bin/evaluation/end_trust_plots.py \\
        --boundaries ${boundaries_tsv} \\
        --summary ${summary_tsv} \\
        --outdir .
    """
}


// =========================================================================
// Weight Sweep processes
// =========================================================================

/*
 * FlairWeightSweepRun
 * --------------------
 * Runs FLAIR transcriptome with a specific (profile, alpha) combination.
 * Used by the weight_sweep track to generate isoforms under each config.
 */
process FlairWeightSweepRun {
    tag "${profile}_alpha${alpha}"
    label 'process_medium'

    conda "${params.conda_flair ?: 'flair-dev'}"

    input:
    tuple val(profile), val(alpha),
          path(bam), path(bai),
          path(genome_fa), path(genome_fai),
          path(annotation_gtf),
          val(extra_args)
    path tss_model, stageAs: 'tss_model.pkl'
    path tts_model, stageAs: 'tts_model.pkl'

    output:
    tuple val(profile), val(alpha),
          path("sweep_${profile}_a${alpha}.isoforms.bed"),
          path("sweep_${profile}_a${alpha}.isoforms.gtf"),
          emit: sweep_isoforms
    tuple val(profile), val(alpha),
          path("sweep_${profile}_a${alpha}.firstpass.unfiltered.bed"),
          optional: true,
          emit: sweep_firstpass

    script:
    def alpha_arg = alpha > 0 ? "--end_scoring_alpha ${alpha} --library_type ${profile}" : ""
    def model_args = ""
    if (tss_model.name != 'NO_TSS_MODEL') {
        model_args += " --tss_model ${tss_model}"
    }
    if (tts_model.name != 'NO_TTS_MODEL') {
        model_args += " --tts_model ${tts_model}"
    }
    """
    flair transcriptome \\
        -b ${bam} \\
        -g ${genome_fa} \\
        -f ${annotation_gtf} \\
        -o sweep_${profile}_a${alpha} \\
        -t ${task.cpus} \\
        --keep_intermediate \\
        ${alpha_arg} \\
        ${model_args} \\
        ${extra_args}
    """
}


/*
 * WeightSweepEval
 * ----------------
 * Runs per-SQANTI-category precision/recall and end-redundancy analysis
 * on a single (profile, alpha) FLAIR output.
 */
process WeightSweepEval {
    tag "${profile}_alpha${alpha}"
    label 'process_low'

    conda "${params.conda_flair ?: 'flair-dev'}"

    publishDir "${params.outdir}/weight_sweep/per_run/${profile}_a${alpha}", mode: 'copy', pattern: '*.png'

    input:
    tuple val(profile), val(alpha),
          path(isoforms_bed),
          path(annotation_gtf),
          path(peaks_5prime),
          path(peaks_3prime)

    output:
    tuple val(profile), val(alpha),
          path("${profile}_a${alpha}_pr.tsv"),
          path("${profile}_a${alpha}_redund.tsv"),
          emit: sweep_metrics
    path "*.png", emit: sweep_plots

    script:
    def peaks5_arg = peaks_5prime.name != 'NO_PEAKS_5' ? "--peaks-5prime ${peaks_5prime}" : ""
    def peaks3_arg = peaks_3prime.name != 'NO_PEAKS_3' ? "--peaks-3prime ${peaks_3prime}" : ""
    """
    python ${projectDir}/bin/evaluation/end_scoring_precision_recall.py \\
        --isoforms-bed ${isoforms_bed} \\
        --gtf ${annotation_gtf} \\
        --library-type ${profile} \\
        --window 50 \\
        --outdir . \\
        ${peaks5_arg} \\
        ${peaks3_arg}

    # The script now writes precision_recall_by_category.${profile}.tsv and
    # end_redundancy.${profile}.tsv with library_type column already present.
    # Rename with alpha tag and inject profile+alpha columns for sweep merge.
    python3 -c "
import csv
for src, dst_name in [('precision_recall_by_category.${profile}.tsv', '${profile}_a${alpha}_pr.tsv'),
                       ('end_redundancy.${profile}.tsv', '${profile}_a${alpha}_redund.tsv')]:
    try:
        with open(src) as f:
            reader = csv.DictReader(f, delimiter='\\t')
            rows = list(reader)
            fieldnames = ['profile', 'alpha'] + [c for c in reader.fieldnames if c not in ('profile', 'alpha')]
        with open(dst_name, 'w', newline='') as out:
            writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\\t')
            writer.writeheader()
            for row in rows:
                row['profile'] = '${profile}'
                row['alpha'] = '${alpha}'
                writer.writerow(row)
    except FileNotFoundError:
        with open(dst_name, 'w') as out:
            out.write('profile\\talpha\\terror\\n')
            out.write('${profile}\\t${alpha}\\tno_data\\n')
"
    """
}


/*
 * CollectSweepResults
 * --------------------
 * Merges per-(profile, alpha) precision/recall and redundancy TSVs
 * into a single combined summary for sweep-level plotting.
 */
process CollectSweepResults {
    tag "collect_sweep"
    label 'process_low'

    conda "${params.conda_flair ?: 'flair-dev'}"

    publishDir "${params.outdir}/weight_sweep", mode: 'copy'

    input:
    path pr_tsvs       // *_pr.tsv files with profile/alpha columns pre-injected
    path redund_tsvs   // *_redund.tsv files with profile/alpha columns pre-injected

    output:
    path "weight_sweep_summary.tsv", emit: sweep_summary

    script:
    """
    python3 -c "
import csv, glob

# Each TSV already has 'profile' and 'alpha' columns injected by WeightSweepEval.
# Merge all precision/recall TSVs into one summary, then append redundancy stats.
pr_files = sorted(glob.glob('*_pr.tsv'))
redund_files = sorted(glob.glob('*_redund.tsv'))

# ─── Build per-(profile, alpha) redundancy lookup ───
redund_lookup = {}
for f in redund_files:
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter='\\t')
        for row in reader:
            key = (row.get('profile', ''), row.get('alpha', ''))
            redund_lookup[key] = row

# ─── Merge precision/recall rows, enrich with redundancy ───
rows = []
for f in pr_files:
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter='\\t')
        for row in reader:
            key = (row.get('profile', ''), row.get('alpha', ''))
            rd = redund_lookup.get(key, {})
            row['total_groups'] = rd.get('total_groups', '')
            row['groups_with_redundancy'] = rd.get('groups_with_redundancy', '')
            row['redundant_tss'] = rd.get('redundant_tss', '')
            row['redundant_tts'] = rd.get('redundant_tts', '')
            rows.append(row)

if rows:
    fieldnames = list(rows[0].keys())
    with open('weight_sweep_summary.tsv', 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\\t')
        writer.writeheader()
        writer.writerows(rows)
else:
    with open('weight_sweep_summary.tsv', 'w') as out:
        out.write('profile\\talpha\\terror\\n')
        out.write('none\\tnone\\tno_input_files\\n')
"
    """
}


/*
 * PlotWeightSweep
 * ----------------
 * Generates pub-quality sweep summary plots from the combined TSV:
 *   1. Heatmap: alpha × profile → mean F1
 *   2. Heatmap: alpha × profile → end redundancy rate
 *   3. Pareto frontier: F1 vs redundancy
 *   4. Per-category sensitivity across the sweep
 */
process PlotWeightSweep {
    tag "plot_sweep"
    label 'process_low'

    conda "${params.conda_analysis ?: 'nextflow_env'}"

    publishDir "${params.outdir}/weight_sweep/plots", mode: 'copy'

    input:
    path sweep_summary

    output:
    path "*.png", emit: plots

    script:
    """
    python3 ${projectDir}/bin/evaluation/plot_weight_sweep.py \\
        --summary ${sweep_summary} \\
        --outdir .
    """
}
