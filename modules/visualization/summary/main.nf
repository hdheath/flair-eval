// Module: SummaryPlots + PeakReasonHeatmap + TpOverlap + IsoformsPerGeneHist
//         + JaccardHeatmap + TotalIsoforms + EndSignalScatter + CumulativeSignal
//         + CombineEvaluationTSVs + InternalPrimingAnalysis + SignalReadSupport + PrSignalBalance
//         + SqantiPrecision + DepthCalibration
// Cross-assembler comparison plots generated after all evaluations complete.

/*
 * CombineEvaluationTSVs
 * ---------------------
 * Concatenates all per-method evaluation TSVs for a test into one file.
 * Header is taken from the first file; subsequent files contribute data rows only.
 */
process CombineEvaluationTSVs {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary", mode: 'copy'
    tag "${test_name}"

    input:
    tuple val(test_name), path(evaluation_files)

    output:
    path "${test_name}_combined_evaluation.tsv", emit: combined_tsv

    script:
    """
    python3 -c "
import csv, sys
from collections import OrderedDict

files = '${evaluation_files}'.split()
rows = []
all_cols = OrderedDict()
for f in files:
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for col in reader.fieldnames:
            all_cols[col] = None
        for row in reader:
            rows.append(row)

cols = list(all_cols.keys())
with open('${test_name}_combined_evaluation.tsv', 'w', newline='') as out:
    w = csv.writer(out, delimiter='\t')
    w.writerow(cols)
    for row in rows:
        w.writerow([row.get(c, '') for c in cols])
"
    """
}

/*
 * CombinePrecisionRecall
 * ----------------------
 * Concatenates all per-method precision_recall_summary.tsv files for a sample
 * into a single TSV (header once, one data row per method).
 */
process CombinePrecisionRecall {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/end_accuracy/precision_recall", mode: 'copy'
    tag "${test_name}"

    input:
    tuple val(test_name), path(pr_tsvs)
    val(suffix)   // "ortho" or "gtf" — distinguishes the two combined TSV files

    output:
    tuple val(test_name), path("${test_name}_precision_recall_${suffix}.tsv"), emit: combined_pr

    script:
    """
    python3 -c "
import sys, csv
files = '${pr_tsvs}'.split()
header = None
rows = []
for f in files:
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter='\\t')
        if header is None:
            header = reader.fieldnames
        for row in reader:
            rows.append(row)
rows.sort(key=lambda r: r.get('transcriptome_mode', ''))
with open('${test_name}_precision_recall_${suffix}.tsv', 'w', newline='') as out:
    writer = csv.DictWriter(out, fieldnames=header, delimiter='\\t')
    writer.writeheader()
    writer.writerows(rows)
"
    """
}

/*
 * PrecisionRecallPlot
 * -------------------
 * Generates 5'/3' P/R scatter and F1 bar plots from JC-deduplicated
 * precision_recall_summary.tsv files produced by TedEndPrecision.
 * Uses interval-edge matching with a 50bp window — no per-isoform naive counting.
 */
process PrecisionRecallPlot {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/end_accuracy/precision_recall", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(precision_recall_tsvs), path(gtf_precision_recall_tsvs)

    output:
    path "*.png", emit: precision_recall_plot, optional: true

    script:
    """
    # v5: add paired P/R scatter plot
    python ${projectDir}/bin/evaluation/precision_recall_plot.py \\
        --input ${precision_recall_tsvs} \\
        --gtf-input ${gtf_precision_recall_tsvs} \\
        --output . \\
        --verbose
    """
}

process SummaryPlots {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/comparison", mode: 'copy'
    tag "${test_name}"

    input:
    tuple val(test_name), path(evaluation_files)

    output:
    path "concordance/*.png", emit: concordance_plot, optional: true
    path "signal_support/*.png", emit: signal_support_plot, optional: true
    path "landscape/*.png", emit: landscape_plot, optional: true

    script:
    """
    python ${projectDir}/bin/evaluation/concordance_plots.py \\
        --input ${evaluation_files} \\
        --output concordance \\
        --verbose || true

    python ${projectDir}/bin/evaluation/signal_support_summary_plot.py \\
        --input ${evaluation_files} \\
        --output signal_support \\
        --verbose || true

    python ${projectDir}/bin/evaluation/transcriptome_landscape_plot.py \\
        --input ${evaluation_files} \\
        --output landscape \\
        --verbose || true
    """
}

process SignalSupportDashboard {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/comparison", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(reason_tsvs)

    output:
    path "${test_name}_signal_support_dashboard.png", emit: dashboard_plot, optional: true

    script:
    """
    cage_tsvs=\$(ls *cage_peak_reasons* 2>/dev/null || true)
    drna_tsvs=\$(ls *drna_peak_reasons* 2>/dev/null || true)

    cage_args=""
    drna_args=""
    [ -n "\$cage_tsvs" ] && cage_args="--cage-tsvs \$cage_tsvs"
    [ -n "\$drna_tsvs" ] && drna_args="--drna-tsvs \$drna_tsvs"

    python ${projectDir}/bin/evaluation/signal_support_dashboard.py \\
        \$cage_args \\
        \$drna_args \\
        --output ${test_name}_signal_support_dashboard.png \\
        --title-prefix "${test_name}: " \\
        --verbose || true
    """
}

process PeakReasonHeatmap {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/ted_diagnostics/peak_reason_heatmaps", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(reason_tsvs)

    output:
    path "*_peak_reason_heatmap.png", optional: true

    script:
    """
    cage_tsvs=\$(ls *cage_peak_reasons* 2>/dev/null || true)
    drna_tsvs=\$(ls *drna_peak_reasons* 2>/dev/null || true)

    cage_args=""
    drna_args=""
    [ -n "\$cage_tsvs" ] && cage_args="--cage-tsvs \$cage_tsvs"
    [ -n "\$drna_tsvs" ] && drna_args="--drna-tsvs \$drna_tsvs"

    python ${projectDir}/bin/evaluation/peak_reason_heatmap.py \\
        \$cage_args \\
        \$drna_args \\
        --output-prefix ${test_name} \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// Isoform diversity parallel-coordinates plot
// Compares reference + FLAIR output GTFs across multi-region partitions.
// Skipped automatically when only one region is present.
// -----------------------------------------------------------------
process IsoformDiversity {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/isoform_structure", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(ref_gtf), val(regions), path(flair_gtfs), val(flair_labels)

    output:
    path "${test_name}_isoform_diversity.html", emit: diversity_html, optional: true
    path "${test_name}_isoform_diversity_bars.png", emit: diversity_bars, optional: true

    script:
    def region_args = regions.collect { "--regions ${it}" }.join(' ')
    // Build --gtfs and --labels arguments: reference first, then each FLAIR mode
    def all_gtfs = ["${ref_gtf}"] + flair_gtfs.collect { it.name }
    def all_labels = ["Reference"] + flair_labels
    def gtf_args = all_gtfs.collect { "\"${it}\"" }.join(' ')
    def label_args = all_labels.collect { "\"${it}\"" }.join(' ')
    """
    python ${projectDir}/bin/evaluation/isoform_diversity_plot.py \\
        --gtfs ${gtf_args} \\
        --labels ${label_args} \\
        --regions ${regions.join(' ')} \\
        --output ${test_name}_isoform_diversity.png \\
        --title-prefix "${test_name}: " \\
        --bars --verbose || true
    """
}

// -----------------------------------------------------------------
// Rescue impact dashboard
// -----------------------------------------------------------------
// TP-overlap analysis: pairwise TP set comparisons vs baseline
// Requires per-peak reason TSVs (CAGE + dRNA) with label:path pairs.
// -----------------------------------------------------------------
process TpOverlapPlot {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/end_accuracy/tp_overlap", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(reason_tsvs)

    output:
    path "tp_overlap/*.png", optional: true
    path "tp_overlap/*.tsv", optional: true

    script:
    // Build label:path pairs from staged peak-reason TSV filenames.
    // Filenames follow: {mode}_transcriptome_{cage|drna}_peak_reasons.tsv
    """
    cage_args=""
    drna_args=""
    for f in *_cage_peak_reasons.tsv; do
        [ -f "\$f" ] || continue
        label=\$(echo "\$f" | sed 's/_transcriptome_cage_peak_reasons\\.tsv//' | sed 's/_cage_peak_reasons\\.tsv//')
        cage_args="\$cage_args \$label:\$f"
    done
    for f in *_drna_peak_reasons.tsv; do
        [ -f "\$f" ] || continue
        label=\$(echo "\$f" | sed 's/_transcriptome_drna_peak_reasons\\.tsv//' | sed 's/_drna_peak_reasons\\.tsv//')
        drna_args="\$drna_args \$label:\$f"
    done

    [ -n "\$cage_args" ] && cage_args="--cage \$cage_args"
    [ -n "\$drna_args" ] && drna_args="--drna \$drna_args"

    python ${projectDir}/bin/evaluation/tp_overlap_plot.py \\
        \$cage_args \\
        \$drna_args \\
        --output tp_overlap \\
        --baseline-label baseline \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// Isoforms-per-gene frequency histogram
// Requires BED12 isoform files with label:path pairs.
// -----------------------------------------------------------------
process IsoformsPerGeneHist {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/isoform_structure", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files)

    output:
    path "isoforms_per_gene/*.png", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/isoforms_per_gene_hist.py \\
        --bed ${bed_args.join(' ')} \\
        --output isoforms_per_gene \\
        --box \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// Jaccard heatmaps: splice-junction + transcript-end Jaccard
// Requires BED12 isoform files with label:path pairs.
// -----------------------------------------------------------------
process JaccardHeatmapPlot {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/isoform_structure", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files)

    output:
    path "jaccard_heatmaps/*.png", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/jaccard_heatmap_plot.py \\
        --bed ${bed_args.join(' ')} \\
        --output jaccard_heatmaps \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// Gene-level isoform variation proportions.
// Classifies genes by variation type (alt ends, alt splicing, both).
// Requires BED12 isoform files with label:path pairs.
// -----------------------------------------------------------------
process GeneVariationPlot {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/isoform_structure", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files)

    output:
    path "gene_variation/*.png", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/gene_variation_plot.py \\
        --bed ${bed_args.join(' ')} \\
        --output gene_variation \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// SJC end distance plots: pairwise distances within splice junction
// chain groups (alternative promoters / polyadenylation sites).
// -----------------------------------------------------------------
process SjcEndDistancePlot {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/isoform_structure", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files)

    output:
    path "sjc_end_distances/*.png", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/sjc_end_distance_plot.py \\
        --bed ${bed_args.join(' ')} \\
        --output sjc_end_distances || true
    """
}

// -----------------------------------------------------------------
// Total isoforms horizontal bar chart.
// Requires evaluation TSV files.
// -----------------------------------------------------------------
process TotalIsoformsPlot {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/isoform_structure", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(eval_tsvs)

    output:
    path "total_isoforms/*.png", optional: true

    script:
    """
    python ${projectDir}/bin/evaluation/total_isoforms_plot.py \\
        --eval ${eval_tsvs} \\
        --output total_isoforms \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// End-signal density scatter (per-isoform KDE-coloured TSS vs TTS).
// Requires BED12 files + signal bedGraph tracks.
// -----------------------------------------------------------------
process EndSignalScatterPlot {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/signal", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files),
          val(cage_signal_plus), val(cage_signal_minus),
          val(drna_signal_plus), val(drna_signal_minus)

    output:
    path "end_signal_scatter/*.png", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/end_signal_scatter_plot.py \\
        --bed ${bed_args.join(' ')} \\
        --cage-plus ${cage_signal_plus} --cage-minus ${cage_signal_minus} \\
        --qs-plus ${drna_signal_plus} --qs-minus ${drna_signal_minus} \\
        --output end_signal_scatter \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// TED score vs boundary signal: scatter + correlation bar chart.
// Compares per-isoform TED scores (depth, model, annot, reality)
// to orthogonal CAGE/dRNA signal.  Non-TED BEDs auto-skipped.
// -----------------------------------------------------------------
process TedScoreVsSignal {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/signal", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files),
          val(cage_signal_plus), val(cage_signal_minus),
          val(drna_signal_plus), val(drna_signal_minus)

    output:
    path "ted_score_signal/*.png", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/ted_score_vs_signal.py \\
        --bed ${bed_args.join(' ')} \\
        --cage-plus ${cage_signal_plus} --cage-minus ${cage_signal_minus} \\
        --qs-plus ${drna_signal_plus} --qs-minus ${drna_signal_minus} \\
        --output ted_score_signal \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// Cumulative orthogonal signal curves.
// Requires BED12 files + signal bedGraph tracks.
// -----------------------------------------------------------------
process CumulativeSignalPlot {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/signal", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files), path(read_map_files),
          val(cage_signal_plus), val(cage_signal_minus),
          val(drna_signal_plus), val(drna_signal_minus)

    output:
    path "cumulative_signal/*.png", optional: true

    script:
    def bed_args = []
    def map_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
        map_args << "${bed_labels[i]}:${read_map_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/cumulative_signal_plot.py \\
        --bed ${bed_args.join(' ')} \\
        --read-map ${map_args.join(' ')} \\
        --cage-plus ${cage_signal_plus} --cage-minus ${cage_signal_minus} \\
        --qs-plus ${drna_signal_plus} --qs-minus ${drna_signal_minus} \\
        --output cumulative_signal \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// SJC alt-end analysis: TP/FP breakdown and boundary signal for
// alternative ends within splice-junction-chain groups.
// Requires BED12 files + read maps + peaks + signal bedGraph tracks.
// -----------------------------------------------------------------
process SjcAltEndAnalysis {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/ted_diagnostics", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files),
          val(cage_peaks), val(drna_peaks),
          val(cage_signal_plus), val(cage_signal_minus),
          val(drna_signal_plus), val(drna_signal_minus)

    output:
    path "sjc_alt_end/*.png", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/sjc_alt_end_analysis.py \\
        --bed ${bed_args.join(' ')} \\
        --cage-peaks ${cage_peaks} --qs-peaks ${drna_peaks} \\
        --cage-plus ${cage_signal_plus} --cage-minus ${cage_signal_minus} \\
        --qs-plus ${drna_signal_plus} --qs-minus ${drna_signal_minus} \\
        --output sjc_alt_end \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// TED component diagnostic: ROC curves, violins, joint heatmaps,
// signal×score scatter, weight sweep, and marginal model value.
// Requires BED12+TED files + peaks + signal tracks.
// -----------------------------------------------------------------
process TedComponentDiagnostic {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/ted_diagnostics", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files),
          val(cage_peaks), val(drna_peaks),
          val(cage_signal_plus), val(cage_signal_minus),
          val(drna_signal_plus), val(drna_signal_minus)

    output:
    path "ted_component_diagnostic/*.{png,tsv}", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    """
    # v2: add joint threshold sweep panel
    python ${projectDir}/bin/evaluation/ted_component_diagnostic.py \\
        --bed ${bed_args.join(' ')} \\
        --cage-peaks ${cage_peaks} --qs-peaks ${drna_peaks} \\
        --cage-plus ${cage_signal_plus} --cage-minus ${cage_signal_minus} \\
        --qs-plus ${drna_signal_plus} --qs-minus ${drna_signal_minus} \\
        --output ted_component_diagnostic \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// Signal-stratified peak recovery + ROC curves (per sample).
// Outputs (per end type, 5prime and 3prime):
//   signal_recall_{end}.png              raw signal x-axis
//   signal_recall_{end}_norm.png         normalised [0,1] x-axis
//   signal_roc_{end}.png                 TPR-FPR ROC curve
//   signal_recall_5v3_overlay.png        5' vs 3' per-tool comparison
// -----------------------------------------------------------------
process PeakRocCurves {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/end_accuracy/peak_roc", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), path(reason_tsvs)

    output:
    path "signal_recall_*.png", optional: true
    path "signal_roc_*.png",    optional: true

    script:
    """
    python ${projectDir}/bin/evaluation/peak_roc_curves.py \\
        --input ${reason_tsvs} \\
        --output . \\
        --dataset "${test_name}" \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// Read end-signal scatter: per-sample KDE-coloured scatter of TSS
// (CAGE) vs TTS (dRNA) signal at raw aligned read ends.
// Shows what the input data looks like before assembly.
// -----------------------------------------------------------------
process ReadEndSignalScatter {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/signal", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(read_labels), path(read_files),
          val(cage_signal_plus), val(cage_signal_minus),
          val(drna_signal_plus), val(drna_signal_minus)

    output:
    path "read_end_signal/*.png", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < read_labels.size(); i++) {
        bed_args << "${read_labels[i]}:${read_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/read_end_signal_scatter.py \\
        --bed ${bed_args.join(' ')} \\
        --cage-plus ${cage_signal_plus} --cage-minus ${cage_signal_minus} \\
        --qs-plus ${drna_signal_plus} --qs-minus ${drna_signal_minus} \\
        --output read_end_signal \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// Internal priming analysis: A-content at TTS, cross-tool comparison,
// and proximal APA distance scatter.
// Requires BED12 isoform files + read maps + genome FASTA + GTF.
// -----------------------------------------------------------------
process InternalPrimingAnalysis {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/isoform_structure/internal_priming", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files), path(read_map_files),
          path(genome), path(gtf)

    output:
    path "*.png", optional: true

    script:
    def bed_args = []
    def map_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
        map_args << "${bed_labels[i]}:${read_map_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/internal_priming.py \\
        --bed ${bed_args.join(' ')} \\
        --read-map ${map_args.join(' ')} \\
        --genome ${genome} \\
        --gtf ${gtf} \\
        --output . \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// Signal distribution at isoform ends across assembly modes.
// Asks: do different modes place isoform ends on stronger signal peaks?
// Requires BED12 files + CAGE/dRNA signal bedGraph tracks (no read maps needed).
// -----------------------------------------------------------------
process SignalReadSupport {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/signal/signal_distribution", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files), path(read_map_files),
          val(cage_signal_plus), val(cage_signal_minus),
          val(drna_signal_plus), val(drna_signal_minus)

    output:
    path "*.png", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/signal_read_support.py \\
        --bed ${bed_args.join(' ')} \\
        --cage-plus ${cage_signal_plus} --cage-minus ${cage_signal_minus} \\
        --qs-plus ${drna_signal_plus} --qs-minus ${drna_signal_minus} \\
        --output . \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// P/R × boundary-signal balance: three-axis scatter, Pareto frontier,
// and dead-zone bar + F1 overlay.
// Requires combined evaluation TSV(s) + BED12 files + signal tracks.
// -----------------------------------------------------------------
process PrSignalBalance {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/signal/pr_signal_balance", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files), path(read_map_files),
          val(cage_signal_plus), val(cage_signal_minus),
          val(drna_signal_plus), val(drna_signal_minus),
          path(combined_tsv)

    output:
    path "*.png", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/pr_signal_balance.py \\
        --tsv ${test_name}:${combined_tsv} \\
        --bed ${bed_args.join(' ')} \\
        --cage-plus ${cage_signal_plus} --cage-minus ${cage_signal_minus} \\
        --qs-plus ${drna_signal_plus} --qs-minus ${drna_signal_minus} \\
        --output . \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// SQANTI-stratified end precision: per-category (FSM/ISM/NIC/NNC) TSS
// and TTS precision across modes. Classifies isoforms from BED12 files
// against the reference GTF and computes end precision within each category.
// -----------------------------------------------------------------
process SqantiPrecision {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/end_accuracy/sqanti_precision", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files), path(gtf), path(category_tsvs),
          path(cage_peaks), path(drna_peaks)

    output:
    path "*.png", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    def cat_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        cat_args << "${bed_labels[i]}:${category_tsvs[i]}"
    }
    def cage_arg = cage_peaks.name != 'NO_CAGE' ? "--cage-peaks ${cage_peaks}" : ""
    def drna_arg = drna_peaks.name != 'NO_DRNA' ? "--drna-peaks ${drna_peaks}" : ""
    """
    # v8: add SEM/SEN single-exon categories
    python ${projectDir}/bin/evaluation/sqanti_precision.py \\
        --bed ${bed_args.join(' ')} \\
        --categories-tsv ${cat_args.join(' ')} \\
        --gtf ${gtf} \\
        ${cage_arg} \\
        ${drna_arg} \\
        --output . \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// Depth calibration: TED acceptance rate, depth-score violin, and
// ECDF of n_reads for accepted vs rejected clusters.
// Requires TED log TSV files (produced when --ted_log is passed).
// Non-TED modes produce empty/absent logs and are skipped automatically.
// -----------------------------------------------------------------
process DepthCalibration {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/ted_diagnostics/depth_calibration", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(ted_log_labels), path(ted_log_files)

    output:
    path "*.png", optional: true

    script:
    def log_args = []
    for (int i = 0; i < ted_log_labels.size(); i++) {
        log_args << "${ted_log_labels[i]}:${ted_log_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/depth_calibration.py \\
        --ted-log ${log_args.join(' ')} \\
        --output . \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// Cluster spread plots (D2): TSS/TTS IQR by pass/reject, spread vs
// signal, spread vs n_reads.  Requires TED log files + signal tracks.
// -----------------------------------------------------------------
process ClusterSpreadPlots {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/ted_diagnostics/cluster_spread", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(ted_log_labels), path(ted_log_files),
          val(cage_signal_plus), val(cage_signal_minus),
          val(drna_signal_plus), val(drna_signal_minus)

    output:
    path "*.png", optional: true

    script:
    def log_args = []
    for (int i = 0; i < ted_log_labels.size(); i++) {
        log_args << "${ted_log_labels[i]}:${ted_log_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/cluster_spread_plots.py \\
        --ted-log ${log_args.join(' ')} \\
        --cage-plus ${cage_signal_plus} --cage-minus ${cage_signal_minus} \\
        --qs-plus ${drna_signal_plus} --qs-minus ${drna_signal_minus} \\
        --output . \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// End-signal meta-profile: average CAGE/dRNA profile centred on
// called TSS and TTS positions.  All methods overlaid on shared axes.
// Requires BED12/GTF isoform files + CAGE/dRNA signal bedGraph tracks.
// -----------------------------------------------------------------
process EndSignalMetaplot {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/signal/end_signal_profiles", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files),
          val(cage_signal_plus), val(cage_signal_minus),
          val(drna_signal_plus), val(drna_signal_minus)

    output:
    path "metaplot_*.png",         optional: true
    path "metaplot_metrics.tsv",   optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/end_signal_metaplot.py \\
        --bed ${bed_args.join(' ')} \\
        --cage-plus ${cage_signal_plus} --cage-minus ${cage_signal_minus} \\
        --qs-plus ${drna_signal_plus} --qs-minus ${drna_signal_minus} \\
        --output . \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// End-signal heatmap (smarca4-style): per-isoform signal density
// rows sorted by read support, columns = bp offset from called end.
// Requires BED12/GTF isoform files + read maps + CAGE/dRNA bedGraphs.
// -----------------------------------------------------------------
process EndSignalHeatmap {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/signal/end_signal_profiles", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files), path(read_map_files),
          val(cage_signal_plus), val(cage_signal_minus),
          val(drna_signal_plus), val(drna_signal_minus)

    output:
    path "heatmap_*.png", optional: true

    script:
    def bed_args = []
    def map_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
        map_args << "${bed_labels[i]}:${read_map_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/end_signal_heatmap.py \\
        --bed ${bed_args.join(' ')} \\
        --read-map ${map_args.join(' ')} \\
        --cage-plus ${cage_signal_plus} --cage-minus ${cage_signal_minus} \\
        --qs-plus ${drna_signal_plus} --qs-minus ${drna_signal_minus} \\
        --output . \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// Read-end heatmap: where assigned reads' ends land relative to the
// called isoform end.  Long-read analogue of EndSignalHeatmap with
// no orthogonal signal — pure read distribution.
// -----------------------------------------------------------------
process ReadEndHeatmap {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/signal/end_signal_profiles", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files), path(read_map_files),
          path(reads_bed)

    output:
    path "read_end_heatmap_*.png", optional: true

    script:
    def bed_args = []
    def map_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
        map_args << "${bed_labels[i]}:${read_map_files[i]}"
    }
    """
    python ${projectDir}/bin/evaluation/read_end_heatmap.py \\
        --bed ${bed_args.join(' ')} \\
        --read-map ${map_args.join(' ')} \\
        --reads-bed ${reads_bed} \\
        --output . \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// TED log analysis (E1+E2): TSS vs TTS spread violin + summary
// scatter, CAGE peak width vs cluster IQR, threshold margin violin.
// Requires TED log files + CAGE peaks BED.
// -----------------------------------------------------------------
process TedLogAnalysis {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/ted_diagnostics/ted_log_analysis", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(ted_log_labels), path(ted_log_files),
          val(cage_peaks)

    output:
    path "*.png", optional: true

    script:
    def log_args = []
    for (int i = 0; i < ted_log_labels.size(); i++) {
        log_args << "${ted_log_labels[i]}:${ted_log_files[i]}"
    }
    def cage_arg = (cage_peaks && cage_peaks != 'NO_CAGE') ? "--cage-peaks ${cage_peaks}" : ""
    """
    python ${projectDir}/bin/evaluation/ted_log_analysis.py \\
        --ted-log ${log_args.join(' ')} \\
        ${cage_arg} \\
        --output . \\
        --verbose || true
    """
}

// -----------------------------------------------------------------
// Isoform end AUC-ROC: treats each isoform as a binary prediction
// (TP = end within --window bp of a reference peak), uses orthogonal
// signal at the called end as the classifier score, and plots TPR vs
// FPR curves for all methods overlaid.  JC-deduplication is applied
// so over-segmented assemblers don't inflate their AUC.
// Outputs:
//   roc_5prime.png        TSS ROC curves (CAGE signal)
//   roc_3prime.png        TTS ROC curves (dRNA signal)
//   roc_combined.png      Side-by-side 5′/3′
//   roc_auc_summary.tsv   mode, end, AUC, n_isoforms, n_tp, n_fp
// -----------------------------------------------------------------
process IsoformEndRoc {
    publishDir "${params.outdir}/evaluations/per_sample/${test_name}/summary/end_accuracy/roc", mode: 'copy'
    tag "${test_name}"
    errorStrategy 'ignore'

    input:
    tuple val(test_name), val(bed_labels), path(bed_files),
          val(cage_peaks), val(drna_peaks),
          val(cage_signal_plus), val(cage_signal_minus),
          val(drna_signal_plus), val(drna_signal_minus)

    output:
    path "roc_*.png",           optional: true
    path "roc_auc_summary.tsv", optional: true

    script:
    def bed_args = []
    for (int i = 0; i < bed_labels.size(); i++) {
        bed_args << "${bed_labels[i]}:${bed_files[i]}"
    }
    def cage_pk_arg = (cage_peaks && cage_peaks != 'NO_CAGE') ? "--cage-peaks ${cage_peaks}" : ""
    def drna_pk_arg = (drna_peaks && drna_peaks != 'NO_DRNA') ? "--drna-peaks ${drna_peaks}" : ""
    """
    python ${projectDir}/bin/evaluation/isoform_end_roc.py \\
        --bed ${bed_args.join(' ')} \\
        ${cage_pk_arg} \\
        ${drna_pk_arg} \\
        --cage-plus ${cage_signal_plus} --cage-minus ${cage_signal_minus} \\
        --qs-plus ${drna_signal_plus}   --qs-minus ${drna_signal_minus} \\
        --output . \\
        --verbose || true
    """
}