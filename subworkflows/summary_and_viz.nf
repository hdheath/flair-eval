// =============================================================================
// Subworkflow: SUMMARY_AND_VIZ
// =============================================================================
// Runs all downstream summary and comparison analyses after evaluation:
//   - SummaryPlots: precision/recall, concordance, signal support comparisons
//   - PeakReasonHeatmap: cross-mode CAGE/dRNA peak reason breakdowns
//   - TpOverlapPlot: pairwise TP set comparisons vs baseline
//   - IsoformsPerGeneHist: isoforms-per-gene frequency histogram + box plot
//   - JaccardHeatmapPlot: splice-junction + transcript-end Jaccard heatmaps
//   - TotalIsoformsPlot: total isoform count per method
//   - EndSignalScatterPlot: per-isoform KDE-coloured TSS vs TTS scatter
//   - CumulativeSignalPlot: cumulative orthogonal signal curves
//   - ToolDivergence: Jaccard agreement + motif collapse across assemblers
//   - UTRFeatures: 5'UTR length, GC content, uORF comparison
//
// Inputs:
//   evaluation_results      — Evaluation.out.evaluation_results
//   cage_peak_reason_tsvs   — Evaluation.out.cage_peak_reason_tsvs
//   drna_peak_reason_tsvs — Evaluation.out.drna_peak_reason_tsvs
//   all_eval_inputs         — full joined eval channel (for isoform file extraction)
//   dataset_signal_ch       — per-dataset signal bedGraph paths
// =============================================================================

include { CombineEvaluationTSVs     } from '../modules/visualization/summary/main'
include { CombinePrecisionRecall                      } from '../modules/visualization/summary/main'
include { CombinePrecisionRecall as CombineGtfPrecisionRecall } from '../modules/visualization/summary/main'
include { PrecisionRecallPlot                         } from '../modules/visualization/summary/main'
include { SummaryPlots              } from '../modules/visualization/summary/main'
include { PeakReasonHeatmap         } from '../modules/visualization/summary/main'
include { SignalSupportDashboard    } from '../modules/visualization/summary/main'
include { IsoformDiversity          } from '../modules/visualization/summary/main'
include { TpOverlapPlot             } from '../modules/visualization/summary/main'
include { IsoformsPerGeneHist       } from '../modules/visualization/summary/main'
include { JaccardHeatmapPlot       } from '../modules/visualization/summary/main'
include { GeneVariationPlot        } from '../modules/visualization/summary/main'
include { SjcEndDistancePlot       } from '../modules/visualization/summary/main'
include { TotalIsoformsPlot        } from '../modules/visualization/summary/main'
include { EndSignalScatterPlot     } from '../modules/visualization/summary/main'
include { CumulativeSignalPlot     } from '../modules/visualization/summary/main'
include { TedScoreVsSignal         } from '../modules/visualization/summary/main'
include { SjcAltEndAnalysis        } from '../modules/visualization/summary/main'
include { TedComponentDiagnostic   } from '../modules/visualization/summary/main'
include { ReadEndSignalScatter     } from '../modules/visualization/summary/main'
include { PeakRocCurves            } from '../modules/visualization/summary/main'
include { InternalPrimingAnalysis  } from '../modules/visualization/summary/main'
include { SignalReadSupport        } from '../modules/visualization/summary/main'
include { PrSignalBalance          } from '../modules/visualization/summary/main'
include { SqantiPrecision          } from '../modules/visualization/summary/main'
include { DepthCalibration         } from '../modules/visualization/summary/main'
include { ClusterSpreadPlots       } from '../modules/visualization/summary/main'
include { TedLogAnalysis           } from '../modules/visualization/summary/main'
include { TedRejectionAnalysis     } from '../modules/visualization/summary/main'
include { TedConfusionMatrix       } from '../modules/visualization/summary/main'
include { EndSignalMetaplot        } from '../modules/visualization/summary/main'
include { EndSignalHeatmap         } from '../modules/visualization/summary/main'
include { ReadEndHeatmap           } from '../modules/visualization/summary/main'
include { IsoformEndRoc            } from '../modules/visualization/summary/main'

include { CrossSamplePrecisionRecall} from '../modules/visualization/cross_sample/main'
include { CrossSamplePeakRoc        } from '../modules/visualization/cross_sample/main'
include { CrossSampleConcordance    } from '../modules/visualization/cross_sample/main'
include { CrossSampleSignalSupport  } from '../modules/visualization/cross_sample/main'
include { CrossSampleLandscape      } from '../modules/visualization/cross_sample/main'

include { CrossSampleToolEndAccuracy} from '../modules/visualization/cross_sample/main'
workflow SUMMARY_AND_VIZ {

    take:
        evaluation_results         // tuple [test_name, dataset_name, align_mode, partition_mode, transcriptome_mode, eval_tsv]
        isoform_categories         // tuple [test_name, dataset_name, align_mode, partition_mode, transcriptome_mode, isoform_categories.tsv]
        cage_peak_reason_tsvs      // tuple [test_name, path]
        drna_peak_reason_tsvs      // tuple [test_name, path]
        all_eval_inputs            // full 23-field eval channel (for divergence + UTR)
        dataset_signal_ch          // tuple [test_name, library_type, cage_signal_plus, cage_signal_minus, qs_signal_plus, qs_signal_minus]
        ted_precision_metrics      // tuple [test_name, dataset_name, transcriptome_mode, precision_recall_summary.tsv, per_junction_chain.tsv]
        ted_gtf_precision          // tuple [test_name, dataset_name, transcriptome_mode, gtf_precision_recall_summary.tsv]
        flair_firstpass            // tuple [test_name, dataset_name, align_mode, partition_mode, partition_args, transcriptome_mode, firstpass_bed]

    main:


        // --- Summary plots: group evaluation TSVs by test_name ---
        all_evaluations = evaluation_results
            .map { test_name, dataset_name, align_mode, partition_mode, transcriptome_mode, eval_tsv ->
                [test_name, eval_tsv]
            }
            .groupTuple()

        SummaryPlots(all_evaluations)
        CombineEvaluationTSVs(all_evaluations)

        // --- Precision/Recall: combine per-method TSVs into one file, then plot ---
        all_pr_tsvs = ted_precision_metrics
            .map { test_name, dataset_name, transcriptome_mode, pr_tsv, jc_tsv ->
                [test_name, pr_tsv]
            }
            .groupTuple()

        CombinePrecisionRecall(all_pr_tsvs, "ortho")

        // Combine GTF-only precision TSVs in parallel (same CombinePrecisionRecall process)
        all_gtf_pr_tsvs = ted_gtf_precision
            .map { test_name, dataset_name, transcriptome_mode, gtf_tsv ->
                [test_name, gtf_tsv]
            }
            .groupTuple()

        CombineGtfPrecisionRecall(all_gtf_pr_tsvs, "gtf")

        // Join orthogonal + GTF combined TSVs by test_name for the scatter plot
        pr_plot_inputs = CombinePrecisionRecall.out.combined_pr
            .join(CombineGtfPrecisionRecall.out.combined_pr)
            .map { test_name, ortho_tsv, gtf_tsv -> [test_name, ortho_tsv, gtf_tsv] }

        PrecisionRecallPlot(pr_plot_inputs)

        // --- Peak-reason heatmap: mix CAGE + dRNA, group by test_name ---
        all_reason_tsvs = cage_peak_reason_tsvs
            .mix(drna_peak_reason_tsvs)
            .groupTuple()
            .map { test_name, tsvs -> [test_name, tsvs.flatten()] }

        PeakReasonHeatmap(all_reason_tsvs)

        // --- Signal-stratified peak recovery curves (per dataset) ---
        PeakRocCurves(all_reason_tsvs)

        // --- Signal-vs-support scatter dashboard: tiled per-mode scatters ---
        SignalSupportDashboard(all_reason_tsvs)

        // --- TP overlap analysis: pairwise TP set comparisons vs baseline ---
        TpOverlapPlot(all_reason_tsvs)

        // --- Isoforms-per-gene histogram: from BED12 or GTF isoform files ---
        // Pick whichever isoform file is available for each tool:
        //   FLAIR → BED12 (items[5]), Bambu/IsoQuant → GTF (items[6])
        isoforms_per_gene_ch = all_eval_inputs
            .map { items ->
                // items[0]=test_name, [4]=transcriptome_mode, [5]=isoforms_bed, [6]=isoforms_gtf
                def isoform_file = items[5].name.contains('NO_ISOFORMS_BED') ? items[6] : items[5]
                [items[0], items[4], isoform_file]
            }
            .filter { !it[2].name.contains('NO_ISOFORMS') }
            .groupTuple(by: [0])
            .map { test_name, labels, files ->
                [test_name, labels, files.flatten()]
            }

        IsoformsPerGeneHist(isoforms_per_gene_ch)

        // --- Jaccard heatmaps: splice-junction + transcript-end Jaccard ---
        // Reuses the same BED12 channel as isoforms_per_gene_ch
        JaccardHeatmapPlot(isoforms_per_gene_ch)

        // --- Gene variation proportions: alt ends vs alt splicing ---
        GeneVariationPlot(isoforms_per_gene_ch)

        // --- SJC end distance plots: pairwise distances within SJC groups ---
        // Labels include dataset_name for per-sample faceting.
        sjc_ch = all_eval_inputs
            .map { items ->
                def isoform_file = items[5].name.contains('NO_ISOFORMS_BED') ? items[6] : items[5]
                [items[0], "${items[1]}::${items[4]}", isoform_file]
            }
            .filter { !it[2].name.contains('NO_ISOFORMS') }
            .groupTuple(by: [0])
            .map { test_name, labels, files ->
                [test_name, labels, files.flatten()]
            }

        SjcEndDistancePlot(sjc_ch)

        // --- Total isoforms bar chart: from eval TSVs ---
        total_iso_ch = evaluation_results
            .map { test_name, dataset_name, align_mode, partition_mode, transcriptome_mode, eval_tsv ->
                [test_name, eval_tsv]
            }
            .groupTuple()

        TotalIsoformsPlot(total_iso_ch)

        // --- Signal-dependent plots: end-signal scatter + cumulative signal ---
        // Only run when signal bedGraph tracks are available.
        // Pick BED12 or GTF isoform file for each tool.
        signal_bed_ch = all_eval_inputs
            .map { items ->
                def isoform_file = items[5].name.contains('NO_ISOFORMS_BED') ? items[6] : items[5]
                [items[0], items[4], isoform_file]
            }
            .filter { !it[2].name.contains('NO_ISOFORMS') }
            .groupTuple(by: [0])
            .map { test_name, labels, files ->
                [test_name, labels, files.flatten()]
            }

        signal_plot_inputs = signal_bed_ch
            .join(
                dataset_signal_ch
                    .filter { it[2] && it[3] && it[4] && it[5] }  // all 4 signal tracks present
                    .map { test_name, library_type,
                           cage_signal_plus, cage_signal_minus,
                           drna_signal_plus, drna_signal_minus ->
                        [test_name,
                         cage_signal_plus, cage_signal_minus,
                         drna_signal_plus, drna_signal_minus]
                    }
            )
            .map { test_name, bed_labels, bed_files,
                   cage_signal_plus, cage_signal_minus,
                   drna_signal_plus, drna_signal_minus ->
                [test_name, bed_labels, bed_files,
                 cage_signal_plus, cage_signal_minus,
                 drna_signal_plus, drna_signal_minus]
            }

        EndSignalScatterPlot(signal_plot_inputs)

        // --- TED score vs signal: compare internal scoring to orthogonal signal ---
        TedScoreVsSignal(signal_plot_inputs)

        // --- Cumulative signal plot: uses read-map files for read counts ---
        // Build channel carrying isoform files + read-map files per method.
        cumulative_bed_ch = all_eval_inputs
            .map { items ->
                def isoform_file = items[5].name.contains('NO_ISOFORMS_BED') ? items[6] : items[5]
                [items[0], items[4], isoform_file, items[7]]
            }
            .filter { !it[2].name.contains('NO_ISOFORMS') }
            .groupTuple(by: [0])
            .map { test_name, labels, files, read_maps ->
                [test_name, labels, files.flatten(), read_maps.flatten()]
            }

        cumulative_signal_inputs = cumulative_bed_ch
            .join(
                dataset_signal_ch
                    .filter { it[2] && it[3] && it[4] && it[5] }
                    .map { test_name, library_type,
                           cage_signal_plus, cage_signal_minus,
                           drna_signal_plus, drna_signal_minus ->
                        [test_name,
                         cage_signal_plus, cage_signal_minus,
                         drna_signal_plus, drna_signal_minus]
                    }
            )
            .map { test_name, bed_labels, bed_files, read_maps,
                   cage_signal_plus, cage_signal_minus,
                   drna_signal_plus, drna_signal_minus ->
                [test_name, bed_labels, bed_files, read_maps,
                 cage_signal_plus, cage_signal_minus,
                 drna_signal_plus, drna_signal_minus]
            }

        CumulativeSignalPlot(cumulative_signal_inputs)

        // --- End-signal meta-profile: average CAGE/dRNA centred on called ends ---
        // Reuses signal_plot_inputs (bed files + signal tracks, no read maps needed).
        EndSignalMetaplot(signal_plot_inputs)

        // --- End-signal heatmap (smarca4-style): per-isoform rows, bp-offset columns ---
        // Reuses cumulative_signal_inputs (bed + read maps + signal tracks).
        EndSignalHeatmap(cumulative_signal_inputs)

        // --- Read-end heatmap: where assigned reads land relative to called ends ---
        // Same layout as EndSignalHeatmap but uses read-map + read_audit BED.
        // reads_bed is per-dataset (all modes share the same reads file), so we
        // pass a single deduplicated reads BED alongside the per-method isoform BEDs.
        // The script receives --reads-bed as a single shared file (not per-label).
        read_end_heatmap_inputs = cumulative_bed_ch
            .join(
                all_eval_inputs
                    .map { items -> [items[0], items[11]] }  // test_name, reads_bed
                    .unique { it[0] }                         // one reads_bed per test_name
            )
            .map { test_name, bed_labels, bed_files, read_maps, reads_bed ->
                [test_name, bed_labels, bed_files, read_maps, reads_bed]
            }
        ReadEndHeatmap(read_end_heatmap_inputs)

        // --- SJC alt-end analysis: TP/FP breakdown + boundary signal ---
        // Reuses cumulative_bed_ch (has read_maps) and adds peak files.
        // Peaks are per-dataset but identical within a test_name; take first.
        dataset_peaks_ch = all_eval_inputs
            .map { items -> [items[0], items[14], items[15]] }  // test_name, cage_peaks, qs_peaks
            .unique { it[0] }  // one per test_name

        sjc_alt_end_inputs = cumulative_bed_ch
            .join(dataset_peaks_ch)
            .join(
                dataset_signal_ch
                    .filter { it[2] && it[3] && it[4] && it[5] }
                    .map { test_name, library_type,
                           cage_signal_plus, cage_signal_minus,
                           drna_signal_plus, drna_signal_minus ->
                        [test_name,
                         cage_signal_plus, cage_signal_minus,
                         drna_signal_plus, drna_signal_minus]
                    }
            )
            .map { test_name, bed_labels, bed_files, read_maps,
                   cage_peaks, drna_peaks,
                   cage_signal_plus, cage_signal_minus,
                   drna_signal_plus, drna_signal_minus ->
                // read_maps not needed — analysis is purely peak-based
                [test_name, bed_labels, bed_files,
                 cage_peaks, drna_peaks,
                 cage_signal_plus, cage_signal_minus,
                 drna_signal_plus, drna_signal_minus]
            }

        SjcAltEndAnalysis(sjc_alt_end_inputs)

        // --- TED component diagnostic: ROC, violins, heatmaps, weight sweep ---
        // Uses signal_bed_ch (no read maps) + peaks + signal tracks.
        ted_diag_inputs = signal_bed_ch
            .join(dataset_peaks_ch)
            .join(
                dataset_signal_ch
                    .filter { it[2] && it[3] && it[4] && it[5] }
                    .map { test_name, library_type,
                           cage_signal_plus, cage_signal_minus,
                           drna_signal_plus, drna_signal_minus ->
                        [test_name,
                         cage_signal_plus, cage_signal_minus,
                         drna_signal_plus, drna_signal_minus]
                    }
            )
            .map { test_name, bed_labels, bed_files,
                   cage_peaks, drna_peaks,
                   cage_signal_plus, cage_signal_minus,
                   drna_signal_plus, drna_signal_minus ->
                [test_name, bed_labels, bed_files,
                 cage_peaks, drna_peaks,
                 cage_signal_plus, cage_signal_minus,
                 drna_signal_plus, drna_signal_minus]
            }

        TedComponentDiagnostic(ted_diag_inputs)

        // --- Isoform end AUC-ROC: signal-as-score, JC-dedup TP labels ---
        // Reuses ted_diag_inputs: same bed files + peaks + signal tracks.
        IsoformEndRoc(ted_diag_inputs)

        // --- Read end-signal scatter: raw read ends vs orthogonal signal ---
        // Deduplicate reads BED by (test_name, dataset_name) since all modes
        // share the same reads, then group by test_name.
        read_signal_ch = all_eval_inputs
            .map { items ->
                // items[1]=dataset_name, items[11]=reads_bed
                [items[0], items[1], items[11]]
            }
            .unique { it[0] + '::' + it[1] }
            .groupTuple(by: [0])
            .map { test_name, labels, files ->
                [test_name, labels, files.flatten()]
            }
            .join(
                dataset_signal_ch
                    .filter { it[2] && it[3] && it[4] && it[5] }
                    .map { test_name, library_type,
                           cage_signal_plus, cage_signal_minus,
                           drna_signal_plus, drna_signal_minus ->
                        [test_name,
                         cage_signal_plus, cage_signal_minus,
                         drna_signal_plus, drna_signal_minus]
                    }
            )
            .map { test_name, read_labels, read_files,
                   cage_signal_plus, cage_signal_minus,
                   drna_signal_plus, drna_signal_minus ->
                [test_name, read_labels, read_files,
                 cage_signal_plus, cage_signal_minus,
                 drna_signal_plus, drna_signal_minus]
            }

        ReadEndSignalScatter(read_signal_ch)

        // --- Internal priming analysis: A-content at TTS + cross-tool + APA scatter ---
        // Reuses cumulative_bed_ch (has read maps) + genome and GTF from eval inputs.
        // Genome and GTF are the same for all modes in a test_name; take first.
        genome_gtf_ch = all_eval_inputs
            .map { items ->
                // items[12]=genome, items[13]=gtf
                [items[0], items[12], items[13]]
            }
            .unique { it[0] }

        internal_priming_inputs = cumulative_bed_ch
            .join(genome_gtf_ch)
            .map { test_name, bed_labels, bed_files, read_maps, genome, gtf ->
                [test_name, bed_labels, bed_files, read_maps, genome, gtf]
            }

        InternalPrimingAnalysis(internal_priming_inputs)

        // --- Signal × read support: scatter + zero-signal fraction by bin ---
        // Reuses cumulative_signal_inputs (bed + read maps + signal tracks).
        SignalReadSupport(cumulative_signal_inputs)

        // --- Combined evaluation TSV channel: used by PrSignalBalance for P/R data ---
        sqanti_ch = CombineEvaluationTSVs.out.combined_tsv
            .map { combined_tsv ->
                // CombineEvaluationTSVs emits only the TSV (no test_name in the output).
                // Recover test_name from the filename: {test_name}_combined_evaluation.tsv
                def fname = combined_tsv.name
                def test_name = fname.replace('_combined_evaluation.tsv', '')
                [test_name, combined_tsv]
            }

        // --- SQANTI precision: per-category end precision (TSS/TTS within 50bp of annotated) ---
        // Pre-classified isoform categories (from flair_eval.py) are joined in so
        // sqanti_precision.py can skip the expensive re-parse + re-classification.
        // Build: [test_name, transcriptome_mode, bed_file, categories_tsv] per method,
        // then group by test_name into [test_name, labels, beds, cat_tsvs, gtf].
        sqanti_prec_ch = all_eval_inputs
            .map { items ->
                def isoform_file = items[5].name.contains('NO_ISOFORMS_BED') ? items[6] : items[5]
                // key = [test_name, dataset_name, align_mode, partition_mode, transcriptome_mode]
                [items[0], items[1], items[2], items[3], items[4], isoform_file, items[13]]
            }
            .filter { !it[5].name.contains('NO_ISOFORMS') }
            .join(
                isoform_categories.map { test_name, dataset_name, align_mode, partition_mode,
                                         transcriptome_mode, cat_tsv ->
                    [test_name, dataset_name, align_mode, partition_mode, transcriptome_mode, cat_tsv]
                },
                by: [0, 1, 2, 3, 4]
            )
            // Now: [test_name, dataset_name, align_mode, partition_mode, transcriptome_mode, bed, gtf, cat_tsv]
            .map { test_name, dataset_name, align_mode, partition_mode,
                   transcriptome_mode, bed, gtf, cat_tsv ->
                [test_name, transcriptome_mode, bed, gtf, cat_tsv]
            }
            .groupTuple(by: [0])
            .map { test_name, labels, beds, gtfs, cat_tsvs ->
                [test_name, labels, beds.flatten(), gtfs[0], cat_tsvs.flatten()]
            }
            .join(
                all_eval_inputs
                    .map { items -> [items[0], items[14], items[15]] }
                    .unique { it[0] }
            )
            .map { test_name, labels, beds, gtf, cat_tsvs, cage_peaks, drna_peaks ->
                [test_name, labels, beds, gtf, cat_tsvs, cage_peaks, drna_peaks]
            }

        SqantiPrecision(sqanti_prec_ch)

        // --- P/R × boundary-signal balance: 3-axis scatter, Pareto, dead-zone bar ---
        // Joins cumulative_signal_inputs (bed + signal tracks) with the combined
        // evaluation TSV (for P/R stats). read_map_files carried through but unused.
        pr_signal_inputs = cumulative_signal_inputs
            .join(sqanti_ch.map { test_name, tsv -> [test_name, tsv] })
            .map { test_name, bed_labels, bed_files, read_maps,
                   cage_plus, cage_minus, qs_plus, qs_minus,
                   combined_tsv ->
                [test_name, bed_labels, bed_files, read_maps,
                 cage_plus, cage_minus, qs_plus, qs_minus,
                 combined_tsv]
            }

        PrSignalBalance(pr_signal_inputs)

        // --- Depth calibration: TED log n_reads, acceptance rate, depth-score violin ---
        // Only include modes that produced a real TED log (non-placeholder, non-empty).
        depth_cal_ch = all_eval_inputs
            .map { items ->
                // items[4]=transcriptome_mode, items[8]=ted_log
                def ted_log = items[8]
                def has_log = ted_log.name != 'NO_TED_LOG' && !ted_log.name.startsWith('NO_') && ted_log.size() > 0
                has_log ? [items[0], items[4], ted_log] : null
            }
            .filter { it != null }
            .groupTuple(by: [0])
            .map { test_name, labels, logs ->
                [test_name, labels, logs.flatten()]
            }

        DepthCalibration(depth_cal_ch)

        // --- Cluster spread plots (D2a): IQR violin by pass/reject status ---
        ClusterSpreadPlots(depth_cal_ch)

        // --- TED log analysis (E1+E2): end spread + CAGE peak width vs cluster IQR ---
        // Reuses depth_cal_ch + CAGE peaks BED (from dataset_peaks_ch).
        ted_log_analysis_inputs = depth_cal_ch
            .join(dataset_peaks_ch)
            .map { test_name, ted_log_labels, ted_log_files, cage_peaks, _drna_peaks ->
                [test_name, ted_log_labels, ted_log_files, cage_peaks]
            }

        TedLogAnalysis(ted_log_analysis_inputs)

        // --- TED rejection analysis: drop_reason breakdown + spliced-length stratification ---
        // Joins TED log channel with firstpass BED (for spliced length) + cage/drna peaks.
        firstpass_keyed = flair_firstpass
            .map { test_name, dataset_name, align_mode, partition_mode, partition_args,
                   transcriptome_mode, firstpass_bed ->
                def has_bed = firstpass_bed != null && firstpass_bed.name != 'NO_FIRSTPASS' &&
                              !firstpass_bed.name.startsWith('NO_') && firstpass_bed.size() > 0
                has_bed ? [test_name, transcriptome_mode, firstpass_bed] : null
            }
            .filter { it != null }
            .groupTuple(by: [0])
            .map { test_name, labels, beds ->
                [test_name, labels, beds.flatten()]
            }

        ted_rejection_inputs = depth_cal_ch
            .join(firstpass_keyed, remainder: true)
            .join(dataset_peaks_ch)
            .map { test_name, ted_log_labels, ted_log_files,
                   fp_labels, fp_files,
                   cage_peaks, drna_peaks ->
                if (fp_labels == null) return null
                // Intersect labels: only configs that have both a ted log and a firstpass BED
                def common = ted_log_labels.intersect(fp_labels)
                if (common.isEmpty()) return null
                def log_map = [ted_log_labels, ted_log_files].transpose().collectEntries()
                def bed_map = [fp_labels, fp_files].transpose().collectEntries()
                def paired_logs = common.collect { log_map[it] }
                def paired_beds = common.collect { bed_map[it] }
                [test_name, common, paired_logs, paired_beds, cage_peaks, drna_peaks]
            }
            .filter { it != null }

        TedRejectionAnalysis(ted_rejection_inputs)

        // --- TED confusion matrix: pass/reject × joint TP/FP per config ---
        // Uses TED log + orthogonal CAGE/dRNA peaks.  No firstpass BED needed
        // (we score the TED log positions directly).
        ted_confusion_inputs = depth_cal_ch
            .join(dataset_peaks_ch)
            .map { test_name, ted_log_labels, ted_log_files, cage_peaks, drna_peaks ->
                [test_name, ted_log_labels, ted_log_files, cage_peaks, drna_peaks]
            }

        TedConfusionMatrix(ted_confusion_inputs)

        // --- Isoform diversity: parallel-coordinates across multi-region partitions ---
        // Extract regions from partition_mode args.  Only runs when ≥2 regions detected.
        diversity_ch = all_eval_inputs
            .map { test_name, dataset_name, align_mode, partition_mode, transcriptome_mode,
                   isoforms_bed, isoforms_gtf, isoform_read_map, ted_log,
                   bam, bai, reads_bed, genome, gtf,
                   cage_peaks, drna_peaks, ref_tss, ref_tts,
                   library_type,
                   cage_signal_plus, cage_signal_minus, drna_signal_plus, drna_signal_minus ->
                def iso_gtf = isoforms_gtf.name.contains('NO_ISOFORMS_GTF') ? null : isoforms_gtf
                [test_name, gtf, partition_mode, iso_gtf, transcriptome_mode]
            }
            .filter { it[3] != null }  // skip entries with no output GTF
            .groupTuple(by: [0])
            .map { test_name, gtfs, partition_modes, flair_gtfs, tx_modes ->
                // Parse regions from partition_mode name (stored in params JSON)
                // Regions are extracted by reading the pipeline config
                def ref_gtf = gtfs[0]  // reference GTF is the same for all entries
                [test_name, ref_gtf, flair_gtfs.flatten(), tx_modes]
            }

        // Read regions from params JSON config and inject into the channel
        def jsonSlurperDiv = new groovy.json.JsonSlurper()
        def divConfig = jsonSlurperDiv.parse(new File(file(params.params_file).toString()))
        def partitionConfig = divConfig.partition ?: [:]
        def allRegions = []
        partitionConfig.each { mode_name, mode_args ->
            def m = (mode_args =~ /--region\s+(.+?)(?:\s+--|$)/)
            if (m.find()) {
                m.group(1).trim().split(/\s+/).each { r -> allRegions << r }
            }
        }
        allRegions = allRegions.unique()

        if (allRegions.size() >= 2) {
            diversity_input = diversity_ch.map { test_name, ref_gtf, flair_gtfs, tx_modes ->
                [test_name, ref_gtf, allRegions, flair_gtfs, tx_modes]
            }
            IsoformDiversity(diversity_input)
        }



        // -----------------------------------------------------------------
        // Cross-sample summary: combined precision/recall across ALL samples
        // -----------------------------------------------------------------
        cross_sample_evals = evaluation_results
            .map { test_name, dataset_name, align_mode, partition_mode, transcriptome_mode, eval_tsv ->
                [params.test_name, eval_tsv]
            }
            .groupTuple()

        // CrossSamplePrecisionRecall uses the per-sample combined P/R TSVs
        cross_sample_pr_tsvs = CombinePrecisionRecall.out.combined_pr
            .map { test_name, pr_tsv -> [params.test_name, pr_tsv] }
            .groupTuple()

        CrossSamplePrecisionRecall(cross_sample_pr_tsvs)
        CrossSampleConcordance(cross_sample_evals)
        CrossSampleSignalSupport(cross_sample_evals)
        CrossSampleLandscape(cross_sample_evals)
        CrossSampleToolEndAccuracy(cross_sample_evals)

        // --- Cross-sample signal-stratified peak recovery curves ---
        cross_sample_reasons = cage_peak_reason_tsvs
            .mix(drna_peak_reason_tsvs)
            .map { test_name, tsv -> [params.test_name, tsv] }
            .groupTuple()
            .map { test_name, tsvs -> [test_name, tsvs.flatten()] }

        CrossSamplePeakRoc(cross_sample_reasons)



    emit:
        precision_recall_plot = PrecisionRecallPlot.out.precision_recall_plot
        cross_sample_pr_plot  = CrossSamplePrecisionRecall.out.cross_sample_pr_plot
        combined_evaluation   = CombineEvaluationTSVs.out.combined_tsv
}
