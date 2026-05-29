// =============================================================================
// Subworkflow: EVALUATE  (STAGE 2 core)
// =============================================================================
// Runs PrepareReferencePeaks + Evaluation + TedEndPrecision +
// FirstpassComparison against FIXED assemblies. Driven by NORMALIZED channels
// that BOTH the in-process monolith path (ASSEMBLE.out) and the split path
// (parsed manifest) can produce. No assembler / align / partition runs here.
//
// Inputs (all normalized — identical shapes from either driver):
//   isoform_tuples    — [test_name, dataset_name, align_mode, partition_mode,
//                        transcriptome_mode, isoforms_bed, isoforms_gtf,
//                        isoform_read_map, ted_log, isoform_counts]
//   eval_context_seed — [test_name, dataset_name, align_mode, partition_mode,
//                        bam, bai, reads_bed, genome(val), gtf, cage_peaks,
//                        drna_peaks, partition_args]
//   flair_firstpass   — [test_name, dataset_name, align_mode, partition_mode,
//                        partition_args, transcriptome_mode, firstpass_bed]
//   dataset_signal_ch — [test_name, library_type, cage_signal_plus,
//                        cage_signal_minus, drna_signal_plus, drna_signal_minus]
//   placeholders      — val map (NO_CAGE, NO_DRNA, ...)
//
// Emits the same surface ASSEMBLE_AND_EVAL did, so SUMMARY_AND_VIZ is unchanged.
// =============================================================================

include { PrepareReferencePeaks } from '../modules/evaluation/prepare_ref_peaks/main'
include { Evaluation            } from '../modules/evaluation/main'
include { TedEndPrecision       } from '../modules/evaluation/main'
include { FirstpassComparison   } from '../modules/evaluation/main'

workflow EVALUATE {

    take:
        isoform_tuples     // see header
        eval_context_seed  // see header
        flair_firstpass    // see header
        dataset_signal_ch  // see header
        placeholders       // val map

    main:
        // ------------------------------------------------------------------
        // Reference peaks (TSS/TTS) — derived from the partition annotation GTF.
        // Keyed on the 4-tuple; the GTF comes from eval_context_seed.
        // ------------------------------------------------------------------
        ref_peak_inputs = eval_context_seed.map {
            test_name, dataset_name, align_mode, partition_mode,
            bam, bai, reads_bed, genome, gtf, cage_peaks, drna_peaks, partition_args ->
            [test_name, dataset_name, align_mode, partition_mode, gtf]
        }
        PrepareReferencePeaks(ref_peak_inputs)

        // ------------------------------------------------------------------
        // Full eval context = partition outputs + ref peaks + signals.
        // Shape matches the old ASSEMBLE_AND_EVAL eval_context_ch exactly so
        // Evaluation's input tuple is unchanged.
        // ------------------------------------------------------------------
        eval_context_ch = eval_context_seed
            .combine(PrepareReferencePeaks.out.ref_peaks, by: [0, 1, 2, 3])
            .combine(dataset_signal_ch, by: [0])
            .map { test_name, dataset_name, align_mode, partition_mode,
                   bam, bai, reads_bed, genome, gtf, cage_peaks, drna_peaks, partition_args,
                   ref_tss, ref_tts,
                   library_type,
                   cage_signal_plus, cage_signal_minus, drna_signal_plus, drna_signal_minus ->
                [test_name, dataset_name, align_mode, partition_mode,
                 bam, bai, reads_bed, genome, gtf,
                 cage_peaks ?: placeholders.NO_CAGE, drna_peaks ?: placeholders.NO_DRNA,
                 ref_tss, ref_tts,
                 library_type,
                 cage_signal_plus, cage_signal_minus, drna_signal_plus, drna_signal_minus]
            }

        // partition_args recovered separately (keyed) for TedEndPrecision/Firstpass
        partition_args_ch = eval_context_seed.map {
            test_name, dataset_name, align_mode, partition_mode,
            bam, bai, reads_bed, genome, gtf, cage_peaks, drna_peaks, partition_args ->
            [test_name, dataset_name, align_mode, partition_mode, partition_args]
        }

        // ------------------------------------------------------------------
        // Evaluation
        // ------------------------------------------------------------------
        all_eval_inputs = isoform_tuples
            .combine(eval_context_ch, by: [0, 1, 2, 3])

        Evaluation(all_eval_inputs)

        // ------------------------------------------------------------------
        // TED End Precision
        // ------------------------------------------------------------------
        ted_precision_inputs = all_eval_inputs
            .combine(partition_args_ch, by: [0, 1, 2, 3])
            .filter { items ->
                def cage = items[15]
                def drna = items[16]
                def has_cage = cage.name != 'NO_CAGE' && cage.size() > 0
                def has_drna = drna.name != 'NO_DRNA' && drna.size() > 0
                has_cage || has_drna
            }
            .map { test_name, dataset_name, align_mode, partition_mode, transcriptome_mode,
                   isoforms_bed, isoforms_gtf, isoform_read_map, ted_log, isoform_counts,
                   bam, bai, reads_bed, genome, gtf,
                   cage_peaks, drna_peaks, ref_tss, ref_tts, library_type,
                   cage_signal_plus, cage_signal_minus, drna_signal_plus, drna_signal_minus,
                   partition_args ->
                [test_name, dataset_name, transcriptome_mode,
                 isoforms_bed, isoforms_gtf, gtf, isoform_counts,
                 cage_peaks, drna_peaks, partition_args]
            }
        TedEndPrecision(ted_precision_inputs)

        // ------------------------------------------------------------------
        // Firstpass vs Final comparison (FLAIR only)
        // ------------------------------------------------------------------
        // Build a final-BED channel keyed by 5-tuple from isoform_tuples (FLAIR
        // rows carry the real BED; others carry NO_ISOFORMS_BED and are filtered).
        final_bed_keyed = isoform_tuples.map {
            test_name, dataset_name, align_mode, partition_mode, transcriptome_mode,
            isoforms_bed, isoforms_gtf, isoform_read_map, ted_log, isoform_counts ->
            [[test_name, dataset_name, align_mode, partition_mode, transcriptome_mode], isoforms_bed]
        }

        firstpass_keyed = flair_firstpass.map {
            test_name, dataset_name, align_mode, partition_mode, partition_args, transcriptome_mode, firstpass_bed ->
            [[test_name, dataset_name, align_mode, partition_mode, transcriptome_mode], firstpass_bed]
        }

        firstpass_comparison_inputs = firstpass_keyed
            .join(final_bed_keyed)
            .map { key5, firstpass_bed, final_bed ->
                [key5[0], key5[1], key5[2], key5[3], key5[4], firstpass_bed, final_bed]
            }
            .combine(eval_context_ch, by: [0, 1, 2, 3])
            .combine(partition_args_ch, by: [0, 1, 2, 3])
            .filter { items ->
                // items: 0-3 keys, 4=mode, 5=firstpass, 6=final_bed,
                //   7=bam,8=bai,9=reads_bed,10=genome,11=gtf,12=cage,13=drna,14=ref_tss,...
                def firstpass = items[5]
                firstpass != null && firstpass.name != 'NO_FIRSTPASS' &&
                    !firstpass.name.startsWith('NO_') && firstpass.size() > 0
            }
            .filter { items ->
                def cage = items[12]
                def drna = items[13]
                def has_cage = cage.name != 'NO_CAGE' && cage.size() > 0
                def has_drna = drna.name != 'NO_DRNA' && drna.size() > 0
                has_cage || has_drna
            }
            .map { test_name, dataset_name, align_mode, partition_mode, transcriptome_mode,
                   firstpass_bed, final_bed,
                   bam, bai, reads_bed, genome, gtf, cage_peaks, drna_peaks,
                   ref_tss, ref_tts, library_type,
                   cage_signal_plus, cage_signal_minus, drna_signal_plus, drna_signal_minus,
                   partition_args ->
                [test_name, dataset_name, transcriptome_mode,
                 firstpass_bed, final_bed, gtf, cage_peaks, drna_peaks, partition_args]
            }
        if (params.run_firstpass_comparison) {
            FirstpassComparison(firstpass_comparison_inputs)
        }

    emit:
        evaluation_results      = Evaluation.out.evaluation_results
        isoform_categories      = Evaluation.out.isoform_categories
        cage_peak_reason_tsvs   = Evaluation.out.cage_peak_reason_tsvs
        drna_peak_reason_tsvs   = Evaluation.out.drna_peak_reason_tsvs
        ted_precision_metrics   = TedEndPrecision.out.metrics
        ted_gtf_precision       = TedEndPrecision.out.gtf_metrics
        all_eval_inputs         = all_eval_inputs
}
