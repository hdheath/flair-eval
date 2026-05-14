// =============================================================================
// Subworkflow: ASSEMBLE_AND_EVAL
// =============================================================================
// Runs all transcriptome assemblers (FLAIR, Bambu, IsoQuant, IsoSeq, FLAMES,
// StringTie2) on partitioned data, prepares reference peaks, and runs the
// unified Evaluation process.
//
// Inputs:
//   partitioned_ch    — FlairPartition.out.partitioned
//   datasets_ch       — master datasets channel (for mode maps + junction_tab)
//   dataset_signal_ch — per-dataset signal bedGraph paths
//   placeholders      — map of placeholder file paths (NO_ISOFORMS_BED, etc.)
//
// Emits:
//   evaluation_results     — tupled Evaluation TSVs
//   cage_peak_reason_tsvs  — per-peak CAGE reason TSVs
//   drna_peak_reason_tsvs — per-peak dRNA reason TSVs
//   all_eval_inputs        — full joined eval channel (for divergence/UTR downstream)
//   flair_transcriptome    — FlairTranscriptome.out.transcriptome (for PlotIsoforms)
// =============================================================================

include { FlairTranscriptome   } from '../modules/transcriptome/flair/main'
include { BambuAssembly        } from '../modules/transcriptome/bambu/main'
include { IsoQuantAssembly     } from '../modules/transcriptome/isoquant/main'
include { IsoSeqAssembly       } from '../modules/transcriptome/isoseq/main'
include { FlamesAssembly       } from '../modules/transcriptome/flames/main'
include { StringTie2Assembly    } from '../modules/transcriptome/stringtie2/main'
include { PrepareReferencePeaks} from '../modules/evaluation/prepare_ref_peaks/main'
include { Evaluation           } from '../modules/evaluation/main'
include { TedEndPrecision       } from '../modules/evaluation/main'
include { FirstpassComparison   } from '../modules/evaluation/main'

workflow ASSEMBLE_AND_EVAL {

    take:
        partitioned_ch     // [test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks]
        datasets_ch        // [test_name, dataset, align_modes, partition_modes, transcriptome_modes, bambu_modes, isoquant_modes, isoseq_modes, flames_modes, stringtie2_modes]
        dataset_signal_ch  // [test_name, library_type, cage_signal_plus, cage_signal_minus, drna_signal_plus, drna_signal_minus]
        placeholders       // val: map with keys NO_ISOFORMS_BED, NO_ISOFORMS_GTF, NO_JUNCTION_TAB, NO_CAGE, NO_DRNA

    main:
        // Cache the join of partitioned outputs × dataset modes
        partitioned_with_modes = partitioned_ch.combine(datasets_ch, by: 0)

        // --- FLAIR ---
        // Filter transcriptome modes by library_type, mirroring the
        // isoquant pattern below. Modes whose names start with
        // "lib_<library_prefix>_" run only on samples with the matching
        // library_type prefix:
        //   "lib_pacbio_cDNA_*"  → only pacbio_cDNA samples
        //   "lib_ont_cDNA_*"     → only ont_cDNA samples
        //   "lib_ont_dRNA_*"     → only ont_dRNA samples
        //   "lib_pacbio_dRNA_*"  → only pacbio_dRNA samples
        //   "lib_dRNA_*"         → any *_dRNA library
        //   "lib_cDNA_*"         → any *_cDNA library
        // All other mode names run unconditionally (back-compat with
        // older configs that don't use the prefix convention).
        transcriptome_inputs = partitioned_with_modes.flatMap {
            test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks,
            dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes ->
            def partition_args = ds_partition_modes[partition_mode] ?: ''
            def junction_tab_file = dataset.junction_tab ? file(dataset.junction_tab) : placeholders.NO_JUNCTION_TAB
            def lib = dataset.library_type ?: 'unknown'
            // List of recognized exact library_type prefixes — used to
            // detect when a mode tag is targeting a SPECIFIC library
            // and shouldn't fall through to broader wildcards.
            //
            // IMPORTANT: Utils.sanitizeModeName (lib/Utils.groovy) replaces
            // underscores with dashes BEFORE these mode names land in
            // ds_transcriptome_modes. So a JSON key like
            // "lib_pacbio_cDNA_TED_global_M85_..." becomes
            // "lib-pacbio-cDNA-TED-global-M85-...". The filter must match
            // on the dash-prefixed form. We also dash-ify the library_type
            // so prefix matching works (e.g. "pacbio-cDNA-").
            def known_libs = ['pacbio-cDNA', 'pacbio-dRNA', 'ont-cDNA', 'ont-dRNA']
            def lib_dashed = lib.replace('_', '-')
            ds_transcriptome_modes.findAll { transcriptome_mode, transcriptome_args ->
                // Single keep boolean — one assignment per branch — so the
                // closure has a single explicit return point at the bottom.
                def keep
                if (!transcriptome_mode.startsWith('lib-')) {
                    keep = true                                  // unrestricted
                } else if (lib == 'unknown') {
                    keep = true                                  // can't filter
                } else {
                    def tag = transcriptome_mode.substring(4)    // drop "lib-"
                    // Exact-match library_type wins.
                    if (tag.startsWith("${lib_dashed}-")) {
                        keep = true
                    } else if (known_libs.any { it != lib_dashed && tag.startsWith("${it}-") }) {
                        // Tag targets a DIFFERENT specific library: reject.
                        keep = false
                    } else if (lib.endsWith('_dRNA') && tag.startsWith('dRNA-')) {
                        keep = true                              // dRNA-family wildcard
                    } else if (lib.endsWith('_cDNA') && tag.startsWith('cDNA-')) {
                        keep = true                              // cDNA-family wildcard
                    } else if (lib.startsWith('pacbio') && tag.startsWith('pacbio-')) {
                        keep = true                              // pacbio-platform wildcard
                    } else if (lib.startsWith('ont') && tag.startsWith('ont-')) {
                        keep = true                              // ont-platform wildcard
                    } else {
                        keep = false
                    }
                }
                return keep
            }.collect { transcriptome_mode, transcriptome_args ->
                [test_name, dataset_name, align_mode, partition_mode, partition_args, bam, bai, genome, gtf,
                 transcriptome_mode, transcriptome_args, junction_tab_file]
            }
        }
        FlairTranscriptome(transcriptome_inputs)

        // --- Bambu ---
        bambu_inputs = partitioned_with_modes.flatMap {
            test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks,
            dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes ->
            if (ds_bambu_modes.isEmpty()) return []
            ds_bambu_modes.collect { bambu_mode, bambu_args ->
                [test_name, dataset_name, align_mode, partition_mode, bam, bai, genome, gtf,
                 bambu_mode, bambu_args]
            }
        }
        BambuAssembly(bambu_inputs)

        // --- IsoQuant ---
        // Filter isoquant modes by library_type: pacbio samples → isoquant_pacbio only,
        // ont samples → isoquant_ont only, unknown → run all modes.
        // Special mode "auto": automatically injects --data_type based on library_type,
        // so all samples share the same mode name in results.
        isoquant_inputs = partitioned_with_modes.flatMap {
            test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks,
            dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes ->
            if (ds_isoquant_modes.isEmpty()) return []
            def lib = dataset.library_type ?: 'unknown'
            ds_isoquant_modes.findAll { isoquant_mode, isoquant_args ->
                if (isoquant_mode == 'auto') return true
                if (lib == 'unknown') return true
                if (lib.startsWith('pacbio') && isoquant_mode == 'pacbio') return true
                if (lib.startsWith('ont') && isoquant_mode == 'ont') return true
                return false
            }.collect { isoquant_mode, isoquant_args ->
                def resolved_args = isoquant_args
                if (isoquant_mode == 'auto') {
                    def auto_data_type = lib.startsWith('pacbio') ? 'pacbio_ccs' : 'nanopore'
                    resolved_args = "--data_type ${auto_data_type}" + (isoquant_args ? " ${isoquant_args}" : "")
                }
                [test_name, dataset_name, align_mode, partition_mode, bam, bai, genome, gtf,
                 isoquant_mode, resolved_args]
            }
        }
        IsoQuantAssembly(isoquant_inputs)

        // --- IsoSeq ---
        // IsoSeq is designed for PacBio CCS/HiFi data; only run on pacbio samples.
        isoseq_inputs = partitioned_with_modes.flatMap {
            test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks,
            dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes ->
            if (ds_isoseq_modes.isEmpty()) return []
            def lib = dataset.library_type ?: 'unknown'
            // Only run IsoSeq on PacBio data or unknown library type
            if (lib != 'unknown' && !lib.startsWith('pacbio')) return []
            ds_isoseq_modes.collect { isoseq_mode, isoseq_args ->
                [test_name, dataset_name, align_mode, partition_mode, bam, bai, genome, gtf,
                 isoseq_mode, isoseq_args]
            }
        }
        IsoSeqAssembly(isoseq_inputs)

        // --- FLAMES ---
        flames_inputs = partitioned_with_modes.flatMap {
            test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks,
            dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes ->
            if (ds_flames_modes.isEmpty()) return []
            ds_flames_modes.collect { flames_mode, flames_args ->
                [test_name, dataset_name, align_mode, partition_mode, bam, bai, genome, gtf,
                 flames_mode, flames_args]
            }
        }
        FlamesAssembly(flames_inputs)

        // --- StringTie2 ---
        stringtie2_inputs = partitioned_with_modes.flatMap {
            test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks,
            dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes ->
            if (ds_stringtie2_modes.isEmpty()) return []
            ds_stringtie2_modes.collect { stringtie2_mode, stringtie2_args ->
                [test_name, dataset_name, align_mode, partition_mode, bam, bai, genome, gtf,
                 stringtie2_mode, stringtie2_args]
            }
        }
        StringTie2Assembly(stringtie2_inputs)

        // --- Reference peaks ---
        ref_peak_inputs = partitioned_ch.map {
            test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks ->
            [test_name, dataset_name, align_mode, partition_mode, gtf]
        }
        PrepareReferencePeaks(ref_peak_inputs)

        // --- Evaluation context: partition + ref peaks + signals + metadata ---
        eval_context_ch = partitioned_ch
            .combine(PrepareReferencePeaks.out.ref_peaks, by: [0, 1, 2, 3])
            .combine(dataset_signal_ch, by: [0])
            .map { test_name, dataset_name, align_mode, partition_mode,
                   bam, bai, bed, genome, gtf, cage_peaks, drna_peaks,
                   ref_tss, ref_tts,
                   library_type,
                   cage_signal_plus, cage_signal_minus, drna_signal_plus, drna_signal_minus ->
                [test_name, dataset_name, align_mode, partition_mode,
                 bam, bai, bed, genome, gtf,
                 cage_peaks ?: placeholders.NO_CAGE, drna_peaks ?: placeholders.NO_DRNA,
                 ref_tss, ref_tts,
                 library_type,
                 cage_signal_plus, cage_signal_minus, drna_signal_plus, drna_signal_minus]
            }

        // --- Normalize assembler outputs to common shape ---
        flair_isoform_ch = FlairTranscriptome.out.transcriptome.map {
            test_name, dataset_name, align_mode, partition_mode, partition_args, transcriptome_mode,
            isoforms_bed, isoforms_gtf, isoforms_fa, isoform_counts, isoform_read_map, ted_log ->
            [test_name, dataset_name, align_mode, partition_mode, transcriptome_mode,
             isoforms_bed, placeholders.NO_ISOFORMS_GTF, isoform_read_map, ted_log]
        }

        // NO_TED_LOG placeholder for non-FLAIR assemblers
        def NO_TED_LOG = file("${workflow.workDir}/NO_TED_LOG")
        if (!NO_TED_LOG.exists()) { NO_TED_LOG.text = '' }

        bambu_isoform_ch = BambuAssembly.out.bambu_gtf.map {
            test_name, dataset_name, align_mode, partition_mode, bambu_mode, isoforms_gtf, isoform_read_map ->
            [test_name, dataset_name, align_mode, partition_mode, "bambu_${bambu_mode}",
             placeholders.NO_ISOFORMS_BED, isoforms_gtf, isoform_read_map, NO_TED_LOG]
        }

        isoquant_isoform_ch = IsoQuantAssembly.out.isoquant_gtf.map {
            test_name, dataset_name, align_mode, partition_mode, isoquant_mode, isoforms_gtf, isoform_read_map ->
            [test_name, dataset_name, align_mode, partition_mode, "isoquant_${isoquant_mode}",
             placeholders.NO_ISOFORMS_BED, isoforms_gtf, isoform_read_map, NO_TED_LOG]
        }

        isoseq_isoform_ch = IsoSeqAssembly.out.isoseq_gff.map {
            test_name, dataset_name, align_mode, partition_mode, isoseq_mode, isoforms_gff, isoform_read_map ->
            [test_name, dataset_name, align_mode, partition_mode, "isoseq_${isoseq_mode}",
             placeholders.NO_ISOFORMS_BED, isoforms_gff, isoform_read_map, NO_TED_LOG]
        }

        flames_isoform_ch = FlamesAssembly.out.flames_gtf.map {
            test_name, dataset_name, align_mode, partition_mode, flames_mode, isoforms_gtf, isoform_read_map ->
            [test_name, dataset_name, align_mode, partition_mode, "flames_${flames_mode}",
             placeholders.NO_ISOFORMS_BED, isoforms_gtf, isoform_read_map, NO_TED_LOG]
        }

        stringtie2_isoform_ch = StringTie2Assembly.out.stringtie2_gtf.map {
            test_name, dataset_name, align_mode, partition_mode, stringtie2_mode, isoforms_gtf, isoform_read_map ->
            [test_name, dataset_name, align_mode, partition_mode, "stringtie2_${stringtie2_mode}",
             placeholders.NO_ISOFORMS_BED, isoforms_gtf, isoform_read_map, NO_TED_LOG]
        }

        // Single combine: normalize → mix → join context
        // Using .mix() instead of .concat() so evaluations can start as
        // soon as any assembler finishes, rather than waiting in order.
        all_eval_inputs = flair_isoform_ch
            .mix(bambu_isoform_ch, isoquant_isoform_ch,
                 isoseq_isoform_ch, flames_isoform_ch, stringtie2_isoform_ch)
            .combine(eval_context_ch, by: [0, 1, 2, 3])

        Evaluation(all_eval_inputs)

        // --- TED End Precision (all assemblers with orthogonal peaks) ---
        // Runs for any mode where at least one peak file exists.
        // FLAIR → uses isoforms_bed; GTF-only assemblers → uses isoforms_gtf.
        // all_eval_inputs shape: [test_name, dataset_name, align_mode, partition_mode, transcriptome_mode,
        //   isoforms_bed, isoforms_gtf, isoform_read_map, ted_log, bam, bai, reads_bed, genome, gtf,
        //   cage_peaks, drna_peaks, ref_tss, ref_tts, library_type, ...]
        // We also need partition_args — retrieve from FLAIR transcriptome for FLAIR modes,
        // and from the partitioned_with_modes for non-FLAIR modes.
        // Simplest: re-derive partition_args from eval_context_ch (it doesn't carry it).
        // Instead, build a separate channel from all_eval_inputs + the partitioned_with_modes
        // to recover partition_args for each entry.

        partition_args_ch = partitioned_with_modes
            .map { test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks,
                   dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes ->
                def partition_args = ds_partition_modes[partition_mode] ?: ''
                [test_name, dataset_name, align_mode, partition_mode, partition_args]
            }

        ted_precision_inputs = all_eval_inputs
            .combine(partition_args_ch, by: [0, 1, 2, 3])
            .filter { items ->
                // items: [0-4]=keys+transcriptome_mode, [5]=isoforms_bed, [6]=isoforms_gtf,
                //   [7]=read_map, [8]=ted_log, [9]=bam, [10]=bai, [11]=reads_bed,
                //   [12]=genome, [13]=gtf, [14]=cage_peaks, [15]=drna_peaks, ..., last=partition_args
                def cage = items[14]
                def drna = items[15]
                def has_cage = cage.name != 'NO_CAGE' && cage.size() > 0
                def has_drna = drna.name != 'NO_DRNA' && drna.size() > 0
                has_cage || has_drna
            }
            .map { test_name, dataset_name, align_mode, partition_mode, transcriptome_mode,
                   isoforms_bed, isoforms_gtf, isoform_read_map, ted_log,
                   bam, bai, reads_bed, genome, gtf,
                   cage_peaks, drna_peaks, ref_tss, ref_tts, library_type,
                   cage_signal_plus, cage_signal_minus, drna_signal_plus, drna_signal_minus,
                   partition_args ->
                [test_name, dataset_name, transcriptome_mode,
                 isoforms_bed, isoforms_gtf, gtf, cage_peaks, drna_peaks, partition_args]
            }
        TedEndPrecision(ted_precision_inputs)

        // --- Firstpass vs Final comparison ---
        // Joins firstpass BED with final isoforms BED + eval context to compare
        // pre-TED and post-TED end precision/recall.
        firstpass_comparison_inputs = FlairTranscriptome.out.firstpass
            .combine(FlairTranscriptome.out.transcriptome, by: [0, 1, 2, 3, 4, 5])
            .combine(eval_context_ch, by: [0, 1, 2, 3])
            .filter { items ->
                // firstpass BED must exist (non-empty)
                def firstpass = items[4 + 2]  // firstpass_bed after 4 keys + partition_args + mode
                firstpass != null && firstpass.size() > 0
            }
            .filter { items ->
                // Need at least one orthogonal peak file
                // After combine: keys(4) + partition_args(1) + mode(1) + firstpass(1) + transcriptome_tuple(6) + eval_context(...)
                // Layout: 0-3=keys, 4=partition_args, 5=mode, 6=firstpass_bed,
                //         7=isoforms_bed, 8=gtf_out, 9=fa, 10=counts, 11=read_map, 12=ted_log,
                //         13=bam, 14=bai, 15=bed, 16=genome, 17=gtf, 18=cage, 19=drna, ...
                def cage = items[18]
                def drna = items[19]
                def has_cage = cage.name != 'NO_CAGE' && cage.size() > 0
                def has_drna = drna.name != 'NO_DRNA' && drna.size() > 0
                has_cage || has_drna
            }
            .map { test_name, dataset_name, align_mode, partition_mode, partition_args, transcriptome_mode,
                   firstpass_bed,
                   isoforms_bed, _isoforms_gtf, _isoforms_fa, _isoform_counts, _read_map, _ted_log,
                   bam, bai, bed, genome, gtf, cage_peaks, drna_peaks,
                   ref_tss, ref_tts, library_type,
                   _cage_plus, _cage_minus, _qs_plus, _qs_minus ->
                [test_name, dataset_name, transcriptome_mode,
                 firstpass_bed, isoforms_bed, gtf, cage_peaks, drna_peaks, partition_args]
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
        flair_transcriptome     = FlairTranscriptome.out.transcriptome
        flair_firstpass         = FlairTranscriptome.out.firstpass
        // CDS-aware BEDs from --predict_cds (predictProductivity output).
        // Consumed by AltEndCDSAnalysis / CDSCoverageRecovery in
        // SUMMARY_AND_VIZ. Tuple shape:
        //   [test_name, dataset_name, align_mode, partition_mode,
        //    partition_args, transcriptome_mode, isoforms_CDS.bed, isoforms_CDS.info.tsv]
        flair_cds_bed           = FlairTranscriptome.out.cds_bed
}
