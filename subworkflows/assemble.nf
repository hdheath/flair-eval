// =============================================================================
// Subworkflow: ASSEMBLE  (STAGE 1)
// =============================================================================
// Runs ALL transcriptome assemblers (FLAIR/TED, Bambu, IsoQuant, IsoSeq,
// StringTie2, LRAA; FLAMES is intentionally dropped) on partitioned data and
// emits:
//   * per-tool normalized isoform tuples (for the in-process monolith path)
//   * FLAIR firstpass + CDS BEDs
//   * a MANIFEST channel: one row per produced assembly, carrying ABSOLUTE
//     published paths + the matched orthogonal/reference paths, so STAGE 2
//     (EVALUATE) can resolve every input without re-running assembly.
//
// This is the SLOW, SLURM-heavy stage. Settle it once, then iterate EVALUATE.
//
// Inputs:
//   partitioned_ch    — FlairPartition.out.partitioned
//                       [test_name, dataset_name, align_mode, partition_mode,
//                        bam, bai, bed, genome(val), gtf, cage_peaks, drna_peaks]
//   datasets_ch       — master datasets channel (mode maps + junction_tab + signals)
//   dataset_signal_ch — per-dataset signal bedGraph paths + library_type
//   placeholders      — map of placeholder file paths
//
// Emits:
//   isoform_tuples    — mixed per-tool normalized isoform tuples (monolith path)
//   eval_context_seed — [test_name, dataset_name, align_mode, partition_mode,
//                        bam, bai, reads_bed, genome, gtf, cage_peaks, drna_peaks,
//                        partition_args] (monolith convenience; ref peaks added in EVALUATE)
//   flair_firstpass   — FlairTranscriptome.out.firstpass
//   flair_cds_bed     — Channel.empty() (separate CDS BED obsolete; inline now)
//   manifest_rows     — per-assembly CSV row strings (header-less) for the manifest
// =============================================================================

include { FlairTranscriptome   } from '../modules/transcriptome/flair/main'
include { BambuAssembly        } from '../modules/transcriptome/bambu/main'
include { IsoQuantAssembly     } from '../modules/transcriptome/isoquant/main'
include { IsoSeqAssembly       } from '../modules/transcriptome/isoseq/main'
include { StringTie2Assembly   } from '../modules/transcriptome/stringtie2/main'
include { LraaAssembly         } from '../modules/transcriptome/lraa/main'
// NOTE: FLAMES is intentionally NOT included — dropped per the 2026-05 decision.

workflow ASSEMBLE {

    take:
        partitioned_ch     // see header
        datasets_ch        // [test_name, dataset, align_modes, partition_modes, transcriptome_modes, bambu, isoquant, isoseq, flames, stringtie2, lraa]
        dataset_signal_ch  // [test_name, library_type, cage_signal_plus, cage_signal_minus, drna_signal_plus, drna_signal_minus]
        placeholders       // val map: NO_ISOFORMS_BED, NO_ISOFORMS_GTF, NO_JUNCTION_TAB, NO_CAGE, NO_DRNA

    main:
        // Join partitioned outputs × dataset modes (single combine, reused below)
        partitioned_with_modes = partitioned_ch.combine(datasets_ch, by: 0)

        // ------------------------------------------------------------------
        // FLAIR / TED
        // ------------------------------------------------------------------
        // (Library-type mode filtering preserved verbatim from ASSEMBLE_AND_EVAL.)
        transcriptome_inputs = partitioned_with_modes.flatMap {
            test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks,
            dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes, ds_lraa_modes ->
            def partition_args = ds_partition_modes[partition_mode] ?: ''
            def junction_tab_file = dataset.junction_tab ? file(dataset.junction_tab) : placeholders.NO_JUNCTION_TAB
            def lib = dataset.library_type ?: 'unknown'
            def known_libs = ['pacbio-cDNA', 'pacbio-dRNA', 'ont-cDNA', 'ont-dRNA']
            def lib_dashed = lib.replace('_', '-')
            ds_transcriptome_modes.findAll { transcriptome_mode, transcriptome_args ->
                def keep
                if (!transcriptome_mode.startsWith('lib-')) {
                    keep = true
                } else if (lib == 'unknown') {
                    keep = true
                } else {
                    def tag = transcriptome_mode.substring(4)
                    if (tag.startsWith("${lib_dashed}-")) {
                        keep = true
                    } else if (known_libs.any { it != lib_dashed && tag.startsWith("${it}-") }) {
                        keep = false
                    } else if (lib.endsWith('_dRNA') && tag.startsWith('dRNA-')) {
                        keep = true
                    } else if (lib.endsWith('_cDNA') && tag.startsWith('cDNA-')) {
                        keep = true
                    } else if (lib.startsWith('pacbio') && tag.startsWith('pacbio-')) {
                        keep = true
                    } else if (lib.startsWith('ont') && tag.startsWith('ont-')) {
                        keep = true
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

        // ------------------------------------------------------------------
        // Bambu
        // ------------------------------------------------------------------
        bambu_inputs = partitioned_with_modes.flatMap {
            test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks,
            dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes, ds_lraa_modes ->
            if (ds_bambu_modes.isEmpty()) return []
            ds_bambu_modes.collect { bambu_mode, bambu_args ->
                [test_name, dataset_name, align_mode, partition_mode, bam, bai, genome, gtf,
                 bambu_mode, bambu_args]
            }
        }
        BambuAssembly(bambu_inputs)

        // ------------------------------------------------------------------
        // IsoQuant
        // ------------------------------------------------------------------
        isoquant_inputs = partitioned_with_modes.flatMap {
            test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks,
            dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes, ds_lraa_modes ->
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

        // ------------------------------------------------------------------
        // IsoSeq (PacBio only)
        // ------------------------------------------------------------------
        isoseq_inputs = partitioned_with_modes.flatMap {
            test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks,
            dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes, ds_lraa_modes ->
            if (ds_isoseq_modes.isEmpty()) return []
            def lib = dataset.library_type ?: 'unknown'
            if (lib != 'unknown' && !lib.startsWith('pacbio')) return []
            ds_isoseq_modes.collect { isoseq_mode, isoseq_args ->
                [test_name, dataset_name, align_mode, partition_mode, bam, bai, genome, gtf,
                 isoseq_mode, isoseq_args]
            }
        }
        IsoSeqAssembly(isoseq_inputs)

        // ------------------------------------------------------------------
        // StringTie2
        // ------------------------------------------------------------------
        stringtie2_inputs = partitioned_with_modes.flatMap {
            test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks,
            dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes, ds_lraa_modes ->
            if (ds_stringtie2_modes.isEmpty()) return []
            ds_stringtie2_modes.collect { stringtie2_mode, stringtie2_args ->
                [test_name, dataset_name, align_mode, partition_mode, bam, bai, genome, gtf,
                 stringtie2_mode, stringtie2_args]
            }
        }
        StringTie2Assembly(stringtie2_inputs)

        // ------------------------------------------------------------------
        // LRAA (--HiFi auto-injected for pacbio)
        // ------------------------------------------------------------------
        lraa_inputs = partitioned_with_modes.flatMap {
            test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks,
            dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes, ds_lraa_modes ->
            if (ds_lraa_modes.isEmpty()) return []
            def lib = dataset.library_type ?: 'unknown'
            ds_lraa_modes.collect { lraa_mode, lraa_args ->
                def resolved = (lraa_args ?: '').toString()
                if (lib.startsWith('pacbio') && !(resolved =~ /(?<![\w-])--HiFi(?![\w-])/)) {
                    resolved = "--HiFi ${resolved}".trim()
                }
                [test_name, dataset_name, align_mode, partition_mode, bam, bai, genome, gtf,
                 lraa_mode, resolved]
            }
        }
        LraaAssembly(lraa_inputs)

        // ==================================================================
        // Normalize assembler outputs to the common isoform tuple shape:
        //   [test_name, dataset_name, align_mode, partition_mode,
        //    transcriptome_mode, isoforms_bed, isoforms_gtf, isoform_read_map,
        //    ted_log, isoform_counts]
        // ==================================================================
        def NO_TED_LOG = file("${workflow.workDir}/NO_TED_LOG")
        if (!NO_TED_LOG.exists()) { NO_TED_LOG.text = '' }
        def NO_COUNTS = file("${workflow.workDir}/NO_COUNTS")
        if (!NO_COUNTS.exists()) { NO_COUNTS.text = '' }

        flair_isoform_ch = FlairTranscriptome.out.transcriptome.map {
            test_name, dataset_name, align_mode, partition_mode, partition_args, transcriptome_mode,
            isoforms_bed, isoforms_gtf, isoforms_fa, isoform_counts, isoform_read_map, ted_log ->
            [test_name, dataset_name, align_mode, partition_mode, transcriptome_mode,
             isoforms_bed, placeholders.NO_ISOFORMS_GTF, isoform_read_map, ted_log,
             isoform_counts ?: NO_COUNTS]
        }

        bambu_isoform_ch = BambuAssembly.out.bambu_gtf.map {
            test_name, dataset_name, align_mode, partition_mode, bambu_mode, isoforms_gtf, isoform_read_map, isoform_counts ->
            [test_name, dataset_name, align_mode, partition_mode, "bambu_${bambu_mode}",
             placeholders.NO_ISOFORMS_BED, isoforms_gtf, isoform_read_map, NO_TED_LOG,
             isoform_counts ?: NO_COUNTS]
        }

        isoquant_isoform_ch = IsoQuantAssembly.out.isoquant_gtf.map {
            test_name, dataset_name, align_mode, partition_mode, isoquant_mode, isoforms_gtf, isoform_read_map, isoform_counts ->
            [test_name, dataset_name, align_mode, partition_mode, "isoquant_${isoquant_mode}",
             placeholders.NO_ISOFORMS_BED, isoforms_gtf, isoform_read_map, NO_TED_LOG,
             isoform_counts ?: NO_COUNTS]
        }

        isoseq_isoform_ch = IsoSeqAssembly.out.isoseq_gff.map {
            test_name, dataset_name, align_mode, partition_mode, isoseq_mode, isoforms_gff, isoform_read_map, isoform_counts ->
            [test_name, dataset_name, align_mode, partition_mode, "isoseq_${isoseq_mode}",
             placeholders.NO_ISOFORMS_BED, isoforms_gff, isoform_read_map, NO_TED_LOG,
             isoform_counts ?: NO_COUNTS]
        }

        stringtie2_isoform_ch = StringTie2Assembly.out.stringtie2_gtf.map {
            test_name, dataset_name, align_mode, partition_mode, stringtie2_mode, isoforms_gtf, isoform_read_map ->
            [test_name, dataset_name, align_mode, partition_mode, "stringtie2_${stringtie2_mode}",
             placeholders.NO_ISOFORMS_BED, isoforms_gtf, isoform_read_map, NO_TED_LOG,
             NO_COUNTS]
        }

        lraa_isoform_ch = LraaAssembly.out.lraa_gtf.map {
            test_name, dataset_name, align_mode, partition_mode, lraa_mode, isoforms_gtf, isoform_read_map, isoform_counts ->
            [test_name, dataset_name, align_mode, partition_mode, "lraa_${lraa_mode}",
             placeholders.NO_ISOFORMS_BED, isoforms_gtf, isoform_read_map, NO_TED_LOG,
             isoform_counts ?: NO_COUNTS]
        }

        isoform_tuples = flair_isoform_ch
            .mix(bambu_isoform_ch, isoquant_isoform_ch,
                 isoseq_isoform_ch, stringtie2_isoform_ch, lraa_isoform_ch)

        // ==================================================================
        // Eval-context seed: partition outputs + partition_args, keyed for
        // the monolith path. PrepareReferencePeaks (ref_tss/ref_tts) is added
        // inside EVALUATE so the cut is clean.
        // ==================================================================
        partition_args_ch = partitioned_with_modes
            .map { test_name, dataset_name, align_mode, partition_mode, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks,
                   dataset, ds_align_modes, ds_partition_modes, ds_transcriptome_modes, ds_bambu_modes, ds_isoquant_modes, ds_isoseq_modes, ds_flames_modes, ds_stringtie2_modes, ds_lraa_modes ->
                def partition_args = ds_partition_modes[partition_mode] ?: ''
                [test_name, dataset_name, align_mode, partition_mode, partition_args]
            }

        eval_context_seed = partitioned_ch
            .combine(partition_args_ch, by: [0, 1, 2, 3])
            .map { test_name, dataset_name, align_mode, partition_mode,
                   bam, bai, bed, genome, gtf, cage_peaks, drna_peaks, partition_args ->
                [test_name, dataset_name, align_mode, partition_mode,
                 bam, bai, bed, genome, gtf,
                 cage_peaks ?: placeholders.NO_CAGE, drna_peaks ?: placeholders.NO_DRNA,
                 partition_args]
            }

        // ==================================================================
        // MANIFEST: one CSV row per produced assembly.
        // We start from the published per-tool isoform tuple (so a row exists
        // ONLY for assemblies that actually produced isoforms), then join:
        //   * eval_context_seed   → partition bam/bai/reads_bed/genome/gtf/peaks + partition_args
        //   * dataset_signal_ch   → library_type + signal tracks
        //   * partition_region    → derived from partition_args (--region ...)
        //   * firstpass / cds      → FLAIR only (placeholder otherwise)
        // All file paths are emitted as ABSOLUTE (.toRealPath where possible,
        // else .toAbsolutePath) so STAGE 2 resolves them from any launch dir.
        // ==================================================================
        def NO_FIRSTPASS = file("${workflow.workDir}/NO_FIRSTPASS")
        if (!NO_FIRSTPASS.exists()) { NO_FIRSTPASS.text = '' }
        def NO_CDS = file("${workflow.workDir}/NO_CDS")
        if (!NO_CDS.exists()) { NO_CDS.text = '' }

        // FLAIR firstpass keyed by the 5-tuple (firstpass is FLAIR-only)
        flair_firstpass_keyed = FlairTranscriptome.out.firstpass.map {
            test_name, dataset_name, align_mode, partition_mode, partition_args, transcriptome_mode, firstpass_bed ->
            [[test_name, dataset_name, align_mode, partition_mode, transcriptome_mode],
             firstpass_bed ?: NO_FIRSTPASS]
        }
        // The separate CDS BED is OBSOLETE under the new dev FLAIR (productivity
        // is inline in .isoforms.bed). FlairTranscriptome no longer emits cds_bed,
        // so this keyed channel is empty: the remainder:true join below fills the
        // manifest cds_bed/cds_info columns with the NO_CDS placeholder for every
        // row. Column stays in the 28-col schema for back-compat.
        flair_cds_keyed = Channel.empty()

        // signal map keyed by test_name
        signal_keyed = dataset_signal_ch.map {
            test_name, library_type, cage_signal_plus, cage_signal_minus, drna_signal_plus, drna_signal_minus ->
            [test_name, library_type, cage_signal_plus, cage_signal_minus, drna_signal_plus, drna_signal_minus]
        }

        manifest_rows = isoform_tuples
            .map { test_name, dataset_name, align_mode, partition_mode, transcriptome_mode,
                   isoforms_bed, isoforms_gtf, isoform_read_map, ted_log, isoform_counts ->
                [[test_name, dataset_name, align_mode, partition_mode],
                 [test_name, dataset_name, align_mode, partition_mode, transcriptome_mode,
                  isoforms_bed, isoforms_gtf, isoform_read_map, ted_log, isoform_counts]]
            }
            // attach partition context (eval_context_seed) by the 4-key
            .combine(
                eval_context_seed.map { test_name, dataset_name, align_mode, partition_mode,
                                        bam, bai, bed, genome, gtf, cage_peaks, drna_peaks, partition_args ->
                    [[test_name, dataset_name, align_mode, partition_mode],
                     bam, bai, bed, genome, gtf, cage_peaks, drna_peaks, partition_args]
                },
                by: 0
            )
            .map { key4, iso, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks, partition_args ->
                // re-key by the full 5-tuple to attach firstpass/cds
                def key5 = [iso[0], iso[1], iso[2], iso[3], iso[4]]
                [key5, iso, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks, partition_args]
            }
            .map { key5, iso, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks, partition_args ->
                [key5, [iso, bam, bai, bed, genome, gtf, cage_peaks, drna_peaks, partition_args]]
            }
            .join(flair_firstpass_keyed, remainder: true)
            .map { key5, payload, firstpass ->
                // remainder:true → non-FLAIR rows get firstpass == null
                payload != null ? [key5, payload, firstpass ?: NO_FIRSTPASS] : null
            }
            .filter { it != null }
            .join(flair_cds_keyed.map { k, cds_bed, cds_info -> [k, cds_bed, cds_info] }, remainder: true)
            .map { key5, payload, firstpass, cds_bed, cds_info ->
                payload != null ? [payload, firstpass, cds_bed ?: NO_CDS, cds_info ?: NO_CDS] : null
            }
            .filter { it != null }
            // attach signals by test_name
            .map { payload, firstpass, cds_bed, cds_info ->
                def test_name = payload[0][0]
                [test_name, payload, firstpass, cds_bed, cds_info]
            }
            .combine(signal_keyed, by: 0)
            .map { test_name, payload, firstpass, cds_bed, cds_info,
                   library_type, cage_signal_plus, cage_signal_minus, drna_signal_plus, drna_signal_minus ->
                def iso          = payload[0]
                def bam          = payload[1]
                def bai          = payload[2]
                def reads_bed    = payload[3]
                def genome       = payload[4]   // val string (absolute canonical path)
                def part_gtf     = payload[5]
                def cage_peaks   = payload[6]
                def drna_peaks   = payload[7]
                def partition_args = payload[8]

                // Derive partition_region from partition_args (--region a:1-2 c:3-4 ...)
                def region = ''
                def m = (partition_args =~ /--region\s+(.+?)(?:\s+--|$)/)
                if (m.find()) { region = m.group(1).trim().replaceAll(/\s+/, ';') }

                // Resolve to absolute paths. Published outputs are symlinks in
                // the workdir; toRealPath() follows to the published/real file
                // so STAGE 2 can read them regardless of launch dir.
                def abspath = { f ->
                    if (f == null) return ''
                    def p = file(f.toString())
                    try { return p.toRealPath().toString() }
                    catch (Exception e) { return p.toAbsolutePath().toString() }
                }

                // tool = transcriptome_mode prefix (flair gets bare mode; others carry tool_ prefix already)
                def tmode = iso[4]
                def tool
                if (tmode.startsWith('bambu_'))        tool = 'bambu'
                else if (tmode.startsWith('isoquant_')) tool = 'isoquant'
                else if (tmode.startsWith('isoseq_'))   tool = 'isoseq'
                else if (tmode.startsWith('stringtie2_')) tool = 'stringtie2'
                else if (tmode.startsWith('lraa_'))     tool = 'lraa'
                else                                     tool = 'flair'

                // CSV row — column order MUST match Utils.MANIFEST_COLUMNS
                [iso[0],              // test_name
                 iso[1],              // dataset_name
                 iso[2],              // align_mode
                 iso[3],              // partition_mode
                 tmode,               // transcriptome_mode
                 tool,                // tool
                 abspath(iso[5]),     // isoforms_bed
                 abspath(iso[6]),     // isoforms_gtf
                 abspath(iso[7]),     // isoform_read_map
                 abspath(iso[8]),     // ted_log
                 abspath(iso[9]),     // isoform_counts
                 abspath(firstpass),  // firstpass
                 abspath(cds_bed),    // cds_bed
                 abspath(cds_info),   // cds_info
                 abspath(bam),        // bam
                 abspath(bai),        // bai
                 abspath(reads_bed),  // reads_bed (partition .bed)
                 (genome ?: '').toString(),  // genome (canonical absolute path, val)
                 abspath(part_gtf),   // gtf (partition annotation.gtf)
                 abspath(cage_peaks), // cage_peaks (partition cage.bed)
                 abspath(drna_peaks), // drna_peaks (partition drna.bed)
                 (cage_signal_plus  ?: '').toString(),
                 (cage_signal_minus ?: '').toString(),
                 (drna_signal_plus  ?: '').toString(),
                 (drna_signal_minus ?: '').toString(),
                 (library_type ?: 'unknown').toString(),
                 region,              // partition_region
                 (partition_args ?: '').toString()
                ].collect { v -> def s = v.toString(); s.contains(',') ? "\"${s}\"" : s }.join(',')
            }

    emit:
        isoform_tuples    = isoform_tuples
        eval_context_seed = eval_context_seed
        flair_firstpass   = FlairTranscriptome.out.firstpass
        // CDS BED is obsolete (inline productivity now) — emit empty so
        // downstream CDS analyses are gracefully fed nothing (no-op).
        flair_cds_bed     = Channel.empty()
        manifest_rows     = manifest_rows
}
