// Module: BambuAssembly
// Bambu: R-based transcript discovery and quantification.
//
// IMPORTANT: XGBoost version compatibility issue
// Bambu's pre-trained XGBoost model is incompatible with XGBoost >= 2.1.0.
// Options: downgrade to r-xgboost=1.7.6, use discovery=FALSE, or update Bambu to 3.5.1+.

process BambuAssembly {
    publishDir "${params.outdir}/assemblers/bambu/${test_name}", mode: 'symlink'
    publishDir "${params.outdir}/logs/${test_name}", mode: 'copy', pattern: '.command.{log,err}', saveAs: { "${dataset_name}_${align_mode}_${partition_mode}_bambu_${bambu_mode}_${it}" }
    errorStrategy 'terminate'
    tag "${dataset_name}_${align_mode}_${partition_mode}_bambu_${bambu_mode}"

    input:
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          path(bam), path(bai), path(genome), path(gtf),
          val(bambu_mode), val(bambu_args)

    output:
    // counts_transcript.txt is included in the main bambu_gtf tuple so the
    // evaluation subworkflow can forward it to ted.py's --counts arg.
    // Bambu, like IsoQuant, emits reference transcripts in the output GTF
    // even when no reads support them in the sample, so without this filter
    // isoforms_observed is inflated by zero-support reference carry-throughs.
    tuple val(test_name), val(dataset_name), val(align_mode), val(partition_mode),
          val(bambu_mode), path("${dataset_name}_${align_mode}_${partition_mode}_bambu_${bambu_mode}.gtf"),
          path("${dataset_name}_${align_mode}_${partition_mode}_bambu_${bambu_mode}_read_map.txt"),
          path("${dataset_name}_${align_mode}_${partition_mode}_bambu_${bambu_mode}_counts_transcript.txt", optional: true),
          emit: bambu_gtf
    path "${dataset_name}_${align_mode}_${partition_mode}_bambu_${bambu_mode}_counts_transcript.txt", optional: true, emit: transcript_counts
    path "${dataset_name}_${align_mode}_${partition_mode}_bambu_${bambu_mode}_counts_gene.txt", optional: true, emit: gene_counts

    script:
    def output_prefix = "${dataset_name}_${align_mode}_${partition_mode}_bambu_${bambu_mode}"
    def bambu_params = bambu_args ?: ""
    """
    #!/private/home/hdheath/miniforge3/envs/bambu/bin/Rscript
    # Absolute path to Rscript (2026-05-14): the conda activate emitted by
    # Nextflow's process directive works on the login node but on some SLURM
    # compute nodes does NOT prepend the env bin to PATH, causing
    # `#!/usr/bin/env Rscript` to find the system R (which has no `bambu`
    # package). Pinning to the env's Rscript bypasses the activation entirely.
    # v2: wrap read-map extraction in tryCatch — on full-genome runs the
    # rbind of per-read data.frames can fail at R level; don't let that
    # error kill the job since the GTF (the main deliverable) is already written.

    library(bambu)

    annotations <- prepareAnnotations("${gtf}")

    se <- bambu(
        reads = "${bam}",
        annotations = annotations,
        genome = "${genome}",
        ncore = ${task.cpus},
        trackReads = TRUE${bambu_params ? ", ${bambu_params}" : ""}
    )

    writeBambuOutput(se, path = ".", prefix = "${output_prefix}")
    file.rename("${output_prefix}extended_annotations.gtf", "${output_prefix}.gtf")
    # Bambu writes counts files as `<prefix>counts_transcript.txt` (no
    # separator); the pipeline expects `<prefix>_counts_transcript.txt`.
    # Also: on some bambu versions / single-sample inputs, writeBambuOutput
    # silently skips the counts file. To make the pipeline's --counts filter
    # work reliably for Bambu, we write the counts file ourselves directly
    # from rowData(se)\$readCount (raw integer read count per transcript).
    # This gives ted.py what it needs: tab-separated `transcript_id\\tcount`
    # which it parses to filter isoforms_observed to count >= 1.
    counts_outfile <- "${output_prefix}_counts_transcript.txt"
    if (file.exists("${output_prefix}counts_transcript.txt")) {
        file.rename("${output_prefix}counts_transcript.txt", counts_outfile)
    } else if ("readCount" %in% colnames(rowData(se))) {
        cat("writeBambuOutput skipped counts file; writing from rowData(se)\$readCount\\n")
        rd <- as.data.frame(rowData(se))
        counts_df <- data.frame(
            feature_id = rd\$TXNAME,
            count      = rd\$readCount,
            stringsAsFactors = FALSE
        )
        # ted.py expects `#feature_id\\tcount` (IsoQuant convention) or
        # whitespace-separated. Use the IsoQuant header format for symmetry.
        cat("#feature_id\\tcount\\n", file = counts_outfile)
        write.table(counts_df, file = counts_outfile, sep = "\\t",
                    quote = FALSE, row.names = FALSE, col.names = FALSE,
                    append = TRUE)
        cat(sprintf("Wrote %d transcript counts to %s\\n", nrow(counts_df), counts_outfile))
    } else {
        cat("WARNING: no counts data available; touching empty placeholder\\n")
        file.create(counts_outfile)
    }
    if (file.exists("${output_prefix}counts_gene.txt")) {
        file.rename("${output_prefix}counts_gene.txt", "${output_prefix}_counts_gene.txt")
    } else {
        file.create("${output_prefix}_counts_gene.txt")
    }

    read_map_ok <- tryCatch({
        read_maps <- metadata(se)\$readToTranscriptMaps
        if (!is.null(read_maps) && length(read_maps) > 0) {
            map_df <- read_maps[[1]]
            tx_names <- rowData(se)\$TXNAME

            # Vectorized extraction — collect lists of (read_id, transcript_id) pairs
            result_list <- vector("list", nrow(map_df))
            for (i in seq_len(nrow(map_df))) {
                rid <- map_df\$readId[i]
                eq_idx <- map_df\$equalMatches[[i]]
                if (length(eq_idx) > 0 && !all(is.na(eq_idx))) {
                    valid <- eq_idx[eq_idx >= 1 & eq_idx <= length(tx_names)]
                    if (length(valid) > 0) {
                        result_list[[i]] <- data.frame(read_id = rid, transcript_id = tx_names[valid], stringsAsFactors = FALSE)
                    }
                } else {
                    compat_idx <- map_df\$compatibleMatches[[i]]
                    if (length(compat_idx) > 0 && !all(is.na(compat_idx))) {
                        valid <- compat_idx[compat_idx >= 1 & compat_idx <= length(tx_names)]
                        if (length(valid) > 0) {
                            result_list[[i]] <- data.frame(read_id = rid, transcript_id = tx_names[valid], stringsAsFactors = FALSE)
                        }
                    }
                }
            }
            assignments <- do.call(rbind, result_list[!sapply(result_list, is.null)])
            if (is.null(assignments)) assignments <- data.frame(read_id = character(), transcript_id = character(), stringsAsFactors = FALSE)

            write.table(assignments, file = "${output_prefix}_bambu_reads.tsv", sep = "\\t",
                         row.names = FALSE, col.names = TRUE, quote = FALSE)
            cat(sprintf("Exported %d read-transcript assignments from Bambu\\n", nrow(assignments)))
        } else {
            writeLines(character(0), "${output_prefix}_bambu_reads.tsv")
            cat("WARNING: No readToTranscriptMaps in Bambu output\\n")
        }
        TRUE
    }, error = function(e) {
        cat(sprintf("WARNING: Bambu read-map extraction failed: %s\\n", conditionMessage(e)))
        cat("Writing empty read-map placeholder; GTF output is unaffected.\\n")
        writeLines(character(0), "${output_prefix}_bambu_reads.tsv")
        FALSE
    })

    system2("python", c("${projectDir}/bin/convert_read_map.py",
                        "--bambu-read-map", "${output_prefix}_bambu_reads.tsv",
                        "--output", "${output_prefix}_read_map.txt",
                        "--verbose"))
    """
}
