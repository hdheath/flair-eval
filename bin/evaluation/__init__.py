"""
Evaluation package for FLAIR isoform analysis.

This package provides modular components for:
- TED (Transcript End Distance) metrics
- FLAIR structural evaluation
- Result synthesis

Main entry points remain in bin/ted.py, bin/flair_eval.py, and bin/synthesize_evaluations.py
"""

from .utils import (
    timed_section,
    which,
    run,
    count_lines,
    safe_f1,
    get_timing_data,
    get_logger,
    write_timing_report,
)

from .bed_utils import (
    read_bed6,
    prepare_bed6_sorted,
    extract_tss_tts_bed,
    run_closest,
    extract_distance_and_peak,
    extract_signed_distances,
    vectorized_overlap_counts,
    parse_isoform_lengths,
)

from .read_analysis import (
    parse_map_file_comprehensive,
    get_assigned_read_ids,
    count_reads_per_isoform,
    count_assigned_reads_by_type,
    count_total_alignments_bam,
    parse_isoform_ends,
    parse_read_map,
    parse_reads_bed_ends,
)

from .entropy import (
    shannon_entropy,
    compute_read_end_entropy,
)

from .truncation import (
    characterize_truncation_pattern,
)

from .peak_analysis import (
    find_recoverable_peaks,
    classify_isoform_recoverability,
    extract_read_end_positions,
    write_recoverable_peaks_bed,
    find_captured_peaks,
    build_peak_to_isoform_map,
    classify_missed_peak_reads,
    get_reads_to_isoforms,
    analyze_missed_peaks_comprehensive,
    analyze_all_recoverable_peaks_truncation,
    write_annotated_peaks_bed,
    write_troubled_regions_tsv,
    parse_read_sj_chains,
    classify_read_sj_support,
    build_isoform_peak_distance_map,
    compute_alt_end_success_metrics,
    write_well_captured_alt_regions,
    write_captured_peaks_annotated,
    classify_false_positive_endpoints,
    write_false_positive_endpoints_tsv,
)

from .motif import (
    extract_sequences_batch,
    extract_sequence_context,
    compute_position_frequency_matrix,
    compute_information_content,
    analyze_motifs_at_ends,
)

from .plots import (
    HAS_MATPLOTLIB,
    configure_plotting,
    save_figure,
    plot_distance_histogram,
    plot_distance_histogram_colored,
    plot_read_end_entropy,
    plot_read_support_distribution,
    plot_read_classification_summary,
    plot_truncation_patterns,
    plot_all_truncation_patterns,
    plot_sequence_logo,
    plot_transcript_classification,
    plot_splice_junction_support,
    plot_missed_peak_sj_support,
    plot_peak_recovery_by_expression,
    plot_peak_recovery_by_width,
    plot_peak_recovery_by_isoform_length,
    plot_read_end_frequency_at_peaks,
    plot_read_end_frequency_stratified_by_width,
    plot_peak_width_histogram,
    plot_internal_priming_summary,
    plot_junction_chain_end_variation_histogram,
    plot_gene_variation_proportions,
    plot_end_distance_dashboard,
    plot_peak_recovery_dashboard,
    plot_missed_peak_diagnostics_dashboard,
    plot_proximal_apa_distance_support,
    plot_single_exon_peaks_support,
)

from .end_variation import (
    summarize_junction_chain_end_variation,
    classify_genes_by_variation,
)

from .dexseq_ends import (
    build_gene_end_bins,
    cluster_end_positions,
    estimate_dispersion_per_gene,
    fit_mean_dispersion_trend,
    shrink_dispersions,
    this_vs_others_test,
    test_all_bins_per_gene,
    benjamini_hochberg,
    compute_end_confidence_scores,
    score_novel_ends,
    run_dexseq_end_analysis,
)

from .concordance_metrics import (
    compute_concordance_metrics,
)

from .signal_profile import (
    compute_meta_profile,
    compute_meta_profile_metrics,
    compute_boundary_scores,
    summarize_boundary_scores,
    write_boundary_scores_tsv,
    compute_signal_profile_metrics,
    plot_meta_profiles,
    plot_boundary_score_distributions,
    plot_boundary_score_cdf,
)

from .ted_core import (
    tss_tts_metrics,
    calculate_ted_metrics,
)

from .training_data import (
    cluster_read_ends,
    compute_cluster_features,
    compute_sequence_features,
    compute_read_to_isoform_features,
    parse_gtf_ends,
    label_clusters_with_peaks,
    generate_training_data,
    write_training_data,
    generate_and_write_training_data,
    check_internal_priming,
    compute_internal_priming_features,
)

from .isoform_end_training import (
    compute_isoform_peak_assignment,
    compute_isoform_end_features,
    compute_isoform_sequence_features,
    generate_isoform_training_data,
    write_isoform_training_data,
    generate_and_write_isoform_training_data,
)

from .flair_structural import (
    get_chromtoint,
    get_regions,
    get_intersect_count,
    extract_sj_info,
    parse_gtf_transcripts,
    build_reference_structures,
    classify_transcripts,
    SINGLE_EXON_END_WINDOW,
)

from .synthesize import (
    parse_filename_metadata,
    parse_tsv_file,
    process_evaluation_files,
    write_tsv_merged,
    METADATA_FIELDS,
)

from .tool_divergence import (
    compute_pairwise_jaccard,
    compute_motif_collapse,
    MOTIF_LIBRARY,
)

# precision_recall_plot requires pandas/matplotlib - import conditionally
try:
    from .precision_recall_plot import (
        identify_assembler,
        load_evaluation_files,
        create_precision_recall_plot,
        ASSEMBLER_COLORS,
        ASSEMBLER_NAMES,
    )
    HAS_PRECISION_RECALL_PLOT = True
except ImportError:
    HAS_PRECISION_RECALL_PLOT = False

__all__ = [
    # utils
    'timed_section',
    'which',
    'run',
    'count_lines',
    'safe_f1',
    'get_timing_data',
    'get_logger',
    'write_timing_report',
    # bed_utils
    'read_bed6',
    'prepare_bed6_sorted',
    'extract_tss_tts_bed',
    'run_closest',
    'extract_distance_and_peak',
    'extract_signed_distances',
    'vectorized_overlap_counts',
    # read_analysis
    'parse_map_file_comprehensive',
    'get_assigned_read_ids',
    'count_reads_per_isoform',
    'count_assigned_reads_by_type',
    'count_total_alignments_bam',
    'parse_isoform_ends',
    'parse_read_map',
    'parse_reads_bed_ends',
    # entropy
    'shannon_entropy',
    'compute_read_end_entropy',
    # truncation
    'characterize_truncation_pattern',
    # peak_analysis
    'find_recoverable_peaks',
    'classify_isoform_recoverability',
    'extract_read_end_positions',
    'write_recoverable_peaks_bed',
    'find_captured_peaks',
    'classify_missed_peak_reads',
    'get_reads_to_isoforms',
    'analyze_missed_peaks_comprehensive',
    'analyze_all_recoverable_peaks_truncation',
    'write_annotated_peaks_bed',
    'write_troubled_regions_tsv',
    'parse_read_sj_chains',
    'classify_read_sj_support',
    'build_isoform_peak_distance_map',
    'compute_alt_end_success_metrics',
    'write_well_captured_alt_regions',
    'write_captured_peaks_annotated',
    # motif
    'extract_sequences_batch',
    'extract_sequence_context',
    'compute_position_frequency_matrix',
    'compute_information_content',
    'analyze_motifs_at_ends',
    # plots
    'HAS_MATPLOTLIB',
    'configure_plotting',
    'save_figure',
    'plot_distance_histogram',
    'plot_distance_histogram_colored',
    'plot_read_end_entropy',
    'plot_read_support_distribution',
    'plot_read_classification_summary',
    'plot_truncation_patterns',
    'plot_all_truncation_patterns',
    'plot_sequence_logo',
    'plot_transcript_classification',
    'plot_splice_junction_support',
    'plot_missed_peak_sj_support',
    'plot_peak_recovery_by_expression',
    'plot_peak_recovery_by_width',
    'plot_read_end_frequency_at_peaks',
    'plot_read_end_frequency_stratified_by_width',
    'plot_peak_width_histogram',
    'plot_internal_priming_summary',
    'plot_junction_chain_end_variation_histogram',
    'plot_gene_variation_proportions',
    'plot_end_distance_dashboard',
    'plot_peak_recovery_dashboard',
    'plot_missed_peak_diagnostics_dashboard',
    # end variation
    'summarize_junction_chain_end_variation',
    'classify_genes_by_variation',
    # concordance_metrics
    'compute_concordance_metrics',
    # signal_profile
    'compute_meta_profile',
    'compute_meta_profile_metrics',
    'compute_boundary_scores',
    'summarize_boundary_scores',
    'write_boundary_scores_tsv',
    'compute_signal_profile_metrics',
    'plot_meta_profiles',
    'plot_boundary_score_distributions',
    'plot_boundary_score_cdf',
    # dexseq_ends
    'build_gene_end_bins',
    'cluster_end_positions',
    'estimate_dispersion_per_gene',
    'fit_mean_dispersion_trend',
    'shrink_dispersions',
    'this_vs_others_test',
    'test_all_bins_per_gene',
    'benjamini_hochberg',
    'compute_end_confidence_scores',
    'score_novel_ends',
    'run_dexseq_end_analysis',
    # ted_core
    'tss_tts_metrics',
    'calculate_ted_metrics',
    # training_data
    'cluster_read_ends',
    'compute_cluster_features',
    'compute_sequence_features',
    'compute_read_to_isoform_features',
    'parse_gtf_ends',
    'label_clusters_with_peaks',
    'generate_training_data',
    'write_training_data',
    'generate_and_write_training_data',
    'check_internal_priming',
    'compute_internal_priming_features',
    # flair_structural
    'get_chromtoint',
    'get_regions',
    'get_intersect_count',
    'extract_sj_info',
    'parse_gtf_transcripts',
    'build_reference_structures',
    'classify_transcripts',
    'SINGLE_EXON_END_WINDOW',
    # synthesize
    'parse_filename_metadata',
    'parse_tsv_file',
    'process_evaluation_files',
    'write_tsv_merged',
    'METADATA_FIELDS',
    # tool_divergence
    'compute_pairwise_jaccard',
    'compute_motif_collapse',
    'MOTIF_LIBRARY',
    # precision_recall_plot availability flag
    'HAS_PRECISION_RECALL_PLOT',
]

# Conditionally add precision_recall_plot exports if available
if HAS_PRECISION_RECALL_PLOT:
    __all__.extend([
        'identify_assembler',
        'load_evaluation_files',
        'create_precision_recall_plot',
        'ASSEMBLER_COLORS',
        'ASSEMBLER_NAMES',
    ])
