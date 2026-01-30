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
    classify_missed_peak_reads,
    get_reads_to_isoforms,
    analyze_missed_peaks_comprehensive,
    analyze_all_recoverable_peaks_truncation,
    write_annotated_peaks_bed,
    write_troubled_regions_tsv,
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
    plot_distance_histogram,
    plot_distance_histogram_colored,
    plot_read_end_entropy,
    plot_read_support_distribution,
    plot_read_classification_summary,
    plot_truncation_patterns,
    plot_all_truncation_patterns,
    plot_sequence_logo,
)

from .ted_core import (
    tss_tts_metrics,
    calculate_ted_metrics,
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
    # motif
    'extract_sequences_batch',
    'extract_sequence_context',
    'compute_position_frequency_matrix',
    'compute_information_content',
    'analyze_motifs_at_ends',
    # plots
    'HAS_MATPLOTLIB',
    'plot_distance_histogram',
    'plot_distance_histogram_colored',
    'plot_read_end_entropy',
    'plot_read_support_distribution',
    'plot_read_classification_summary',
    'plot_truncation_patterns',
    'plot_all_truncation_patterns',
    'plot_sequence_logo',
    # ted_core
    'tss_tts_metrics',
    'calculate_ted_metrics',
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
]
