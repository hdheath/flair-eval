/*
 * PREP_REFERENCE: with --regions, cut each distinct (genome, GTF) down to the chromosomes in the
 * regions BED (whole chromosomes, so coordinates are unchanged). ASSEMBLE, QUANTIFY,
 * PRODUCTIVITY and EVALUATE_STRUCTURE use the cut copy; alignment keeps the full genome so
 * reads from elsewhere are not forced onto the kept chromosomes. bundle.json "inputs" still
 * records the original paths.
 */
process PREP_REFERENCE {
    tag "${key}"
    label 'process_single'

    input:
    tuple val(key), path(genome), path(gtf)
    path(regions)

    output:
    tuple val(key), path('ref/genome.fa'), path('ref/genome.fa.fai'), path('ref/annotation.gtf'), emit: ref

    script:
    """
    flair_eval_subset_reference.py --genome ${genome} --gtf ${gtf} --regions ${regions} --out ref
    """

    stub:
    """
    flair_eval_subset_reference.py --genome ${genome} --gtf ${gtf} --regions ${regions} --out ref
    """
}
