/*
 * EVALUATE_STRUCTURE: isoform structure vs the reference GTF (bin/structure_metrics.py, owned by
 * the report stream). Writes bundle.json + tables/ into eval/ (CONTRACTS section 3).
 */
process EVALUATE_STRUCTURE {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(assemble), path(gtf), path(env_json)
    path(regions)

    output:
    tuple val(meta), path('eval', type: 'dir'), emit: dir

    script:
    """
    cat > run.json <<'JSON'
${groovy.json.JsonOutput.toJson(meta.run)}
JSON
    structure_metrics.py \\
        --isoforms-bed ${assemble}/isoforms.bed \\
        --isoforms-gtf ${assemble}/isoforms.gtf \\
        --gtf ${gtf} \\
        --sample ${meta.sample} \\
        --run-json run.json \\
        --env-json ${env_json} \\
        ${regions ? "--regions ${regions}" : ''} \\
        --out eval
    """

    stub:
    def bundle = [
        schema           : 'flair-eval/bundle/1',
        sample           : meta.sample,
        run_id           : meta.run_id,
        label            : meta.run.label,
        tags             : meta.run.tags,
        flair            : [source: 'release', requested: 'stub', version: 'stub', commit: null],
        summary          : stubZeros(['isoforms', 'genes', 'loci', 'mono_exon', 'multi_exon', 'ref_transcripts_total', 'ref_transcripts_fsm', 'ref_transcripts_cim', 'ref_genes_detected']) + [isoforms: 1, genes: 1, loci: 1, multi_exon: 1, ref_transcripts_total: 1],
        categories       : stubZeros(['FSM', 'ISM', 'NIC', 'NNC', 'genic', 'genic_intron', 'antisense', 'intergenic', 'fusion']) + [intergenic: 1],
        fsm_subcategories: stubZeros(['CIM', 'FSM5', 'FSM3', 'FSMTE']),
        junctions        : stubZeros(['total', 'known', 'novel_canonical', 'novel_noncanonical']),
        histograms       : [
            isoforms_per_gene: stubHist([1, 2, 3, 4, 5, 6, 8, 11, 16, 21]),
            transcript_length: stubHist([0, 250, 500, 750, 1000, 1500, 2000, 3000, 4000, 6000, 10000]),
            exon_length      : stubHist([0, 50, 100, 150, 200, 300, 500, 1000, 2000, 5000]),
            exons_per_isoform: stubHist([1, 2, 3, 4, 5, 6, 8, 11, 16, 21, 31]),
        ],
        sections         : [:],
    ]
    """
    structure_metrics.py --help > /dev/null
    mkdir -p eval/tables
    cat > eval/bundle.json <<'JSON'
${groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(bundle))}
JSON
    printf 'isoform_id\\tgene_id\\tchrom\\tstart\\tend\\tstrand\\tn_exons\\tlength\\tcategory\\tsubcategory\\tref_transcript\\tref_gene\\n' > eval/tables/isoforms.tsv
    printf 'iso1_stub\\tnovel_gene_1\\tchr1\\t100\\t900\\t+\\t2\\t500\\tintergenic\\t\\t\\t\\n' >> eval/tables/isoforms.tsv
    printf 'category\\tsubcategory\\tn\\nintergenic\\t\\t1\\n' > eval/tables/categories.tsv
    """
}

def stubZeros(List keys) {
    return keys.collectEntries { k -> [(k): 0] }
}

def stubHist(List edges) {
    return [edges: edges, counts: edges.collect { _e -> 0 }]
}
