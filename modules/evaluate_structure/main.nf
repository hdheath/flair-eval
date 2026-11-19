/*
 * EVALUATE_STRUCTURE: isoform structure vs the reference GTF (bin/structure_metrics.py, owned by
 * the report stream). Writes bundle.json + tables/ into eval/ (CONTRACTS section 3).
 * Always gets the genome (junction motifs, A13); reads/BAM paths are only recorded in inputs.
 * Round 3 (A26-A32): the end references (refTSS, PolyA_DB BED6) and the TUSCO list are optional;
 * pass [] for a missing one and only its section is left out. Writes tables/ends.tsv and
 * sections.ends / tusco / priming / variation.
 * Section 13: sections.structure_qc always; the sample's `junctions` file (STAR
 * SJ.out.tab or BED) is optional and adds structure_qc.junction_support. A45: the sample's
 * `sr_junctions` file (short reads, scoring only) is optional and adds sr_junction_support.
 */
process EVALUATE_STRUCTURE {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(assemble), path(gtf), path(genome), path(genome_fai), path(env_json), path(ref_tss, stageAs: 'refs/tss/*'), path(ref_pas, stageAs: 'refs/pas/*'), path(tusco, stageAs: 'refs/tusco/*'), path(junction_file, stageAs: 'refs/sj/*'), path(sr_junction_file, stageAs: 'refs/srsj/*'), val(inputs)
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
        --genome ${genome} \\
        ${inputs.reads ? "--reads '${inputs.reads}'" : ''} \\
        ${inputs.bam ? "--bam '${inputs.bam}'" : ''} \\
        ${regions ? "--regions ${regions}" : ''} \\
        ${ref_tss ? "--ref-tss ${ref_tss}" : ''} \\
        ${ref_pas ? "--ref-pas ${ref_pas}" : ''} \\
        ${tusco ? "--tusco ${tusco}" : ''} \\
        ${junction_file ? "--junction-file ${junction_file}" : ''} \\
        ${sr_junction_file ? "--sr-junction-file ${sr_junction_file}" : ''} \\
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
        sections         : stubSections(ref_tss, ref_pas, tusco, junction_file, sr_junction_file),
    ]
    """
    structure_metrics.py --help > /dev/null
    mkdir -p eval/tables
    cat > eval/bundle.json <<'JSON'
${groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(bundle))}
JSON
    printf 'isoform_id\\tgene_id\\tchrom\\tstart\\tend\\tstrand\\tn_exons\\tlength\\tcategory\\tsubcategory\\tref_transcript\\tref_gene\\n' > eval/tables/isoforms.tsv
    printf 'iso1_stub\\tnovel_gene_1\\tchr1\\t100\\t900\\t+\\t2\\t500\\tintergenic\\t\\t\\t\\n' >> eval/tables/isoforms.tsv
    printf 'iso2_stub\\tg1\\tchr1\\t1\\t20\\t+\\t1\\t20\\tFSM\\tCIM\\tt1\\tg1\\n' >> eval/tables/isoforms.tsv
    printf 'gene_id\\tgene_name\\ng1\\tGENE1\\n' > eval/tables/genes.tsv
    printf 'isoform_id\\tcategory\\tsubcategory\\tstrand\\ttss\\ttes\\td_reftss\\td_pas\\td_ref_tss\\td_ref_tes\\tpolya_motif\\tdown_a_frac\\n' > eval/tables/ends.tsv
    printf 'iso1_stub\\tintergenic\\t\\t+\\t100\\t899\\t\\t\\t\\t\\t\\t\\n' >> eval/tables/ends.tsv
    printf 'iso2_stub\\tFSM\\tCIM\\t+\\t1\\t19\\t\\t\\t0\\t0\\t\\t\\n' >> eval/tables/ends.tsv
    printf 'category\\tsubcategory\\tn\\nintergenic\\t\\t1\\n' > eval/tables/categories.tsv
    """
}

def stubZeros(List keys) {
    return keys.collectEntries { k -> [(k): 0] }
}

def stubHist(List edges) {
    return [edges: edges, counts: edges.collect { _e -> 0 }]
}

// Round-3 sections with zeros and nulls (A29-A32), so stub bundles have the real shape. ends and
// tusco appear only when their files were staged, as in the real script.
// Parser v2 cannot call a closure held in a local variable as name(), so these are functions.
def stubDistHist() {
    def distEdges = [-1000, -200, -50, -10, 0, 10, 50, 200, 1000]
    return [edges: distEdges, counts: (0..distEdges.size()).collect { _i -> 0 }]
}

def stubEndHalf() {
    return [n: 0, within: [0, 0, 0, 0, 0, 0, 0, 0], dist_hist: stubDistHist(), no_site: 0]
}

def stubRefMeta(Object f) {
    return [file: f.toString().tokenize('/')[-1], md5: null, sites: 0]
}

def stubSections(Object refTss, Object refPas, Object tusco, Object junctionFile, Object srJunctionFile) {
    def cats = ['FSM', 'ISM', 'NIC', 'NNC', 'genic', 'genic_intron', 'antisense', 'intergenic', 'fusion']
    def sec = [:]
    if (refTss || refPas) {
        sec.ends = [
            refs        : [tss: refTss ? stubRefMeta(refTss) : null, pas: refPas ? stubRefMeta(refPas) : null],
            windows     : [0, 10, 25, 50, 100, 200, 300, 500],
            tss         : refTss ? stubEndHalf() : null,
            tes         : refPas ? stubEndHalf() : null,
            fsm_tss_hist: stubDistHist(),
            fsm_tes_hist: stubDistHist(),
            by_category : cats.collectEntries { c -> [(c): [n: 0, tss_within_50: refTss ? 0 : null, tes_within_50: refPas ? 0 : null]] },
            polya_motif : [n_with_motif: 0, n: 0],
        ]
    }
    if (tusco) {
        sec.tusco = [file: tusco.toString().tokenize('/')[-1], md5: null, id_type: 'ensembl', genes_total: 0, genes_in_regions: 0,
                     TP: 0, PTP: 0, FP: 0, FN: 0, sensitivity: null, precision: null, precision_redundant: null,
                     positive_detection_rate: null, transcripts: [n: 0, TP: 0, PTP: 0, FP: 0], per_gene: [:], gene_names: [:]]
    }
    sec.priming = [n: 0, risk: 0, down_a_hist: [edges: (0..10).collect { i -> i / 10 }, counts: (0..10).collect { _i -> 0 }]]
    sec.variation = [genes_multi: 0, ends_only: 0, splice_only: 0, both: 0, neither: 0]
    sec.structure_qc = stubStructureQc(cats, junctionFile, srJunctionFile)
    return sec
}

// Section 13 with zeros (bin/flair_eval_lib/structure_qc.py); the genome is always staged, so
// canonical, rts and priming are always there; junction_support only with the sample's file,
// sr_junction_support only with its sr_junctions file (A45).
def stubZeroList(List edges) {
    return edges.collect { _e -> 0 }
}

def stubStructureQc(List cats, Object junctionFile, Object srJunctionFile) {
    def lenEdges = [0, 250, 500, 750, 1000, 1500, 2000, 3000, 4000, 6000, 10000]
    def exEdges = [1, 2, 3, 4, 5, 6, 8, 11, 16, 21, 31]
    def redEdges = [1, 2, 3, 4, 5, 6, 8, 11]
    def jt = ['known_canonical', 'known_noncanonical', 'novel_canonical', 'novel_noncanonical']
    def qc = [
        subcategories: [
            ISM: stubZeros(['3prime_fragment', '5prime_fragment', 'internal_fragment', 'intron_retention', 'mono-exon']),
            NIC: stubZeros(['combination_of_known_junctions', 'combination_of_known_splicesites', 'intron_retention', 'mono-exon_by_intron_retention']),
            NNC: stubZeros(['at_least_one_novel_splicesite', 'intron_retention']),
        ],
        redundancy   : stubZeros(['fsm_isoforms', 'fsm_refs', 'fsm_refs_multi', 'fsm_ism_isoforms', 'fsm_ism_refs', 'fsm_ism_refs_multi'])
            + [fsm_hist: [edges: redEdges, counts: stubZeroList(redEdges)], fsm_ism_hist: [edges: redEdges, counts: stubZeroList(redEdges)]],
        by_category  : [length: [edges: lenEdges, counts: cats.collectEntries { c -> [(c): stubZeroList(lenEdges)] }],
                        exons : [edges: exEdges, counts: cats.collectEntries { c -> [(c): stubZeroList(exEdges)] }]],
        ref_recovery : [length: [edges: lenEdges, all: stubZeroList(lenEdges), fsm: stubZeroList(lenEdges)],
                        exons : [edges: exEdges, all: stubZeroList(exEdges), fsm: stubZeroList(exEdges)]],
        canonical    : [multi_exon: 0, all_canonical: 0, by_category: [:]],
        rts          : [junctions: 0, rts: 0, by_type: jt.collectEntries { t -> [(t): [n: 0, rts: 0]] }, isoforms: 0, isoforms_with_rts: 0,
                        params: [pattern: 10, wiggle: 1, min_match: 8, mismatches: 1]],
        priming      : [by_category: [:], by_exons: [mono: [n: 0, risk: 0], multi: [n: 0, risk: 0]]],
    ]
    if (junctionFile) {
        qc.junction_support = stubJunctionSupport(junctionFile)
    }
    if (srJunctionFile) {
        qc.sr_junction_support = stubJunctionSupport(srJunctionFile)
    }
    return qc
}

def stubJunctionSupport(Object f) {
    return [file: f.toString().tokenize('/')[-1], md5: null, file_junctions: 0, junctions: 0, supported: 0,
            known: [n: 0, supported: 0], novel: [n: 0, supported: 0], multi_exon: 0, all_supported: 0, by_category: [:]]
}
