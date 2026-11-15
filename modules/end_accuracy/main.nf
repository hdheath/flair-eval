/*
 * END_ACCURACY: one FLAIR run's isoform ends against its sample's end sets (CONTRACTS section 12).
 *
 * End sets = the sample's truth sets from --truth (kind truth) plus the end references
 * --ref_tss / --ref_pas or the genome preset's (kind reference). bin/end_accuracy.py runs
 * flair4_ted's judge, ted_end_precision.py, vendored into resources/usr/bin/ with its helpers and
 * extended with per-chain columns after the legacy ones (provenance: README.md here), and
 * writes end_accuracy/end_accuracy.json (sections.end_accuracy, merged by BUNDLE) and
 * end_accuracy/tables/end_accuracy.tsv.
 *
 * Runs only for jobs whose sample has at least one end set. Pass [] for a missing reference.
 */
process END_ACCURACY {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(assemble), path(gtf), path(truth, stageAs: 'truth/?/*'), val(info), path(ref_tss, stageAs: 'refs/tss/*'), path(ref_pas, stageAs: 'refs/pas/*')
    path(regions)

    output:
    tuple val(meta), path('end_accuracy', type: 'dir'), emit: dir

    script:
    def sets = endSets(truth, info, ref_tss, ref_pas)
    """
    cat > sets.json <<'JSON'
${groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(sets))}
JSON
    end_accuracy.py \\
        --isoforms-bed ${assemble}/isoforms.bed \\
        --gtf ${gtf} \\
        --sets sets.json \\
        ${regions ? "--regions ${regions}" : ''} \\
        --out end_accuracy
    """

    stub:
    def sets = endSets(truth, info, ref_tss, ref_pas)
    def section = [
        scorer     : [script: 'ted_end_precision.py', source: 'stub', md5: null, window: 50, regions: null],
        provisional: sets.any { s -> s.kind == 'truth' && s.provisional },
        sets       : sets.collect { s -> [id: s.id, kind: s.kind, end: s.end, assay: s.assay, file: s.source.tokenize('/')[-1], path: s.source, md5: null, provisional: s.provisional] },
        ends       : sets.collect { s ->
            [set: s.id, kind: s.kind, end: s.end, assay: s.assay, isoforms: 0, isoform_ends: 0, sites: 0,
             sites_hit: s.kind == 'truth' ? 0 : null, precision: null, naive_precision: null, dedup_tp: 0, redundant_calls: 0,
             recall: null, f1: null]
        },
        pairs      : stubPairs(sets),
    ]
    """
    end_accuracy.py --help > /dev/null
    mkdir -p end_accuracy/tables
    cat > end_accuracy/end_accuracy.json <<'JSON'
${groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(section))}
JSON
    printf 'kind\\tend\\tset\\tassay\\tpartner\\tisoforms\\tisoform_ends\\tsites\\tsites_hit\\tprecision\\tnaive_precision\\tdedup_tp\\tredundant_calls\\trecall\\tf1\\tunique_pairs\\tboth_hit\\tboth_hit_raw\\tpair_precision\\tpair_naive_precision\\n' > end_accuracy/tables/end_accuracy.tsv
    """
}

/*
 * The sets.json list for bin/end_accuracy.py: references first (ids refTSS / PolyA_DB), then the
 * truth sets in sheet order. `truth` is the staged file list (a single path when there is one),
 * in the same order as info.truth ([assay, end, source, provisional] maps); info.ref_tss /
 * info.ref_pas are the references' original paths, recorded in the section.
 */
def endSets(Object truth, Map info, Object refTss, Object refPas) {
    def staged = truth instanceof List ? truth : (truth ? [truth] : [])
    def sets = []
    if (refTss) {
        sets << [id: 'refTSS', kind: 'reference', end: 'tss', assay: 'refTSS', path: refTss.toString(), source: info.ref_tss ?: refTss.toString(), provisional: false]
    }
    if (refPas) {
        sets << [id: 'PolyA_DB', kind: 'reference', end: 'tes', assay: 'PolyA_DB', path: refPas.toString(), source: info.ref_pas ?: refPas.toString(), provisional: false]
    }
    (info.truth ?: []).eachWithIndex { t, i ->
        sets << [id: t.assay, kind: 'truth', end: t.end, assay: t.assay, path: staged[i].toString(), source: t.source, provisional: t.provisional]
    }
    return sets
}

def stubPairs(List sets) {
    def out = []
    ['reference', 'truth'].each { k ->
        def t5 = sets.findAll { s -> s.kind == k && s.end == 'tss' }
        def t3 = sets.findAll { s -> s.kind == k && s.end == 'tes' }
        t5.each { a ->
            t3.each { b ->
                out << [kind: k, tss_set: a.id, tes_set: b.id, isoforms: 0, unique_pairs: 0, both_hit: 0, both_hit_raw: 0, precision: null, naive_precision: null]
            }
        }
    }
    return out
}
