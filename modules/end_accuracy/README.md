# END_ACCURACY

Scores one FLAIR run's isoform 5′ and 3′ ends against its sample's end sets: the truth sets
from `--truth`, plus the refTSS / PolyA_DB end references (CONTRACTS section 12, A36–A40).

`bin/end_accuracy.py` scores nothing itself. It plans the judge calls, runs them, reads each
call's one-row `precision_recall_summary.tsv`, and shapes `sections.end_accuracy`. Every number
comes from the judge below.

## The judge (vendored, extended)

`resources/usr/bin/` holds flair4_ted's end judge and its two helpers, copied from
`legacy/bin/evaluation/` (flair4_ted `bin/evaluation/`, commit `40fdd21`). The helpers are
unchanged. The judge was extended here for the TED end clustering effort (ted-clustering map
tickets 05 and 12, 2026-11-15): it gained per-chain columns after its legacy ones, and the
`--read-ends` input. flair4_ted's own copy stays frozen. `tests/python/test_end_accuracy.py`
checks that the helpers still match `legacy/`, and that the judge's legacy columns and
`per_junction_chain.tsv` are exactly the legacy judge's on its hand-built case.

| file | md5 |
|---|---|
| `ted_end_precision.py` | `ecb9f1366c49c1fb04d73091daa56330` (legacy `71b8835195ec7d38ea0a1e9b45f378d2`) |
| `peak_matching.py` | `655a40d17bd8ac7dac23c5ce48695436` |
| `signal_utils.py` | `a9466e1b59b564cccefe0742cf279723` |
| `build_truth_union.py` | new here (see below) |

Change the legacy columns only by re-vendoring from flair4_ted (then update this table, the
test and `JUDGE_SOURCE` in `bin/flair_eval_lib/end_accuracy.py`).

## The per-chain columns (ted-clustering tickets 05 and 12)

Appended to `precision_recall_summary.tsv` after the legacy columns (`CHAIN_FIELDS` in the
judge), for each axis (`5prime`, `3prime`) that has truth peaks; blank otherwise.
`bin/end_accuracy.py` does not read them yet: they are for direct judge calls (the clustering
map's offline arms, whose `<prefix>.variants.bed` goes in as `--isoforms-bed`).

- **One-to-one within each junction chain** (the headline). A chain is a group of
  `per_junction_chain.tsv` (multi-exon isoforms by junction chain; single-exon ones by overlap).
  Its called ends are the distinct end positions of its isoforms. They pair one-to-one with
  truth calls (same strand, 50 bp, the judge's interval predicate). Two chains that share a
  truth call both get credit; two calls of one chain near one truth call count once.
  - `*_chain_ends`, `*_chain_tp`, `*_chain_precision` = tp / ends
  - `*_chain_sites_hit` (distinct truth calls hit by any chain), `*_chain_recall_all_calls`
    = sites hit / every truth call in scope
- **Read-supported recall** (needs `--read-ends`). A truth call is read-supported when at least
  one corrected spliced read of the library ends within 50 bp of it (same strand).
  `*_read_supported_sites`, `*_chain_sites_hit_read_supported`,
  `*_chain_recall_read_supported` = read-supported calls hit / read-supported calls.
  `n_read_ends` is the number of reads read.
- **Alternative-end pairs**, in two bins of distance: `51_200` (51–200 bp) and `over_200`.
  - Truth side (needs `--read-ends`): two read-supported calls on one axis, at least 51 bp
    apart, each with at least one read of the same junction chain ending at it. A read ends at
    one call, its nearest within 50 bp (ties to the lower coordinate), so one read between two
    calls 51–100 bp apart cannot make a pair on its own. `*_altpair_<bin>_truth_pairs`;
    `*_pairs_hit` = pairs where that chain's one-to-one matching hits both calls; `*_recall` =
    hit / pairs.
  - Distances are between summits, which for 1-bp truth calls are the reference points (the
    ortho-truth BEDs, as ticket 05 assumes). On a wide-peak BED a pair 51–200 bp apart by gap
    can be over 200 bp apart by summit and land in `over_200`.
  - Called side: a chain's called ends other than its main end, binned by distance to it.
    The main end has the most reads (summed BED score of the chain's isoforms at that end; 1
    per isoform from a GTF); ties go to the 5′-most end (lowest on +, highest on −), on both
    axes. `*_alt_ends`, `*_alt_ends_hit` (hit in the chain's one-to-one matching),
    `*_precision` = hit / alt ends. The reads are the input BED's scores, not `--read-ends`:
    for the harness's `variants.bed` they are each variant's reads, so the main end is the
    variant end with the most reads, not the 5′ cluster with the most reads.

`--read-ends` takes a genome-record dump (ted-clustering ticket 01's format): a directory of
region TSVs, a `.tgz` of one, or one `.tsv[.gz]` such as the joined `*.genome_records.tsv.gz`.
It uses the columns `chrom`, `strand`, `introns`, `five_prime`, `three_prime`; reads of unknown
strand are skipped, and a dump with no reads left is an error. Every arm of a library is
scored with the same reads.

## Truth unions and atlas sets

`build_truth_union.py` writes one line's 5′ truth union as a truth BED, so a union is one more
`--truth` sheet entry:

    build_truth_union.py --out A549.union5.truth.bed CAGE=A549.CAGE.truth.bed PROcap=A549.PROcap.truth.bed

Sets go in tier order: CAGE, then PRO-cap, then RAMPAGE, for each line with two or more of
them. Every first-tier call is kept at its reference point; each next set's calls pair
one-to-one (sweep) with the calls kept so far within 50 bp, same strand. A paired call is the
same site and is dropped, so a site seen by both sets counts once; an unpaired call is added.
Pairing is one-to-one, so where two next-tier calls lie within 50 bp of one kept call, one is
dropped and the other is added.
Inputs must be 1-bp reference points (the ortho-truth BEDs). There is no 3′ union, and the
atlases never join one.

The atlases (refTSS v4.1, PolyA_DB v4.1) are scored on their own: here they are the
`reference` kind (precision only in the section). A dRNA library is never scored against its
own line's dRNA truth (circular): leave that set off the library's sheet.

## The rule (as the judge applies it)

- **Hit.** An isoform end is a hit when it lies on the same strand within 50 bp of any part of
  a site's interval (`[start − 50, end − 1 + 50]`). Only chr1–22, X and Y are scored. A truth BED
  holds each call's 1-bp mode point, so on those files this means "within 50 bp of the mode".
- **Precision** = isoform ends that are hits / isoform ends.
- **Recall** = sites with a hit near them / all sites in scope. This is shown for truth sets
  only. The databases pool many tissues, so recall against them means nothing.
- Single-end matching is many-to-one: several isoform ends may share one site.
- **Both ends** (`pairs`) = isoforms whose 5′ and 3′ ends both hit, divided by all isoforms.
  Within a junction chain, isoforms on the same (5′ site, 3′ site) pair count once.

## Known quirks (documented, not patched)

- **End position is one past the end.** `tss_tts()` takes the raw BED start/end. So the 5′ end of
  a minus-strand isoform, and the 3′ end of a plus-strand one, are one past the last transcribed
  base. With a 50 bp window this moves almost no calls.
- **Imports.** `signal_utils.py` imports numpy and `scipy.stats.gaussian_kde` at the top, so the
  judge needs scipy even though the end scoring never uses it. That is why the module has its
  own environment (`environment.yml`).
- **Region check.** With `--region`, the overlap test is inclusive at both ends.
- **Recall denominator.** It counts every site in scope, including sites of genes the run's reads
  never covered. Recall is therefore a lower bound for a shallow library. The per-chain
  `*_chain_recall_read_supported` column divides by read-supported sites instead.
- **`dedup_tp` / `redundant_calls`.** These come from the junction-chain grouping. They are not
  the numerator of the headline precision (which is many-to-one). They are kept in the table as
  diagnostics only.
- **Duplicate database rows.** PolyA_DB has duplicate rows that inflate its site count. This
  matters only for recall, which is not shown for the databases.
- **Without a 3′ set.** The judge then reports GTF-based 3′ numbers. `end_accuracy.py` reads only
  the ends that each call was given sites for.
- **Unstable ids.** `per_junction_chain.tsv` keys chains with Python `hash()` of the junction
  tuple. The file is kept for the record, not for joins.
