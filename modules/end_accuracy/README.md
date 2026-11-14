# END_ACCURACY

Scores one FLAIR run's isoform 5′ and 3′ ends against its sample's end sets: the truth sets
from `--truth`, plus the refTSS / PolyA_DB end references (CONTRACTS section 12, A36–A40).

`bin/end_accuracy.py` scores nothing itself. It plans the judge calls, runs them, reads each
call's one-row `precision_recall_summary.tsv`, and shapes `sections.end_accuracy`. Every number
comes from the judge below.

## The judge (vendored, unchanged)

`resources/usr/bin/` holds flair4_ted's end judge and its two helpers. They were copied byte for
byte from `legacy/bin/evaluation/` (flair4_ted `bin/evaluation/`, commit `40fdd21`). The files
are not changed, and `tests/python/test_end_accuracy.py` checks that they still match `legacy/`.

| file | md5 |
|---|---|
| `ted_end_precision.py` | `71b8835195ec7d38ea0a1e9b45f378d2` |
| `peak_matching.py` | `655a40d17bd8ac7dac23c5ce48695436` |
| `signal_utils.py` | `a9466e1b59b564cccefe0742cf279723` |

To update them, copy the new files from flair4_ted, then update this table, the md5s in the test,
and `JUDGE_SOURCE` in `bin/flair_eval_lib/end_accuracy.py`. Never edit them here.

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
  never covered. Recall is therefore a lower bound for a shallow library.
- **`dedup_tp` / `redundant_calls`.** These come from the junction-chain grouping. They are not
  the numerator of the headline precision (which is many-to-one). They are kept in the table as
  diagnostics only.
- **Duplicate database rows.** PolyA_DB has duplicate rows that inflate its site count. This
  matters only for recall, which is not shown for the databases.
- **Without a 3′ set.** The judge then reports GTF-based 3′ numbers. `end_accuracy.py` reads only
  the ends that each call was given sites for.
- **Unstable ids.** `per_junction_chain.tsv` keys chains with Python `hash()`, which changes from
  run to run. The file is kept for the record, not for joins.
