#!/usr/bin/env bash
# Builds the truth sets for the -profile test dataset (CONTRACTS section 12, A36): WTC11
# orthogonal end calls cut to the test window (hg38 chr19:1-2,000,000; the test genome keeps hg38
# coordinates, so no offset is applied), plus the truth sheet tests/data/test_truth.csv. The test
# reads are WTC11 (LRGASP), so these are the test samples' own cell line. build_test_data.sh
# calls it; it can also be run on its own from the repo root:
#
#   bash tests/data/build_test_truth.sh [OUT_DIR]      # OUT_DIR defaults to tests/data
#
# Sources: TRUTH_DIR/<assay>/WTC11.<assay>.truth.bed, provisional truth sets (no final quality
# gate yet) built from public data:
#   CAGE    GEO GSE185917, biological replicate 1 (SRA SRR16340036-SRR16340038)
#   PROcap  ENCODE ENCSR178XKL (PINTS calls of the ENCODE PRO-cap pipeline)
#   dRNA    ENCODE LRGASP WTC11 direct RNA (ONT RNA002)
# Each truth BED is one row per call: the call's 1-bp mode point in columns 2-3, strand in
# column 6, a "#" header line. Calls are kept when they lie wholly inside the window; rows are
# copied unchanged. Output is deterministic (sorted, no timestamps).
set -euo pipefail

OUT=${1:-tests/data}
TRUTH_DIR=${TRUTH_DIR:-/private/warm-archive/brookslab/hdheath/ted/23_10_26_ortho_truth/work/WTC11}
ASSAYS=${ASSAYS:-"CAGE:tss PROcap:tss dRNA:tes"}
HUMAN_CHROM=${HUMAN_CHROM:-chr19}
HUMAN_END=${HUMAN_END:-2000000}

H=$OUT/human_region
mkdir -p "$H"
SHEET=$OUT/test_truth.csv
echo "sample,end,assay,path,provisional" > "$SHEET"
for ae in $ASSAYS; do
  a=${ae%%:*} e=${ae##*:}
  src=$TRUTH_DIR/$a/WTC11.$a.truth.bed
  [ -s "$src" ] || { echo "ERROR: missing source: $src" >&2; exit 1; }
  dst=$H/human_region.$a.truth.bed
  { grep -m1 '^#' "$src"
    awk -F'\t' -v OFS='\t' -v c="$HUMAN_CHROM" -v e="$HUMAN_END" '$1 == c && $2 >= 0 && $3 <= e' "$src" \
      | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k3,3n -k6,6 -k4,4
  } > "$dst"
  echo "  $(basename "$dst"): $(grep -vc '^#' "$dst") calls"
  for s in human_region_pacbio human_region_ont; do
    echo "$s,$e,$a,human_region/$(basename "$dst"),true" >> "$SHEET"
  done
done
echo "  $(basename "$SHEET"): $(($(wc -l < "$SHEET") - 1)) rows"
