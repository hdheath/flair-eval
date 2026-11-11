#!/usr/bin/env bash
# Builds the end references and the TUSCO list for the -profile test dataset (CONTRACTS A26):
# the lab's hg38 reference files cut to the test window (hg38 chr19:1-2,000,000; the test
# genome keeps hg38 coordinates, so no offset is applied). build_test_data.sh calls it; it can
# also be run on its own from the repo root:
#
#   bash tests/data/build_test_refs.sh [OUT_DIR]      # OUT_DIR defaults to tests/data
#
# Sources (override each with an env var of the same name):
#   REFTSS   refTSS v4.1 human hg38 (BED6)
#   POLYADB  PolyA_DB v4.1 HumanPas main sites, hg38, BED6 [pos-1,pos)
#   TUSCO    SQANTI3 tusco_human.tsv (commit 589dfd4), copied whole (46 genes; none lies in the
#            window, so sections.tusco has genes_in_regions = 0 on the test data)
# Each source folder has a SOURCE.tsv with provenance. Sites are kept when they lie wholly
# inside the window. Output is deterministic (sorted, no timestamps).
set -euo pipefail

OUT=${1:-tests/data}
REF=/private/groups/brookslab/flair-eval/reference
REFTSS=${REFTSS:-$REF/end_refs/refTSS_v4.1.hg38.bed}
POLYADB=${POLYADB:-$REF/end_refs/PolyA_DB_v4.1.main.hg38.bed}
TUSCO=${TUSCO:-$REF/tusco/tusco_human.tsv}
HUMAN_CHROM=${HUMAN_CHROM:-chr19}
HUMAN_END=${HUMAN_END:-2000000}

for f in "$REFTSS" "$POLYADB" "$TUSCO"; do
  [ -s "$f" ] || { echo "ERROR: missing source: $f" >&2; exit 1; }
done
H=$OUT/human_region
mkdir -p "$H"

cut_window(){  # $1 source BED, $2 output BED
  awk -F'\t' -v OFS='\t' -v c="$HUMAN_CHROM" -v e="$HUMAN_END" '$1 == c && $2 >= 0 && $3 <= e' "$1" \
    | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k3,3n -k6,6 -k4,4 > "$2"
  echo "  $(basename "$2"): $(wc -l < "$2") sites"
}
cut_window "$REFTSS"  "$H/human_region.refTSS.bed"
cut_window "$POLYADB" "$H/human_region.polyA_DB.bed"
cp "$TUSCO" "$H/human_region.tusco.tsv"
echo "  human_region.tusco.tsv: $(grep -vc '^#' "$H/human_region.tusco.tsv") genes"
