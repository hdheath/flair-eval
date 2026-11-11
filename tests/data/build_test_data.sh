#!/usr/bin/env bash
# Rebuilds the -profile test dataset (human region) in tests/data/ from local copies of public data.
# Provenance and licences: docs/test_data.md. You only need this to regenerate the data;
# the outputs are committed. Run from the repo root:
#
#   bash tests/data/build_test_data.sh [OUT_DIR]      # OUT_DIR defaults to tests/data
#
# Needs samtools (>=1.12, for -e endpos), python3 and gzip. Every source path can be
# overridden with an env var of the same name. Output is deterministic: reads are picked by
# md5(SEED:read_name), and gzip runs with -n so md5s don't change between rebuilds.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 XGBOOST_NTHREAD=1

OUT=${1:-tests/data}
N_READS=${N_READS:-20000}
SEED=${SEED:-flair-eval-testdata-v1}

# --- sources (local copies of public files; see docs/test_data.md) -------------------------
LRGASP=${LRGASP:-/private/warm-archive/brookslab/hdheath/04_08_26_lrgasp_wtc11_ont}
GRCH38=${GRCH38:-/private/groups/brookslab/reference_sequence/GRCh38.primary_assembly.genome.fa}
GENCODE=${GENCODE:-/private/groups/brookslab/reference_annotations/gencode.v49.basic.annotation.gtf}
# ONT: ENCODE ENCSR539ZXJ rep1 (ENCFF263YFG), minimap2 -ax splice to GRCh38 primary
ONT_HUMAN_BAM=${ONT_HUMAN_BAM:-$LRGASP/aligned/WTC11_ONT_cDNA_rep1.bam}
# PacBio: ENCODE ENCSR507JOF reps 1-3 (ENCFF563QZR, ENCFF370NFS, ENCFF245IPA; read names
# carry the ENCFF prefix), minimap2 -ax splice to GRCh38 primary (+ spike-in contigs, unused here)
PB_BAM=${PB_BAM:-$LRGASP/spikein/aligned/pacbio.spikeingenome.bam}

HUMAN_CHROM=chr19
HUMAN_END=2000000                       # hg38 chr19:1-2,000,000, kept in hg38 coordinates

die(){ echo "ERROR: $*" >&2; exit 1; }
log(){ echo "[$(date '+%H:%M:%S')] $*"; }
for f in "$GRCH38" "$GRCH38.fai" "$GENCODE" "$ONT_HUMAN_BAM" "$PB_BAM"; do
  [ -s "$f" ] || die "missing source: $f"
done
command -v samtools >/dev/null || die "samtools not on PATH"

H=$OUT/human_region
mkdir -p "$H"
TMP=$(mktemp -d "${TMPDIR:-/tmp}/flair_eval_testdata.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

# pick N names deterministically: rank by md5(SEED:name), keep the N smallest
pick_names(){  # stdin: read names; $1: N
  python3 -I -c '
import hashlib, sys
seed, n = sys.argv[1], int(sys.argv[2])
names = sorted(set(l.strip() for l in sys.stdin if l.strip()))
key = lambda x: hashlib.md5(f"{seed}:{x}".encode()).hexdigest()
sys.stdout.write("".join(x + "\n" for x in sorted(sorted(names, key=key)[:n])))
' "$SEED" "$1"
}

# subsample primary alignments (-F 0x904) from BAM $1 over regions $3.. into gz FASTA $2.
# Reads keep their BAM (coordinate) order; FASTA has no quality strings by design.
subsample(){
  local bam=$1 out=$2; shift 2
  local filt=(-F 0x904)
  [ -n "${ENDPOS_MAX:-}" ] && filt+=(-e "endpos <= $ENDPOS_MAX")
  samtools view "${filt[@]}" "$bam" "$@" | cut -f1 > "$TMP/all.txt"
  local avail; avail=$(sort -u "$TMP/all.txt" | wc -l)
  pick_names "$N_READS" < "$TMP/all.txt" > "$TMP/keep.txt"
  samtools view -b "${filt[@]}" -N "$TMP/keep.txt" "$bam" "$@" \
    | samtools fasta - 2>/dev/null | gzip -n -9 > "$out"
  log "  $(basename "$out"): $(zcat "$out" | grep -c '^>') reads (of $avail available)"
}

# --- human region ------------------------------------------------------------------------
log "human region $HUMAN_CHROM:1-$HUMAN_END"
samtools faidx "$GRCH38" "$HUMAN_CHROM:1-$HUMAN_END" | sed "1s/.*/>$HUMAN_CHROM/" > "$H/human_region.genome.fa"
samtools faidx "$H/human_region.genome.fa"
printf '%s\t0\t%s\n' "$HUMAN_CHROM" "$HUMAN_END" > "$H/human_region.regions.bed"

# every line of every GENCODE gene that lies wholly inside the window
awk -F'\t' -v c="$HUMAN_CHROM" -v e="$HUMAN_END" '
  /^#/ { print; next }
  $1 != c { next }
  { match($9, /gene_id "[^"]+"/); g = substr($9, RSTART, RLENGTH) }
  $3 == "gene" { keep[g] = ($4 >= 1 && $5 <= e) }
  keep[g] { print }' "$GENCODE" > "$H/human_region.gtf"

ENDPOS_MAX=$HUMAN_END subsample "$ONT_HUMAN_BAM" "$H/human_region.ont.fa.gz" "$HUMAN_CHROM:1-$HUMAN_END"
ENDPOS_MAX=$HUMAN_END subsample "$PB_BAM"        "$H/human_region.pacbio.fa.gz" "$HUMAN_CHROM:1-$HUMAN_END"

# --- end references + TUSCO (CONTRACTS A26), cut to the same window ---------------------
HUMAN_CHROM=$HUMAN_CHROM HUMAN_END=$HUMAN_END bash "$(dirname "$0")/build_test_refs.sh" "$OUT"

# --- checksums -----------------------------------------------------------------------------
log "md5sums"
( cd "$OUT" && find human_region -type f ! -name '*.md5' | LC_ALL=C sort | xargs md5sum > md5sums.txt )
log "done -> $OUT"
ls -l "$H"
