#!/bin/bash
#SBATCH --job-name=ted_pileup_tracks
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --output=/private/home/hdheath/projects/flair4_ted/bin/debug_work/pileup_tracks_%j.log
#SBATCH --error=/private/home/hdheath/projects/flair4_ted/bin/debug_work/pileup_tracks_%j.err
# ===========================================================================
# Genome-wide TSS/TTS read-end pileup tracks (IGV bigWig) for the 3 A549 BAMs,
# reproducing TED's global-peak-caller pileup (_build_pileup_coverage):
#   TSS = where transcripts START, TTS = where they END, strand-aware, RAW counts.
# Per BAM -> 4 bigWigs: {sample}.{tss,tts}_{fwd,rev}.bw
#
# Single per-chrom streaming pass over each BAM (index-based fetch; the dRNA BAM
# has a .csi — pysam handles it). bedGraph + sort on node-local /data/tmp (fast);
# only the small bigWigs land on shared storage. m6A MM/ML tags are never touched
# (we read alignment coords only) so the 23 GB dRNA pass stays cheap.
# ===========================================================================
set -euo pipefail
PY=/private/home/hdheath/miniforge3/envs/flair-ted-52826/bin/python
SCRIPT=/private/home/hdheath/projects/flair4_ted/bin/ted_pileup_tracks.py
FAI=/private/groups/brookslab/reference_sequence/GRCh38.primary_assembly.genome.fa.fai
# save to a plain groups dir (NOT prism-published copyable_data — the smoothed
# genome-wide tracks can be large); access via the groups mount / scp to IGV.
OUT=/private/groups/brookslab/hdheath/ted_pileup_tracks
mkdir -p "$OUT" /private/home/hdheath/projects/flair4_ted/bin/debug_work

LOCAL=/data/tmp/${USER}_${SLURM_JOB_ID:-manual}; mkdir -p "$LOCAL"
trap 'rm -rf "$LOCAL"' EXIT

run() {  # $1=sample  $2=bam
  local sample=$1 bam=$2
  echo "[$(date +%T)] === $sample ===  $bam"
  "$PY" "$SCRIPT" --bam "$bam" --fai "$FAI" --prefix "$LOCAL/$sample" \
    --smooth --sigma-tss 75 --sigma-tts 15
  cp "$LOCAL/$sample".*.bw "$OUT/"
  echo "[$(date +%T)] $sample bigWigs -> $OUT/"
}

run A549_pacbio_cDNA /private/groups/brookslab/hdheath/projects/flair-eval/data/D_2411_PB_CL_Lu_N_A549.pbmm2.sorted.bam
run A549_ONT_cDNA    /private/groups/brookslab/data.rep/A549_PCS114.sorted.bam
run A549_dRNA        /private/groups/brookslab/brett/260127_A549_dRNA/A549_dRNA_m6A_aln.bam

echo "=== done $(date '+%F %H:%M') ==="
ls -la "$OUT"/*.bw 2>/dev/null
echo "tracks saved to (groups; scp to IGV or load from the mount): $OUT/"
echo "  RAW    = .{tss,tts}_{fwd,rev}.bw           (1-bp read-end counts)"
echo "  SMOOTH = .{tss,tts}_{fwd,rev}.smooth.bw    (gaussian σ_tss=75 / σ_tts=15, the find_peaks input)"
