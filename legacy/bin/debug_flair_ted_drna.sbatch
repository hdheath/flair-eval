#!/bin/bash
#SBATCH --job-name=flair_ted_drna_debug
#SBATCH --partition=medium
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --output=/private/home/hdheath/projects/flair4_ted/bin/logs/flair_ted_drna_debug_%j.log
#SBATCH --error=/private/home/hdheath/projects/flair4_ted/bin/logs/flair_ted_drna_debug_%j.err
#SBATCH --no-requeue
#
# ============================================================================
# Reproduce the A549_dRNA TED FlairTranscriptome failure with --log-debug
# so we get a proper traceback (the Nextflow run exited with a bare
# FileNotFoundError-looking message and "Specify --log-debug for details").
#
# Pipeline task reproduced:
#   ASSEMBLE_AND_EVAL:FlairTranscriptome (A549_dRNA_pre-aligned_All_TED_transcriptome)
#   Original work dir:
#     /private/groups/brookslab/hdheath/results/flair4/work/63/172150cb237dd0ef9fa7b44048fb26
#   Original SLURM job: 32514117 (failed after 2h28m, log stopped writing at 24 min)
#
# Inputs are symlinked in from the cached FlairPartition output dir so we
# don't rebuild the partitioned BAM/GTF/genome.
#
# Outputs + SLURM logs go to:
#   /private/home/hdheath/projects/flair4_ted/bin/logs/
#   /private/home/hdheath/projects/flair4_ted/bin/debug_work/
# ============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
# shellcheck disable=SC1091
source /private/home/hdheath/miniforge3/etc/profile.d/conda.sh
conda activate /private/home/hdheath/miniforge3/envs/flair4

export OMP_NUM_THREADS=1

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PART_DIR=/private/groups/brookslab/hdheath/results/flair4/work/45/3cd0ffd97422133ce7995077073daf
JUNCTIONS=/private/groups/brookslab/hdheath/tools/intronProspector_results/A549_dRNA_m6A_aln.junctions.bed

TSS_MODEL=/private/groups/brookslab/hdheath/projects/flAIr/model_weights/tss/xgboost.pkl
TTS_MODEL=/private/groups/brookslab/hdheath/projects/flAIr/model_weights/tts/xgboost.pkl

SCRATCH=/data/tmp/flair_ted_drna_debug_${SLURM_JOB_ID:-manual}
OUT_DIR=/private/home/hdheath/projects/flair4_ted/bin/debug_work/drna_ted_${SLURM_JOB_ID:-manual}

mkdir -p "$SCRATCH" "$OUT_DIR" /private/home/hdheath/projects/flair4_ted/bin/logs
cd "$SCRATCH"

# ---------------------------------------------------------------------------
# Stage inputs as symlinks (same files Nextflow staged)
# ---------------------------------------------------------------------------
ln -sf "$PART_DIR/A549_dRNA_pre-aligned_All.bam"                A549_dRNA_pre-aligned_All.bam
ln -sf "$PART_DIR/A549_dRNA_pre-aligned_All.bam.bai"            A549_dRNA_pre-aligned_All.bam.bai
ln -sf "$PART_DIR/A549_dRNA_pre-aligned_All_genome.fa"          A549_dRNA_pre-aligned_All_genome.fa
ln -sf "$PART_DIR/A549_dRNA_pre-aligned_All_annotation.gtf"     A549_dRNA_pre-aligned_All_annotation.gtf
ln -sf "$JUNCTIONS"                                             A549_dRNA_m6A_aln.junctions.bed

# ---------------------------------------------------------------------------
# Run FLAIR TED with --log-debug (global option before the module name)
# ---------------------------------------------------------------------------
echo "### Starting at $(date) on $(hostname) in $(pwd)"

flair --log-debug transcriptome \
    -b A549_dRNA_pre-aligned_All.bam \
    --genome A549_dRNA_pre-aligned_All_genome.fa \
    -f A549_dRNA_pre-aligned_All_annotation.gtf \
    -t ${SLURM_CPUS_PER_TASK:-12} \
    --keep_intermediate \
    --junction_bed A549_dRNA_m6A_aln.junctions.bed \
    --ted_log \
    --ted \
    --tss_model "$TSS_MODEL" \
    --tts_model "$TTS_MODEL" \
    --ted_w_depth_tss 0.18 --ted_w_model_tss 0.44 --ted_w_annot_tss 0.38 \
    --ted_threshold_tss 0.65 \
    --ted_w_depth_tts 0.42 --ted_w_model_tts 0.21 --ted_w_annot_tts 0.38 \
    --ted_threshold_tts 0.55 \
    --ted_tss_slack 1.5 \
    --trust_ends \
    --frac_support 0.00 \
    --ted_tss_snap 15 \
    -o A549_dRNA_pre-aligned_All_TED_transcriptome
rc=$?
echo "### FLAIR exited with $rc at $(date)"

# ---------------------------------------------------------------------------
# Copy back outputs so they persist after the scratch cleanup
# ---------------------------------------------------------------------------
cp -r A549_dRNA_pre-aligned_All_TED_transcriptome* "$OUT_DIR/" 2>/dev/null || true
echo "### Outputs copied to $OUT_DIR"

exit $rc
