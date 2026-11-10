#!/bin/bash
#SBATCH -J nf-driver-a549_full_ted_0606
#SBATCH -o /private/home/hdheath/projects/flair4_ted/logs/nf-driver-%j.log
#SBATCH -c 2
#SBATCH --mem 8G
#SBATCH -t 3-00:00:00
#SBATCH --partition=long
#SBATCH --no-requeue
#
# Nextflow DRIVER as a SLURM job.
# Why: the driver must NOT run in the foreground on the login node (mustard).
# A login-node session reap (systemd-logind sends SIGTERM on logout/disconnect)
# kills the foreground driver, which then scancels all its child tasks -- which
# is exactly what aborted run `fabulous_wilson` (jobs 34089841/44/47) at
# 2026-06-07 00:48 with failedCount=0 (no task crashed; the driver was killed).
# Running the driver itself under SLURM puts it on a compute node, immune to that.
# The driver submits the per-sample FlairTranscriptome tasks via the slurm executor.

set -euo pipefail

source /private/home/hdheath/miniforge3/etc/profile.d/conda.sh
conda activate nextflow_env

cd /private/home/hdheath/projects/flair4_ted

nextflow run /private/home/hdheath/projects/flair4_ted/main.nf \
    -entry ASSEMBLY -profile standard -resume -ansi-log false \
    -w /private/groups/brookslab/hdheath/results/flair4/work \
    --input /private/home/hdheath/projects/flair4_ted/data/a549_full_2026_06_05/a549_full.csv \
    --params_file /private/home/hdheath/projects/flair4_ted/data/a549_full_2026_06_05/a549_full_ted.json \
    --test_name a549_full_ted_0606
