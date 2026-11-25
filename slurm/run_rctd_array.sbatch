#!/bin/bash
# Template: one array task per sample. Set --array to 1-<number of samples>, edit the account and environment lines, then
#   sbatch slurm/run_rctd_array.sbatch my_rctd_config.yaml my_samples.tsv
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --job-name=rctd
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-4
#SBATCH --output=logs/rctd_%A_%a.out
#SBATCH --error=logs/rctd_%A_%a.err
set -euo pipefail
# module load R/4.5.2            # or: conda activate <env>
CONFIG="${1:?config yaml}"; SAMPLES="${2:?sample sheet}"
ID=$(awk -F'\t' -v n="$SLURM_ARRAY_TASK_ID" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "sample_id") c = i; next } NR == n + 1 { print $c }' "$SAMPLES")
Rscript run_rctd.R --config "$CONFIG" --samples "$SAMPLES" --sample "$ID"
