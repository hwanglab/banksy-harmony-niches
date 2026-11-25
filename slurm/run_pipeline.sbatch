#!/bin/bash
# Template: edit the account, the environment lines and the config path, then `sbatch slurm/run_pipeline.sbatch`.
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --job-name=banksy_harmony
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=logs/banksy_harmony_%j.out
#SBATCH --error=logs/banksy_harmony_%j.err
set -euo pipefail
# module load R/4.5.2            # or: conda activate <env with the packages listed in the README>
Rscript run_pipeline.R --config "${1:-config/example_config.yaml}"
