#!/bin/bash
#SBATCH --job-name=mt_mod_pipe
#SBATCH --output=log/mt_mod_pipe_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=0-158

set -euo pipefail

mkdir -p log

CONFIG=${CONFIG:-config.example.ini}
PIPE_PY=${PIPE_PY:-mt_pipeline.py}
CONDA_ENV=${CONDA_ENV:-mt_match_env}

module load miniconda/24.11.3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

python3 "$PIPE_PY" run-array --config "$CONFIG"
