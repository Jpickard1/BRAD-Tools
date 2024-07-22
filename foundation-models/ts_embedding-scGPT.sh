#!/bin/bash

# SBATCH --array=0-30
# SBATCH --job-name=scGPT
# SBATCH --mail-type=BEGIN,END
# SBATCH --mail-user=jpic@umich.edu
# SBATCH --nodes=1
# SBATCH --mem-per-gpu=180g
# SBATCH --time=10:00:00
# SBATCH --account=indikar0
# SBATCH --partition=gpu,spgpu,gpu_mig40
# SBATCH --gpus=1
# SBATCH --output=/home/jpic/BRAD-Tools/foundation-models/logs

# Load conda
source ~/.bashrc
source $(conda info --base)/etc/profile.d/conda.sh

# Activate environment
conda activate scgpt3

# run the script
python /home/jpic/BRAD-Tools/foundation-models/ts_embedding-scGPT.py $SLURM_ARRAY_TASK_ID

# kill the environment
conda deactivate
