#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=END
#SBATCH --output=run_SMASH_humanDLPFC.o
#SBATCH --error=run_SMASH_humanDLPFC.e
#SBATCH --job-name=SMASH_humanDLPFC

# conda create -n SMASH_env
# conda activate SMASH_env
# downladed .zip file from https://github.com/sealx017/SMASH-package and unzipped in /envs folder
# pip install matplotlib
# pip install matplotlib-venn
# pip install blosc

module load conda
conda activate SMASH_env
python run_SMASH_humanDLPFC.py