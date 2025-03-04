#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=END
#SBATCH --output=run_SpaGFT_humanDLPFC.o
#SBATCH --error=run_SpaGFT_humanDLPFC.e
#SBATCH --job-name=SpaGFT_humanDLPFC

# create conda environment
# conda env create -f ../../../envs/SpaGFT/environment.yml


module load conda_R/4.4.x
Rscript 01_prep_humanDLPFC.R

module load conda
conda activate spagft_env
python 02_run_SpaGFT_humanDLPFC.py

Rscript 03_finish_SpaGFT_humanDLPFC.R