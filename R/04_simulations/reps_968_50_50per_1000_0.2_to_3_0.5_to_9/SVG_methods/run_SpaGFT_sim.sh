#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=END
#SBATCH --output=run_SpaGFT_sim.o
#SBATCH --error=run_SpaGFT_sim.e
#SBATCH --job-name=SpaGFT_sim

# create conda environment
# conda env create -f ../../../envs/SpaGFT/environment.yml


module load conda_R/4.4.x
Rscript 01_prep_sim.R

module load conda
conda activate spagft_env
python 02_run_SpaGFT_sim.py

Rscript 03_finish_SpaGFT_sim.R