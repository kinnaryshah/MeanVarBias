#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=END
#SBATCH --output=run_SpatialDE2_sim.o
#SBATCH --error=run_SpatialDE2_sim.e
#SBATCH --job-name=SpatialDE2_sim

# create conda environment
# conda env create -f ../../../envs/SpatialDE2/environment.yml

module load conda_R/4.4.x
Rscript 01_prep_sim.R

module load conda
conda activate spatialde2_env
python 02_run_SpatialDE2_sim.py

Rscript 03_finish_SpatialDE2_sim.R