#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=END
#SBATCH --output=run_SpatialDE2_humanDLPFC.o
#SBATCH --error=run_SpatialDE2_humanDLPFC.e
#SBATCH --job-name=SpatialDE2_humanDLPFC

module load conda_R/4.4.x
Rscript 01_prep_humanDLPFC.R

module load conda
conda activate spatialde2_env
python 02_run_SpatialDE2_humanDLPFC.py