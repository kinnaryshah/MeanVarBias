#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --mail-user=kinnaryshahh@gmail.com

#module load conda
#conda activate simstpy_env
#python large_scale_sims.py

module load conda_R
Rscript large_scale_sims.R

