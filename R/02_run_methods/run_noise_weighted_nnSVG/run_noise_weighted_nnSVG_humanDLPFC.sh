#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --output=run_noise_weighted_nnSVG_humanDLPFC.o
#SBATCH --error=run_noise_weighted_nnSVG_humanDLPFC.e
#SBATCH --job-name=noise_weighted_nnSVG_humanDLPFC

module load conda_R/4.4.x
Rscript run_noise_weighted_nnSVG_humanDLPFC.R

