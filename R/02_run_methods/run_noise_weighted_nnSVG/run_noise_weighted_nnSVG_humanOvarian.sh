#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --output=run_noise_weighted_nnSVG_humanOvarian.o
#SBATCH --error=run_noise_weighted_nnSVG_humanOvarian.e
#SBATCH --job-name=noise_weighted_nnSVG_humanOvarian

module load conda_R/4.4.x
Rscript run_noise_weighted_nnSVG_humanOvarian.R

