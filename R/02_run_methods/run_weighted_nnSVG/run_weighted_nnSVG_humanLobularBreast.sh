#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --time=2-0:00:00
#SBATCH --output=run_weighted_nnSVG_humanLobularBreast.o
#SBATCH --error=run_weighted_nnSVG_humanLobularBreast.e
#SBATCH --job-name=weighted_nnSVG_humanLobularBreast

module load conda_R/4.4.x
Rscript run_weighted_nnSVG_humanLobularBreast.R

