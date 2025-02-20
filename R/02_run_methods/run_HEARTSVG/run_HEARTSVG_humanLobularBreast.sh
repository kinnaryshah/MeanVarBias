#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60G
#SBATCH --mail-type=END
#SBATCH --output=run_HEARTSVG_humanLobularBreast.o
#SBATCH --error=run_HEARTSVG_humanLobularBreast.e
#SBATCH --job-name=HEARTSVG_humanLobularBreast

module load conda_R/4.4.x
Rscript run_HEARTSVG_humanLobularBreast.R
