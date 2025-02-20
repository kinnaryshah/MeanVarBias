#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --mail-type=END
#SBATCH --output=run_HEARTSVG_humanDLPFC.o
#SBATCH --error=run_HEARTSVG_humanDLPFC.e
#SBATCH --job-name=HEARTSVG_humanDLPFC

module load conda_R/4.4.x
Rscript run_HEARTSVG_humanDLPFC.R