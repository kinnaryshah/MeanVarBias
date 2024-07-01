#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --output=run_HVG_humanDLPFC.o
#SBATCH --error=run_HVG_humanDLPFC.e
#SBATCH --job-name=HVG_humanDLPFC

module load conda_R/4.4.x
Rscript run_HVG_humanDLPFC.R