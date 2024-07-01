#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --output=run_MoransI_humanDLPFC.o
#SBATCH --error=run_MoransI_humanDLPFC.e
#SBATCH --job-name=MoransI_humanDLPFC

module load conda_R/4.4.x
Rscript run_MoransI_humanDLPFC.R