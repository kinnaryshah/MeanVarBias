#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=END
#SBATCH --output=run_SPARKX_humanDLPFC.o
#SBATCH --error=run_SPARKX_humanDLPFC.e
#SBATCH --job-name=SPARKX_humanDLPFC

module load conda_R/4.4.x
Rscript run_SPARKX_humanDLPFC.R