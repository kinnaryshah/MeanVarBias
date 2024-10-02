#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=END
#SBATCH --output=run_SpatialDE_humanDLPFC.o
#SBATCH --error=run_SpatialDE_humanDLPFC.e
#SBATCH --job-name=SpatialDE_humanDLPFC

module load conda_R/4.4.x
Rscript run_SpatialDE_humanDLPFC.R