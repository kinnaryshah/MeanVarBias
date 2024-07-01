#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --output=intermediate_weights_humanBreast.o
#SBATCH --error=intermediate_weights_humanBreast.e
#SBATCH --job-name=intermediate_weights_humanBreast

module load conda_R/4.4.x
Rscript intermediate_weights_humanBreast.R

