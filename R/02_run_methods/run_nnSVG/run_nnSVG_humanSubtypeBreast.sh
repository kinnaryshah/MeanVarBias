#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --output=run_nnSVG_humanSubtypeBreast.o
#SBATCH --error=run_nnSVG_humanSubtypeBreast.e
#SBATCH --job-name=nnSVG_humanSubtypeBreast

module load conda_R/4.4.x
Rscript run_nnSVG_humanSubtypeBreast.R

