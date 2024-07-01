#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --output=run_nnSVG_humanDLPFC_layers.o
#SBATCH --error=run_nnSVG_humanDLPFC_layers.e
#SBATCH --job-name=nnSVG_humanDLPFC_layers

module load conda_R/4.4.x
Rscript run_nnSVG_humanDLPFC_layers.R

