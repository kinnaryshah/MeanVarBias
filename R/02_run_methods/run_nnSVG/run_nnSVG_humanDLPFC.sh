#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END

module load conda_R/4.4.x
Rscript run_nnSVG_humanDLPFC.R

