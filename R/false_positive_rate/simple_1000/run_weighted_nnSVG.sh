#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --mail-user=kinnaryshahh@gmail.com

module load conda_R
Rscript run_weighted_nnSVG.R
