#!/bin/bash
#SBATCH --job-name=sparkx
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=END
#SBATCH --mail-user=kinnaryshahh@gmail.com
module load conda_R
Rscript sparkx.R
