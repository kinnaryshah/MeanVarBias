#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --mail-user=kinnaryshahh@gmail.com
#SBATCH --time=4-00:00:00

module load conda_R
Rscript run_fdr.R
