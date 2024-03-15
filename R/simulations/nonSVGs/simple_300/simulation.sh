#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=25G
#SBATCH --mail-type=END
#SBATCH --mail-user=kinnaryshahh@gmail.com
#SBATCH --time=10-00:00:00
module load conda_R
Rscript simulation.R
