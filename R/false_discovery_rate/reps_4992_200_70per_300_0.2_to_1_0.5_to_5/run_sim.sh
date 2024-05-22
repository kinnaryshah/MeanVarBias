#!/bin/bash
#SBATCH --job-name=high_beta
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-5
#SBATCH --mail-type=END
#SBATCH --output=sim.%a.txt
#SBATCH --error=sim.%a.txt
#SBATCH --mail-user=kinnaryshahh@gmail.com
#SBATCH --time=4-00:00:00

module load conda_R
Rscript run_sim.R
