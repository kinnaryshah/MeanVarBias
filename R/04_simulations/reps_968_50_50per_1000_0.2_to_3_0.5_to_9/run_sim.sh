#!/bin/bash
#SBATCH --job-name=run_sim
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-50
#SBATCH --mail-type=END
#SBATCH --output=sim.%a.o
#SBATCH --error=sim.%a.e
#SBATCH --time=4-00:00:00

module load conda_R/4.4.x
Rscript run_sim.R
