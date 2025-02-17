#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=END
#SBATCH --output=run_BOOST-GP_humanDLPFC.o
#SBATCH --error=run_BOOST-GP_humanDLPFC.e
#SBATCH --job-name=BOOST-GP_humanDLPFC

module load conda_R/4.4.x
Rscript run_BOOST-GP_humanDLPFC.R