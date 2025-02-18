#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --output=run_voom_nnSVG_sim.o
#SBATCH --error=run_voom_nnSVG_sim.e
#SBATCH --job-name=voom_nnSVG_sim

module load conda_R/4.4.x
Rscript run_voom_nnSVG_sim.R