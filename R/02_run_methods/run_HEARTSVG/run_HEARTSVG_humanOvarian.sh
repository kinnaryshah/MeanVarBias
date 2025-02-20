#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=70G
#SBATCH --mail-type=END
#SBATCH --output=run_HEARTSVG_humanOvarian.o
#SBATCH --error=run_HEARTSVG_humanOvarian.e
#SBATCH --job-name=HEARTSVG_humanOvarian

module load conda_R/4.4.x
Rscript run_HEARTSVG_humanOvarian.R
