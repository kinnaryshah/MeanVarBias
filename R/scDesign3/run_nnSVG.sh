#!/bin/bash
#SBATCH --job-name=nnSVG
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=7G
#SBATCH --array=50,100,150,200
#SBATCH --mail-type=END
#SBATCH --output=run.%a.txt
#SBATCH --error=run.%a.txt
#SBATCH --mail-user=kinnaryshahh@gmail.com
#SBATCH --time=2-00:00:00

module load conda_R/4.3.x
Rscript run_nnSVG.R
