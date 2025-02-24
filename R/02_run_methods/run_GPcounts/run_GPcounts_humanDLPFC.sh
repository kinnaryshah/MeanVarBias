#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --mail-type=END
#SBATCH --output=run_GPcounts_humanDLPFC.o
#SBATCH --error=run_GPcounts_humanDLPFC.e
#SBATCH --job-name=GPcounts_humanDLPFC

# create conda environment
# conda create -n gpcounts
# conda activate gpcounts
# git clone https://github.com/ManchesterBioinference/GPcounts.git
# cd GPcounts
# pip install --user -r requirements.txt
# cd ..
# git clone https://github.com/markvdw/RobustGP
# cd RobustGP
# pip install --user .
# pip install --upgrade tensorflow-probability
# cd ..
# cd GPcounts
# pip install -U tensorflow
# pip install --user .

module load conda
conda activate gpcounts

# already run prep_GPcounts_humanDLPFC.py
#python prep_GPcounts_humanDLPFC.py

python run_GPcounts_humanDLPFC.py

