
# conda create -n SMASH_env
# conda activate SMASH_env
# downladed .zip file from https://github.com/sealx017/SMASH-package and unzipped in /envs folder
# pip install matplotlib
# pip install matplotlib-venn
# pip install blosc

# conda activate SMASH_env

import os
path = "/users/kshah/MeanVarBias/envs/SMASH-package-main/"
os.chdir(path + "/SMASH")

import SMASH 
import pickle
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from matplotlib import pyplot as plt
from matplotlib_venn import venn2_unweighted
pd.set_option('display.float_format', '{:.2g}'.format)
from pathlib import Path

# Set up the base directory (typically the project root)
base_dir = path + "../../"
base_dir = base_dir / ".." / ".." / ".." 
print(base_dir)

# ---------
# load data
# ---------

data_dir = base_dir + "outputs/preprocessed/"
counts_dir = data_dir + "counts_humanLobularBreast_preprocessed.csv"
coords_dir = data_dir + "coords_humanLobularBreast_preprocessed.csv"
print("File path:", counts_dir)

# ---------
# run SMASH
# ---------

# need to transpose counts matrix for AnnData object
count_mtx = pd.read_csv(counts_dir,index_col=0).transpose()
coord_mtx = pd.read_csv(coords_dir,index_col=0)
coord_mtx = coord_mtx.rename(columns={'pxl_col_in_fullres': 'array_col', 'pxl_row_in_fullres': 'array_row'})

Run_method = SMASH.SMASH(count_mtx, coord_mtx)

SMASH_pvals = Run_method['SMASH']
alpha_cutoff = 0.01
SMASH_pvals['Adjusted p-val'] = multipletests(SMASH_pvals['p-val'], alpha=alpha_cutoff, method='fdr_by')[1]
SVGS_SMASH = SMASH_pvals[SMASH_pvals['Adjusted p-val'] < alpha_cutoff]
print(SVGS_SMASH)
print({'# detected SVGs ': SVGS_SMASH.shape[0]})

# save df as csv to open in R
df_dir = base_dir + "outputs/results/SMASH_results_humanLobularBreast.csv"
SMASH_pvals.to_csv(df_dir)
