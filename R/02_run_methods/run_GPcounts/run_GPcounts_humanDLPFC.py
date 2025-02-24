import pandas as pd 
import numpy as np 
import gpflow
import tensorflow as tf
from GPcounts.RNA_seq_GP import rna_seq_gp
from GPcounts.GPcounts_Module import Fit_GPcounts
from matplotlib import pyplot as plt
from pathlib import Path

# Set up the base directory (typically the project root)
base_dir = Path(__file__).resolve().parent
base_dir = base_dir / ".." / ".." / ".." 
print(base_dir)

# ---------
# load data
# ---------

data_dir = base_dir / "outputs" / "preprocessed"
counts_dir = data_dir / "counts_humanDLPFC_preprocessed.csv"
coords_dir = data_dir / "coords_humanDLPFC_preprocessed.csv"
print("File path:", counts_dir)

# ---------
# run GPcounts
# ---------

# need to transpose counts matrix for AnnData object
count_mtx = pd.read_csv(counts_dir,index_col=0).transpose()
coord_mtx = pd.read_csv(coords_dir,index_col=0)
coord_mtx = coord_mtx.rename(columns={'pxl_col_in_fullres': 'x', 'pxl_row_in_fullres': 'y'})

count_mtx = count_mtx.T

gene_name = []
gene_name = count_mtx.index
likelihood = 'Negative_binomial' 
gp_counts = rna_seq_gp(coord_mtx,count_mtx.loc[gene_name], safe_mode=False)

log_likelihood_ratio = gp_counts.One_sample_test(likelihood)
log_likelihood_ratio

results = gp_counts.calculate_FDR(log_likelihood_ratio)
results

# save df as csv to open in R
df_dir = base_dir / "outputs" / "results" / "GPcounts_results_humanDLPFC.csv"
results.to_csv(df_dir)
