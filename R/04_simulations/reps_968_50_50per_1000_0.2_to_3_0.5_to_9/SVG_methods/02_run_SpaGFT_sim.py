import SpaGFT as spg
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os
from pathlib import Path

# Set up the base directory (typically the project root)
base_dir = Path(__file__).resolve().parent
base_dir = base_dir / ".." / ".." / ".." 
print(base_dir)

# ---------
# load data
# ---------

data_dir = base_dir / "outputs" / "preprocessed" / "simulations" / "reps_968_50_50per_1000_0.2_to_3_0.5_to_9"
counts_dir = data_dir / "counts_preprocessed.csv"
coords_dir = data_dir / "coords_preprocessed.csv"
print("File path:", counts_dir)

# ---------
# run SpaGFT
# ---------

# need to transpose counts matrix for AnnData object
count_mtx = pd.read_csv(counts_dir,index_col=0).transpose()
coord_mtx = pd.read_csv(coords_dir,index_col=0)
coord_mtx = coord_mtx.rename(columns={'pxl_col_in_fullres': 'array_col', 'pxl_row_in_fullres': 'array_row'})
adata = sc.AnnData(count_mtx)
adata.obs.loc[:, ['array_col', 'array_row']] = coord_mtx
adata.var_names_make_unique()
adata.raw = adata

# QC
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

# determine the number of low-frequency FMs and high-frequency FMs
(ratio_low, ratio_high) = spg.gft.determine_frequency_ratio(adata,
                                                            ratio_neighbors=1)

gene_df = spg.detect_svg(adata,
                         spatial_info=['array_row', 'array_col'],
                         ratio_low_freq=ratio_low,
                         ratio_high_freq=ratio_high,
                         ratio_neighbors=1,
                         filter_peaks=True,
                         S=6)
# S determines the  sensitivity of kneedle algorithm

# save gene_df as csv to open in R
df_dir = base_dir / "outputs" / "results" / "simulations" / "reps_968_50_50per_1000_0.2_to_3_0.5_to_9" / "SpaGFT_results.csv"
gene_df.to_csv(df_dir)
