# installation following
# https://spagft.readthedocs.io/en/latest/Installation.html

# created csv files of counts mat and coords mat to read in 

import SpaGFT as spg
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os

sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Load data

# need to transpose counts matrix for AnnData object
count_mtx = pd.read_csv('/users/kshah/MeanVarBias/R/compare_methods/DLPFC_data/DLPFC_151507/counts.csv',index_col=0).transpose()
coord_mtx = pd.read_csv('/users/kshah/MeanVarBias/R/compare_methods/DLPFC_data/DLPFC_151507/coords.csv',index_col=0)
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
gene_df.to_csv('/users/kshah/MeanVarBias/R/compare_methods/DLPFC_data/DLPFC_151507/SpaGFT_results.csv')
