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
count_mtx = pd.read_csv('/users/kshah/MeanVarBias/R/simulations/sample_means_300/counts.csv')
coord_mtx = pd.read_csv('/users/kshah/MeanVarBias/R/simulations/sample_means_300/coords.csv')
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
# extract spaitally variable genes
svg_list = gene_df[gene_df.cutoff_gft_score][gene_df.qvalue < 0.05].index.tolist()
print("The number of SVGs: ", len(svg_list))
# the top 20 SVGs
print(svg_list[:20])

plot_svgs = ['CD3E', 'IL7R', 'CCR7', 'PCNA', 'CDK1', 'CDC20', 'CD19', 'CD79B']
sc.pl.spatial(adata, color=plot_svgs, size=1.6, cmap='magma', use_raw=False)

# save gene_df to csv file

if(np.isinf(adata.X).any()):
    print("The Array contain NaN values")
else:
    print("The Array does not contain NaN values")
