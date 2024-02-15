# installed in virtual environment spatialde2_env
# pip install git+https://github.com/PMBio/SpatialDE

import pandas as pd
import numpy as np
import scanpy as sc
import SpatialDE as sd
from tqdm import tqdm
from functools import partialmethod

# need to transpose counts matrix for AnnData object
count_mtx = pd.read_csv('/users/kshah/MeanVarBias/R/compare_methods/DLPFC_data/DLPFC_151507/counts.csv',index_col=0).transpose()
coord_mtx = pd.read_csv('/users/kshah/MeanVarBias/R/compare_methods/DLPFC_data/DLPFC_151507/coords.csv',index_col=0)
coord_mtx = coord_mtx.rename(columns={'pxl_col_in_fullres': 'array_col', 'pxl_row_in_fullres': 'array_row'})
# fit() looks for coordinates stored in Key in ``adata.obsm`
adata = sc.AnnData(count_mtx, obsm={"spatial": coord_mtx})
adata.raw = adata

sc.pp.normalize_total(adata, target_sum=1e4)
df = sd.fit(adata, normalized=True, control=None)

df.set_index("gene", inplace=True)
df.sort_values("FSV", inplace=True, ascending=False)

# save df as csv to open in R
df.to_csv('/users/kshah/MeanVarBias/R/compare_methods/DLPFC_data/DLPFC_151507/SpatialDE2_results.csv')
