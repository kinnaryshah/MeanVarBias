import pandas as pd
import numpy as np
import scanpy as sc
import SpatialDE as sd
from tqdm import tqdm
from functools import partialmethod
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
# run SpatialDE2
# ---------

# need to transpose counts matrix for AnnData object
count_mtx = pd.read_csv(counts_dir,index_col=0).transpose()
coord_mtx = pd.read_csv(coords_dir,index_col=0)
coord_mtx = coord_mtx.rename(columns={'pxl_col_in_fullres': 'array_col', 'pxl_row_in_fullres': 'array_row'})
# fit() looks for coordinates stored in Key in `adata.obsm`
adata = sc.AnnData(count_mtx, obsm={"spatial": coord_mtx})
adata.raw = adata

sc.pp.normalize_total(adata, target_sum=1e4)
df = sd.fit(adata, normalized=True, control=None)

df.set_index("gene", inplace=True)
df.sort_values("FSV", inplace=True, ascending=False)

# save df as csv to open in R
df_dir = base_dir / "outputs" / "results" / "simulations" / "reps_968_50_50per_1000_0.2_to_3_0.5_to_9" / "SpatialDE2_results.csv"
df.to_csv(df_dir)
