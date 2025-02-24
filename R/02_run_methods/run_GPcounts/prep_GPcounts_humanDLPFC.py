import pandas as pd 
import numpy as np 
import gpflow
import tensorflow as tf
from GPcounts import NegativeBinomialLikelihood
from GPcounts.RNA_seq_GP import rna_seq_gp
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
import statsmodels.formula.api as smf 
import statsmodels.api as sm

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
# scale data
# ---------

# need to transpose counts matrix for AnnData object
count_mtx = pd.read_csv(counts_dir,index_col=0).transpose()
coord_mtx = pd.read_csv(coords_dir,index_col=0)
coord_mtx = coord_mtx.rename(columns={'pxl_col_in_fullres': 'x', 'pxl_row_in_fullres': 'y'})

coord_mtx['total_counts'] = count_mtx.sum(1)

Y = count_mtx
spatial_locations = coord_mtx

scales = []
for i in range(0,len(Y.columns)):
    model=smf.glm(formula = "Y.iloc[:,i]~0+spatial_locations['total_counts']", data=Y,family=sm.families.NegativeBinomial(sm.families.links.identity())).fit()
    res = model.params[0]*spatial_locations['total_counts']
    scales.append(res)
scalesdf=pd.DataFrame(scales)
scalesdf=scalesdf.T

# save df as csv to open in R
df_dir = base_dir / "outputs" / "results" / "GPcounts_scales_humanDLPFC.csv"
scalesdf.to_csv(df_dir)
