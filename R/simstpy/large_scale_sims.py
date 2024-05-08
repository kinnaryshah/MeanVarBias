import numpy as np
import pandas as pd
import scipy as sp
import scanpy as sc
import squidpy as sq
import anndata as ad
import seaborn as sns
import simstpy as sim
import matplotlib.pyplot as plt

n_SVGs = 25*5
n_nonSVGs = 75*5
# simulate SVGs with low sigma
adata1 = sim.rna.sim_svgs(n_svgs=n_SVGs, n_non_svgs=0, alpha=0.1, sigma=2, library_size=100)

adata1.layers['counts'] = adata1.X.copy()
gene_names = ["SVG_2_" + str(i) for i in range(1, n_SVGs+1)]
adata1.var_names = gene_names

mean_expression_across_genes = np.mean(adata1.X, axis=0)
min_value = np.min(mean_expression_across_genes)
max_value = np.max(mean_expression_across_genes)

print("Minimum:", min_value)
print("Maximum:", max_value)
# how many genes have means less than 10
print("Number of genes with mean less than 10:", np.sum(mean_expression_across_genes < 10))

# simulate SVGs with high sigma
adata2 = sim.rna.sim_svgs(n_svgs=n_SVGs, n_non_svgs=0, alpha=0.1, sigma=3, library_size=100)

adata2.layers['counts'] = adata2.X.copy()
gene_names = ["SVG_3_" + str(i) for i in range(1, n_SVGs+1)]
adata2.var_names = gene_names

mean_expression_across_genes = np.mean(adata2.X, axis=0)
min_value = np.min(mean_expression_across_genes)
max_value = np.max(mean_expression_across_genes)

print("Minimum:", min_value)
print("Maximum:", max_value)

# how many genes have means less than 10
print("Number of genes with mean less than 10:", np.sum(mean_expression_across_genes < 10))

init_adata = ad.concat([adata1, adata2], axis=1)
init_adata.write_h5ad("large_scale_sims_init.h5ad")

# simulate nonSVGs with low sigma
adata3 = sim.rna.sim_svgs(n_svgs=0, n_non_svgs=n_nonSVGs, alpha=0.1, sigma=2, library_size=100)

adata3.layers['counts'] = adata3.X.copy()
gene_names = ["nonSVG_2_" + str(i) for i in range(1, n_nonSVGs+1)]
adata3.var_names = gene_names

fold_change_factors = np.arange(1, 11)
fold_change_factors = np.tile(fold_change_factors, round(n_nonSVGs/5))[:n_nonSVGs]
inflated_counts_matrix = adata3.X.multiply(fold_change_factors)
adata3.X = inflated_counts_matrix

mean_expression_across_genes = np.mean(adata3.X, axis=0)
min_value = np.min(mean_expression_across_genes)
max_value = np.max(mean_expression_across_genes)

print("Minimum:", min_value)
print("Maximum:", max_value)
# how many genes have means less than 10
print("Number of genes with mean less than 10:", np.sum(mean_expression_across_genes < 10))

# simulate nonSVGs with high sigma
adata4 = sim.rna.sim_svgs(n_svgs=0, n_non_svgs=n_nonSVGs, alpha=0.1, sigma=3, library_size=100)

adata4.layers['counts'] = adata4.X.copy()
gene_names = ["nonSVG_3_" + str(i) for i in range(1, n_nonSVGs+1)]
adata4.var_names = gene_names

fold_change_factors = np.arange(1, 11)
fold_change_factors = np.tile(fold_change_factors, round(n_nonSVGs/5))[:n_nonSVGs]
inflated_counts_matrix = adata4.X.multiply(fold_change_factors)
adata4.X = inflated_counts_matrix

mean_expression_across_genes = np.mean(adata4.X, axis=0)
min_value = np.min(mean_expression_across_genes)
max_value = np.max(mean_expression_across_genes)

print("Minimum:", min_value)
print("Maximum:", max_value)
# how many genes have means less than 10
print("Number of genes with mean less than 10:", np.sum(mean_expression_across_genes < 10))

combined_adata = ad.concat([adata1, adata2], axis=1)
combined_adata = ad.concat([combined_adata, adata3], axis=1)
combined_adata = ad.concat([combined_adata, adata4], axis=1)

# Extract spatial coordinates from adata3.obsm
spatial_coordinates = adata3.obsm["spatial"]

# Add spatial coordinates to observation metadata
combined_adata.obs["X_coord"] = spatial_coordinates[:, 0]
combined_adata.obs["Y_coord"] = spatial_coordinates[:, 1]

combined_adata.write_h5ad("large_scale_sims.h5ad")

