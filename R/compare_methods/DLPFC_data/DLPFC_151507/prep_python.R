library(SpatialExperiment)
library(scran)
library(spatialLIBD)
library(nnSVG)

#already contains logcounts
spe <- fetch_data(type = "spe")

spe <- spe[, spe$sample_id == "151507"]
spe <- spe[, !is.na(unfactor(colData(spe)$spatialLIBD))]

## Remove genes without enough data
spe <- filter_genes(spe, 2, 0.2)

# sparse to dense matrix
dense_counts_mat <- as.matrix(counts(spe))

# save csv with 5 genes
write.csv(dense_counts_mat[c(1:5),],
          '/users/kshah/MeanVarBias/R/compare_methods/DLPFC_data/DLPFC_151507/counts_5.csv',
          row.names = T)

# save csv with all genes after filter
write.csv(dense_counts_mat,
          '/users/kshah/MeanVarBias/R/compare_methods/DLPFC_data/DLPFC_151507/counts.csv',
          row.names = T)

# save coords csv
write.csv(spatialCoords(spe),
          '/users/kshah/MeanVarBias/R/compare_methods/DLPFC_data/DLPFC_151507/coords.csv',
          row.names = T)

