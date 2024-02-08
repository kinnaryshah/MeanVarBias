library(scran)
library(spatialDE)
library(SpatialExperiment)
library(reticulate)

use_virtualenv("~/myenv")

fn <- ("../../../simulations/sample_means_300/spe_simulation.rds") 

spe <- readRDS(fn)

dim(spe)

out <- spatialDE(spe, assay_type = "counts", verbose = FALSE)
head(out[order(out$qval), ])

rowData(spe)$spatialDE_rank <- rank(-1 * out$qval, ties.method = "first")

file = "spe_spatialDE.rds"
saveRDS(spe, file=file)

