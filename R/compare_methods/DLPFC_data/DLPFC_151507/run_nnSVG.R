library(SpatialExperiment)
library(nnSVG)
library(scran)
library(spatialLIBD)

#already contains logcounts
spe <- fetch_data(type = "spe")

spe <- spe[, spe$sample_id == "151507"]
spe <- spe[, !is.na(unfactor(colData(spe)$spatialLIBD))]
dim(spe)

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) <=5)
length(no_expr)
length(no_expr) / nrow(spe) * 100
spe <- spe[-no_expr, ]
dim(spe)

set.seed(123)
runtime <- system.time({
  spe <- nnSVG(
    spe,
    X = NULL,
    assay_name = "logcounts",
    n_neighbors = 10,
    order = "AMMD",
    n_threads = 10,
    verbose = FALSE
  )
})

file = "spe_nnSVG.rds"
saveRDS(spe, file = file)
