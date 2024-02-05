library(SpatialExperiment)
library(nnSVG)
library(scran)

fn <- ("../../../simulations/sample_means_300/spe_simulation.rds") 

spe <- readRDS(fn)

dim(spe)

# note that this needs logcounts
spe <- logNormCounts(spe)

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
