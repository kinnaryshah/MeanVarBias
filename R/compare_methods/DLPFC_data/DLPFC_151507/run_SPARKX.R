library(SpatialExperiment)
library(SPARK)

fn <- ("../../../simulations/sample_means_300/spe_simulation.rds") 

spe <- readRDS(fn)

dim(spe)

set.seed(123)
runtime <- system.time({
  sparkx_out <- sparkx(
    count_in = counts(spe),
    locus_in = spatialCoords(spe),
    X_in = NULL,
    numCores = 1,
    option = "mixture",
    verbose = TRUE
  )
})

# results for combined kernels
head(sparkx_out$res_mtest)

# store results in SPE object
stopifnot(all(rownames(sparkx_out$res_mtest) == rowData(spe)$gene_name))

rowData(spe) <- cbind(rowData(spe), sparkx_out$res_mtest)

# calculate ranks
rowData(spe)$SPARKX_rank <- rank(rowData(spe)$combinedPval, ties.method = "first")

file <- "spe_SPARKX.rds"
saveRDS(spe, file=file)
