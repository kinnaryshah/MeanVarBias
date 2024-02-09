library(SpatialExperiment)
library(nnSVG)
library(scran)
library(spatialLIBD)

#already contains logcounts
spe <- fetch_data(type = "spe")

spe <- spe[, spe$sample_id == "151507"]
spe <- spe[, !is.na(unfactor(colData(spe)$spatialLIBD))]
dim(spe)

## Remove genes without enough data
spe <- filter_genes(spe, 2, 0.2)

dim(spe)

set.seed(123)
runtime <- system.time({
  tryCatch({
    spe <- nnSVG(
      spe,
      X = NULL,
      assay_name = "logcounts",
      n_neighbors = 10,
      order = "AMMD",
      n_threads = 10,
      verbose = FALSE
    )
  }, error = function(e) {
    print(paste("Error:", e$message))
  })
})

file = "spe_nnSVG.rds"
saveRDS(spe, file = file)
