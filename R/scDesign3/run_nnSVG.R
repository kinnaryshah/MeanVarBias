library(SpatialExperiment)
library(MASS)
library(scuttle)

K <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(K)

spe <- readRDS(paste0("spe_",K,".rds"))

### RUN NNSVG

library(nnSVG)
library(SpatialExperiment)
library(BRISC)
library(BiocParallel)
library(scuttle)
library(Matrix)
library(spoon)

print(dim(spe))

spe <- logNormCounts(spe)

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

file = paste0("spe_nnSVG_",K,".rds")
saveRDS(spe, file = file)

spe_unweighted <- spe

### RUN WEIGHTED NNSVG

spe <- readRDS(paste0("spe_",K,".rds"))

spe <- logNormCounts(spe)

weights <- generate_weights(input = spe,
                            stabilize = TRUE,
                            BPPARAM = MulticoreParam(workers = 10,
                                                     RNGseed = 4))

spe <- weighted_nnSVG(input = spe,
                      w = weights,
                      BPPARAM = MulticoreParam(workers = 10, RNGseed = 5))

file = paste0("spe_weighted_nnSVG_",K,".rds")
saveRDS(spe, file = file)

