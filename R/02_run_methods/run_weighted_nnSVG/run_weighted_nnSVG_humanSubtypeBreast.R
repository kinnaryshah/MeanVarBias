library(SpatialExperiment)
library(nnSVG)
library(scran)
library(here)
library(spoon)
library(Matrix)
library(scuttle)
library(BRISC)
library(BiocParallel)

# ---------
# load data
# ---------

fn <- here("outputs", "preprocessed", "spe_humanSubtypeBreast_preprocessed.rds")
spe <- readRDS(file=fn)

# ---------
# run weighted nnSVG
# ---------

weights <- generate_weights(
  input = spe,
  stabilize = TRUE,
  BPPARAM = MulticoreParam(workers = 10, RNGseed = 4))

spe <- weighted_nnSVG(
  input = spe,
  w = weights,
  BPPARAM = MulticoreParam(workers = 10, RNGseed = 5))

print(spe)

# -----------
# save objects
# -----------

file = here("outputs","results","spe_humanSubtypeBreast_weighted_nnSVG.rds")
saveRDS(spe, file = file)

file = here("outputs","results","spe_humanSubtypeBreast_weights.rds")
saveRDS(weights, file = file)

# -----------
# session information
# -----------

sessionInfo()


