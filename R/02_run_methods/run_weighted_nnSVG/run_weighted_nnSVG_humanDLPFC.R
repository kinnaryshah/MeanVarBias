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

fn <- here("outputs", "preprocessed", "spe_humanDLPFC_preprocessed.rds")
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

# -----------
# save objects
# -----------

file = here("outputs","results","spe_humanDLPFC_weighted_nnSVG.RData")
save(spe, weights, file = file)

# -----------
# session information
# -----------

sessionInfo()


