library(SpatialExperiment)
library(nnSVG)
library(scran)
library(here)
library(limma)
library(Matrix)
library(scuttle)
library(BiocParallel)

# ---------
# load data
# ---------

fn <- here("outputs", "preprocessed", "spe_humanDLPFC_preprocessed.rds")
spe <- readRDS(file=fn)

# ---------
# run weighted nnSVG
# ---------

weights <- voom(
  counts(spe))$weights

spe <- weighted_nnSVG(
  input = spe,
  w = weights,
  BPPARAM = MulticoreParam(workers = 10, RNGseed = 5))

print(spe)

# -----------
# save objects
# -----------

file = here("outputs","results","spe_humanDLPFC_weighted_nnSVG.rds")
saveRDS(spe, file = file)

file = here("outputs","results","spe_humanDLPFC_weights.rds")
saveRDS(weights, file = file)

# -----------
# session information
# -----------

sessionInfo()


