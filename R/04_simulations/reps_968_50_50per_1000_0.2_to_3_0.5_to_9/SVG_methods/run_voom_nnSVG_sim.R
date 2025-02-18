library(SpatialExperiment)
library(SPARK)
library(here)
library(scran)
library(limma)
library(spoon)
library(BiocParallel)

# ---------
# load data
# ---------

sim_dir <- "reps_968_50_50per_1000_0.2_to_3_0.5_to_9"
K=1
fn = here("outputs", "simulations", sim_dir, paste0("spe_simulation_",K,".rds"))
spe <- readRDS(file=fn)
spe <- logNormCounts(spe)

# ---------
# run voom nnSVG
# ---------

voom_list <- voom(
  counts(spe))

weights <- voom_list$weights
weights <- t(weights)

spe <- weighted_nnSVG(input = spe,
                      w = weights,
                      BPPARAM = MulticoreParam(workers = 10, RNGseed = 5))

file = here("outputs", "simulations", sim_dir, paste0("voom_weights_",K,".rds"))
saveRDS(weights, file = file)

file = here("outputs", "simulations", sim_dir, paste0("spe_voom_weighted_nnSVG_",K,".rds"))
saveRDS(spe, file = file)
