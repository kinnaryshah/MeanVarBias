
# devtools::install_github('xzhoulab/SPARK')
library(SpatialExperiment)
library(SPARK)
library(here)
library(scran)

# ---------
# load data
# ---------

sim_dir <- "reps_968_50_50per_1000_0.2_to_3_0.5_to_9"
K=1
fn = here("outputs", "simulations", sim_dir, paste0("spe_simulation_",K,".rds"))
spe <- readRDS(file=fn)
spe <- logNormCounts(spe)

# ---------
# run SPARKX
# ---------

set.seed(123)
sparkx_out <- sparkx(
  count_in = counts(spe),
  locus_in = spatialCoords(spe),
  X_in = NULL,
  numCores = 1,
  option = "mixture",
  verbose = TRUE
)


# store results in spe object
stopifnot(all(rownames(sparkx_out$res_mtest) == rowData(spe)$gene_id))

rowData(spe) <- cbind(rowData(spe), sparkx_out$res_mtest)

# calculate ranks
rowData(spe)$SPARKX_rank <- rank(rowData(spe)$combinedPval, ties.method = "first")

# -----------
# save object
# -----------

file = here("outputs","results","simulations",sim_dir,"spe_SPARKX.rds")
saveRDS(spe, file = file)

# -----------
# session information
# -----------

sessionInfo()