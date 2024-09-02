library(SpatialExperiment)
library(scran)
library(here)

# ---------
# load data
# ---------

sim_dir <- "reps_968_50_50per_1000_0.2_to_3_0.5_to_9"
K=1
fn = here("outputs", "simulations", sim_dir, paste0("spe_simulation_",K,".rds"))
spe <- readRDS(file=fn)
spe <- logNormCounts(spe)

# ---------
# run HVG
# ---------

set.seed(123)
dec <- modelGeneVar(spe)

# store results in spe object
rowData(spe) <- cbind(rowData(spe), dec)

# calculate ranks
rowData(spe)$HVG_rank <- rank(-1 * rowData(spe)$bio, ties.method = "first")

# -----------
# save object
# -----------

file = here("outputs","results","simulations",sim_dir,"spe_HVG.rds")
saveRDS(spe, file = file)

# -----------
# session information
# -----------

sessionInfo()