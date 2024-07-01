library(SpatialExperiment)
library(scran)
library(here)

# ---------
# load data
# ---------

fn <- here("outputs", "preprocessed", "spe_humanDLPFC_preprocessed.rds")
spe <- readRDS(file=fn)

# ---------
# run HVG
# ---------

set.seed(123)
dec <- modelGeneVar(spe)

res <- unlist(res)

# store results in spe object
rowData(spe) <- cbind(rowData(spe), dec)

# calculate ranks
rowData(spe)$HVG_rank <- rank(-1 * rowData(spe)$bio, ties.method = "first")

# -----------
# save object
# -----------

file = here("outputs","results","spe_humanDLPFC_HVG.rds")
saveRDS(spe, file = file)

# -----------
# session information
# -----------

sessionInfo()