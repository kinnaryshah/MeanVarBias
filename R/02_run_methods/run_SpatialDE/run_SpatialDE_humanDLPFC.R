library(SpatialExperiment)
library(spatialDE)
library(here)

# ---------
# load data
# ---------

fn <- here("outputs", "preprocessed", "spe_humanDLPFC_preprocessed.rds")
spe <- readRDS(file=fn)

# ---------
# run SpatialDE
# ---------

set.seed(123)
spatialDE_out <- spatialDE(spe)

# match rownames with g in spatialDE_out
rownames(spatialDE_out) <- spatialDE_out$g
df_merge <- merge(rowData(spe), spatialDE_out, by=0)
row.names(df_merge) <- df_merge$Row.names
df_merge$Row.names <- NULL

# rank largest fraction of spatial variances on top
df_merge$SpatialDE2_rank <- rank(-df_merge$FSV, ties.method = "first")

# add to spe
rowData(spe) <- df_merge

# -----------
# save object
# -----------

file = here("outputs","results","spe_humanDLPFC_SpatialDE.rds")
saveRDS(spe, file = file)

# -----------
# session information
# -----------

sessionInfo()
