# run in envs folder
# git clone https://github.com/Minzhe/BOOST-GP
library(SpatialExperiment)
library(here)
setwd(here("envs","BOOST-GP"))
source("R/boost.gp.R")

# ---------
# load data
# ---------

fn <- here("outputs", "preprocessed", "spe_humanDLPFC_preprocessed.rds")
spe <- readRDS(file=fn)

# ---------
# run BOOST-GP
# ---------

set.seed(123)
boostgp_out <- boost.gp(Y = counts(spe), loc = spatialCoords(spe))

saveRDS(boostgp_out, file = here("outputs","results","intermediate_humanDLPFC_BOOST-GP.rds"))

# store results in spe object
stopifnot(all(rownames(boostgp_out) == rowData(spe)$gene_id))

rowData(spe) <- cbind(rowData(spe), boostgp_out)

# calculate ranks
rowData(spe)$BOOSTGP_rank <- rank(rowData(spe)$pval, ties.method = "first")

# -----------
# save object
# -----------

file = here("outputs","results","spe_humanDLPFC_BOOST-GP.rds")
saveRDS(spe, file = file)

# -----------
# session information
# -----------

sessionInfo()