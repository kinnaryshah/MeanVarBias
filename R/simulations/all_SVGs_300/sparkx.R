library(SPARK)
library(SpatialExperiment)
spe <- readRDS(file = "spe_simulation.rds")

sp_count <- counts(spe)
location <- spatialCoords(spe)
sparkX_unweighted <- sparkx(sp_count,location,numCores=1,option="mixture")

load(file = "weights.Rdata")
weighted_sp_count <- t(weights)*sp_count
sparkX_weighted <- sparkx(weighted_sp_count,location,X_in=weights,numCores=2,option="mixture")

save(sparkX_unweighted, sparkX_weighted, "spe_simulation_weighted_SPARKX.Rdata")