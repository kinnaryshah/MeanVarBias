library(here)
library(SpatialExperiment)
library(nnSVG)
library(scran)

#load in spe object
load(file = "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/02_build_spe/spe_transform.Rdata")

#subset to one sample
spe <- spe[, spe$sample_id == "V12D07-335_D1"]

spe <- filter_genes(spe)

spe <- logNormCounts(spe)

fn <- here("outputs", "preprocessed", "spe_humanHPC_V12D07-335_D1_preprocessed.rds")
saveRDS(spe, file = fn)
