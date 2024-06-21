library(here)
library(SpatialExperiment)
library(nnSVG)
library(scran)

#load in spe object
load(file = "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/02_build_spe/spe_transform.Rdata")

#subset to one sample
spe <- spe[, spe$sample_id == "V10B01-085_A1"]

spe <- filter_genes(spe)

spe <- logNormCounts(spe)

fn <- here("outputs", "preprocessed", "spe_humanHPC_preprocessed.rds")
saveRDS(spe, file = fn)



