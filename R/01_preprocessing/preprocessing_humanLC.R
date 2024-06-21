library(SpatialExperiment)
library(ggspavis)
library(WeberDivechaLCdata)
library(here)

spe <- WeberDivechaLCdata_Visium()
#subset one sample
spe <- spe[, spe$sample_id == "Br2701_LC_round2"]

#https://github.com/lmweber/locus-c/blob/main/code/analyses/Visium/06_LC_gene_filtering.R
# filter out genes with extremely low expression
# using simple threshold on total UMI counts summed across all spots
n_umis <- 80
ix_low_genes <- rowSums(counts(spe)) < n_umis
table(ix_low_genes)

spe <- spe[!ix_low_genes, ]
dim(spe)

fn <- here("outputs", "preprocessed", "spe_humanLC_preprocessed.rds")
saveRDS(spe, file = fn)
