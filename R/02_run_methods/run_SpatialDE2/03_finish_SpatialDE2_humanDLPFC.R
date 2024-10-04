library(SpatialExperiment)
library(here)

# ---------
# load data
# ---------

fn <- here("outputs", "preprocessed", "spe_humanDLPFC_preprocessed.rds")
spe <- readRDS(file=fn)

fn <- here("outputs", "results", "SpatialDE2_results_humanDLPFC.csv")
results <- read.csv(file=fn, header = T, sep = ",", row.names = 1)

# ---------
# merge results into spe
# ---------

df_merge <- merge(rowData(spe), results, by=0)
row.names(df_merge) <- df_merge$Row.names
df_merge$Row.names <- NULL

# rank largest fraction of spatial variances on top
df_merge$SpatialDE2_rank <- rank(-df_merge$FSV, ties.method = "first")

# add to spe
df_merge <- df_merge[match(rownames(spe), rownames(df_merge)), ]
rowData(spe)$SpatialDE2_rank <- df_merge$SpatialDE2_rank

# ---------
# save spe
# ---------

fn <- here("outputs", "results", "spe_humanDLPFC_SpatialDE2.rds")
saveRDS(spe, file=fn)
