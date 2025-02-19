library(SpatialExperiment)
library(here)

# ---------
# load data
# ---------

fn <- here("outputs", "preprocessed", "spe_humanDLPFC_preprocessed.rds")
spe <- readRDS(file=fn)

fn <- here("outputs", "results", "SMASH_results_humanDLPFC.csv")
results <- read.csv(file=fn, header = T, sep = ",", row.names = 1)

# ---------
# merge results into spe
# ---------

df_merge <- merge(rowData(spe), results, by.x=0, by.y=1)
row.names(df_merge) <- df_merge$Row.names
df_merge$Row.names <- NULL

# add to spe
df_merge <- df_merge[match(rownames(spe), rownames(df_merge)), ]
rowData(spe) <- df_merge
rowData(spe)$SMASH_rank <- rank(rowData(spe)$Adjusted.p.val, ties.method = "first")

# ---------
# save spe
# ---------

fn <- here("outputs", "results", "spe_humanDLPFC_SMASH.rds")
saveRDS(spe, file=fn)
