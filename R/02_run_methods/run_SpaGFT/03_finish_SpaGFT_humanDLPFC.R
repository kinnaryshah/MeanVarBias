library(SpatialExperiment)
library(here)

# ---------
# load data
# ---------

fn <- here("outputs", "preprocessed", "spe_humanDLPFC_preprocessed.rds")
spe <- readRDS(file=fn)

fn <- here("outputs", "results", "SpaGFT_results_humanDLPFC.csv")
results <- read.csv(file=fn, header = T, sep = ",", row.names = 1)

# ---------
# merge results into spe
# ---------

df_merge <- merge(rowData(spe), results, by=0)
row.names(df_merge) <- df_merge$Row.names
df_merge$Row.names <- NULL

# rename rank column
df_merge$SpaGFT_rank <- df_merge$svg_rank

# add to spe
df_merge <- df_merge[match(rownames(spe), rownames(df_merge)), ]
rowData(spe) <- df_merge

# ---------
# save spe
# ---------

fn <- here("outputs", "results", "spe_humanDLPFC_SpaGFT.rds")
saveRDS(spe, file=fn)