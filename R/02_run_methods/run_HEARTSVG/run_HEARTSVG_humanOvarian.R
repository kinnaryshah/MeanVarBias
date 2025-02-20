library(spatstat);library(Seurat);library(plot3D);library(plotly);library(TTR);
library(stringr);library(fpp2);library(seasonal);library(cluster);library(poolr);
library(reshape2);library(data.table);library(plyr);library(dplyr);library(ggplot2);
library(MASS);library(EnvStats);library(BreakPoints);library(gprofiler2)

library(HEARTSVG)
library(SpatialExperiment)
library(here)

# ---------
# load data
# ---------

fn <- here("outputs", "preprocessed", "spe_humanOvarian_preprocessed.rds")
spe <- readRDS(file=fn)

# ---------
# run HEARTSVG
# ---------

df <- spatialCoords(spe)
colnames(df) <- c("col", "row")
counts_mat <- t(counts(spe))
counts_mat <- as(counts_mat, "dgCMatrix")
df <- cbind(df, counts_mat)

results <- heartsvg(df,scale=T)
head(results)

fn <- here("outputs", "results", "intermediate_humanOvarian_HEARTSVG.rds")
saveRDS(results, file=fn)

# ---------
# merge results into spe
# ---------

df_merge <- merge(rowData(spe), results, by.x=0, by.y=1)
row.names(df_merge) <- df_merge$Row.names
df_merge$Row.names <- NULL

# add to spe
df_merge <- df_merge[match(rownames(spe), rownames(df_merge)), ]
rowData(spe) <- df_merge
rowData(spe)$HEARTSVG_rank <- rowData(spe)$rank

# ---------
# save spe
# ---------

fn <- here("outputs", "results", "spe_humanOvarian_HEARTSVG.rds")
saveRDS(spe, file=fn)
