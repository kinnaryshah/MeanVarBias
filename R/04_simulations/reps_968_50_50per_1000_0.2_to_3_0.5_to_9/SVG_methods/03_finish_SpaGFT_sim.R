library(SpatialExperiment)
library(here)

# ---------
# load data
# ---------

sim_dir <- "reps_968_50_50per_1000_0.2_to_3_0.5_to_9"
K=1
fn = here("outputs", "simulations", sim_dir, paste0("spe_simulation_",K,".rds"))
spe <- readRDS(file=fn)

fn <- here("outputs", "results", "simulations", sim_dir, "SpaGFT_results.csv")
results <- read.csv(file=fn, header = T, sep = ",", row.names = 1)

# ---------
# merge results into spe
# ---------

df_merge <- merge(rowData(spe), results, by=0)
row.names(df_merge) <- df_merge$Row.names
df_merge$Row.names <- NULL

# rank by qvalues
df_merge$SpaGFT_rank <- rank(df_merge$qvalue, ties.method = "first")

# add to spe
rowData(spe) <- df_merge

# ---------
# save spe
# ---------

fn = here("outputs","results","simulations",sim_dir,"spe_SpaGFT.rds")
saveRDS(spe, file=fn)