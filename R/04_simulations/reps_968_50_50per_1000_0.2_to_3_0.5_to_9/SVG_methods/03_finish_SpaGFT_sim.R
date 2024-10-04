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

rownames(spe) <- c(1:1000)
df_merge <- merge(rowData(spe), results, by=0)
row.names(df_merge) <- df_merge$Row.names
df_merge$Row.names <- NULL

# re-name
df_merge$SpaGFT_rank <- df_merge$svg_rank

# add to spe
df_merge <- df_merge[match(rownames(spe), rownames(df_merge)), ]
rowData(spe) <- df_merge

# ---------
# save spe
# ---------

fn = here("outputs","results","simulations",sim_dir,"spe_SpaGFT.rds")
saveRDS(spe, file=fn)