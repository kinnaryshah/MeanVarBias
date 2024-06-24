###################################
# dataset: simulated observations from a negative binomial distribution
###################################
library(splatter)
library(here)

# use default parameters from splatter: 10k genes, 100 cells
set.seed(100)
params <- newSplatParams()
sim <- splatSimulate(params)

# -----------
# save object
# -----------

fn <- here("outputs", "preprocessed", "sce_simulated_Gamma-Poisson_preprocessed.rds")
saveRDS(sim, file = fn)
