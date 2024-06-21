library(splatter)
library(here)

set.seed(100)

# use default parameters from splatter: 10k genes, 100 cells
params <- newSplatParams()
sim <- splatSimulate(params)

fn <- here("outputs", "preprocessed", "sce_simulated_Gamma-Poisson_preprocessed.rds")
saveRDS(sim, file = fn)
