library(splatter)
library(here)

params <- newSplatParams()
sim <- splatSimulate(params)

fn <- here("outputs", "preprocessed", "sce_simulated_Gamma-Poisson_preprocessed.rds")
saveRDS(sim, file = fn)
