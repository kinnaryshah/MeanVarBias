library(here)

###############################
##### Top-level folders #######
###############################

# Create top-level outputs folder
if(!file.exists(here("outputs"))){
  dir.create(here("outputs"))
}

# Create top-level plots folder
if(!file.exists(here("plots"))){
  dir.create(here("plots"))
}

###############################
##### 2nd-level folders #######
###############################

# create `/outputs/raw` folder
if(!file.exists(here("outputs", "raw"))){
  dir.create(here("outputs", "raw"))
}

# create `/outputs/preprocessed` folder
if(!file.exists(here("outputs", "preprocessed"))){
  dir.create(here("outputs", "preprocessed"))
}

# create `/outputs/results` folder
if(!file.exists(here("outputs", "results"))){
  dir.create(here("outputs", "results"))
}

# create `/outputs/simulations` folder
if(!file.exists(here("outputs", "simulations"))){
  dir.create(here("outputs", "simulations"))
}

# create `/plots/supplementary` folder
if(!file.exists(here("plots", "supplementary"))){
  dir.create(here("plots", "supplementary"))
}

# create `/plots/main` folder
if(!file.exists(here("plots", "main"))){
  dir.create(here("plots", "main"))
}

