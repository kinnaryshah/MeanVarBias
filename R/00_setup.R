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

# create `/plots/supplementary` folder
if(!file.exists(here("plots", "supplementary"))){
  dir.create(here("plots", "supplementary"))
}

# create `/outputs/logs` folder
if(!file.exists(here("outputs", "logs"))){
  dir.create(here("outputs", "logs"))
}

###############################
##### 3rd-level folders #######
###############################

# create `/outputs/logs/02_run_methods` folder
if(!file.exists(here("outputs", "logs", "02_run_methods"))){
  dir.create(here("outputs", "logs", "02_run_methods"))
}