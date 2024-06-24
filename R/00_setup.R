
# Create top-level outputs folder
if(!file.exists(here("outputs"))){
  dir.create(here("outputs"))
}

# create `/outputs/raw` folder
if(!file.exists(here("outputs", "raw"))){
  dir.create(here("outputs", "raw"))
}

# create `/outputs/preprocessed` folder
if(!file.exists(here("outputs", "preprocessed"))){
  dir.create(here("outputs", "preprocessed"))
}
