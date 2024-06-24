
# Create output folders for the preprocessed data from 01_preprocessing
if(!file.exists(here("outputs"))){
  dir.create(here("outputs"))
}

if(!file.exists(here("outputs", "preprocessed"))){
  dir.create(here("outputs", "preprocessed"))
}
