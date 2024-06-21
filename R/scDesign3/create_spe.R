library(SpatialExperiment)

for(num_de in c(50, 100, 150, 200)){
  
  print(num_de)

  loc = read.csv(file = "location.csv")
  example_newcount = read.csv(file = glue::glue('counts_svg_{num_de}.csv'))
  de_idx = read.csv(file = glue::glue('.labels_svg_{num_de}.csv'))
  non_de_idx <- read.csv(file = glue::glue('labels_non_svg_{num_de}.csv'))
  
  rownames(example_newcount) = example_newcount$X
  example_newcount$X = NULL
  
  rownames(loc) = loc$X
  loc$X = NULL
  loc$x = loc$spatial1
  loc$y = loc$spatial2
  loc$spatial1 = NULL
  loc$spatial2 = NULL
  
  # create SpatialExperiment object
  spe <- SpatialExperiment(
    assays = list(counts=example_newcount),
    spatialCoords = as.matrix(loc))
  
  de_indices = which(rownames(spe) %in% de_idx$x)
  # make true for de_indices and false for other indices
  rowData(spe)$de = rep(FALSE, nrow(spe))
  rowData(spe)$de[de_indices] = TRUE
  
  non_de_indices = which(rownames(spe) %in% non_de_idx$x)
  # make true for non_de_indices and false for other indices
  rowData(spe)$non_de = rep(FALSE, nrow(spe))
  rowData(spe)$non_de[non_de_indices] = TRUE
  
  saveRDS(spe, file = glue::glue('spe_{num_de}.rds'))
  
  }

