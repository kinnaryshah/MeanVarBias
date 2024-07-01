generate_weights <- function(input, spatial_coords = NULL,
                             assay_name = "counts",
                             stabilize = TRUE,
                             n_threads = 1, BPPARAM = NULL){
  
  if (is(input, "SpatialExperiment")) {
    spe <- input
    stopifnot(assay_name %in% assayNames(spe))
  }
  
  if (is.null(BPPARAM)) {
    BPPARAM <- MulticoreParam(workers = n_threads)
  }
  
  #calculate weights
  
  # Count Matrix, transpose so each row is a spot, and each column is a gene
  if (is(input, "SpatialExperiment")) {
    r <- t(as.matrix(assays(spe)[[assay_name]]))
    coords <- spatialCoords(spe)
    
  } else {
    stopifnot(dim(spatial_coords)[1] == dim(input)[2])
    r <- t(as.matrix(input))
    coords <- spatial_coords
  }
  
  n <- dim(r)[1] # Number of Cells
  G <- dim(r)[2] # Number of Genes
  
  # Sample-specific Library Size
  R <- rowSums(r)
  stopifnot(length(R)==n)
  
  # Temporary Matrix, replicate library size for each row
  tmp_R_mat <- matrix(
    rep(R, each = G),
    byrow = TRUE, nrow = n, ncol = G
  )
  
  # logCPM
  y <- log2(r+0.5) - log2(tmp_R_mat+1) + log2(10^6)
  
  # scale coordinates proportionally
  range_all <- max(apply(coords, 2, function(col) diff(range(col))))
  coords <- apply(coords, 2, function(col) (col - min(col)) / range_all)
  
  # calculate ordering of coordinates
  order_brisc <- BRISC_order(coords, order = "AMMD", verbose = FALSE)
  
  # calculate nearest neighbors
  nn_brisc <- BRISC_neighbor(coords, n.neighbors = 10, n_omp = 1,
                             search.type = "tree", ordering = order_brisc,
                             verbose = FALSE)
  
  # run BRISC using parallelization
  # run BRISC by column of y so BRISC is run per gene
  ix <- seq_len(ncol(y))
  out_brisc <- bplapply(ix, function(i) {
    # fit model (intercept-only model if x is NULL)
    y_i <- y[ ,i]
    suppressWarnings({
      runtime <- system.time({
        out_i <- BRISC_estimation(coords = coords, y = y_i, x = NULL,
                                  cov.model = "exponential",
                                  ordering = order_brisc, neighbor = nn_brisc,
                                  verbose = FALSE)
      })
    })
    
    
    pred_i <- BRISC_prediction(out_i, coords.0 = coords, X.0 = NULL,
                               verbose = FALSE)
    residual_i <- y_i - pred_i$prediction
    
    return(c(pred_i$prediction, residual_i))
  }, BPPARAM = BPPARAM)
  
  # collapse output list into matrix
  mat_brisc <- do.call("rbind", out_brisc)
  
  # *Voom Variance Modelling -------------------------------------------------
  
  mu_hat <- t(mat_brisc[seq_len(G),seq_len(n), drop=FALSE])
  stopifnot(dim(mu_hat) == c(n, G))
  
  s_g <- as.data.frame(t(mat_brisc[seq_len(G),n+seq_len(n), drop=FALSE])) |>
    apply(MARGIN = 2,  # Column wise
          FUN = sd)
  stopifnot(length(s_g) == G)
  
  y_bar <- colMeans(mu_hat)
  stopifnot(length(y_bar) == G)
  
  # Geometric Mean
  R_tilda <- exp(mean(log(R)))
  # The reason of calculating log is to avoid integer overflow
  
  # Log2 Counts
  # Note: slight notation abuse. Prev r denotes read counts
  r_tilda <- y_bar + log2(R_tilda) - log2(10^6)
  stopifnot(length(r_tilda)==G)
  
  return(list(s_g = s_g, r_tilda = r_tilda))
}

library(SpatialExperiment)
library(nnSVG)
library(scran)
library(here)
library(spoon)
library(Matrix)
library(scuttle)
library(BRISC)
library(BiocParallel)

# ---------
# load data
# ---------

fn <- here("outputs", "preprocessed", "spe_humanBreast_preprocessed.rds")
spe <- readRDS(file=fn)

# ---------
# calculate sg and r_tilda
# ---------

output_list <- generate_weights(
  input = spe,
  stabilize = TRUE,
  BPPARAM = MulticoreParam(workers = 10, RNGseed = 4))

# ---------
# save output
# ---------

fn <- here("outputs", "results", "weights_output_list_humanBreast.rds")
saveRDS(output_list, file=fn)

