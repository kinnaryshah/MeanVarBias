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

fn <- here("outputs", "preprocessed", "spe_humanSubtypeBreast_preprocessed.rds")
spe <- readRDS(file=fn)

# new noise weighted function

generate_weights_noise <- function(input, spatial_coords = NULL,
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

  # run BRISC using parallezation
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

  # *PREDICT MODEL -----------------------------------------------------------
  stopifnot(dim(mu_hat)==dim(tmp_R_mat))
  lambda_hat <- mu_hat + log2(tmp_R_mat+1) - log2(10^6)

  spline_fit <- smooth.spline(y=sqrt(s_g), x=r_tilda)

  # NOTE: It is possible that lambda is out of range of r_tilda
  # which will produce NA predicted values due to extrapolation
  tmp_pred_sqrt_sg <- predict(
    spline_fit,
    x = c(lambda_hat)
  )$y |>
    matrix(
      nrow = n, ncol = G
    )

  w <- tmp_pred_sqrt_sg^(-4)

  if(stabilize){
    y_bar <- r_tilda

    max_ybar <- max(y_bar)
    min_ybar <- min(y_bar)

    s_g_max_ybar <- predict(spline_fit, x=max_ybar)$y
    s_g_min_ybar <- predict(spline_fit, x=min_ybar)$y

    #this matrix has same dimensions of lambda_hat
    tmp_pred_sqrt_sg <- predict(
      spline_fit,
      x = c(lambda_hat)
    )$y |>
      matrix(
        nrow = n, ncol = G
      )

    #constrain observation weights with more extreme lambda hat than r_tilda
    count_changes <- 0
    for (i in seq_len(nrow(lambda_hat))) {
      for (j in seq_len(ncol(lambda_hat))) {
        #if observation is greater than max_ybar, change the weight matrix
        if(lambda_hat[i,j] > max_ybar){
          count_changes <- count_changes + 1
          tmp_pred_sqrt_sg[i,j] <- s_g_max_ybar + runif(1,0,0.2)
        }
        #if observation is less than min_ybar, change the weight matrix
        if(lambda_hat[i,j] < min_ybar){
          count_changes <- count_changes + 1
          tmp_pred_sqrt_sg[i,j] <- s_g_min_ybar + runif(1,0,0.2)
        }
      }
    }

    message(count_changes/(n*G)*100,
            "% of observations had their weight stabilized")

    w <- tmp_pred_sqrt_sg^(-4)

  }

  return(w)
}

# ---------
# run weighted nnSVG
# ---------

weights <- generate_weights_noise(
  input = spe,
  stabilize = TRUE,
  BPPARAM = MulticoreParam(workers = 10, RNGseed = 4))

spe <- weighted_nnSVG(
  input = spe,
  w = weights,
  BPPARAM = MulticoreParam(workers = 10, RNGseed = 5))

print(spe)

# -----------
# save objects
# -----------

file = here("outputs","results","spe_humanSubtypeBreast_noise_weighted_nnSVG.rds")
saveRDS(spe, file = file)

file = here("outputs","results","spe_humanSubtypeBreast_noise_weights.rds")
saveRDS(weights, file = file)

# -----------
# session information
# -----------

sessionInfo()


