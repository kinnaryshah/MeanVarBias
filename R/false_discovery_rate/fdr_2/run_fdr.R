library(SpatialExperiment)
library(STexampleData)
library(MASS)
library(scuttle)

set.seed(12)

### SIMULATE DATA
# beta range 0.5 to 1
# sigma.sq range 1 to 3

n_genes <- 300
fraction <- 0.5
max_sigma.sq <- 3
low_range_beta <- c(0.5,1)

#check if integer
stopifnot(n_genes*fraction*0.5 == round(n_genes*fraction*0.5))

#some genes have some nonzero sigma.sq
sigma.sq <- runif(n_genes, 1, max_sigma.sq)
#some genes have zero sigma.sq
sigma.sq[sample(1:n_genes, round(n_genes*fraction*0.5))] <- 0
ground_truth_rank <- rank(-sigma.sq)

#all genes will have nonzero beta values
beta <- runif(n_genes, log(low_range_beta[1]), log(low_range_beta[2]))

#choose fixed length scale parameter (~medium from nnSVG paper)

scale_length <- 200

params <- data.frame(sigma.sq, beta)

#sampling from a poisson distribution - mean controls variance, so we don't specify tau.sq:
#step 1: use ST example distance matrix instead of creating a new one (Euclidean distance)

spe_demo <- Visium_humanDLPFC()
points_coord <- spatialCoords(spe_demo)
n_points <- nrow(points_coord)

pair.points <- cbind(
  matrix( rep(points_coord, each = n_points), ncol = 2, byrow = FALSE),
  rep(1, times = n_points) %x% points_coord # Creating the combinations using kronecker product.
) |> data.frame()
colnames(pair.points) <- c("si.x", "si.y", "sj.x", "sj.y")

#step 2: calculate gaussian process/kernel 

kernel.fun <- function(si.x, si.y, sj.x, sj.y,  l = 0.2){
  exp(-1*sqrt(((si.x-sj.x)^2+(si.y-sj.y)^2))/l)
}

C_theta <- with(pair.points, kernel.fun(si.x, si.y, sj.x, sj.y, l = scale_length)) |> 
  matrix(nrow = n_points, ncol = n_points)

counts <- matrix(NA, nrow = n_genes, ncol = n_points)
eta_list <- list()

for (i in c(1:n_genes)) {
  
  print(i)
  sigma.sq_i <- sigma.sq[i]
  beta_i <- beta[i]
  print(paste0(sigma.sq_i, " ", beta_i))
  
  #step 3: simulate gaussian process per gene
  
  gp_dat <- mvrnorm(n = 1, rep(0,n_points), sigma.sq_i* C_theta) 
  
  #step 4: calculate lambda = exp(beta + gaussian process) per gene
  
  #eta_i <- mean(gp_dat + beta_i)
  #eta_list <- append(eta_list, eta_i)
  #lambda_i <- exp(eta_i)
  lambda_i <- exp(gp_dat + beta_i)
  
  #step 5: use rpois() to simulate 4992 values per gene
  
  counts_i <- rpois(n = n_points, lambda_i)
  
  #put all counts in matrix 
  #orientation: genes x spots
  
  counts[i,] <- counts_i
}

#create spe using counts and distance matrix

spe <- SpatialExperiment(
  assays = list(counts = counts),
  spatialCoords = points_coord)

rowData(spe)$ground_truth_rank <- ground_truth_rank
rowData(spe)$ground_truth_sigma.sq <- sigma.sq
rowData(spe)$ground_truth_beta <- beta

saveRDS(spe, file = "spe_simulation.rds")

### RUN NNSVG

library(nnSVG)
library(SpatialExperiment)
library(BRISC)
library(BiocParallel)
library(scuttle)
library(Matrix)
library(spoon)

print(dim(spe))

spe <- logNormCounts(spe)

runtime <- system.time({
  tryCatch({
    spe <- nnSVG(
      spe,
      X = NULL,
      assay_name = "logcounts",
      n_neighbors = 10,
      order = "AMMD",
      n_threads = 10,
      verbose = FALSE
    )
  }, error = function(e) {
    print(paste("Error:", e$message))
  })
})

file = "spe_nnSVG.rds"
saveRDS(spe, file = file)

spe_unweighted <- spe

### RUN WEIGHTED NNSVG

fn <- ("spe_simulation.rds") 

spe <- readRDS(fn)

spe <- logNormCounts(spe)

weights <- generate_weights(input = spe,
                            stabilize = TRUE,
                            BPPARAM = MulticoreParam(workers = 10,
                                                     RNGseed = 4))

spe <- weighted_nnSVG(input = spe,
                      w = weights,
                      BPPARAM = MulticoreParam(workers = 10, RNGseed = 5))

file = "spe_weighted_nnSVG.rds"
saveRDS(spe, file = file)

spe_weighted <- spe

### CREATE FDR PLOTS

alpha_levels <- seq(0.01, 0.1, by = 0.01)
fdr_unweighted <- numeric(length(alpha_levels))
fdr_weighted <- numeric(length(alpha_levels))

for (i in 1:length(alpha_levels)) {
  alpha <- alpha_levels[i]
  
  # UNWEIGHTED
  
  rowData(spe_unweighted)$pos_test = rowData(spe_unweighted)$padj < alpha
  rowData(spe_unweighted)$SVG = rowData(spe_unweighted)$ground_truth_sigma.sq != 0
  
  total_genes_unweighted <- dim(spe_unweighted)[1]
  false_positives_unweighted <- sum(rowData(spe_unweighted)$pos_test & !rowData(spe_unweighted)$SVG)
  
  fdr_unweighted[i] <- false_positives_unweighted / sum(rowData(spe_unweighted)$pos_test)
  
  # WEIGHTED
  
  pvalues_weighted <- 1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df = 2)
  pvalues_weighted <- p.adjust(pvalues_weighted, method = "BH")
  rowData(spe_weighted)$pos_test = pvalues_weighted < alpha
  rowData(spe_weighted)$SVG = rowData(spe_weighted)$ground_truth_sigma.sq != 0
  
  total_genes_weighted <- dim(spe_weighted)[1]
  false_positives_weighted <- sum(rowData(spe_weighted)$pos_test & !rowData(spe_weighted)$SVG)
  
  fdr_weighted[i] <- false_positives_weighted / sum(rowData(spe_weighted)$pos_test)
}

# plot using ggplot2
df <- data.frame(alpha = alpha_levels, fdr_unweighted = fdr_unweighted, fdr_weighted = fdr_weighted)
df <- tidyr::gather(df, model, FDR, -alpha)

ggplot(df, aes(x = alpha, y = FDR, color = model)) +
  geom_line() +
  geom_point() +
  labs(title = "FDR vs. Alpha Levels",
       x = "Alpha Level",
       y = "False Discovery Rate",
       color = "Model") +
  theme_bw()

ggsave("FDR_vs_alpha.png")
