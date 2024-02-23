library(SpatialExperiment)
library(STexampleData)
library(MASS)
library(scuttle)

spe_demo <- Visium_humanDLPFC()

set.seed(1)

### GENERATING DATA

#4992 spots and 300 genes

n_genes <- 300
fraction <- 0.6
max_sigma.sq <- 1

#some genes have some nonzero sigma.sq
#most genes have zero sigma.sq
sigma.sq <- runif(n_genes, 0.2, max_sigma.sq)
sigma.sq[sample(1:n_genes, round(n_genes*fraction))] <- 0
ground_truth_rank <- rank(-sigma.sq)

#all genes will have nonzero beta values
#sampled from DLPFC dataset to get more realistic distribution of marginal means

#normalizing (not log) transformation
spe_demo <- logNormCounts(spe_demo, transform = "none")

#calculate the marginal means of all genes
marginal_means <- rowMeans(normcounts(spe_demo))

#remove outliers
marginal_means <- marginal_means[rowMeans(normcounts(spe_demo)) < 10]

beta <- sample(log(marginal_means+1), n_genes)

#choose fixed length scale parameter (~medium from nnSVG paper)

scale_length <- 200

params <- data.frame(sigma.sq, beta)

plot(beta, sigma.sq)

#sampling from a poisson distribution - mean controls variance, so we don't specify tau.sq:
#step 1: use ST example distance matrix instead of creating a new one (Euclidean distance)

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
