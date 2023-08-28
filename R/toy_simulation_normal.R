set.seed(1)


library(SpatialExperiment)
library(STexampleData)
library(MASS)
library(scuttle)

### GENERATING DATA

#4992 spots and 2 genes

n_genes <- 4

sigma.sq <- c(2, 1.2, 0, 1.5)
tau.sq <- c(0.2, 0.15, 0.25, 0.1)

ground_truth <- c(T, T, F, T)
ground_truth_rank <- rank(-sigma.sq)

beta <- c(9, 10, 8, 7)

#choose fixed length scale parameter (~medium from nnSVG paper)

scale_length <- 200

params <- data.frame(sigma.sq, beta)

#plot(beta, sigma.sq)

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

kernel.fun <- function(si.x, si.y, sj.x, sj.y, l){
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
  lambda_i <- gp_dat + beta_i
  
  #step 5: use rpois() to simulate 4992 values per gene
  
  #counts_i <- rpois(n = n_points, lambda_i)
  counts_i <- rnorm(n = n_points, lambda_i, tau.sq[i]) 
  
  #put all counts in matrix 
  #orientation: genes x spots
  
  counts[i,] <- counts_i
}

#pdf("beta_eta_simulation.pdf", width = 21, height = 10)
#plot(beta, unlist(eta_list))
#dev.off()

#create spe using counts and distance matrix

spe <- SpatialExperiment(
  assays = list(counts = counts),
  spatialCoords = points_coord)

rowData(spe)$ground_truth <- ground_truth
rowData(spe)$ground_truth_rank <- ground_truth_rank
rowData(spe)$ground_truth_sigma.sq <- sigma.sq
rowData(spe)$ground_truth_beta <- beta



library(nnSVG)
library(SpatialExperiment)
library(scuttle)

#run nnSVG on the logcounts matrix

spe <- nnSVG(spe, assay_name = "counts")



#~~~~~~~~~~~~~~~~~
# checking if we can recover true sigma.sq value using BRISC
#~~~~~~~~~~~~~~~~~

library(SpatialExperiment)
library(STexampleData)
library(MASS)
library(scuttle)
library(BRISC)

n_genes <- 1

sigma.sq <- 3
tau.sq <- 0.1
beta <- 0
scale_length <- 200

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

kernel.fun <- function(si.x, si.y, sj.x, sj.y,  l){
  exp(-1*sqrt(((si.x-sj.x)^2+(si.y-sj.y)^2))/l)
}

C_theta <- with(pair.points, kernel.fun(si.x, si.y, sj.x, sj.y, l = scale_length)) |> 
  matrix(nrow = n_points, ncol = n_points)

#step 3: simulate gaussian process per gene

gp_dat <- mvrnorm(n = 1, rep(0,n_points), sigma.sq* C_theta) 

#step 4: calculate lambda = exp(beta + gaussian process) per gene

eta <- mean(gp_dat + beta)
#lambda <- exp(eta)
lambda <- gp_dat + beta

#step 5: use rpois() to simulate 4992 values per gene

#counts <- rpois(n = n_points, lambda)
counts <- rnorm(n = n_points, lambda, tau.sq)

estimation_result <- BRISC_estimation(coords = points_coord, y = counts)

estimation_result$Theta ##Estimates of covariance model parameters.
estimation_result$Beta ##Estimates of Beta
