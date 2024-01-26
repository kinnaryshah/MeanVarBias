library(SpatialExperiment)
library(STexampleData)
library(MASS)
library(scuttle)

set.seed(1)

### GENERATING DATA

#4992 spots and 300 genes

n_genes <- 300
fraction <- 0.5
max_sigma.sq <- 1
low_range_beta <- c(0.5,1)

#check if integer
stopifnot(n_genes*fraction*0.5 == round(n_genes*fraction*0.5))

#all genes have some nonzero sigma.sq
sigma.sq <- runif(n_genes, 0.2, max_sigma.sq)
ground_truth_rank <- rank(-sigma.sq)

#all genes will have nonzero beta values
beta <- runif(n_genes, log(low_range_beta[1]), log(low_range_beta[2]))

#choose fixed length scale parameter (~medium from nnSVG paper)

scale_length <- 200

params <- data.frame(sigma.sq, beta)

plot(beta, sigma.sq)

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

pdf("beta_eta_simulation.pdf", width = 21, height = 10)
plot(beta, unlist(eta_list))
dev.off()

#create spe using counts and distance matrix

spe <- SpatialExperiment(
    assays = list(counts = counts),
    spatialCoords = points_coord)

rowData(spe)$ground_truth_rank <- ground_truth_rank
rowData(spe)$ground_truth_sigma.sq <- sigma.sq
rowData(spe)$ground_truth_beta <- beta

saveRDS(spe, file = "spe_simulation.rds")
