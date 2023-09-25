library(scuttle)
library(ggplot2)

#load spline_fit, w matrix, and spe
load(file = "mean_var_project/simulation/poisson_3/simulation_BRISC_estimation_spline.rds")
spe <- readRDS(file = "mean_var_project/simulation/poisson_3/spe_simulation.rds")

spe <- spe[, colSums(counts(spe)) > 0]
dim(spe)

spe <- logNormCounts(spe)

#plot log counts vs weights
data.frame(x = y_bar,
           y = rowMeans(w)) |>
  ggplot() +
  aes(x=x, y=y) + 
  scale_y_continuous(trans = "log10") +
  geom_line(aes(x=x, y=y)) +
  geom_point() +
  labs(x="logcounts", y="mean weight")
  

fivenum(w)
#[1] 9.639114e-03 1.254207e+00 2.158240e+00 3.157108e+00 5.869133e+17

fivenum(logcounts(spe))
#[1] 0.0000000 0.9467601 1.9281053 2.6885182 5.8868834

#get mean of logcounts per gene
y_bar <- colMeans(logcounts(spe))

fivenum(y_bar)
#[1] 1.575665 1.668835 1.688914 1.707375 1.802108

max_ybar <- max(y_bar)
min_ybar <- min(y_bar)

max_ybar; min_ybar
#[1] 1.802108
#[1] 1.575665

w_max_ybar <- predict(spline_fit, x=max_ybar)$y
w_min_ybar <-predict(spline_fit, x=min_ybar)$y

w_max_ybar; w_min_ybar
#[1] 1.026443
#[1] 0.9774143

#constrain individual observations that are greater than the mean to the weight of the mean 
count_changes <- 0
for (i in 1:ncol(w)) {
  print(i)
  for (j in nrow(w)) {
    print(j)
    #if this observation is greater than the max_ybar mean, change the weight matrix
    if(logcounts(spe)[i,j] > max_ybar){
      count_changes <- count_changes + 1
      w[j,i] <- w_max_ybar
    } 
    #if this observation is less than the min_ybar mean, change the weight matrix
    if(logcounts(spe)[i,j] < min_ybar){
      count_changes <- count_changes + 1
      w[j,i] <- w_min_ybar
    }
  }
}

count_changes
#180

fivenum(w)
#[1] 1.282067e-02 1.360334e+00 2.183129e+00 3.128545e+00 5.183556e+14

#plot log counts vs weights after constraining
