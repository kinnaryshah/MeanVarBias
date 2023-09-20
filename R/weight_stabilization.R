w <- readRDS(file = "simulation_BRISC_estimation_spline.rds")

fivenum(w)
#[1] 9.639114e-03 1.254207e+00 2.158240e+00 3.157108e+00 5.869133e+17

#get mean of weights per gene
y_bar <- colMeans(w)

fivenum(y_bar)
#[1] 1.541238e+00 2.052111e+00 2.768669e+00 1.568808e+01 1.176866e+14

#constrain individual observation weights that are greater than the mean to the mean
count_changes <- 0
for (i in 1:ncol(w)) {
  print(i)
  for (j in nrow(w)) {
    print(j)
    #if the weight for this observation is greater than the mean for this gene, change
    if(w[j,i] > y_bar[i]){
      count_changes <- count_changes + 1
      w[j,i] <- y_bar[i]
    }
  }
}

fivenum(w)
#[1] 9.639114e-03 1.254207e+00 2.158190e+00 3.157056e+00 5.869133e+17