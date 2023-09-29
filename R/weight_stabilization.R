library(ggplot2)
library(ggformula)
library(purrr)
#load spline_fit, w matrix, r_tilda, s_g, lambda_hat
load(file = "mean_var_project/simulation/poisson_3/simulation_BRISC_estimation_spline.rds")

#get dimensions of spe
spe <- readRDS(file = "mean_var_project/simulation/poisson_3/spe_simulation.rds")
n <- dim(spe)[2] # Number of Cells
G <- dim(spe)[1] # Number of Genes

#plot r_tilda vs sqrt(s_g)
p <- data.frame(x = r_tilda, y = sqrt(s_g)) |>
  ggplot() +
  stat_spline(aes(x = x, y = y), color = "blue", linewidth = 3) +
  labs(
    x = "log2(count size)",
    y = "Sqrt(s_g)"
  )

plot(new_y ~ new_x, xlim = c(0,2), ylim = c(0,5))
plot.new()
plot(sqrt(s_g) ~ r_tilda)
lines(spline_fit, lty = 2, col = "red")

fivenum(w)
#[1] 6.458559e-03 1.189445e+00 2.189256e+00 3.267646e+00 1.581570e+19

fivenum(r_tilda)
#[1] -0.08877544  0.11008756  1.26336581  2.38983768  2.95709919

#find min and max 
y_bar <- r_tilda

max_ybar <- max(y_bar)
min_ybar <- min(y_bar)

max_ybar; min_ybar
#2.957099
#-0.08877544

s_g_max_ybar <- predict(spline_fit, x=max_ybar)$y
s_g_min_ybar <- predict(spline_fit, x=min_ybar)$y

s_g_max_ybar; s_g_min_ybar #these values are very close to the actual sqrt(s_g) values corresponding to the max and min r_tilda values
#0.6745045
#0.8839731

#this matrix has same dimensions of lambda_hat
tmp_pred_sqrt_sg <- predict(
  spline_fit, 
  x = c(lambda_hat)
)$y |> 
  matrix(
    nrow = n, ncol = G
  )

#plot log counts vs s_g before constraining
df <- data.frame(y = exec(c, tmp_pred_sqrt_sg, use.names = FALSE),
                 x = exec(c, lambda_hat, use.names = FALSE))
p + geom_point(data = df, aes(x, y), size = 0.2)

#constrain individual observation weights that have lambda hat more extreme than range of r_tilda 
count_changes <- 0
for (i in 1:nrow(lambda_hat)) {
  print(i)
  for (j in 1:ncol(lambda_hat)) {
    print(j)
    #if this observation is greater than the max_ybar, change the weight matrix
    if(lambda_hat[i,j] > max_ybar){
      count_changes <- count_changes + 1
      tmp_pred_sqrt_sg[i,j] <- s_g_max_ybar
    } 
    #if this observation is less than the min_ybar, change the weight matrix
    if(lambda_hat[i,j] < min_ybar){
      count_changes <- count_changes + 1
      tmp_pred_sqrt_sg[i,j] <- s_g_min_ybar
    }
  }
}

count_changes/(n*G)
#0.2540655

w_new <- tmp_pred_sqrt_sg^(-4) 
fivenum(w_new)
#[1] 0.7213785  1.6377358  2.1905490  3.2739918 32.2346907

#plot log counts vs s_g after constraining
#make matrix into one long vector to plot each observation as a point
df <- data.frame(y = exec(c, tmp_pred_sqrt_sg, use.names = FALSE),
                 x = exec(c, lambda_hat, use.names = FALSE))
p + geom_point(data = df, aes(x, y), size = 0.2)
