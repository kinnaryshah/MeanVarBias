n <- 200
x <- runif(n, 0, 1)
error <- rnorm(n, 0, 0.1)
y <- x^2 + error
#plot(y ~ x, main = "weight stabilization simuation")

spline_fit <- smooth.spline(y=y, x=x, df = 5)

new_x <- runif(50, 1, 2)
new_error <- rnorm(50, 0, 0.1)
new_y <- new_x^2 + new_error
plot(new_y ~ new_x, xlim = c(0,2), ylim = c(0,5))
lines(spline_fit, lty = 2, col = "red")

count_changes <- 0
for (i in c(1:length(new_x))) {
  print(i)
  #if this observation is greater than the max x, change the new_y value 
  if(new_x[i] > max(x)){
    count_changes <- count_changes + 1
    new_y[i] <- predict(spline_fit, x=max(x))$y
  } 
}

count_changes
plot(new_y ~ new_x, xlim = c(0,2), ylim = c(0,5))
lines(spline_fit, lty = 2, col = "red")

