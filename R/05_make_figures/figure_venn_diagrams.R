# Load the required libraries
library(ggplot2)
library(eulerr)
library(gridExtra)

# Define the first set of scenarios
scenarios_row1 <- list(
  c(A = 5, "A&B" = 2),  # 1: 7 genes in List 1, 2 overlaps
  c(A = 20, "A&B" = 11), # 2: 31 genes in List 1, 11 overlaps
  c(A = 18, "A&B" = 8),   # 3: 26 genes in List 1, 8 overlaps
  c(A = 1, "A&B" = 0)      # 4: 1 gene in List 1, 0 overlaps
)

# Define the second set of scenarios
scenarios_row2 <- list(
  c(A = 69, "A&B" = 8),   # 1: 77 genes in List 1, 8 overlaps
  c(A = 65, "A&B" = 19),  # 2: 84 genes in List 1, 19 overlaps
  c(A = 59, "A&B" = 19),  # 3: 78 genes in List 1, 19 overlaps
  c(A = 43, "A&B" = 16)   # 4: 59 genes in List 1, 16 overlaps
)

# Create each plot for the first row without set labels
plots_row1 <- lapply(scenarios_row1, function(scenario) {
  fit <- euler(scenario)
  plot(fit, 
       quantities = TRUE, 
       labels = FALSE, 
       fills = c("red", "green"), 
       alpha = 0.5,
       labels_cex = 0.9)
})

# Create each plot for the second row without set labels
plots_row2 <- lapply(scenarios_row2, function(scenario) {
  fit <- euler(scenario)
  plot(fit, 
       quantities = TRUE, 
       labels = FALSE, 
       fills = c("red", "green"), 
       alpha = 0.5,
       labels_cex = 0.9)
})

# Arrange the plots in a 2x4 grid layout
do.call(grid.arrange, c(plots_row1, plots_row2, nrow = 2))


do.call(grid.arrange, c(plots_row1,nrow = 1))
do.call(grid.arrange, c(plots_row2,nrow = 1))
