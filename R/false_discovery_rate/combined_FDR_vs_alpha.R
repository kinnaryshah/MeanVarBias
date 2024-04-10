library(SpatialExperiment)
library(dplyr)
library(ggplot2)

# make a function with input of spe_unweighted and spe_weighted, return a data frame of fdr_unweighted and fdr_weighted 
create_FDR_df <- function(spe_unweighted, spe_weighted) {
  
  # UNWEIGHTED
  
  # Sort unadjusted p-values and calculate FDR
  sorted_p_values_unweighted <- sort(rowData(spe_unweighted)$pval)
  
  # WEIGHTED
  
  # Sort unadjusted p-values and calculate FDR
  pvalues_weighted <- 1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df = 2)
  sorted_p_values_weighted <- sort(pvalues_weighted, na.last = TRUE)
  
  num_tests <- length(sorted_p_values_unweighted)
  fdr <- (1:num_tests) * 0.05 / num_tests
  
  df <- data.frame(p_vals_unweighted = sorted_p_values_unweighted,
                   p_vals_weighted = sorted_p_values_weighted,
                   fdr = fdr)
  
  return(df)
}


# use the function to create the data frame for the files from all the folders from fdr_1 to fdr_8
# create a column for the fdr_1 to fdr_8 

df <- data.frame()
for (i in 1:8) {
  print(i)
  spe_unweighted <- readRDS(paste0("fdr_", i, "/spe_nnSVG.rds"))
  spe_weighted <- readRDS(paste0("fdr_", i, "/spe_weighted_nnSVG.rds"))
  
  df_temp <- create_FDR_df(spe_unweighted, spe_weighted)
  df_temp$sim <- i+1
  
  df <- rbind(df, df_temp)
}

# add another df_temp to df using i=1 from the folder "../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/"
spe_unweighted <- readRDS("../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/spe_nnSVG.rds")
spe_weighted <- readRDS("../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/spe_weighted_nnSVG.rds")

df_temp <- create_FDR_df(spe_unweighted, spe_weighted)
df_temp$sim <- 1

df <- rbind(df, df_temp)

# plot the FDRs for unweighted and weighted models
# separate plots for each fdr 1 to 9
# use faceting from ggplot2
# x axis is the significance level
# y axis is the FDR
# make same y-axis for all the plots

p <- ggplot(df, aes(x = p_vals_unweighted, y = fdr)) +
  geom_line(aes(color = "Unweighted")) +
  geom_line(aes(x = p_vals_weighted, y = fdr, color = "Weighted")) +
  facet_wrap(~sim, scales = "free_y") +
  labs(x = "Significance level", y = "FDR") +
  theme_bw()

# save the plot
ggsave("combined_FDR_vs_significance.png", p, width = 10, height = 10)
