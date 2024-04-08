library(SpatialExperiment)
library(dplyr)
library(ggplot2)

# make a function with input of spe_unweighted and spe_weighted, return a data frame of fdr_unweighted and fdr_weighted 
create_FDR_df <- function(spe_unweighted, spe_weighted) {
  
  # create a loop from the FDR calculation for varying alpha levels from 0.01 to 0.1
  # then plot the FDRs for unweighted and weighted models
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
    false_positives_weighted <- sum(rowData(spe_weighted)$pos_test & !rowData(spe_weighted)$SVG, na.rm = TRUE)
    
    fdr_weighted[i] <- false_positives_weighted / sum(rowData(spe_weighted)$pos_test, na.rm = TRUE)
  }
  
  df <- data.frame(alpha = alpha_levels, fdr_unweighted = fdr_unweighted, fdr_weighted = fdr_weighted)
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
  df_temp$fdr <- i+1
  
  df <- rbind(df, df_temp)
}

# add another df_temp to df using i=1 from the folder "../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/"
spe_unweighted <- readRDS("../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/spe_nnSVG.rds")
spe_weighted <- readRDS("../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/spe_weighted_nnSVG.rds")

df_temp <- create_FDR_df(spe_unweighted, spe_weighted)
df_temp$fdr <- 1

df <- rbind(df, df_temp)

# plot the FDRs for unweighted and weighted models
# separate plots for each fdr 1 to 9
# use faceting from ggplot2
# make same y-axis for all the plots

p <- ggplot(df, aes(x = alpha)) +
  geom_line(aes(y = fdr_unweighted, color = "Unweighted")) +
  geom_line(aes(y = fdr_weighted, color = "Weighted")) +
  labs(title = "FDR vs alpha",
       x = "alpha",
       y = "FDR") +
  facet_wrap(~fdr) +
  theme_minimal() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  scale_y_continuous(limits = c(0, 0.2)) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "increasing beta range -->", breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "increasing sigma.sq range -->", breaks = NULL, labels = NULL))

# save the plot
ggsave("combined_FDR_vs_alpha.png", p, width = 10, height = 10)
