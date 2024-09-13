# from https://github.com/pinellolab/SVG_Benchmarking/blob/main/simulate_data/07_scDesign3.ipynb

library(scDesign3)
library(scales)
library(dplyr)
library(ggplot2)
library(cowplot)


example_sce <- readRDS((url("https://figshare.com/ndownloader/files/40582019")))
print(example_sce)


mt_idx<- grep("mt-",rownames(example_sce))
if(length(mt_idx)!=0){
  example_sce <- example_sce[-mt_idx,]
}


set.seed(1)
example_data <- construct_data(
  sce = example_sce,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = NULL,
  spatial = c("spatial1", "spatial2"),
  other_covariates = NULL,
  corr_by = "1"
)


example_marginal <- fit_marginal(
  data = example_data,
  predictor = "gene",
  mu_formula = "s(spatial1, spatial2, bs = 'gp', k=50)",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 10,
  usebam = FALSE
)


set.seed(1)
example_copula <- fit_copula(
  sce = example_sce,
  assay_use = "counts",
  marginal_list = example_marginal,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 10,
  input_data = example_data$dat
)


example_para <- extract_para(
  sce = example_sce,
  marginal_list = example_marginal,
  n_cores = 10,
  family_use = "nb",
  new_covariate = example_data$new_covariate,
  data = example_data$dat
)


dev_explain <- sapply(example_marginal, function(x){
  sum = summary(x$fit)
  return(sum$dev.expl)
})
dev_ordered <- order(dev_explain, decreasing = TRUE)


df_dev_explain <- as.data.frame(dev_explain)
df_dev_explain$gene <- rownames(df_dev_explain)
df_dev_explain$ranking <- rank(-df_dev_explain$dev_explain)

write.csv(df_dev_explain, 'dev_explain.csv')
head(df_dev_explain)


options(repr.plot.width = 8, repr.plot.height = 5)

p <- ggplot(data = df_dev_explain, aes(x = ranking, y = dev_explain)) +
  geom_point() +
  geom_vline(xintercept = 50) +
  geom_vline(xintercept = 100) +
  geom_vline(xintercept = 150) +
  geom_vline(xintercept = 200) +
  theme_cowplot() +
  xlab("Genes") + ylab("deviation explained by GP") +
  theme(axis.text.x = element_blank())

head(df_dev_explain)

for(num_de in c(50, 100, 150, 200)){
  ordered <- dev_explain[dev_ordered]
  de_idx <- names(ordered)[1:num_de]
  non_de_idx <- names(ordered)[-(1:num_de)]
  non_de_mat <- apply(example_para$mean_mat[, non_de_idx], 2, function(x){
    avg <- (max(x)+min(x))/2
    new_mean <- rep(avg, length(x))
    return(new_mean)
  })
  example_para$mean_mat[, non_de_idx] <- non_de_mat
  
  set.seed(1)
  example_newcount <- simu_new(
    sce = example_sce,
    mean_mat = example_para$mean_mat,
    sigma_mat = example_para$sigma_mat,
    zero_mat = example_para$zero_mat,
    quantile_mat = NULL,
    copula_list = example_copula$copula_list,
    n_cores = 10,
    family_use = "nb",
    input_data = example_data$dat,
    new_covariate = example_data$newCovariate,
    important_feature = rep(TRUE, dim(example_sce)[1]),
    filtered_gene = NULL
  )
  
  simu_sce <- SingleCellExperiment(list(counts =example_newcount), colData = example_data$newCovariate)
  loc = colData(simu_sce)[,c("spatial1","spatial2")]
  
  write.csv(loc, file = "location.csv")
  write.csv(example_newcount, file = glue::glue('counts_svg_{num_de}.csv'))
  write.csv(de_idx, file = glue::glue('.labels_svg_{num_de}.csv'))
  write.csv(non_de_idx, file = glue::glue('labels_non_svg_{num_de}.csv'))
  
}

