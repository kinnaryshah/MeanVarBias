library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(viridis)
library(ggpubr)
library(purrr)
library(patchwork)

df_fxn <- function(file_name, dataset) {
  
  spe <- readRDS(file = file_name)
  
  spe_df <- rowData(spe)
  
  df_nnSVG <- 
    as.data.frame(spe_df)
  
  df_effect <- 
    as.data.frame(df_nnSVG) %>% 
    mutate(l = 1 / phi) %>% 
    filter(rank <= 1000)
  
  df_effect$LR_stat_scaled <- df_effect$LR_stat / max(df_effect$LR_stat)
  
  #add col of dataset name
  df_effect$dataset_val <- rep(dataset,1000)
  
  #keep only nnSVG output because breast cancer dataset has different cols
  df_effect <- df_effect %>%
    select(mean,var,LR_stat,LR_stat_scaled,dataset_val,sigma.sq,tau.sq,prop_sv)
  
  return(df_effect)
}

file_list <- c(here("outputs", "results", "spe_humanHPC_nnSVG.rds"),
               here("outputs", "results", "spe_humanBreast_nnSVG.rds"),
               here("outputs", "results", "spe_humanDLPFC_nnSVG.rds"),
               here("outputs", "results", "spe_humanLC_nnSVG.rds"),
               here("outputs", "results", "spe_humanOvarian_nnSVG.rds"))

dataset_list <- c("HPC", "Breast","DLPFC","LC", "Ovarian")

df <- map2(.x=file_list, .y=dataset_list, .f=~df_fxn(.x,.y))
reduce_df <- Reduce(rbind,df)

# variance vs. mean
var <- ggplot(reduce_df, 
              aes(x = mean, y = var, color = LR_stat_scaled)) + 
  geom_point(size = 1) + 
  geom_smooth(method="loess", color="black", linewidth=0.5) +
  facet_grid(dataset_val~., switch = "y") +
  scale_color_viridis(trans = "log10") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "mean logcounts", 
       y = "", 
       color = "LR stat",
       title = expression(sigma^{2}+tau^{2})) + 
  theme_bw() + 
  theme(strip.text.y.left = element_text(angle = 0)) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="black", fill="white") ) +
  theme(axis.title.x = element_text(size = 7))


# spatial variance vs. mean
spat_var <- ggplot(reduce_df, 
                   aes(x = mean, y = sigma.sq, color = LR_stat_scaled)) + 
  geom_point(size = 1) + 
  facet_grid(dataset_val~., switch = "y") +
  scale_color_viridis(trans = "log10") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "mean logcounts", 
       y = "", 
       color = "LR stat",
       title = expression(sigma^{2})) + 
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text.y = element_blank()) +
  theme(axis.title.x = element_text(size = 7))

nonspat_var <- ggplot(reduce_df, 
                      aes(x = mean, y = tau.sq, color = LR_stat_scaled)) + 
  geom_point(size = 1) + 
  facet_grid(dataset_val~., switch = "y") +
  scale_color_viridis(trans = "log10") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "mean logcounts", 
       y = "", 
       color = "LR stat",
       title = expression(tau^{2})) + 
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text.y = element_blank()) +
  theme(axis.title.x = element_text(size = 7))

prop <- ggplot(reduce_df, 
               aes(x = mean, y = prop_sv, color = LR_stat_scaled)) + 
  geom_point(size = 1) + 
  facet_grid(dataset_val~., switch = "y") +
  scale_color_viridis(trans = "log10") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "mean logcounts", 
       y = "", 
       color = "LR stat",
       title = expression(sigma^{2}/(sigma^{2}+tau^{2}))) + 
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text.y = element_blank()) +
  theme(axis.title.x = element_text(size = 7))

ggsave(here("plots", "main", "motivation.png"),
       wrap_plots(var, spat_var, nonspat_var, prop, guides="collect",
           ncol=4, nrow=1))

