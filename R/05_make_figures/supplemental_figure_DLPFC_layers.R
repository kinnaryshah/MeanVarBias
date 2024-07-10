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

df_fxn <- function(layer) {
  
  file = here("outputs","results",paste0("spe_humanDLPFC_layer_", layer, "_nnSVG.rds"))
  spe <- readRDS(file = file)
  
  spe_df <- rowData(spe)
  
  df_nnSVG <- 
    as.data.frame(spe_df)
  
  df_effect <- 
    as.data.frame(df_nnSVG) %>% 
    mutate(l = 1 / phi) %>% 
    filter(rank <= 1000)
  
  df_effect$LR_stat_scaled <- df_effect$LR_stat / max(df_effect$LR_stat)
  
  #add col of layer name
  df_effect$layer_val <- rep(layer,1000)
  return(df_effect)
}

layer_list <- c("L1", "L2", "L3", "L4", "L5", "L6", "WM")
#create dfs of all layers
layer_df <- map(.x=layer_list, .f=~df_fxn(.x))
#convert to df and rbind, will facet on layer column which was added in df_fxn
all_layers_df <- Reduce(rbind,layer_df)

# variance vs. mean
var <- ggplot(all_layers_df, 
              aes(x = mean, y = var, color = LR_stat_scaled)) + 
  geom_point(size = 1) + 
  geom_smooth(method="loess", color="black", linewidth=0.5) +
  facet_grid(layer_val~., switch = "y") +
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
spat_var <- ggplot(all_layers_df, 
                   aes(x = mean, y = sigma.sq, color = LR_stat_scaled)) + 
  geom_point(size = 1) + 
  facet_grid(layer_val~., switch = "y") +
  scale_color_viridis(trans = "log10") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "mean logcounts", 
       y = "", 
       color = "LR stat",
       title = expression(sigma^{2})) + 
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text.y = element_blank()) +
  theme(axis.title.x = element_text(size = 7))

nonspat_var <- ggplot(all_layers_df, 
                      aes(x = mean, y = tau.sq, color = LR_stat_scaled)) + 
  geom_point(size = 1) + 
  facet_grid(layer_val~., switch = "y") +
  scale_color_viridis(trans = "log10") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "mean logcounts", 
       y = "", 
       color = "LR stat",
       title =expression(tau^{2})) + 
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text.y = element_blank()) +
  theme(axis.title.x = element_text(size = 7))

prop <- ggplot(all_layers_df, 
               aes(x = mean, y = prop_sv, color = LR_stat_scaled)) + 
  geom_point(size = 1) + 
  facet_grid(layer_val~., switch = "y") +
  scale_color_viridis(trans = "log10") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "mean logcounts", 
       y = "", 
       color = "LR stat",
       title = expression(sigma^{2}/(sigma^{2}+tau^{2}))) + 
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text.y = element_blank()) +
  theme(axis.title.x = element_text(size = 7))

ggsave(here("plots", "supplementary", "DLPFC_layers.png"),
       wrap_plots(var, spat_var, nonspat_var, prop, guides="collect",
           ncol=4, nrow=1))
