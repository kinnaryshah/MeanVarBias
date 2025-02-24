library(ggplot2)
library(SpatialExperiment)
library(patchwork)
library(GGally)
library(dplyr)
library(ggridges)
library(here)


My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  legend.title = element_text(size=20),
  legend.text = element_text(size=20),
  title = element_text(size = 20))

# plot each gene's average weight on each axis, colored by mean expression level

# download files
file <- here("outputs","results","spe_humanDLPFC_voom_weights.rds")
voom_DLPFC <- readRDS(file = file)
voom_DLPFC <- t(voom_DLPFC)

file <- here("outputs","results","spe_humanDLPFC_weights.rds")
spoon_DLPFC <- readRDS(file = file)
spoon_DLPFC <- t(spoon_DLPFC)

# get average weight for each gene
voom_weight <- rowMeans(voom_DLPFC)
spoon_weight <- rowMeans(spoon_DLPFC)

# get mean expression for each gene from spe 
file = here("outputs","results","spe_humanDLPFC_nnSVG.rds")
spe <- readRDS(file = file)
mean_expr <- rowData(spe)$mean

# make plot
df <- data.frame(
  voom_weight = voom_weight,
  spoon_weight = spoon_weight,
  mean_expr = mean_expr
)

p <- ggplot(data = df, aes(x=voom_weight, y=spoon_weight, color = mean_expr)) +
  geom_point(size=1, alpha = 0.3) +
  scale_colour_gradient(low = "green", high = "purple") +
  labs(
    y = "Average spoon Weight",
    x = "Average voom Weight",
    title = "Comparing voom and spoon weights",
    color = "Mean Expression"
  ) +
  theme_bw() +
  My_Theme

# could plot each observation's weight on each axis, colored by expression level of the gene it comes from?

ggsave(here("plots", "supplementary", "comparing_voom_spoon.png"),
       wrap_plots(p,
                  nrow = 1,
                  guides = "collect"),
       width = 15, height = 15)
