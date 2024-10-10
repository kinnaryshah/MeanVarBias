library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(scales)
library(ggpubr)
library(purrr)
library(patchwork)

file_list <- c(here("outputs", "results", "spe_humanHPC_V12D07-335_D1_nnSVG.rds"),
               here("outputs", "results", "spe_humanBreast_nnSVG.rds"),
               here("outputs", "results", "spe_humanLC_nnSVG.rds"),
               here("outputs", "results", "spe_humanOvarian_nnSVG.rds"))

dataset_list <- c("HPC", "Breast","LC", "Ovarian")

spe_HPC <- readRDS(file_list[1])
spe_Breast <- readRDS(file_list[2])
spe_LC <- readRDS(file_list[3])
spe_Ovarian <- readRDS(file_list[4])

# HPC
lengthscale_data <- data.frame(lengthscale = 1 / rowData(spe_HPC)$phi)
lengthscale_data <- lengthscale_data[lengthscale_data$lengthscale < 1, ]

HPC_hist <- ggplot(as.data.frame(lengthscale_data), aes(x = lengthscale_data)) + 
  geom_density(color = "black", fill = "#804585") + 
  xlim(c(0, 1)) + 
  xlab("Estimated Length Scale") + 
  ylab("Density") + 
  ggtitle("HPC") + 
  theme_bw()

# Breast
lengthscale_data <- data.frame(lengthscale = 1 / rowData(spe_Breast)$phi)
lengthscale_data <- lengthscale_data[lengthscale_data$lengthscale < 1, ]

Breast_hist <- ggplot(as.data.frame(lengthscale_data), aes(x = lengthscale_data)) + 
  geom_density(color = "black", fill = "#804585") + 
  xlim(c(0, 1)) + 
  xlab("Estimated Length Scale") + 
  ylab("Density") + 
  ggtitle("Ductal Breast") + 
  theme_bw()

# LC
lengthscale_data <- data.frame(lengthscale = 1 / rowData(spe_LC)$phi)
lengthscale_data <- lengthscale_data[lengthscale_data$lengthscale < 1, ]

LC_hist <- ggplot(as.data.frame(lengthscale_data), aes(x = lengthscale_data)) + 
  geom_density(color = "black", fill = "#804585") + 
  xlim(c(0, 1)) + 
  xlab("Estimated Length Scale") + 
  ylab("Density") + 
  ggtitle("LC") + 
  theme_bw()

# Ovarian
lengthscale_data <- data.frame(lengthscale = 1 / rowData(spe_Ovarian)$phi)
lengthscale_data <- lengthscale_data[lengthscale_data$lengthscale < 1, ]

Ovarian_hist <- ggplot(as.data.frame(lengthscale_data), aes(x = lengthscale_data)) + 
  geom_density(color = "black", fill = "#804585") + 
  xlim(c(0, 1)) + 
  xlab("Estimated Length Scale") + 
  ylab("Density") + 
  ggtitle("Ovarian") + 
  theme_bw()

plot <- wrap_plots(
  list(HPC_hist, Breast_hist, LC_hist, Ovarian_hist),
  ncol = 2, nrow = 2,
  axis_titles = "collect"
) +
  plot_annotation(
    tag_levels = 'A',
  theme = theme(
    plot.tag = element_text(face = "bold", size = 14))
  )

ggsave(here("plots", "supplementary", "lengthscale_distributions.png"),
       plot = plot,
       width = 20, height = 8)
