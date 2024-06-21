library(splatter)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scuttle)
library(here)

# load data: simulated counts from splatter (10k genes, 100 cells, default parameters)
fn <- here("outputs", "preprocessed", "sce_simulated_Gamma-Poisson_preprocessed.rds")
sim <- readRDS(fn)

### calculate sample mean and variance 
emp_mean <- rowMeans(counts(sim))
emp_var <- rowVars(counts(sim))

p1 <- data.frame(x=emp_mean, y=emp_var) %>%
  ggplot() +
  geom_point(aes(x=x, y=y)) +
  theme_bw() +
  labs(x = "sample mean", y="sample variance", 
       title = "Using raw counts")

# voom plot
sim <- logNormCounts(sim)
emp_var_log <- rowVars(logcounts(sim)) # sample variance of log-counts
emp_mean_log <- rowMeans(logcounts(sim)) # sample mean of log-counts

p2 <- data.frame(x=emp_mean_log,
                 y=emp_var_log) %>%
  ggplot() +
  geom_point(aes(x=x, y=y)) +
  theme_bw() +
  labs(x = "sample mean", y="standard variance", 
       title = "Using log-counts")

fig1 <- wrap_plots(p1, p2, nrow=1)
ggsave(here("plots", "supplementary", "scales.png"), fig1)
