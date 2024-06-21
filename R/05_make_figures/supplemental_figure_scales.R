library(splatter)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scuttle)
library(here)

fn <- here("outputs", "preprocessed", "sce_simulated_Gamma-Poisson_preprocessed.rds")
sim <- readRDS(fn)

### calculate empirical mean and variance 
emp_mean <- rowMeans(counts(sim))
emp_var <- rowVars(counts(sim))

p1 <- data.frame(x=emp_mean, y=emp_var) %>%
  ggplot() +
  geom_point(aes(x=x, y=y)) +
  theme_bw() +
  labs(x = "empirical mean", y="empirical variance", title = "Empirical scale")

# voom plot
sim <- logNormCounts(sim)
emp_var_log <- rowVars(logcounts(sim))

p2 <- data.frame(x=log2(emp_mean + 0.5), y=sqrt(emp_var_log)) %>%
  ggplot() +
  geom_point(aes(x=x, y=y)) +
  theme_bw() +
  labs(x = "log2(empirical mean)", y="std. dev. logcounts", title = "log2 scale")

fig1 <- wrap_plots(p1, p2, nrow=1)
ggsave(here("plots", "supplementary", "scales.png"), fig1)
