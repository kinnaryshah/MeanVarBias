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
  labs(x = "Sample Mean", y="Sample Variance", 
       title = "Using raw counts")

# voom plot
sim <- logNormCounts(sim, pseudo.count = 1)
emp_var_log <- rowVars(logcounts(sim)) # sample variance of log-counts
emp_mean_log <- rowMeans(logcounts(sim)) # sample mean of log-counts

p2 <- data.frame(x=emp_mean_log,
                 y=emp_var_log) %>%
  ggplot() +
  geom_point(aes(x=x, y=y)) +
  theme_bw() +
  labs(x = "Sample Mean", y="Sample Variance", 
       title = "Using log2(counts + 1)")

fig1 <- wrap_plots(A = p1, B = p2, nrow=1, axis_titles = "collect") + 
          plot_annotation(tag_levels = 'A') 
ggsave(filename=here("plots", "supplementary", "scales.png"), 
       plot = fig1, 
       width = 12, 
       height = 5,
       units = "in")
