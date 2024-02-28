library(ggplot2)
library(SpatialExperiment)
library(patchwork)

spe_HPC <- readRDS("HPC_data/HPC_V10B01_085_A1/spe_nnSVG.rds")
spe_ovarian <- readRDS("ovarian_data/spe_nnSVG.rds")

# use ggplot to plot densities of the marginal distributions of the logcounts
# for the HPC and ovarian datasets on the same plot
df <- data.frame(
  means <- c(rowData(spe_HPC)$mean, rowData(spe_ovarian)$mean),
  dataset <- c(rep("HPC", dim(spe_HPC)[1]), rep("ovarian", dim(spe_ovarian)[1]))
)

p1 <- ggplot(df, aes(x = means, fill = dataset)) +
  geom_density(alpha = 0.3) +
  theme_minimal() +
  labs(
    title = "Marginal distributions of logcounts",
    x = "mean of logcounts",
    y = "density"
  )

# create a new plot that zooms in to the bottom left corner of the previous plot
p2 <- ggplot(df, aes(x = means, fill = dataset)) +
  geom_density(alpha = 0.3) +
  theme_minimal() +
  labs(
    title = "Marginal distributions of logcounts",
    x = "mean of logcounts",
    y = "density"
  ) +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 4)) +
  #add caption
  annotate("text", x = 3.5, y = 3.5, label = "Zoomed in view of the bottom left corner")

pdf("vis_marginal_distributions.pdf", width = 21, height = 20)
p1
p2
dev.off()

