library(ggplot2)
library(SpatialExperiment)
library(patchwork)
library(spatialLIBD)
library(nnSVG)
library(scuttle)

spe_sim <- readRDS("spe_simulation.rds")

#already contains logcounts
spe <- fetch_data(type = "spe")

spe <- spe[, spe$sample_id == "151507"]
spe <- spe[, !is.na(unfactor(colData(spe)$spatialLIBD))]

## Remove genes without enough data
spe <- filter_genes(spe, 2, 0.2)

spe_sim <- logNormCounts(spe_sim, transform = "none")
spe_sim <- logNormCounts(spe_sim, transform = "log")
rowData(spe_sim)$norm_mean <- rowMeans(normcounts(spe_sim))
rowData(spe_sim)$log_mean <- rowMeans(logcounts(spe_sim))

spe <- logNormCounts(spe, transform = "none")
spe <- logNormCounts(spe, transform = "log")
rowData(spe)$norm_mean <- rowMeans(normcounts(spe))
rowData(spe)$log_mean <- rowMeans(logcounts(spe))

# use ggplot to plot densities of the marginal distributions of the logcounts
# for the HPC and ovarian datasets on the same plot
df <- data.frame(
  means <- c(rowData(spe)$norm_mean, rowData(spe_sim)$norm_mean),
  dataset <- c(rep("DLPFC", dim(spe)[1]), rep("simulated", dim(spe_sim)[1]))
)

p1 <- ggplot(df, aes(x = means, fill = dataset)) +
  geom_density(alpha = 0.3) +
  theme_minimal() +
  labs(
    title = "Marginal distributions of norm counts",
    x = "mean of logcounts",
    y = "density"
  )

# create similar plot for logcounts
df <- data.frame(
  means <- c(rowData(spe)$log_mean, rowData(spe_sim)$log_mean),
  dataset <- c(rep("DLPFC", dim(spe)[1]), rep("simulated", dim(spe_sim)[1]))
)

p2 <- ggplot(df, aes(x = means, fill = dataset)) +
  geom_density(alpha = 0.3) +
  theme_minimal() +
  labs(
    title = "Marginal distributions of logcounts",
    x = "mean of logcounts",
    y = "density"
  )

# find index of gene from spe_sim with median logcounts close to the median logcounts of spe

spe_sim_median <- median(rowData(spe_sim)$log_mean)
spe_median <- median(rowData(spe)$log_mean)

index <- which.min(abs(rowData(spe_sim)$log_mean - spe_sim_median))

df <- as.data.frame(cbind(spatialCoords(spe_sim), expr = counts(spe_sim)[index, ]))
p3 <- ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = expr)) + 
  geom_point(size = 0.8) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_gradient(low = "gray90", high = "blue", 
                       trans = "sqrt", breaks = range(df$expr), 
                       name = "counts") + 
  ggtitle("example of length scale = 200") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "italic"), 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

df <- as.data.frame(cbind(spatialCoords(spe_sim), expr = logcounts(spe_sim)[index, ]))
p4 <- ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                     color = expr)) + 
  geom_point(size = 0.8) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_gradient(low = "gray90", high = "blue", 
                       trans = "sqrt", breaks = range(df$expr), 
                       name = "counts") + 
  ggtitle("example of length scale = 200") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "italic"), 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

# plot the gene expression of the gene in spe_sim and spe
p3 <- plot_gene_expression(spe_sim, gene, "simulated")

pdf("vis_marginal_distributions.pdf", width = 21, height = 20)
p1
p2
p3
p4
dev.off()

