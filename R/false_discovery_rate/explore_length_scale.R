library(SpatialExperiment)
library(spatialLIBD)
library(ggplot2)

spe_unweighted <- readRDS("spe_nnSVG_5.rds")

my_list <- 1/rowData(spe_unweighted)$phi
my_list <- my_list[my_list <= 1]


# make ggplot histogram of my list and include in the caption the number of elements in the list
p = ggplot(data = data.frame(x = my_list), aes(x = x)) +
  geom_histogram(bins = 30, fill = "blue", color = "black") +
  labs(title = paste("Histogram of 1/rowData(spe_unweighted)$phi with lengthscale","10"),
       caption = paste("Number of elements in the list:", length(my_list), "out of", length(rowData(spe_unweighted)$phi), "total elements")) +

  theme_bw()

# create some spot plots to visualize the differences 
# sample 10 genes that are true positives
# all true positives are getting captured by both simulations
set.seed(1)
tp_unweighted <- which(rowData(spe_unweighted)$ground_truth_sigma.sq != 0 & rowData(spe_unweighted)$padj <= 0.05)
ixs <- sample(tp_unweighted, 10)

tn_unweighted <- rowData(spe_unweighted)$ground_truth_sigma.sq == 0 & rowData(spe_unweighted)$padj > 0.05
ixs_tn <- sample(which(tn_unweighted), 10)

pdf("explore_length_scale.pdf")
print(p)

for (ix in ixs) {
  true_sigma.sq <- rowData(spe_unweighted)$ground_truth_sigma.sq[ix]
  est_sigma.sq <- rowData(spe_unweighted)$sigma.sq[ix]
  est_length_scale <- 1/rowData(spe_unweighted)$phi[ix]
  df <- as.data.frame(cbind(spatialCoords(spe_unweighted), expr = logcounts(spe_unweighted)[ix, ]))
  print(ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                       color = expr)) + 
          geom_point(size = 2) + 
          coord_fixed() + 
          scale_y_reverse() + 
          scale_color_viridis_c(name = "logcounts") + 
          labs(title=paste("unw true positive, true sigma.sq=", true_sigma.sq),
               caption=paste("est sigma.sq=", est_sigma.sq, ", est length scale=", est_length_scale)) + 
          theme_bw() + 
          theme(plot.title = element_text(face = "italic"), 
                panel.grid = element_blank(), 
                axis.title = element_blank(), 
                axis.text = element_blank(), 
                axis.ticks = element_blank())
  )
}

for (ix in ixs_tn) {
  true_sigma.sq <- rowData(spe_unweighted)$ground_truth_sigma.sq[ix]
  est_sigma.sq <- rowData(spe_unweighted)$sigma.sq[ix]
  est_length_scale <- 1/rowData(spe_unweighted)$phi[ix]
  df <- as.data.frame(cbind(spatialCoords(spe_unweighted), expr = logcounts(spe_unweighted)[ix, ]))
  print(ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                       color = expr)) + 
          geom_point(size = 2) + 
          coord_fixed() + 
          scale_y_reverse() + 
          scale_color_viridis_c(name = "logcounts") + 
          labs(title=paste("unw true negative, true sigma.sq=", true_sigma.sq),
               caption=paste("est sigma.sq=", est_sigma.sq, ", est length scale=", est_length_scale)) + 
          theme_bw() + 
          theme(plot.title = element_text(face = "italic"), 
                panel.grid = element_blank(), 
                axis.title = element_blank(), 
                axis.text = element_blank(), 
                axis.ticks = element_blank())
  )
}

dev.off()


