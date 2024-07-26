###################################
# dataset: Visium human breast cancer 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9044823/
###################################
library(scran)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(nnSVG)
# BiocManager::install("DropletUtils") if not installed


# ---------
# load data
# ---------

spe <-  read10xVisium(here("outputs", "raw", "humanSubtypeBreast", "outs"),
                      type = "sparse",   # use sparse (not HDF5) format
                      data = "raw",     
                      images = "lowres", # specify which image(s) to include
                      load = TRUE)      # specify whether or not to load image(s)

# -------------
# preprocessing
# -------------

#row data formatting
rowData(spe)$gene_id <- rownames(spe)
rowData(spe)$gene_name <- rowData(spe)$symbol
rowData(spe)$symbol <- NULL

# get meta data to remove "artefact" spots
meta <- read.csv(here("outputs", "raw", "humanSubtypeBreast", "metadata", "CID4290_metadata.csv"),
                 row.names = 1)
meta <- meta %>%
  filter(Classification != "Artefact")

intersect.spots <- intersect(colnames(spe), rownames(meta))
spe <- spe[, intersect.spots]

## remove genes without enough data
spe <- filter_genes(spe, 2, 0.2)
dim(spe)
# [1] 12325  2419

# create logcounts matrix
spe <- logNormCounts(spe)

# -----------
# save object
# -----------

fn <- here("outputs", "preprocessed", "spe_humanSubtypeBreast_preprocessed.rds")
saveRDS(spe, file = fn)

# -----------
# session information
# -----------

sessionInfo()

# R version 4.4.1 Patched (2024-07-24 r86925)
# Platform: x86_64-conda-linux-gnu
# Running under: Rocky Linux 9.4 (Blue Onyx)
# 
# Matrix products: default
# BLAS:   /jhpce/shared/community/core/conda_R/4.4.x/R/lib64/R/lib/libRblas.so 
# LAPACK: /jhpce/shared/community/core/conda_R/4.4.x/R/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: US/Eastern
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices datasets  utils     methods  
# [8] base     
# 
# other attached packages:
#   [1] nnSVG_1.9.0                 lubridate_1.9.3            
# [3] forcats_1.0.0               stringr_1.5.1              
# [5] dplyr_1.1.4                 purrr_1.0.2                
# [7] readr_2.1.5                 tidyr_1.3.1                
# [9] tibble_3.2.1                ggplot2_3.5.1              
# [11] tidyverse_2.0.0             here_1.0.1                 
# [13] SpatialExperiment_1.15.1    scran_1.33.1               
# [15] scuttle_1.15.2              SingleCellExperiment_1.27.2
# [17] SummarizedExperiment_1.35.1 Biobase_2.65.0             
# [19] GenomicRanges_1.57.1        GenomeInfoDb_1.41.1        
# [21] IRanges_2.39.2              S4Vectors_0.43.2           
# [23] BiocGenerics_0.51.0         MatrixGenerics_1.17.0      
# [25] matrixStats_1.3.0          
# 
# loaded via a namespace (and not attached):
#   [1] pbapply_1.7-2             rlang_1.1.4              
# [3] magrittr_2.0.3            compiler_4.4.1           
# [5] DelayedMatrixStats_1.27.2 vctrs_0.6.5              
# [7] rdist_0.0.5               pkgconfig_2.0.3          
# [9] crayon_1.5.3              magick_2.8.4             
# [11] XVector_0.45.0            utf8_1.2.4               
# [13] BRISC_1.0.5               tzdb_0.4.0               
# [15] UCSC.utils_1.1.0          bluster_1.15.0           
# [17] zlibbioc_1.51.1           beachmat_2.21.4          
# [19] jsonlite_1.8.8            rhdf5filters_1.17.0      
# [21] DelayedArray_0.31.9       Rhdf5lib_1.27.0          
# [23] BiocParallel_1.39.0       irlba_2.3.5.1            
# [25] parallel_4.4.1            cluster_2.1.6            
# [27] R6_2.5.1                  stringi_1.8.4            
# [29] limma_3.61.5              Rcpp_1.0.13              
# [31] R.utils_2.12.3            Matrix_1.7-0             
# [33] igraph_2.0.3              timechange_0.3.0         
# [35] tidyselect_1.2.1          abind_1.4-5              
# [37] codetools_0.2-20          lattice_0.22-6           
# [39] withr_3.0.0               pillar_1.9.0             
# [41] generics_0.1.3            rprojroot_2.0.4          
# [43] hms_1.1.3                 sparseMatrixStats_1.17.2 
# [45] munsell_0.5.1             scales_1.3.0             
# [47] glue_1.7.0                metapod_1.13.0           
# [49] tools_4.4.1               BiocNeighbors_1.23.0     
# [51] ScaledMatrix_1.13.0       locfit_1.5-9.10          
# [53] RANN_2.6.1                rhdf5_2.49.0             
# [55] grid_4.4.1                DropletUtils_1.25.2      
# [57] edgeR_4.3.5               colorspace_2.1-0         
# [59] GenomeInfoDbData_1.2.12   BiocSingular_1.21.2      
# [61] HDF5Array_1.33.5          cli_3.6.3                
# [63] rsvd_1.0.5                fansi_1.0.6              
# [65] S4Arrays_1.5.5            gtable_0.3.5             
# [67] R.methodsS3_1.8.2         SparseArray_1.5.25       
# [69] dqrng_0.4.1               rjson_0.2.21             
# [71] R.oo_1.26.0               lifecycle_1.0.4          
# [73] httr_1.4.7                statmod_1.5.0  