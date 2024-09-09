###################################
# dataset: Visium human hippocampus
###################################
library(here)
library(SpatialExperiment)
library(nnSVG)
library(scran)
# BiocManager::install("DropletUtils") if not installed

# ---------
# load data
# ---------

spe <- read10xVisium(here("outputs/raw/humanHPC_V12D07-335_D1/outs/"),
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

spe <- filter_genes(spe)

spe <- logNormCounts(spe)


# -----------
# save object
# -----------

fn <- here("outputs", "preprocessed", "spe_humanHPC_V12D07-335_D1_preprocessed.rds")
saveRDS(spe, file = fn)

# -----------
# session information
# -----------

sessionInfo()

# R version 4.4.1 Patched (2024-08-21 r87049)
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
#   [1] scran_1.33.2                scuttle_1.15.4             
# [3] nnSVG_1.9.0                 SpatialExperiment_1.15.1   
# [5] SingleCellExperiment_1.27.2 SummarizedExperiment_1.35.1
# [7] Biobase_2.65.1              GenomicRanges_1.57.1       
# [9] GenomeInfoDb_1.41.1         IRanges_2.39.2             
# [11] S4Vectors_0.43.2            BiocGenerics_0.51.1        
# [13] MatrixGenerics_1.17.0       matrixStats_1.4.1          
# [15] here_1.0.1                 
# 
# loaded via a namespace (and not attached):
#   [1] rjson_0.2.22              rhdf5_2.49.0             
# [3] lattice_0.22-6            rhdf5filters_1.17.0      
# [5] tools_4.4.1               parallel_4.4.1           
# [7] cluster_2.1.6             R.oo_1.26.0              
# [9] pkgconfig_2.0.3           BiocNeighbors_1.99.0     
# [11] Matrix_1.7-0              sparseMatrixStats_1.17.2 
# [13] dqrng_0.4.1               lifecycle_1.0.4          
# [15] GenomeInfoDbData_1.2.12   compiler_4.4.1           
# [17] statmod_1.5.0             bluster_1.15.1           
# [19] codetools_0.2-20          crayon_1.5.3             
# [21] R.utils_2.12.3            BiocParallel_1.39.0      
# [23] DelayedArray_0.31.11      limma_3.61.9             
# [25] BRISC_1.0.5               magick_2.8.4             
# [27] abind_1.4-5               rsvd_1.0.5               
# [29] locfit_1.5-9.10           metapod_1.13.0           
# [31] BiocSingular_1.21.2       rprojroot_2.0.4          
# [33] grid_4.4.1                cli_3.6.3                
# [35] SparseArray_1.5.31        magrittr_2.0.3           
# [37] S4Arrays_1.5.5            edgeR_4.3.14             
# [39] DelayedMatrixStats_1.27.3 UCSC.utils_1.1.0         
# [41] XVector_0.45.0            httr_1.4.7               
# [43] DropletUtils_1.25.2       igraph_2.0.3             
# [45] RANN_2.6.2                R.methodsS3_1.8.2        
# [47] ScaledMatrix_1.13.0       beachmat_2.21.6          
# [49] HDF5Array_1.33.5          pbapply_1.7-2            
# [51] irlba_2.3.5.1             rlang_1.1.4              
# [53] Rcpp_1.0.13               rdist_0.0.5              
# [55] jsonlite_1.8.8            Rhdf5lib_1.27.0          
# [57] R6_2.5.1                  zlibbioc_1.51.1   