library(SpatialExperiment)
library(here)
library(scran)

# ---------
# load data
# ---------

sim_dir <- "reps_968_50_50per_1000_0.2_to_3_0.5_to_9"
K=1
fn = here("outputs", "simulations", sim_dir, paste0("spe_simulation_",K,".rds"))
spe <- readRDS(file=fn)
spe <- logNormCounts(spe)

# ---------
# save in .csv files
# ---------

fn <- here("outputs", "preprocessed", "simulations", sim_dir, "counts_preprocessed.csv")
dense_matrix <- as.matrix(counts(spe))
df <- as.data.frame(dense_matrix)
write.csv(df, file=fn, row.names=TRUE)

fn <- here("outputs", "preprocessed", "simulations", sim_dir, "coords_preprocessed.csv")
df <- as.data.frame(spatialCoords(spe))
write.csv(df, file=fn, row.names=TRUE)

# -----------
# session information
# -----------

sessionInfo()

# R version 4.4.1 Patched (2024-06-25 r86854)
# Platform: x86_64-conda-linux-gnu
# Running under: Rocky Linux 9.2 (Blue Onyx)
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
#   [1] here_1.0.1                  SpatialExperiment_1.15.1   
# [3] SingleCellExperiment_1.27.2 SummarizedExperiment_1.35.1
# [5] Biobase_2.65.0              GenomicRanges_1.57.1       
# [7] GenomeInfoDb_1.41.1         IRanges_2.39.0             
# [9] S4Vectors_0.43.0            BiocGenerics_0.51.0        
# [11] MatrixGenerics_1.17.0       matrixStats_1.3.0          
# 
# loaded via a namespace (and not attached):
#   [1] crayon_1.5.3            httr_1.4.7              magick_2.8.3           
# [4] UCSC.utils_1.1.0        jsonlite_1.8.8          DelayedArray_0.31.3    
# [7] rjson_0.2.21            rprojroot_2.0.4         grid_4.4.1             
# [10] abind_1.4-5             compiler_4.4.1          Rcpp_1.0.12            
# [13] XVector_0.45.0          lattice_0.22-6          R6_2.5.1               
# [16] SparseArray_1.5.10      GenomeInfoDbData_1.2.12 magrittr_2.0.3         
# [19] Matrix_1.7-0            tools_4.4.1             zlibbioc_1.51.1        
# [22] S4Arrays_1.5.1         
