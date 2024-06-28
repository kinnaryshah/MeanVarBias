# Overview of `02_run_methods`

The folders contain the scripts to run the individual SVG methods for each dataset in the `spoon` paper. 

Before running the scripts below, please generate the RDS objects  the `/R/01_preprocessing` folder. 

The SVG methods considered are: 

- highly variable genes (HVGs)
  - Paper: 
  - scripts: 
- Moran's I
  - Paper: 
  - scripts: 
- SpaGFT
  - Paper: 
  - scripts: 
- Spark-X
  - Paper: 
  - scripts: 
- SpatialDE
  - Paper: 
  - scripts: 
- nnSVG
  - Paper: [Weber et al. 2023](https://www.nature.com/articles/s41467-023-39748-z)
  - scripts: `/R/02_run_methods/run_nnSVG`
- weighted nnSVG (using weights from `spoon`)
  - scripts: `/R/02_run_methods/run_weighted_nnSVG`

Output RDS files are saved to the `/outputs/results` folder (not tracked on github).
