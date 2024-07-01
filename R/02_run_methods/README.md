# Overview of `02_run_methods`

The folders contain the scripts to run the individual SVG methods for each dataset in the `spoon` paper. 

Before running the scripts below, please generate the RDS objects  the `/R/01_preprocessing` folder. 

The SVG methods considered are: 

- highly variable genes (HVGs)
  - Paper: [Lun et al. 2016](https://doi.org/10.12688/f1000research.9501.2)
  - scripts: `/R/02_run_methods/run_HVGs`
- Moran's I
  - Paper: [Papadakis et al. 2023](https://cran.r-project.org/web/packages/Rfast2/Rfast2.pdf)
  - scripts: `/R/02_run_methods/run_MoransI`
- SpaGFT
  - Paper: [Chang et al. 2022](https://www.biorxiv.org/content/10.1101/2022.12.10.519929v3)
  - scripts: `/R/02_run_methods/run_SpaGFT`
    - `01_prep_humanDLPFC.R` is used to prepare the human DLPFC data for use in Python
    - `02_run_SpaGFT_humanDLPFC.py` is used to run SpaGFT in Python
    - `03_finish_SpaGFT_humanDLPFC.R` is used to add the SpaGFT results to an R object
- SPARK-X
  - Paper: [Zhu et al. 2021](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02404-0)
  - scripts: `/R/02_run_methods/run_SPARKX`
- SpatialDE2
  - Paper: [Kats et al. 2021](https://doi.org/10.1101/2021.10.27.466045)
  - scripts: `/R/02_run_methods/run_SpatialDE2`
    - `01_prep_humanDLPFC.R` is used to prepare the human DLPFC data for use in Python
    - `02_run_SpatialDE2_humanDLPFC.py` is used to run SpatialDE2 in Python
    - `03_finish_SpatialDE2_humanDLPFC.R` is used to add the SpatialDE2 results to an R object
- nnSVG
  - Paper: [Weber et al. 2023](https://www.nature.com/articles/s41467-023-39748-z)
  - scripts: `/R/02_run_methods/run_nnSVG`
- weighted nnSVG (using weights from `spoon`)
  - scripts: `/R/02_run_methods/run_weighted_nnSVG`

Output RDS files are saved to the `/outputs/results` folder (not tracked on github).
