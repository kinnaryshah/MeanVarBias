# code for spoon manuscript

This repository contains the code to reproduce the analyses and figures generated for the spoon manuscript.

Inside the R folder, there are many subfolders for the analyses and figures.

-   compare_methods: visualizing the mean-rank relationship after using various SVG detection methods

    -   DLPFC dataset: did not add spe objects due to GitHub size constraints

    -   simulated datasets are described in a README.md file inside the /simulations folder (did not do the calculations for all datasets)

-   false_discovery_rate: many iterations of simulations with varying parameters (described further in a README.md file inside the folder)

-   false_positive rate: calculating the false positive rate for two null simulations, one simulation using uniform marginal means and the other simulation using DLPFC-based marginal means

    -   uses simulated data in the /simulations/nonSVGs subfolder

-   old: contains older code for various exploratory ideas

-   simstpy: trying out the Python package from the Pinello benchmarking paper to simulate data

    -   the actual data files are in JHPCE - can be added to GitHub if we want to look more into them!

-   simulations: contains the code to create various simulations and the simulated data itself (described further in a README.md file inside the folder)

    -   either all SVGs or all nonSVGs

-   weight_effect: using nnSVG and weighted nnSVG in various real and simulated datasets and creating various visualizations

    -   spe objects for real datasets were too large to add to GitHub

    -   simulated datasets are described in a README.md file inside the /simulations folder (did not do the calculations for all datasets)
