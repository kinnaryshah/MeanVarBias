scripts:

-   compare_factors.R: create grid of plots (ridge plots, avg FDR, avg TPR, avg TNR) to compare different simulations

-   combined_FDR_vs_alpha.R: create p-value distribution, avg FDR, avg TPR, avg TNR plots

-   explore_sims\_\*\_percent.Rmd: explore distributions of means/weights and spot plots for chosen simulation to investigate how to improve the simulations

folders: I have not added most spe objects because they are quite large for some of the simulations.

1.  first term in file name

    1.  Folders that start with "reps" are run with arrays to create multiple iterations of the same simulation parameters with different random seeds - this allows us to average over multiple iterations. The names of results files correspond to the number of the iteration, so p_value_distribution_1.pdf indicates the first iteration.

    2.  Folders that start with "single" are run with arrays to vary the percentages of SVGs with otherwise the same parameters. The names of results files correspond to the fraction of nonSVGs, so p_value_distribution_3.pdf indicates the 30% nonSVGs.

2.  number of spots

3.  length scale parameter

4.  is it a specific percentage of nonSVGs or varying? if "varying", this means that there are subfolders that contain a few different percentages of nonSVGs with the parameters outlined in the folder name

5.  number of genes

6.  range of $\sigma^2$

7.  range of $\beta$
