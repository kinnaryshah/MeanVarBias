-   all_SVGs_300: 300 true positives, uniform marginal means with $\beta \in (0.5,1)$

-   all_SVGs_1000: 1000 true positives, uniform marginal means with $\beta \in (0.5,1)$

-   nonSVGs: (the analyses for these are in the false_positive_rate folder)

    -   simple_1000: 1000 true negatives, marginal means sampled from DLPFC dataset

    -   simple_1000_unif_means: 1000 true negatives, uniform marginal means with $\beta \in (0.5,1)$

-   sample_means_300: 300 true positives, marginal means sampled from DLPFC dataset

-   sample_means_300_log2: 300 true positives, marginal means sampled from DLPFC dataset

    -   When sampling from the DLPFC dataset, I calculate the gene-wise means from normalized counts and then use a log transformation when sampling these means. In this simulation, I used the log2 scale for this transformation instead of the ln scale in the other simulations.

        -   beta \<- sample(log2(marginal_means+1), n_genes)

-   sample_means_300_log2_threshold_3: 300 true positives, marginal means sampled from DLPFC dataset

    -   same as log2 explanation above

    -   I used 10 as the threshold for outliers (removing gene-wise means of normalized counts below 10 before sampling) in the other simulations. In this one, I used 3 as the threshold.

-   sample_means_1000: 1000 true positives, marginal means sampled from DLPFC dataset
