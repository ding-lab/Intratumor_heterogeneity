# Analysis of Intra-Tumor Heterogeneity programs (also, called as MPs, Meta Programs).

## Analysis steps:

1. ```Cluster_DEGs.All_samples.v2.20210116.R``` -- Calculates DEGs for each cluster vs other clusters for each sample separately. It is required to have >1 cluster for this and >50 cells overall.

2. ```Calculate_DAMs_perCluster_perSample.20220130.R``` -- Calculates DAMs for each cluster vs other clusters for each sample separately. It is required to have >1 cluster for this and >50 cells overall.