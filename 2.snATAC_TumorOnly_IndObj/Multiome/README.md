# Subsetting of tumor cells from multiome objects and further re-clustering, dim. reduction.

## Analyses steps:

1. ```Subset_TumorOnly_reCluster.20211228.R``` -- Make tumor only snATAC objects.

   * Extracts for each cohort set of snATAC samples ids that need to be processed. Individual objects paths are taken from the catalog.

   * Subsetting of tumor cells from snATAC objects (cases with paired snRNA/snATAC) and further re-clustering, dim. reduction. 


2. ```Get_embeddings.20211229.R``` -- Extracts UMAP coordinates from the tumor only objects and saves them. This step is needed to make work with objects easier, so that DimPlots can be made quicker with no need for loading RDS objects each time. Those tables will be used for manual cluster re-annotation.
