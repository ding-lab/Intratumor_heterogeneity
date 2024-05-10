# Subsetting of tumor cells from snRNA objects (cases with paired snRNA/snATAC) and further re-clustering, dim. reduction.

All scripts are in the ```snRNA``` folder.


## Analysis steps:


1. ```Identify_samples.20211219.R``` -- Identify and make a table of cases, that have both snRNA and snATAC samples. The table should have ids that can be used to extract the paths for the objects from the catalog.


2. ```Subset_TumorOnly_snRNAreCluster.20211208.R``` -- Subset tumor only cells, and do further re-clustering and dim reduction.

   * Extracts for each cohort set of snRNA samples ids that need to be processed. Individual objects paths are taken from the catalog.

   * Susets tumor cells with further clustering, dim. reduction (N tumor cells should be > 50). Saved DimPlots and RDS objects.


3. ```Get_embeddings.20211219.R``` -- Extracts UMAP coordinates from the tumor only objects and saves them. This needed to make work with objects easier, so that DimPlots can be made quicker with no need for loading RDS objects each time. Those tables will be used for manual cluster re-annotation.