# Make tumor only snATAC objects and do snRNA cluster labels transfer.

## Analyses steps:

1. ```Subset_TumorOnly_reCluster.20211126.R```

   * Extracts for each cohort set of snATAC samples ids that need to be processed. Individual objects paths are taken from the catalog.

   * Subsets tumor cells with further clustering, dim. reduction (N tumor cells should be > 50). Save DimPlots and RDS objects.


2. ```Cluster_label_transfer``` -- Transfer snRNA clusters labels to tumor cells in snATAC objects.

   * Note: go to the above dir to make further analysis. This analysis requires prior manual cluster annotation of tumor only snRNA samples, and saved meta.data -- this is done on the local PC. 

   * ```TumorOnly_Integrate_withRNA.20211219.R``` -- Transfer snRNA clusters labels to tumor cells in snATAC objects. It adds predictions to snATAC objects meta.data and saves them in ```TumorOnly_RNAClusters_mapped/metadata```. RDS are not saved, because annotation in the meta.data is enough.

   * ```Re_assign_RNA_clusters.20211228.R``` -- Re-assigning of cluster ids using predicted annotation from label transfer (previous step). It goes over all seurat clusters in ATAC, and checks if there are >50% of cells annotated with any cluster id. If yes, the cluster will be annotated with the respective cluster label. 

     	+ Note: meta.data and DimPlots are saved in ```TumorOnly.RDS.dimplots.RNAClustersReassigned```.

   * ```Combine_tables.20220101.R``` -- Combine all manual cluster annotations from separate tables into one master table, that will be used in further analysis steps.