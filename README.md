# Script for tumor heterogeneity analysis.

Note: the procedure is a bit different for cases which are multiome, and the ones with paired snRNA/snATAC.

## Key differences:

   * Annotation used is based on cohort-level merged objects:

     + For multiome objects -- They have tumor cells annotated same in RNA/ATAC assays, same barcodes, so we use only Tumor cells that are annotated as Tumor in both RNA and ATAC assays.

     + For paired snRNA/snATAC -- We subtract separately tumor cells using snRNA and snATAC annotations respectively, because different barcodes in two objects.

   * Identification of manual clusters (need to have same cluster labels between snRNA/snATAC assays):

     + For multiome -- They are identified using joint wnn clustering, that is based on both snRNA and snATAC assays.

     + For separate snRNA/snATAC - We re-assign clusters manually in snRNA object, and then perform label transfer procedure of cluster labels, and assign clusters in snATAC based on snRNA cluster labels if >X % of cells were labeled of any particular cluster.



## Steps of the analysis:

1. ```1.snRNA_TumorOnly_IndObj``` -- Subsetting of tumor cells from snRNA objects (cases with paired snRNA/snATAC) and further re-clustering, dim. reduction.


2. ```2.snATAC_TumorOnly_IndObj```:

   * Subsetting of tumor cells from snATAC objects (cases with paired snRNA/snATAC) and further re-clustering, dim. reduction. Also, transfer of snRNA cluster labels to snATAC tumor cells.

   * Subsetting of tumor cells from multiome objects and further re-clustering, dim. reduction
  

3. ```3.DAM_analysis``` -- Calculate motif scores per identified tumor cluster.


4. ```4.ITH_programs``` -- Analysis of Intra-Tumor Heterogeneity programs (also, called as MPs, Meta Programs).


## These analyses were not included in the final set of analyses.