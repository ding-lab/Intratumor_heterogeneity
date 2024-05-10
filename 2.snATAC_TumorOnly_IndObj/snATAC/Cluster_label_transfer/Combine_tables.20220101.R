library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)

library(plyr)
library(dplyr)
library(tibble)
library(reshape)
library(RColorBrewer)

library(EnsDb.Hsapiens.v86)

all_meta=NULL
cancers=c('ccRCC','BRCA','CRC','OV','CESC','OV','UCEC','MM')

tab=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/8.Tumor_heterogeneity/',
'1.snRNA_TumorOnly_IndObj/snRNA/snRNA_Sample_IDs.20211219.tsv',sep=''),header=T,sep='\t')

for (disease in cancers){
tab_s=tab[tab$Disease==disease,]
samples_atac=tab_s$sample_atac
if (disease=='ccRCC'){
samples_atac=samples_atac[1:18]
}
#piece_ids_rna=tab_s$piece_id_rna

for (sample_atac in samples_atac){
meta=read.table(paste('TumorOnly_RNAClusters_mapped/metadata/',disease,'/',sample_atac,'_TumorOnly.snATAC.meta.data',sep=''),
sep='\t',header=T)
meta= meta %>%dplyr::select('orig.ident','Cluster_ID','Piece_ID_rna')
meta$Sample_ATAC=sample_atac
meta$Disease=disease
all_meta=rbind(all_meta,meta)
}
}

all_meta=as.data.frame(all_meta)
all_meta$Barcode=rownames(all_meta)

write.table(all_meta,'Combined_cluster_IDs.20211220.tsv',sep='\t',quote=F,row.names=F)