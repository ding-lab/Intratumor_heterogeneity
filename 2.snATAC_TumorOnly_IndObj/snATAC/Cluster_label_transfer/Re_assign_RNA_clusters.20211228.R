library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)

library(plyr)
library(dplyr)
library(tibble)
library(reshape)
library(reshape2)
library(RColorBrewer)
library(EnsDb.Hsapiens.v86)


date='20211228'

con<-file(paste('logs/Clusters_reassigned/All_canceres_LabelTransfer_TumorOnly.log.',date,'.txt',sep=''))
sink(con,append=TRUE)
sink(con,type='message',append=TRUE)

cancers=c('MM','GBM','ccRCC','OV','BRCA','CRC','UCEC','CESC')

annot_atac=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/DATA/',
'All_149_samples_metadata_data_freeze_v4.1.tsv',sep=''),sep='\t',header=T)

tab=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/8.Tumor_heterogeneity/',
'1.snRNA_TumorOnly_IndObj/snRNA/snRNA_Sample_IDs.20211219.tsv',sep=''),header=T,sep='\t')

all_tab=NULL
all_res_sel=NULL
for (disease in cancers){

tab_s=tab[tab$Disease==disease,]
piece_ids_rna=tab_s$piece_id_rna

for (piece_id_rna in piece_ids_rna){
print(piece_id_rna)
rna_path=paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/8.Tumor_heterogeneity/',
'1.snRNA_TumorOnly_IndObj/snRNA/out/TumorOnly.RDS/',disease,'/',piece_id_rna,'_TumorOnly.20211216.rds',sep='')

RNA=readRDS(rna_path)
sample_atac=tab_s$sample_atac[tab_s$piece_id_rna==piece_id_rna]

###Check that we have data for that Piece:
sample_atac_2=paste(disease,sample_atac,sep="_")
n_cells=nrow(annot_atac[annot_atac$Sample==sample_atac_2 & annot_atac$cell_type.harmonized.cancer=='Tumor',])
piece_id_atac=unique(annot_atac$Piece_ID[annot_atac$Sample==sample_atac_2])
piece_id_atac_2=paste(disease,piece_id_atac,sep='_')


if (n_cells > 50){
ATAC=readRDS(paste("../out/TumorOnly.RDS/",disease,"/",sample_atac_2,"_TumorOnly.20211129.rds",sep=""))
clusters=read.table(paste('TumorOnly_RNAClusters_mapped/metadata/',disease,'/',sample_atac,
'_TumorOnly.snATAC.meta.data',sep=''),sep='\t',header=T)


##############
x=as.data.frame(table(clusters %>% dplyr::select ('seurat_clusters','Cluster_ID')))
x1=aggregate(x$Freq, by=list(x$seurat_clusters),FUN='sum')
colnames(x1)=c('seurat_clusters','Sum')

res=merge(x,x1,all.x=T)
res$Fraction=res$Freq/res$Sum
res_sel=res[res$Fraction>0.5,]

orig_1=as.data.frame(ATAC$seurat_clusters)
colnames(orig_1)='seurat_clusters'
orig_1$Barcode=rownames(orig_1)
meta=orig_1
meta=merge(meta,res_sel,all.x=T)
rownames(meta)=meta$Barcode
meta=meta[rownames(orig_1),]
ATAC$RNA_Cluster_ID=meta$Cluster_ID

orig_1=as.data.frame(ATAC$RNA_Cluster_ID)
colnames(orig_1)='cell_group'
orig_1$Barcode=rownames(orig_1)
orig_1$Piece_ID_ATAC=piece_id_atac_2
orig_1$Sample_ATAC=sample_atac_2
orig_1$Disease=disease

all_tab=rbind(all_tab,orig_1)



#####Labelling cell-types for the snRNA

meta_sample=read.table(paste('DATA_snRNA_Manual_clusters/out.meta/',disease,'/',piece_id_rna,
'.meta_data.cluster_manual.20211219.tsv',sep=''),sep='\t',header=T)
rownames(meta_sample)=meta_sample$Barcode

orig_1=as.data.frame(RNA$orig.ident)
meta_sample=meta_sample[rownames(orig_1),]

RNA$RNA_Cluster_ID=meta_sample$RNA_cluster_manual

plot1 <- DimPlot(
  object = RNA,
  group.by = 'RNA_Cluster_ID',
  label = TRUE,
  repel = TRUE,
  reduction='rna.umap') + ggtitle(paste(piece_id_atac_2,' snRNA-seq',sep=""))

plot2 <- DimPlot(
  object = ATAC,
  group.by = 'RNA_Cluster_ID',
  label = TRUE,
  repel = TRUE,
  reduction='atac.umap') + ggtitle(paste(piece_id_atac_2,' snATAC-seq',sep=""))

pdf(paste("TumorOnly.RDS.dimplots.RNAClustersReassigned/plots/",disease,"/",sample_atac,"_",piece_id_atac_2,
"_snRNA_ClusterIDs_mapped_to_snATAC.pdf",sep=""),height=6,width=16)
p=CombinePlots(list(plot1,plot2), ncol = 2)
print(p)
dev.off()

res_sel$Piece_ID_ATAC=piece_id_atac_2
res_sel$Sample_ATAC=sample_atac_2
res_sel$Disease=disease

all_res_sel=rbind(all_res_sel,res_sel)
print(paste(piece_id_atac_2,disease,sep=' '))
}
}
}

write.table(all_tab,paste("TumorOnly.RDS.dimplots.RNAClustersReassigned/metadata/",
"RNA_clusters_TumorOnly.snATAC.tsv",sep=""),sep='\t',quote=FALSE)
write.table(all_res_sel,paste("TumorOnly.RDS.dimplots.RNAClustersReassigned/metadata/RNA_clusters_Fractions_",
"inTumorOnly_snATAC_Clusters.tsv",sep=""),sep='\t',quote=FALSE)
