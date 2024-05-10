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

args = commandArgs(trailingOnly=TRUE)
disease=args[1]
date='20211228'

con<-file(paste('logs/LabelTransfer_TumorOnly/',disease,'_LabelTransfer_TumorOnly.log.',date,'.txt',sep=''))
sink(con,append=TRUE)
sink(con,type='message',append=TRUE)

tab=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/8.Tumor_heterogeneity/',
'1.snRNA_TumorOnly_IndObj/snRNA/snRNA_Sample_IDs.20211219.tsv',sep=''),header=T,sep='\t')
tab_s=tab[tab$Disease==disease,]
piece_ids_rna=tab_s$piece_id_rna

for (piece_id_rna in piece_ids_rna){
print(piece_id_rna)
rna_path=paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/8.Tumor_heterogeneity/',
'1.snRNA_TumorOnly_IndObj/snRNA/out/TumorOnly.RDS/',disease,'/',piece_id_rna,'_TumorOnly.20211216.rds',sep='')

RNA=readRDS(rna_path)
sample_atac=tab_s$sample_atac[tab_s$piece_id_rna==piece_id_rna]

#Check that we have data for that Piece id:
annot=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/DATA/',
'All_149_samples_metadata_data_freeze_v4.1.tsv',sep=''),sep='\t',header=T)
sample_atac_2=paste(disease,sample_atac,sep="_")
n_cells=nrow(annot[annot$Sample==sample_atac_2 & annot$cell_type.harmonized.cancer=='Tumor',])

if (n_cells > 50){
ATAC=readRDS(paste("/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/8.Tumor_heterogeneity/2.",
"snATAC_TumorOnly_IndObj/snATAC/out/TumorOnly.RDS/",disease,"/",disease,"_",sample_atac,
"_TumorOnly.20211129.rds",sep=""))

DefaultAssay(ATAC)='X500peaksMACS2'

#Labelling cell-types for the snRNA

meta=read.table(paste('DATA_snRNA_Manual_clusters/out.meta/',disease,'/',piece_id_rna,
'.meta_data.cluster_manual.20211219.tsv',sep=''),sep='\t',header=T)
rownames(meta)=meta$Barcode
orig_1=as.data.frame(RNA$orig.ident)
meta=meta[rownames(orig_1),]

RNA$Cluster_ID=meta$RNA_cluster_manual


DefaultAssay(ATAC) <- 'ATACGeneActivity'

ATAC <- NormalizeData(
  object = ATAC,
  assay = 'ATACGeneActivity',
  normalization.method = 'LogNormalize',
  scale.factor = median(ATAC$nCount_ATACGeneActivity)
)

DefaultAssay(RNA)='RNA'
RNA <- NormalizeData(RNA, normalization.method = "LogNormalize", scale.factor = 10000)
RNA <- FindVariableFeatures(RNA, selection.method = "vst", nfeatures = 3000)

transfer.anchors <- FindTransferAnchors(
  reference = RNA,
  query = ATAC,
  reduction = 'cca'
#  features=intersect(rownames(RNA),rownames(ATAC))
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = RNA$Cluster_ID,
  weight.reduction = ATAC[['lsi']],
  dims=2:30
)

ATAC <- AddMetaData(object = ATAC, metadata = predicted.labels)
ATAC$Cluster_ID=ATAC$predicted.id
ATAC$Piece_ID_rna=piece_id_rna

#Making plots:

plot1 <- DimPlot(
  object = RNA,
  group.by = 'Cluster_ID',
  label = TRUE,
  reduction='rna.umap',
  repel = TRUE) + ggtitle(paste(piece_id_rna,' snRNA-seq',sep=""))

plot2 <- DimPlot(
  object = ATAC,
  group.by = 'Cluster_ID',
  label = TRUE,
  repel = TRUE) + ggtitle(paste(piece_id_rna,' snATAC-seq',sep=""))

pdf(paste("TumorOnly_RNAClusters_mapped/plots/",disease,"/",sample_atac,"_",piece_id_rna,
"_snRNA_ClusterIDs_mapped_to_snATAC.pdf",sep=""),height=6,width=16)
p=CombinePlots(list(plot1,plot2), ncol = 2)
print(p)
dev.off()

atac_meta=ATAC@meta.data
write.table(atac_meta,paste("TumorOnly_RNAClusters_mapped/metadata/",disease,"/",sample_atac,
"_TumorOnly.snATAC.meta.data",sep=""),sep='\t',quote=FALSE)
}
}

