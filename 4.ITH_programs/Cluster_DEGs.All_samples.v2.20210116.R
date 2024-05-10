#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(ggplot2)
RhpcBLASctl::blas_set_num_threads(50)
library(future)
plan("multiprocess", workers = 40)
options(future.globals.maxSize = 50 * 1024^3)



###Sample ids are different for snRNA/multiome, so process them separately
###1.snRNA:
data='snRNA'
meta=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/10.snRNA/3.Markers_expression/',
'Tumor_ManualClusters_snRNA_multiome.20220102.tsv',sep=''),header=T,sep='\t')
meta=meta[meta$Data_type=='snRNA',]
cancers=unique(meta$Disease)

for (can in cancers){
meta_s=meta[meta$Disease==can,]
samples=unique(meta_s$Piece_ID)

for (sample in samples){
data=unique(meta_s$Data_type[meta_s$Piece_ID==sample])
meta_sample=meta_s[meta_s$Piece_ID==sample,]
sample=gsub(paste(can,'_',sep=''),'',sample)

rna_dir=paste("/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/8.Tumor_heterogeneity/",
"1.snRNA_TumorOnly_IndObj/snRNA/out/TumorOnly.RDS",sep="")
rna=readRDS(paste(rna_dir,"/",can,"/",sample,"_TumorOnly.20211216.rds",sep=""))
rownames(meta_sample)=meta_sample$Barcode
orig_1=as.data.frame(rna$orig.ident)
meta_sample=meta_sample[rownames(orig_1),]

all_idents=table(meta_sample$cell_group)
idents=names(all_idents[all_idents>50])

rna$Cluster_ID=meta_sample$cell_group

Idents(rna)<-'Cluster_ID'
p=DimPlot(rna,reduction='rna.umap')
pdf(paste('plots/',can,'/',can,'_',sample,'_clusters.Dimplot.pdf',sep=''),width=6,height=5)
print(p)
dev.off()


if(length(idents)>1){
DefaultAssay(rna)='SCT'
all_degs=NULL

for (ident in idents){
degs <- FindMarkers(
  object = rna,
  ident.1 = ident,
#  ident.2='',
  only.pos = FALSE,
  min.pct = 0.1,
  min.diff.pct=0,
  assay='SCT',
  logfc.threshold=0.1
)

degs$Ident=ident
degs$Gene=rownames(degs)
degs$Piece_ID=paste(can,sample,sep='_')
all_degs=rbind(all_degs,degs)
print(paste(sample,ident,sep=' '))
}

write.table(all_degs, paste("out/",can,"/degs_",can,"_",sample,".AllTumorClusters.20220116.tsv",sep=""),
sep="\t",quote=FALSE,row.names=FALSE)
}
}
}

#########################
###2.Multiome samples:###
#########################
data='Multiome'

meta=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/10.snRNA/3.Markers_expression/',
'Tumor_ManualClusters_snRNA_multiome.20220102.tsv',sep=''),header=T,sep='\t')
meta=meta[meta$Data_type==data,]
cancers=unique(meta$Disease)

for (can in cancers){
meta_s=meta[meta$Disease==can,]
samples=unique(meta_s$Sample)

for (sample in samples){
meta_sample=meta_s[meta_s$Sample==sample,]
piece_id=unique(meta_s$Piece_ID[meta_s$Sample==sample])
piece_id=gsub(paste(can,'_',sep=''),'',piece_id)


mult_dir=paste("/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/8.Tumor_heterogeneity/",
"2.snATAC_TumorOnly_IndObj/Multiome/out/TumorOnly.RDS",sep="")
rna=readRDS(paste(mult_dir,"/",can,"/",sample,"_TumorOnly.20211129.rds",sep=""))
rownames(meta_sample)=meta_sample$Barcode
orig_1=as.data.frame(rna$orig.ident)
meta_sample=meta_sample[rownames(orig_1),]

all_idents=table(meta_sample$cell_group)
idents=names(all_idents[all_idents>50])

rna$Cluster_ID=meta_sample$cell_group

Idents(rna)<-'Cluster_ID'
p=DimPlot(rna,reduction='wnn.umap')
pdf(paste('plots/',can,'/',piece_id,'_clusters.Dimplot.pdf',sep=''),width=6,height=5)
print(p)
dev.off()

if(length(idents)>1){
DefaultAssay(rna)='SCT'
all_degs=NULL

for (ident in idents){
degs <- FindMarkers(
  object = rna,
  ident.1 = ident,
#  ident.2='',
  only.pos = FALSE,
  min.pct = 0.1,
  min.diff.pct=0,
  assay='SCT',
  logfc.threshold=0.1
)

degs$Ident=ident
degs$Gene=rownames(degs)
degs$Piece_ID=paste(can,piece_id,sep='_')
all_degs=rbind(all_degs,degs)
print(paste(piece_id,ident,sep=' '))
}

write.table(all_degs, paste("out/",can,"/degs_",can,"_",piece_id,".AllTumorClusters.20220116.tsv",sep=""),
sep="\t",quote=FALSE,row.names=FALSE)
}
}
}

