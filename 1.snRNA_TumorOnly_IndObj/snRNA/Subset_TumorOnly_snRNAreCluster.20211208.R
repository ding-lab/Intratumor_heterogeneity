#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp

library(Signac)
library(Seurat)
library(tidyverse)
library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 100 * 1024 ^ 3)
RhpcBLASctl::blas_set_num_threads(30)

args = commandArgs(trailingOnly=TRUE)
disease=args[1]
date='20211216'



con<-file(paste('logs/Subset_TumorOnly_reCluster/',disease,'_SubsetTumor_Recluster.log.',date,'.txt',sep=''))
sink(con,append=TRUE)
sink(con,type='message',append=TRUE)

#we don't include NAT here:
data_dir="/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/8.Tumor_heterogeneity/Data"

tab=read.table(paste(data_dir,'/ATAC_samples.FinalList.v2.20211020.txt',sep=''),sep='\t',header=TRUE)
tab=tab[tab$Sample.Type!='Normal',]

if (disease=='ccRCC'){
tab=tab[tab$Include.in.the.downstream.analysis=='TRUE' &
!is.na(tab$Include.in.the.downstream.analysis=='TRUE') & tab$Disease.Type %in% c(disease,'PKD'),]
}else{
tab=tab[tab$Include.in.the.downstream.analysis=='TRUE' &
!is.na(tab$Include.in.the.downstream.analysis=='TRUE') & tab$Disease.Type %in% c(disease),]
}

tab_atac=tab[tab$Data.Type %in% c('snATAC','scATAC'),]
samples_atac=paste(tab_atac$Sample.ID,sep='_')
piece_ids_atac=paste(tab_atac$Disease.Type,tab_atac$Piece_ID,sep='_')

tab_multiome=tab[tab$Data.Type=='10x_SC_Multi_ATAC_SEQ',]
samples_multiome=paste(tab_multiome$Sample.ID,sep='_')
piece_ids_multiome=paste(tab_multiome$Disease.Type,tab_multiome$Piece_ID,sep='_')


paths_snRNA=read_delim(paste(data_dir,'/snRNA_samples_included_inDownstream.v2.20211209.txt',sep=''),
delim='\t')
paths_snRNA=as.data.frame(paths_snRNA)
colnames(paths_snRNA)=gsub(' ','_',colnames(paths_snRNA))
colnames(paths_snRNA)=gsub('_\\(127_total\\)','',colnames(paths_snRNA))
paths_snRNA=paths_snRNA[paths_snRNA$Include_in_snRNA_analysis==TRUE,]


anno_c=read.table(paste(data_dir,'/cell_type_snRNA_merged/',disease,'/',disease,
'_cell_type.harmonized.cancer.meta.data',sep=''),sep='\t',header=TRUE)
anno_c=anno_c[anno_c$Sample_ATAC %in% samples_atac,]

piece_ids_rna=unique(anno_c$Piece_ID)
#samples_rna=unique(anno_c$Sample_RNA)

for (piece_id_rna in piece_ids_rna){

path_obj=paths_snRNA$Processed_Object_RDS[paths_snRNA$Piece_ID==piece_id_rna]
obj=readRDS(path_obj)
anno_c_sample=anno_c[anno_c$Piece_ID==piece_id_rna,]

rownames(anno_c_sample)=gsub('(.*)_(.*)_(.*)','\\3',anno_c_sample[,1])
obj_barc=as.data.frame(obj$orig.ident)
anno_c_sample=anno_c_sample[rownames(obj_barc),]

obj$cell_type=anno_c_sample$cell_type.harmonized.cancer
n_cells=nrow(anno_c_sample[anno_c_sample$cell_type.harmonized.cancer=='Tumor' & 
!is.na(anno_c_sample$cell_type.harmonized.cancer),])

if (n_cells>50){
obj_s=subset(obj,cell_type %in% c('Tumor') )

###Now re-cluster:
DefaultAssay(obj_s) <- "RNA"
obj_s <- SCTransform(obj_s, vars.to.regress = c("nCount_RNA","percent.mito"),return.only.var.genes = F)
obj_s <- RunPCA(obj_s, npcs=30,verbose = FALSE)
obj_s = RunUMAP(
obj_s, 
reduction = 'pca',
dims = 1:30,
reduction.name = 'rna.umap',
reduction.key = 'rnaUMAP_'
)
obj_s <- FindNeighbors(obj_s, reduction = "pca", dims = 1:30)
###this helps:
#options(future.globals.maxSize= 891289600)
obj_s <- FindClusters(obj_s)

umap_data <- Seurat::FetchData(object = obj_s, vars = c("UMAP_1","UMAP_2",'seurat_clusters'))
umap_data$Barcode=rownames(umap_data)

#write.table(umap_data, paste('out/Meta.data/',disease,'/',piece_id_rna,'.meta.data.',date,'.tsv',sep=''),
#sep='\t',quote=F,row.names=F)
plot1=DimPlot(obj_s, group.by='seurat_clusters',label=T,reduction='rna.umap')
plot2=DimPlot(obj_s, group.by='cell_type',label=T,reduction='rna.umap')
p=plot1+plot2

pdf(paste('out/TumorOnly.RDS.dimplots/',disease,'/',piece_id_rna,'_TumorOnly.',date,'.pdf',sep=''),width=12,
height=5)
print(p)
dev.off()

saveRDS(obj_s,paste('out/TumorOnly.RDS/',disease,'/',piece_id_rna,'_TumorOnly.',date,'.rds',sep=''))
print(sample)
}
}
