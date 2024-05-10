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

obj=readRDS(paste('out/TumorOnly.RDS/',disease,'/',piece_id_rna,'_TumorOnly.20211216.rds',sep=''))
emb=Embeddings(object = obj, reduction = "rna.umap")
emb=as.data.frame(emb)
emb$Barcode=rownames(emb)
cl_data <- Seurat::FetchData(object = obj, vars = c('seurat_clusters'))
cl_data$Barcode=rownames(cl_data)
res_data=merge(emb,cl_data)

write.table(res_data, paste('out/Meta.data/',disease,'/',piece_id_rna,'.meta.data.',date,'.tsv',sep=''),
sep='\t',quote=F,row.names=F)
print(piece_id_rna)
}
