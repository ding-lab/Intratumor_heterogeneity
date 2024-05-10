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
date='20211229'

#we don't include NAT here:
data_dir="/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/8.Tumor_heterogeneity/Data"


tab=read.table('SampleList_multiome.20211229.txt',sep='\t',header=F)
tab=tab[,1:2]
colnames(tab)=c('Disease','sample_file')
tab=tab[tab$Disease==disease,]

for (sample_file in tab$sample_file){

sample=gsub('(.*)_TumorOnly.20211129.rds','\\1',sample_file)

obj=readRDS(paste('out/TumorOnly.RDS/',disease,'/',sample,'_TumorOnly.20211129.rds',sep=''))
emb=Embeddings(object = obj, reduction = "wnn.umap")
emb=as.data.frame(emb)
emb$Barcode=rownames(emb)
cl_data <- Seurat::FetchData(object = obj, vars = c('seurat_clusters'))
cl_data$Barcode=rownames(cl_data)
res_data=merge(emb,cl_data)

write.table(res_data, paste('out/Meta.data/',disease,'/',sample,'.meta.data.',date,'.tsv',sep=''),
sep='\t',quote=F,row.names=F)
print(sample)

}
