#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp

library(Signac)
library(Seurat)
library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 100 * 1024 ^ 2)

args = commandArgs(trailingOnly=TRUE)
disease=args[1]
date='20211129'

con<-file(paste('logs/Subset_TumorOnly_reCluster/',disease,'_SubsetTumor_Recluster.log.',date,'.txt',sep=''))
sink(con,append=TRUE)
sink(con,type='message',append=TRUE)

data_dir="/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/8.Tumor_heterogeneity/Data"

#we don't include NAT here:
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
samples_atac=paste(tab_atac$Disease.Type,tab_atac$Sample.ID,sep='_')
piece_ids_atac=paste(tab_atac$Disease.Type,tab_atac$Piece_ID,sep='_')

tab_multiome=tab[tab$Data.Type=='10x_SC_Multi_ATAC_SEQ',]
samples_multiome=paste(tab_multiome$Disease.Type,tab_multiome$Sample.ID,sep='_')
piece_ids_multiome=paste(tab_multiome$Disease.Type,tab_multiome$Piece_ID,sep='_')


anno=read.table(paste(data_dir,'/All_149_samples_metadata_data_freeze_v4.1.tsv',sep=''),sep='\t',header=TRUE)
anno_c=anno[anno$Cancer==disease,]

samples=samples_atac

for (sample in samples){
obj=readRDS(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/0.Individual_objects/',
'snATAC/out/',disease,'/',sample,'/',sample,'_processed_atac.rds',sep=''))
anno_c_sample=anno_c[anno_c$Sample==sample,]

rownames(anno_c_sample)=gsub('(.*)_(.*)_(.*)','\\3',anno_c_sample$V1)
obj_barc=as.data.frame(obj$orig.ident)
anno_c_sample=anno_c_sample[rownames(obj_barc),]

obj$cell_type=anno_c_sample$cell_type.harmonized.cancer
n_cells=nrow(anno_c_sample[anno_c_sample$cell_type.harmonized.cancer=='Tumor' & 
!is.na(anno_c_sample$cell_type.harmonized.cancer),])

if (n_cells>50){
obj_s=subset(obj,cell_type %in% c('Tumor') )

###Now re-cluster:
obj_s <- FindTopFeatures(obj_s, min.cutoff = 'q0')
obj_s <- RunTFIDF(obj_s)
obj_s <- RunSVD(obj_s)

obj_s <- RunUMAP(object = obj_s, reduction = 'lsi', dims = 2:30,reduction.name = "atac.umap")
obj_s <- FindNeighbors(object = obj_s, reduction = 'lsi', dims = 2:30)


###this helps:
options(future.globals.maxSize= 891289600)
obj_s <- FindClusters(object = obj_s, verbose = FALSE, algorithm = 3)

plot1=DimPlot(obj_s, group.by='seurat_clusters',label=T)
plot2=DimPlot(obj_s, group.by='cell_type',label=T)
p=plot1+plot2

pdf(paste('out/TumorOnly.RDS.dimplots/',disease,'/',sample,'_TumorOnly.20211129.pdf',sep=''),width=12,height=5)
print(p)
dev.off()

saveRDS(obj_s,paste('out/TumorOnly.RDS/',disease,'/',sample,'_TumorOnly.20211129.rds',sep=''))
print(sample)
}
}
