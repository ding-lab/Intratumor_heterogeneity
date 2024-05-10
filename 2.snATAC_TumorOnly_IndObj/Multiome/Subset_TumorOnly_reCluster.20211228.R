###Make cutoff of 200 cells, so that FindMultiModalNeighbors will be able to run: https://github.com/satijalab/seurat/issues/3646

#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp

library(Signac)
library(Seurat)
library(future)
library(ggplot2)
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 100 * 1024 ^ 3)

args = commandArgs(trailingOnly=TRUE)
disease=args[1]
date='20211228'

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

anno_rna=read.table(paste(data_dir,'/132_snRNA_samples_MetaData.20211215.tsv',sep=''),sep='\t',header=T)
anno_rna_c=anno_rna[anno_rna$Disease==disease,]


for (sample in samples_multiome){
obj=readRDS(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/0.Individual_objects/',
'Multiome/out/',disease,'/',sample,'/',sample,'_processed_multiome.rds',sep=''))

anno_c_sample=anno_c[anno_c$Sample==sample,]
piece_id_atac=unique(anno_c_sample$Piece_ID)

rownames(anno_c_sample)=gsub('(.*)_(.*)_(.*)','\\3',anno_c_sample$V1)
obj_barc=as.data.frame(obj$orig.ident)
anno_c_sample=anno_c_sample[rownames(obj_barc),]
obj$cell_type_ATAC=anno_c_sample$cell_type.harmonized.cancer

piece_id_rna=piece_id_atac
anno_rna_c_sample=anno_rna_c[anno_rna_c$Piece_ID==piece_id_rna,]
rownames(anno_rna_c_sample)=gsub('(.*)_(.*)_(.*)','\\3',anno_rna_c_sample$Barcode)
obj_barc=as.data.frame(obj$orig.ident)
anno_rna_c_sample=anno_rna_c_sample[rownames(obj_barc),]
obj$cell_type_RNA=anno_rna_c_sample$cell_type.harmonized.cancer


n_cells=length(obj$atac_barcode[obj$cell_type_RNA=='Tumor' & obj$cell_type_ATAC=='Tumor' & 
!is.na(obj$cell_type_RNA) & !is.na(obj$cell_type_RNA)])

if (n_cells>200){
obj_s=subset(obj,cell_type_RNA %in% c('Tumor') & cell_type_ATAC %in% c('Tumor') )

##################################################
########Gene expression data processing###########
##################################################

#Use params as in the Dan's script:
DefaultAssay(obj_s) <- "RNA"
obj_s <- SCTransform(obj_s, vars.to.regress = c("nCount_RNA","percent.mito"),return.only.var.genes = F)
obj_s <- RunPCA(obj_s, npcs=30,verbose = FALSE)
obj_s = RunUMAP(obj_s, reduction = 'pca',dims = 1:30,reduction.name = 'rna.umap',reduction.key = 'rnaUMAP_')


##################################################
########DNA accessibility data processing#########
##################################################

###Now re-cluster, ATAC:
DefaultAssay(obj_s) = "X500peaksMACS2"
obj_s <- FindTopFeatures(obj_s, min.cutoff = 'q0')
obj_s <- RunTFIDF(obj_s)
obj_s <- RunSVD(obj_s)

obj_s <- RunUMAP(object = obj_s, reduction = 'lsi', dims = 2:30,reduction.name = "atac.umap",
reduction.key = "atacUMAP_")


# build a joint neighbor graph using both assays
obj_s = FindMultiModalNeighbors(object = obj_s,reduction.list = list("pca", "lsi"),
dims.list = list(1:30, 2:30),verbose = TRUE)

# build a joint UMAP visualization
obj_s = RunUMAP(object = obj_s, nn.name = "weighted.nn",reduction.name = "wnn.umap",reduction.key = "wnnUMAP_",
verbose = TRUE)

obj_s = FindClusters(obj_s,graph.name = "wsnn",algorithm = 3,verbose = TRUE)

p1 = DimPlot(obj_s, label = TRUE, repel = TRUE, reduction = "wnn.umap") +
NoLegend() +ggtitle('Joint WNN')
p2 = DimPlot(obj_s, reduction = 'rna.umap', label = TRUE,
repel = TRUE, label.size = 2.5) +NoLegend() + ggtitle('RNA')
p3 = DimPlot(obj_s, reduction = 'atac.umap', label = TRUE,
repel = TRUE, label.size = 2.5) +NoLegend() +ggtitle('ATAC')
p = p1 + p2 + p3

pdf(paste('out/TumorOnly.RDS.dimplots/',disease,'/',sample,'_TumorOnly.20211129.pdf',sep=""),height=4,width=12)
print(p)
dev.off()

saveRDS(obj_s,paste('out/TumorOnly.RDS/',disease,'/',sample,'_TumorOnly.20211129.rds',sep=''))
print(sample)

}
}
