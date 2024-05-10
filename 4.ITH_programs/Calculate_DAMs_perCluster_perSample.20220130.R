#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(reshape)
library(reshape2)
library(EnsDb.Hsapiens.v86)
library(plyr)
library(dplyr)

atac=readRDS(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/3.Merge_snATAC/',
'Merge.vers.20211020/PanCan_object/Tumor_Normal.v.20211026/149_snATAC_148K_peaks_TumorNormal.motifsAdded.',
'chromvar.20211029.rds.gz',sep=''))

###Use latest cell type annotation

annot=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/3.Merge_snATAC/Merge.',
'vers.20211020/All_149_samples_metadata_data_freeze_v4.1.tsv',sep=''),sep='\t',header=T)
orig_1=as.data.frame(atac$Piece_ID)
rownames(annot)=annot[,1]
annot=annot[rownames(orig_1),]
atac$cell_type.harmonized.cancer=annot$cell_type.harmonized.cancer
atac$Disease=ifelse(atac$Piece_ID %in% basal_piece_ids,, "BRCA_Basal",
atac$Disease)

meta_snatac=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/',
'8.Tumor_heterogeneity/2.snATAC_TumorOnly_IndObj/snATAC/Cluster_label_transfer/TumorOnly.RDS.dimplots.',
'RNAClustersReassigned/metadata/RNA_clusters_TumorOnly.snATAC.tsv',sep=''),sep='\t',header=T)
meta_mult=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/',
'8.Tumor_heterogeneity/2.snATAC_TumorOnly_IndObj/snATAC/Cluster_label_transfer/Multiome_clustersManual_',
'TumorOnly.snATAC.20220101.tsv',sep=''),sep='\t',header=T)
meta_snatac$Data_type='snATAC'
meta_mult$Data_type='Multiome'

meta_tab=rbind(meta_snatac, meta_mult)
write.table(meta_tab, "Tumor_ManualClusters_snATAC_multiome.20220101.tsv",sep='\t',quote=F,row.names=F)

meta_tab=read.table("../../3.DAM_analysis/Tumor_ManualClusters_snATAC_multiome.20220101.tsv",sep='\t',header=T)
orig_1=as.data.frame(atac$Piece_ID)
rownames(meta_tab)=paste(meta_tab$Piece_ID_ATAC, meta_tab$Barcode,sep='_')
meta_tab=meta_tab[rownames(orig_1),]

ATAC=atac
ATAC$Cluster_ID=meta_tab$cell_group

ATAC$Cluster_ID_2=ifelse(!is.na(ATAC$Cluster_ID),paste(ATAC$Piece_ID, ATAC$Cluster_ID,sep='_'),'Other')

ATAC$test=ATAC$Cluster_ID_2

DefaultAssay(ATAC) <- 'chromvar'

chromv= GetAssayData(object = ATAC)

jaspar=read.table('/diskmnt/Projects/HTAN_analysis/snATAC/Signac/CCRCC/JASPAR2020_motifs.txt',sep='\t',
header=TRUE)

mtx0=chromv
res=merge(mtx0,jaspar,by=0,all.x=TRUE)
rownames(res)=res$motif.name
res=res[,-1]
res=res[,1:(ncol(res)-2)]


ann_col0 = data.frame('cell_types' = ATAC$test,'Piece_ID'=ATAC$Piece_ID)
ann_col0$cel_2=ann_col0$cell_types
ann_col0=ann_col0[order(ann_col0$cell_types),]
ann_col1=data.frame("cell_types"=ann_col0$cell_types, "Piece_ID"=ann_col0$Piece_ID)
rownames(ann_col1)=rownames(ann_col0)

res=res[,rownames(ann_col0)]

final_wilcoxon_stat=NULL
final_wilcoxon_stat_all=NULL

piece_ids=unique(ATAC$Piece_ID)

for (piece_id in piece_ids){
    cell_types=unique(ATAC$test[ATAC$Piece_ID==piece_id])
    cell_types=cell_types[cell_types!='Other']
    if(length(cell_types)>1){
    for (cell_type in cell_types){
    if(length(rownames(ann_col1)[ann_col1$cell_types==cell_type])>=50){
    res_1=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_type]]
    res_2=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types!=cell_type & 
ann_col1$Piece_ID==piece_id]]
    all_wilcoxon_stat=NULL
        for (motif in 1:nrow(res)){
           mean_score1=mean(as.numeric(as.character(unlist(res_1[motif,]))),na.rm=TRUE)
           mean_score2=mean(as.numeric(as.character(unlist(res_2[motif,]))),na.rm=TRUE)
           w_test=wilcox.test(as.numeric(as.character(unlist(res_1[motif,]))),
as.numeric(as.character(unlist(res_2[motif,]))))
           stat=cbind(cell_type,rownames(res)[motif],mean_score1,mean_score2,w_test$p.value,piece_id)
           all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
        all_wilcoxon_stat$V5=as.numeric(as.character(unlist(all_wilcoxon_stat$V5)))
	all_wilcoxon_stat=all_wilcoxon_stat[order(all_wilcoxon_stat$V5),]
        all_wilcoxon_stat$mean_score1=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score1)))
        all_wilcoxon_stat$mean_score2=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score2)))
        all_wilcoxon_stat$diff=all_wilcoxon_stat$mean_score1-all_wilcoxon_stat$mean_score2
        all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$diff),]

colnames(all_wilcoxon_stat)[1:2]=c('cell_type','TF')
colnames(all_wilcoxon_stat)[5]=c('pvalue')
all_wilcoxon_stat$FDR=p.adjust(all_wilcoxon_stat$pvalue, method="fdr")
all_wilcoxon_stat$cell_t2='Other'

final_wilcoxon_stat=rbind(final_wilcoxon_stat,all_wilcoxon_stat)
}
}
print(paste(piece_id,cell_type,sep=' '))
}
}
final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
colnames(final_wilcoxon_stat)[1:2]=c('cell_t1','TF_Name')



write.table(final_wilcoxon_stat,paste("out/DAMs_perCluster_perSample.ITH.20220131.tsv"),
quote=FALSE,sep="\t",row.names=FALSE)



