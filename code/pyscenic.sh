##R：生成loom文件
#track_genes from monocle3.R
rna_filter_sub <- readRDS("rna_filter_sub.rds")
rna_filter_sub_trackgene <- subset(x=rna_filter_sub, features = rownames(track_genes))
write.csv(as.matrix(rna_filter_sub_trackgene@meta.data),file = "metadata_rna_filter_sub.csv")
library(SeuratDisk)
rna_filter_sub_trackgene.loom <- as.loom(rna_filter_sub_trackgene,filename = "rna_filter_sub_trackgene.loom",assay = "RNA")
rna_filter_sub_trackgene.loom$close_all()

##shell：python
pyscenic grn \
         --num_workers 20 \
         --output adj.TF_genetop200_5837.tsv \
         --method grnboost2 \
         --cell_id_attribute combined \
         rna_filter_sub_trackgenetop200_5837.loom \
         mm_mgi_tfs.txt \
         >pyscenic.grntop200_5837.out \
         2>pyscenic.grntop200_5837.err
#pyscenic 的3个步骤之 cistarget
pyscenic ctx \
          adj.TF_genetop200_5837.tsv \
          mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
          --annotations_fname motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
          --expression_mtx_fname rna_filter_sub_trackgenetop200_5837.loom \
          --cell_id_attribute combined \
          --mode "dask_multiprocessing" \
          --output reg_motifenrichtop200_5837.csv \
          --num_workers 10 \
          --mask_dropouts \
          >pyscenic.cistargettop200_5837.out \
          2>pyscenic.cistargettop200_5837.err
#pyscenic 的3个步骤之 AUCell
pyscenic aucell \
          rna_filter_sub_trackgenetop200_5837.loom \
          reg_motifenrichtop200_5837.csv \
          --cell_id_attribute combined \
          --output rna_filter_sub_trackgenetop200_5837_SCENIC.loom \
          --num_workers 10 \
          >pyscenic.AUCelltop200_5837.out \
          2>pyscenic.AUCelltop200.err

##R：RSS与可视化
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(pheatmap)
library(data.table)
library(SCENIC)

scenicLoomPath="~/rna_filter_sub_trackgene_SCENIC.loom"
loom <- open_loom(scenicLoomPath) 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)

meta <- read.csv("c",header=T,stringsAsFactor=F,row.names=1)
cellinfo <- meta[,c("combined","nFeature_RNA","nCount_RNA")]
cellinfo <- rna_filter_sub@meta.data[,c("combined83","nFeature_RNA","nCount_RNA")]
colnames(cellinfo)=c('celltype', 'nGene' ,'nUMI')
cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"

sub_regulonAUC <- regulonAUC
#SCENIC包自己提供了一个 calcRSS函数，帮助我们来挑选各个单细胞亚群特异性的转录因子，全称是：Calculates the regulon specificity score
rss <- calcRSS(AUC=getAUC(sub_regulonAUC),
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])
rss=na.omit(rss) 

rssPlot <- plotRSS(rss, zThreshold=0.8,cluster_columns = TRUE)
#提取具体数据
rssDf <- rssPlot$df
rssDf <- reshape2::dcast(rssDf,Topic~rssDf$cellType,value.var='Z')
rownames(rssDf) <- rssDf[,1]
rssDf <- rssDf[,-1]
rssDf_matrix <- as.matrix(rssDf)




