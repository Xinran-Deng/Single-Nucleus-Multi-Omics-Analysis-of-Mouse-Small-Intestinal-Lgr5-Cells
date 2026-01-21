library(CytoTRACE2)
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(hdf5r)
library(biovizBase)
library(ggsci)
library(RColorBrewer)
library(ggpubr)
library(forcats)
rna_obj_filter <- readRDS("rna_obj_filter.rds")
rna_filter_sub <- readRDS("rna_filter_sub.rds")

cytotrace2_result <- cytotrace2(rna_obj_filter,
                  species = "mouse",
                  is_seurat = TRUE,
                  slot_type = "counts",
                  #full_model = TRUE,
                  batch_size = NULL,
                  smooth_batch_size = 1000,
                  parallelize_models = TRUE,
                  parallelize_smoothing = TRUE,
                  ncores = NULL,
                  #max_pcs = 25,
                  seed = 777)
cytotrace2_result$combined <- Idents(cytotrace2_result)  
pdf("CytoTrace2_combinedISC.pdf",width=10,height=5)
ggboxplot(cytotrace2_result@meta.data, x="ISCcombined", y="CytoTRACE2_Score", width = 0.6, 
                color = "black",#轮廓颜色
                fill="ISCcombined",#填充
                palette = "simpsons",
                xlab = T, #不显示x轴的标签
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.5, #误差条大小
                outlier.shape=NA, #不显示outlier
                legend = "right")+
                scale_fill_manual(values=c("#c5ac84",
                              "#197EC0FF","#F05C3BFF","#46732EFF","#71D0F5FF","#370335FF","#075149FF"))+
                stat_compare_means()
dev.off()

rna_obj_filter_sub <- subset(rna_obj_filter,idents=c("Goblet","Enteroendocrine","Tuft"),invert=TRUE)
cytotrace2_result115 <- cytotrace2(rna_obj_filter_sub,
                  species = "mouse",
                  is_seurat = TRUE,
                  slot_type = "counts",
                  #full_model = TRUE,
                  batch_size = NULL,
                  smooth_batch_size = 3000,
                  parallelize_models = TRUE,
                  parallelize_smoothing = TRUE,
                  ncores = NULL,
                  #max_pcs = 25,
                  seed = 115)
cytotrace2_result_sub <- cytotrace2_result115@meta.data[cytotrace2_result115@meta.data$combined %in% levels(rna_obj_filter$combined)[1:6],] 
cytotrace2_result_sub <- subset(cytotrace2_result115,combined %in% levels(rna_obj_filter$combined)[1:6]) 
Idents(cytotrace2_result115) <- fct_relevel(Idents(cytotrace2_result),"Lgr5-low SC","Top2a+ TA","CBC","Hmgb2+ TA","Lgr5+ LRC","AP","pre-Paneth","Paneth","Goblet","Enteroendocrine","Tuft","Enterocyte")  

pdf("CytoTrace2_all_subISC_test_order.pdf",width=10,height=4)
Cytotrace2_p <- ggboxplot(cytotrace2_result_sub@meta.data, x="combined", y="CytoTRACE2_Score", width = 0.6, 
                color = "black",#轮廓颜色
                fill="combined",#填充
                palette = "simpsons",
                xlab = F, #不显示x轴的标签
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.5, #误差条大小
                outlier.shape=NA, #不显示outlier
                legend = "right")+ylim(c(0.5,0.66))+
                scale_fill_manual(values=c("#fed439",
                              "#8a9197","#fd7446","#d2af81","#d5e4a2","#709ae1"))+
                stat_compare_means()
Cytotrace2_p
dev.off() 

##mean
df <- cytotrace2_result_sub
df <- aggregate(df$CytoTRACE2_Score,by=list(type=df$combined),mean)
df <- df[order(df$x),]
Cytotrace2_df <- df %>% cbind(ggplot_build(Cytotrace2_p)$data[[1]])
colnames(Cytotrace2_df)[2] <- "mean"
pdf("CytoTrace2_all_subISC_test_meanorder.pdf",width=10,height=5)
Cytotrace2_p+
  geom_segment(data=Cytotrace2_df,
               aes(x=xmin,xend=xmax,
                   y=mean,
                   yend=mean),
               color="red")
dev.off()