library(ggpointdensity)
library(ggpubr)
library(viridis)
library(ggsci)
rna_filter_sub <- readRDS("rna_filter_sub.rds")
markers_motifs_sub_50kb <- presto:::wilcoxauc.Seurat(X = rna_filter_sub, group_by = 'celltype', assay = 'data', seurat_assay = 'chromvar')
##DEG&motif
t1 <- intersect(markers_motifs_sub_50kb_up[markers_motifs_sub_50kb_up$group=='Lgr5-low SC','gene'],markers_motifs_sub_50kb_down[markers_motifs_sub_50kb_down$group=='Lgr5+ LRC','gene'])
t2 <- intersect(markers_motifs_sub_50kb_up[markers_motifs_sub_50kb_up$group=='Lgr5-low SC','gene'],markers_motifs_sub_50kb_down[markers_motifs_sub_50kb_down$group=='AP','gene'])
t3 <- intersect(markers_motifs_sub_50kb_up[markers_motifs_sub_50kb_up$group=='Lgr5+ LRC','gene'],markers_motifs_sub_50kb_down[markers_motifs_sub_50kb_down$group=='Lgr5-low SC','gene'])
t4 <- intersect(markers_motifs_sub_50kb_up[markers_motifs_sub_50kb_up$group=='AP','gene'],markers_motifs_sub_50kb_down[markers_motifs_sub_50kb_down$group=='Lgr5-low SC','gene'])
test <- c(t1,t2,t3,t4)

###AUC motif
regulonAUC <- read.csv("/bios-store5/data_mouse_Intestine_scRNA-co-scATAC/TFfind/regulonAUC.csv",row.names=1)
lm_df1 <- t(rna_filter_sub@assays$chromvar[motif_TF[motif_TF$Transcription_factor %in% test,'Model'],])
lm_df1 <- cbind(lm_df1,t(regulonAUC[paste0(test,'(+)'),]))
lm_df1 <- as.data.frame(lm_df1)
colnames(lm_df1) <- gsub("\\(.*?\\)","",colnames(lm_df1))
colnames(lm_df1) <- gsub("-",".",colnames(lm_df1))
lm_df1$celltype <- rna_filter_sub@meta.data[rownames(lm_df1),'celltype']

###mean(AUC) mean(motif)
lm_df2 <- aggregate(lm_df1[,1:11], by=list(type=lm_df1$celltype),mean)
t <- aggregate(lm_df1[,12:20], by=list(type=lm_df1$celltype),mean)
lm_df2 <- cbind(lm_df2,t[,-1])
colnames(lm_df2) <- gsub("\\(.*?\\)","",colnames(lm_df2))

###RSS mean(motif)
tmp <- rssPlot$df
tmp <- reshape2::dcast(tmp,Topic~tmp$cellType,value.var='RSS')
rownames(tmp) <- tmp[,1]
tmp <- tmp[,-1]
tmp <- as.matrix(tmp)
lm_df2 <- aggregate(lm_df1[,1:8], by=list(type=lm_df1$celltype),mean)
lm_df2 <- cbind(lm_df2,t(tmp[paste0(test,"(+)"),]))
colnames(lm_df2) <- gsub("\\(.*?\\)","",colnames(lm_df2))

###correlation
for (i in test){
    motif <- gsub("-",".",motif_TF[motif_TF$Transcription_factor==i,'Model'])
    for(j in motif){
        #整体的绘制
        p1 <- ggplot(data = lm_df1, mapping = aes(x=lm_df1[,i],y=lm_df1[,j]))+
               geom_pointdensity()+
               scale_color_viridis()+ #viridis (version 0.5.1)
               xlab(paste0("regulonAUC of ",i))+
               ylab(j)+
               theme_classic()+
               geom_smooth(method = 'lm', color = 'black', fill = 'lightgray')+
               stat_cor(method = 'pearson',label.x.npc = "left",label.y.npc = "top")
        #分组绘制
        p2 <- ggplot(data = lm_df1, mapping = aes(x=lm_df1[,i],y=lm_df1[,j],color=celltype))+
               geom_pointdensity()+
               scale_color_simpsons()+
               xlab(paste0("regulonAUC of ",i))+
               ylab(j)+
               theme_classic()+
               geom_smooth(method = 'lm', color = 'black', fill = 'lightgray')+
               stat_cor(method = 'pearson',label.x.npc = "left",label.y.npc = "top")
        p3 <- ggplot(data = lm_df1, mapping = aes(x=lm_df1[,i],y=lm_df1[,j],color=celltype))+
               geom_pointdensity()+
               scale_color_simpsons()+ #scale_color_manual(values = pal_simpsons()(6)[c(1,3,5,4,6,2)]) +
               xlab(paste0("regulonAUC of ",i))+
               ylab(j)+
               theme_classic()+
               geom_smooth(method = 'lm', color = 'black', fill = 'lightgray')+
               stat_cor(method = 'pearson',label.x = 0,label.y = 6)+
               facet_wrap(~celltype,nrow = 2)
               
        print(paste0("TF:",i,"    motif:",j))
        lm_result1 <- ggplot_build(p1)$data[[3]]
        lm_result2 <- ggplot_build(p2)$data[[3]]
        cat("Total  R:",lm_result1$r,"  p.vaule:",lm_result1$p,'\n')
        cat("  ",paste0("cluster",c(5,7,1,8,2,4,3)),'\n')
        cat("R:",lm_result2$r,'\n')
        cat("p:",lm_result2$p,"\n")
        
        #if(all(lm_result1$r>0,lm_result2$r>0,lm_result1$p < 0.05,lm_result2$p < 0.05)){
        if(all(lm_result1$r>0,lm_result1$p < 0.05)){
            #pdf(paste0(i,"-",j,".pdf"),height=12,width=16)
            #p1|p2
            #p3
            ggarrange(ggarrange(p1,p2,ncol = 2),p3,nrow = 2,labels=i,label.x=0.45)
            ggsave(paste0(i,"-",j,".pdf"),height=16,width=16)
            #dev.off() 
        }
        p4 <- ggscatter(data = lm_df2, x=i,y=j,color="type",
              add = "reg.line", conf.int = TRUE,
              add.params = list(color = "red", fill = 'lightgray'),
              cor.coef = TRUE,
              cor.coeff.args = list(method = 'pearson',label.x.npc = "left",label.y.npc = "top",label.sep = "\n"))+
              scale_color_simpsons()+
              xlab(paste0(i,"(+)[z.score]"))+
              ylab(paste0(j,"[mean]"))
        lm_result4 <- ggplot_build(p4)$data[[3]]
        cat("==============================",'\n')
        cat("cluster R:",lm_result4$r,'\n')
        cat("cluster p:",lm_result4$p,"\n")
        if(all(lm_result4$r>0,lm_result4$p < 0.05)){
            #pdf(paste0(i,"-",j,"_clusterCorrelation.pdf"),height=12,width=16)
            p4
            ggsave(paste0(i,"-",j,"_clusterCorrelation.pdf"))
            #dev.off()  
        }
    }
}
