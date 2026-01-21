library(Seurat)
library(monocle3)
library(ggplot2)
rna_obj_filter <- readRDS("rna_obj_filter.rds")
rna_filter_sub <- readRDS("rna_filter_sub.rds")
get_earliest_principal_node <- function(cds,type){
  cell_ids <- which(colData(cds)[, "combined"] == type)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  #igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex))))]
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

##overall
exp <- as.sparse(rna_obj_filter@assays$RNA@counts)
cell_metadata <-rna_obj_filter@meta.data
gene_annotation <- data.frame(gene_short_name = row.names(exp),row.names = row.names(exp))
cds <- new_cell_data_set(exp,
			     cell_metadata = cell_metadata,
			     gene_metadata = gene_annotation)
cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds)
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(rna_obj_filter,reduction = "wnn.umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed 
cds <- cluster_cells(cds)
cds@clusters$UMAP$clusters <- Idents(rna_obj_filter)[rownames(colData(cds))]
cds <- learn_graph(cds,use_partition =FALSE,close_loop=FALSE,learn_graph_control=list(rann.k=NULL,ncenter=700))
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds,'5'))

##ISC
exp_sub <- as.sparse(rna_filter_sub@assays$RNA@counts)
cell_metadata_sub <- rna_filter_sub@meta.data
gene_annotation_sub <- data.frame(gene_short_name = row.names(exp_sub),row.names = row.names(exp_sub))
cds_sub <- new_cell_data_set(exp_sub,
			     cell_metadata = cell_metadata_sub,
			     gene_metadata = gene_annotation_sub)
cds_sub <- preprocess_cds(cds_sub)
cds_sub <- reduce_dimension(cds_sub)
cds_sub.embed <- cds_sub@int_colData$reducedDims$UMAP
int.embed <- Embeddings(rna_filter_sub,reduction = "subwnn.umap")
int.embed <- int.embed[rownames(cds_sub.embed),]
cds_sub@int_colData$reducedDims$UMAP <- int.embed 
cds_sub <- cluster_cells(cds_sub)
cds_sub@clusters$UMAP$clusters <- Idents(rna_filter_sub)[rownames(colData(cds_sub))]
cds_sub <- learn_graph(cds_sub,close_loop=FALSE,learn_graph_control=list(rann.k=NULL,ncenter=700))
cds_sub <- order_cells(cds_sub, root_pr_nodes=get_earliest_principal_node(cds_sub,'5'))

#寻找拟时轨迹差异基因
track_genes <- graph_test(cds_sub, neighbor_graph="principal_graph",cores=6)
track_genes <- track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3) 
track_genes <- track_genes %>% arrange(desc(abs(morans_I)))
track_genes_sig <- track_genes %>% top_n(n=10,morans_I) %>% pull(gene_short_name) %>% as.character()
plot_genes_in_pseudotime(cds[track_genes_sig,],color_cells_by="celltype",min_expr=0.5)
plot_cells(cds, genes=track_genes_sig,
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
###collect the trajectory-variable genes into modules
pr_deg_ids <- row.names(subset(track_genes, morans_I>0.1&q_value == 0))
gene_module_df_sub1 <- find_gene_modules(cds_sub[pr_deg_ids,], resolution=c(10^seq(-6,-1)))
cell_group_df_sub1 <- tibble::tibble(cell=row.names(colData(cds_sub)), 
                                cell_group=colData(cds_sub)$combined_split8)
agg_mat_sub1 <- aggregate_gene_expression(cds_sub, gene_module_df_sub1, cell_group_df_sub1)
row.names(agg_mat_sub1) <- stringr::str_c("Module ", row.names(agg_mat_sub1))
colnames(agg_mat_sub1) <- c("CBC","AP","Top2a+ TA","Hmgb2+ TA","SP","Lgr5+ LRC")
pdf("heatmap_monocle3_trackgenemodule.pdf")
pheatmap::pheatmap(agg_mat_sub1,
                   scale="column", clustering_method="ward.D2")
dev.off()
