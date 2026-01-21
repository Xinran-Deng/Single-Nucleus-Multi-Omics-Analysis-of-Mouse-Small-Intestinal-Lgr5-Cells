library(chromVAR)
library(JASPAR2022)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Seurat)
library(Signac)
library(dplyr)
library(universalmotif)
library(rlist)

#RNA analysis
DefaultAssay(rna_obj_filter) <- "RNA"
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
library(stringr)
s.genes <- str_to_title(s.genes)
g2m.genes <- str_to_title(g2m.genes)
rna_obj_filter <- SCTransform(rna_obj_filter, vars.to.regress = c("percent.mt","nFeature_RNA","nCount_RNA"),assay = 'RNA', new.assay.name = 'SCT')
rna_obj_filter <- CellCycleScoring(rna_obj_filter, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
rna_obj_filter$CC.Difference <- rna_obj_filter$S.Score - rna_obj_filter$G2M.Score
rna_obj_filter <- SCTransform(rna_obj_filter, vars.to.regress = "CC.Difference", verbose = T,assay = 'RNA', new.assay.name = 'SCT')
DimPlot(rna_obj_filter, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE,group.by='Phase')
rna_obj_filter <- RunPCA(rna_obj_filter)
pdf("ElbowPlot.pdf")
ElbowPlot(rna_obj_filter,ndims=50)
DimHeatmap(rna_obj_filter,dims=1:20,cells = 500,balanced=TRUE)
dev.off()
rna_obj_filter <- RunUMAP(rna_obj_filter,dims = 1:20, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
DefaultAssay(rna_obj_filter) <- "ATAC"
rna_obj_filter <- RunTFIDF(rna_obj_filter)
rna_obj_filter <- FindTopFeatures(rna_obj_filter, min.cutoff = 'q0')
rna_obj_filter <- RunSVD(rna_obj_filter)
# We exclude the first dimension as this is typically correlated with sequencing depth
rna_obj_filter <- RunUMAP(rna_obj_filter, reduction = 'lsi', dims = 2:20, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#WNN
rna_obj_filter <- FindMultiModalNeighbors(rna_obj_filter, reduction.list = list("pca", "lsi"), dims.list = list(1:20, 2:20))
rna_obj_filter <- RunUMAP(rna_obj_filter, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
DimPlot(rna_obj_filter, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE,group.by='Phase')
rna_obj_filter <- FindClusters(rna_obj_filter, graph.name = "wsnn", algorithm = 4, verbose = FALSE,resolution = 1)

# call peaks using MACS2
peaks <- CallPeaks(rna_obj_filter)
peaks <- keepStandardChromosomes(peaks, species="Mus_musculus", pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(rna_obj_filter),
  features = peaks,
  cells = colnames(rna_obj_filter)
)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "~/M3F-001_result_cellranger-arc_count/outs/atac_fragments.tsv.gz"
peakchromassay <- CreateChromatinAssay(
  counts = macs2_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  min.cells = 10,
  fragments = frag.file,
  annotation = annotations
)
rna_obj_filter[["peaks"]] <- peakchromassay

#motif score
DefaultAssay(rna_obj_filter) <- "peaks"
##HOCOMOCO V11
pathToPWMs = "~/pwm/"
pwmfiles = list.files(pathToPWMs)
empty.list <- list()
for (x in 1:length(pwmfiles)) {
 
    test1 = read_matrix(paste0(pathToPWMs,pwmfiles[x]), 
                      headers = ">", sep = "", positions = "rows")
    lel = convert_motifs(test1, class = "TFBSTools-PFMatrix")
  
    empty.list <- c(empty.list, lel)
  
  }
pfm.list <- do.call(PFMatrixList, empty.list)
motif.matrix <- CreateMotifMatrix(
      features = StringToGRanges(rownames(rna_obj_filter), sep = c("-", "-")),
      pwm = pfm.list,
      genome = 'mm10',
      sep = c("-", "-"),
      use.counts = FALSE
    )
colnames(motif.matrix) = TFBSTools::name(pfm.list)
motif.object <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm.list
)
rna_obj_filter <- SetAssayData(rna_obj_filter, assay = 'peaks', slot = 'motifs', new.data = motif.object)
DefaultAssay(rna_obj_filter) <- "peaks"
rna_obj_filter <- RunChromVAR(
  object = rna_obj_filter,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  assay = "peaks"
)
saveRDS(rna_obj_filter,"rna_obj_filter.rds")

