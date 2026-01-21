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

#ISC cell
rna_obj_filter <- readRDS("rna_obj_filter.rds")
cells.sub <- subset(rna_obj_filter@meta.data, wsnn_res.0.5%in%c(1 ,2 ,5 ,6 ,7 ,8 ,9 ,10))
#重新导入数据
inputdata.10x <- Read10X_h5("~/M3F-001_result_cellranger-arc_count/outs/filtered_feature_bc_matrix.h5")
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
rna_counts <- rna_counts[,row.names(cells.sub)]
rna_filter_sub <- CreateSeuratObject(counts = rna_counts)
rna_filter_sub[["percent.mt"]] <- PercentageFeatureSet(rna_filter_sub, pattern = "^mt-")

#RNA analysis
DefaultAssay(rna_filter_sub) <- "RNA"
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
library(stringr)
s.genes <- str_to_title(s.genes)
g2m.genes <- str_to_title(g2m.genes)
rna_filter_sub <- SCTransform(rna_filter_sub, vars.to.regress = c("percent.mt","nFeature_RNA","nCount_RNA"),assay = 'RNA', new.assay.name = 'SCT')
rna_filter_sub <- CellCycleScoring(rna_filter_sub, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
rna_filter_sub$CC.Difference <- rna_filter_sub$S.Score - rna_filter_sub$G2M.Score
rna_filter_sub <- SCTransform(rna_filter_sub, vars.to.regress = "CC.Difference", verbose = T,assay = 'RNA', new.assay.name = 'SCT')
rna_filter_sub <- RunPCA(rna_filter_sub)
rna_filter_sub <- RunUMAP(rna_filter_sub,dims = 1:30, reduction.name = 'subumap.rna', reduction.key = 'subrnaUMAP_')

# ATAC analysis
atac_counts <- inputdata.10x$Peaks
atac_counts <- atac_counts[,row.names(cells.sub)]

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

frag.file <- "/bios-store5/data_mouse_Intestine_scRNA-co-scATAC/M3F-001_result_cellranger-arc_count/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
rna_filter_sub[["ATAC"]] <- chrom_assay

DefaultAssay(rna_filter_sub) <- "ATAC"
rna_filter_sub <- RunTFIDF(rna_filter_sub)
rna_filter_sub <- FindTopFeatures(rna_filter_sub, min.cutoff = 'q0')
rna_filter_sub <- RunSVD(rna_filter_sub)
# We exclude the first dimension as this is typically correlated with sequencing depth
rna_filter_sub <- RunUMAP(rna_filter_sub, reduction = 'lsi', dims = 2:30, reduction.name = "subumap.atac", reduction.key = "subatacUMAP_")

#WNN
rna_filter_sub <- FindMultiModalNeighbors(rna_filter_sub, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30),k.nn=5)
rna_filter_sub <- RunUMAP(rna_filter_sub, nn.name = "weighted.nn", reduction.name = "subwnn.umap", reduction.key = "subwnnUMAP_")
rna_filter_sub <- FindClusters(rna_filter_sub, graph.name = "wsnn", algorithm = 4, verbose = FALSE,resolution = 1)
rna_filter_sub <- FindSubCluster(rna_filter_sub, cluster = 5, graph.name = "wsnn", algorithm = 4,resolution = 0.5)
Idents(rna_filter_sub) <- "sub.cluster"

# call peaks using MACS2
peaks <- CallPeaks(rna_filter_sub)
peaks <- keepStandardChromosomes(peaks, species="Mus_musculus", pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(rna_filter_sub),
  features = peaks,
  cells = colnames(rna_filter_sub)
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
rna_filter_sub[["peaks"]] <- peakchromassay

#motif score
DefaultAssay(rna_filter_sub) <- "peaks"
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
      features = StringToGRanges(rownames(rna_filter_sub), sep = c("-", "-")),
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
rna_filter_sub <- SetAssayData(rna_filter_sub, assay = 'peaks', slot = 'motifs', new.data = motif.object)
DefaultAssay(rna_filter_sub) <- "peaks"
rna_filter_sub <- RunChromVAR(
  object = rna_filter_sub,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  assay = "peaks"
)


saveRDS(rna_filter_sub,"rna_filter_sub.rds")
