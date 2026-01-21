library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(ggplot2)
library(hdf5r)
library(biovizBase)
library(ggsci)
library(RColorBrewer)

inputdata.10x <- Read10X_h5("~/M3F-001_result_cellranger-arc_count/outs/filtered_feature_bc_matrix.h5")
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

#去除核糖体基因
rb.genes <- rownames(rna_counts)[grep("^Rp[sl]",rownames(rna_counts))]
rna_counts <- rna_counts[-which(rownames(rna_counts) %in% rb.genes),]

# Create Seurat object
rna_obj <- CreateSeuratObject(counts = rna_counts)
rna_obj[["percent.mt"]] <- PercentageFeatureSet(rna_obj, pattern = "^mt-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
frag.file <- "~/M3F-001_result_cellranger-arc_count/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'mm10',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
rna_obj[["ATAC"]] <- chrom_assay

DefaultAssay(rna_obj) <- "ATAC"
rna_obj <- NucleosomeSignal(rna_obj)
rna_obj <- TSSEnrichment(rna_obj)

#filter low quality cell
pdf("QC_atac.pdf")
VlnPlot(
  object = rna_obj,
  features = c("nCount_RNA","nFeature_RNA","nCount_ATAC","nFeature_ATAC","percent.mt", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
dev.off()

rna_obj_filter <- subset(
  x = rna_obj,
  subset = nCount_ATAC < 1e5 &
    nCount_ATAC > 1e3 &
    nCount_RNA < 50000 &
    nCount_RNA > 2000 &
    percent.mt < 20 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2
)