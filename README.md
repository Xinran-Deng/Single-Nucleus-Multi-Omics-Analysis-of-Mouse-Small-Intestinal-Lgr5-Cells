# Single-Nucleus-Multi-Omics-Analysis-of-Mouse-Small-Intestinal-Lgr5-Cells
This is the code collection for analysing single-cell multi-omics data in the article entitled "Single-nucleus multi-omics analysis of mouse small-intestinal Lgr5+ cell populations reveals Foxa3-induced Paneth cell-lineage differentiation."


1.load_QC.R: Import RNA+ATAC dual-modal single-cell sequencing data and perform quality control.

2.WNN.R: Perform dimensionality reduction and clustering using the Weighted Nearest Neighbor (WNN) algorithm.

3.WNN_ISCsub.R: Extract and reanalyse ISC cells.

4.CytoTrace2.R: Evaluate stemness of each cell subtype using CytoTrace2.

5.monocle3.R: Conduct pseudo-time analysis with Monocle3.

6.pyscenic.sh: Infer subtype-specific transcription factors for ISCs using Pyscenic.

7.TF_correlation.R: Conduct consistency tests for gene expression and chromatin accessibility of candidate transcription factors for ISC subtypes.
