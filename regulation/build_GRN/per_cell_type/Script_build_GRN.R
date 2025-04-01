# import packages
library(Seurat)
library(Signac)
library(tidyverse)
library(Pando)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(doParallel)
registerDoParallel(10)
library(dplyr)
library(ggplot2)

# set working folder
setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/GRN_per_cell_type")

# load data
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
ct <- args[1]
print(paste0("Read Res_chip_infection_multiome_", ct, ".rds"))
combined <- readRDS(paste0("Res_chip_infection_multiome_", ct, ".rds"))
# load motif and TF info
data('motifs')
data('motif2tf')

# infer gene regulatory network associated with infection
# get pre-selected peaks that are both conserved or overlap with ENCODE regions
data('phastConsElements20Mammals.UCSC.hg38')
data('SCREEN.ccRE.UCSC.hg38')
gr_preselected_regions <- union(phastConsElements20Mammals.UCSC.hg38, SCREEN.ccRE.UCSC.hg38)
## run LinkPeaks, high resolution clustering, hvg identification for each sample
DefaultAssay(combined) <- "peaks_merged"
combined <- RunTFIDF(combined)
# perform high resolution clustering on the original scATAC-seq data 
print("perform high resolution clustering on the original scATAC-seq data")
det_rates <- rowMeans(combined@assays$peaks_merged@counts > 0) 
top <- names(which(det_rates>0.01))
combined <- FindTopFeatures(combined, min.cutoff = 'q0')
VariableFeatures(combined) <- top
combined <- RunSVD(combined)
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:50)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3, resolution = 10)
saveRDS(combined, file=paste0("Res_chip_infection_",ct,"_with_ATAC_HighResClustering.rds"))

#  link peaks to nearby genes based on correlation between gene expression and peak detection rate
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(object = combined, normalization.method = "LogNormalize", scale.factor = 1e4)
combined <- FindVariableFeatures(object = combined, selection.method = "vst", nfeatures = 5000)
# load DEGs induced upon infection
combined_sig_res <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/Res_Salmonella_infection_on_colon_epi_goblet_and_s2E_DEG_sig_res.rds")
hvg <- union(VariableFeatures(combined), combined_sig_res$feature)

# first compute the GC content for each peak
DefaultAssay(combined) <- "peaks_merged"
combined <- RegionStats(combined, genome = BSgenome.Hsapiens.UCSC.hg38)
# link peaks to genes
print("link peaks to genes")
combined <- LinkPeaks(
  object=combined,
  peak.assay="peaks_merged",
  expression.assay="RNA",
  genes.use = hvg
)
saveRDS(combined, file=paste0("Res_chip_infection_",ct,"_with_LinkPeaks.rds"))

# only use the peaks with significant correlation with nearby genes
print("initiate gene regulatory network inference")
combined <- initiate_grn(
  combined,
  rna_assay = 'RNA',
  peak_assay = 'ATAC_binarized',
  regions = intersect(gr_preselected_regions, combined@assays$peaks_merged@links)
)

print("find TFBS motifs")
combined <- find_motifs(
  combined, 
  pfm = motifs, 
  motif_tfs = motif2tf,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

print("infer gene regulatory network")
combined <- infer_grn(
  combined,
  peak_to_gene_method = 'Signac',
  genes = hvg,
  parallel = T,
  aggregate_peaks_col = "peaks_merged_snn_res.10"
  
)

print("find modules")
combined <- find_modules(
  combined, 
  p_thresh = 0.05,
  rsq_thresh = 0.1,
  nvar_thresh = 10,
  min_genes_per_module = 5
)

combined <- get_network_graph(combined, graph="module_graph")
saveRDS(combined, file=paste0("Res_chip_infection_",ct,"_with_GRN.rds"))