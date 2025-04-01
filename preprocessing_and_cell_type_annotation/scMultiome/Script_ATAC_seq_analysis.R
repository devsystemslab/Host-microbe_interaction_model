setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/")
library(Seurat)
library(Signac)
library(dplyr)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")

data_folder <- "/home/yuq22/Bacteria_TRM/salmonella/data/count_res/filtered_feature_bc_matrix"
sample_id <- list.files(data_folder)
seu_obj_list <- list()
for(x in sample_id){
  path <- file.path(data_folder, x, "matrix.h5")
  count_list <- Read10X_h5(path)
  rna_count <- count_list$`Gene Expression`
  id <- strsplit(x, split="_")[[1]][3]
  #colnames(rna_count) <- paste(id, colnames(rna_count), sep="_")
  
  # Create Seurat object
  seu_obj <- CreateSeuratObject(counts = rna_count)
  group <- substr(id,1,nchar(id)-1)
  seu_obj$condition <- group

  atac_count <- count_list$Peaks
  #colnames(atac_count) <- paste(id, colnames(atac_count), sep="_")
  # Now add in the ATAC-seq data
  # we'll only use peaks in standard chromosomes
  grange.counts <- StringToGRanges(rownames(atac_count), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_count <- atac_count[as.vector(grange.use), ]
  
  frag.file <- file.path("/home/yuq22/Bacteria_TRM/salmonella/data/count_res/atac", paste0(x, ".tsv.gz"))
  chrom_assay <- CreateChromatinAssay(
    counts = atac_count,
    sep = c(":", "-"),
    genome = "hg38",
    fragments = frag.file,
    min.cells = 1
  )
  seu_obj[["ATAC"]] <- chrom_assay  
  DefaultAssay(seu_obj) <- "ATAC"
  seu_obj_list[[id]] <- seu_obj
} 


# call peaks in each sample separately and get a unified peak list
# requantify the peaks in each sample
seurat_atac_aggr <- merge(
  x = seu_obj_list[[1]],
  y = seu_obj_list[-1],
  add.cell.ids = names(seu_obj_list),
)
saveRDS(seurat_atac_aggr, file="Res_seurat_atac_aggr.rds")

# call peaks in each individual sample
seurat_atac_aggr <- readRDS("Res_seurat_atac_aggr.rds")
seurat_atac_aggr$orig.ident <- sapply(colnames(seurat_atac_aggr), function(x) strsplit(x, split="_")[[1]][1])
peaks <- CallPeaks(
  object = seurat_atac_aggr,
  group.by="orig.ident",
  macs2.path = "/home/yuq22/miniconda3/envs/r-kernel/bin/macs2"
)
saveRDS(peaks, file="Res_human_peaks_called_in_each_sample.rds")

# quantify counts of each peaks of the union list
counts_atac_aggr <- FeatureMatrix(seurat_atac_aggr@assays$ATAC@fragments,
                                  features = peaks,
                                  cells = colnames(seurat_atac_aggr))
assay_peaks_merged_atac <- CreateChromatinAssay(counts_atac_aggr, fragments = seurat_atac_aggr@assays$ATAC@fragments)

# get annotation file necessary for ATAC samples
annotations <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/EnsDB_for_ATAC/data.ens98_annot_for_atac.rds")

Annotation(assay_peaks_merged_atac) <- annotations
seurat_atac_aggr[['peaks_merged']] <- assay_peaks_merged_atac
saveRDS(seurat_atac_aggr, file="Res_seurat_atac_aggr_with_peaks_called_in_individual_sample.rds")

# calculate nucleosome signal and TSS enrichment
seurat_atac_aggr <- readRDS("Res_seurat_atac_aggr_with_peaks_called_in_individual_sample.rds")
Idents(seurat_atac_aggr) <- seurat_atac_aggr$orig.ident
DefaultAssay(seurat_atac_aggr) <- "peaks_merged"
seurat_atac_aggr <- NucleosomeSignal(seurat_atac_aggr)
seurat_atac_aggr <- TSSEnrichment(seurat_atac_aggr)
DefaultAssay(seurat_atac_aggr) <- "RNA"
# calculate Mt%
seurat_atac_aggr[["percent.mt"]] <- PercentageFeatureSet(seurat_atac_aggr, pattern = "^MT-")
saveRDS(seurat_atac_aggr, file="Res_seurat_atac_aggr_with_peaks_called_in_individual_sample_before_filtering.rds")

# Visualize QC metrics as a violin plot
p1 <- SCpubr::do_ViolinPlot(seurat_atac_aggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "nCount_peaks_merged", "nFeature_peaks_merged", "TSS.enrichment", "nucleosome_signal"), ncol = 4)
pdf("Plot_ViolinPlot_QC_metrics_per_sample.pdf" , width = 20, height = 15)
p1
dev.off()

# filter out low quality cells
combined <- subset(
  x = seurat_atac_aggr,
  subset = nCount_peaks_merged < 50000 &
    nCount_RNA < 30000 &
    nCount_peaks_merged > 500 &
    nFeature_RNA > 200 &
    nucleosome_signal < 1 &
    TSS.enrichment > 1
)
saveRDS(combined, file="Dat_filtered_colon_epithelium_infection_multiome.rds")
rm(seurat_atac_aggr)
combined <- readRDS("Dat_filtered_colon_epithelium_infection_multiome.rds")

# load data
combined <- readRDS("Dat_filtered_colon_epithelium_infection_multiome.rds")
# add blacklist ratio
combined$blacklist_ratio <- FractionCountsInRegion(
  object = combined, 
  assay = 'peaks_merged',
  regions = blacklist_hg38_unified
)
bl_cells <- colnames(combined)[which(combined$blacklist_ratio>0.005)]
combined <- subset(
  x = combined,
  subset = blacklist_ratio < 0.005
)
saveRDS(combined, file="Dat_filtered_colon_epithelium_infection_multiome.rds")

## peak annotation
# annotate peaks with ChIPseeker
setwd("/Volumes/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/")
library(Signac)
gr_atac <- readRDS("Res_human_peaks_called_in_each_sample.rds")
library(ChIPseeker) # need to install the most updated ChIPseeker from github, otherwise there would be errors
#library(Signac)
#library(Seurat)
library(AnnotationHub)
library(GenomeInfoDb)
ah <- AnnotationHub()
ensdb.list <- query(ah, c("EnsDb.Hsapiens"))
mat <- cbind(ensdb.list$ah_id, ensdb.list$title)
ens.version <- 98 
id <- mat[which(mat[,2]==paste("Ensembl", ens.version, "EnsDb for Homo sapiens")),1]
ensdb <- ensdb.list[[id]] 
seqlevelsStyle(ensdb) <- 'UCSC'

peakAnno <- annotatePeak(gr_atac,
                         tssRegion=c(-2000, 2000),
                         TxDb = ensdb,
                         annoDb=NULL,
                         overlap="all")
peakAnno@anno$annotation_category <- gsub(" \\(.+\\)$", "", peakAnno@anno$annotation)
anno_df <- as.data.frame(peakAnno@anno)
saveRDS(peakAnno, file="Res_peakAnno_by_ChIPseeker.rds")
anno_df$name <- paste(anno_df$seqnames, anno_df$start, anno_df$end, sep="-")
saveRDS(anno_df, file="Data_frame_peakAnno_by_ChIPseeker.rds")

peak_gene <- anno_df[,c("annotation", "geneId", "name", "annotation_category")]
idx <- which(peak_gene$annotation_category %in% c("Exon", "Intron"))
for(i in which(peak_gene$annotation_category =="Exon")){
  peak_gene$geneId[i] <- substr(peak_gene$annotation[i],23,37)
}
for(i in which(peak_gene$annotation_category =="Intron")){
  peak_gene$geneId[i] <- substr(peak_gene$annotation[i],25,39)
}

# download v98 gene ID and gene name pairs
library(biomaRt)
listEnsemblArchives()
mart <- useMart("ensembl", host="https://sep2019.archive.ensembl.org", dataset="hsapiens_gene_ensembl")
att <- listAttributes(mart)
selected_att <- c("ensembl_gene_id",
                  "external_gene_name")
att_mat <- getBM(attributes=selected_att, mart=mart,
                 filters = "external_gene_name",
                 values = "LGR5")
att_mat <- getBM(attributes=selected_att, mart=mart)
saveRDS(att_mat, file="/Volumes/pred/ihb-intestine-evo/Annotation/Ensembl/Human/Table_v98_ensembl_ID_and_gene_symbol.rds")

# combine two data frames by matching the column "geneId" in the dataframe "peak_gene" with the column "ensembl_gene_id" in the dataframe "att_mat"
library(tidyverse)
combined_df <- peak_gene %>%
  left_join(att_mat, by = c("geneId" = "ensembl_gene_id"))
saveRDS(combined_df, file="Data_frame_peak_gene_pairs.rds")


combined_rna <- combined
combined_rna[["ATAC"]] <- NULL
combined_rna[["peaks_merged"]] <- NULL
DefaultAssay(combined_rna) <- "RNA"
combined_rna <- NormalizeData(object = combined_rna, normalization.method = "LogNormalize", scale.factor = 1e4)
combined_rna <- FindVariableFeatures(combined_rna, selection.method = "vst", nfeatures = 3000)
combined_rna <- ScaleData(object = combined_rna, verbose = T)
combined_rna <- RunPCA(object = combined_rna, features = VariableFeatures(combined_rna), verbose = F, npcs = 50)
usefulPCs <- 1:50
combined_rna <- FindNeighbors(object = combined_rna, dims = usefulPCs)
combined_rna <- FindClusters(object = combined_rna, resolution = 1)
combined_rna <- RunUMAP(object = combined_rna, dims = usefulPCs)
SCpubr::do_DimPlot(combined_rna, reduction = "umap", split.by = "orig.ident")
SCpubr::do_FeaturePlot(combined_rna, reduction = "umap", feature="CCL20", split.by = "orig.ident")
SCpubr::do_FeaturePlot(combined_rna, reduction = "umap", feature="BEST4", order=T)
SCpubr::do_DimPlot(combined_rna, reduction = "umap", group.by = "RNA_snn_res.1")
saveRDS(combined_rna, file="Dat_filtered_colon_epithelium_infection_RNA_preprocessed.rds")

p1 <- SCpubr::do_DimPlot(combined_rna, split.by = "orig.ident", pt.size=6, font.size=100, border.size=1.5, legend.icon.size=40)
png("Plot_UMAP_RNA_colored_by_sample.png", width = 2000*3, height = 2000*2)
p1
dev.off()

p2 <- SCpubr::do_DimPlot(combined_rna, group.by = "RNA_snn_res.1", pt.size=6, font.size=100, border.size=1.5, legend.icon.size=40, label=T, label.size=20)+NoLegend()
p3 <- SCpubr::do_DimPlot(combined_rna, group.by = "condition", pt.size=6, font.size=100, border.size=1.5, legend.icon.size=40)
png("Plot_UMAP_RNA_colored_by_cluster.png", width = 2000, height = 2000*2)
p2/p3
dev.off()

p2 <- SCpubr::do_FeaturePlot(combined_rna, feature="CCL20", split.by = "orig.ident", pt.size=6, font.size=100, border.size=1.5)
png("Plot_UMAP_RNA_CCL20_split_by_sample.png", width = 2000*2, height = 2000*2)
p2
dev.off()


# ATAC analysis
DefaultAssay(combined) <- "peaks_merged"
# only keep peaks in standard chromosomes
peaks.keep <- seqnames(granges(combined)) %in% standardChromosomes(granges(combined))
combined <- combined[as.vector(peaks.keep), ]
combined[["ATAC"]] <- NULL

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 'q0')
det_rate <- rowMeans(combined@assays$peaks_merged@counts>0)
VariableFeatures(combined) <- names(det_rate)[which(det_rate>0.01)]
combined <- RunSVD(combined)
combined <- RunUMAP(combined, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:50)
combined <- FindClusters(object = combined, algorithm = 3, resolution=0.5)
saveRDS(combined, file="Dat_filtered_colon_epithelium_infection_ATAC_preprocessed.rds")
p1 <- SCpubr::do_DimPlot(combined, reduction="umap.atac", split.by = "orig.ident", pt.size=6, font.size=100, border.size=1.5, legend.icon.size=40)
png("Plot_UMAP_ATAC_colored_by_sample.png", width = 2000*3, height = 2000*2)
p1
dev.off()


p3 <- SCpubr::do_DimPlot(combined, reduction="umap.atac", group.by = "peaks_merged_snn_res.0.5", pt.size=6, border.size=1.5, label=T, label.size=20)+NoLegend()
png("Plot_UMAP_ATAC_colored_by_cluster.png", width = 2000, height = 2000)
p3
dev.off()

# coarse cell type annotation based on scRNA-seq data
combined_rna <- readRDS("Dat_filtered_colon_epithelium_infection_multiome_RNA_preprocessed.rds")
genes <- c("ALDOB", "CEACAM7", "LGR5", "OLFM4", "ASCL2","SLC12A2","CHGA", "CHGB", "ISL1", "FCGBP", "MUC2", "MEIS1","CA7","BEST4","MKI67")
p2 <- SCpubr::do_DotPlot(sample = combined_rna, 
             features = genes,
             scale=F,
             group.by="RNA_snn_res.1",
             cluster=T)
p2
ggsave(p2, filename="Plot_dotplot_adult_colon_chip_Salmonella_infection_cell_type_markers.pdf", width=5, height=7)

genes <- c("ALDOB", "CEACAM7", "LGR5", "OLFM4", "ASCL2","SLC12A2")
plotFeature(seu.obj=combined_rna, dr="umap", col.num=3, genes.to.plot=genes)

combined_rna$Coarse_cell_type <- "Colonocyte"
ct_anno <- list(
  "SC"=c(20,8),
  "Goblet"=c(6,10),
  "EEC"=23,
  "BEST4+ cell"=22,
  "TA"=c(14,19)
)
for(x in names(ct_anno)){
  combined_rna$Coarse_cell_type[combined_rna$RNA_snn_res.1%in%ct_anno[[x]]] <- x
}
saveRDS(combined_rna, file="Dat_filtered_colon_epithelium_infection_multiome_RNA_preprocessed.rds")
saveRDS(combined, file="Dat_filtered_colon_epithelium_infection_ATAC_preprocessed.rds")

combined$Coarse_cell_type <- combined_rna$Coarse_cell_type
p1 <- SCpubr::do_DimPlot(combined_rna, group.by = "Coarse_cell_type", pt.size=6, font.size=100, border.size=1.5, legend.icon.size=40, label=T, label.size=20)+NoLegend()
p2 <- SCpubr::do_DimPlot(combined, reduction="umap.atac", group.by = "Coarse_cell_type", pt.size=6, font.size=100, border.size=1.5, legend.icon.size=40, label=T, label.size=20)+NoLegend()
cells <- colnames(combined)[which(combined$peaks_merged_snn_res.0.5%in%c(8,9))]
p3 <- SCpubr::do_DimPlot(combined_rna, reduction="umap", cells.highlight=cells, pt.size=6, border.size=1.5)
p4 <- SCpubr::do_DimPlot(combined, reduction="umap.atac", cells.highlight=cells, pt.size=6, border.size=1.5)
png("Plot_UMAP_RNA_and_ATAC_cell_type.png", width = 2000*2, height = 2000*2)
(p1+p2)/(p3+p4)
dev.off()

# run chromVAR
# build GRN
# DAR identification
# ATAC-seq analysis
# there are colonocyte subset specific to Salmonella infection
# subset to the colonocytes
coloncyte_atac <- subset(combined, subset = Coarse_cell_type == "Colonocyte")
saveRDS(coloncyte_atac, file="Dat_filtered_colonocyte_infection_dataset_ATAC_preprocessed.rds")
SCpubr::do_DimPlot(coloncyte_atac, reduction="umap.atac")
# identify marker regions
de_res <- presto::wilcoxauc(coloncyte_atac, group_by="peaks_merged_snn_res.0.5", seurat_assay="peaks_merged")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
top_res <- sig_res %>% group_by(group) %>% top_n(200, wt=logFC)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_colonocyte_infection_dataset_cluster_marker_regions.rds")

# run DAR analysis per cell type
## for colonocyte, exclude C8 and C9
## run test per cell type
setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/")
library(Seurat)
library(Signac)
library(dplyr)
library(doParallel)
registerDoParallel(20)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")
combined <- readRDS("Dat_filtered_colon_epithelium_infection_ATAC_preprocessed.rds")
combined$Coarse_cell_type_per_condition <- paste(combined$Coarse_cell_type, combined$condition, sep=":")
# load scRNA-seq data
combined_rna <- readRDS("Dat_filtered_colon_epithelium_infection_multiome_RNA_preprocessed.rds")
combined_rna$Coarse_cell_type_per_condition <- paste(combined_rna$Coarse_cell_type, combined_rna$condition, sep=":")
combined_rna[["peaks_merged"]] <- combined[["peaks_merged"]]
cells <- colnames(combined)[which(!combined$peaks_merged_snn_res.0.5%in%c(8,9))]
seu_obj <- subset(combined, cells=cells)
# get per cell type detection rate
cell_type_det_rate <- getExpressedProp(seu_obj, feature.to.calc = "Coarse_cell_type_per_condition", assay.type="peaks_merged", colname.prefix =NULL)
saveRDS(cell_type_det_rate, file="Res_cell_type_detection_rate_per_condition_ex_C8_C9.rds")
library(doParallel)
ddp_by_cell_type_top_res <- list()
ddp_by_cell_type_sig_res <- list()
group <- sort(unique(seu_obj$condition))
for(x in sort(unique(seu_obj$Coarse_cell_type))){
  print(paste(x, "start"))
  
  seu_obj_x <- subset(seu_obj, subset = Coarse_cell_type == x)
  de_res <- presto::wilcoxauc(seu_obj_x, group_by="Coarse_cell_type_per_condition", seurat_assay="peaks_merged")
  de_res$pct_diff <- de_res$pct_in - de_res$pct_out
  sig_res <- de_res %>% filter(padj<0.05 & pct_diff>5 & logFC>0.01)
  top_res <- sig_res %>% group_by(group) %>% top_n(200, wt=logFC)
  deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
  saveRDS(deg_res, file=paste0("Res_",x,"_infection_related_DAR.rds"))

  ddp_by_cell_type_top_res[[x]] <- top_res
  ddp_by_cell_type_sig_res[[x]] <- sig_res
}
ddp_by_cell_type <- list("top_res"=do.call('rbind', ddp_by_cell_type_top_res), "sig_res"=do.call('rbind', ddp_by_cell_type_sig_res))
saveRDS(ddp_by_cell_type, file="Res_infection_related_DAR_by_cell_type.rds")

sig_res <- ddp_by_cell_type$sig_res
pdf("Plot_barplot_DAR_per_cell_type.pdf")
par(mar=c(10,5,5,5))
barplot(table(sig_res$group), las=2, ylab="# Region with higher accessibility")
dev.off()

combined_rna$peaks_merged_snn_res.0.5 <- combined$peaks_merged_snn_res.0.5
atac_coor <- Embeddings(combined, reduction="umap.atac")
combined_rna[["umap.atac"]] <- CreateDimReducObject(embeddings = atac_coor, key="ATACUMAP_", assay="peaks_merged")
SCpubr::do_DimPlot(combined_rna, reduction="umap.atac", group.by="peaks_merged_snn_res.0.5", pt.size=6, font.size=100, border.size=1.5, legend.icon.size=40, label=T, label.size=20)+NoLegend()
saveRDS(combined_rna, file="Dat_filtered_colon_epithelium_infection_RNA_ATAC_preprocessed_ctByCondition.rds")

# generate coverage plot
setwd("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/coverage_plot")
library(Seurat)
library(Signac)
library(IRanges)
library(GenomicRanges)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
# load evolutionary signature list
evo_list <- readRDS("/projects/site/pred/ihb-intestine-evo/evo_signature/summary_Feb_2024/Dat_SNC_GWAS_HAR_PS_HAQER.rds")

# as a sanity check, check whether regions associated with selected genes also show differential accessibility
genes <- c("CCL20", "TNFAIP3", "IL1B", "IRAK2")
peak_anno <- readRDS("Data_frame_peak_gene_pairs.rds")
peaks <- peak_anno$name[which(peak_anno$external_gene_name%in%genes & peak_anno$name%in%sig_res$feature)]
peaks <- "chr2-227812975-227814642" # promoter of CCL20
colors <- readRDS()

source("/home/yuq22/ihb-intestine-evo/common_script/plotting/Script_plotting_functions.R")
p_list <- coverage_plot_with_evo_sig(seu_obj=combined_rna, 
                                     peak_assay="peaks_merged", 
                                     rna_assay="RNA", 
                                     group="Coarse_cell_type_per_condition", 
                                     colors.use=NULL, 
                                     peaks=peaks, 
                                     evo_list_path="/projects/site/pred/ihb-intestine-evo/evo_signature/summary_Feb_2024/Dat_SNC_GWAS_HAR_PS_HAQER.rds", 
                                     peak_anno_path="Data_frame_peak_gene_pairs.rds", 
                                     peak_anno_peak_col="name",
                                     peak_anno_gene_col="external_gene_name",
                                     extend_up=10000,
                                     extend_down=10000,
                                     window_size=500,
                                     color="#303030",
                                     gene_mode="ChIPseeker", 
                                     g_vec=NULL,
                                     do_plot_combined=TRUE,
                                     do_plot_individual=TRUE,
                                     combined_plot_name=NULL,
                                     plot_suffix="-2")


# specifically focusing on C8 and C9 which are the infection-specific clusters
library(Pando)
data('motif2tf')
peak_gene <- readRDS("Data_frame_peak_gene_pairs.rds")
# generate coverage plot for top 10 markers of C8 and C9
deg_res <- readRDS("Res_colonocyte_infection_dataset_cluster_marker_regions.rds")
sig_res <- deg_res$sig_res
cl=8
df <- sig_res[which(sig_res$group==cl), ]
peaks <- df$feature[order(df$logFC, decreasing=T)[1:12]]
p_list <- coverage_plot_with_evo_sig(seu_obj=combined_rna, 
                                     peak_assay="peaks_merged", 
                                     rna_assay="RNA", 
                                     group="peaks_merged_snn_res.0.5", 
                                     colors.use=NULL, 
                                     peaks=peaks, 
                                     evo_list_path="/projects/site/pred/ihb-intestine-evo/evo_signature/summary_Feb_2024/Dat_SNC_GWAS_HAR_PS_HAQER.rds", 
                                     peak_anno_path="/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/Data_frame_peak_gene_pairs.rds", 
                                     peak_anno_peak_col="name",
                                     peak_anno_gene_col="external_gene_name",
                                     extend_up=10000,
                                     extend_down=10000,
                                     window_size=500,
                                     color="#303030",
                                     gene_mode="ChIPseeker", 
                                     g_vec=NULL,
                                     do_plot_combined=TRUE,
                                     do_plot_individual=TRUE,
                                     combined_plot_name=NULL,
                                     plot_suffix="-2")

p <- plot_grid(plotlist = p_list, align = "hv", ncol = 4)
nrow=3
ncol=4
plot_name <- "Plot_coveragePlot_for_hg38_selected_region.pdf"
ggsave(p, filename=plot_name, width=10*ncol, height=15*nrow)
# get condition proportion per scATAC-seq cell cluster
n1 <- sapply(sort(unique(combined_rna$peaks_merged_snn_res.0.5)), function(x){
  sapply(sort(unique(combined_rna$condition)), function(y){
    sum(combined_rna$condition==y & combined_rna$peaks_merged_snn_res.0.5==x)
  })
})
colnames(n1) <- paste0("C",sort(unique(combined_rna$peaks_merged_snn_res.0.5)))
