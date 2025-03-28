# this script is to integrate public adult intestine epithelial scRNA-seq data with the cell type annotation per paper
setwd("/home/yuq22/ihb-intestine-evo/adult_primate_intestine_atlas/primate_duo/analysis/public_human/intestine_multi_region_scRNA-seq_atlas/")
library(Seurat)
library(dplyr)
# dataset-1
#library(SeuratDisk)

# to get to know how to get dataset-1 to -3 data, please refer to the following script: "/home/yuq22/ihb-intestine-evo/adult_primate_intestine_atlas/primate_duo/analysis/public_human/Script_get_public_adult_human_duodenum_epithelium.R"
# dataset-1
#ds1 <- readRDS("../Dat_gut_cell_atlas_full_dataset_raw_count.rds")
### subset to adult gut cells (not including the mesenteric lymph nodes (mLN))
#ds1_adult <- subset(ds1, Age_group=="Adult")
#saveRDS(ds1_adult, file="Dat_gut_cell_atlas_healthy_adult_gut_raw_count.rds")
ds1 <- readRDS("Dat_gut_cell_atlas_healthy_adult_gut_raw_count.rds")
ds1$paper <- "gut_cell_atlas"
# extract epithelial cells
ds1_epi <- subset(ds1, category=="Epithelial")
ds1_epi$Cell_type <- ds1_epi$Integrated_05
saveRDS(ds1_epi, file="Dat_gut_cell_atlas_healthy_adult_gut_epi_raw_count.rds")

## dataset-2, Yu-2021, hg38-based, so not yet have the cell type annotation
ds2 <- readRDS("../Dat_guttuber_adult_duo_raw_count.rds")
ds2$paper <- "guttuber"
#seu_obj <- readRDS("~/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/used_seurat_objects/Res_tHIO_fetal_and_adult_duodenum_epi_integrated_with_CSS.rds") 
seu_obj <- readRDS("/Volumes/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/used_seurat_objects/Res_tHIO_fetal_and_adult_duodenum_epi_integrated_with_CSS.rds") 
# extract the adult cells 
# add cell type annotation to the dataset
hg19 <- subset(seu_obj, Stage=="adult")
ct_vec <- setNames(hg19$Cell_type, colnames(hg19))
bc <- do.call('rbind', strsplit(names(ct_vec), split="_"))[,2]
bc2 <- paste0(bc, "-1_", hg19$Individual)
names(ct_vec) <- bc2
shared_cells <- intersect(names(ct_vec), colnames(ds2))
ds2_anno <- subset(ds2, cells=shared_cells)
ds2_anno$Cell_type <- ct_vec[colnames(ds2_anno)]
saveRDS(ds2_anno, file="Dat_guttuber_adult_duo_epi_raw_count_with_anno.rds")

## dataset-3
ds3 <- readRDS("../Dat_Burclaff_et_all_adult_human_multi_region_full_dataset_raw_counts.rds")
ds3$paper <- "Burclaff"
ds3$Cell_type <- ds3$type
seu_obj_list <- list(ds1_epi, ds2_anno, ds3)
combined <- merge(x=seu_obj_list[[1]],
                  y=seu_obj_list[-1])
combined <- NormalizeData(object = combined, normalization.method = "LogNormalize", scale.factor = 1e4)
saveRDS(combined, file="Dat_combined_adult_human_epi_normalized_data.rds")

# identify hvg per paper
papers <- sort(unique(combined$paper))
hvg_idx <- matrix(FALSE, nrow=nrow(combined), ncol=length(papers))
colnames(hvg_idx) <- papers
rownames(hvg_idx) <- rownames(combined)
for(x in papers){
  seu_obj <- subset(combined, paper==x)
  seu_obj <- FindVariableFeatures(object = seu_obj, selection.method = "vst", nfeatures = 3000)
  hvg <- VariableFeatures(seu_obj)
  hvg_idx[hvg,x] <- TRUE
}
freq <- rowSums(hvg_idx)
rownames(hvg_idx)[which(rowSums(hvg_idx)>1)] # 2297 hvg shred by at least 2 papers

seu_obj_list <- SplitObject(combined, split.by = "paper")
selected_hvg <- SelectIntegrationFeatures(
  seu_obj_list,
  nfeatures = 3000,
  fvf.nfeatures = 3000
)
combined <- FindVariableFeatures(object = combined, selection.method = "vst", nfeatures = 3000)
VariableFeatures(combined) <- selected_hvg
combined <- ScaleData(object = combined, verbose = T)
combined <- RunPCA(object = combined, features = VariableFeatures(combined), verbose = F, npcs = 50)
saveRDS(combined, file="Dat_combined_adult_human_epi_normalized_data_with_PCA.rds")
usefulPCs <- 1:20
combined <- FindNeighbors(object = combined, dims = usefulPCs)
combined <- FindClusters(object = combined, resolution = 1)
combined <- RunUMAP(object = combined, dims = usefulPCs)
combined$RNA_PCA_snn_res.1 <- combined$RNA_snn_res.1
saveRDS(combined, file="Res_combined_adult_human_epi_no_batch_correction.rds")

p1 <- SCpubr::do_DimPlot(combined, reduction = "umap", group.by = "paper", pt.size=4, border.size = 2, label=T, label.size=3)
png("Plot_UMAP_scRNA_seq_adult_human_epi_no_batch_correction_colored_by_paper.png", width=2000, height=2000)
p1
dev.off()

library(ggplot2)
library(cowplot)
plot_list <- list()
for(x in sort(unique(combined$paper))){
  seu_obj <- subset(combined, paper==x)
  p <- SCpubr::do_DimPlot(sample=seu_obj, reduction = "umap", group.by = "Cell_type", pt.size=4, border.size = 2, label=T, label.size=6)
  plot_list[[x]] <- p
}
p1 <- plot_grid(plotlist = plot_list, nrow = 1, ncol = 3)
png("Plot_UMAP_scRNA_seq_adult_human_epi_no_batch_correction.png", width=2000*3, height=2000)
p1
dev.off()
combined$orig.ident[which(combined$paper=="gut_cell_atlas")] <- combined$sample.name[which(combined$paper=="gut_cell_atlas")]
saveRDS(combined, file="Res_combined_adult_human_epi_no_batch_correction.rds")

# harmony integration - stratify by orig.ident
library(harmony)
DefaultAssay(combined) <- "RNA"

aa <- RunHarmony(combined,
                 group.by.vars = "orig.ident",
                 dims.use = 1:20)
aa <- RunUMAP(aa,
                 reduction = "harmony",
                 dims = 1:20,
                 reduction.key = "UMAPHARMONY_",
                 reduction.name="umap_harmony")
saveRDS(aa, file="Res_combined_adult_human_epi_with_harmony.rds")

library(ggplot2)
library(cowplot)
plot_list <- list()
p1 <- SCpubr::do_DimPlot(aa, reduction = "umap_harmony", group.by = "paper", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10)
plot_list[[1]] <- p1

for(x in sort(unique(aa$paper))){
  seu_obj <- subset(aa, paper==x)
  p <- SCpubr::do_DimPlot(sample=seu_obj, reduction = "umap_harmony", group.by = "Cell_type", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10)
  plot_list[[x]] <- p
}
p1 <- plot_grid(plotlist = plot_list, nrow = 2, ncol = 2)
png("Plot_UMAP_scRNA_seq_adult_human_epi_with_harmony.png", width=2000*2, height=2000*2)
p1
dev.off()


# CSS integration
css <- simspec::cluster_sim_spectrum(combined, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
css <- RunUMAP(css, reduction = "css", dims = 1:ncol(css@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
saveRDS(css, file="Res_combined_adult_human_epi_with_CSS.rds")
css <- FindNeighbors(object = css, reduction = "css", dims = 1:ncol(css@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 1)
css[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- css[[paste0("RNA_snn_res.", 0.2*5)]]
saveRDS(css, file="Res_combined_adult_human_epi_with_CSS.rds")

SCpubr::do_DimPlot(css, reduction = "umap_css", group.by = "RNA_CSS_snn_res.1", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10)

plot_list <- list()
p1 <- SCpubr::do_DimPlot(css, reduction = "umap_css", group.by = "paper", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10)
plot_list[[1]] <- p1

for(x in sort(unique(css$paper))){
  seu_obj <- subset(css, paper==x)
  p <- SCpubr::do_DimPlot(sample=seu_obj, reduction = "umap_css", group.by = "Cell_type", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10)
  plot_list[[x]] <- p
}
p1 <- plot_grid(plotlist = plot_list, nrow = 2, ncol = 2)
png("Plot_UMAP_scRNA_seq_adult_human_epi_with_CSS.png", width=2000*2, height=2000*2)
p1
dev.off()

# identify cluster markers of epi ref
de_res <- presto::wilcoxauc(css, group_by="RNA_CSS_snn_res.1")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(auc > 0.6 & padj<0.05 & pct_diff>0.1 & pct_in>0.1 & logFC>0.1) 
top_res <- sig_res %>% group_by(group) %>% top_n(10, logFC)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_CSS_integrated_adult_human_epi_cluster_markers.rds")

selected_groups <- 16
genes <- unique(c("DEFA5", "DEFA6", "SEMA3F", "KIF19", "EPHB3", "NOTUM", top_res$feature[which(top_res$group %in% selected_groups)]))
genes <- c("DEFA5", "DEFA6", "SEMA3F", "KIF19", "EPHB3", "NOTUM")

plotFeature(seu.obj=css, genes.to.plot=genes, dr="umap_css", plot.name="Plot_UMAP_CSS_adult_epi_ref_paneth_c16_cluster_markers.png", do.plot=TRUE, col.num=5)
# c16, c12 - paneth cell subtypes
# sample info
table(css$Region[which(css$RNA_snn_res.1==12)])

# unify the region info of reference data
# for Burclaff et al., the regions are:
#duodenum, jejunum, ileum, and ascending, transverse, and descending colon
# for gut cell atlas, the regions are:
#duodenum (two locations, DUO1 and DUO2, which were pooled in the analysis), jejunum (JEJ), ileum (two locations, ILE1 and ILE2, which were pooled in the analysis), 
#appendix (APD), caecum (CAE), ascending colon (ACL), transverse colon (TCL), descending colon (DCL), sigmoid colon (SCL), rectum (REC) and mesenteric lymph nodes (mLN)
css$Unified_region_code <- css$Region.code
css$Unified_region_code <- toupper(css$Unified_region_code)
css$Unified_region_code[which(css$Unified_region_code%in%c("TC", "AC", "DC"))] <- paste0(css$Unified_region_code[which(css$Unified_region_code%in%c("TC", "AC", "DC"))], "L")
css$Unified_region_code[grep("ILE", css$Unified_region_code)] <- "ILE"
css$Region[which(css$Unified_region_code%in%c("DUO", "JEJ", "ILE"))] <- "SmallInt"
css$Region[which(css$Unified_region_code%in%c("ACL", "TCL", "DCL"))] <- "LargeInt"
saveRDS(css, file="Res_combined_adult_human_epi_with_CSS.rds")
p1 <- SCpubr::do_DimPlot(sample=css, reduction = "umap_css", group.by = "Region", pt.size=2, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10)
p1

# unify the cell type annotation 
ds_burclaff <- subset(css, paper=="Burclaff")
p1 <- SCpubr::do_DimPlot(sample=ds_burclaff , reduction = "umap_css", group.by = "lineage", pt.size=2, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10)
p1 
png("Plot_UMAP_CSS_adult_epi_ref_Burclaff_lineage.png", width=2000, height=2000)
p1  
dev.off()

cell_type_anno <- list(
  "Absorptive" = c("absorptive", "Enterocyte",  "Enterocyte_precursor", "Colonocyte"),
  "Goblet cell" = c("goblet", "Goblet", "BEST2+ Goblet cell", "Goblet cell"),
  "BEST4+ cell" = c("BEST4+", "BEST4+ epithelial", "M_cell"), 
  "EEC" = c("Progenitor (NEUROG3+)", "D cells (SST+)", "EC cells (TAC1+)", "EEC", "EECs", "Enteroendocrine", "I cells (CCK+)","K cells (GIP+)", "L cells (PYY+)","M/X cells (MLN/GHRL+)", "N cells (NTS+)" ),
  "FAE" = "FAE",
  "Stem cell" = c("ISC", "Stem cells", "Stem_cell"),
  "M cell"=c("Microfold cell", "SPIB-/GP2+/MUC1_M-cell"),
  "Paneth cell"= c("paneth", "Paneth"),
  "Secretory prog."="secretory_prog",
  "TA"="TA",
  "Tuft cell"=c("tuft", "Tuft")
)
css$Unified_coarse_cell_type <- NA
for(x in names(cell_type_anno)){
  css$Unified_coarse_cell_type[which(css$Cell_type %in% cell_type_anno[[x]])] <- x
}
p1 <- SCpubr::do_DimPlot(sample=css, reduction = "umap_css", group.by="Unified_coarse_cell_type", split.by = "Unified_coarse_cell_type", pt.size=1, border.size = 1, label=F, font.size=10)
p1
# for aborptive cells and BEST4+ cell, further separate them based on the region
css$Unified_coarse_cell_type[which(css$Unified_coarse_cell_type=="BEST4+ cell" & css$Region=="SmallInt")] <- "SmallInt BEST4+"
css$Unified_coarse_cell_type[which(css$Unified_coarse_cell_type=="BEST4+ cell" & css$Region!="SmallInt")] <- "non-SmallInt BEST4+"
css$Unified_coarse_cell_type[which((css$Cell_type%in%c("Enterocyte",  "Enterocyte_precursor")) | (css$Cell_type=="absorptive" & css$Region=="SmallInt"))] <- "Enterocyte"
css$Unified_coarse_cell_type[which((css$Cell_type%in%c("Colonocyte")) | (css$Cell_type=="absorptive" & css$Region!="SmallInt"))] <- "Colonocyte"
css$Unified_coarse_cell_type[which(css$Unified_coarse_cell_type=="Stem cell" & css$Region=="SmallInt")] <- "SmallInt stem cell"
css$Unified_coarse_cell_type[which(css$Unified_coarse_cell_type=="Stem cell" & css$Region!="SmallInt")] <- "non-SmallInt stem cell"
css$Unified_coarse_cell_type[which(css$Unified_coarse_cell_type=="TA" & css$Region=="SmallInt")] <- "SmallInt TA"
css$Unified_coarse_cell_type[which(css$Unified_coarse_cell_type=="TA" & css$Region!="SmallInt")] <- "non-SmallInt TA"
css$Unified_coarse_cell_type[which(css$Unified_coarse_cell_type=="Secretory prog." & css$Region=="SmallInt")] <- "SmallInt secretory prog."
css$Unified_coarse_cell_type[which(css$Unified_coarse_cell_type=="Secretory prog." & css$Region!="SmallInt")] <- "non-SmallInt secretory prog."
saveRDS(css, file="Res_combined_adult_human_epi_with_CSS_ct_anno_updated.rds")

# load Yu-2021 data and refine the cell type annotation
setwd("/home/yuq22/Bacteria_TRM/intestine_multi_region_scRNA-seq_atlas/update_yu_adult_data")
seu_obj <- readRDS("/projects/site/pred/ihb-intestine-evo/USERS/Qianhui_Yu/Endoderm/used_seurat_objects/Res_tHIO_fetal_and_adult_duodenum_epi_integrated_with_CSS.rds")
# subset to adult cells
adult_obj <- subset(seu_obj, Stage=="adult")
SCpubr::do_DimPlot(sample=adult_obj, reduction = "umap_css", group.by = "Cell_type", pt.size=2, border.size = 2, label=F, label.size=10, font.size=30, legend.icon.size=10)
genes <-  c("PGC", "TFF2", "MUC1", "BEST4", "CFTR", "KRT20", "ADGRR4")
plotFeature(seu.obj=adult_obj, dr="umap_css", genes.to.plot=genes, col.num=6, plot.name="Plot_UMAPCSS_adult_M_and_BEST4_cell_marker.png", nCols=beach.col, per.plot.size=2000, cex=5)
genes <- c("DEFA5", "DEFA6", "RGS13", "AVIL", "ALOX5")
plotFeature(seu.obj=adult_obj, dr="umap_css", genes.to.plot=genes, col.num=6, plot.name="Plot_UMAPCSS_adult_tuft_and_paneth_cell_marker.png", nCols=beach.col, per.plot.size=2000, cex=5)
# "SPIB-/GP2+/MUC1_M-cell" is M cells, previously annotated M_cells are BEST4+ cells
# annotated tuft cells and paneth cells
anchor_genes <- c("DEFA5","DEFA6")
seu_obj <- adult_obj
seu_obj <- FindVariableFeatures(object = seu_obj, selection.method = "vst", nfeatures = 3000)
ref_pattern <- colSums(t(scale(t(as.matrix(seu_obj@assays$RNA@data[anchor_genes,])))))
cor_vec <- cor(t(as.matrix(seu_obj@assays$RNA@data[VariableFeatures(seu_obj),])), ref_pattern)[,1]
genes <- names(sort(cor_vec, decreasing = T)[1:5])
plotFeature(seu.obj=seu_obj, dr="umap_css", genes.to.plot=genes, col.num=5, do.plot=T, plot.name="Plot_UMAPCSS_paneth_cell_marker_feature_plot.png", nCols=beach.col, per.plot.size=2000, cex=10)
seu_obj$expr_based_paneth_score <- colSums(t(scale(t(as.matrix(seu_obj@assays$RNA@data[genes,])))))
SCpubr::do_FeaturePlot(seu_obj, reduction="umap_css", feature="expr_based_paneth_score", pt.size=5, order=T)
score2 <- seu_obj$expr_based_paneth_score
score2 <- score2[which(score2>0)]
plot(sort(score2), seq(length(score2)), pch=16)
cutoff2 <- quantile(score2, 0.975)
abline(v=cutoff2, col="red")
paneth_cells <- colnames(seu_obj)[which(seu_obj$expr_based_paneth_score>cutoff2)]
SCpubr::do_DimPlot(seu_obj, reduction = "umap_css", cells.highlight=paneth_cells)

anchor_genes <- c("RGS13", "AVIL", "ALOX5")
ref_pattern <- colSums(t(scale(t(as.matrix(seu_obj@assays$RNA@data[anchor_genes,])))))
cor_vec <- cor(t(as.matrix(seu_obj@assays$RNA@data[VariableFeatures(seu_obj),])), ref_pattern)[,1]
genes <- names(sort(cor_vec, decreasing = T)[1:10])
plotFeature(seu.obj=seu_obj, dr="umap_css", genes.to.plot=genes, col.num=5, do.plot=T, plot.name="Plot_UMAPCSS_tuft_cell_marker_feature_plot.png", nCols=beach.col, per.plot.size=2000, cex=10)
seu_obj$expr_based_tuft_score <- colSums(t(scale(t(as.matrix(seu_obj@assays$RNA@data[genes,])))))
SCpubr::do_FeaturePlot(seu_obj, reduction="umap_css", feature="expr_based_tuft_score", pt.size=5, order=T)
score2 <- seu_obj$expr_based_tuft_score
score2 <- score2[which(score2>0)]
plot(sort(score2), seq(length(score2)), pch=16)
cutoff2 <- quantile(score2, 0.95)
abline(v=cutoff2, col="red")
tuft_cells <- colnames(seu_obj)[which(seu_obj$expr_based_tuft_score>cutoff2)]
SCpubr::do_DimPlot(seu_obj, reduction = "umap_css", cells.highlight=tuft_cells, pt.size=3)

cell <- intersect(paneth_cells, tuft_cells)
sum(seu_obj$expr_based_tuft_score>seu_obj@meta.data[cell, "expr_based_tuft_score"])
sum(seu_obj$expr_based_paneth_score>seu_obj@meta.data[cell, "expr_based_paneth_score"])
# tuft cell score rank higher than paneth cell score, so assign to tuft cells
seu_obj$updated_cell_type <- seu_obj$Cell_type
seu_obj@meta.data[paneth_cells, "updated_cell_type"] <- "Paneth cell"
seu_obj@meta.data[tuft_cells, "updated_cell_type"] <- "Tuft cell"
seu_obj$updated_cell_type[which(seu_obj$updated_cell_type=="SPIB-/GP2+/MUC1_M-cell")] <- "M cell"
seu_obj$updated_cell_type[which(seu_obj$updated_cell_type=="M_cell")] <- "BEST4+ cell"
table(seu_obj$updated_cell_type)
name_pair <- list(
"Enterocyte" = c("Enterocyte", "Enterocyte_precursor"),
"SmallInt stem cell" = "Stem_cell",
"SmallInt BEST4+" ="BEST4+ cell",
"SmallInt goblet cell" = "Goblet",
"EEC"="Enteroendocrine"
)
for(x in names(name_pair)){
  seu_obj$updated_cell_type[which(seu_obj$updated_cell_type %in% name_pair[[x]])] <- x
}
saveRDS(seu_obj, file="Res_Yu_et_al_adult_duodenum_epi_with_updated_cell_type.rds")
SCpubr::do_DimPlot(seu_obj, reduction="umap_css", group.by="updated_cell_type", pt.size=2, border.size = 2, label=F, label.size=10, font.size=30, legend.icon.size=10)
# update tissue data annotation
ct_vec <- setNames(seu_obj$updated_cell_type, colnames(seu_obj))
bc <- do.call('rbind', strsplit(names(ct_vec), split="_"))[,2]
bc2 <- paste0(bc, "-1_", seu_obj$Individual)
names(ct_vec) <- bc2
shared_cells <- intersect(names(ct_vec), colnames(tissue_data))
tissue_data@meta.data[shared_cells, "Unified_coarse_cell_type"] <- ct_vec[shared_cells]
saveRDS(tissue_data, file="/home/yuq22/Bacteria_TRM/used_object/Res_combined_adult_human_epi_with_CSS_ct_anno_updated.rds")

p1 <- SCpubr::do_DimPlot(sample=css, reduction = "umap_css", group.by="Unified_coarse_cell_type", split.by = "Unified_coarse_cell_type", pt.size=1, border.size = 1, label=F, font.size=10)+NoLegend()
png("Plot_UMAP_CSS_adult_epi_unified_coarse_cell_type_annotation.png", width=2000, height=2000)
p1  
dev.off()

setwd("/home/yuq22/ihb-intestine-evo/adult_primate_intestine_atlas/primate_duo/analysis/public_human/intestine_multi_region_scRNA-seq_atlas/")
css <- readRDS("Res_combined_adult_human_epi_with_CSS_ct_anno_updated.rds")
saveRDS(css, file="/home/yuq22/Bacteria_TRM/used_object/Res_combined_adult_human_epi_with_CSS_ct_anno_updated.rds")
# further classify SI and non-SI goblet cells
css$Unified_coarse_cell_type[which(css$Unified_coarse_cell_type=="Goblet cell" & css$Region=="SmallInt")] <- "SmallInt goblet cell"
css$Unified_coarse_cell_type[which(css$Unified_coarse_cell_type=="Goblet cell" & css$Region!="SmallInt")] <- "non-SmallInt goblet cell"
saveRDS(css, file="Res_combined_adult_human_epi_with_CSS_ct_anno_updated.rds")
p1 <- SCpubr::do_DimPlot(css, reduction = "umap_css", group.by = "Unified_region_code", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=12)
p2 <- SCpubr::do_DimPlot(css, reduction = "umap_css", group.by = "Unified_coarse_cell_type", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=12)
p3 <- SCpubr::do_DimPlot(css, reduction = "umap_css", group.by = "RNA_CSS_snn_res.1", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=12)
p4 <- SCpubr::do_DimPlot(css, reduction = "umap_css", group.by = "paper", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=12)
png("Plot_UMAP_CSS_adult_intestine_tissue_epi.png", width=2000*2, height=2000*2)
(p1+p2)/(p3+p4)
dev.off()

# trajectory reconstruction
# extract the SI continuum
setwd("/home/yuq22/Bacteria_TRM/intestine_multi_region_scRNA-seq_atlas/SI_subset")
selected_cell_types <- c("SmallInt TA", "Enterocyte", "Paneth cell", "SmallInt stem cell", "SmallInt BEST4+",  "SmallInt secretory prog.", "FAE", "SmallInt goblet cell")
si <- subset(css, cells=colnames(css)[which(css$Unified_coarse_cell_type %in% selected_cell_types)])
si <- simspec::cluster_sim_spectrum(si, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
si <- RunUMAP(si, reduction = "css", dims = 1:ncol(si@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
saveRDS(si, file="Res_selected_SI_cell_types_with_CSS.rds")
p1 <- SCpubr::do_DimPlot(si, reduction = "umap_css", group.by = "Unified_coarse_cell_type", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=12)
p2 <- SCpubr::do_DimPlot(si, reduction = "umap_css", group.by = "Unified_region_code", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=12)
png("Plot_UMAP_CSS_adult_intestine_SI_non_Tuft_cell_and_EEC.png", width=2000*2, height=2000)
p1+p2
dev.off()

# trajectory reconstruction


# extract the non-SI continuum
setwd("/home/yuq22/Bacteria_TRM/intestine_multi_region_scRNA-seq_atlas/non_SI_subset")
selected_cell_types <- c("non-SmallInt TA", "non-SmallInt BEST4+",  "Colonocyte",  "non-SmallInt goblet cell",  "non-SmallInt stem cell", "SPIB-/GP2+/MUC1_M-cell",  "non-SmallInt secretory prog.")  
li <- subset(css, cells=colnames(css)[which(css$Unified_coarse_cell_type %in% selected_cell_types)])
li <- simspec::cluster_sim_spectrum(li, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
li <- RunUMAP(li, reduction = "css", dims = 1:ncol(li@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
saveRDS(li, file="Res_selected_non-SI_cell_types_with_CSS.rds")
p1 <- SCpubr::do_DimPlot(li, reduction = "umap_css", group.by = "Unified_coarse_cell_type", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=12)
p2 <- SCpubr::do_DimPlot(li, reduction = "umap_css", group.by = "Unified_region_code", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=12)
png("Plot_UMAP_CSS_adult_intestine_non-SI_non_Tuft_cell_and_EEC.png", width=2000*2, height=2000)
p1+p2
dev.off()




# subtypes of Paneth cells
# follicle associated epithelium (FAE) cells
Idents(css) <- setNames(css$Cell_type, colnames(css))
selected_cell_types <- c("Paneth", "paneth")
cells <- colnames(css)[which(css$Cell_type %in% selected_cell_types & css$paper=="Burclaff")]
SCpubr::do_DimPlot(css, reduction = "umap_css", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10, idents.highlight=selected_cell_types)
SCpubr::do_DimPlot(css, reduction = "umap_css", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10, cells.highlight=cells)
p1 <- SCpubr::do_DimPlot(css, reduction = "umap_css", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10, idents.highlight="Paneth")
p2 <- SCpubr::do_DimPlot(css, reduction = "umap_css", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10, idents.highlight="paneth")
p1+p2
# load adult reference data and refine the cell type annotation
css <- readRDS("Res_combined_adult_human_epi_with_CSS.rds")
p3 <- SCpubr::do_DimPlot(sample=css, reduction = "umap_css", group.by = "Cell_type", pt.size=2, border.size = 2, label=T, label.size=3, font.size=5, legend.icon.size=3)
p4 <- SCpubr::do_DimPlot(sample=css, reduction = "umap_css", group.by = "RNA_CSS_snn_res.1", pt.size=2, border.size = 2, label=T, label.size=3, font.size=5, legend.icon.size=3)
(p1+p2)/(p3+p4)
sort(table(css$RNA_snn_res.1[which(css$Cell_type %in% c("Paneth", "paneth"))]))


# extract paneth cells and perform subclustering
cells <- colnames(css)[which(css$Cell_type %in% selected_cell_types)]
paneth_obj <- subset(css, cells=cells)
paneth_obj <- FindVariableFeatures(object = paneth_obj, selection.method = "vst", nfeatures = 3000)
paneth_obj <- ScaleData(object = paneth_obj, verbose = T)
paneth_obj <- RunPCA(object = paneth_obj, features = VariableFeatures(paneth_obj), verbose = F, npcs = 20)
paneth_obj <- simspec::cluster_sim_spectrum(paneth_obj, label_tag = "paper", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
paneth_obj <- RunUMAP(paneth_obj, reduction = "css", dims = 1:ncol(paneth_obj@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
paneth_obj <- FindNeighbors(object = paneth_obj, reduction = "css", dims = 1:ncol(paneth_obj@reductions$css@cell.embeddings), force.recalc = T) %>%
  FindClusters(resolution = 0.4)
paneth_obj[[paste0("RNA_CSS_snn_res.", 0.2*2)]] <- paneth_obj[[paste0("RNA_snn_res.", 0.2*2)]]
saveRDS(paneth_obj, file="Res_combined_adult_human_paneth_cells_with_CSS.rds")
p1 <- SCpubr::do_DimPlot(sample=paneth_obj, reduction = "umap_css", group.by = "RNA_CSS_snn_res.0.4", pt.size=2, border.size = 2, label=T, label.size=8, font.size=15, legend.icon.size=3)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")
plotFeature(seu.obj=paneth_obj, genes.to.plot=c("DEFA5", "DEFA6", "SEMA3F", "KIF19", "EPHB3", "NOTUM"), dr="umap_css", plot.name="Plot_UMAP_CSS_paneth_cells_selected_gene_feature_plot.png")
genes <- c("EPCAM", "PTPRC", "ALDOB")
plotFeature(seu.obj=paneth_obj, genes.to.plot=genes, dr="umap_css", col.num=3, do.plot=FALSE)
SCpubr::do_FeaturePlot(sample=paneth_obj, reduction = "umap_css", features = c("DEFA5", "DEFA6", "SEMA3F", "KIF19", "EPHB3", "NOTUM"), pt.size=2, order=TRUE)
# identify cluster markers for paneth cells
de_res <- presto::wilcoxauc(paneth_obj, group.by="RNA_CSS_snn_res.0.4")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
library(dplyr)
sig_res <- de_res %>% filter(auc > 0.6 & padj<0.05 & pct_diff>0.1 & pct_in>0.1 & logFC>0.1) 
top_res <- sig_res %>% group_by(group) %>% top_n(10, logFC)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_paneth_cell_subcluster_markers.rds")
# cluster 8, 7 and 5 are depleted of DEFA5 and DEFA6, plot their positive markers
selected_groups <- c(0,1,2,3,6)
genes <- unique(top_res$feature[which(top_res$group %in% selected_groups)])
plotFeature(seu.obj=paneth_obj, genes.to.plot=genes, dr="umap_css", plot.name="Plot_UMAP_CSS_paneth_cells_de_novo_markers_of_DEFA_pos_big_clusters.png")

cells <- colnames(css)[which(css$Cell_type=="Paneth" & css$RNA_snn_res.1==16)]
SCpubr::do_DimPlot(paneth_obj, reduction = "umap_css", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10, cells.highlight=cells)

# try CSS stratified on region per individual # for some reason, with the CSS integration, the guttuber data is not well integrated
#setwd("/home/yuq22/Bacteria_TRM/intestine_multi_region_scRNA-seq_atlas/css_region_by_individual")
#library(Seurat)
#library(dplyr)
#library(ggplot2)
#
## Load the data
#combined <- readRDS("/home/yuq22/ihb-intestine-evo/adult_primate_intestine_atlas/primate_duo/analysis/public_human/intestine_multi_region_scRNA-seq_atlas/Res_combined_adult_human_epi_with_CSS_ct_anno_updated.rds")
#combined$Region_per_individual <- combined$orig.ident
#combined$Region_per_individual[which(combined$paper=="gut_cell_atlas")] <- paste(combined$Sample.name[which(combined$paper=="gut_cell_atlas")], combined$Unified_region_code[which(combined$paper=="gut_cell_atlas")], sep="_")
#combined$Region_per_individual[which(combined$paper=="Burclaff")] <- combined$orig.ident[which(combined$paper=="Burclaff")]
#combined$Region_per_individual[which(combined$paper=="guttuber")] <- paste0(combined$Region_per_individual[which(combined$paper=="guttuber")], "_DUO")
#saveRDS(combined, file="Res_combined_adult_enteroid_and_tissue_data_before_integration.rds")
#
#sample_size <- sort(table(combined$Region_per_individual))
## filter out samples with less than 30 cells (requried for CCA integration)
#selected_samples <- names(sample_size)[which(sample_size>=30)]
#combined_sub <- subset(combined, cells=colnames(combined)[which(combined$Region_per_individual%in%selected_samples)])
#
#bk <- combined 
#combined <- combined_sub
#seu_obj_list <- SplitObject(combined, split.by = "Region_per_individual")
#selected_hvg <- SelectIntegrationFeatures(
#  seu_obj_list,
#  nfeatures = 3000,
#  fvf.nfeatures = 3000
#)
#combined <- FindVariableFeatures(object = combined, selection.method = "vst", nfeatures = 3000)
#VariableFeatures(combined) <- selected_hvg
#combined <- ScaleData(object = combined, verbose = T)
#combined <- RunPCA(object = combined, features = VariableFeatures(combined), verbose = F, npcs = 50)
#usefulPCs <- 1:20
#combined <- FindNeighbors(object = combined, dims = usefulPCs)
#combined <- FindClusters(object = combined, resolution = 1)
#combined <- RunUMAP(object = combined, dims = usefulPCs)
#combined$RNA_PCA_snn_res.1 <- combined$RNA_snn_res.1
#saveRDS(combined, file="Res_combined_adult_tissue_no_batch_effect_correction.rds")
#
## CSS integration
#combined <- readRDS("Res_combined_adult_tissue_no_batch_effect_correction.rds")
#sample_size <- sort(table(combined$Region_per_individual))
#css <- simspec::cluster_sim_spectrum(combined, label_tag = "Region_per_individual", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
#css <- RunUMAP(css, reduction = "css", dims = 1:ncol(css@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
#css <- FindNeighbors(object = css, reduction = "css", dims = 1:ncol(css@reductions$css@cell.embeddings), force.recalc = T) %>%
#  FindClusters(resolution = 1)
#css[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- css[[paste0("RNA_snn_res.", 0.2*5)]]
#saveRDS(css, file="Res_combined_adult_human_epi_with_CSS.rds")
#
#p1 <- SCpubr::do_DimPlot(css, reduction = "umap_css", group.by = "Unified_region_code", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10)
#p2 <- SCpubr::do_DimPlot(css, reduction = "umap_css", group.by = "Unified_coarse_cell_type", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10)
#p3 <- SCpubr::do_DimPlot(css, reduction = "umap_css", group.by = "RNA_CSS_snn_res.1", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10)
#p4 <- SCpubr::do_DimPlot(css, reduction = "umap_css", group.by = "paper", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10)
#png("Plot_UMAP_CSS_adult_intestine_tissue_epi.png", width=2000*2, height=2000*2)
#(p1+p2)/(p3+p4)
#dev.off()
#
## 
