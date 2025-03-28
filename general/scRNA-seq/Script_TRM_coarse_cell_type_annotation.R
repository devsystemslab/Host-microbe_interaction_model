setwd("/home/yuq22/group_share/Bacteria_TRM/TRM/TRM_from_coculture/cell_type_annotation/coarse_cell_type_annotation")
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(dplyr)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")

# load TRM data
trm <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/TRM/TRM_from_coculture/Dat_merged_filtered_TRM_from_coculture_no_epi.rds")
# load Epi from coculture data
epi <- readRDS("/home/yuq22/Bacteria_TRM/TRM/Dat_merged_filtered_ileum_epithelium_TRM_coculture_noImmune_noLowQuality_clustered.rds")

# coarse grain cell type annotation on TRM data
SCpubr::do_DimPlot(trm, group.by="RNA_snn_res.1", label=T)
# based on the similarity to gut cell atlas immune cell types, merged C4/5/8, C1/2/3/6/7, repectively
cell_type_anno <- list(
    "B" = 9,
    "ILC3-like CD4 T"=0,
    "SELL+ CD4 T"=10,
    "CD4 type 1 T"=c(1,2,3,6,7),
    "CD8 type 1 T"=c(4,5,8),
    "CD8 Tmem/NK"=11,
    "TRGV3/4/5 CD8 T"=12,
    "ILC3"=13,
    "G2M CD4"=14
)
trm$Cell_type <- NA
for(x in names(cell_type_anno)){
    trm$Cell_type[which(trm$RNA_snn_res.1 %in% cell_type_anno[[x]])] <- x
}
p1 <- SCpubr::do_DimPlot(trm, group.by="Cell_type", label=T)+NoLegend()
ggsave(p1, filename="Plot_UMAP_TRM_with_coarse_grain_cell_type_anno.png", width=6, height=6)
saveRDS(trm, file="Dat_TRM_with_coarse_grain_cell_type_anno.rds")

trm <- readRDS("Dat_TRM_with_coarse_grain_cell_type_anno.rds")
p1 <- SCpubr::do_DimPlot(trm, group.by="Cell_type", label=T)+NoLegend()
ggsave(p1, filename="Plot_UMAP_TRM_with_coarse_grain_cell_type_anno.png", width=6, height=6)
# load tissue immune cell data
tissue <- readRDS("/home/yuq22/Bacteria_TRM/intestine_multi_region_scRNA-seq_atlas/immune_cells/Res_gut_cell_atlas_adult_human_immune_no_batch_correction.rds")
# subset to T cells
tissue_t_cells <- subset(tissue, category=="T cells")
tissue_t_cells$Cell_type <- as.character(tissue_t_cells$Integrated_05)
size <- sort(table(tissue_t_cells$Cell_type))
saveRDS(tissue_t_cells, file="Dat_gut_cell_atlas_adult_human_immune_T_cells.rds")
trm_t_cells <- subset(trm, Cell_type!="B")
saveRDS(trm_t_cells, file="Dat_TRM_T_cells.rds")

tissue_t_cells <- readRDS("Dat_gut_cell_atlas_adult_human_immune_T_cells.rds")
trm_t_cells <- readRDS("Dat_TRM_T_cells.rds")
# identify cell type markers in tissue and in vitro data separately
de_res_tissue <- presto::wilcoxauc(tissue_t_cells, group_by="Cell_type")
de_res_tissue$pct_diff <- de_res_tissue$pct_in - de_res_tissue$pct_out  
sig_res_tissue <- de_res_tissue %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
sig_res_tissue$Cell_type <- paste("tissue", as.character(sig_res_tissue$group), sep=":")

de_res_trm <- presto::wilcoxauc(trm_t_cells, group_by="Cell_type")
de_res_trm$pct_diff <- de_res_trm$pct_in - de_res_trm$pct_out
sig_res_trm <- de_res_trm %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
sig_res_trm$Cell_type <- paste("inVitro", as.character(sig_res_trm$group), sep=":")

sig_res <- rbind(sig_res_tissue, sig_res_trm)
saveRDS(sig_res, file="Res_tissue_and_inVitro_T_cell_subtype_markers.rds")

# get the overlap of cell type markers between tissue and in vitro data
tissue_ct <- sort(unique(grep("tissue:", sig_res$Cell_type, value=T)))
vitro_ct <- sort(unique(grep("inVitro:", sig_res$Cell_type, value=T)))
n1 <- sapply(tissue_ct, function(ct1){
    sapply(vitro_ct, function(ct2){
        tissue_markers <- sig_res$feature[which(sig_res$Cell_type==ct1)]
        vitro_markers <- sig_res$feature[which(sig_res$Cell_type==ct2)]
        length(intersect(tissue_markers, vitro_markers))
    })
})
# test the overlap significance
union_markers <- sort(unique(sig_res$feature))
pval <- sapply(tissue_ct, function(ct1){
    sapply(vitro_ct, function(ct2){
        tissue_ct1_markers <- sig_res$feature[which(sig_res$Cell_type==ct1)]
        vitro_ct2_markers <- sig_res$feature[which(sig_res$Cell_type==ct2)]
        overlap_markers <- intersect(tissue_ct1_markers, vitro_ct2_markers)
        ct1_non_ct2_markers <- setdiff(tissue_ct1_markers, vitro_ct2_markers)
        ct2_non_ct1_markers <- setdiff(vitro_ct2_markers, tissue_ct1_markers)
        non_ct1_ct2_markers <- setdiff(union_markers, c(tissue_ct1_markers, vitro_ct2_markers))
        a <- length(overlap_markers)
        b <- length(ct1_non_ct2_markers)
        c <- length(ct2_non_ct1_markers)
        d <- length(non_ct1_ct2_markers)
        fisher.test(matrix(c(a, b, c, d), nrow=2))$p.value
    })
})
saveRDS(pval, file="Res_cell_type_marker_overlap_significance_pval.rds")
pval <- readRDS("Res_cell_type_marker_overlap_significance_pval.rds")

res <- -log10(pval)-(-log10(0.05))
rownames(res) <- sub("inVitro:", "", rownames(res))
colnames(res) <- sub("tissue:", "", colnames(res))
pdf("Plot_heatmap_T_cell_type_marker_overlap_between_chip_and_tissue.pdf", height=10, width=10)
gplots::heatmap.2(res, scale="row", cellnote=round(res,1), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10), cexRow=1, cexCol=1, xlab="Tissue", ylab="inVitro")
dev.off()

# for each invitro cell type, get the top 3 tissue cell type with the strongest overlap, unless the difference between the top 1 and top 2/3 is larger than 10
# do not include in vitro G2M CD4 for this analysis
vitro_ct <- sub("inVitro:", "", setdiff(all_vitro_ct, "inVitro:G2M CD4"))
matched_tissue_cts <- lapply(vitro_ct, function(ct){
    vec <- res[ct, ]
    top3 <- order(vec, decreasing=T)[1:3]
    if(vec[top3[1]]-vec[top3[2]]>10){
        idx <- top3[1]
    }else if(vec[top3[1]]-vec[top3[3]]>10){
        idx <- top3[1:2]
    }else{
        idx <- top3
    }
    names(vec)[idx]

})
names(matched_tissue_cts) <- vitro_ct
# CD4 cells won't be matched to any CD8 cells, so take top 1,2,4
matched_tissue_cts[["CD4 type 1 T"]] <- c(matched_tissue_cts[["CD4 type 1 T"]][1:2],"Treg")


# calculate the similarity of T cell subtypes between tissue and in vitro
tissue_expr <- getAveExpr(seu.obj=tissue_t_cells, feature.to.calc="Cell_type", colname.prefix="tissue:")
vitro_expr <- getAveExpr(seu.obj=trm_t_cells, feature.to.calc="Cell_type", colname.prefix="inVitro:")
saveRDS(tissue_expr, file="Res_tissue_T_cell_subtype_expression.rds")
saveRDS(vitro_expr, file="Res_inVitro_T_cell_subtype_expression.rds")

# get markers for each cell type
sig_res <- readRDS("Res_tissue_and_inVitro_T_cell_subtype_markers.rds")
vitro_marker_genes <- sort(unique(sig_res$feature[grep("inVitro", sig_res$Cell_type)]))
tissue_marker_genes <- sort(unique(sig_res$feature[grep("tissue", sig_res$Cell_type)]))
both_marker_genes <- intersect(vitro_marker_genes, tissue_marker_genes)

scc_mat <- cor(vitro_expr[both_marker_genes, ], tissue_expr[both_marker_genes, ], method="spearman")
res <- scc_mat
rownames(res) <- sub("inVitro:", "", rownames(res))
colnames(res) <- sub("tissue:", "", colnames(res))
pdf("Plot_heatmap_T_cell_type_SCC_between_chip_and_tissue.pdf", height=10, width=10)
gplots::heatmap.2(res, scale="row", cellnote=round(res,2), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10), cexRow=1, cexCol=1, xlab="Tissue", ylab="inVitro")
dev.off()

# for each invitro cell type, get the top 3 most similar tissue cell type, unless the difference between the top 1 and top 2/3 is larger than 0.1
# do not include in vitro G2M CD4 for this analysis
vitro_ct <- sub("inVitro:", "", setdiff(all_vitro_ct, "inVitro:G2M CD4"))
scc_matched_tissue_cts <- lapply(vitro_ct, function(ct){
    vec <- res[ct, ]
    top3 <- order(vec, decreasing=T)[1:3]
    if(vec[top3[1]]-vec[top3[2]]>0.1){
        idx <- top3[1]
    }else if(vec[top3[1]]-vec[top3[3]]>0.1){
        idx <- top3[1:2]
    }else{
        idx <- top3
    }
    names(vec)[idx]

})
names(scc_matched_tissue_cts) <- vitro_ct
# CD4 cells won't be matched to any CD8 cells, so take top 1,2,4
scc_matched_tissue_cts[["CD4 type 1 T"]] <- c(scc_matched_tissue_cts[["CD4 type 1 T"]][1:2],"Activated CD4 T")

# get the overlap of overlap-based and similarity-based matched tissue cell types
overlap_matched_tissue_cts <- lapply(vitro_ct, function(ct){
    intersect(matched_tissue_cts[[ct]], scc_matched_tissue_cts[[ct]])
})
names(overlap_matched_tissue_cts) <- vitro_ct
# inVitro CD8 Tmem/NK also show high similarity to tissue NK cells (SCC=0.61, rank 4; rank 1 SCC=0.66)
# so manually match inVitro CD8 Tmem/NK to tissue NK
overlap_matched_tissue_cts[["CD8 Tmem/NK"]] <- "NK cell"

# update cell type annotation
ct1 <- c("CD4 type 1 T", "CD8 Tmem/NK", "CD8 type 1 T", "ILC3-like CD4 T", "SELL+ CD4 T", "TRGV3/4/5 CD8 T")
ct2 <- c("Activated CD4 T", "NK cell", "Activated CD8 T", "MAIT-like cell", "Treg", "TRGV3/4/5 gdT")
trm$Updated_cell_type <- trm$Cell_type
for(x in ct1){
    trm$Updated_cell_type[which(trm$Cell_type==x)] <- ct2[which(ct1==x)]
}
df <- unique(trm@meta.data[,c("Cell_type", "Updated_cell_type")])
trm$Cell_type <- NULL
p1 <- SCpubr::do_DimPlot(trm, group.by="Updated_cell_type", label=T)+NoLegend()
ggsave(p1, filename="Plot_UMAP_TRM_with_updated_coarse_grain_cell_type_anno.png", width=6, height=6)
saveRDS(trm, file="Dat_TRM_with_updated_cell_type_anno.rds")

# update the heatmap of cell type marker overlap
res <- -log10(pval)-(-log10(0.05))
rownames(res) <- sub("inVitro:", "", rownames(res))
colnames(res) <- sub("tissue:", "", colnames(res))
rownames(res) <- df$Updated_cell_type[match(rownames(res), df$Cell_type)]
pdf("Plot_heatmap_T_cell_type_marker_overlap_between_chip_and_tissue_updated_ct.pdf", height=10, width=10)
gplots::heatmap.2(res, scale="row", cellnote=round(res,1), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10), cexRow=1, cexCol=1, xlab="Tissue", ylab="inVitro")
dev.off()

res <- scc_mat
rownames(res) <- sub("inVitro:", "", rownames(res))
colnames(res) <- sub("tissue:", "", colnames(res))
rownames(res) <- df$Updated_cell_type[match(rownames(res), df$Cell_type)]
pdf("Plot_heatmap_T_cell_type_SCC_between_chip_and_tissue_updated.pdf", height=10, width=10)
gplots::heatmap.2(res, scale="row", cellnote=round(res,2), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10), cexRow=1, cexCol=1, xlab="Tissue", ylab="inVitro")
dev.off()

# for each invitro cell type, get the top 3 tissue cell type with the strongest overlap, unless the difference between the top 1 and top 2/3 is larger than 10
# do not include in vitro G2M CD4 for this analysis
vitro_ct <- sub("inVitro:", "", setdiff(all_vitro_ct, "inVitro:G2M CD4"))
matched_tissue_cts <- lapply(vitro_ct, function(ct){
    vec <- res[ct, ]
    top3 <- order(vec, decreasing=T)[1:3]
    if(vec[top3[1]]-vec[top3[2]]>10){
        idx <- top3[1]
    }else if(vec[top3[1]]-vec[top3[3]]>10){
        idx <- top3[1:2]
    }else{
        idx <- top3
    }
    names(vec)[idx]

})
names(matched_tissue_cts) <- vitro_ct
# CD4 cells won't be matched to any CD8 cells, so take top 1,2,4
matched_tissue_cts[["CD4 type 1 T"]] <- c(matched_tissue_cts[["CD4 type 1 T"]][1:2],"Treg")

# update the cell type annotation in the sig_res object
sig_res_2 <- sig_res
for(ct in df$Cell_type){
    sig_res_2$Cell_type[which(sig_res$Cell_type==paste0("inVitro:", ct))] <- paste("inVitro", df$Updated_cell_type[which(df$Cell_type==ct)], sep=":")
}
saveRDS(sig_res_2, file="Res_tissue_and_inVitro_T_cell_subtype_markers_ct_updated.rds")
sig_res <- sig_res_2
top_res <- sig_res %>% group_by(Cell_type) %>% top_n(20, logFC)
all_tissue_ct <- sort(unique(grep("tissue:", sig_res$Cell_type, value=T)))
all_vitro_ct <- sort(unique(grep("inVitro:", sig_res$Cell_type, value=T)))


tissue_ct <- "tissue:Activated CD4 T"
vitro_ct <- "inVitro:Activated CD4 T"
tissue_markers <- top_res$feature[which(top_res$Cell_type%in%tissue_ct)]
vitro_markers <- top_res$feature[which(top_res$Cell_type==vitro_ct)]
genes <- intersect(tissue_markers, vitro_markers)

tissue <- readRDS("/home/yuq22/Bacteria_TRM/intestine_multi_region_scRNA-seq_atlas/immune_cells/Res_gut_cell_atlas_adult_human_immune_no_batch_correction.rds")
# subset to non G2M B cells
cells <- colnames(tissue)[which(tissue$Cell_type%in%c("Memory B", "Naive B"))]
tissue_b_cells <- subset(tissue, cells=cells)
saveRDS(tissue_b_cells, file="Dat_gut_cell_atlas_adult_human_immune_B_cells.rds")
# identify cell type markers in tissue and in vitro data separately
de_res_tissue <- presto::wilcoxauc(tissue_b_cells, group_by="Cell_type")
de_res_tissue$pct_diff <- de_res_tissue$pct_in - de_res_tissue$pct_out  
sig_res_tissue <- de_res_tissue %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
sig_res_tissue$Cell_type <- paste("tissue", as.character(sig_res_tissue$group), sep=":")
top_res <- sig_res_tissue %>% group_by(Cell_type) %>% top_n(10, logFC)
genes <- top_res$feature

# update cell type annotation in invivo data from "TRGV3/4/5 gdT" to "gdT"
trm$Updated_cell_type[which(trm$Updated_cell_type=="TRGV3/4/5 gdT")] <- "gdT"
trm$Updated_cell_type[which(trm$Updated_cell_type=="B")] <- "Memory B"
trm$Updated_cell_type[which(trm$Updated_cell_type=="MAIT-like cell")] <- "MAIT-like CD4"
saveRDS(trm, file="Dat_TRM_with_updated_cell_type_anno.rds")

# update the heatmap of cell type marker overlap
df$Updated_cell_type[which(df$Updated_cell_type=="TRGV3/4/5 gdT")] <- "gdT"
res <- -log10(pval)-(-log10(0.05))
rownames(res) <- sub("inVitro:", "", rownames(res))
colnames(res) <- sub("tissue:", "", colnames(res))
rownames(res) <- df$Updated_cell_type[match(rownames(res), df$Cell_type)]
pdf("Plot_heatmap_T_cell_type_marker_overlap_between_chip_and_tissue_updated_ct.pdf", height=10, width=10)
gplots::heatmap.2(res, scale="row", cellnote=round(res,1), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10), cexRow=1, cexCol=1, xlab="Tissue", ylab="inVitro")
dev.off()

res <- scc_mat
rownames(res) <- sub("inVitro:", "", rownames(res))
colnames(res) <- sub("tissue:", "", colnames(res))
rownames(res) <- df$Updated_cell_type[match(rownames(res), df$Cell_type)]
pdf("Plot_heatmap_T_cell_type_SCC_between_chip_and_tissue_updated.pdf", height=10, width=10)
gplots::heatmap.2(res, scale="row", cellnote=round(res,2), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10), cexRow=1, cexCol=1, xlab="Tissue", ylab="inVitro")
dev.off()

tissue_subset <- subset(tissue, cells=colnames(tissue)[which(tissue$category%in%c("B cells", "T cells"))])
marker_list <- list(
    "G2M CD4"=c("CD4", "MKI67"),
    "Memory B"=c("MZB1", "CD79B", "LINC01781", "CD27", "TNFRSF13B"),
    "ILC3"=c("KIT", "ZFP36L1"),
    "NK"=c("FCER1G", "GZMK",   "CTSW", "CCL4",   "TYROBP", "NKG7"),
    "MAIT"=c("SLC4A10", "KLRB1"),
    "Treg"=c("LTB", "ARID5B", "CCR7"),
    "gdT"=c("CD160", "GNLY", "KLRC2", "KLRD1", "TRG-AS1"),
    "Activated CD8 T"=c("XCL1", "CD8A", "CD8B", "IFNG"),
    "Activated CD4 T"=c("RGCC", "GPR183")
)

p1 <- SCpubr::do_DotPlot(trm, features=marker_list, group.by="Updated_cell_type")
p2 <- SCpubr::do_DotPlot(tissue_subset, features=marker_list, group.by="Cell_type")
p <- p1/p2
ggsave(p, filename="Plot_DotPlot_TRM_vs_tissue_immune_cell_type_markers.pdf", width=10, height=15)

# identify cell type markers
de_res <- presto::wilcoxauc(trm, group_by="Updated_cell_type")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
top_res <- sig_res %>% group_by(group) %>% top_n(20, wt=logFC)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_duo_TRM_immune_cell_type_markers.rds")

# get cluster average expression
cluster_expr <- getAveExpr(seu.obj=trm, feature.to.calc="Updated_cell_type", colname.prefix=NULL)
saveRDS(cluster_expr, file="Dat_duo_TRM_immune_cell_type_expr.rds")

genes <- unique(top_res$feature)
expr_mat <- cluster_expr[genes,]
res <- prepareTreeAndHeatmapInput(expr.mat = expr_mat, hc.method = "ward.D2", norm.method = "quantile")
saveRDS(res, file="Res_duo_TRM_immune_cell_type_marker_expr_heatmap_input.rds")
pdf("Plot_heatmap_duo_TRM_immune_cell_type_marker_expr.pdf", height=10)
gplots::heatmap.2(res$heatmap_input, col = beach.col.heatmap, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 1, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

pdf("Plot_hc_duo_TRM_immune_cell_type_marker_expr.pdf")
plot(res$hc_col, main="Hierarchical clustering of adult tissue subcluster marker expression", hang=-1, cex=0.5)
dev.off()

# add cell type index to the trm data
idx <- setNames(paste(seq(ncol(res$heatmap_input)), colnames(res$heatmap_input), sep=":"), colnames(res$heatmap_input))
trm$Updated_cell_type_index <- idx[trm$Updated_cell_type]
trm$Cell_type <- trm$Updated_cell_type
saveRDS(trm, file="Dat_TRM_with_updated_cell_type_anno.rds")
saveRDS(trm, file="/home/yuq22/group_share/Bacteria_TRM/used_object/Dat_TRM_with_updated_cell_type_anno.rds")
p1 <- SCpubr::do_DimPlot(trm, group.by="Updated_cell_type_index")+NoLegend()
ggsave(p1, filename="Plot_UMAP_TRM_with_updated_coarse_grain_cell_type_anno.png", width=6, height=6)


# subset tissue immune data to those matched with in invitro data
trm$Source <- "inVitro"
selected_cell_type <- c("gdT", "ILC3", "Memory B", "Activated CD4 T", "Activated CD8 T", "Treg", "NK cell", "MAIT cell")
tissue_subset <- subset(tissue, cells=colnames(tissue)[which(tissue$Cell_type%in%selected_cell_type)])
tissue_subset$Source <- "tissue"

ct1 <- c("Activated CD4 T", "Activated CD8 T", "gdT", "ILC3", "MAIT cell", "Memory B", "NK cell", "Treg")  
ct2 <- setNames(c("8:Activated CD4 T", "3:Activated CD8 T", "4:gdT", "5:ILC3", "9:MAIT cell", "1:Memory B", "2:NK cell", "7:Treg"), ct1)
tissue_subset$Updated_cell_type_index <- ct2[as.character(tissue_subset$Cell_type)]
combined <- merge(x=trm, y=tissue_subset)
combined$Cell_type_per_source <- paste(combined$Updated_cell_type_index, combined$Source, sep=":")
combined$Cell_type_per_source <- factor(combined$Cell_type_per_source, levels=sort(unique(combined$Cell_type_per_source)))
Idents(combined) <- combined$Cell_type_per_source
saveRDS(combined, file="Dat_combined_immune_inVitro_and_tissue_data.rds")
saveRDS(combined, file="/home/yuq22/group_share/Bacteria_TRM/used_object/Dat_combined_immune_inVitro_and_tissue_data.rds")
marker_list <- list(
    "Memory B"=c("CD79B", "LINC01781", "TNFRSF13B"),
    "NK"=c("GZMB", "GZMK", "NKG7"),
    "Act. CD8"=c("XCL1", "CD8A", "CD8B"),
    "gdT"=c("CD160", "KLRC2", "KLRD1"),
    "ILC3"=c("KIT", "ZFP36L1"),
    "G2M"=c("CD4", "MKI67"),
    "Treg"=c("LTB", "ARID5B", "CCR7"),
    "Act. CD4"=c("RGCC", "GPR183"),
    "MAIT"=c("SLC4A10", "KLRB1")
)


p2 <- SCpubr::do_DotPlot(sample = combined, 
             features = marker_list,
             scale=T)
p2

ggsave(p2, filename="Plot_dotplot_adult_human_tissue_and_inVitro_immune_cell_type_markers.pdf", width=8, height=6)



# identify cell type markers
trm <- readRDS("/home/yuq22/group_share/Bacteria_TRM/used_object/Dat_TRM_with_updated_cell_type_anno.rds")
# focus on expressed genes not starting with ENSG and not from confounding genes
gene <- rownames(trm)[!grepl("^ENSG", rownames(trm))]
#remove ribosomal genes, mitochondrial genes, and genes located on sex chromosomes before normalization
confound_genes <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/confound_genes/Dat_confound_genes.rds")
genes_to_exclude <- unique(unlist(confound_genes[c(4,5,6)]))
genes <- setdiff(gene, genes_to_exclude)
X <- trm@assays$RNA@data[genes, ]
y <- trm$Updated_cell_type_index
de_res <- presto::wilcoxauc(X=X, y=y)
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
top_res <- sig_res %>% group_by(group) %>% top_n(20, wt=logFC)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_duo_TRM_immune_cell_type_markers.rds")
write.table(top_res, file="Res_duo_TRM_immune_cell_type_markers.txt", sep="\t", quote=F, row.names=F)



