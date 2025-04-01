setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data")
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
# load data
combined <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/Res_chip_infection_multiome.rds")
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
combined <- FindClusters(object = combined, algorithm = 3, resolution=0.6)
saveRDS(combined, file="Dat_filtered_colon_epithelium_infection_ATAC_preprocessed.rds")

combined <- readRDS("Dat_filtered_colon_epithelium_infection_ATAC_preprocessed.rds")
p1 <- SCpubr::do_DimPlot(combined, reduction="umap.atac", group.by="RNA_snn_res.1", label=TRUE)+NoLegend()
p2 <- SCpubr::do_DimPlot(combined, reduction="umap.atac", group.by="condition")
p3 <- SCpubr::do_DimPlot(combined, reduction="umap.atac", group.by="peaks_merged_snn_res.0.6", label=TRUE)
p4 <- SCpubr::do_DimPlot(combined, reduction="umap.atac", group.by="Subtype", label=TRUE)
p <- cowplot::plot_grid(plotlist=list(p1, p2, p3, p4), nrow=2)
ggsave(p, filename="Plot_UMAP_ATAC.png", width=12, height=12)

# check the cell type distribution on the ATAC cluster level
n1 <- sapply(sort(unique(combined$peaks_merged_snn_res.0.6)), function(cl){
  sapply(sort(unique(combined$Subtype)), function(ct){
    sum(combined$Subtype==ct & combined$peaks_merged_snn_res.0.6==cl)
  })
})
colnames(n1) <- paste0("C", sort(unique(combined$peaks_merged_snn_res.0.6)))
p1 <- t(t(n1)/colSums(n1))
p_diff <- apply(p1, 2, function(vec){
    sorted_vec <- sort(vec, decreasing=TRUE)
    sorted_vec[1] - sorted_vec[2]
})
barplot(p_diff)
# only EEC, Goblet, BEST4+ and the infection-specific colonocyte populations show clear enrichment
# so perform DAR analysis a broader level based on caterogization on RNA
# identify marker regions
de_res <- presto::wilcoxauc(combined, group_by="peaks_merged_snn_res.0.6", seurat_assay="peaks_merged")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
top_res <- sig_res %>% group_by(group) %>% top_n(200, wt=logFC)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_colon_epithelium_Salmonella_infection_dataset_cluster_marker_regions.rds")

# check the correspondence between RNA clusters and ATAC clusters
n1 <- sapply(sort(unique(combined$peaks_merged_snn_res.0.6)), function(cl_atac){
  sapply(sort(unique(combined$RNA_snn_res.1)), function(cl_rna){
    sum(combined$RNA_snn_res.1==cl_rna & combined$peaks_merged_snn_res.0.6==cl_atac)
  })
})
colnames(n1) <- paste0("ATAC_", sort(unique(combined$peaks_merged_snn_res.0.6)))
rownames(n1) <- paste0("RNA_", sort(unique(combined$RNA_snn_res.1)))


# get the condition composition of each cluster
n1 <- sapply(sort(unique(combined$RNA_snn_res.1)), function(cl){
  sapply(sort(unique(combined$condition)), function(condition){
    sum(combined$condition==condition & combined$RNA_snn_res.1==cl)
  })
})
colnames(n1) <- paste0("C", sort(unique(combined$RNA_snn_res.1)))
infection_specific_clusters <- names(which(n1["Control",]<20))

# check condition distribution of ATAC clusters
n1 <- sapply(sort(unique(combined$peaks_merged_snn_res.0.6)), function(cl){
  sapply(sort(unique(combined$condition)), function(condition){
    sum(combined$condition==condition & combined$peaks_merged_snn_res.0.6==cl)
  })
})
colnames(n1) <- paste0("C", sort(unique(combined$peaks_merged_snn_res.0.6)))
infection_specific_clusters <- names(which(n1["Control",]<20))
selected_cl <- setdiff(combined$peaks_merged_snn_res.0.6, infection_specific_clusters)
df_out <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/Res_Salmonella_and_Control_nearest_neighbors_PCA_based.rds")
# for infection specific clusters, identify DAR between infection cluster and matched cells
# for non-infection-specific clusters, identify DAR between co-cluster infected and control cells
cluster_list <- list(
    "infection_specific_clusters"=sub("C", "", infection_specific_clusters),
    "non_infection_specific_cl"=selected_cl
)
# identify DARs of goblet, BEST4+ cells, EECs, Colon_SC, 
# infection-speicfic colonocyte subtype 1 (RNA-based cluster 14 and 17, which mainly correspond to ATAC-based cluster 9),
# infection-specific colonocyte subtype 2 (RNA-based cluster 16, which mainly correspond to ATAC-based cluster 10),
# and non-infection-specific coloncytes (other RNA-based colonocyte clusters)
group_vec <- combined$Subtype
group_vec[which(combined$RNA_snn_res.1%in%c(14, 17))] <- "Infection_specific_colonocyte_subtype_1"
group_vec[which(combined$RNA_snn_res.1==16)] <- "Infection_specific_colonocyte_subtype_2"
group_vec[which(!combined$RNA_snn_res.1%in%c(14,16,17) & grepl("Colonocyte_", combined$Subtype))] <- "Non_infection_specific_colonocyte"
combined$RNA_coarse_cell_type <- group_vec
df <- unique(combined@meta.data[,c("Subtype", "RNA_coarse_cell_type", "RNA_snn_res.1")])
rownames(df) <- NULL
df <- df[order(df$RNA_coarse_cell_type),]
saveRDS(combined, file="Dat_filtered_colon_epithelium_infection_ATAC_preprocessed.rds")

combined <- readRDS("Dat_filtered_colon_epithelium_infection_ATAC_preprocessed.rds")
# DAR on coarser cell type level
de_res_list <- lapply(sort(unique(combined$RNA_coarse_cell_type)), function(x){
  salmonella_cells <- colnames(combined)[which(combined$RNA_coarse_cell_type==x & combined$condition=="Salmonella")]
  if(grepl("Infection_specific", x)){
      matched_control_cells <- unique(df_out$Control_Cell[which(df_out$Salmonella_Cell%in%salmonella_cells)])
  }else{
      rna_clusters <- unique(combined$RNA_snn_res.1[which(combined$RNA_coarse_cell_type==x)])
      matched_control_cells <- unique(colnames(combined)[which(combined$RNA_snn_res.1%in%rna_clusters & combined$condition=="Control")])
  }

  cells <- c(salmonella_cells, matched_control_cells)
  seu_obj <- subset(combined, cells=cells)
  de_res <- presto::wilcoxauc(seu_obj, group_by="condition", seurat_assay="peaks_merged")
  de_res$group <- paste(x, de_res$group, sep="_")
  de_res$pct_diff <- de_res$pct_in - de_res$pct_out
  return(de_res)
})
names(de_res_list) <- sort(unique(combined$RNA_coarse_cell_type))
combined_de_res <- do.call('rbind', de_res_list)
combined_sig_res <- combined_de_res %>% filter(padj<0.05 & pct_diff>5 & logFC>0.01)
saveRDS(combined_de_res, file="Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_combined_control_cells_DAR_all_res_on_RNA_cell_type_level.rds")
saveRDS(combined_sig_res, file="Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_combined_control_cells_DAR_sig_res_on_RNA_cell_type_level.rds")
pdf("Plot_DAR_per_cell_type.pdf", width=10)
par(mar=c(4, 22, 4, 2))
barplot(table(combined_sig_res$group), las=2, xlab="Number of DARs", horiz=TRUE)
dev.off()
union_dar <- sort(unique(combined_sig_res$feature))
# method to be improved

## DAR on ATAC cluster level (not recommended as the number of DARs per condition is still not balanced and it is harder to explain)
#de_res_list_2 <- lapply(sort(unique(combined$peaks_merged_snn_res.0.6)), function(x){
#  salmonella_cells <- colnames(combined)[which(combined$peaks_merged_snn_res.0.6==x & combined$condition=="Salmonella")]
#  if(x %in% cluster_list[["infection_specific_clusters"]]){
#      matched_control_cells <- unique(df_out$Control_Cell[which(df_out$Salmonella_Cell%in%salmonella_cells)])
#  }else{
#      matched_control_cells <- unique(colnames(combined)[which(combined$peaks_merged_snn_res.0.6==x & combined$condition=="Control")])
#  }
#
#  cells <- c(salmonella_cells, matched_control_cells)
#  seu_obj <- subset(combined, cells=cells)
#  de_res <- presto::wilcoxauc(seu_obj, group_by="condition", seurat_assay="peaks_merged")
#  de_res$group <- paste(x, de_res$group, sep="_")
#  de_res$pct_diff <- de_res$pct_in - de_res$pct_out
#  return(de_res)
#})
#names(de_res_list_2) <- paste0("C", sort(unique(combined$peaks_merged_snn_res.0.6)))
#combined_de_res_2 <- do.call('rbind', de_res_list_2)
#combined_sig_res_2 <- combined_de_res_2 %>% filter(padj<0.05 & pct_diff>5 & logFC>0.01)
#saveRDS(combined_de_res_2, file="Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_combined_control_cells_DAR_all_res_on_atac_cluster_level.rds")
#saveRDS(combined_sig_res_2, file="Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_combined_control_cells_DAR_sig_res_on_atac_cluster_level.rds")
#pdf("Plot_DAR_per_ATAC_cluster.pdf", width=10)
#par(mar=c(4, 22, 4, 2))
#barplot(table(combined_sig_res_2$group), las=2, xlab="Number of DARs", horiz=TRUE)
#dev.off()
#union_dar_2 <- sort(unique(combined_sig_res_2$feature))

library(Pando)
seurat_grn <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/Res_chip_infection_with_GRN.rds")
modules <- NetworkModules(seurat_grn) 
saveRDS(modules, file="Res_chip_infection_GRN_modules.rds")

module_list <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/Analysis/Res_salmonella_infection_GRN_module_list.rds")
selected_modules <- unique(unlist(module_list))
peak_list_pos <- modules@features$peaks_pos
peak_list_neg <- modules@features$peaks_neg
names(peak_list_pos) <- paste0("pos_", names(peak_list_pos))
names(peak_list_neg) <- paste0("neg_", names(peak_list_neg))
peak_list_grn <- c(peak_list_pos, peak_list_neg)
grn_regions <- unique(unlist(peak_list_grn))
selected_grn_regions <- unique(unlist(peak_list_grn[selected_modules]))
grn_dar <- intersect(selected_grn_regions, union_dar)

cluster_marker_res <- deg_res$sig_res
infection_specific_cluster_marker_regions <- unique(cluster_marker_res$feature[which(cluster_marker_res$group%in%c(9,10))])

grn_marker_dar <- intersect(grn_dar, infection_specific_cluster_marker_regions)
length(grn_marker_dar)
saveRDS(grn_marker_dar, file="Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_combined_control_cells_DAR_GRN_marker_regions.rds")
grn_marker_dar <- readRDS("Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_combined_control_cells_DAR_GRN_marker_regions.rds")

# load scRNA-seq data
# create bi-modal object
combined_rna <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/Dat_salmonella_colon_epithelium_scRNA-seq_data_with_subtypes.rds")
combined_rna$Coarse_cell_type_per_condition <- paste(combined_rna$Subtype, combined_rna$condition, sep=":")
combined <- readRDS("Dat_filtered_colon_epithelium_infection_ATAC_preprocessed.rds")
combined_rna$RNA_subtype_group <- combined$RNA_coarse_cell_type
combined_rna$RNA_subtype_group_per_condition <- paste(combined_rna$RNA_subtype_group, combined_rna$condition, sep=":")
combined_rna[["peaks_merged"]] <- combined[["peaks_merged"]]
combined_rna$peaks_merged_snn_res.0.6 <- combined$peaks_merged_snn_res.0.6
atac_coor <- Embeddings(combined, reduction="umap.atac")
combined_rna[["umap.atac"]] <- CreateDimReducObject(embeddings = atac_coor, key="ATACUMAP_", assay="peaks_merged")
SCpubr::do_DimPlot(combined_rna, reduction="umap.atac", group.by="peaks_merged_snn_res.0.6", pt.size=6, font.size=100, border.size=1.5, legend.icon.size=40, label=T, label.size=20)+NoLegend()
SCpubr::do_DimPlot(combined_rna, reduction="umap", group.by="peaks_merged_snn_res.0.6", pt.size=6, font.size=100, border.size=1.5, legend.icon.size=40, label=T, label.size=20)+NoLegend()
saveRDS(combined_rna, file="Dat_filtered_colon_epithelium_infection_RNA_ATAC_preprocessed_ctByCondition.rds")

combined_rna <- readRDS("Dat_filtered_colon_epithelium_infection_RNA_ATAC_preprocessed_ctByCondition.rds")
# generate heatmap of regions with top20 logFC between conditions per cell type
# exclude TA cels
combined_de_res <- readRDS("Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_combined_control_cells_DAR_all_res_on_RNA_cell_type_level.rds")
# rank by logFC
combined_de_res[grepl("Salmonella", combined_de_res$group) & !grepl("TA", combined_de_res$group),] %>% group_by(group) %>% top_n(20, wt=logFC) %>% select(feature) %>% pull() %>% unique()-> sal_up_regions
combined_de_res[grepl("Control", combined_de_res$group) & !grepl("TA", combined_de_res$group),] %>% group_by(group) %>% top_n(20, wt=logFC) %>% select(feature) %>% pull() %>% unique()-> con_up_regions
# rank by pct_diff # the heatmap is more messy than the logFC-based one, not recommended
#combined_de_res[grepl("Salmonella", combined_de_res$group) & !grepl("TA", combined_de_res$group),] %>% group_by(group) %>% top_n(20, wt=pct_diff) %>% select(feature) %>% pull() %>% unique()-> sal_up_regions
#combined_de_res[grepl("Control", combined_de_res$group) & !grepl("TA", combined_de_res$group),] %>% group_by(group) %>% top_n(20, wt=pct_diff) %>% select(feature) %>% pull() %>% unique()-> con_up_regions

# get average accessibility per cell type per condition
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")
ave_acc <- getAveExpr(seu.obj=combined_rna, assay.type="peaks_merged", feature.to.calc="RNA_subtype_group_per_condition", colname.prefix=NULL, size.cutoff=20)
saveRDS(ave_acc, file="Dat_average_accessibility_per_RNA_subtype_group_per_condition.rds")
# exclude TA cells
ave_acc <- ave_acc[,!grepl("TA", colnames(ave_acc))]
saveRDS(ave_acc, file="Dat_average_accessibility_per_RNA_subtype_group_per_condition_noTA.rds")

top_dar <- unique(c(sal_up_regions, con_up_regions))
input <- ave_acc[top_dar,]
res <- prepareTreeAndHeatmapInput(expr.mat = input, hc.method = "ward.D2", norm.method = "quantile")
saveRDS(res, file="Res_salmonella_up_and_down_region_heatmap_input.rds")
pdf("Plot_heatmap_salmonella_up_and_down_region_accessibility.pdf", height=10)
gplots::heatmap.2(res$heatmap_input, col = beach.col.heatmap, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 1, cexCol=0.5, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

pdf("Plot_hc_salmonella_up_and_down_regio_accessibility.pdf")
plot(res$hc_col, main="Hierarchical clustering of adult tissue subcluster marker expression", hang=-1, cex=0.5)
dev.off()

# highlight the infection upregulated regions
setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/")
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")
combined_de_res <- readRDS("Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_combined_control_cells_DAR_all_res_on_RNA_cell_type_level.rds")
# rank by logFC
# exclude BEST4+ cell, EEC and TA
combined_de_res[grepl("Salmonella", combined_de_res$group) & !grepl("TA|BEST4|EEC", combined_de_res$group),] %>% group_by(group) %>% top_n(20, wt=logFC) %>% select(feature) %>% pull() %>% unique()-> sal_up_regions
combined_de_res[grepl("Control", combined_de_res$group) & !grepl("TA|BEST4|EEC", combined_de_res$group),] %>% group_by(group) %>% top_n(20, wt=logFC) %>% select(feature) %>% pull() %>% unique()-> con_up_regions
ave_acc <- readRDS("Dat_average_accessibility_per_RNA_subtype_group_per_condition_noTA.rds")
top_dar <- union(sal_up_regions, con_up_regions)
col_orders <- c(
"Colon_SC:Salmonella", 
"Colon_SC:Control",                                  
"Infection_specific_colonocyte_subtype_1:Salmonella",
"Infection_specific_colonocyte_subtype_2:Salmonella", 
"Non_infection_specific_colonocyte:Salmonella", 
"Non_infection_specific_colonocyte:Control",                        
"Goblet:Salmonella", 
"Goblet:Control"                                    
)                                
input <- ave_acc[top_dar,col_orders]
res <- prepareTreeAndHeatmapInput(expr.mat = input, hc.method = "ward.D2", norm.method = "quantile", column.reorder = FALSE)
saveRDS(res, file="Res_GC_s2E_salmonella_up_and_down_region_heatmap_input.rds")
pdf("Plot_heatmap_GC_s2E_salmonella_up_and_down_region_accessibility.pdf", height=10)
gplots::heatmap.2(res$heatmap_input, col = pink.col.heatmap, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 1, cexCol=0.5, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()


# perform GREAT analysis on salmonella infection up-regulated regions in each of non-TA cell types
sal_up_200 <- combined_de_res[grepl("Salmonella", combined_de_res$group) & !grepl("TA", combined_de_res$group),] %>% filter(logFC>0) %>% group_by(group) %>% top_n(200, wt=logFC)
# perform GREAT enrichment analysis on the selected peaks lists
library(rGREAT)
library(Signac)

great_res_list <- list()
for(ct in unique(sal_up_200$group)){
  peaks <- sal_up_200$feature[which(sal_up_200$group==ct)]
  great_input <- StringToGRanges(peaks)
  job = submitGreatJob(great_input, species = "hg38")
  tbl = getEnrichmentTables(job)
  res <- tbl[["GO Biological Process"]]
  res$sig_idx <- res$Binom_Adjp_BH<0.05
  saveRDS(res, file=paste0("Res_salmonella_infection_on_", ct, "_enriched_peaks_functional_BP.rds"))
  great_res_list[[ct]] <- res
}
saveRDS(great_res_list, file="Res_salmonella_infection_on_non_TA_cell_types_enriched_peaks_functional_BP.rds")

great_res_list <- readRDS("Res_salmonella_infection_on_non_TA_cell_types_enriched_peaks_functional_BP.rds")
terms <- great_res_list[[1]]$name

# hit number matrix
hit_mat <- sapply(names(great_res_list), function(ct){
  great_res_list[[ct]]$Binom_Observed_Region_Hits
})
rownames(hit_mat) <- terms
saveRDS(hit_mat, file="Dat_salmonella_infection_on_non_TA_cell_types_enriched_peaks_functional_BP_hit_mat.rds")

# nominal binomial p-value matrix
pval_mat <- sapply(names(great_res_list), function(ct){
  great_res_list[[ct]]$Binom_Raw_PValue
})
rownames(pval_mat) <- terms
saveRDS(pval_mat, file="Dat_salmonella_infection_on_non_TA_cell_types_enriched_peaks_functional_BP_pval_mat.rds")

# BH-adjusted binomial p-value matrix
padj_mat <- sapply(names(great_res_list), function(ct){
  great_res_list[[ct]]$Binom_Adjp_BH
})
rownames(padj_mat) <- terms
saveRDS(padj_mat, file="Dat_salmonella_infection_on_non_TA_cell_types_enriched_peaks_functional_BP_padj_mat.rds")

fold_mat <- sapply(names(great_res_list), function(ct){
  great_res_list[[ct]]$Binom_Fold_Enrichment
})
rownames(fold_mat) <- terms
saveRDS(fold_mat, file="Dat_salmonella_infection_on_non_TA_cell_types_enriched_peaks_functional_BP_binom_fold_enrichment_mat.rds")

pval_mat <- readRDS("Dat_salmonella_infection_on_non_TA_cell_types_enriched_peaks_functional_BP_pval_mat.rds")
fold_mat <- readRDS("Dat_salmonella_infection_on_non_TA_cell_types_enriched_peaks_functional_BP_binom_fold_enrichment_mat.rds")
enrichment_mat <- -log10(pval_mat)
enrichment_diff_mat <- enrichment_mat - rowMeans(enrichment_mat)
res <- lapply(seq(ncol(enrichment_diff_mat)), function(j){
  idx <- which(padj_mat[,j]<0.05 & enrichment_diff_mat[,j]>0)
  N <- min(c(10, length(idx)))
  values <- setNames(enrichment_diff_mat[idx,j], terms[idx])
  names(values)[order(values, decreasing=T)[1:N]]
})
names(res) <- colnames(pval_mat) 
union_terms <- setdiff(unique(unlist(res)), NA)
union_terms <- c(
  "positive regulation of T cell mediated cytotoxicity",
  "cell-cell junction maintenance",
  "negative regulation of entry of bacterium into host cell",
  "beta-catenin-TCF complex assembly",
  "positive regulation of cation transmembrane transport",
  "positive regulation of toll-like receptor signaling pathway",
  "positive regulation of cellular response to hypoxia",
  "endocytosis",
  "cellular response to carbon monoxide",
  "positive regulation of sodium:potassium-exchanging ATPase activity",
  "positive regulation of T-helper 2 cell cytokine production",
  "interleukin-10 secretion",
  "negative regulation of chemokine secretion",
  "positive regulation of neutrophil mediated killing of gram-negative bacterium",
  "regulation of cell death",
  "positive regulation of apoptotic signaling pathway",
  "protein modification by small protein conjugation or removal",
  "regulation of apoptotic process",
  "leukocyte activation",
  "positive regulation of cellular protein metabolic process",
  "negative regulation of lipopolysaccharide-mediated signaling pathway"
)

# visualize the enrichment of selected terms using dot plot
used_pval_mat <- -log10(padj_mat[union_terms, ])
used_fold_mat <- t(apply(fold_mat[union_terms, ],1,function(vec){(vec-min(vec))/(max(vec)-min(vec))}))
used_hit_mat <- hit_mat[union_terms, ]
df <- data.frame('BP_term'=rep(rownames(used_fold_mat), ncol(used_fold_mat)),
                 'cell_type'=rep(colnames(used_fold_mat), each=nrow(used_fold_mat)),
                 'pval'=as.vector(used_pval_mat),
                 'fold'=as.vector(used_fold_mat),
                 'hit'=as.vector(used_hit_mat),
                 stringsAsFactors = F)
saveRDS(df, file="Dat_term_padj_per_ct.rds")

library(ggplot2)
# order the condition
res <- prepareTreeAndHeatmapInput(expr.mat = used_fold_mat, hc.method = "ward.D2", norm.method = "quantile")
term_order <- rownames(res$heatmap_input)
cell_type_order <- colnames(res$heatmap_input)
df$BP_term <- factor(df$BP_term, levels=term_order)
df$cell_type <- factor(df$cell_type, levels=cell_type_order)

p1 <- ggplot(df, aes(x=cell_type, y=BP_term, size=fold, fill=pval))+
    geom_point(shape=21, col="#202020")+
    scale_size(range=c(1,10))+
    scale_fill_gradientn(colors=beach.col)+
    theme(text = element_text(size = 30))+
    theme_light()+
    guides(x =  guide_axis(angle = 45))+
    labs(x="Cell type", y="Biological process", fill="-log10[BH-P]", size="Binom fold enrichment")
p1
ggsave(p1, filename="Plot_dot_selected_term_enrichment_2.pdf", width=10, height=10)

term <- "regulation of apoptotic process"
ct <- "EEC_Salmonella"
great_res_list[[ct]] -> res
res[which(res$name==term),]
ct2 <- "Colon_SC_Salmonella" 
great_res_list[[ct2]] -> res2
res2[which(res2$name==term),]

# check motif enrichment 


# generate coverage plot
setwd("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/coverage_plot/")
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

peak_anno <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/peak_annotation/Data_frame_peak_gene_pairs.rds")
peaks <- grn_marker_dar
colors <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/cols.rds")
group <- sort(unique(combined_rna$Coarse_cell_type_per_condition))
mat <- do.call('rbind', strsplit(group, split=":"))
ct <- mat[,1]
names(ct) <- group
group_cols <- setNames(colors$cell_type[ct], group)
setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/coverage_plot")
source("/home/yuq22/ihb-intestine-evo/common_script/plotting/Script_plotting_functions.R")
p_list <- coverage_plot_with_evo_sig(seu_obj=combined_rna, 
                                     peak_assay="peaks_merged", 
                                     rna_assay="RNA", 
                                     group="Coarse_cell_type_per_condition", 
                                     colors.use=group_cols, 
                                     peaks=peaks, 
                                     evo_list_path="/projects/site/pred/ihb-intestine-evo/evo_signature/summary_Feb_2024/Dat_SNC_GWAS_HAR_PS_HAQER.rds", 
                                     peak_anno_path="/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/peak_annotation/Data_frame_peak_gene_pairs.rds", 
                                     peak_anno_peak_col="name",
                                     peak_anno_gene_col="external_gene_name",
                                     extend_up=10000,
                                     extend_down=10000,
                                     window_size=500,
                                     color="#303030",
                                     gene_mode="ChIPseeker", 
                                     g_vec=NULL,
                                     do_plot_combined=FALSE,
                                     do_plot_individual=TRUE,
                                     combined_plot_name=NULL,
                                     plot_suffix="-2")

# get coarse cell type marker region
# generate coverage plot on ATAC cluster levels
setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/RNA_coarse_cell_type_marker_regions")
# identify infection-specific colonocyte subtype markers only using the salmonella infected data, not the control data
sal_data <- subset(combined, subset=condition=="Salmonella")
de_res <- presto::wilcoxauc(sal_data, group_by="RNA_coarse_cell_type", seurat_assay="peaks_merged")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
top_res <- sig_res %>% group_by(group) %>% top_n(200, wt=logFC)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="../Res_colon_epithelium_Salmonella_infection_only_dataset_coarse_cell_type_marker_regions.rds")
peaks <- sig_res %>% filter(group%in%c("Infection_specific_colonocyte_subtype_1", "Infection_specific_colonocyte_subtype_2", "Non_infection_specific_colonocyte")) %>% group_by(group) %>% top_n(5, wt=pct_diff) %>% pull(feature) %>% unique()

# generate coverage plot of canonical cell type marker regulatory regions
genes <- c("MUC2", "LGR5", "CHGA", "MEIS1", "CEACAM7", "CA1", "SMOC2", "ASCL2")
# specific genes
genes <- c("LTB", "IL1B", "CXCL8")
# gene associated with specific pathways
kegg_gene_list <- readRDS("/home/yuq22/Bacteria_TRM/Annotation/KEGG/Res_KEGG_gene_list.rds")
selected_pathways <- c("Antigen processing and presentation", "Salmonella infection", "Bacterial invasion of epithelial cells", "Toll-like receptor signaling pathway", "NF-kappa B signaling pathway", "Regulation of actin cytoskeleton", "Tight junction")
pathway_genes <- unique(unlist(kegg_gene_list[selected_pathways]))
deg_sig_res <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_matched_control_cells_and_cocluster_control_cells_DEG_sig_res.rds")
deg_top_res <- deg_sig_res %>% group_by(group) %>% top_n(200, wt=pct_diff)
genes <- intersect(pathway_genes, deg_top_res$feature)
length(genes)

# intersect with significant DARs
dar_sig_res <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_combined_control_cells_DAR_sig_res_on_RNA_cell_type_level.rds")
dar_top_res <- dar_sig_res %>% group_by(group) %>% top_n(200, wt=logFC)
union_dar <- sort(unique(dar_top_res$feature))
peaks <- intersect(unique(peak_anno$name[which(peak_anno$external_gene_name%in%genes)]), union_dar)
length(peaks)
#peaks <- intersect(peak_anno$name[which(peak_anno$external_gene_name%in%genes)], sig_res$feature)
#top_res <- sig_res %>% filter(group%in%c("Infection_specific_colonocyte_subtype_1", "Infection_specific_colonocyte_subtype_2", "Non_infection_specific_colonocyte")) %>% group_by(group) %>% top_n(5, wt=logFC)
#peaks <-  unique(top_res$feature)
sal_rna <- subset(combined_rna, subset=condition=="Salmonella")
p_list <- coverage_plot_with_evo_sig(#seu_obj=combined_rna, 
                                     seu_obj=sal_rna,
                                     peak_assay="peaks_merged", 
                                     rna_assay="RNA", 
                                     #group="Coarse_cell_type_per_condition", 
                                     group="RNA_subtype_group_per_condition",
                                     #colors.use=group_cols, 
                                     colors.use=NULL,
                                     peaks=peaks, 
                                     evo_list_path="/projects/site/pred/ihb-intestine-evo/evo_signature/summary_Feb_2024/Dat_SNC_GWAS_HAR_PS_HAQER.rds", 
                                     peak_anno_path="/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/peak_annotation/Data_frame_peak_gene_pairs.rds", 
                                     peak_anno_peak_col="name",
                                     peak_anno_gene_col="external_gene_name",
                                     extend_up=10000,
                                     extend_down=10000,
                                     window_size=500,
                                     color="#303030",
                                     gene_mode="ChIPseeker", 
                                     g_vec=NULL,
                                     do_plot_combined=FALSE,
                                     do_plot_individual=TRUE,
                                     combined_plot_name=NULL,
                                     plot_suffix=NULL)
 