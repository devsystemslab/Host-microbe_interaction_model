setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data")
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(rGREAT)
library(Signac)
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")

data("motifs")

combined_rna <- readRDS("Dat_filtered_colon_epithelium_infection_RNA_ATAC_preprocessed_ctByCondition.rds")
combined_de_res <- readRDS("Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_combined_control_cells_DAR_all_res_on_RNA_cell_type_level.rds")
sal_up_200 <- combined_de_res[grepl("Salmonella", combined_de_res$group) & !grepl("TA", combined_de_res$group),] %>% filter(logFC>0) %>% group_by(group) %>% top_n(200, wt=logFC)

setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/motif_analysis")
pdf("Plot_condition_logFC_per_cell_type.pdf", height=10, width=20)
par(mfrow=c(2,4))
for(x in unique(sal_up_200$group)){
    hist(sal_up_200$logFC[sal_up_200$group==x], main=x, xlab="logFC")
}
dev.off()

DefaultAssay(combined_rna) <- "peaks_merged"
combined_rna <- AddMotifs(
  object = combined_rna,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = motifs
)
saveRDS(combined_rna, file="Res_salmonella_scMultiome_with_motif_res.rds")

# takes larger than 100GB memory, so submit a job
combined_rna <- RunChromVAR(
  object = combined_rna,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
saveRDS(combined_rna, file="Res_salmonella_scMultiome_with_chromVAR_and_motif_res.rds")

setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/motif_analysis")
combined_rna <- readRDS("Res_salmonella_scMultiome_with_chromVAR_and_motif_res.rds")

module_score_data <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/NFKB_barrier_innate_sensing/Dat_salmonella_dataset_all_epi_cell_type_GRN_module_score.rds")
combined_rna$cluster_group_per_condition <- module_score_data$cluster_group_per_condition
# get cell type average chromVAR score per condition
X <- combined_rna@assays$chromvar@data
y <- combined_rna$cluster_group_per_condition
ave_score <- getAveExpr(input.type="matrix", X=X, y=y, size.cutoff=20, colname.prefix=NULL)
saveRDS(ave_score, file="Dat_salmonella_dataset_all_epi_cell_type_chromVAR_ave_score_per_condition.rds")

idx <- setdiff(intersect(grep("Colonocyte", colnames(ave_score)), grep("Salmonella", colnames(ave_score))), grep("infection_specific", colnames(ave_score)))
colonocyte_sal_max_score <- apply(ave_score[,idx], 1, max)

idx <- grep("infection_specific", colnames(ave_score))
colonocyte_sal_specific_max_score <- apply(ave_score[,idx], 1, max)

idx <- intersect(grep("Colonocyte", colnames(ave_score)), grep("Control", colnames(ave_score)))
colonocyte_con_max_score <- apply(ave_score[,idx], 1, max)
ct_ave_score <- cbind(ave_score[,!grepl("Colonocyte|TA", colnames(ave_score))], colonocyte_con_max_score, colonocyte_sal_max_score, colonocyte_sal_specific_max_score)
colnames(ct_ave_score)[c((ncol(ct_ave_score)-2):ncol(ct_ave_score))] <- c("Colonocyte_Control", "Colonocyte_Salmonella", "Colonocyte_Salmonella_specific")
saveRDS(ct_ave_score, file="Dat_salmonella_dataset_all_epi_cell_type_max_chromVAR_score_per_condition.rds")

# convert the rownames from motif ID to TF name
library(Pando)
data('motif2tf')
rownames(ct_ave_score) <- sub("-", "_", rownames(ct_ave_score))
pairs <- sapply(unique(motif2tf$motif), function(x) {
    paste(motif2tf$tf[motif2tf$motif==x], collapse=",")
})
rownames(ct_ave_score) <- pairs[rownames(ct_ave_score)]
rownames(X) <- sub("-", "_", rownames(X))
rownames(X) <- pairs[rownames(X)]
saveRDS(X, file="Dat_salmonella_dataset_all_epi_cell_type_chromVAR_score_per_cell_with_TF_name.rds")
saveRDS(ct_ave_score, file="Dat_salmonella_dataset_all_epi_cell_type_max_chromVAR_score_per_condition_with_TF_name.rds")
rownames(ave_score) <- sub("-", "_", rownames(ave_score))
rownames(ave_score) <- pairs[rownames(ave_score)]
saveRDS(ave_score, file="Dat_salmonella_dataset_all_epi_cell_type_chromVAR_ave_score_per_condition_with_TF_names.rds")

# get motifs with significant variance explained condition on top of cell type
res <- t(sapply(combined_rna$cluster_group_per_condition, function(x){
  vec <- strsplit(x, "_")[[1]]
  c(paste(vec[1:(length(vec)-1)], collapse="_"), vec[length(vec)])
}))
module_score_data$cluster_group <- res[,1]
module_score_data$condition <- res[,2]
saveRDS(module_score_data, file="/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/NFKB_barrier_innate_sensing/Dat_salmonella_dataset_all_epi_cell_type_GRN_module_score.rds")

library(doParallel)
registerDoParallel(20)
ct_vec <- as.factor(res[,1])
condition_vec <- as.factor(res[,2])
p.value <- foreach(gene = rownames(X), .combine = c) %dopar% {
            e <- as.vector(as.matrix(X[gene, ]))
            m0 <- lm(e ~ ct_vec)
            m1 <- lm(e ~ ct_vec+condition_vec)
            a0 <- anova(m0)
            a1 <- anova(m1)
            p_anova <- anova(m1,m0)$Pr[2]
            return(p_anova)
        }
names(p.value) <- rownames(X)
stopImplicitCluster()
saveRDS(p.value, file="Res_condition_dependent_chromVAR_score_p_value.rds")

padj <- p.adjust(p.value, method="BH")

# intersect with DE-TFs
epi_sig_res <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/Res_Salmonella_infection_on_colon_epi_goblet_and_s2E_DEG_sig_res.rds")
de_tfs <- sort(intersect(motif2tf$tf, epi_sig_res$feature[which(epi_sig_res$condition=="Salmonella" & epi_sig_res$logFC>0.15)]))
de_tf_chromvar_sig_diff <- intersect(de_tfs, unlist(strsplit(split=",", names(padj)[padj<0.05])))
saveRDS(de_tf_chromvar_sig_diff, file="Res_DE_TFs_with_sig_chromVAR_score_diff_between_condition.rds")
# generate heatmap showing chromVAR score overlapping with the DE-TFs
idx <- sapply(rownames(ave_score), function(x){
  length(intersect(strsplit(x, split=",")[[1]], de_tf_chromvar_sig_diff))>0 
})
names(idx) <- rownames(ave_score)

input <- ct_ave_score[idx,]
saveRDS(input, file="Dat_cell_type_condition_average_chromVAR_score_per_DE_TF_per_motif.rds")

res <- prepareTreeAndHeatmapInput(expr.mat=input, hc.method="ward.D2")
saveRDS(res, file="Res_chromVAR_score_heatmap_input.rds")

pdf("Plot_chromVAR_cell_type_score_heatmap_2.pdf", width=7, height=10)
gplots::heatmap.2(res$heatmap_input, col=beach.col.heatmap, scale="none", trace="none", Rowv=FALSE, Colv=FALSE, cexRow=0.5, cexCol=1, key=0.5, dendrogram="none", margin=c(10,5))
dev.off()

# calculate the correlation between expression logFC and chromVAR score difference induced upon Salmonella infection 
# get chromVAR score difference
chromvar_diff <- sapply(sort(setdiff(unique(module_score_data$cluster_group), "TA")), function(ct){
  if(grepl("infection_specific", ct)){
    ave_score[,paste(ct,"Salmonella",sep="_")]
  }else{
    ave_score[,paste(ct,"Salmonella",sep="_")] - ave_score[,paste(ct,"Control",sep="_")]
  }
})
chromvar_diff <- data.frame(chromvar_diff, gene_symbol=rownames(chromvar_diff))
chromvar_diff$feature <- rownames(chromvar_diff)
chromvar_diff$motif_ID <- names(pairs)
saveRDS(chromvar_diff, file="Dat_chromVAR_score_diff_per_gene_per_motif.rds")

# get regulomes with significant variance explained condition on top of cell type
regulome_data <-  t(module_score_data[,grep("pos|neg",colnames(module_score_data))])
grn_p_value <- foreach(gene = rownames(regulome_data), .combine = c) %dopar% {
            e <- as.vector(as.matrix(regulome_data[gene, ]))
            m0 <- lm(e ~ ct_vec)
            m1 <- lm(e ~ ct_vec+condition_vec)
            a0 <- anova(m0)
            a1 <- anova(m1)
            p_anova <- anova(m1,m0)$Pr[2]
            return(p_anova)
        }
names(grn_p_value) <- rownames(regulome_data)
stopImplicitCluster()
saveRDS(grn_p_value,file="Res_condition_dependent_GRN_module_score_p_value.rds")

# perform motif enrichment analysis on top DARs
setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/motif_analysis")
combined_de_res <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_combined_control_cells_DAR_all_res_on_RNA_cell_type_level.rds")
sal_up_200 <- combined_de_res[grepl("Salmonella", combined_de_res$group) & !grepl("TA", combined_de_res$group),] %>% filter(logFC>0) %>% group_by(group) %>% top_n(200, wt=logFC)
combined_rna <- readRDS("Res_salmonella_scMultiome_with_chromVAR_and_motif_res.rds")

# find peaks open in each cell type
Idents(combined_rna) <- combined_rna$RNA_subtype_group
open_peak_list <- lapply(sort(unique(combined_rna$RNA_subtype_group)), function(x){
  AccessiblePeaks(combined_rna, idents = x)
})
names(open_peak_list) <- sort(unique(combined_rna$RNA_subtype_group))
saveRDS(open_peak_list, file="Dat_open_peaks_per_cell_type.rds")
# match the overall GC content in the peak set
meta_feature <- combined_rna@assays$peaks_merged@meta.features
matched_peak_list <- lapply(sub("_Salmonella", "", sort(unique(sal_up_200$group))), function(x){
  top_da_peaks <- sal_up_200$feature[grep(x, sal_up_200$group, fixed=T)]
  peaks_matched <- MatchRegionStats(
    meta.feature = meta_feature[open_peak_list[[x]], ],
    query.feature = meta_feature[top_da_peaks, ],
    n = 30000
  )
  return(peaks_matched)
})
names(matched_peak_list) <- sort(unique(sal_up_200$group))
saveRDS(matched_peak_list, file="Dat_GC_content_matched_with_top_sal_up_open_peaks_per_cell_type.rds")

enriched_motif_res <- lapply(names(matched_peak_list), function(x){
  top_da_peaks <- sal_up_200$feature[which(sal_up_200$group==x)]
  enriched_motifs <- FindMotifs(
    object = combined_rna,
    features = top_da_peaks,
    assay="peaks_merged",
    background = matched_peak_list[[x]],
  )
  return(enriched_motifs)
})
names(enriched_motif_res) <- names(matched_peak_list)
saveRDS(enriched_motif_res, file="Res_sal_up_DAR_enriched_motifs_per_cell_type.rds")

# generate scatter plot for each cell type
# X axis showing overrepresenting of motifs in the top DARs
# Y axis showing the chromVAR score difference between Salmonella and control condition
# size represents logFC of corresponding TF (if one motif is associated with multiple TFs, the logFC with the largest absolute value is used)
# color represents condition change direction - red for up-regulated and blue for down-regulated upon Salmonella infection


# get cell type average chromVAR score per condition
X <- combined_rna@assays$chromvar@data
y <- combined_rna$RNA_subtype_group_per_condition
ave_score <- getAveExpr(input.type="matrix", X=X, y=y, size.cutoff=20, colname.prefix=NULL)
saveRDS(ave_score, file="Dat_salmonella_dataset_all_epi_RNA_subtype_group_chromVAR_ave_score_per_condition.rds")

# get chromVAR score difference
cell_types <- sort(setdiff(unique(combined_rna$RNA_subtype_group), "TA"))
chromvar_diff <- sapply(cell_types, function(ct){
  if(grepl("Infection_specific", ct)){
    ave_score[,paste(ct,"Salmonella",sep=":")]
  }else{
    ave_score[,paste(ct,"Salmonella",sep=":")] - ave_score[,paste(ct,"Control",sep=":")]
  }
})
rownames(chromvar_diff) <- sub("-", "_", rownames(chromvar_diff))
chromvar_diff <- data.frame(chromvar_diff, gene_symbol=as.character(pairs), motif_ID=rownames(chromvar_diff))
chromvar_diff$feature <- rownames(chromvar_diff)
colnames(chromvar_diff)[1] <- "BEST4+_cell"
saveRDS(chromvar_diff, file="Dat_chromVAR_score_diff_per_RNA_subtype_group_per_gene_per_motif.rds")

# get gene expression at the RNA subtype group level
ave_expr <- getAveExpr(seu.obj=combined_rna, feature.to.calc="RNA_subtype_group_per_condition", size.cutoff=20, colname.prefix=NULL)
saveRDS(ave_expr, file="Dat_salmonella_dataset_gene_expression_per_epi_RNA_subtype_group_per_condition.rds")

expr_diff <- sapply(cell_types, function(ct){
  if(grepl("Infection_specific", ct)){
    ave_expr[,paste(ct,"Salmonella",sep=":")]
  }else{
    ave_expr[,paste(ct,"Salmonella",sep=":")] - ave_expr[,paste(ct,"Control",sep=":")]
  }
})
saveRDS(expr_diff, file="Dat_gene_expr_logFC_per_RNA_subtype_group.rds")


input_list <- lapply(cell_types, function(ct){
  res1<- enriched_motif_res[[paste(ct, "Salmonella", sep="_")]]
  logp <- -log10(res1[rownames(chromvar_diff),"pvalue"])
  diff <- chromvar_diff[,ct]
  expr_fc <- sapply(chromvar_diff$gene_symbol, function(gene){
    genes <- intersect(unlist(strsplit(gene, split=",")), rownames(expr_diff))
    if(length(genes)==0){
      return(0)
    }else{
      vec <- expr_diff[genes, ct]
      vec[which.max(abs(vec))]
    }
  })

  df <- data.frame(
    'motif'=rownames(chromvar_diff), 
    'gene_symbol'=chromvar_diff$gene_symbol,
    'logp'=logp, 
    'chromvar_diff'=diff,
    'expr_logFC'=expr_fc,
    'abs_expr_logFC'=abs(expr_fc),
    stringsAsFactors=FALSE)
  return(df)
})
names(input_list) <- cell_types
saveRDS(input_list, file="Dat_motif_chromVAR_score_diff_TF_expr_logFC_DAR_motif_enrichment_per_RNA_subtype_group.rds")

input_list <- readRDS("Dat_motif_chromVAR_score_diff_TF_expr_logFC_DAR_motif_enrichment_per_RNA_subtype_group.rds")

library(ggrepel)
library(cowplot)
plot_list <- lapply(names(input_list), function(ct){
  df <- input_list[[ct]]
  # highlight motifs with significant motif enrichment and top expression logFC
  sig_df <- df[df$expr_logFC>0,]
  top_motif <- sig_df$motif[order(sig_df$logp, decreasing=TRUE)[1:5]]
  df$highlight <- FALSE
  df$highlight[match(top_motif, df$motif)] <- TRUE
  #df$highlight[grep("NR1D1|NFIL3|BHLHE40", df$gene_symbol)] <- TRUE
  ggplot(df, aes(x=logp, y=chromvar_diff, size=abs_expr_logFC, fill=expr_logFC)) + 
    geom_point(shape=21, color="#303030") + 
    scale_fill_gradient2(low="#053061", mid="#F7F7F7", high="#67001F", midpoint=0)+
    ggrepel::geom_text_repel(data=df[df$highlight,], aes(label=gene_symbol), size=4) +
    theme_minimal() + 
    labs(title=ct, x="-log10(Infection-up region motif enrichment p-value)", y="ChromVAR score difference", size="TF expression logFC", color="TF expression logFC")
})
p <- plot_grid(plotlist=plot_list, ncol=4, nrow=2)
ggsave(p, filename="Plot_motif_chromVAR_score_diff_expr_logFC_per_RNA_subtype_group.pdf", width=25, height=10)

ct <- "Infection_specific_colonocyte_subtype_2"
df <- input_list[[ct]]
# highlight motifs associated with selected TFs
selected_tfs <- c("STAT3", "BHLHE40")
df$highlight <- grepl(paste(selected_tfs, collapse="|"), df$gene_symbol)
#df$highlight[grep("NR1D1|NFIL3|BHLHE40", df$gene_symbol)] <- TRUE
p1 <- ggplot(df, aes(x=logp, y=chromvar_diff, size=abs_expr_logFC, fill=expr_logFC)) + 
  geom_point(shape=21, color="#303030") + 
  scale_fill_gradient2(low="#053061", mid="#F7F7F7", high="#67001F", midpoint=0)+
  ggrepel::geom_text_repel(data=df[df$highlight,], aes(label=gene_symbol), size=4) +
  theme_minimal() + 
  labs(title=ct, x="-log10(Infection-up region motif enrichment p-value)", y="ChromVAR score difference", size="TF expression logFC", color="TF expression logFC")
ggsave(p1, filename="Plot_DAR_motif_enrichment_highlight_STAT3.pdf")

enriched_motif_res <- readRDS("Res_sal_up_DAR_enriched_motifs_per_cell_type.rds")
enriched_motif_res_ct <- enriched_motif_res[["Infection_specific_colonocyte_subtype_2_Salmonella"]]
