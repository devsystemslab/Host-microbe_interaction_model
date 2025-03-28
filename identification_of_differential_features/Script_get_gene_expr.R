setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/C16/gene_expr")
library(Seurat)
library(Signac)
library(dplyr)  
library(dplyr)
library(ggrepel)
library(Pando)
data(motif2tf)
source("~/ihb-intestine-evo/common_script/Script_functions.R")
cols <- readRDS("~/Bacteria_TRM/used_object/cols.rds")


combined_rna <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/motif_analysis/Res_salmonella_scMultiome_with_chromVAR_and_motif_res.rds")
all <- combined_rna
size <- table(combined_rna$RNA_subtype_group_per_condition)
combined_subset <- subset(combined_rna, subset = RNA_subtype_group_per_condition %in% names(size)[which(size>20)] & Coarse_cell_type != "TA")
combined_rna <- combined_subset

# plot gene expression logFC in C16 and C14/17
# load infection vs control DEG test result
combined_de_res <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/Res_Salmonella_infection_on_colon_epi_goblet_and_s2E_all_DEG_res.rds")

# get the fold change of infection - control
all_fc_data <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/NFKB_barrier_innate_sensing/combined_with_Ruben_gene_set/Dat_expr_logFC_infection_control.rds")
all_pval_data <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/NFKB_barrier_innate_sensing/combined_with_Ruben_gene_set/Dat_expr_pval_infection_control.rds")

groups <- grep("infection_specific", colnames(all_fc_data), value=T)
pairs <- combn(groups, 2)

plot_list <- list()
for(i in seq(ncol(pairs))){
    df <- data.frame(
        "gene"=rownames(all_fc_data),
        all_fc_data[,pairs[,i]], 
        stringsAsFactors = F
    )
    plot_list[[i]] <- ggplot(df, aes_string(x=pairs[1,i], y=pairs[2,i])) + 
        geom_point(pch=16) + 
        geom_abline(intercept=0, slope=1) + 
        theme_bw() + 
        ggtitle(paste(pairs[,i], collapse=" vs "))
}
p_combined <- cowplot::plot_grid(plotlist=plot_list, ncol=3)
ggsave(p_combined, file="plot_gene_expr_logFC_C16_vs_C14_C17.png", width=30, height=10, dpi=300)
# C14 correlates with C17, C16 is distinct
all_cell_types <- sapply(colnames(all_fc_data), function(x){
    vec <- strsplit(x, "_")[[1]]
    paste(vec[-c(length(vec)-1, length(vec))], collapse="_")
})
fc_by_cell_type <- sapply(unique(all_cell_types), function(ct){
    idx <- which(all_cell_types==ct)
    if(length(idx)>1){
        return(apply(all_fc_data[,idx], 1, function(vec){ vec[which.max(abs(vec))] }))
    }else{
        return(all_fc_data[,idx])
    }
})
saveRDS(fc_by_cell_type, file="Dat_gene_expr_logFC_per_cell_type_data_only.rds")
df <- data.frame(
    'gene'=rownames(all_fc_data),
    fc_by_cell_type,
    'tf'=rownames(all_fc_data)%in%motif2tf$tf,
    stringsAsFactors = F
)
saveRDS(df, file="Dat_gene_expr_logFC_per_cell_type.rds")


# cluster the genes based on the average expression level per cell type per condition
combined_rna <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/motif_analysis/Res_salmonella_scMultiome_with_chromVAR_and_motif_res.rds")
DefaultAssay(combined_rna) <- "RNA"
combined_rna$RNA_cluster_group <- combined_rna$Subtype
for(i in c(14, 16, 17)){
    combined_rna$RNA_cluster_group[which(combined_rna$RNA_snn_res.1==i)] <- paste("Colonocyte_2_infection_specific_C", i, sep="")
}
combined_rna$RNA_cluster_group_per_condition <- paste(combined_rna$RNA_cluster_group, combined_rna$condition, sep=":")
combined_rna$RNA_subtype_group_per_condition <- paste(combined_rna$RNA_subtype_group, combined_rna$condition, sep=":")
saveRDS(combined_rna, file="/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/motif_analysis/Res_salmonella_scMultiome_with_chromVAR_and_motif_res_v2.rds")
# only include goblet cells, stem cells and colonocytes
combined_sub <- subset(combined_rna, subset = Coarse_cell_type %in% c("Goblet", "Colon_SC", "Colonocyte"))
saveRDS(combined_sub, file="/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/motif_analysis/Res_salmonella_scMultiome_with_chromVAR_and_motif_res_v2_SC_GC_Colonocyte.rds")
ave_expr <- getAveExpr(seu.obj=combined_sub, feature.to.calc="RNA_cluster_group_per_condition", colname.prefix=NULL, size.cutoff=20)
saveRDS(ave_expr, file="Dat_salmonella_dataset_GC_SC_Colonocyte_RNA_cluster_group_per_condition_ave_expr.rds")

combined_sig_de_res <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/Res_Salmonella_infection_on_colon_epi_goblet_and_s2E_DEG_sig_res.rds") 
union_deg <- unique(combined_sig_de_res$feature)
# cluster based on the logFC so that the clustering is not driven by cell type differentiation, but plot the average expression level per cell type per condition to ease the interpretation
fc_by_cell_type <- readRDS("Dat_gene_expr_logFC_per_cell_type_data_only.rds")
input_expr_mat <- fc_by_cell_type[union_deg,]
hc <- hclust(as.dist(1-cor(t(input_expr_mat))), method="ward.D2")
plot(hc, hang=-1)
h=10
abline(h=h)
g <- cutree(hc, h=h)
saveRDS(g, file="Res_Salmonella_infection_DEG_cluster.rds")

#scaled_fc_data <- t(scale(t(input_expr_mat)))
#saveRDS(scaled_fc_data, file="Dat_DEG_scaled_expr_logFC_per_cell_type.rds")
#scaled_fc_data <- readRDS("Dat_DEG_scaled_expr_logFC_per_cell_type.rds")
#g <- readRDS("Res_Salmonella_infection_DEG_cluster.rds")
col_num=4
row_num=ceiling(length(unique(g))/col_num)

ave_expr <- readRDS("Dat_salmonella_dataset_GC_SC_Colonocyte_RNA_cluster_group_per_condition_ave_expr.rds")
mat <- do.call('rbind', strsplit(split=":", colnames(ave_expr)))
group_cols <- setNames(c("C0", "00"), c("Salmonella", "Control"))

cols$cell_type["Colonocyte_1"] <- "#ffeda0"
box_cols <- paste0(cols$cell_type[mat[,1]], group_cols[mat[,2]])
border_cols <- cols$cell_type[mat[,1]]
scaled_expr_data <- t(scale(t(ave_expr)))
#colnames(scaled_fc_data)[grep("specific", colnames(scaled_fc_data))] <- sub("Colonocyte_2_infection_specific", "Infection", colnames(scaled_fc_data)[grep("specific", colnames(scaled_fc_data))])
pdf("Plot_boxplot_gene_expr_per_cell_type_per_condition.pdf", width=6*col_num, height=5*row_num)
par(mfrow=c(row_num,col_num), mar=c(10,5,5,5))
for(i in seq(length(unique(g)))){
    boxplot(scaled_expr_data[names(g)[which(g==i)],], las=2, main=paste("Module", i, "(", sum(g==i), ")"), border=border_cols, col=box_cols, outline=F, lwd=2)
}
dev.off()

########
# get the marker gene overlap with feature gene sets (KEGG, HALLMARK)
# marker genes overlap with gene sets
path <- "/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/more_DEG_characterization/"
gene_list <- readRDS(paste(path, "Res_KEGG_and_HALLMARK_gene_list.rds", sep="/"))

all_genes <- sort(unique(unlist(gene_list)))
gene_mat <- sapply(names(gene_list), function(gs){
    gene_set <- gene_list[[gs]]
    all_genes %in% gene_set
})
rownames(gene_mat) <- all_genes
saveRDS(gene_mat, file="Dat_KEGG_HALLMARK_gene_annotation.rds")
saveRDS(gene_mat, )
all_expressed_genes <- readRDS(paste(path, "Dat_s2E_and_goblet_all_expressed_genes.rds", sep="/"))
all_overlap_genes <- intersect(unlist(gene_list), all_expressed_genes)

g <- readRDS("Res_Salmonella_infection_DEG_cluster.rds")
# check enrichment on infection induced DEGs on stem cell, colonocyte and goblet cell
pval <- sapply(sort(unique(g)), function(i){
    sapply(names(gene_list), function(gs){
        gene_set <- gene_list[[gs]]
        marker_genes <- names(g)[which(g==i)]
        overlap_genes <- intersect(intersect(gene_set, marker_genes), all_overlap_genes)
        marker_non_gs_genes <- intersect(setdiff(marker_genes, gene_set), all_overlap_genes)
        gs_non_marker_genes <- intersect(setdiff(gene_set, marker_genes), all_overlap_genes)
        non_gs_non_marker_genes <- setdiff(all_overlap_genes, union(gene_set, marker_genes))
        a <- length(overlap_genes)
        b <- length(marker_non_gs_genes)
        c <- length(gs_non_marker_genes)
        d <- length(non_gs_non_marker_genes)
        fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="g")$p.value
    })
})
colnames(pval) <- paste0("M", sort(unique(g)))
saveRDS(pval, file="Res_DEG_cluster_and_gene_set_overlap_enrichment_pval.rds")
write.table(pval, file="Table_Salmonella_DEG_cluster_and_gene_set_overlap_enrichment_nominal_pval.txt", sep="\t", quote=F)

# result visualization
# manually seleted terms
manual_terms <- c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
    "HALLMARK_HYPOXIA", 
    "Adherens junction", 
    "HALLMARK_GLYCOLYSIS", 
    "Bacterial invasion of epithelial cells", 
    "Tight junction", 
    "Fructose and mannose metabolism", 
    "HALLMARK_TGF_BETA_SIGNALING", 
    "HIF-1 signaling pathway", 
    "HALLMARK_APOPTOSIS", 
    "Pathogenic Escherichia coli infection", 
    "Leukocyte transendothelial migration", 
    "NF-kappa B signaling pathway", 
    "Protein processing in endoplasmic reticulum", 
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", 
    "HALLMARK_MTORC1_SIGNALING", 
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", 
    "Shigellosis")
# de novo top terms of M1, M3 and M6 where genes exhibit highest gene expression levels in the Salmonella specific clusters 
selected_clusters <- paste0("M", c(1, 3, 6))
# exclude the cancer and disease related terms
data <- -log10(pval[!grepl("cancer|carcino|disease|muscle|Drug|cardio|neurode|Amyotrophic", rownames(pval)), selected_clusters])
top_terms <- lapply(colnames(data), function(x){
    intersect(rownames(data)[order(data[,x], decreasing=T)[1:10]], rownames(data)[data[,x]> -log10(0.05)])
})
names(top_terms) <- colnames(data)
saveRDS(top_terms, file="Res_DEG_KEGG_hypergeometric_test_top_enriched_terms_for_Sal_specific_cluster_high_genes.rds")

terms <- sort(union(unique(unlist(top_terms)), manual_terms))
res <- prepareTreeAndHeatmapInput(expr.mat = -log10(pval[terms,]), hc.method = "ward.D2", norm.method = "quantile")
saveRDS(res, file="Res_DEG_KEGG_hypergeometric_test_heatmap_input.rds")

black_cols <- c("#f7f7f7", "#cccccc", "#969696", "#636363", "#252525")
black_cols_heatmap <- colorRampPalette(black_cols)(30) 
pdf("Plot_heatmap_kegg_pathway_enrichment.pdf", height=10)
gplots::heatmap.2(res$heatmap_input, col = black_cols_heatmap, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 1, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

########
# check which genes are involved in Ferroptosis
# check whether they are target genes of STAT3 or SOCS3

####
## visualize the functional enrichment of the DEGs 
pval <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/more_DEG_characterization/Res_cell_group_DEG_and_gene_set_overlap_enrichment_pval.rds")
# focus on salmonella infection upregulated genes
input <- -log10(pval[manual_terms,grep("Salmonella",colnames(pval))])
# generate barplot to show infection upregulated genes enriched terms
enrichment <- apply(input, 1, max)
pdf("Plot_barplot_selected_term_enrichment_v2.pdf", height=5)
par(mar=c(5,20,2,2))
barplot(sqrt(sort(enrichment)), horiz=T, las=2, xlab="sqrt(-log10(p-value))")
abline(v=sqrt(-log10(0.05)), lty=2)
dev.off()













combined_sig_de_res <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/Res_Salmonella_infection_on_colon_epi_goblet_and_s2E_DEG_sig_res.rds") 
union_deg <- unique(combined_sig_de_res$feature)
hc <- hclust(as.dist(1-cor(fc_by_cell_type[union_deg, ])), method="ward.D2")
plot(hc, hang=-1)

pairs <- combn(colnames(fc_by_cell_type), 2)
cor_vec <- sapply(1:ncol(pairs), function(i){
    cor(df[,pairs[1,i]], df[,pairs[2,i]], use="complete.obs")
})
names(cor_vec) <- apply(pairs, 2, paste, collapse=" vs ")
plot(seq(length(cor_vec)), sort(cor_vec), pch=16)

plot_list <- list()
for(i in seq(ncol(pairs))){
    df_i <- df[,c("gene", pairs[,i])]
    plot_list[[i]] <- ggplot(df_i, aes_string(x=pairs[1,i], y=pairs[2,i])) + 
        geom_point(pch=16) + 
        geom_abline(intercept=0, slope=1) + 
        theme_bw() + 
        ggtitle(paste(pairs[,i], collapse=" vs "))
}
p_combined <- cowplot::plot_grid(plotlist=plot_list, ncol=6)
ggsave(p_combined, file="plot_gene_expr_logFC_per_cell_type_pairs.png", width=5*6, height=5*6, dpi=300)


# get the genes with the largest logFC
t1_genes_to_highlight<- df$gene[order(df$subtype1_logFC, decreasing=T)[1:10]]
t2_genes_to_highlight<- df$gene[order(df$subtype2_logFC, decreasing=T)[1:10]]
genes_to_highlight <- union(t1_genes_to_highlight, t2_genes_to_highlight)
p1 <- ggplot(df, aes(x=subtype1_logFC, y=subtype2_logFC)) + 
        geom_point(pch=16) + 
        geom_text_repel(data=df[which(df$gene %in% genes_to_highlight),], aes(label=gene), col="red", size=2, box.padding=0.5)+
        geom_hline(yintercept=0, linetype="dashed") +
        geom_vline(xintercept=0, linetype="dashed") +
        theme_bw()
ggsave(p1, file="plot_TF_gene_expr_logFC_C16_vs_C14_C17_highlighted_genes.png", width=5, height=5, dpi=300)



df <- readRDS("Dat_gene_expr_logFC_C16_vs_C14_C17.rds")
DefaultAssay(combined_rna) <- "RNA"
sum(df[,"subtype2_logFC"]>df["PLA2G2A", "subtype2_logFC"])
SCpubr::do_BoxPlot(combined_rna, feature="PLA2G2A", group.by="RNA_subtype_group_per_condition")

par(mar=c(25,5,5,5))
g <- "SOCS3"
barplot(ave_expr[g,], las=2, main=g, ylab="Normalized expression")
genes_to_highlight <- c("PLA2G2A", "DMBT1", "STAT3", "SOCS3", "BHLHE40", "NDRG1", "DUOX2", "NFKB1", "NKFB2")
df_tf <- df[which(df$tf),]

# plot the logFC of DE-TFs
combined_sig_de_res <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/Res_Salmonella_infection_on_colon_epi_goblet_and_s2E_DEG_sig_res.rds") 
union_deg <- unique(combined_sig_de_res$feature)
union_de_tf <- intersect(union_deg, motif2tf$tf)
all_fc_data <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/NFKB_barrier_innate_sensing/combined_with_Ruben_gene_set/Dat_expr_logFC_infection_control.rds")
tf_fc_data <- all_fc_data[union_de_tf,]
# get the strongest logFC and the cell type with the strongest logFC
res <- t(apply(tf_fc_data, 1, function(vec){
    idx <- which.max(abs(vec))
    return(c(vec[idx], colnames(tf_fc_data)[idx]))
}))
df <- data.frame(
    "gene"=rownames(tf_fc_data),
    "logFC"=as.numeric(res[,1]),
    "cell_cluster"=res[,2],
    stringsAsFactors = F
)
df$cell_type <- sapply(df$cell_cluster, function(x){
    vec <- strsplit(x, "_")[[1]]
    paste(vec[-c(length(vec)-1, length(vec))], collapse="_")
})
df <- df[order(df$logFC),]
df$gene <- factor(df$gene, levels=df$gene)
saveRDS(df, file="Dat_DE_TF_max_logFC.rds")


ggplot(data=df, aes(x = gene, y = logFC, fill=cell_type))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = cols$cell_type)+
    theme_minimal() +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



# get SOCS3+ cells %
c16_cells <- colnames(combined_rna)[which(combined_rna$RNA_snn_res.1==16)]
all_cell_types <- sort(unique(combined_rna$RNA_subtype_group_per_condition))
prop_mat <- sapply(0:5, function(cutoff){
    socs3_pos_cells <- colnames(combined_rna)[which(combined_rna@assays$RNA@counts["SOCS3", ]>cutoff)]
    prop <- sapply(all_cell_types, function(ct){
        cell0 <- colnames(combined_rna)[which(combined_rna$RNA_subtype_group_per_condition==ct)]
        cell1 <- intersect(cell0, socs3_pos_cells)
        length(cell1)/length(cell0)
    })
    return(prop)
})

plot(seq(nrow(prop_mat)), prop_mat[,1], type = "n",  xaxt = "n", xlab="", ylab="Proportion of SOCS3+ cells")
axis(1, at = seq(nrow(prop_mat)), labels =rownames(prop_mat), las=2)
for(j in seq(ncol(prop_mat))){
    points(seq(nrow(prop_mat)), prop_mat[,j], col=rainbow(ncol(prop_mat))[j], pch=16)
    lines(seq(nrow(prop_mat)), prop_mat[,j], col=rainbow(ncol(prop_mat))[j])
}
legend("topleft", legend=paste("count cutoff=", 0:5), col=rainbow(ncol(prop_mat)), pch=16)
