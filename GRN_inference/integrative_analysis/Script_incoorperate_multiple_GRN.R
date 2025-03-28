setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/incoorperate_multiple_GRN/")

# load all epi inferred GRN
#epi_modules <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/Res_chip_infection_GRN_modules.rds")
#epi_modules_df <- data.frame(epi_modules@meta)

# load s2E inferred GRN
s2e_modules <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/GRN_on_subset/Analysis/Res_chip_infection_GRN_modules.rds")
s2e_modules_df <- data.frame(s2e_modules@meta)

# load per cell type inferred GRN
folder <- "/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/GRN_per_cell_type"
files <- list.files(folder, pattern="Res_chip_infection_.*_with_GRN.rds")
module_df_list <- list()
for(file in files){
    id <- sub("_with_GRN.rds", "", sub("Res_chip_infection_","",file))
    seu_obj <- readRDS(paste0(folder, "/", file))
    modules <- NetworkModules(seu_obj) 
    modules_df <- data.frame(modules@meta)
    module_df_list[[id]] <- modules_df
}

module_df_list[["epi"]] <- epi_modules_df
module_df_list[["s2e"]] <- s2e_modules_df

module_combined <- do.call(rbind, module_df_list)
saveRDS(module_combined, file="Res_combined_GRN_modules.rds")

ruben_genes <- list(
    "Cytokines"=c("CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL8", "CXCL16", "CCL20", "CCL28", "IL1A", "IL1B", "IL32"),
    "Cytokine_receptors"=c("CXCR4", "IL1RAP", "IL2RG", "IL4R", "IL6R", "IL6ST", "IL10RB", "IL17RA", "IL22RA1"),
    "NFKB"=c("REL", "RELB", "NFKB1", "NFKB2", "BIRC3", "LTB"),
    "Interferon_family"=c("IRF1", "IRF3", "IFNGR1", "IFNGR2"),
    "TNF"=c("TNF", "TNFAIP3", "TNIP1", "TNFSF15"),
    "Toll_like_receptor_pathway"=c("IRAK2", "TRAF6", "TLR2", "TLR3", "TLR4","TLR5"),
    "Antimicrobials_ROS"=c("LCN2", "SAA1", "SAA2", "DUOXA2", "DUOX2", "GPX2", "PIGR"),
    "Mucins_TFFs_glycosylation"=c("TFF1", "TFF2", "MUC1", "MUC13", "ST3GAL1", "ST3GAL4"),
    "Junctions_barrier_membrane"=c("ACTB", "ACTN4", "CLDN1", "CLDN23", "TJP1", "CEACAM1", "CEACAM5", "CEACAM7"),
    "Hypoxia_angiogenesis"=c("VEGFA","HIF1A","NDRG1"),
    "Others"=c("PLAUR", "PI3"),
    "Antigen_presentation"=c("CD74", "HLA-A", "B2M", "HLA-E", "TAPBP"),
    "Cytoskeleton"=c("ACTB", "CDC42", "RAC1")
)
genes <- unique(unlist(ruben_genes))

# exclude the modules coming from epi all cell types
module_combined <- do.call(rbind, module_df_list[-length(module_df_list)])
saveRDS(module_combined, file="Res_combined_GRN_modules_no_all_epi_cell_type_module.rds")

module_combined <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/incoorperate_multiple_GRN/Res_combined_GRN_modules_no_all_epi_cell_type_module.rds")
sig_module <- unique(module_combined[which(module_combined$padj<0.05 & module_combined$target%in%genes), c("tf", "target")])
# get target genes of STAT3 and BHLHE40
selected_sig_module <- unique(module_combined[which(module_combined$padj<0.05 & module_combined$tf %in% c("STAT3", "BHLHE40")), c("tf", "target", "estimate", "pval")])
#stat3_targets <- module_combined$target[which(module_combined$padj<0.05 & module_combined$tf=="STAT3")]



rna_sub <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/Dat_salmonella_colon_SC_colonocyte_goblet_scRNA-seq_data.rds")
p1 <- SCpubr::do_DotPlot(rna_sub, feature=c(stat1_targets,"STAT1"), group.by="Coarse_cell_type_per_condition", scale=T)
p1
filename="Plot_STAT1_target_gene_expression_per_cell_type_per_condition.pdf"
ggsave(p1, filename=filename, width=5, height=7)

gs <- c("ACTB", "ACTG1", "RAC1", "CDC42", "RHOA", "WAS", "ARPC1", "ARPC2", "ARPC3", "ARPC4", "ARPC5", "ACTR2", "ACTR3", "CFL1", "CFL2", "GSN", "PFN1")
p1 <- SCpubr::do_DotPlot(rna_sub, feature=gs, group.by="Coarse_cell_type_per_condition", scale=T)
p1
filename="Plot_selected_gene_expression_per_cell_type_per_condition.pdf"
ggsave(p1, filename=filename, width=5, height=7)

# for each target genes get the top-3 TFs showing the strongest correlation 
genes <- unique(unlist(ruben_genes))
sig_module <- unique(module_combined[which(module_combined$padj<0.05 & module_combined$target%in%genes), c("tf", "target")])
saveRDS(sig_module, file="Data_significant_TF-target_pairs_involving_ruben_genes.rds")

rna_sub <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/Dat_salmonella_colon_SC_colonocyte_goblet_scRNA-seq_data.rds")
ave_expr <- getAveExpr(seu.obj=rna_sub, feature.to.calc="cluster_group_per_condition", size.cutoff=20, colname.prefix=NULL, assay.type="RNA")
saveRDS(ave_expr, file="Data_ave_expr_per_cluster_group_per_condition.rds")


cor_vec <- sapply(seq(nrow(sig_module)), function(i){
    tf <- sig_module$tf[i]
    target <- sig_module$target[i]
    cor(ave_expr[tf,], ave_expr[target,])
})
sig_module$pcc <- cor_vec
saveRDS(sig_module, file="Data_significant_TF-target_pairs_involving_ruben_genes.rds")

target_genes <- sort(unique(sig_module$target))
selected_idx <- rep(F, length=nrow(sig_module))
for(g in target_genes){
    idx <- which(sig_module$target==g)
    if(length(idx) < 3){
        selected_idx[idx] <- T
    }else{
        sig_module_sub <- sig_module[idx,]
        selected_idx[idx[order(abs(sig_module_sub$pcc), decreasing=T)[1:3]]] <- T
    }
}
sig_module$selected <- selected_idx
saveRDS(sig_module, file="Data_significant_TF-target_pairs_involving_ruben_genes.rds")

# update TF-pairs after adding cytoskeleton genes
module_combined <- readRDS("Res_combined_GRN_modules_no_all_epi_cell_type_module.rds")
sig_module <- unique(module_combined[which(module_combined$padj<0.05 & module_combined$target%in%genes), c("tf", "target")])
saveRDS(sig_module, file="Data_significant_TF-target_pairs_involving_ruben_genes_v2_no_cor_cutoff.rds")
cor_vec <- sapply(seq(nrow(sig_module)), function(i){
    tf <- sig_module$tf[i]
    target <- sig_module$target[i]
    cor(ave_expr[tf,], ave_expr[target,])
})
sig_module$pcc <- cor_vec

input_module <- sig_module
input_gene <- union(input_module$target, input_module$tf)
# construct co-expression network of the genes
ave_expr <- readRDS("Data_ave_expr_per_cluster_group_per_condition.rds")
pcc_mat <- cor(t(ave_expr[input_gene,]), method="pearson")
input <- pcc_mat
pca_mat <- irlba::prcomp_irlba(input, n=20)$x
rownames(pca_mat) <- rownames(input)
umap_layout <- uwot::umap(as.matrix(pca_mat))
saveRDS(umap_layout, file="Res_salmonella_infection_GRN_umap_layout.rds")

# load infection vs control DEG test result
combined_de_res <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/Res_Salmonella_infection_on_colon_epi_goblet_and_s2E_all_DEG_res.rds")
# get the fold change of infection - control
fc_data <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/NFKB_barrier_innate_sensing/combined_with_Ruben_gene_set/Dat_expr_logFC_infection_control.rds")
pval_data <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/NFKB_barrier_innate_sensing/combined_with_Ruben_gene_set/Dat_expr_pval_infection_control.rds")
# for each expressed genes, get the minimum p-value and the maximum absolute fold change
min_pval <- setNames(apply(pval_data, 1, min),rownames(pval_data))
max_abs_fc <- setNames(apply(fc_data, 1, function(vec){vec[which.max(abs(vec))]}), rownames(fc_data))

umap_layout <- readRDS("Res_salmonella_infection_GRN_umap_layout.rds")
umap_genes <- intersect(rownames(umap_layout), names(min_pval))
node_df <- data.frame(
    "name"=umap_genes,
    "X"=umap_layout[umap_genes,1],
    "Y"=umap_layout[umap_genes,2],
    "min_pval"=min_pval[umap_genes],
    "logp"=-log10(min_pval[umap_genes]),
    "max_abs_fc"=max_abs_fc[umap_genes],
    "direction"=ifelse(max_abs_fc[umap_genes]>0, 1, -1),
    "tf"=ifelse(umap_genes %in% input_module$tf, 1, 0),
    "target"=ifelse(umap_genes %in% input_module$target, 1, 0),
    stringsAsFactors = F
)
node_df$logp[is.infinite(node_df$logp)] <- ceiling(node_df$logp[!is.infinite(node_df$logp)])
node_df$logp_transformed <- node_df$logp^0.25
# highlight the TFs with the largest logFC
tf_df <- node_df[node_df$tf==1,]
top <- tf_df[order(tf_df$max_abs_fc, decreasing=T)[1:10], "name"]
bottom <- tf_df[order(tf_df$max_abs_fc, decreasing=F)[1:10], "name"]
node_df$highlight <- ifelse(node_df$name %in% c(top, bottom), 1, 0)
saveRDS(node_df, file="Res_salmonella_infection_GRN_coex_network_node_df.rds")

# get edge list
edge_df <- input_module[which(input_module$tf %in% node_df$name & input_module$target %in% node_df$name), ]
edge_df$direction <- ifelse(edge_df[,3]>0, 1, -1)
saveRDS(edge_df, file="Res_salmonella_infection_GRN_coex_network_edge_df.rds")

# generate a graph to represent the similarity of functional enrichment between modules
# Create a graph object using tidygraph
library(tidygraph)
library(ggraph)
# Create a graph object from the edge list
ggraph <- tbl_graph(nodes = node_df, edges = edge_df, directed = FALSE)
saveRDS(ggraph, file="Res_salmonella_infection_GRN_module_coex_graph.rds")


p1 <- ggraph(ggraph, x=X, y=Y) +
    # color the edges based on the direction of the correlation
    geom_edge_diagonal(aes(color=as.factor(direction)),
                     width=0.5,
                     alpha=0.1)+
    scale_edge_color_manual(values = c("1"="#67001F", "-1"="#053061")) +
    geom_node_point(aes(size = logp, fill = max_abs_fc),
                  shape = 21, color='darkgrey') +
    scale_fill_gradient2(low="#053061", mid="#F7F7F7", high="#67001F", midpoint=0)+
    scale_size_continuous(range = c(1,10)) +
    geom_node_label(aes(label=name, filter = node_df$highlight == 1, color=as.factor(direction)), repel=T, max.overlaps = 999) +
    scale_color_manual(values = c("1"="#67001F", "-1"="#053061")) +
    #labs(title=x) +
    theme_void()
p1  
ggsave(p1, filename="Plot_ggraph_regulome_coexpression_network.pdf", height=5, width=5)

write.table(edge_df, file="Dat_inferred_TF_target_pairs.txt", sep="\t", row.names=F, quote=F)


setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/incoorperate_multiple_GRN/C16")
module_combined <- readRDS("Res_combined_GRN_modules_no_all_epi_cell_type_module.rds")

gene_x_module <- module_combined[which(module_combined$tf=="SOCS3"),]

sig_module <- unique(module_combined[which(module_combined$padj<0.05 & module_combined$target%in%genes), c("tf", "target")])
# get target genes of STAT3 and BHLHE40
selected_sig_module <- unique(module_combined[which(module_combined$padj<0.05 & module_combined$tf %in% c("STAT3", "BHLHE40")), c("tf", "target", "estimate", "pval")])
#stat3_targets <- module_combined$target[which(module_combined$padj<0.05 & module_combined$tf=="STAT3")]
# get the fold change of infection - control
fc_data <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/NFKB_barrier_innate_sensing/combined_with_Ruben_gene_set/Dat_expr_logFC_infection_control.rds")
pval_data <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/NFKB_barrier_innate_sensing/combined_with_Ruben_gene_set/Dat_expr_pval_infection_control.rds")
# focus on C16 (a Salmonella specific cluster)
fc_vec <- setNames(fc_data[,"Colonocyte_2_infection_specific_C16_C16_Salmonella"], rownames(fc_data))
pval_vec <- setNames(pval_data[,"Colonocyte_2_infection_specific_C16_C16_Salmonella"], rownames(pval_data))
selected_sig_module$target_fc <- fc_vec[selected_sig_module$target]
selected_sig_module$target_pval <- pval_vec[selected_sig_module$target]
selected_sig_module$target_pval_logp <- -log10(selected_sig_module$target_pval)
selected_sig_module$target_infection_vs_control_fc <- selected_sig_module$target_fc
selected_sig_module$target_control_vs_infection_fc <- -selected_sig_module$target_fc
selected_sig_module$target_fc_abs <- abs(fc_vec[selected_sig_module$target])
saveRDS(selected_sig_module, file="Data_significant_TF-target_pairs_by_STAT3_BHLHE40_C16.rds")
library(dplyr)
top_module <- selected_sig_module %>% group_by (tf) %>% top_n(30, target_fc_abs)
write.table(top_module, file="Table_significant_TF-target_pairs_by_STAT3_BHLHE40_top30_ranked_by_fc_in_C16.txt", sep="\t", row.names=F, quote=F)

# get related regions
module_combined <- readRDS("Res_combined_GRN_modules_no_all_epi_cell_type_module.rds")

# generate coverage plot
setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/incoorperate_multiple_GRN/coverage_plot")
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

colors <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/cols.rds")
group <- sort(unique(combined_rna$RNA_subtype_group_per_condition))
mat <- do.call('rbind', strsplit(group, split=":"))
ct <- mat[,1]
names(ct) <- group
group_cols <- setNames(colors$cell_type[ct], group)
setwd("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/NFKB_barrier_innate_sensing/start_from_featrue_gene_set/coverage_plot")
source("/home/yuq22/ihb-intestine-evo/common_script/plotting/Script_plotting_functions.R")
p_list <- coverage_plot_with_evo_sig(seu_obj=combined_rna, 
                                     peak_assay="peaks_merged", 
                                     rna_assay="RNA", 
                                     group="RNA_subtype_group_per_condition", 
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
                                     plot_suffix="-3")


