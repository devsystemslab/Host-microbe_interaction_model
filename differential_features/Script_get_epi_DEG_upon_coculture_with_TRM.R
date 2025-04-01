setwd("/home/yuq22/group_share/Bacteria_TRM/TRM/epi_DEG")
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(ggrepel)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")

# load Epi from coculture data  --- need to determine whether how high the resolution should be
epi <- readRDS("/home/yuq22/Bacteria_TRM/TRM/Dat_merged_filtered_ileum_epithelium_TRM_coculture_noImmune_noLowQuality_clustered.rds")
cols <- readRDS("/home/yuq22/group_share/Bacteria_TRM/used_object/cols.rds")
ct1 <- c("Enterocyte", "TA", "SC",  "EEC", "GC", "BEST4+")
ct2 <- setNames(c("Prox_SI_Ent", "TA", "Prox_SI_SC", "EEC", "Goblet", "BEST4+_cell"), ct1)
epi$Coarse_cell_type <- ct2[epi$anno_v1_coarse]
epi$Cluster_per_condition <- paste(epi$anno_v1, epi$orig.ident, sep=":")
saveRDS(epi, file="/home/yuq22/Bacteria_TRM/TRM/Dat_merged_filtered_ileum_epithelium_TRM_coculture_noImmune_noLowQuality_clustered_updated_v2.rds")
epi <- readRDS("/home/yuq22/Bacteria_TRM/TRM/Dat_merged_filtered_ileum_epithelium_TRM_coculture_noImmune_noLowQuality_clustered_updated_v2.rds")
p1 <- SCpubr::do_DimPlot(epi, group.by="Coarse_cell_type", colors.use=cols$cell_type, pt.size=1.5)+NoLegend()
p2 <- SCpubr::do_DimPlot(epi, group.by="orig.ident", colors.use=c("withTRM"="#303030", "noTRM"="#cccccc"),  pt.size=1.5)+NoLegend()
p <- p1+p2
p
ggsave(p, filename="Plot_UMAP_duo_epi_coculture_with_TRM.png", width=12, height=6)

# high resolution DEGs on Epi to identify co-culture induced DEGs
# perform DE between each withTRM cluster and its matched noTRM cells
# clean up the genes without gene symbol
detected_genes <- rownames(epi)[which(rowSums(epi@assays$RNA@data)>0)]
gene <- detected_genes[!grepl("^ENSG", detected_genes)]
#remove ribosomal genes, mitochondrial genes, and genes located on sex chromosomes before normalization
confound_genes <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/confound_genes/Dat_confound_genes.rds")
genes_to_exclude <- unique(unlist(confound_genes[c(4,5,6)]))
genes <- setdiff(gene, genes_to_exclude)
saveRDS(genes, file="List_input_genes_for_epi_DEG.rds")

# DEG on cluster level
selected_cl <- sort(unique(epi$anno_v1))
de_res_list <- list()
for(i in selected_cl){
    print(i)
    withTRM_cells <- colnames(epi)[which(epi$anno_v1==i & epi$orig.ident=="withTRM")]
    matched_noTRM_cells <- colnames(epi)[which(epi$anno_v1==i & epi$orig.ident=="noTRM")]
    cells <- c(withTRM_cells, matched_noTRM_cells)
    X <- epi@assays$RNA@data[genes, cells]
    y <- epi@meta.data[cells, "orig.ident"]
    de_res <- presto::wilcoxauc(X=X, y=y)
    de_res$group <- paste(i, de_res$group, sep=":")
    de_res$pct_diff <- de_res$pct_in - de_res$pct_out
    de_res$pval_log <- -log10(de_res$pval)
    de_res$pval_log[is.infinite(de_res$pval_log)] <- max(de_res$pval_log[!is.infinite(de_res$pval_log)])+1
    de_res$pval_transformed <- de_res$pval_log^0.25
    de_res$auc_sym <- abs(de_res$auc-0.5)
    de_res_list[[i]] <- de_res
    
}
cocluster_de_res <- do.call('rbind', de_res_list)
cocluster_sig_res <- cocluster_de_res %>% filter(padj<0.05 & pct_in>10 & logFC>0.1)
saveRDS(cocluster_de_res, file="Res_withTRM_vs_noTRM_cells_DEG_all_res.rds")
saveRDS(cocluster_sig_res, file="Res_withTRM_vs_noTRM_cells_DEG_sig_res.rds")

cocluster_sig_res <- readRDS("Res_withTRM_vs_noTRM_cells_DEG_sig_res.rds")
deg_num <- sort(table(cocluster_sig_res$group))
length(unique(cocluster_sig_res$feature)) # 1014
# summarize the DEG number into coarser cell types
sort(unique(sapply(unique(epi$anno_v1), function(x){strsplit(x, "-")[[1]][1]}))) -> cell_type
deg_per_cell_type <- list()
for(ct in cell_type){
    subtypes <- sort(unique(grep(paste(paste0(ct,"-"),paste0(ct,":"),sep="|"), cocluster_sig_res$group, value=T)))
    coculture_up_groups <- grep("withTRM", subtypes, value=T)
    control_up_groups <- grep("noTRM", subtypes, value=T)
    deg_per_cell_type[[paste0(ct,":coculture_up")]] <- sort(unique(cocluster_sig_res$feature[which(cocluster_sig_res$group%in%coculture_up_groups)]))
    deg_per_cell_type[[paste0(ct,":control_up")]] <- sort(unique(cocluster_sig_res$feature[which(cocluster_sig_res$group%in%control_up_groups)]))
}
saveRDS(deg_per_cell_type, file="Dat_epi_plus_minus_TRM_condition_DEG_per_cell_type.rds")
deg_num <- sapply(deg_per_cell_type, length)
n1 <- deg_num[grep("E:|SC:|TA:|GC", names(deg_num))]
df <- data.frame("cell_type"=rep(c("Prox_SI_Ent", "Goblet", "Prox_SI_SC", "TA"), each=2), "up_condition"=rep(c("withTRM","noTRM"), 4), "num_DEG"=n1)
df$num_DEG[df$up_condition=="noTRM"] <- -df$num_DEG[df$up_condition=="noTRM"]
df$abs_num_DEG <- abs(df$num_DEG)
df$cell_type <- factor(df$cell_type, levels=c("Prox_SI_Ent", "Prox_SI_SC", "TA", "Goblet"))
cols <- readRDS("/home/yuq22/group_share/Bacteria_TRM/used_object/cols.rds")
p1 <- ggplot(df, aes(x=cell_type, y=num_DEG, fill=cell_type))+
    geom_bar(stat="identity", position="dodge", col="#303030")+
    scale_fill_manual(values=cols$cell_type)+
    geom_text(aes(label=abs_num_DEG), position=position_dodge(width=0.9), vjust=-0.25)+
    theme_bw()
p1
ggsave(p1, filename="Plot_DEG_number_per_cell_type_with_GC.pdf", width=5, height=5)

# generate heatmap to show the DEG expression across cell types
ave_expr <- getAveExpr(seu.obj=epi, feature.to.calc="Cell_type_per_condition", colname.prefix=NULL)
saveRDS(ave_expr, file="Dat_epi_plus_minus_TRM_cell_type_per_condition_ave_expr.rds")
# exclude EECs and BEST4+ cells
ave_expr_sub <- ave_expr[,!grepl("BEST4|EEC", colnames(ave_expr))]
saveRDS(ave_expr_sub, file="Dat_epi_plus_minus_TRM_cell_type_per_condition_ave_expr_no_BEST4_EEC.rds")
ave_expr <- readRDS("Dat_epi_plus_minus_TRM_cell_type_per_condition_ave_expr_no_BEST4_EEC.rds")
# load DEG results
deg_per_cell_type <- readRDS("Dat_epi_plus_minus_TRM_condition_DEG_per_cell_type.rds")
union_deg <- unique(unlist(deg_per_cell_type[names(deg_per_cell_type)[!grepl("BEST4|EEC", names(deg_per_cell_type))]]))

expr_mat <- ave_expr[union_deg,]
res <- prepareTreeAndHeatmapInput(expr.mat = expr_mat, hc.method = "ward.D2",  norm.method = "quantile")
saveRDS(res, file="Res_adult_chip_epithelium_plus_minus_TRM_DEG_expr_heatmap_input.rds")
pdf("Plot_heatmap_adult_chip_epithelium_plus_minus_TRM_DEG_expr.pdf", height=10)
gplots::heatmap.2(res$heatmap_input, col = beach.col.heatmap, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 1, cexCol=0.5, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

pdf("Plot_hc_adult_chip_epithelium_plus_minus_TRM_DEG_expr.pdf")
plot(res$hc_col, main="Hierarchical clustering of adult tissue subcluster marker expression", hang=-1, cex=0.5)
dev.off()

# get top 10 genes per cluster and summarize the result on the coarse cell type level
cocluster_sig_res <- readRDS("Res_withTRM_vs_noTRM_cells_DEG_sig_res.rds")
cocluster_top_res <- cocluster_sig_res %>% group_by(group) %>% top_n(20, wt=logFC) 
saveRDS(cocluster_top_res, file="Res_withTRM_vs_noTRM_cells_DEG_top20_per_cluster.rds")
# summarize the DEG number into coarser cell types
sort(unique(sapply(unique(epi$anno_v1), function(x){strsplit(x, "-")[[1]][1]}))) -> cell_type
top_deg_per_cell_type <- list()
for(ct in cell_type){
    subtypes <- sort(unique(grep(paste(paste0(ct,"-"),paste0(ct,":"),sep="|"), cocluster_sig_res$group, value=T)))
    coculture_up_groups <- grep("withTRM", subtypes, value=T)
    control_up_groups <- grep("noTRM", subtypes, value=T)
    top_deg_per_cell_type[[paste0(ct,":coculture_up")]] <- sort(unique(cocluster_top_res$feature[which(cocluster_top_res$group%in%coculture_up_groups)]))
    top_deg_per_cell_type[[paste0(ct,":control_up")]] <- sort(unique(cocluster_top_res$feature[which(cocluster_top_res$group%in%control_up_groups)]))
}
saveRDS(top_deg_per_cell_type, file="Dat_epi_plus_minus_TRM_condition_top_DEG_per_cell_type.rds")


# check DEG overlap between cell types
n1 <- sapply(cell_type, function(ct1){
    sapply(cell_type, function(ct2){
        length(intersect(deg_per_cell_type[[ct1]], deg_per_cell_type[[ct2]]))
    })
})
# get number of DEGs that are specific to each cell type
n2 <- sapply(cell_type, function(ct){
    length(setdiff(deg_per_cell_type[[ct]], unlist(deg_per_cell_type[-which(cell_type==ct)])))
})
df <- data.frame("cell_type"=cell_type, "num_DEG"=sapply(deg_per_cell_type, length), "num_DEG_specific"=n2)
write.table(df, file="Res_withTRM_vs_noTRM_cells_DEG_summary.txt", sep="\t", quote=F, row.names=F)
df <- df[which(df$num_DEG>100),]
saveRDS(df, file="Res_withTRM_vs_noTRM_cells_DEG_number_subset.rds")
df <- readRDS("Res_withTRM_vs_noTRM_cells_DEG_number_subset.rds")
df$cell_type <- c("Prox_SI_Ent", "Prox_SI_SC", "TA")
# generate barplot to show the total numbers of DEGs and cell type specific DEGs
cols <- readRDS("/home/yuq22/group_share/Bacteria_TRM/used_object/cols.rds")
pdf("Plot_DEG_number_barplot.pdf", width=5, height=5)
barplot(df$num_DEG, names.arg=c("Enterocyte", "SC", "TA"), col=cols$cell_type[df$cell_type],  ylim=c(0, 1000), ylab="Number of DEGs")
barplot(df$num_DEG_specific, names="", add=T, density=20, angle=45, col="#202020", ylim=c(0, 1000), xaxt="n", yaxt="n")
text(df$num_DEG, labels=df$num_DEG, pos=3, cex=0.8, adj=0.5)
text(df$num_DEG_specific, labels=df$num_DEG_specific, pos=3, cex=0.8, adj=0.5)
dev.off()

# KEGG pathway enrichment on the identified DEGs
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_statiscal_test.R")
genes <- readRDS("List_input_genes_for_epi_DEG.rds")
background_genes <- genes
cocluster_sig_res <- readRDS("Res_withTRM_vs_noTRM_cells_DEG_sig_res.rds")

# focus on top DEGs
#epi_degs <- sort(unique(cocluster_sig_res$feature))
#write.table(epi_degs, file ="/home/yuq22/group_share/Bacteria_TRM/TRM/epi_immune_interaction/data/epi_degs.tsv")
# all cell types
input_res <- cocluster_sig_res
union_coculture_up_genes <- sort(unique(input_res$feature[grep(":withTRM", input_res$group)]))
union_control_up_genes <- sort(unique(input_res$feature[grep(":noTRM", input_res$group)]))
shared_genes <- intersect(union_coculture_up_genes, union_control_up_genes) # only 7 genes that show differential change direction depending on the cell types, so it is fine to include them in both list
# exclude the shared genes
gene_list <- list("withTRM"=setdiff(union_coculture_up_genes, shared_genes), "noTRM"=setdiff(union_control_up_genes, shared_genes))
saveRDS(gene_list, file="Res_epi_coculture_DEG_gene_list_exclusive.rds")

# focus on enterocytes
ent_sig_res <- cocluster_sig_res[grep("E-", cocluster_sig_res$group),]
# for all DEGs in enterocytes
input_res <- ent_sig_res
union_coculture_up_genes <- sort(unique(input_res$feature[grep(":withTRM", input_res$group)]))
union_control_up_genes <- sort(unique(input_res$feature[grep(":noTRM", input_res$group)]))
shared_genes <- intersect(union_coculture_up_genes, union_control_up_genes) # only 7 genes that show differential change direction depending on the cell types, so it is fine to include them in both list
# exclude the shared genes
gene_list <- list("withTRM"=setdiff(union_coculture_up_genes, shared_genes), "noTRM"=setdiff(union_control_up_genes, shared_genes))
saveRDS(gene_list, file="Res_epi_coculture_enterocyte_DEG_gene_list_exclusive.rds") # genes with different change directions in different clusters are included when counting DEG number, but excluded when performing functional enrichment test
# prepare input for KEGG pathway enrichment
setwd("/home/yuq22/Bacteria_TRM/TRM/epi_DEG/DAVID")
output_list <- list()
output_list[["withTRM"]] <- gene_list[["withTRM"]]
output_list[["noTRM"]] <- gene_list[["noTRM"]]
output_list[["background"]] <- background_genes

ID_gene_pairs <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/Ensembl/Human/Table_v112_ensembl_ID_and_gene_symbol.rds")
for(x in names(output_list)){
    g1 <- output_list[[x]]
    # translate the gene symbol to ensembl ID
    id1 <- unique(ID_gene_pairs$ensembl_gene_id[which(ID_gene_pairs$external_gene_name %in% g1)])
    writeLines(id1, con=paste0("List_epi_coculture_DEG_",x,".txt"))
}

# for top DEGs in enterocytes
genes <- sort(unique(ent_sig_res$feature))
max_fc_vec <- sapply(genes, function(g){
    max(ent_sig_res$logFC[which(ent_sig_res$feature==g)])
})
min_fc_vec <- sapply(genes, function(g){
    min(ent_sig_res$logFC[which(ent_sig_res$feature==g)])
})
df <- data.frame("gene"=genes, "max_fc"=max_fc_vec, "min_fc"=min_fc_vec)
saveRDS(df, file="Dat_enterocyte_DEG_max_min_fc.rds")

kegg_gene_list <- readRDS("/home/yuq22/Bacteria_TRM/Annotation/KEGG/Res_KEGG_gene_list.rds")
kegg_gene_df <- readRDS("/home/yuq22/Bacteria_TRM/Annotation/KEGG/Res_KEGG_gene_list_df.rds")
pathway_enrichment_pval <- enrichment_test(query_gene_list=gene_list, anno_gene_list=kegg_gene_list, all_genes=background_genes)
saveRDS(pathway_enrichment_pval, file="Res_immune_coculture_DEG_KEGG_hypergeometric_test_nominal_pval_enterocyte_only.rds")

pathway_enrichment_padj <- apply(pathway_enrichment_pval, 2, p.adjust, method="BH")
saveRDS(pathway_enrichment_padj, file="Res_immune_coculture_DEG_KEGG_hypergeometric_test_BH-corrected_pval_enterocyte_only.rds")

padj_data <- -log10(pathway_enrichment_padj)
# exclude the cancer and disease related terms
input <- padj_data[!grepl("disease|cancer|carcino|muscle", rownames(padj_data)),]
top_terms <- lapply(colnames(padj_data), function(x){
    intersect(rownames(input)[order(input[,x], decreasing=T)[1:10]], rownames(pathway_enrichment_padj)[pathway_enrichment_padj[,x]<0.05])
})
names(top_terms) <- colnames(padj_data)
saveRDS(top_terms, file="Res_immune_coculture_DEG_KEGG_hypergeometric_test_top_enriched_terms_enterocyte_only.rds")

terms <- c(top_terms[["withTRM"]], top_terms[["noTRM"]])
val <- c(padj_data[top_terms[["withTRM"]],"withTRM"], padj_data[top_terms[["noTRM"]],"noTRM"])
df <- data.frame('pathway'=terms, 'padj'=val, 'group'=rep(c("withTRM", "noTRM"), each=10))
library(ggpubr)
p <- ggbarplot(df, x="pathway", y="padj", orientation = "horiz", fill="#303030")+labs('title'=x, 'x'="-log10[BH-corrected P]", 'y'="Top 10 enriched KEGG pathways")+facet_wrap(~group, scales="free_y")+theme_pubr(base_size=8)
ggsave(p, filename="Plot_barplot_kegg_pathway_enrichment_gene_set_exclusive_enterocyte_only.pdf", height=5, width=10)

plot_list <- list()
for(x in c("withTRM", "noTRM")){
    df <- data.frame('pathway'=top_terms[[x]], 'padj'=padj_data[top_terms[[x]],x])
    plot_list[[x]] <- ggbarplot(df, x="pathway", y="padj", orientation = "horiz", fill="#303030")+labs('title'=x, 'x'="-log10[BH-corrected P]", 'y'="Top 10 enriched KEGG pathways")
}
p <- plot_grid(plotlist = plot_list, align = "hv",nrow = 1, axis = "l")
ggsave(p, filename="Plot_barplot_kegg_pathway_enrichment-2_gene_set_exclusive_enterocyte_only.pdf", height=5, width=10)


enriched_pathways <- lapply(seq(ncol(pathway_enrichment_padj)), function(i){
    rownames(pathway_enrichment_padj)[pathway_enrichment_padj[,i]<0.05]
})
names(enriched_pathways) <- colnames(pathway_enrichment_padj)


gs_to_check <- c("Salmonella infection", "Bacterial invasion of epithelial cells") 
gs_deg <- list()
for(gs in gs_to_check){
    for(x in names(gene_list)){
        gs_deg[[paste(gs,x,"high",sep="_")]] <- intersect(gene_list[[x]], kegg_gene_list[[gs]])
    } 
}

unique(epi$Cluster_per_condition)

# exclude "TA-1" which has lower gene numbers - low quality cells

# visualize the epithelial DEG results +/- TRM coculture
# plot vocalno plot of enterocyte DEGs upon coculture with TRM
epi <- readRDS("/home/yuq22/Bacteria_TRM/TRM/Dat_merged_filtered_ileum_epithelium_TRM_coculture_noImmune_noLowQuality_clustered.rds")

# subset to enterocyte
genes <- readRDS("List_input_genes_for_epi_DEG.rds")
i <- "Enterocyte"
withTRM_cells <- colnames(epi)[which(epi$anno_v1_coarse==i & epi$orig.ident=="withTRM")]
matched_noTRM_cells <- colnames(epi)[which(epi$anno_v1_coarse==i & epi$orig.ident=="noTRM")]
cells <- c(withTRM_cells, matched_noTRM_cells)
X <- epi@assays$RNA@data[genes, cells]
y <- epi@meta.data[cells, "orig.ident"]
de_res <- presto::wilcoxauc(X=X, y=y)
de_res$group <- paste(i, de_res$group, sep=":")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
de_res$pval_log <- -log10(de_res$pval)
de_res$pval_log[is.infinite(de_res$pval_log)] <- max(de_res$pval_log[!is.infinite(de_res$pval_log)])+1
de_res$pval_transformed <- de_res$pval_log^0.25
de_res$auc_sym <- abs(de_res$auc-0.5)
saveRDS(de_res, file="Res_enterocyte_withTRM_vs_noTRM_cells_DEG_all_res.rds")

de_res <- readRDS("Res_enterocyte_withTRM_vs_noTRM_cells_DEG_all_res.rds")
# plot vocalno plot
# subset to "Enterocyte:withTRM"
de_res <- de_res[which(de_res$group=="Enterocyte:withTRM"),]
# highlight the top 5 genes with the smallest p-value from both the upregualtion and downregulation
pos_res <- de_res[which(de_res$logFC>0),]
neg_res <- de_res[which(de_res$logFC<0),]
up_genes_to_highlight<- pos_res$feature[order(pos_res$pval_log, decreasing=T)[1:20]]
down_genes_to_highlight<- neg_res$feature[order(neg_res$pval_log, decreasing=T)[1:20]]
de_res$de_group <- rep("All", nrow(de_res))
de_res$de_group[which(de_res$feature %in% up_genes_to_highlight)] <- "Up"
de_res$de_group[which(de_res$feature %in% down_genes_to_highlight)] <- "Down" 
saveRDS(de_res, file="Res_enterocyte_withTRM_vs_noTRM_TFs_DEG_all_res_vocalno_plot_input_top20.rds")
p1 <- ggplot(de_res, aes(x=logFC, y=pval_log, color=de_group))+
    geom_point(size=1.5, aes(color=de_group))+
    geom_text_repel(data=de_res[which(de_res$de_group %in% c("Up", "Down")),], aes(label=feature, color=de_group), size=2, box.padding=0.5)+
    scale_color_manual(values=c("Up"="#B2182B", "Down"="#2166AC", "All"="#969696"))+
    theme_bw()+
    labs(x="Expr. logFC", y="-log10(P-value)")+
    ggtitle("Enterocyte DEGs upon coculture with TRM")+
    theme(plot.title = element_text(hjust = 0.5))
ggsave(p1, filename="Plot_volcano_enterocyte_DEGs_withTRM_vs_noTRM_all_gene_top20.pdf", width=6, height=5)


# highlight the top 5 genes with the smallest p-value from both the upregualtion and downregulation
pos_tf_res <- de_tf_res[which(de_tf_res$logFC>0),]
neg_tf_res <- de_tf_res[which(de_tf_res$logFC<0),]
up_genes_to_highlight<- pos_tf_res$feature[order(pos_tf_res$pval_log, decreasing=T)[1:5]]
down_genes_to_highlight<- neg_tf_res$feature[order(neg_tf_res$pval_log, decreasing=T)[1:5]]
#genes_to_highlight <- c(up_genes_to_highlight, down_genes_to_highlight)
#de_tf_res$highlight <- de_tf_res$feature %in% genes_to_highlight
de_tf_res$de_group <- rep("All", nrow(de_tf_res))
de_tf_res$de_group[which(de_tf_res$feature %in% up_genes_to_highlight)] <- "Up"
de_tf_res$de_group[which(de_tf_res$feature %in% down_genes_to_highlight)] <- "Down" 
saveRDS(de_tf_res, file="Res_enterocyte_withTRM_vs_noTRM_TFs_DEG_TF_res.rds")
p1 <- ggplot(de_tf_res, aes(x=logFC, y=pval_log, color=de_group))+
    geom_point(size=1.5, aes(color=de_group))+
    #scale_color_manual(values=c("Up"="#B2182B", "Down"="#2166AC", "All"="#969696"))+
    geom_text_repel(data=de_tf_res[which(de_tf_res$de_group %in% c("Up", "Down")),], aes(label=feature, color=de_group), size=2, box.padding=0.5)+
    scale_color_manual(values=c("Up"="#B2182B", "Down"="#2166AC", "All"="#969696"))+
    theme_bw()+
    labs(x="Expr. logFC", y="-log10(P-value)")+
    ggtitle("Enterocyte DE-TFs upon coculture with TRM")+
    theme(plot.title = element_text(hjust = 0.5))
ggsave(p1, filename="Plot_volcano_enterocyte_DEGs_withTRM_vs_noTRM_TFs.pdf", width=6, height=5)

# subset to TFs
library(Pando)
data(motif2tf)
tf_genes <- unique(motif2tf$tf)
de_tf_res <- de_res[which(de_res$feature%in%tf_genes),]

# highlight the top 5 genes with the smallest p-value from both the upregualtion and downregulation
pos_tf_res <- de_tf_res[which(de_tf_res$logFC>0),]
neg_tf_res <- de_tf_res[which(de_tf_res$logFC<0),]
up_genes_to_highlight<- pos_tf_res$feature[order(pos_tf_res$pval_log, decreasing=T)[1:5]]
down_genes_to_highlight<- neg_tf_res$feature[order(neg_tf_res$pval_log, decreasing=T)[1:5]]
#genes_to_highlight <- c(up_genes_to_highlight, down_genes_to_highlight)
#de_tf_res$highlight <- de_tf_res$feature %in% genes_to_highlight
de_tf_res$de_group <- rep("All", nrow(de_tf_res))
de_tf_res$de_group[which(de_tf_res$feature %in% up_genes_to_highlight)] <- "Up"
de_tf_res$de_group[which(de_tf_res$feature %in% down_genes_to_highlight)] <- "Down" 
saveRDS(de_tf_res, file="Res_enterocyte_withTRM_vs_noTRM_TFs_DEG_TF_res.rds")
p1 <- ggplot(de_tf_res, aes(x=logFC, y=pval_log, color=de_group))+
    geom_point(size=1.5, aes(color=de_group))+
    #scale_color_manual(values=c("Up"="#B2182B", "Down"="#2166AC", "All"="#969696"))+
    geom_text_repel(data=de_tf_res[which(de_tf_res$de_group %in% c("Up", "Down")),], aes(label=feature, color=de_group), size=2, box.padding=0.5)+
    scale_color_manual(values=c("Up"="#B2182B", "Down"="#2166AC", "All"="#969696"))+
    theme_bw()+
    labs(x="Expr. logFC", y="-log10(P-value)")+
    ggtitle("Enterocyte DE-TFs upon coculture with TRM")+
    theme(plot.title = element_text(hjust = 0.5))
ggsave(p1, filename="Plot_volcano_enterocyte_DEGs_withTRM_vs_noTRM_TFs.pdf", width=6, height=5)


SCpubr::do_DotPlot(epi, feature="ANXA2", group.by="Cell_type_per_condition")

# load and present DAVID results
up_res <- read.table("/home/yuq22/Bacteria_TRM/TRM/epi_DEG/DAVID/DAVID_enterocyte_up.txt", header=T, sep="\t")
down_res <- read.table("/home/yuq22/Bacteria_TRM/TRM/epi_DEG/DAVID/DAVID_enterocyte_down.txt", header=T, sep="\t")
up_res$Term[order(up_res$Benjamini)[1:25]]

selected_terms <- c(
"GO:0070062~extracellular exosome",                              
"hsa04612:Antigen processing and presentation",                  
"GO:0042605~peptide antigen binding",                            
"DOMAIN:Ig-like C1-type",                                        
"GO:0042613~MHC class II protein complex",                       
"GO:0042612~MHC class I protein complex",                                                                          
"KW-0443~Lipid metabolism",                                                                          
"KW-0256~Endoplasmic reticulum" ,                                
"GO:0001916~positive regulation of T cell mediated cytotoxicity",
"GO:0016324~apical plasma membrane"                             
)

pdf("Plot_barplot_DAVID_enterocyte_up.pdf", width=5, height=5)
barplot(-log10(up_res$Benjamini[match(selected_terms, up_res$Term)]), names=selected_terms, las=2)
dev.off()

up_res[match(selected_terms, up_res$Term), c("Benjamini", "Term")]
