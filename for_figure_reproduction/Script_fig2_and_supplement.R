# perform GO enrichemnt on the non-proliferative goblet cell cluster markers
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_statiscal_test.R")
expressed_genes <- rownames(chip_gc)[rowSums(chip_gc@assays$RNA@data)>0]
GO_anno <- read_GO_anno(expressed_genes=expressed_genes)
# not separating clusters
genes <- unique(sig_res$feature)
res <- GO_enrichment_test(deg=genes, GO_anno=GO_anno)
saveRDS(res, file="Res_GO_enrichment_adult_chip_non_proli_GC_cluster_markers_merged.rds")

# separating into the clusters
GO_res <- list()
for(i in sort(unique(sig_res$group))){
    genes <- sig_res$feature[which(sig_res$group==i)]
    res <- GO_enrichment_test(deg=genes, GO_anno=GO_anno)
    GO_res[[paste0("C",i)]] <- res
}
saveRDS(GO_res, file="Res_GO_enrichment_adult_chip_non_proli_GC_cluster_markers.rds")
# combine adjusted p values of all clusters
padj_mat <- sapply(seq(length(GO_res)), function(i){
    GO_res[[i]]$"BH_corrected_P"
})
rownames(padj_mat) <- rownames(GO_res[[1]])
colnames(padj_mat) <- names(GO_res)
saveRDS(padj_mat, file="Res_GO_enrichment_adult_chip_non_proli_GC_cluster_markers_padj_mat.rds")
sig_padj_mat <- padj_mat[which(rowSums(padj_mat<0.05)>0),]
sig_terms_merged <- rownames(sig_padj_mat)
gc_cl_markers_merged <- unique(sig_res$feature)
# get GO involved cluster markers and plot the heatmap to show the expression pattern across clusters

sapply(sig_terms_merged, function(term){
    genes <- GO_anno$HGNC.symbol[which(GO_anno$GO.term.name==term)], gc_cl_markers_merged
    plotFeature(seu.obj=chip_gc, dr="umap", genes.to.plot=genes, col.num=6, plot.name=paste0("Plot_UMAP_adult_chip_non_proli_GC_cluster_markers_",term,".png"), nCols=beach.col, per.plot.size=2000, cex=5)
})

sig_term_per_cluster <- lapply(seq(ncol(sig_padj_mat)), function(j){
    rownames(sig_padj_mat)[which(sig_padj_mat[,j]<0.05)]
})
names(sig_term_per_cluster) <- colnames(sig_padj_mat)

sig_res$group <- paste0("C", sig_res$group)
score_mat <- c()
id_vec <- c()
for(x in names(sig_term_per_cluster)){
    for(term in sig_term_per_cluster[[x]]){
        genes <- intersect(GO_anno$HGNC.symbol[which(GO_anno$GO.term.name==term)], sig_res$feature[which(sig_res$group==x)])
        id_vec <- c(id_vec, paste0(x, "_", term, " (", length(genes), " genes)"))
        score <- colSums(t(scale(t(cluster_expr[genes,]))))
        score_mat <- rbind(score_mat, score)

    }
}
rownames(score_mat) <- id_vec
saveRDS(score_mat, file="Dat_adult_chip_non_proli_GC_cluster_markers_enriched_GO_term_sum_scaled_expr_mat.rds")
pdf("Plot_heatmap_non_proli_goblet_cell_cluster_marker_enriched_GO_term_score.pdf", height=20, width=15)
gplots::heatmap.2(score_mat, col = beach.col.heatmap, trace="none", scale="row", density.info = "none", margins = c(15,10), cexRow = 0.2, Colv = FALSE, dendrogram = "row")
dev.off()
score_diff_vec <- sapply(rownames(score_mat), function(x){
    cl_id <- strsplit(x, "_")[[1]][1]
    cl_id <- sub("C", "Cluster_", cl_id)
    score_mat[x, cl_id]-mean(score_mat[x, setdiff(colnames(score_mat), cl_id)])
})
p_diff_vec <- c()
log_sig_padj_mat <- -log10(sig_padj_mat)
for(x in names(sig_term_per_cluster)){
    for(term in sig_term_per_cluster[[x]]){
        p_diff <- log_sig_padj_mat[term, x]-mean(log_sig_padj_mat[term, setdiff(colnames(log_sig_padj_mat), x)])
        p_diff_vec <- c(p_diff_vec, p_diff)

    }
}
df <- data.frame("feature"=id_vec, "p_diff"=p_diff_vec, "expr_score_diff"=score_diff_vec ,stringsAsFactors=F)

count_vec <- sapply(df$feature, function(x){
    id <- as.numeric(sub(" genes)", "", strsplit(x, split="(", fixed=T)[[1]][2]))
})
df$gene_count <- count_vec
df$cluster <- sapply(df$feature, function(x){
    strsplit(x, "_")[[1]][1]
})
df$GO_term <- sapply(df$feature, function(x){
    strsplit(strsplit(x, "_")[[1]][2], " (", fixed=T)[[1]][1]
})
saveRDS(df, file="Res_non_proli_GC_cluster_markers_enriched_GO_term_score_diff.rds")

df <- readRDS("Res_non_proli_GC_cluster_markers_enriched_GO_term_score_diff.rds")

library(ggrepel)
for(cl in sort(unique(df$cluster))){
    df_sub <- df[which(df$cluster==cl),]
    highlighted_terms <- union(df_sub$GO_term[order(df_sub$p_diff, decreasing=T)[1:10]], df_sub$GO_term[order(df_sub$expr_score_diff, decreasing=T)[1:10]])
    p1 <- ggplot(df_sub, aes(x=p_diff, y=expr_score_diff))+
        geom_point(aes(size=gene_count), pch=16, col="#303030")+
        geom_text_repel(data=df_sub[df_sub$GO_term%in%highlighted_terms,], aes(label=GO_term), box.padding=2, size=3, max.overlaps=999)+
        theme_bw()+
        xlab("BH-adjusted P diff")+
        ylab("Scaled expr. idff")+
        #scale_color_manual(values=c("grey", "dark red"))+
        theme(axis.text=element_text(size=10), axis.title=element_text(size=10), axis.ticks.length=unit(0.3, "cm"), axis.ticks = element_line(colour = "black", linewidth = 0.5))
    pdf(paste0("Plot_",cl, "_marker_enriched_GO_term_score_diff.pdf"), height=10, width=10)
    print(p1)
    dev.off()
}
