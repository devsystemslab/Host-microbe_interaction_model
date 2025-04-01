setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings")
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")

# load data
combined <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/Dat_filtered_colon_epithelium_infection_RNA_preprocessed.rds")
# filter out low quality cells
combined <- subset(
  x = combined,
  subset = nFeature_RNA < 5000 & 
    nFeature_RNA > 500 &
    nCount_RNA < 20000 & 
    percent.mt < 50
)
p2 <- SCpubr::do_ViolinPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=4, group.by="orig.ident")&NoLegend()
pdf("Plot_quality_check_filtered_data.pdf", height=5, width=15)
p2
dev.off()
# clean up the genes without gene symbol
#count <- combined@assays$RNA@counts
#gene <- rownames(count)[!grepl("^ENSG", rownames(count))]
#remove ribosomal genes, mitochondrial genes, and genes located on sex chromosomes before normalization
#confound_genes <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/confound_genes/Dat_confound_genes.rds")
#genes_to_exclude <- unique(unlist(confound_genes[c(4,5,6)]))
#genes <- setdiff(gene, genes_to_exclude)
#count <- count[genes,]
#meta <- combined@meta.data
#combined <- CreateSeuratObject(counts=count, meta.data=meta)
combined <- NormalizeData(object = combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000)
combined <- ScaleData(object = combined, verbose = T)
combined <- RunPCA(object = combined, features = VariableFeatures(combined), verbose = F, npcs = 50)
usefulPCs <- 1:20
combined <- FindNeighbors(object = combined, dims = usefulPCs)
combined <- FindClusters(object = combined, resolution = 1)
combined <- RunUMAP(object = combined, dims = usefulPCs)
SCpubr::do_DimPlot(combined, reduction="umap", group.by="orig.ident", label=T)
saveRDS(combined, file="Dat_salmonella_colon_epithelium_scRNA-seq_data.rds")

#combined <- simspec::cluster_sim_spectrum(combined, label_tag = "condition", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, return_seuratObj = TRUE)
#combined <- RunUMAP(combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
#combined <- FindNeighbors(object = combined, reduction = "css", dims = 1:ncol(combined@reductions$css@cell.embeddings), force.recalc = T) %>%
#  FindClusters(resolution = 1)
#combined[[paste0("RNA_CSS_snn_res.", 0.2*5)]] <- combined[[paste0("RNA_snn_res.", 0.2*5)]]
#saveRDS(combined, file="Dat_salmonella_colon_epithelium_scRNA-seq_data_with_CSS_integration.rds")
## the CSS integration result less resolved the cell types, so stick with the merged data

combined <- readRDS("Dat_salmonella_colon_epithelium_scRNA-seq_data.rds")
# get the condition distribution per cluster
n1 <- sapply(sort(unique(combined$RNA_snn_res.1)), function(i){
    sapply(sort(unique(combined$condition)), function(j){
        sum(combined$condition[combined$RNA_snn_res.1==i]==j)
    })
})
colnames(n1) <- paste0("C", sort(unique(combined$RNA_snn_res.1)))
prop1 <- t(t(n1)/colSums(n1))
prop_diff <- prop1[2,]-prop1[1,]

p1 <- SCpubr::do_DimPlot(combined, reduction="umap", group.by="RNA_snn_res.1", label=T)+NoLegend()
p2 <- SCpubr::do_DimPlot(combined, reduction="umap", group.by="condition")
p <- plot_grid(plotlist = list(p1, p2), align = "hv",nrow = 1, axis = "l")
ggsave(p, filename="Plot_UMAP_combined_salmonella.png", width=5*2, height=5)

de_res <- presto::wilcoxauc(combined, group_by="condition", seurat_assay="RNA")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_in>10 & logFC>0.1)
table(sig_res$group)
top_res <- sig_res %>% group_by(group) %>% top_n(100, wt=logFC)
deg_res <- list("top_res"=top_res, "sig_res"=sig_res)
saveRDS(deg_res, file="Res_Salmonella_infection_on_colon_epi_global_DEG.rds")

condition_deg_res <- readRDS("Res_Salmonella_infection_on_colon_epi_global_DEG.rds")
condition_sig_deg <- condition_deg_res$sig_res
condition_sig_deg[which(condition_sig_deg$feature%in%c("NEDD4L", "SHARPIN")),]
plotFeature(combined, genes.to.plot = c("HILPDA","ID2","PLCG2","HES1","CSKMT","SLC7A11","TMC5","GFPT1","SLC26A2","CLCA4","CA1","CA2","CES2"), col.num=4, plot.name="Plot_UMAP_combined_salmonella_feature_plot_2.png")
# coarse cell type annotation
combined$Coarse_cell_type <- "Colonocyte"
ct <- list(
    "SC"=c(6,13),
    "BEST4+_cell"=19,
    "TA"=8,
    "GC"=c(5,9),
    "EEC"=20
)
for(x in names(ct)){
    combined$Coarse_cell_type[combined$RNA_snn_res.1%in%ct[[x]]] <- x
}
SCpubr::do_DimPlot(combined, reduction="umap", group.by="Coarse_cell_type", label=T)+NoLegend()
combined$Ct_per_condition <- paste(combined$Coarse_cell_type, combined$condition, sep=":")
saveRDS(combined, file="Dat_salmonella_colon_epithelium_scRNA-seq_data.rds")

p1 <- SCpubr::do_FeaturePlot(combined, reduction="umap", features="nFeature_RNA", order=T)
ggsave(p1, filename="Plot_UMAP_RNA_colored_by_gene_number.png", width=5, height=5)


p1 <- SCpubr::do_DimPlot(combined, reduction="umap", group.by="RNA_snn_res.1", split.by="RNA_snn_res.1")&NoLegend()
p1

# identify cluster markers
de_res <- presto::wilcoxauc(combined, group_by="RNA_snn_res.1", seurat_assay="RNA")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_in>10 & logFC>0.1)
table(sig_res$group)
top_res <- sig_res %>% group_by(group) %>% top_n(10, wt=logFC)
deg_res <- list("top_res"=top_res, "sig_res"=sig_res)
saveRDS(deg_res, file="Res_Salmonella_infection_on_colon_epi_cluster_markers.rds")

for(cl in c(10,16,17,18)){
    genes <- top_res$feature[top_res$group==cl]
    plotFeature(combined, genes.to.plot = genes, col.num=5, plot.name=paste0("Plot_UMAP_combined_salmonella_feature_plot_cluster_",cl,".png"))
}

# generate a scatter plot where x and y axis represent the specificity of logFC and maximal logFC between condition, with size proportional to reference number with the gene and Salmonella cooccurrence, and color represent the cell type with the maximal logFC

# get the genes showing top correlation with SAA1 and SAA2 gene expression in colonocytes
colonocyte <- subset(combined, subset=Coarse_cell_type=="Colonocyte")
colonocyte <- FindVariableFeatures(colonocyte, selection.method = "vst", nfeatures = 2000)
expr_mat <- as.matrix(colonocyte@assays$RNA@data[VariableFeatures(colonocyte),])
sd_vec <- apply(expr_mat, 1, sd)
ref_pattern <- colSums(expr_mat[c("SAA1","SAA2"),])
cor_vec <- cor(t(expr_mat), ref_pattern)[,1]
genes <- names(cor_vec)[order(cor_vec, decreasing=T)[1:12]]
plotFeature(seu.obj=colonocyte, genes.to.plot=genes, plot.name="Plot_colonocyte_SAA1_SAA2_correlated_genes.png")
png("Plot_colonocyte_SAA1_SAA2_correlated_genes_featurePlot_by_condition.png", width=2000*6, height=2000*4)
par(mfrow=c(4,6), mar=c(10,10,10,10))
for(g in genes){
    for(x in unique(colonocyte$condition)){
        plotFeature2(coor=Embeddings(colonocyte, reduction="umap"), values=colonocyte@assays$RNA@data[g,], emphasize=colonocyte$condition==x, main=paste(g, x, sep=":"), cex.main=10, cex=8, point.order="sorted")
    }
}
dev.off()

p1 <- SCpubr::do_ViolinPlot(colonocyte, features=genes, ncol=4, group.by="condition")&NoLegend()
ggsave(p1, filename="Plot_colonocyte_SAA1_SAA2_correlated_genes_violinPlot.pdf", width=5*4, height=5*3)

combined$Ct_per_condition <- paste(combined$Coarse_cell_type, combined$condition, sep=":")
saveRDS(combined, file="Dat_salmonella_colon_epithelium_scRNA-seq_data.rds")

combined <- readRDS("Dat_salmonella_colon_epithelium_scRNA-seq_data.rds")
ave_expr <- getAveExpr(seu.obj=combined, feature.to.calc="Ct_per_condition", colname.prefix=NULL)
saveRDS(ave_expr, file="Res_cell_type_expression_per_condition.rds")
expr_diff <- ave_expr[,grep("Salmonella",colnames(ave_expr))] - ave_expr[,grep("Control",colnames(ave_expr))]
saveRDS(expr_diff, file="Res_cell_type_expression_diff_Salmonella_over_Control_per_condition.rds")


tissue_s2e <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/RNA_per_condition/tissue_colon/Res_colon_s2e_adult_human_epi_with_CSS.rds")
# subset to colonocytes
tissue_colonocyte <- subset(tissue_s2e, cells=colnames(tissue_s2e)[which(!tissue_s2e$RNA_CSS_snn_res.0.4%in%c(0,2,8))])
colonocyte <- subset(combined, subset=Coarse_cell_type=="Colonocyte")
p1 <- SCpubr::do_DimPlot(tissue_colonocyte, reduction="umap_css", group.by="RNA_CSS_snn_res.0.4", label=T)
p2 <- SCpubr::do_DimPlot(colonocyte, reduction="umap", group.by="RNA_snn_res.1", label=T)
p1/p2

de_res <- presto::wilcoxauc(tissue_colonocyte, group_by="RNA_CSS_snn_res.0.4", seurat_assay="RNA")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_in>10 & logFC>0.1)
table(sig_res$group)
top_res <- sig_res %>% group_by(group) %>% top_n(100, wt=logFC)
genes <- intersect(top_res$feature, rownames(combined))

tissue_expr <- getAveExpr(seu.obj=tissue_colonocyte, feature.to.calc="RNA_CSS_snn_res.0.4", colname.prefix="tisssue")
chip_expr <-   getAveExpr(seu.obj=colonocyte, feature.to.calc="RNA_snn_res.1", colname.prefix="chip")
cor_mat <- cor(tissue_expr[genes,], chip_expr[genes,], method="spearman")
pdf("Plot_heatmap_SCC_s2e_cluster_between_chip_and_tissue.pdf", height=10, width=10)
gplots::heatmap.2(cor_mat, scale="column", cellnote=round(cor_mat,2), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10))
dev.off()

chip_colonocyte_3 <- c(10, 2)
chip_colonocyte_2 <- c(1,3,7,18)
chip_colonocyte_1 <- setdiff(colonocyte$RNA_snn_res.1, c(chip_colonocyte_3, chip_colonocyte_2))
combined$Coarse_cell_type[which(combined$Coarse_cell_type=="SC")] <- "Colon_SC"
combined$Coarse_cell_type[which(combined$Coarse_cell_type=="GC")] <- "Goblet"
combined$Ct_per_condition <- paste(combined$Coarse_cell_type, combined$condition, sep=":")
combined$Cell_type_with_colonocyte_subtypes <- combined$Coarse_cell_type
combined$Cell_type_with_colonocyte_subtypes[combined$RNA_snn_res.1%in%chip_colonocyte_3] <- "Colonocyte_3"
combined$Cell_type_with_colonocyte_subtypes[combined$RNA_snn_res.1%in%chip_colonocyte_2] <- "Colonocyte_2"
combined$Cell_type_with_colonocyte_subtypes[combined$RNA_snn_res.1%in%chip_colonocyte_1] <- "Colonocyte_1"
SCpubr::do_DimPlot(combined, reduction="umap", group.by="Cell_type_with_colonocyte_subtypes", label=T)+NoLegend()
n1 <- sapply(sort(unique(combined$Cell_type_with_colonocyte_subtypes)), function(i){
    sapply(sort(unique(combined$condition)), function(j){
        sum(combined$condition[combined$Cell_type_with_colonocyte_subtypes==i]==j)
    })
})
prop1 <- n1/rowSums(n1)

combined$Subtype_per_condition <- paste(combined$Cell_type_with_colonocyte_subtypes, combined$condition, sep=":")
saveRDS(combined, file="Dat_salmonella_colon_epithelium_scRNA-seq_data.rds")

# generate stacked bar plot showing the distribution of cell types in each condition
# create a dataset
df <- data.frame("condition"=rep(rownames(n1), ncol(n1)), "cell_type"=rep(colnames(n1), each=nrow(n1)), "cell_number"=as.vector(n1))
saveRDS(df, file="Dat_cell_type_composition_per_condition.rds")

# Stacked + percent
cols <- readRDS("~/Bacteria_TRM/used_object/cols.rds")
cols$cell_type["BEST4+_cell"] <- "#C51B8A"
cols$cell_type <- c(cols$cell_type, setNames(c("#FFFFB2", "#FECC5C", "#FD8D3C"), paste("Colonocyte", 1:3, sep="_")))
saveRDS(cols, file="~/Bacteria_TRM/used_object/cols.rds")
p1 <- ggplot(df, aes(fill=cell_type, y=cell_number, x=condition)) + 
    geom_bar(position="fill", stat="identity")+
    scale_fill_manual(values=cols$cell_type)+
    labs(y="Proportion", x="Condition", fill="Cell type")+
    theme_minimal()

ggsave(p1, filename="Plot_cell_type_composition_per_condition_stacked_bar_plot.pdf", width=5, height=5)

# separate the SC subtypes
combined <- readRDS("Dat_salmonella_colon_epithelium_scRNA-seq_data.rds")
p1 <- SCpubr::do_DimPlot(combined, reduction="umap", group.by="Cell_type_with_colonocyte_subtypes", label=T)+NoLegend()
p2 <- SCpubr::do_DimPlot(combined, reduction="umap", group.by="RNA_snn_res.1", label=T)+NoLegend()


tissue_s2e <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/RNA_per_condition/tissue_colon/Res_colon_s2e_adult_human_epi_with_CSS.rds")
tissue_sc <- subset(tissue_s2e, cells=colnames(tissue_s2e)[which(tissue_s2e$RNA_CSS_snn_res.0.4%in%c(2,8))])
sc <- subset(combined, subset=Coarse_cell_type=="Colon_SC")
p1 <- SCpubr::do_DimPlot(tissue_sc, reduction="umap_css", group.by="RNA_CSS_snn_res.0.4", label=T)
p2 <- SCpubr::do_DimPlot(sc, reduction="umap", group.by="RNA_snn_res.1", label=T)
p1/p2

de_res <- presto::wilcoxauc(tissue_sc, group_by="RNA_CSS_snn_res.0.4", seurat_assay="RNA")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_in>10 & logFC>0.1)
table(sig_res$group)
top_res <- sig_res %>% group_by(group) %>% top_n(100, wt=logFC)
genes <- intersect(top_res$feature, rownames(combined))

tissue_expr <- getAveExpr(seu.obj=tissue_sc, feature.to.calc="RNA_CSS_snn_res.0.4", colname.prefix="tisssue")
chip_expr <-   getAveExpr(seu.obj=sc, feature.to.calc="RNA_snn_res.1", colname.prefix="chip")
cor_mat <- cor(tissue_expr[genes,], chip_expr[genes,], method="spearman")

# while chip C13 is much more similar to the SMOC+/LGR5+ stem cell cluster in tissue, C6 shows comparable similarity to the SMOC+/LGR5+ stem cell cluster and the SMOC-/LGR5-/OLFM4+ stem cell cluster in tissue
combined$Subtype <- combined$Cell_type_with_colonocyte_subtypes
combined$Subtype[which(combined$RNA_snn_res.1%in%13)] <- "Colon_SC_1"
combined$Subtype[which(combined$RNA_snn_res.1%in%6)] <- "Colon_SC_2"
combined$Colonocyte_subtype_per_condition <- combined$Subtype_per_condition
combined$Subtype_per_condition <- paste(combined$Subtype, combined$condition, sep=":")
table(combined$Subtype_per_condition)
saveRDS(combined, file="Dat_salmonella_colon_epithelium_scRNA-seq_data_with_subtypes.rds")

n1 <- sapply(sort(unique(combined$Subtype)), function(i){
    sapply(sort(unique(combined$condition)), function(j){
        sum(combined$condition[combined$Subtype==i]==j)
    })
})
# generate stacked bar plot showing the distribution of cell types in each condition
# create a dataset
df <- data.frame("condition"=rep(rownames(n1), ncol(n1)), "cell_type"=rep(colnames(n1), each=nrow(n1)), "cell_number"=as.vector(n1))
saveRDS(df, file="Dat_cell_type_composition_per_condition_2.rds")

# Stacked + percent
cols <- readRDS("~/Bacteria_TRM/used_object/cols.rds")
cols$cell_type <- c(cols$cell_type, setNames(c(cols$cell_type["Colon_SC"], "#9b59b6"), paste("Colon_SC", 1:2, sep="_")))
saveRDS(cols, file="~/Bacteria_TRM/used_object/cols.rds")
p1 <- ggplot(df, aes(fill=cell_type, y=cell_number, x=condition)) + 
    geom_bar(position="fill", stat="identity")+
    scale_fill_manual(values=cols$cell_type)+
    labs(y="Proportion", x="Condition", fill="Cell type")+
    theme_minimal()

ggsave(p1, filename="Plot_cell_type_composition_per_condition_stacked_bar_plot_2.pdf", width=5, height=5)


tissue_s2e <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/RNA_per_condition/tissue_colon/Res_colon_s2e_adult_human_epi_with_CSS.rds")
tissue_s2e <- subset(tissue_s2e, subset=RNA_CSS_snn_res.0.4!=0)
chip_s2e <- subset(combined, subset=Coarse_cell_type%in%c("Colon_SC", "Colonocyte"))
p1 <- SCpubr::do_DimPlot(tissue_s2e, reduction="umap_css", group.by="RNA_CSS_snn_res.0.4", label=T)
p2 <- SCpubr::do_DimPlot(chip_s2e, reduction="umap", group.by="RNA_snn_res.1", label=T)
p1/p2

de_res <- presto::wilcoxauc(tissue_s2e, group_by="RNA_CSS_snn_res.0.4", seurat_assay="RNA")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_in>10 & logFC>0.1)
table(sig_res$group)
top_res <- sig_res %>% group_by(group) %>% top_n(100, wt=logFC)
genes <- intersect(top_res$feature, rownames(combined))

tissue_expr <- getAveExpr(seu.obj=tissue_s2e, feature.to.calc="RNA_CSS_snn_res.0.4", colname.prefix="tisssue")
chip_expr <-   getAveExpr(seu.obj=chip_s2e, feature.to.calc="RNA_snn_res.1", colname.prefix="chip")
cor_mat <- cor(tissue_expr[genes,], chip_expr[genes,], method="spearman")
res <- gplots::heatmap.2(cor_mat, scale="column", cellnote=round(cor_mat,2), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10))
saveRDS(res, file="Res_heatmap_input_SCC_s2e_cluster_between_chip_and_tissue.rds")
pdf("Plot_heatmap_SCC_s2e_cluster_between_chip_and_tissue_2.pdf", height=10, width=10)
gplots::heatmap.2(cor_mat, scale="column", cellnote=round(cor_mat,2), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10))
dev.off()

hc <- as.hclust(res$colDendrogram)
plot(hc, hang=-1)
g <- cutree(hc, h=0.3)
ct <- c("Colonocyte_2","Colonocyte_3","Colonocyte_4","Colonocyte_1","Colon_SC")
for(i in unique(g)){
    g[which(g==i)] <- ct[i]
}
combined$Coarse_cell_type[which(combined$RNA_snn_res.1==6)] <- "Colonocyte"
combined$Ct_per_condition <- NULL
combined$Cell_type_with_colonocyte_subtypes <- NULL
combined$Subtype_per_condition <- NULL
combined$Subtype <- NULL
combined$Colonocyte_subtype_per_condition <- NULL
combined$Subtype <- combined$Coarse_cell_type
for(x in unique(g)){
    combined$Subtype[which(combined$RNA_snn_res.1%in%sub("chip", "", names(g)[which(g==x)]))] <- x
}
df <- unique(combined@meta.data[,c("Subtype", "RNA_snn_res.1")])
df <- df[order(df$Subtype),]
SCpubr::do_DimPlot(combined, reduction="umap", group.by="Subtype", label=T)+NoLegend()
combined$Subtype_per_condition <- paste(combined$Subtype, combined$condition, sep=":")
saveRDS(combined, file="Dat_salmonella_colon_epithelium_scRNA-seq_data_with_subtypes.rds")


combined <- readRDS("Dat_salmonella_colon_epithelium_scRNA-seq_data_with_subtypes.rds")
n1 <- sapply(sort(unique(combined$Subtype)), function(i){
    sapply(sort(unique(combined$condition)), function(j){
        sum(combined$condition[combined$Subtype==i]==j)
    })
})
# generate stacked bar plot showing the distribution of cell types in each condition
# create a dataset
df <- data.frame("condition"=rep(rownames(n1), ncol(n1)), "cell_type"=rep(colnames(n1), each=nrow(n1)), "cell_number"=as.vector(n1))
saveRDS(df, file="Dat_cell_type_composition_per_condition_2.rds")

# Stacked + percent
cols <- readRDS("~/Bacteria_TRM/used_object/cols.rds")
cols$cell_type <- c(cols$cell_type[1:10], setNames(c("#FFFFB2", "#FECC5C", "#FD8D3C", "#D35400"), paste("Colonocyte", 1:4, sep="_")))
saveRDS(cols, file="~/Bacteria_TRM/used_object/cols.rds")
p1 <- ggplot(df, aes(fill=cell_type, y=cell_number, x=condition)) + 
    geom_bar(position="fill", stat="identity")+
    scale_fill_manual(values=cols$cell_type)+
    labs(y="Proportion", x="Condition", fill="Cell type")+
    theme_minimal()
ggsave(p1, filename="Plot_cell_type_composition_per_condition_stacked_bar_plot_2.pdf", width=5, height=5)

p1 <- SCpubr::do_DimPlot(combined, reduction="umap", group.by="Subtype", label=T, colors.use=cols$cell_type)+NoLegend()
ggsave(p1, filename="Plot_UMAP_salmonella_dataset_subtypes.png", width=5, height=5)

# perform DEG analysis on the cell types
# generate volcano plot per cell type
# only perform DE test on selected genes
# clean up the genes without gene symbol
gene <- rownames(combined)[!grepl("^ENSG", rownames(combined))]
#remove ribosomal genes, mitochondrial genes, and genes located on sex chromosomes before normalization
confound_genes <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/confound_genes/Dat_confound_genes.rds")
genes_to_exclude <- unique(unlist(confound_genes[c(4,5,6)]))
genes <- setdiff(gene, genes_to_exclude)

# generate vocalno plot per cell type
sig_deg_res_per_ct <- list()
input_list <- list()
for(ct in sort(unique(combined$Subtype))){
    
    print(paste(ct, 'start'))
    idx <- which(combined$Subtype==ct)
    X <- combined@assays$RNA@data[genes,idx]
    y <- combined$Subtype_per_condition[idx]

    de_res <- presto::wilcoxauc(X=X, y=y)
    de_res$Subtype <- ct
    de_res$pct_diff <- de_res$pct_in - de_res$pct_out
    de_res$pval_log <- -log10(de_res$pval)
    de_res$pval_log[is.infinite(de_res$pval_log)] <- max(de_res$pval_log[!is.infinite(de_res$pval_log)])+1
    de_res$pval_transformed <- de_res$pval_log^0.25
    de_res$auc_sym <- abs(de_res$auc-0.5)
    sig_res <- de_res %>% filter(padj<0.05 & pct_in>5 & logFC>0.05)
    sig_deg_res_per_ct[[ct]] <- sig_res
    sig_genes <- sig_res$feature

    input <- de_res[grep("Salmonella", de_res$group),]
    
    pos_genes <- input$feature[order(input$logFC, decreasing=T)[1:10]]
    neg_genes <- input$feature[order(input$logFC, decreasing=F)[1:10]]
    input$group <- "Background"
    input$group[input$feature%in%sig_genes] <- "DEG"
    input$group[input$feature%in%pos_genes] <- "Top Salmonella high"
    input$group[input$feature%in%neg_genes] <- "Top control high"
    input_list[[ct]] <- input
 
}
sig_deg_res_per_ct_combined <- do.call('rbind', sig_deg_res_per_ct)
saveRDS(sig_deg_res_per_ct_combined, file="Res_Salmonella_infection_on_colon_epi_subtype_DEG.rds")
top_deg_res_per_ct_combined <- sig_deg_res_per_ct_combined %>% group_by(group) %>% top_n(50, wt=logFC)
saveRDS(top_deg_res_per_ct_combined, file="Res_Salmonella_infection_on_colon_epi_subtype_top_DEG.rds")
union_deg <- unique(top_deg_res_per_ct_combined$feature)

library(cowplot)
library(grid)
library(gridExtra)
library(cowplot)
library(ggrepel)

group_cols <- setNames(c("#a0a0a0", "#696969", "#922b21", "#1f618d"), c("Background", "DEG", "Top Salmonella high", "Top control high"))
input <- do.call('rbind', input_list)
saveRDS(input, file="Res_Salmonella_infection_on_colon_epi_subtype_DEG_all_res.rds")
p1 <- ggplot(input, aes(x=logFC, y=pval_transformed, size=auc_sym, fill=group))+
    geom_point(shape=21, color="#202020")+
    theme_bw()+
    theme(legend.position="bottom")+
    scale_fill_manual(values=group_cols, guide = guide_legend(override.aes = list(label = "")))+
    labs(x="logFC (Infection - Control)", size="|AUC-0.5|", y="-log10(P)^0.25")+
    #geom_vline(xintercept=0, color="black", linetype="dashed")+
    geom_text_repel(data=input[input$group=="Top Salmonella high",], aes(label=feature), color="#922b21", size=4, point.padding=1, segment.color='grey50', min.segment.length=0, box.padding=0.5, max.overlaps=999)+
    geom_text_repel(data=input[input$group=="Top control high",], aes(label=feature), color="#1f618d", size=4, point.padding=1, segment.color='grey50', min.segment.length=0, box.padding=0.5, max.overlaps=999)
    #scale_color_manual(values=cols$cell_type, guide = "none")
p1_split <- p1 + facet_wrap(~ Subtype, scales = "free")
#p <- plot_grid(plotlist = plot_list, align = "hv", nrow = 2, axis = "l")
ggsave(p1_split, filename="Plot_Salmonella_infection_on_colon_epi_subtype_DEG_volcano_plot.png", width=10*3, height=10*ceiling(length(sort(unique(combined$Subtype)))/3))


df <- input[which(input$Subtype=="Colon_SC"),]
p1 <- ggplot(df, aes(x=logFC, y=pval_transformed, size=auc_sym, fill=group))+
    geom_point(shape=21, color="#202020")+
    theme_bw()+
    theme(legend.position="bottom")+
    scale_fill_manual(values=group_cols, guide = guide_legend(override.aes = list(label = "")))+
    labs(x="logFC (Infection - Control)", size="|AUC-0.5|", y="-log10(P)^0.25")+
    geom_text_repel(data=df[df$feature=="LTB",], aes(label=feature), color="#922b21", size=4, point.padding=1, segment.color='grey50', min.segment.length=0, box.padding=0.5, max.overlaps=999)
p1

p2 <- ggplot(input, aes(x=Subtype, y=auc))+
      geom_boxplot()

SCpubr::do_FeaturePlot(combined, feature="NFKBIA", split.by="orig.ident")

selected_cell_type <- c("Colonocyte_1", "Colonocyte_2", "Colonocyte_3", "Colonocyte_4", "Colon_SC", "Goblet")
sig_deg_res_subset <- sig_deg_res_per_ct_combined[sig_deg_res_per_ct_combined$Subtype%in%selected_cell_type,]
deg_list <- tapply(sig_deg_res_subset$feature, sig_deg_res_subset$group, list)
# perform gene set enrichment analysis on the DEGs
# KEGG pathway enrichment test
library(msigdbr)
all_gene_sets <- msigdbr(species = "Homo sapiens")
h_genes <- msigdbr(species = "Homo sapiens", category = "H")
h_gene_list <- tapply(h_genes$gene_symbol, h_genes$gs_name, list)
kegg_genes <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
kegg_gene_list <- tapply(kegg_genes$gene_symbol, kegg_genes$gs_name, list)
wikipathway_genes <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS")
wikipathway_gene_list <- tapply(wikipathway_genes$gene_symbol, wikipathway_genes$gs_name, list)
merged_list <- c(h_gene_list, kegg_gene_list, wikipathway_gene_list)
saveRDS(merged_list , file="Res_KEGG_wikiPathway_hallMark_merged_annotated_gene_list.rds")

source("/projects/site/pred/ihb-intestine-evo/common_script/Script_statiscal_test.R")
pathway_enrichment_pval <- enrichment_test(query_gene_list=deg_list, anno_gene_list=merged_list, all_genes=rownames(combined))
saveRDS(pathway_enrichment_pval, file="Res_salmonella_infection_DEG_KEGG_hallmark_wikiPathway_hypergeometric_test_nominal_pval.rds")

# perform p-value correction per gene list
idx_kegg <- grep("KEGG", rownames(pathway_enrichment_pval))
idx_wp <- grep("WP", rownames(pathway_enrichment_pval))
idx_hallmark <- grep("HALLMARK", rownames(pathway_enrichment_pval))
padj_list <- lapply(c("HALLMARK", "KEGG", "WP"), function(x){
    idx <- grep(x, rownames(pathway_enrichment_pval))
    pval <- pathway_enrichment_pval[idx,]
    apply(pval, 2, p.adjust, method="BH")
})
padj_mat <- do.call('rbind', padj_list)
saveRDS(padj_mat, file="Res_salmonella_infection_DEG_KEGG_hallmark_wikiPathway_hypergeometric_test_BH_padj.rds")

combined <- readRDS("Dat_salmonella_colon_epithelium_scRNA-seq_data_with_subtypes.rds")
s2e <- subset(combined, subset=Coarse_cell_type%in%c("Colonocyte", "Colon_SC"))
SCpubr::do_DimPlot(s2e, reduction="umap", group.by="Subtype", label=T)+NoLegend()
dm <- destiny::DiffusionMap(Embeddings(s2e, "pca")[,1:20])
# plot the DC1 and DC2 of the colonocytes
dc <- dm@eigenvectors
s2e@meta.data[,colnames(dc)] <- dc
# use rank of DC1 as pseudotime
dpt <- destiny::DPT(dm)
s2e$dpt <- rank(dpt$dpt)
s2e$dc1_rank <- rank(-dc[,1])
coor <- dc[,c(1,2)]
colnames(coor) <- c("DC1", "DC2")
s2e[["dm"]] <- CreateDimReducObject(embeddings = coor, key="DC_", assay="RNA")
saveRDS(s2e, file="Res_salmonella_data_s2e_with_DM.rds")
cols <- readRDS("~/Bacteria_TRM/used_object/cols.rds")
median_rank_per_cluster <- sapply(sort(unique(s2e$RNA_snn_res.1)), function(i){
    median(s2e$dc1_rank[which(s2e$RNA_snn_res.1==i)])
})
names(median_rank_per_cluster) <- paste0("C", sort(unique(s2e$RNA_snn_res.1)))
s2e$RNA_snn_res.1 <- factor(s2e$RNA_snn_res.1, levels=sub("C", "", names(sort(median_rank_per_cluster))))
saveRDS(s2e, file="Res_salmonella_data_s2e_with_DM.rds")
p1 <- SCpubr::do_RidgePlot(s2e, feature="dc1_rank", group.by="RNA_snn_res.1")
p1
ggsave(p1, filenam="Plot_cluster_DC1_rank.pdf", width=5, height=10)
p1 <- SCpubr::do_DimPlot(s2e, reduction="dm", group.by="Subtype", pt.size=2.5, border.size=1.5, colors.use=cols$cell_type, legend.nrow=2)
p2 <- SCpubr::do_DimPlot(s2e, reduction="dm", group.by="RNA_snn_res.1", pt.size=2.5, border.size=1.5, label=T)+NoLegend()
p3 <- SCpubr::do_DimPlot(s2e, reduction="dm", group.by="condition", pt.size=2.5, border.size=1.5)
p4 <- SCpubr::do_DimPlot(s2e, reduction="umap", group.by="RNA_snn_res.1", pt.size=2.5, border.size=1.5, label=T)+NoLegend()
p5 <- SCpubr::do_FeaturePlot(s2e, reduction="umap", feature="dpt", pt.size=2.5, border.size=1.5)
p6 <- SCpubr::do_FeaturePlot(s2e, reduction="umap", feature="dc1_rank", pt.size=2.5, border.size=1.5)
p <- plot_grid(plotlist = list(p1, p3, p2, p4, p5, p6), align = "hv", nrow = 3, axis = "l")
ggsave(p, filename="Plot_DC1_2_salmonella_data_s2E.png", width=10, height=15)

# establish a framework to perform high resolution DE analysis
# and visualize the DEG expression pattern on a differentiation related single dimension
s2e <- readRDS("Res_salmonella_data_s2e_with_DM.rds")
# for each cell in the Salmonella infected sample, get its control cell counterpart in the top 20 PCs space
library(FNN)
# Get the top PC values for Salmonella and Control samples
salmonella_pca <- Embeddings(s2e, "pca")[which(s2e$condition=="Salmonella"),1:20]
control_pca <- Embeddings(s2e, "pca")[which(s2e$condition=="Control"),1:20]
# Find the nearest neighbors in the control sample for each Salmonella cell
nn <- get.knnx(data = control_pca, query = salmonella_pca, k = 1)
# Create a data frame to store the results
nn_df <- data.frame(
    Salmonella_Cell = rownames(salmonella_pca),
    Control_Cell = rownames(control_pca)[nn$nn.index],
    Distance = nn$nn.dist
)
# Find the nearest neighbors in the infected sample for each control cell
nn_2 <- get.knnx(data = salmonella_pca, query = control_pca, k = 1)
# Create a data frame to store the results
nn_df_2 <- data.frame(
    Salmonella_Cell = rownames(salmonella_pca)[nn_2$nn.index],
    Control_Cell = rownames(control_pca),
    Distance = nn_2$nn.dist
)
df_out<- data.frame(unique(rbind(nn_df, nn_df_2)))
# Save the results
saveRDS(df_out, file = "Res_Salmonella_and_Control_nearest_neighbors_PCA_based.rds")

# get the cluster information of the matched cells
n1_pair <- sapply(sort(unique(s2e$RNA_snn_res.1)), function(control_cl){
    sapply(sort(unique(s2e$RNA_snn_res.1)), function(salmonella_cl){
        sum(df_out$Salmonella_Cell%in%colnames(s2e)[which(s2e$RNA_snn_res.1==salmonella_cl & s2e$condition=="Salmonella")] & df_out$Control_Cell%in%colnames(s2e)[which(s2e$RNA_snn_res.1==control_cl & s2e$condition=="Control")])
    })
})
colnames(n1_pair) <- paste0("C", sort(unique(s2e$RNA_snn_res.1)))
rownames(n1_pair) <- paste0("S", sort(unique(s2e$RNA_snn_res.1)))


n1_salmonella_cells <- sapply(sort(unique(s2e$RNA_snn_res.1)), function(control_cl){
    sapply(sort(unique(s2e$RNA_snn_res.1)), function(salmonella_cl){
        idx <- which(df_out$Salmonella_Cell%in%colnames(s2e)[which(s2e$RNA_snn_res.1==salmonella_cl & s2e$condition=="Salmonella")] & df_out$Control_Cell%in%colnames(s2e)[which(s2e$RNA_snn_res.1==control_cl & s2e$condition=="Control")])
        length(unique(df_out$Salmonella_Cell[idx]))
    })
})
colnames(n1_salmonella_cells) <- paste0("C", sort(unique(s2e$RNA_snn_res.1)))
rownames(n1_salmonella_cells) <- paste0("S", sort(unique(s2e$RNA_snn_res.1)))


n1_control_cells <- sapply(sort(unique(s2e$RNA_snn_res.1)), function(control_cl){
    sapply(sort(unique(s2e$RNA_snn_res.1)), function(salmonella_cl){
        idx <- which(df_out$Salmonella_Cell%in%colnames(s2e)[which(s2e$RNA_snn_res.1==salmonella_cl & s2e$condition=="Salmonella")] & df_out$Control_Cell%in%colnames(s2e)[which(s2e$RNA_snn_res.1==control_cl & s2e$condition=="Control")])
        length(unique(df_out$Control_Cell[idx]))
    })
})
colnames(n1_control_cells) <- paste0("C", sort(unique(s2e$RNA_snn_res.1)))
rownames(n1_control_cells) <- paste0("S", sort(unique(s2e$RNA_snn_res.1)))


c7_cells <- colnames(s2e)[which(s2e$RNA_snn_res.1==7 & s2e$condition=="Control")]
s7_cells <- colnames(s2e)[which(s2e$RNA_snn_res.1==7 & s2e$condition=="Salmonella")]

pdf("Plot_heatmap_matched_s2e_cell_pairs_between_condition_pca_based.pdf", height=7, width=7)
gplots::heatmap.2(n1_pair, scale="row", cellnote=round(n1), notecol="#202020", trace="none", Rowv = FALSE, Colv = FALSE, dendrogram="none", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10))
dev.off()

# for each of the salmonella sample cluster, get the DEGs between the matched control cells
# get the condition composition per cluster
n1 <- sapply(sort(unique(s2e$RNA_snn_res.1)), function(i){
    sapply(sort(unique(s2e$condition)), function(j){
        sum(s2e$condition[s2e$RNA_snn_res.1==i]==j)
    })
})
colnames(n1) <- paste0("C", sort(unique(s2e$RNA_snn_res.1)))

# check how many control cells are mapped to each of the Samonella cluster
matched_cell_num_per_sal_cluster <- t(sapply(sort(unique(s2e$RNA_snn_res.1)), function(i){
    idx <- which(df_out$Salmonella_Cell%in%colnames(s2e)[which(s2e$RNA_snn_res.1==i)])
    n_control <- length(unique(df_out$Control_Cell[idx]))
    n_salmonella_matched <- length(unique(df_out$Salmonella_Cell[idx]))
    n_salmonella_total <- sum(s2e$RNA_snn_res.1==i & s2e$condition=="Salmonella")
    return(c(n_control, n_salmonella_matched, n_salmonella_total))
}))
colnames(matched_cell_num_per_sal_cluster) <- c("Control", "Salmonella_matched", "Salmonella_total")
rownames(matched_cell_num_per_sal_cluster) <- paste0("C", sort(unique(s2e$RNA_snn_res.1)))
unique(s2e@meta.data[,c("RNA_snn_res.1","Subtype")]) -> df
table(df$Subtype[df$RNA_snn_res.1%in%sub("C", "", rownames(matched_control_cells)[which(matched_control_cells[,1]>0)])])

ratio_cell_based_matching <- matched_cell_num_per_sal_cluster[,1]/matched_cell_num_per_sal_cluster[,3]
ratio_cluster_based_matching <- n1[1,]/n1[2,]
abs_ratio_cell_based_matching <- abs(ratio_cell_based_matching-1)
abs_ratio_cluster_based_matching <- abs(ratio_cluster_based_matching-1)
ratio_diff <- abs_ratio_cell_based_matching-abs_ratio_cluster_based_matching
cl_to_check <- names(which(ratio_diff>0))
ratio_cell_based_matching[cl_to_check]
ratio_cluster_based_matching[cl_to_check]

i <- 14
salmonella_cells <- colnames(s2e)[which(s2e$RNA_snn_res.1==i & s2e$condition=="Salmonella")]
matched_control_cells <- unique(df_out$Control_Cell[which(df_out$Salmonella_Cell%in%salmonella_cells)])
p1 <- SCpubr::do_DimPlot(s2e, reduction="umap", cells.highlight=matched_control_cells)
p2 <- SCpubr::do_DimPlot(s2e, reduction="umap", cells.highlight=salmonella_cells)
p <- p1+p2

# perform DE between each Salmonella cluster and its matched control cells
# clean up the genes without gene symbol
detected_genes <- rownames(s2e)[which(rowSums(s2e@assays$RNA@data)>0)]
gene <- detected_genes[!grepl("^ENSG", detected_genes)]
#remove ribosomal genes, mitochondrial genes, and genes located on sex chromosomes before normalization
confound_genes <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/confound_genes/Dat_confound_genes.rds")
genes_to_exclude <- unique(unlist(confound_genes[c(4,5,6)]))
genes <- setdiff(gene, genes_to_exclude)
saveRDS(genes, file="List_input_genes_for_s2E_DEG.rds")
de_res_list <- list()
for(i in sort(unique(s2e$RNA_snn_res.1))){
    print(i)
    salmonella_cells <- colnames(s2e)[which(s2e$RNA_snn_res.1==i & s2e$condition=="Salmonella")]
    matched_control_cells <- unique(df_out$Control_Cell[which(df_out$Salmonella_Cell%in%salmonella_cells)])
    cells <- c(salmonella_cells, matched_control_cells)
    X <- s2e@assays$RNA@data[genes, cells]
    y <- s2e@meta.data[cells, "condition"]
    de_res <- presto::wilcoxauc(X=X, y=y)
    de_res$group <- paste(paste0("C", i), de_res$group, sep="_")
    de_res$pct_diff <- de_res$pct_in - de_res$pct_out
    de_res$pval_log <- -log10(de_res$pval)
    de_res$pval_log[is.infinite(de_res$pval_log)] <- max(de_res$pval_log[!is.infinite(de_res$pval_log)])+1
    de_res$pval_transformed <- de_res$pval_log^0.25
    de_res$auc_sym <- abs(de_res$auc-0.5)
    de_res_list[[i]] <- de_res
    
}
combined_de_res <- do.call('rbind', de_res_list)
sig_res <- combined_de_res %>% filter(padj<0.05 & pct_in>10 & logFC>0.1)
saveRDS(combined_de_res, file="Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_matched_control_cells_DEG_all_res.rds")
saveRDS(sig_res, file="Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_matched_control_cells_DEG_sig_res.rds")

# for each of the Salmonella cluster with sufficient co-clustered control cells, perform DE analysis
selected_cl <- sub("C", "", names(which(n1[1,]>50)))
de_res_list <- list()
for(i in selected_cl){
    print(i)
    salmonella_cells <- colnames(s2e)[which(s2e$RNA_snn_res.1==i & s2e$condition=="Salmonella")]
    matched_control_cells <- colnames(s2e)[which(s2e$RNA_snn_res.1==i & s2e$condition=="Control")]
    cells <- c(salmonella_cells, matched_control_cells)
    X <- s2e@assays$RNA@data[genes, cells]
    y <- s2e@meta.data[cells, "condition"]
    de_res <- presto::wilcoxauc(X=X, y=y)
    de_res$group <- paste(paste0("C", i), de_res$group, sep="_")
    de_res$pct_diff <- de_res$pct_in - de_res$pct_out
    de_res$pval_log <- -log10(de_res$pval)
    de_res$pval_log[is.infinite(de_res$pval_log)] <- max(de_res$pval_log[!is.infinite(de_res$pval_log)])+1
    de_res$pval_transformed <- de_res$pval_log^0.25
    de_res$auc_sym <- abs(de_res$auc-0.5)
    de_res_list[[i]] <- de_res
    
}
cocluster_de_res <- do.call('rbind', de_res_list)
cocluster_sig_res <- cocluster_de_res %>% filter(padj<0.05 & pct_in>10 & logFC>0.1)
saveRDS(cocluster_de_res, file="Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_cocluster_control_cells_DEG_all_res.rds")
saveRDS(cocluster_sig_res, file="Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_cocluster_control_cells_DEG_sig_res.rds")

cell_based_res <- table(sig_res$group)
cluster_based_res <- table(cocluster_sig_res$group)
cell_based_res_subset <- cell_based_res[names(cluster_based_res)]
df <- data.frame('group'=names(cluster_based_res), "cell_based_res"=as.numeric(cell_based_res_subset), "cluster_based_res"=as.numeric(cluster_based_res))
df$condition <- sapply(df$group, function(x){
    strsplit(x, "_")[[1]][2]
})
p1 <- ggplot(df, aes(x=cell_based_res, y=cluster_based_res, color=condition))+
geom_point()+
geom_abline(intercept=0, slope=1)+
theme_minimal()
ggsave(p1, file="Plot_DEG_comparison_between_cell_based_and_cluster_based.pdf", width=5, height=5)

cell_based_degs <- unique(sig_res$feature)
saveRDS(cell_based_degs, file="Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_matched_control_cells_DEG_sig_genes.rds")


no_control_cl <- paste(names(which(n1[1,]<50)),collapse="|")
selected_cell_based_sig_res <- sig_res[grep(no_control_cl, sig_res$group),]
combined_sig_res <- rbind(cocluster_sig_res, selected_cell_based_sig_res)
combined_degs <- unique(combined_sig_res$feature)
saveRDS(combined_sig_res, file="Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_matched_control_cells_and_cocluster_control_cells_DEG_sig_res.rds")
saveRDS(combined_degs, file="Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_matched_control_cells_and_cocluster_control_cells_DEG_sig_genes.rds")


############################################
cocluster_de_res, file="Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_cocluster_control_cells_DEG_all_res.rds"
combined_de_res, file="Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_matched_control_cells_DEG_all_res.rds"
s2e <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/Res_salmonella_data_s2e_with_DM.rds")
n1 <- sapply(sort(unique(s2e$RNA_snn_res.1)), function(i){
    sapply(sort(unique(s2e$condition)), function(j){
        sum(s2e$condition[s2e$RNA_snn_res.1==i]==j)
    })
})
colnames(n1) <- paste0("C", sort(unique(s2e$RNA_snn_res.1)))
no_control_cl <- paste(names(which(n1[1,]<50)),collapse="|")
selected_cell_based_sig_res <- sig_res[grep(no_control_cl, sig_res$group),]
combined_sig_res <- rbind(cocluster_sig_res, selected_cell_based_sig_res)
############################################

combined_top_res <- combined_sig_res %>% group_by(group) %>% top_n(50, wt=logFC)
combined_top_degs <- unique(combined_top_res$feature)
intersect(combined_top_degs, c("LTB", "CXCL8", "IL1B"))

# get cluster average expression per condition
expr_mat_list <- list()
for(x in sort(unique(s2e$condition))){
    seu_obj <- subset(s2e, subset=condition==x)
    expr_mat <- getAveExpr(seu.obj=seu_obj, feature.to.calc="RNA_snn_res.1", colname.prefix=x, size.cutoff=50)
    expr_mat_list[[x]] <- expr_mat
}
expr_mat_control <- expr_mat_list[["Control"]]
no_control_cl_vec <- names(which(n1[1,]<50))
expr_mat_2 <- sapply(no_control_cl_vec, function(x){
    i <- sub("C", "", x)
    salmonella_cells <- colnames(s2e)[which(s2e$RNA_snn_res.1==i & s2e$condition=="Salmonella")]
    matched_control_cells <- unique(df_out$Control_Cell[which(df_out$Salmonella_Cell%in%salmonella_cells)])
    rowMeans(s2e@assays$RNA@data[, matched_control_cells])
})
colnames(expr_mat_2) <- sub("C", "Control", colnames(expr_mat_2))
expr_mat_control <- cbind(expr_mat_control, expr_mat_2)
expr_mat_salmonella <- expr_mat_list[["Salmonella"]] 
expr_mat_control <- expr_mat_control[,sub("Salmonella", "Control", colnames(expr_mat_salmonella))]   
combined_expr_mat <- cbind(expr_mat_salmonella, expr_mat_control)
group_vec <- rep(c("Salmonella", "Control"), each=ncol(expr_mat_salmonella)) 
time_vec <- rep(seq(ncol(expr_mat_salmonella)), 2)
colnames(combined_expr_mat) <- paste(group_vec, time_vec, sep="_") 
saveRDS(combined_expr_mat, file="Dat_combined_cluster_average_expression_per_condition.rds")

combined_expr_mat <- readRDS("Dat_combined_cluster_average_expression_per_condition.rds")
# identify the DEGs with differential change magnitude along the trajectory
diff_expr_mat <- combined_expr_mat[,grep("Salmonella", colnames(combined_expr_mat))] - combined_expr_mat[,grep("Control", colnames(combined_expr_mat))]
#diff_expr_mat <- expr_mat_salmonella - expr_mat_control
cl_enriched_de_idx <- apply(diff_expr_mat[combined_degs,], 1, function(vec){
    idx <- sapply(seq(length(vec)), function(i){
        vec[i]>mean(vec[-i])+3*sd(vec[-i]) | vec[i]<mean(vec[-i])-3*sd(vec[-i])
    })
    sum(idx>0)
})
genes <- names(cl_enriched_de_idx)[cl_enriched_de_idx>0]
length(genes)


#s2e$Cluster_per_condition <- paste(s2e$RNA_snn_res.1, s2e$condition, sep="_")
#SCpubr::do_BoxPlot(s2e, group.by="Cluster_per_condition", feature="IL1B")
# exclude genes with ".", which are novel transcripts
genes_to_exclude <- grep(".", fixed=TRUE, genes, value=T)
genes_2 <- setdiff(genes, genes_to_exclude)
# get the top 200 DEGs with the largest change magnitude
# change magnitude should be larger than 0.3
max_mag <- apply(abs(diff_expr_mat[genes_2,]), 1, max)
genes_top <- names(max_mag[which(max_mag>0.3)])
s2e_degs <- list("all"=combined_degs, "cl_enriched"=genes, "top"=genes_top)
saveRDS(s2e_degs, file="Res_Salmonella_infection_s2E_condition_dependent_DEG.rds")

s2e_degs <- readRDS("Res_Salmonella_infection_s2E_condition_dependent_DEG.rds")
#scaled_expr_mat <- t(scale(t(combined_expr_mat)))
#diff_scaled_expr_mat <- scaled_expr_mat[,grep("Salmonella", colnames(combined_expr_mat))] - scaled_expr_mat[,grep("Control", colnames(combined_expr_mat))]
#expr_mat <- diff_scaled_expr_mat[genes_top,]
expr_mat <- combined_expr_mat[genes_top,]
saveRDS(expr_mat, file="Dat_salmonella_data_s2e_DEG_expression_change.rds")


highlight_genes <- c("LTB", "CXCL8", "IL1B")

hc.row <- hclust(as.dist(1-cor(t(expr_mat))), method="ward.D2")
row.orders <- hc.row$labels[hc.row$order] 
expr.mat <- expr_mat[row.orders,]
max.ct <- colnames(expr.mat)[apply(expr.mat, 1, which.max)]
idx <- unlist(lapply(colnames(expr.mat), function(x){
  which(max.ct==x)
}))
expr.mat <- expr.mat[idx,]
normed.expr <- t(apply(expr.mat, 1, function(vec){(vec-min(vec))/(max(vec)-min(vec))}))
highlight.idx <- ifelse(rownames(normed.expr)%in%highlight_genes, 1, NA)
input <- cbind(normed.expr, highlight.idx)
highlight.input <- input
rownames(highlight.input)[which(is.na(highlight.idx))] <- NA

pdf("Plot_heatmap_topDEG_expr.pdf", height=10)
gplots::heatmap.2(highlight.input, col = darkBlue2Red.heatmap, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 1, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()


combined_expr_mat <- readRDS("Dat_combined_cluster_average_expression_per_condition.rds")
s2e_degs <- readRDS("Res_Salmonella_infection_s2E_condition_dependent_DEG.rds")
all_degs <- s2e_degs$all

# perform gene set enrichment analysis on the DEGs
library(msigdbr)
h_genes <- msigdbr(species = "Homo sapiens", category = "H")
h_gene_list <- tapply(h_genes$gene_symbol, h_genes$gs_name, list)
kegg_genes <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
kegg_gene_list <- tapply(kegg_genes$gene_symbol, kegg_genes$gs_name, list)
wikipathway_genes <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS")
wikipathway_gene_list <- tapply(wikipathway_genes$gene_symbol, wikipathway_genes$gs_name, list)
merged_list <- c(h_gene_list, kegg_gene_list, wikipathway_gene_list)
saveRDS(merged_list , file="Res_KEGG_wikiPathway_hallMark_merged_annotated_gene_list.rds")

merged_df <- rbind(h_genes, kegg_genes, wikipathway_genes)
saveRDS(merged_df , file="Res_KEGG_wikiPathway_hallMark_merged_annotated_gene_list_combined_df.rds")

# perform enrichment analysis on all DEGs
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_statiscal_test.R")
background_genes <- readRDS("List_input_genes_for_s2E_DEG.rds")
pathway_enrichment_pval <- enrichment_test(query_gene_list=s2e_degs, anno_gene_list=merged_list, all_genes=background_genes)
saveRDS(pathway_enrichment_pval, file="Res_salmonella_infection_union_DEG_KEGG_hallmark_wikiPathway_hypergeometric_test_nominal_pval.rds")

# perform enrichment analysis per cluster, focusing on Salmonella infection upregulated genes
combined_sig_res <- readRDS("Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_matched_control_cells_and_cocluster_control_cells_DEG_sig_res.rds")
combined_up_res <- combined_sig_res[grep("Salmonella", combined_sig_res$group),]
up_gene_list <- tapply(combined_up_res$feature, combined_up_res$group, list)
pathway_enrichment_pval <- enrichment_test(query_gene_list=up_gene_list, anno_gene_list=merged_list, all_genes=background_genes)
saveRDS(pathway_enrichment_pval, file="Res_salmonella_infection_per_cluster_DEG_KEGG_hallmark_wikiPathway_hypergeometric_test_nominal_pval.rds")

# perform p-value correction per gene set
idx_kegg <- grep("KEGG", rownames(pathway_enrichment_pval))
idx_wp <- grep("WP", rownames(pathway_enrichment_pval))
idx_hallmark <- grep("HALLMARK", rownames(pathway_enrichment_pval))
padj_list <- lapply(c("HALLMARK", "KEGG", "WP"), function(x){
    idx <- grep(x, rownames(pathway_enrichment_pval))
    pval <- pathway_enrichment_pval[idx,]
    apply(pval, 2, p.adjust, method="BH")
})
padj_mat <- do.call('rbind', padj_list)
#saveRDS(padj_mat, file="Res_salmonella_infection_union_DEG_KEGG_hallmark_wikiPathway_hypergeometric_test_BH_padj.rds")
saveRDS(padj_mat, file="Res_salmonella_infection_per_cluster_DEG_KEGG_hallmark_wikiPathway_hypergeometric_test_BH_padj.rds")
# terms or genes to ge highlighted
# Chemokines, cytokines, interleukins ligands and receptors
# CXCL1, CXCL2, CXCL3, CXCL4, CXCL8; CCL20, CCL28; IL1A, IL1B 
# CXCR4, IL1R2, IL6R
# NFKB pathway
# TNF gene
# bacterial response
# Toll like receptor pathway
# TLR3
# TFFs, antimicrobial peptides, mucins
# LCN2, SAA1,SAA2, REG - antimicrobial 
# DUOX2 - ROS
# REG1A
# autophage
# immune response to intracellular bacteria - GBP1
#unique(merged_df$gs_name[merged_df$gene_symbol%in%"LTB"])
key_words <- c("SIGNALING", "CYTOKINE", "CHEMOKINE", "INTERLEUKIN", "RECEPTOR", "LIGAND", "NFKB", "TNF", "BACTERIAL", "TOLL", "TFF", "ANTIMICROBIAL", "MUCIN", "AUTOPHAGY", "IMMUNE", "INFECTION")
#res <- padj_mat[grepl(paste(key_words, collapse="|"), rownames(padj_mat)) & padj_mat[, "all"]<0.05 & grepl("KEGG_", rownames(padj_mat)),]
sig_terms <- rownames(padj_mat)[which(rowSums(padj_mat<0.05)>0)]
res <- padj_mat[grepl(paste(key_words, collapse="|"), rownames(padj_mat)) & rownames(padj_mat)%in%sig_terms & grepl("KEGG_", rownames(padj_mat)),]
selected_gs_names <- rownames(res)

# plot logFC selected gene set expression patterns across clusters
combined_expr_mat <- readRDS("Dat_combined_cluster_average_expression_per_condition.rds")
expr_diff_mat <- combined_expr_mat[,grep("Salmonella", colnames(combined_expr_mat))] - combined_expr_mat[,grep("Control", colnames(combined_expr_mat))]
union_up_genes <- unique(combined_up_res$feature)
#score_mat <- t(sapply(selected_gs_names, function(x){
#    term_degs <- intersect(merged_df$gene_symbol[which(merged_df$gs_name==x)], union_up_genes)
#    apply(expr_diff_mat[term_degs,],2,mean)
#    
#    
#}))
#
#kegg_sig_terms <- grep("KEGG", sig_terms, value=T)
selected_gs_names <- sig_terms
score_mat <- t(sapply(selected_gs_names, function(x){
    term_degs <- intersect(merged_df$gene_symbol[which(merged_df$gs_name==x)], union_up_genes)
    #apply(expr_diff_mat[term_degs,],2,mean)
    apply(abs(expr_diff_mat[term_degs,]),2,mean)
}))

# show BH-adjusted p-value of selected terms
#values <- -log10(padj_mat[selected_gs_names, "all"])
values <- -log10(apply(padj_mat[selected_gs_names, ], 1, min))
gs_sorted <- names(values)[order(values)]

df <- unique(s2e@meta.data[,c("RNA_snn_res.1", "Subtype")])
df$new_cl_idx <- rank(df$RNA_snn_res.1)
df_2 <- df[order(df$RNA_snn_res.1),]
df_2$id_2 <- paste("Salmonella", df_2$new_cl_idx, sep="_")
col_sidebar_id <- df_2$Subtype
cols <- readRDS("~/Bacteria_TRM/used_object/cols.rds")
# get the median pt value per cluster
median_pt_per_cluster <- sapply(sort(unique(s2e$RNA_snn_res.1)), function(i){
    median(s2e$dc1_rank[which(s2e$RNA_snn_res.1==i)])
})
colorPal <- grDevices::colorRampPalette(c("#f7f7f7", "#cccccc", "#969696", "#636363", "#252525"))
cellColor <- adjustcolor(colorPal(30), alpha = 0.8)[as.numeric(cut(median_pt_per_cluster, 
            breaks = 30, right = F, include.lowest = T))]
#pdf("Plot_selected_per_cluster_enriched_KEGG_pathway_DEG_logFC_across_clusters.pdf", height=7, width=8)
#gplots::heatmap.2(score_mat[gs_sorted,], scale="row", trace="none", dendrogram="none", col=darkBlue2Red.heatmap, ColSideColors=cellColor, density.info="none", Colv=FALSE, Rowv=FALSE, keysize=0.8, mar=c(5,30))
#dev.off()

pdf("Plot_selected_per_cluster_enriched_KEGG_pathway_DEG_logFC_across_clusters_2.pdf", height=7, width=8)
gplots::heatmap.2(score_mat[gs_sorted,], scale="row", trace="none", dendrogram="row", col=darkBlue2Red.heatmap, ColSideColors=cols$cell_type[df_2$Subtype], density.info="none", Colv=FALSE, Rowv=FALSE, keysize=0.8, mar=c(5,30))
dev.off()

pdf("Plot_all_per_cluster_enriched_pathway_DEG_logFC_across_clusters_2.pdf", height=7, width=8)
gplots::heatmap.2(score_mat[gs_sorted,], scale="row", trace="none", dendrogram="row", col=darkBlue2Red.heatmap, ColSideColors=cols$cell_type[df_2$Subtype], density.info="none", Colv=FALSE, Rowv=TRUE, keysize=0.8, mar=c(5,30))
dev.off()


pdf("Plot_barplot_selected_per_cluster_terms_BH_padj.pdf", height=5, width=10)
barplot(values[gs_sorted], las=2)
dev.off()

combined_expr_mat <- readRDS("Dat_combined_cluster_average_expression_per_condition.rds")
expr_diff_mat <- combined_expr_mat[,grep("Salmonella", colnames(combined_expr_mat))] - combined_expr_mat[,grep("Control", colnames(combined_expr_mat))]
combined_sig_res <- readRDS("Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_matched_control_cells_and_cocluster_control_cells_DEG_sig_res.rds")
combined_up_res <- combined_sig_res[grep("Salmonella", combined_sig_res$group),]
union_up_genes <- unique(combined_up_res$feature)
padj_mat <- readRDS("Res_salmonella_infection_per_cluster_DEG_KEGG_hallmark_wikiPathway_hypergeometric_test_BH_padj.rds")
sig_terms <- rownames(padj_mat)[which(rowSums(padj_mat<0.05)>0)]
selected_gs_names <- sig_terms
merged_df <- readRDS("Res_KEGG_wikiPathway_hallMark_merged_annotated_gene_list_combined_df.rds")
score_mat <- t(sapply(selected_gs_names, function(x){
    term_degs <- intersect(merged_df$gene_symbol[which(merged_df$gs_name==x)], union_up_genes)
    #apply(expr_diff_mat[term_degs,],2,mean)
    apply(abs(expr_diff_mat[term_degs,]),2,mean)
}))
df <- unique(s2e@meta.data[,c("RNA_snn_res.1", "Subtype")])
df$new_cl_idx <- rank(df$RNA_snn_res.1)
df_2 <- df[order(df$RNA_snn_res.1),]
df_2$id_2 <- paste("Salmonella", df_2$new_cl_idx, sep="_")
col_sidebar_id <- df_2$Subtype
cols <- readRDS("~/Bacteria_TRM/used_object/cols.rds")
# summarize to per subtype - for subtype with more than one cluster, take the maximum 
subtype_score_mat <- sapply(c("Colon_SC", paste("Colonocyte", seq(4), sep="_")), function(i){
    idx <- which(df_2$Subtype==i)
    if(length(idx)>1){
        apply(score_mat[,idx], 1, max)
    }else{
        score_mat[,idx]
    }
})
saveRDS(subtype_score_mat, file="Dat_pathway_union_up_DEG_absolute_logFC_per_s2E_subtype.rds")
res <- gplots::heatmap.2(subtype_score_mat, scale="row", trace="none", dendrogram="row", col=darkBlue2Red.heatmap, ColSideColors=cols$cell_type[colnames(subtype_score_mat)], density.info="none", Colv=FALSE, keysize=0.8, mar=c(5,30))
saveRDS(res, file="Res_heatmap_input_pathway_union_up_DEG_absolute_logFC_per_s2E_subtype.rds")
pdf('Plot_heatmap_enriched_pathway_union_upregulated_DEG_absolute_logFC_per_subtype.pdf', height=10, width=10)
gplots::heatmap.2(subtype_score_mat, scale="row", trace="none", dendrogram="row", col=darkBlue2Red.heatmap, ColSideColors=cols$cell_type[colnames(subtype_score_mat)], density.info="none", Colv=FALSE, keysize=0.8, mar=c(10,30))
dev.off()

max_subtype <- setNames(colnames(subtype_score_mat)[apply(subtype_score_mat, 1, which.max)], rownames(subtype_score_mat))
gs_max_logFC_subtype <- lapply(c("Colon_SC", paste("Colonocyte", seq(4), sep="_")), function(i){
    names(max_subtype)[which(max_subtype==i)]
})
names(gs_max_logFC_subtype) <- c("Colon_SC", paste("Colonocyte", seq(4), sep="_"))  
saveRDS(gs_max_logFC_subtype, file="Res_pathway_union_up_DEG_max_logFC_subtype.rds")
id <- grep("CYTOSKELETON", names(max_subtype), value=T)
gplots::heatmap.2(subtype_score_mat[id,], scale="row", trace="none", dendrogram="row", col=darkBlue2Red.heatmap, ColSideColors=cols$cell_type[colnames(subtype_score_mat)], density.info="none", Colv=FALSE, keysize=0.8, mar=c(10,30))


#saveRDS(res, file="Res_heatmap_selected_enriched_KEGG_pathway_DEG_logFC_across_clusters.rds")

# check the change magnitude per cluster
combined_expr_mat <- readRDS("Dat_combined_cluster_average_expression_per_condition.rds")
expr_diff_mat <- combined_expr_mat[,grep("Salmonella", colnames(combined_expr_mat))] - combined_expr_mat[,grep("Control", colnames(combined_expr_mat))]
s2e_degs <- readRDS("Res_Salmonella_infection_s2E_condition_dependent_DEG.rds")
all_degs <- s2e_degs$all
combined_sig_res <- readRDS("Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_matched_control_cells_and_cocluster_control_cells_DEG_sig_res.rds")
up_genes <- unique(combined_sig_res$feature[grep("Salmonella", combined_sig_res$group)])
down_genes <- unique(combined_sig_res$feature[grep("Control", combined_sig_res$group)])

pdf("Plot_change_magnitude_per_cluster_for_union_DEGs.pdf", height=5, width=15)
par(mfrow=c(1,3))
boxplot(sqrt(log(abs(expr_diff_mat[all_degs,])+1)), las=2, main="All DEGs", ylab="Change magnitude", col=cols$cell_type[df_2$Subtype])
boxplot(sqrt(log(abs(expr_diff_mat[up_genes,])+1)), las=2, main="Union infection-up DEGs", ylab="Change magnitude", col=cols$cell_type[df_2$Subtype])
boxplot(sqrt(log(abs(expr_diff_mat[down_genes,])+1)), las=2, main="Union infection-down DEGs", ylab="Change magnitude", col=cols$cell_type[df_2$Subtype])
dev.off()

# get the control cells that are mapped to C17 and C14 Salmonella cells
matched_control_cells <- sapply(c(14, 17), function(i){
    idx <- which(df_out$Salmonella_Cell%in%colnames(s2e)[which(s2e$RNA_snn_res.1==i)])
    control_cells <- unique(df_out$Control_Cell[idx])
    sapply(sort(unique(s2e$RNA_snn_res.1)), function(j){
        sum(control_cells%in%colnames(s2e)[which(s2e$RNA_snn_res.1==j)])
    })
})
colnames(matched_control_cells) <- c("C14", "C17")
rownames(matched_control_cells) <- paste0("C", sort(unique(s2e$RNA_snn_res.1)))

detected_genes <- rownames(s2e)[which(rowSums(s2e@assays$RNA@data)>0)]
pt_expr_list <- list()
s2e$pt_bin_per_condition <- NA
for(x in sort(unique(s2e$condition))){
    idx <- which(s2e$condition==x)
    X <- s2e@assays$RNA@data[detected_genes, idx]
    res <- getExprByPt(pt.vec=s2e$dc1_rank[idx], expr.mat=X, mode="fix.bin.num", bin.num=30, return.idx=T)
    s2e@meta.data[colnames(X), "pt_bin_per_condition"] <- paste(x, res$cell.idx.vec, sep='_')
    expr_mat <- res$expr.mat
    colnames(expr_mat) <- paste(x, seq(ncol(expr_mat)), sep="_")
    pt_expr_list[[x]] <- expr_mat
}
SCpubr::do_BoxPlot(s2e, group.by="pt_bin_per_condition", feature="IL1B")
SCpubr::do_BoxPlot(s2e, group.by="RNA_snn_res.1", feature="IL1B")
saveRDS(pt_expr_list, file="Res_salmonella_data_s2e_Pt_bin_expression_per_condition.rds")
#saveRDS(s2e, file="Res_salmonella_data_s2e_with_DM.rds")
combined_pt_expr <- do.call('cbind', pt_expr_list)


# get the top DEGs of stem cells, colonocyte subtypes and goblet cells
res <- top_deg_res_per_ct_combined[grep("Colonocyte|SC|Goblet", top_deg_res_per_ct_combined$group),]
features <- unique(res$feature)
ave_expr <- getAveExpr(seu.obj=combined, feature.to.calc="Subtype_per_condition", colname.prefix=NULL)
saveRDS(ave_expr, file="Res_salmonella_dataset_subtype_average_expr.rds")
ave_expr_subset <- ave_expr[,colnames(ave_expr)%in%res$group]
expr_mat <- ave_expr[features,colnames(ave_expr)%in%res$group]

res <- prepareTreeAndHeatmapInput(expr.mat = expr_mat, hc.method = "ward.D2", norm.method = "quantile")
pdf("Plot_heatmap_topDEG_expr.pdf", height=10)
gplots::heatmap.2(res$heatmap_input, col = beach.col.heatmap, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 1, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

# generate scatter plot in the combined dataset
# X shows maximal absolute value of infection - control logFC
# Y shows expression specificity across the combined clusters
# color represents cell type with the maximal logFC
# size represents significance of the DEG
# get the DEG results
condition_deg_res <- list()
for(ct in sort(unique(combined$Cell_type_with_colonocyte_subtypes))){
    
    print(paste(ct, 'start'))
    idx <- which(combined$Cell_type_with_colonocyte_subtypes==ct)
    X <- combined@assays$RNA@data[genes,idx]
    y <- combined$Subtype_per_condition[idx]

    de_res <- presto::wilcoxauc(X=X, y=y)
    de_res$pct_diff <- de_res$pct_in - de_res$pct_out
    pval_log <- -log10(de_res$pval)
    pval_log[is.infinite(pval_log)] <- max(pval_log[!is.infinite(pval_log)])+1
    de_res$pval_transformed <- pval_log^0.25
    de_res$auc_sym <- abs(de_res$auc-0.5)
    de_res <- de_res[grep(":Salmonella", de_res$group),]
    condition_deg_res[[ct]] <- de_res
}
df <- data.frame(do.call('cbind', condition_deg_res))
saveRDS(df, file="Res_Salmonella_infection_on_colon_epi_per_cell_type_DEG_all_res_combined.rds")



df_logFC <- df[,grep("logFC",colnames(df))]
max_abs_logFC_idx <- apply(abs(df_logFC), 1, which.max)
max_abs_logFC <- sapply(seq(nrow(df_logFC)), function(i){
    df_logFC[i, max_abs_logFC_idx[i]]
})
max_abs_logFC_ct <- sub(".logFC", "", colnames(df_logFC)[max_abs_logFC_idx])
df_pval <- df[,grep("pval_transformed", colnames(df))]
max_abs_logFC_pval <- sapply(seq(nrow(df_pval)), function(i){
    df_pval[i, max_abs_logFC_idx[i]]
})

cluster_expr <- getAveExpr(seu.obj=combined, feature.to.calc="RNA_snn_res.1")
saveRDS(cluster_expr, file="Res_salmonella_dataset_combined_cluster_expr.rds")
tau_combined <- calc_tau(m=cluster_expr, byRow=T)

input <- data.frame('feature'=df$BEST4._cell.feature, "max_logFC"=max_abs_logFC, "max_logFC_ct"=max_abs_logFC_ct, "max_logFC_pval"=max_abs_logFC_pval, "tau_combined"=tau_combined[df$BEST4._cell.feature])
input$max_logFC_ct[which(input$max_logFC_ct=="BEST4._cell")] <- "BEST4+_cell"
saveRDS(input, file="Res_Salmonella_infection_on_colon_epi_per_cell_type_DEG_max_logFC_and_tau.rds")

p1 <- ggplot(input, aes(x=max_logFC, y=tau_combined, size=max_logFC_pval, fill=max_logFC_ct))+
    geom_point(shape=21, color="#202020")+
    theme_bw()+
    theme(legend.position="bottom")+
    scale_fill_manual(values=cols$cell_type, guide = guide_legend(override.aes = list(label = "")))+
    labs(x="logFC (Infection - Control)", size="-log10(P)^0.25", y="Expression specificity", fill="Cell type")

# Split the plot into multiple plots according to max_logFC_ct
p1_split <- p1 + facet_wrap(~ max_logFC_ct, scales = "free")

p1_split
# quantify the expression enrichment in each cell type in infection condition respectively
idx <- which(combined$condition=="Salmonella")
X <- combined@assays$RNA@data[genes,idx]
y <- combined$Subtype_per_condition[idx]
sal_ct_de_res <- presto::wilcoxauc(X=X, y=y)    
saveRDS(sal_ct_de_res, file="Res_Salmonella_infected_sample_expr_enrichment_per_cell_type.rds")
condition_de_res <- do.call('rbind', condition_deg_res)
colnames(condition_de_res) <- paste("condition", colnames(condition_de_res), sep="_")
colnames(sal_ct_de_res) <- paste("salmonellaCellType", colnames(sal_ct_de_res), sep="_")
df_2 <- data.frame(cbind(condition_de_res, sal_ct_de_res))
df_2$Cell_type <- sub(":Salmonella", "", df_2$condition_group)
saveRDS(df_2, file="Res_per_cell_type_condition_DE_and_salmonella_infected_sample_expr_enrichment.rds")

# for each cell type, highlight the DEGs with either the top 10 logFC between condition or the top 10 expression enrichment in infection condition 
plot_list <- list()
for(x in unique(df_2$Cell_type)){
    print(x)
    idx <- which(df_2$Cell_type==x)
    input <- df_2[idx,]
    de_res <- input[,grep("condition_",colnames(input))]
    sig_res <- de_res %>% filter(condition_padj<0.05 & condition_pct_in>5 & condition_logFC>0.05)
    sig_genes <- sig_res$condition_feature

    pos_genes <- input$condition_feature[order(input$condition_logFC, decreasing=T)[1:10]]
    neg_genes <- input$condition_feature[order(input$condition_logFC, decreasing=F)[1:10]]
    ct_enriched_genes <- input$salmonellaCellType_feature[order(input$salmonellaCellType_logFC, decreasing=T)[1:10]]
    ct_depleted_genes <- input$salmonellaCellType_feature[order(input$salmonellaCellType_logFC, decreasing=F)[1:10]]
    genes <- unique(c(pos_genes, neg_genes, ct_enriched_genes, ct_depleted_genes))
    union_pos_genes <- intersect(genes, input$condition_feature[input$condition_logFC>0])
    union_neg_genes <- intersect(genes, input$condition_feature[input$condition_logFC<0])

    input$group <- "Background"
    input$group[input$condition_feature%in%sig_genes] <- "DEG"
    input$group[input$condition_feature%in%union_pos_genes] <- "Infection up"
    input$group[input$condition_feature%in%union_neg_genes] <- "Infection down"
    
    group_cols <- setNames(c("#a0a0a0", "#696969", "#922b21", "#1f618d"), c("Background", "DEG", "Infection up", "Infection down"))
    p1 <- ggplot(input, aes(x=condition_logFC, y=salmonellaCellType_logFC, size=condition_pval_transformed, alpha=condition_auc_sym, fill=group))+
        geom_point(shape=21, color="#202020")+
        theme_bw()+
        theme(legend.position="bottom")+
        scale_fill_manual(values=group_cols, guide = guide_legend(override.aes = list(label = "")))+
        labs(x="logFC (Infection - Control)", size="-log10(P)^0.25", y="Expression enrichment in infected sample", title=x)+
        geom_text_repel(data=input[input$condition_feature%in%union_pos_genes,], aes(label=condition_feature), color="#922b21", size=4, point.padding=1, segment.color='grey50', min.segment.length=0, box.padding=0.5, max.overlaps=999)+
        geom_text_repel(data=input[input$condition_feature%in%union_neg_genes,], aes(label=condition_feature), color="#1f618d", size=4, point.padding=1, segment.color='grey50', min.segment.length=0, box.padding=0.5, max.overlaps=999)
    plot_list[[x]] <- p1
}

library(cowplot)
library(grid)
library(gridExtra)
library(cowplot)
p <- plot_grid(plotlist = plot_list, align = "hv", nrow = 2, axis = "l")
ggsave(p, filename="Plot_Salmonella_infection_on_colon_epi_subtype_DEG_volcano_plot.png", height=10*2, width=10*length(plot_list)/2)

# generate scatter plot in the combined dataset
# X shows infection - control logFC
# Y shows expression enrichment in infection condition
# size represents significance of the DEG
p1 <- ggplot(df_2, aes(x=condition_logFC, y=salmonellaCellType_logFC, size=condition_pval_transformed, alpha=condition_auc_sym))+
    geom_point()+
    theme_bw()+
    theme(legend.position="bottom")+
    #scale_fill_manual(values=cols$cell_type, guide = guide_legend(override.aes = list(label = "")))+
    labs(x="logFC (Infection - Control)", size="-log10(P)^0.25", y="Expression enrichment in infected sample")

# Split the plot into multiple plots according to max_logFC_ct
p1_split <- p1 + facet_wrap(~ Cell_type, scales = "free")

ggsave(p1_split, filename="Plot_Salmonella_infection_on_colon_epi_per_cell_type_DEG_condition_and_cell_type_logFC.png", width=5*3, height=5*3)




## generate scatter plot per cell type
## with X showing infection - control logFC
## y showing the expression enrichment of logFC in infection condition
# generate vocalno plot per cell type
sig_deg_res_per_ct <- list()
plot_list <- list()
for(ct in sort(unique(combined$Cell_type_with_colonocyte_subtypes))){
    
    print(paste(ct, 'start'))
    idx <- which(combined$Cell_type_with_colonocyte_subtypes==ct)
    X <- combined@assays$RNA@data[genes,idx]
    y <- combined$Ct_per_condition[idx]

    de_res <- presto::wilcoxauc(X=X, y=y)
    de_res$pct_diff <- de_res$pct_in - de_res$pct_out
    de_res$pval_log <- -log10(de_res$pval)
    de_res$pval_log[is.infinite(de_res$pval_log)] <- max(de_res$pval_log[!is.infinite(de_res$pval_log)])+1
    de_res$pval_transformed <- de_res$pval_log^0.25
    de_res$auc_sym <- abs(de_res$auc-0.5)
    sig_res <- de_res %>% filter(padj<0.05 & pct_in>5 & logFC>0.05)
    sig_deg_res_per_ct[[ct]] <- sig_res
    sig_genes <- sig_res$feature

    input <- de_res[grep("Salmonella", de_res$group),]
    
    pos_genes <- input$feature[order(input$logFC, decreasing=T)[1:10]]
    if(ct == "Colonocyte"){
        pos_genes <- c(pos_genes, "IL1B", "CXCL8", "SAA1", "SAA2")
    }
    neg_genes <- input$feature[order(input$logFC, decreasing=F)[1:10]]
    input$group <- "Background"
    input$group[input$feature%in%sig_genes] <- "DEG"
    input$group[input$feature%in%pos_genes] <- "Top Salmonella high"
    input$group[input$feature%in%neg_genes] <- "Top control high"
    group_cols <- setNames(c("#a0a0a0", "#696969", "#922b21", "#1f618d"), c("Background", "DEG", "Top Salmonella high", "Top control high"))
    p1 <- ggplot(input, aes(x=logFC, y=pval_transformed, size=auc_sym, fill=group))+
        geom_point(shape=21, color="#202020")+
        theme_bw()+
        theme(legend.position="bottom")+
        scale_fill_manual(values=group_cols, guide = guide_legend(override.aes = list(label = "")))+
        labs(x="logFC (Infection - Control)", size="|AUC-0.5|", y="-log10(P)^0.25", title=ct)+
        #geom_vline(xintercept=0, color="black", linetype="dashed")+
        geom_text_repel(data=input[input$feature%in%pos_genes,], aes(label=feature), color="#922b21", size=4, point.padding=1, segment.color='grey50', min.segment.length=0, box.padding=0.5, max.overlaps=999)+
        geom_text_repel(data=input[input$feature%in%neg_genes,], aes(label=feature), color="#1f618d", size=4, point.padding=1, segment.color='grey50', min.segment.length=0, box.padding=0.5, max.overlaps=999)+
        theme(text=element_text(family="Arial"))
        #scale_color_manual(values=cols$cell_type, guide = "none")
    plot_list[[ct]] <- p1
 
}
sig_deg_res_per_ct_combined <- do.call('rbind', sig_deg_res_per_ct)
saveRDS(sig_deg_res_per_ct_combined, file="Res_Salmonella_infection_on_colon_epi_subtype_DEG.rds")
top_deg_res_per_ct_combined <- sig_deg_res_per_ct_combined %>% group_by(group) %>% top_n(50, wt=logFC)
saveRDS(top_deg_res_per_ct_combined, file="Res_Salmonella_infection_on_colon_epi_subtype_top_DEG.rds")
union_deg <- unique(top_deg_res_per_ct_combined$feature)


library(cowplot)
library(grid)
library(gridExtra)
library(cowplot)
p <- plot_grid(plotlist = plot_list, align = "hv", nrow = 2, axis = "l")
ggsave(p, filename="Plot_Salmonella_infection_on_colon_epi_subtype_DEG_volcano_plot.png", height=10*2, width=10*length(plot_list)/2)


ave_expr <- readRDS("Res_cell_type_expression_per_condition.rds")
expr_diff <- readRDS("Res_cell_type_expression_diff_Salmonella_over_Control_per_condition.rds")

sd_vec <- apply(expr_diff, 1, sd)
idx <- which(sd_vec>0)
expr_diff <- expr_diff[idx,]
ave_expr <- ave_expr[idx,]
expr_diff_enrichment <- apply(expr_diff, 1, function(vec){
    a <- max(abs(vec))
    b <- mean(vec)
    if(sign(a)==sign(b)){
        enrichment <- abs(a-b)
    }else{
        enrichment <- abs(a)+abs(b)
    }
        
})
max_expr_diff <- apply(expr_diff, 1, function(vec){
    vec[which.max(abs(vec))]
})
max_expr_diff_ct <- sub(":Salmonella", "", colnames(expr_diff)[apply(abs(expr_diff), 1, which.max)])
df <- data.frame("feature"=rownames(expr_diff), "max_logFC"=max_expr_diff, "max_logFC_enrichment"=expr_diff_enrichment, "max_logFC_ct"=max_expr_diff_ct)

# clean up the genes without gene symbol
gene <- df$feature[!grepl("^ENSG", df$feature)]
#remove ribosomal genes, mitochondrial genes, and genes located on sex chromosomes before normalization
confound_genes <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/confound_genes/Dat_confound_genes.rds")
genes_to_exclude <- unique(unlist(confound_genes[c(4,5,6)]))
genes <- setdiff(gene, genes_to_exclude)
df <- df[which(df$feature%in%genes),]

# calculate the cell type specificity expression in control and infection condition respectively
control_tau <- calc_tau(m=ave_expr[,grep("Control",colnames(ave_expr))], byRow=T)
sal_tau <- calc_tau(m=ave_expr[,grep("Salmonella",colnames(ave_expr))], byRow=T)
df$control_tau <- control_tau[df$feature]
df$sal_tau <- sal_tau[df$feature]
df$max_logFC_condition_tau <- ifelse(df$max_logFC>0, df$sal_tau, df$control_tau)
df$max_logFC_ct[which(df$max_logFC_ct=="BEST4+ cell")] <- "BEST4+_cell"
df$max_logFC_ct[which(df$max_logFC_ct=="SC")] <- "Colon_SC"
saveRDS(df, file="Res_Salmonella_infection_on_colon_epi_gene_difference_specificity.rds")

# highlight the top10 DEGs
pos_genes <- df$feature[order(df$max_logFC, decreasing=T)[1:10]]
neg_genes <- df$feature[order(df$max_logFC, decreasing=F)[1:10]]

# generate scatter plot with X and Y axis representing maximal logFC and difference between maximal logFC and average logFC
# with size proportional to the expression specificity in the condition with the maximal logFC 
# and color representing the cell type with the maximal logFC
library(ggrepel)
cols <- readRDS("~/Bacteria_TRM/used_object/cols.rds")
p1 <- ggplot(df, aes(x=max_logFC, y=max_logFC_enrichment, size=max_logFC_condition_tau, fill=max_logFC_ct))+
    geom_point(shape=21, color="#303030")+
    theme_bw()+
    theme(legend.position="bottom")+
    scale_fill_manual(values=cols$cell_type, guide = guide_legend(override.aes = list(label = "")))+
    labs(y="logFC enrichment", x="Maximal logFC (Infection-Control)", size="Cel type specificity", fill="Cell type")+
    geom_vline(xintercept=0, color="black", linetype="dashed")+
    geom_text_repel(data=df[df$feature%in%c(pos_genes,neg_genes),], aes(label=feature, color=max_logFC_ct), size=3, point.padding=1, segment.color='grey50', min.segment.length=0, box.padding=0.3)+
    scale_color_manual(values=cols$cell_type, guide = "none")
p1

ggsave(p1, filename="Plot_Salmonella_infection_on_colon_epi_gene_difference_specificity.pdf", height=10, width=10)



# retrieve the number of scientific publications where a gene and a trait cooccur
ref_num <- setNames(pbapply::pbsapply(union_deg, function(g){
  rentrez::entrez_search(db="pubmed", term=paste(g, "Salmonella", sep=' '))$count
}), union_deg)
df_deg$ref_num_log1p <- log1p(ref_num[df_deg$feature])
saveRDS(df_deg, file="Res_Salmonella_infection_on_colon_epi_topDEG_difference_specificity.rds")

df_deg <- readRDS("Res_Salmonella_infection_on_colon_epi_topDEG_difference_specificity.rds")
df_deg$max_logFC_ct[which(df_deg$max_logFC_ct=="SC")] <- "Colon_SC"
df_deg$max_logFC_ct[which(df_deg$max_logFC_ct=="BEST4+ cell")] <- "BEST4+_cell"
df_deg$ref_num_log1p_sqrt <- sqrt(df_deg$ref_num_log1p)
saveRDS(df_deg, file="Res_Salmonella_infection_on_colon_epi_topDEG_difference_specificity.rds")



# cluster the DEGs
sig_deg_res_per_ct_combined <- readRDS("Res_Salmonella_infection_on_colon_epi_cell_type_DEG.rds")
genes <- unique(sig_deg_res_per_ct_combined$feature)
confound_genes <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/confound_genes/Dat_confound_genes.rds")
genes_to_exclude <- unique(unlist(confound_genes[c(4,5,6)]))
condition_deg <- setdiff(genes, genes_to_exclude)

# cluster the DEG profile
combined$Ct_per_condition <- paste(combined$Coarse_cell_type, combined$condition, sep=":")
saveRDS(combined, file="Dat_salmonella_colon_epithelium_scRNA-seq_data.rds")

ave_expr <- getAveExpr(seu.obj=combined, feature.to.calc="Ct_per_condition", colname.prefix=NULL)
saveRDS(ave_expr, file="Res_cell_type_expression_per_condition.rds")
expr_diff <- ave_expr[,grep("Salmonella",colnames(ave_expr))] - ave_expr[,grep("Control",colnames(ave_expr))]
saveRDS(expr_diff, file="Res_cell_type_expression_diff_Salmonella_over_Control_per_condition.rds")

#expr_mat <- ave_expr[condition_deg,]
#cl_enriched_de_idx <- apply(expr_mat, 1, function(vec){
#    idx <- sapply(seq(length(vec)), function(i){
#        vec[i]>mean(vec[-i])+2*sd(vec[-i]) | vec[i]<mean(vec[-i])-2*sd(vec[-i])
#    })
#    sum(idx>0)
#})
#genes <- names(cl_enriched_de_idx)[cl_enriched_de_idx>0]
#length(genes)
#saveRDS(genes, file="Res_cell_type_enriched_DEG.rds")
#deg_list <- list("all"=condition_deg, "cell_type_enriched"=genes)
#saveRDS(deg_list, file="Res_per_cell_type_wilcox_test_DEG.rds")
input_expr_mat <- ave_expr[condition_deg,]
hc <- hclust(as.dist(1-cor(t(input_expr_mat))), method="ward.D2")
plot(hc, hang=-1)
abline(h=7)
g <- cutree(hc, h=7)
mat <- do.call('rbind', strsplit(split=":", colnames(input_expr_mat)))

res <- plotClusterExprProfile(expr = input_expr_mat, time.vec = mat[,1], time.vec.type="categoric", group.vec = mat[,2], 
 cluster.vec = g, group.cols = c("#31a354", "#3182bd"), return.value = T, to.plot = T, plot.name = "Plot_DEG_cluster_profile.pdf", 
add.legend = T, legend.pos = "topleft", cex.legend = 2, col.num = 3, 
border.do.smooth = F, mean.do.smooth = F, df = 8, xlab = "Cell type", 
ylab = "Relative expr.", sep=":") 
saveRDS(res, file="Res_DEG_cluster_profile.rds")


# generate a scatter plot showing gene set / pathway / GO enrichment 



si <- readRDS("/home/yuq22/Bacteria_TRM/intestine_multi_region_scRNA-seq_atlas/SI_subset/Res_selected_SI_cell_types_with_CSS.rds")

li <- readRDS("/home/yuq22/Bacteria_TRM/intestine_multi_region_scRNA-seq_atlas/non_SI_subset/Res_selected_non-SI_cell_types_with_CSS.rds")


dir.create('goblet')
setwd('/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/goblet')
combined <- readRDS("../Dat_salmonella_colon_epithelium_scRNA-seq_data.rds")
SCpubr::do_FeaturePlot(combined, feature="SOX4", order=T)

goblet <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/RNA_per_condition/tissue_colon/Res_colon_tissue_subset_to_goblet_cells.rds")
# subset to goblet cells of gut cell atlas
goblet_gca <- subset(goblet, subset = paper=="gut_cell_atlas")
SCpubr::do_DimPlot(goblet_gca, reduction="umap_css")

ref_obj <- goblet_gca
que_obj <- subset(combined, subset=Coarse_cell_type=="Goblet")
saveRDS(que_obj, file="Res_chip_goblet_cells.rds")
de_res <- presto::wilcoxauc(ref_obj, group_by="RNA_CSS_snn_res.1", seurat_assay="RNA")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_in>10 & logFC>0.1)
table(sig_res$group)
top_res <- sig_res %>% group_by(group) %>% top_n(100, wt=logFC)
genes <- intersect(top_res$feature, rownames(combined))

tissue_expr <- getAveExpr(seu.obj=ref_obj, feature.to.calc="RNA_CSS_snn_res.1", colname.prefix="tisssue")
chip_expr <-   getAveExpr(seu.obj=que_obj, feature.to.calc="RNA_snn_res.1", colname.prefix="chip")
cor_mat <- cor(tissue_expr[genes,], chip_expr[genes,], method="spearman")
pdf("Plot_heatmap_SCC_goblet_cluster_between_chip_and_tissue.pdf", height=10, width=10)
gplots::heatmap.2(cor_mat, scale="column", cellnote=round(cor_mat,2), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10))
dev.off()
# no big difference between similarity to different tissue goblet cell clusters

# resolve the heterogeneity of goblet cells in the chip data
de_res <- presto::wilcoxauc(que_obj, group_by="RNA_snn_res.1", seurat_assay="RNA")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_in>10 & logFC>0.1)
table(sig_res$group)
top_res <- sig_res %>% group_by(group) %>% top_n(10, wt=logFC)
features <- top_res$feature
p2 <- SCpubr::do_DotPlot(sample = que_obj, 
             features = features,
             scale=F,
             group.by="RNA_snn_res.1", flip=T)

# perform DE test between condition in each of the goblet cell clusters
detected_genes <- rownames(que_obj)[which(rowSums(que_obj@assays$RNA@data)>0)]
gene <- detected_genes[!grepl("^ENSG", detected_genes)]
#remove ribosomal genes, mitochondrial genes, and genes located on sex chromosomes before normalization
confound_genes <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/confound_genes/Dat_confound_genes.rds")
genes_to_exclude <- unique(unlist(confound_genes[c(4,5,6)]))
genes <- setdiff(gene, genes_to_exclude)

que_obj$Cluster_per_condition <- paste(que_obj$condition, que_obj$RNA_snn_res.1, sep="_")
saveRDS(que_obj, file="Res_chip_goblet_cells_with_cluster_per_condition.rds")
deg_res_per_ct <- list()
input_list <- list()
for(ct in sort(unique(que_obj$RNA_snn_res.1))){
    
    print(paste(ct, 'start'))
    idx <- which(que_obj$RNA_snn_res.1==ct)
    X <- que_obj@assays$RNA@data[genes,idx]
    y <- que_obj$Cluster_per_condition[idx]

    de_res <- presto::wilcoxauc(X=X, y=y)
    de_res$pct_diff <- de_res$pct_in - de_res$pct_out
    de_res$pval_log <- -log10(de_res$pval)
    de_res$pval_log[is.infinite(de_res$pval_log)] <- max(de_res$pval_log[!is.infinite(de_res$pval_log)])+1
    de_res$pval_transformed <- de_res$pval_log^0.25
    de_res$auc_sym <- abs(de_res$auc-0.5)
    de_res$Cluster <- paste0("C", ct)
    deg_res_per_ct[[ct]] <- de_res
    sig_res <- de_res %>% filter(padj<0.05 & pct_in>10 & logFC>0.1)
    sig_genes <- sig_res$feature

    input <- de_res[grep("Salmonella", de_res$group),]
    pos_genes <- input$feature[order(input$logFC, decreasing=T)[1:10]]
    neg_genes <- input$feature[order(input$logFC, decreasing=F)[1:10]]
    input$group <- "Background"
    input$group[input$feature%in%sig_genes] <- "DEG"
    input$group[input$feature%in%pos_genes] <- "Top Salmonella high"
    input$group[input$feature%in%neg_genes] <- "Top control high"
    input_list[[ct]] <- input
    
}
deg_res_per_ct_combined <- do.call('rbind', deg_res_per_ct)
saveRDS(deg_res_per_ct_combined, file="Res_Salmonella_infection_on_colon_goblet_cluster_DEG.rds")

library(cowplot)
library(grid)
library(gridExtra)
library(cowplot)
library(ggrepel)

group_cols <- setNames(c("#a0a0a0", "#696969", "#922b21", "#1f618d"), c("Background", "DEG", "Top Salmonella high", "Top control high"))
input <- do.call('rbind', input_list)
saveRDS(input, file="Res_Salmonella_infection_on_colon_goblet_cluster_DEG_all_res.rds")
p1 <- ggplot(input, aes(x=logFC, y=pval_transformed, size=auc_sym, fill=group))+
    geom_point(shape=21, color="#202020")+
    theme_bw()+
    theme(legend.position="bottom")+
    scale_fill_manual(values=group_cols, guide = guide_legend(override.aes = list(label = "")))+
    labs(x="logFC (Infection - Control)", size="|AUC-0.5|", y="-log10(P)^0.25")+
    #geom_vline(xintercept=0, color="black", linetype="dashed")+
    geom_text_repel(data=input[input$group=="Top Salmonella high",], aes(label=feature), color="#922b21", size=4, point.padding=1, segment.color='grey50', min.segment.length=0, box.padding=0.5, max.overlaps=999)+
    geom_text_repel(data=input[input$group=="Top control high",], aes(label=feature), color="#1f618d", size=4, point.padding=1, segment.color='grey50', min.segment.length=0, box.padding=0.5, max.overlaps=999)
    #scale_color_manual(values=cols$cell_type, guide = "none")
p1_split <- p1 + facet_wrap(~ Cluster, scales = "free")
#p <- plot_grid(plotlist = plot_list, align = "hv", nrow = 2, axis = "l")
ggsave(p1_split, filename="Plot_Salmonella_infection_on_colon_goblet_cluster_DEG_volcano_plot.png", width=10*2, height=10)


combined_sig_res <- deg_res_per_ct_combined %>% filter(padj<0.05 & pct_in>10 & logFC>0.1)
# compare the logFC in C5 and C9
ave_expr <- getAveExpr(seu.obj=que_obj, feature.to.calc="Cluster_per_condition", colname.prefix=NULL)
saveRDS(ave_expr, file="Res_salmonella_infection_per_goblet_cluster_expr.rds")
expr_diff <- ave_expr[,grep("Salmonella",colnames(ave_expr))] - ave_expr[,grep("Control",colnames(ave_expr))]
df <- data.frame(feature=rownames(expr_diff), expr_diff)
p1 <- ggplot(df, aes(x=Salmonella_5, y=Salmonella_9))+
    geom_point()
ggsave(p1, filename="Plot_logFC_comparison_between_goblet_clusters.pdf", width=5, height=5)

combined_up_res <- combined_sig_res[grep("Salmonella", combined_sig_res$group),]
up_gene_list <- tapply(combined_up_res$feature, combined_up_res$group, list)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_statiscal_test.R")
merged_list <- readRDS("../Res_KEGG_wikiPathway_hallMark_merged_annotated_gene_list.rds")
pathway_enrichment_pval <- enrichment_test(query_gene_list=up_gene_list, anno_gene_list=merged_list, all_genes=genes)
saveRDS(pathway_enrichment_pval, file="Res_salmonella_infection_per_goblet_cluster_DEG_KEGG_hallmark_wikiPathway_hypergeometric_test_nominal_pval.rds")

# perform p-value correction per gene set
idx_kegg <- grep("KEGG", rownames(pathway_enrichment_pval))
idx_wp <- grep("WP", rownames(pathway_enrichment_pval))
idx_hallmark <- grep("HALLMARK", rownames(pathway_enrichment_pval))
padj_list <- lapply(c("HALLMARK", "KEGG", "WP"), function(x){
    idx <- grep(x, rownames(pathway_enrichment_pval))
    pval <- pathway_enrichment_pval[idx,]
    apply(pval, 2, p.adjust, method="BH")
})
padj_mat <- do.call('rbind', padj_list)
saveRDS(padj_mat, file="Res_salmonella_infection_per_goblet_cluster_DEG_KEGG_hallmark_wikiPathway_hypergeometric_test_BH_padj.rds")

# Pull all pathways for human  
library(KEGGREST)
pathways.list <- keggList("pathway", "hsa")
head(pathways.list)
# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list)) 
genes.by.pathway <- lapply(pathway.codes,
    function(pwid){
        pw <- keggGet(pwid)
        if (is.null(pw[[1]]$GENE)) return(NA)
        pw2 <- pw[[1]]$GENE[c(FALSE,TRUE)] # may need to modify this to c(FALSE, TRUE) for other organisms
        pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
        pw2 <- pw2[!grepl("[", pw2, fixed = T)]
        return(pw2)
    }
)
pwid <- sapply(pathways.list, function(x){strsplit(x, split=" - ")[[1]][1]})
names(genes.by.pathway) <- pwid 
size <- sapply(genes.by.pathway, length)
kegg_gene_list <- genes.by.pathway[which(size>1)]
saveRDS(kegg_gene_list, file="Res_KEGG_gene_list.rds")
saveRDS(kegg_gene_list, file="/home/yuq22/Bacteria_TRM/Annotation/KEGG/Res_KEGG_gene_list.rds")
# merge the gene list with the same pathway name
merged_kegg_gene_list <- lapply(sort(unique(names(kegg_gene_list))), function(x){
    idx <- which(names(kegg_gene_list)==x)
    genes <- sort(unique(unlist(kegg_gene_list[idx])))
    return(genes)
})
names(merged_kegg_gene_list) <- sort(unique(names(kegg_gene_list)))
saveRDS(merged_kegg_gene_list, file="/home/yuq22/Bacteria_TRM/Annotation/KEGG/Res_KEGG_gene_list.rds")

size_2 <- sapply(kegg_gene_list, length)
terms_to_check <- c("Salmonella infection", "Bacterial invasion of epithelial cells", "Epithelial cell signaling in Helicobacter pylori infection", "Vibrio cholerae infection", "Pathogenic Escherichia coli infection")

kegg_gene_df <- data.frame('pathway'=rep(names(kegg_gene_list), size_2), 'gene'=unlist(kegg_gene_list), stringsAsFactors=F)
saveRDS(kegg_gene_df, file="/home/yuq22/Bacteria_TRM/Annotation/KEGG/Res_KEGG_gene_list_df.rds")

# check expression of genes associated with Salmonella infection / bacterial invasion of epithelial cells / epithelial cell signaling in Helicobacter pylori infection / Vibrio cholerae infection / Pathogenic Escherichia coli infection
# across colon tissue s2E clusters
tissue_s2e <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/RNA_per_condition/tissue_colon/Res_colon_s2e_adult_human_epi_with_CSS.rds")
# exclude TA cells from the tissue data
tissue_s2e <- subset(tissue_s2e, subset=RNA_CSS_snn_res.0.4!=0)
SCpubr::do_DimPlot(tissue_s2e, reduction="umap_css", group.by="RNA_CSS_snn_res.0.4", label=T, label.size=10)
tissue_expr <- getAveExpr(seu.obj=tissue_s2e, feature.to.calc="RNA_CSS_snn_res.0.4", colname.prefix="tissue")
tissue_s2e <- AddModuleScore(tissue_s2e, features=kegg_gene_list[terms_to_check], name="KEGG_gene_score", assay="RNA")
idx <- grep("KEGG_gene_score", colnames(tissue_s2e@meta.data))
colnames(tissue_s2e@meta.data)[idx] <- terms_to_check
tissue_s2e$RNA_CSS_snn_res.0.4 <- factor(tissue_s2e$RNA_CSS_snn_res.0.4, levels=c(8,2,3,4,5,6,7,1))
saveRDS(tissue_s2e, file="Res_adult_human_colon_epi_s2E_no_TA_with_CSS_and_KEGG_gene_score.rds")
p1 <- SCpubr::do_ViolinPlot(tissue_s2e, feature=terms_to_check, group.by="RNA_CSS_snn_res.0.4")&NoLegend()
ggsave(p1, filename="Plot_expression_of_genes_associated_with_Salmonella_infection_across_colon_tissue_s2E_clusters.pdf", width=5*3, height=5*2)
# across colon chip s2E clusters, control and infection condition separately and their logFC magnitude
s2e <- readRDS("Res_salmonella_data_s2e_with_DM.rds")
s2e$Cluster_per_condition <- paste(s2e$condition, s2e$RNA_snn_res.1, sep="_")
s2e <- AddModuleScore(s2e, features=kegg_gene_list[terms_to_check], name="KEGG_gene_score", assay="RNA")
idx <- grep("KEGG_gene_score", colnames(s2e@meta.data))
colnames(s2e@meta.data)[idx] <- terms_to_check
for(x in unique(s2e$condition)){
    seu_obj <- subset(s2e, subset=condition==x)
    p2 <- SCpubr::do_ViolinPlot(seu_obj, feature=terms_to_check, group.by="RNA_snn_res.1", ncol=1)&NoLegend()
    ggsave(p2, filename=paste0("Plot_expression_of_genes_associated_with_Salmonella_infection_across_",x,"_colon_chip_s2E_clusters.pdf"), width=7, height=5*length(terms_to_check))
}


combined_expr_mat <- readRDS("Dat_combined_cluster_average_expression_per_condition.rds")
# identify the DEGs with differential change magnitude along the trajectory
diff_expr_mat <- combined_expr_mat[,grep("Salmonella", colnames(combined_expr_mat))] - combined_expr_mat[,grep("Control", colnames(combined_expr_mat))]

selected_gs_names <- sig_terms
score_mat <- t(sapply(selected_gs_names, function(x){
    term_degs <- intersect(merged_df$gene_symbol[which(merged_df$gs_name==x)], union_up_genes)
    #apply(expr_diff_mat[term_degs,],2,mean)
    apply(abs(expr_diff_mat[term_degs,]),2,mean)
}))

# show BH-adjusted p-value of selected terms
#values <- -log10(padj_mat[selected_gs_names, "all"])
values <- -log10(apply(padj_mat[selected_gs_names, ], 1, min))
gs_sorted <- names(values)[order(values)]

df <- unique(s2e@meta.data[,c("RNA_snn_res.1", "Subtype")])
df$new_cl_idx <- rank(df$RNA_snn_res.1)
df_2 <- df[order(df$RNA_snn_res.1),]
df_2$id_2 <- paste("Salmonella", df_2$new_cl_idx, sep="_")
col_sidebar_id <- df_2$Subtype
cols <- readRDS("~/Bacteria_TRM/used_object/cols.rds")
# get the median pt value per cluster
median_pt_per_cluster <- sapply(sort(unique(s2e$RNA_snn_res.1)), function(i){
    median(s2e$dc1_rank[which(s2e$RNA_snn_res.1==i)])
})
colorPal <- grDevices::colorRampPalette(c("#f7f7f7", "#cccccc", "#969696", "#636363", "#252525"))
cellColor <- adjustcolor(colorPal(30), alpha = 0.8)[as.numeric(cut(median_pt_per_cluster, 
            breaks = 30, right = F, include.lowest = T))]
#pdf("Plot_selected_per_cluster_enriched_KEGG_pathway_DEG_logFC_across_clusters.pdf", height=7, width=8)
#gplots::heatmap.2(score_mat[gs_sorted,], scale="row", trace="none", dendrogram="none", col=darkBlue2Red.heatmap, ColSideColors=cellColor, density.info="none", Colv=FALSE, Rowv=FALSE, keysize=0.8, mar=c(5,30))
#dev.off()

pdf("Plot_selected_per_cluster_enriched_KEGG_pathway_DEG_logFC_across_clusters_2.pdf", height=7, width=8)
gplots::heatmap.2(score_mat[gs_sorted,], scale="row", trace="none", dendrogram="row", col=darkBlue2Red.heatmap, ColSideColors=cols$cell_type[df_2$Subtype], density.info="none", Colv=FALSE, Rowv=FALSE, keysize=0.8, mar=c(5,30))
dev.off()

# try to summarize the results for figure
setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/for_figure")
# how many DEGs are detetected in each cluster of stem cells, colonocytes and goblet cells
## for goblet cells
goblet_deg_res_per_ct_combined <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/goblet/Res_Salmonella_infection_on_colon_goblet_cluster_DEG.rds")
goblet_sig_res <- goblet_deg_res_per_ct_combined %>% filter(padj<0.05 & pct_in>10 & logFC>0.1)
goblet_up_genes <- unique(goblet_sig_res$feature[grep("Salmonella", goblet_sig_res$group)])
## for stem cells and colonocytes
s2e_sig_res <- readRDS("../Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_matched_control_cells_DEG_sig_res.rds")
stem_up_genes <- unique(s2e_sig_res$feature[s2e_sig_res$group=="C13_Salmonella"])
colonocyte_up_genes <- unique(s2e_sig_res$feature[grepl("Salmonella", s2e_sig_res$group) & s2e_sig_res$group!="C13_Salmonella"])
sig_res <- list("goblet"=goblet_sig_res, "s2e"=s2e_sig_res)
saveRDS(sig_res, file="Res_Salmonella_infection_on_colon_epi_per_cell_type_sig_DE_test_res.rds")
deg_list <- list("goblet"=goblet_up_genes, "stem"=stem_up_genes, "colonocyte"=colonocyte_up_genes)
saveRDS(deg_list, file="Res_Salmonella_infection_on_colon_epi_per_cell_type_Salmonella_up_genes.rds")

# focusing on up regulated genes
union_up_genes <- unique(unlist(deg_list))
idx <- matrix(F, nrow=length(union_up_genes), ncol=3)
rownames(idx) <- union_up_genes
colnames(idx) <- names(deg_list)
for(x in names(deg_list)){
    idx[deg_list[[x]], x] <- T
}
saveRDS(idx, file="Res_Salmonella_infection_on_colon_epi_per_cell_type_Salmonella_up_genes_sig_idx.rds")

goblet_only <- rownames(idx)[rowSums(idx)==1 & idx[,"goblet"]==1]
stem_only <- rownames(idx)[rowSums(idx)==1 & idx[,"stem"]==1]
colonocyte_only <- rownames(idx)[rowSums(idx)==1 & idx[,"colonocyte"]==1]
goblet_stem_coloncyte <- rownames(idx)[rowSums(idx)==3]
goblet_stem <- rownames(idx)[rowSums(idx)==2 & idx[,"goblet"]==1 & idx[,"stem"]==1]
goblet_colonocyte <- rownames(idx)[rowSums(idx)==2 & idx[,"goblet"]==1 & idx[,"colonocyte"]==1]
stem_colonocyte <- rownames(idx)[rowSums(idx)==2 & idx[,"stem"]==1 & idx[,"colonocyte"]==1]
counts <- setNames(c(length(goblet_only), length(stem_only), length(colonocyte_only), length(goblet_colonocyte), length(stem_colonocyte), length(goblet_stem), length(goblet_stem_coloncyte)), 
                   c("goblet_only", "stem_only", "colonocyte_only", "goblet_colonocyte", "stem_colonocyte", "goblet_stem", "goblet_stem_coloncyte"))
ps <- counts/sum(counts)                   
saveRDS(ps, file="Res_Salmonella_up_genes_overlaps.rds")
vd <-  venneuler::venneuler(c(G=as.numeric(ps["goblet_only"]), S = as.numeric(ps["stem_only"]), C = as.numeric(ps["colonocyte_only"]),
                  "G&S" = as.numeric(ps["goblet_stem"]),
                  "G&C" = as.numeric(ps["goblet_colonocyte"]),
                  "S&C" = as.numeric(ps["stem_colonocyte"]),
                  "G&S&C" = as.numeric(ps["goblet_stem_coloncyte"])))
pdf("Plot_DEG_per_cell_type_overlap.pdf", height=5, width=5)
plot(vd)
dev.off()



s2e <- readRDS("Res_salmonella_data_s2e_with_DM.rds")
SCpubr::do_ViolinPlot(s2e, feature="CDC42", group.by="RNA_snn_res.1", ncol=1)
tissue_s2e <- readRDS("Res_adult_human_colon_epi_s2E_no_TA_with_CSS_and_KEGG_gene_score.rds")
invasion_genes <- c("CDC42", "RAC1", "RHOA", "ARF1", "ARF6")
defense_genes <- c("CDHR2", "CDHR5", "FUT2", "FUT3", "MUC1", "MUC5B", "MUC3A", "MUC3A", "MUC13")
p1 <- SCpubr::do_ViolinPlot(tissue_s2e, feature=defense_genes, group.by="RNA_CSS_snn_res.0.4", ncol=3)&NoLegend()
ggsave(p1, filename="Plot_colon_tissue_s2E_defense_related_genes.pdf", width=5*3, height=5*3)

# cytosekeleton and junction related genes
kegg_gene_list <- readRDS("/home/yuq22/Bacteria_TRM/Annotation/KEGG/Res_KEGG_gene_list.rds")
terms_to_check <- c("Regulation of actin cytoskeleton", "Tight junction", "Gap junction")
tissue_s2e <- AddModuleScore(tissue_s2e, features=kegg_gene_list[terms_to_check], name="KEGG_gene_score", assay="RNA")
idx <- grep("KEGG_gene_score", colnames(tissue_s2e@meta.data))
colnames(tissue_s2e@meta.data)[idx] <- terms_to_check
p1 <- SCpubr::do_ViolinPlot(tissue_s2e, feature=terms_to_check, group.by="RNA_CSS_snn_res.0.4")&NoLegend()
ggsave(p1, filename="Plot_colon_tissue_s2E_actin_and_junction_related_genes.pdf", width=5*3, height=5)


s2e$Cluster_per_condition <- paste(s2e$condition, s2e$RNA_snn_res.1, sep="_")
saveRDS(s2e, file="Res_salmonella_data_s2e_with_DM.rds")
s2e <- AddModuleScore(s2e, features=kegg_gene_list[terms_to_check], name="KEGG_gene_score", assay="RNA")
idx <- grep("KEGG_gene_score", colnames(s2e@meta.data))
colnames(s2e@meta.data)[idx] <- terms_to_check
p2 <- SCpubr::do_ViolinPlot(s2e, feature=terms_to_check, group.by="Subtype_per_condition", ncol=1)&NoLegend()
ggsave(p2, filename="Plot_expression_of_genes_associated_with_cytoskeleton_and_junction_across_colon_chip_s2E_subtypes.pdf", width=10, height=5*length(terms_to_check))


for(x in unique(s2e$condition)){
    seu_obj <- subset(s2e, subset=condition==x)
    p2 <- SCpubr::do_ViolinPlot(seu_obj, feature=terms_to_check, group.by="Subtype", ncol=1)&NoLegend()
    ggsave(p2, filename=paste0("Plot_expression_of_genes_associated_with_cytoskeleton_and_junction_across_",x,"_colon_chip_s2E_subtypes.pdf"), width=7, height=5*length(terms_to_check))
}

cols <- readRDS("/home/yuq22/Bacteria_TRM/used_object/cols.rds")

# compare the condition logFC between cell types
s2e_ave_exp <- readRDS("../Dat_combined_cluster_average_expression_per_condition.rds")
goblet_ave_expr <- readRDS("../goblet/Res_salmonella_infection_per_goblet_cluster_expr.rds")

s2e_logFC <- s2e_ave_exp[,grep("Salmonella", colnames(s2e_ave_exp))] - s2e_ave_exp[,grep("Control", colnames(s2e_ave_exp))]
stem_cell_logFC <- s2e_logFC[,"Salmonella_1"]
colonocyte_logFC <- apply(s2e_logFC[,-1], 1, function(vec){
    idx <- which.max(abs(vec))
    vec[idx]
})
goblet_cl_logFC <- goblet_ave_expr[,grep("Salmonella", colnames(goblet_ave_expr))] - goblet_ave_expr[,grep("Control", colnames(goblet_ave_expr))]
goblet_logFC <- apply(goblet_cl_logFC, 1, function(vec){
    idx <- which.max(abs(vec))
    vec[idx]
})
length(stem_cell_logFC)
ct_logFC <- data.frame('feature'=rownames(s2e_logFC), 'stem_cell'=stem_cell_logFC, 'colonocyte'=colonocyte_logFC, 'goblet'=goblet_logFC)
gene <- ct_logFC$feature[!grepl("^ENSG", ct_logFC$feature)]
#remove ribosomal genes, mitochondrial genes, and genes located on sex chromosomes before normalization
confound_genes <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/confound_genes/Dat_confound_genes.rds")
genes_to_exclude <- unique(unlist(confound_genes[c(4,5,6)]))
s2e_expressed_genes <- rownames(s2e_ave_exp)[which(rowSums(s2e_ave_exp)>0)]
goblet_expressed_genes <- rownames(goblet_ave_expr)[which(rowSums(goblet_ave_expr)>0)]
union_expressed_genes <- union(s2e_expressed_genes, goblet_expressed_genes)
genes <- intersect(setdiff(gene, genes_to_exclude), union_expressed_genes)
ct_logFC <- ct_logFC[ct_logFC$feature%in%genes,]
saveRDS(ct_logFC, file="Dat_max_logFC_per_s2e_and_goblet_cell_type.rds")


highlight_genes <- c("LTB", "CXCL8", "IL1B", "CCL20", "LCN2")
idx <- readRDS("Res_Salmonella_infection_on_colon_epi_per_cell_type_Salmonella_up_genes_sig_idx.rds")
goblet_only <- rownames(idx)[rowSums(idx)==1 & idx[,"goblet"]==1]
stem_only <- rownames(idx)[rowSums(idx)==1 & idx[,"stem"]==1]
colonocyte_only <- rownames(idx)[rowSums(idx)==1 & idx[,"colonocyte"]==1]
goblet_stem_coloncyte <- rownames(idx)[rowSums(idx)==3]
goblet_stem <- rownames(idx)[rowSums(idx)==2 & idx[,"goblet"]==1 & idx[,"stem"]==1]
goblet_colonocyte <- rownames(idx)[rowSums(idx)==2 & idx[,"goblet"]==1 & idx[,"colonocyte"]==1]
stem_colonocyte <- rownames(idx)[rowSums(idx)==2 & idx[,"stem"]==1 & idx[,"colonocyte"]==1]

# get cell type specific DEGs with top10 logFC
gene_list <- list("goblet"=goblet_only, "stem_cell"=stem_only, "colonocyte"=colonocyte_only, "goblet&stem_cell"=goblet_stem, "goblet&colonocyte"=goblet_colonocyte, "stem_cell&colonocyte"=stem_colonocyte, "goblet&stem_cell&colonocyte"=goblet_stem_coloncyte)
gene_df <- data.frame('gene_symbol'=unlist(gene_list), 'group'=rep(names(gene_list), sapply(gene_list, length)), stringsAsFactors=F)
saveRDS(gene_df, file="Res_DEG_per_cell_type.rds")
saveRDS(gene_list, file="Res_DEG_per_cell_type_gene_list.rds")
top_gene_list <- list()
for(x in names(gene_list)){
    g <- gene_list[[x]]
    ct <- strsplit(x, split="&")[[1]]
    if(length(ct)==1){
        vec <- ct_logFC[g, ct]
    }else{
        vec <- apply(ct_logFC[g, ct], 1, function(fc){
            idx <- which.max(abs(fc))
            fc[idx]
        })
    }
    n <- min(c(10, length(g)))
    top_gene_list[[x]] <- g[order(abs(vec), decreasing=T)[1:n]]
}
saveRDS(top_gene_list, file="Res_DEG_per_cell_type_top_gene_list.rds")

ct_logFC$top_genes <- NA
for(x in names(top_gene_list)){
    ct_logFC$top_genes[ct_logFC$feature%in%top_gene_list[[x]]] <- x
}
saveRDS(ct_logFC, file="Dat_max_logFC_per_s2e_and_goblet_cell_type.rds")

setwd("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/for_figure")
ct_logFC <- readRDS("Dat_max_logFC_per_s2e_and_goblet_cell_type.rds")
data <- as.matrix(ct_logFC[,c("goblet", "stem_cell", "colonocyte")])
ct_pairs <- t(combn(c("goblet","stem_cell","colonocyte"),2))

par(mfrow=c(2,2))
for(i in seq(nrow(ct_pairs))){
    ct1 <- ct_pairs[i,1]
    ct2 <- ct_pairs[i,2]
    dat1 <- data[,c(ct1, ct2)]
    plot(dat1)
    abline(h=0, lty=2)
    abline(v=0, lty=2)
}

# get the genes upregulated in goblet cells but downregulated in colonocytes
both_up <- rownames(data)[which(data[,"goblet"]>0.1 & data[,"colonocyte"]>0.1)]
both_down <- rownames(data)[which(data[,"goblet"]< -0.1 & data[,"colonocyte"]< -0.1)]
goblet_up_col_down <- rownames(data)[which(data[,"goblet"]>0.1 & data[,"colonocyte"]< -0.1)]
col_up_goblet_down <- rownames(data)[which(data[,"goblet"]< -0.1 & data[,"colonocyte"]>0.1)]
gene_list <- list("both_up"=both_up, "both_down"=both_down, "goblet_up_col_down"=goblet_up_col_down, "col_up_goblet_down"=col_up_goblet_down)

# kegg pathway enrichment analysis for DEGs of each cell type and take the union of enriched pathways
#deg_idx <- readRDS("Res_Salmonella_infection_on_colon_epi_per_cell_type_Salmonella_up_genes_sig_idx.rds")
# perform enrichment analysis on all DEGs
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_statiscal_test.R")
background_genes <- rownames(data)
kegg_gene_list <- readRDS("/home/yuq22/Bacteria_TRM/Annotation/KEGG/Res_KEGG_gene_list.rds")
pathway_enrichment_pval <- enrichment_test(query_gene_list=gene_list, anno_gene_list=kegg_gene_list, all_genes=background_genes)
saveRDS(pathway_enrichment_pval, file="Res_salmonella_infection_goblet_and_colonocyte_DEG_KEGG_hypergeometric_test_nominal_pval.rds")

pathway_enrichment_padj <- apply(pathway_enrichment_pval, 2, p.adjust, method="BH")
saveRDS(pathway_enrichment_padj, file="Res_salmonella_infection_goblet_and_colonocyte_DEG_KEGG_hypergeometric_test_BH-corrected_pval.rds")

pathway_enrichment_padj <- readRDS("Res_salmonella_infection_goblet_and_colonocyte_DEG_KEGG_hypergeometric_test_BH-corrected_pval.rds")
terms_to_check <- c("Salmonella infection", "Bacterial invasion of epithelial cells", "Regulation of actin cytoskeleton")
padj_data <- -log10(pathway_enrichment_padj)
padj_diff <- padj_data - rowMeans(padj_data)
# exclude the cancer and disease related terms
input <- padj_diff[!grepl("disease|cancer|carcino|muscle", rownames(padj_diff)),]
top_terms <- lapply(colnames(padj_diff), function(x){
    intersect(rownames(input)[order(input[,x], decreasing=T)[1:10]], rownames(pathway_enrichment_padj)[pathway_enrichment_padj[,x]<0.1])
})
names(top_terms) <- colnames(padj_diff)
saveRDS(top_terms, file="Res_goblet_and_colonocyte_DEG_KEGG_hypergeometric_test_top_enriched_terms.rds")

terms <- unique(unlist(top_terms))
res <- prepareTreeAndHeatmapInput(expr.mat = padj_data[terms,], hc.method = "ward.D2", norm.method = "quantile")
saveRDS(res, file="Res_goblet_and_colonocyte_DEG_KEGG_hypergeometric_test_heatmap_input.rds")

black_cols <- c("#f7f7f7", "#cccccc", "#969696", "#636363", "#252525")
black_cols_heatmap <- colorRampPalette(black_cols)(30) 
pdf("Plot_heatmap_kegg_pathway_enrichment.pdf", height=10)
gplots::heatmap.2(res$heatmap_input, col = black_cols_heatmap, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 1, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

num_mat_list <- list()
plot_list <- list()
for(i in seq(nrow(ct_pairs))){
    ct1 <- ct_pairs[i,1]
    ct2 <- ct_pairs[i,2]
    dat1 <- data[,c(ct1, ct2)]
    ct1_pos_genes <- rownames(data)[which(data[,ct1]>0)]
    ct2_pos_genes <- rownames(data)[which(data[,ct2]>0)]
    ct1_neg_genes <- rownames(data)[which(data[,ct1]<0)]
    ct2_neg_genes <- rownames(data)[which(data[,ct2]<0)]
    ct1_pos_ct2_pos_genes <- intersect(ct1_pos_genes, ct2_pos_genes)
    ct1_pos_ct2_neg_genes <- intersect(ct1_pos_genes, ct2_neg_genes)
    ct1_neg_ct2_pos_genes <- intersect(ct1_neg_genes, ct2_pos_genes)
    ct1_neg_ct2_neg_genes <- intersect(ct1_neg_genes, ct2_neg_genes)
    n1 <- length(ct1_pos_ct2_pos_genes)
    n2 <- length(ct1_pos_ct2_neg_genes)
    n3 <- length(ct1_neg_ct2_pos_genes)
    n4 <- length(ct1_neg_ct2_neg_genes)
    n_total <- n1+n2+n3+n4
    num_mat <- matrix(c(n1/n_total, n2/n_total, n3/n_total, n4/n_total), nrow=2, ncol=2)
    rownames(num_mat) <- c(paste0(ct2, "_pos"), paste0(ct2, "_neg"))
    colnames(num_mat) <- c(paste0(ct1, "_pos"), paste0(ct1, "_neg"))
    num_mat_list[[i]] <- num_mat

    #plot_list[[i]] <- ggplot(dat1, aes_string(x=ct1, y=ct2)) +
    #stat_density_2d(geom = "polygon", contour = TRUE,
    #              aes(fill = after_stat(level)), colour = "black",
    #              bins = 10)+
    #geom_hline(yintercept=0, color="black", linetype="dashed")+
    #geom_vline(xintercept=0, color="black", linetype="dashed")+
    #scale_fill_distiller(palette = "Greens", direction = 1) +
    #theme_bw() +
    #labs(x=paste0(ct1, " logFC"), y=paste0(ct2, " logFC"))            
}
p <- plot_grid(plotlist = plot_list, align = "hv",nrow = 1, axis = "l")
p
ggsave(p, filename="Plot_density_plot_logFC_comparison_between_cell_types.png", width=6*3, height=5)

for(i in seq(nrow(ct_pairs))){
    ct1 <- ct_pairs[i,1]
    ct2 <- ct_pairs[i,2]
    highlight_genes <- c(top_gene_list[[ct1]], top_gene_list[[ct2]])
    input <- ct_logFC
    input$group <- "All"
    input$group[input$feature%in%highlight_genes] <- input$top_genes[input$feature%in%highlight_genes]
    p1 <- ggplot(input, aes_string(x=ct1, y=ct2))+
        geom_point(shape=21, color="#303030", aes(fill=group))+
        scale_fill_manual(values=c("All"="grey", "goblet"="#4DA7B0", "stem_cell"="#A84798", "colonocyte"="#FDD884"), guide = guide_legend(override.aes = list(label = "")))+
        geom_abline(slope=1, intercept=0, linetype="dashed")+
        geom_hline(yintercept=0, color="black", linetype="dashed")+
        geom_vline(xintercept=0, color="black", linetype="dashed")+
        theme_bw()+
        labs(x=paste0(ct1, " logFC"), y=paste0(ct2, " logFC"))+
        ggrepel::geom_text_repel(data=input[input$feature%in%highlight_genes,], aes(label=feature, color=group), size=4, point.padding=1, segment.color='grey50', min.segment.length=0, box.padding=0.5, max.overlaps=999)+
        scale_color_manual(values=c("All"="grey", "goblet"="#4DA7B0", "stem_cell"="#A84798", "colonocyte"="#FDD884"), guide = "none")
    
    plot_list[[i]] <- p1
}
p <- plot_grid(plotlist = plot_list, align = "hv",nrow = 1, axis = "l")
ggsave(p, filename="Plot_logFC_comparison_between_cell_types.png", width=6*3, height=5)

#example_genes <- c("LTB", "CXCL8", "IL1B")

# load the ligand receptor pairs
load("/home/yuq22/ihb-intestine-evo/Annotation/CellChat_2021/CellChatDB.human.rda")
complex_members <- setNames(apply(CellChatDB.human$complex, 1, function(vec){
    paste(vec[which(vec!="")], collapse=",")
}), rownames(CellChatDB.human$complex))
# include gene members of the complex
## for ligands
CellChatDB.human$interaction$ligand_genes <- CellChatDB.human$interaction$ligand
id <- intersect(CellChatDB.human$interaction$ligand_genes, names(complex_members))
idx <- which(CellChatDB.human$interaction$ligand_genes%in%id)
CellChatDB.human$interaction$ligand_genes[idx] <- complex_members[CellChatDB.human$interaction$ligand_genes[idx]]
## for receptors
CellChatDB.human$interaction$receptor_genes <- CellChatDB.human$interaction$receptor
id <- intersect(CellChatDB.human$interaction$receptor_genes, names(complex_members))
idx <- which(CellChatDB.human$interaction$receptor_genes%in%id)
CellChatDB.human$interaction$receptor_genes[idx] <- complex_members[CellChatDB.human$interaction$receptor_genes[idx]]
saveRDS(CellChatDB.human, file="/home/yuq22/ihb-intestine-evo/Annotation/CellChat_2021/CellChatDB_human.rds")

interaction_df <- CellChatDB.human$interaction
receptor_gene_vec <- unlist(strsplit(CellChatDB.human$interaction$receptor_genes, split=","))
all_receptor_genes <- unique(receptor_gene_vec)
receptor_num <- sapply(interaction_df$receptor_genes, function(x)length(unlist(strsplit(x, split=","))))
df <- data.frame('ligand'=rep(interaction_df$ligand_genes, receptor_num), 'receptor'=receptor_gene_vec, stringsAsFactors=F)

ligand_gene_vec <- unlist(strsplit(df$ligand, split=","))
all_receptor_genes <- unique(ligand_gene_vec)
ligand_num <- sapply(df$ligand, function(x)length(unlist(strsplit(x, split=","))))
df_2 <- data.frame('ligand'=ligand_gene_vec, 'receptor'=rep(df$receptor, ligand_num), stringsAsFactors=F)
df_3 <- unique(df_2)
saveRDS(df_3, file="/home/yuq22/ihb-intestine-evo/Annotation/CellChat_2021/Res_ligand_receptor_pairs_df.rds")
# get the pairs whose ligands are involved in "NF-kappa B signaling pathway" and expressed in colonocytes/goblet cells and whole paired receptors are expressed in our TRM data (receptors no need to be involved in the pathway)
kegg_gene_df <- readRDS("/home/yuq22/Bacteria_TRM/Annotation/KEGG/Res_KEGG_gene_list_df.rds")
gs_name <- "NF-kappa B signaling pathway"
pathway_genes <- kegg_gene_df$gene[kegg_gene_df$pathway==gs_name]
pathway_ligand_genes <- intersect(pathway_genes, all_ligand_genes)
paired_receptor_genes <- unique(df_3$receptor[df_3$ligand%in%pathway_ligand_genes])
# load TRM data
trm <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/TRM/TRM_from_coculture/Dat_merged_filtered_TRM_from_coculture_no_epi.rds")
trm_cluster_expr <- getAveExpr(seu.obj=trm, feature.to.calc="RNA_snn_res.1", colname.prefix="TRM")
trm_expressed_genes <- rownames(trm_cluster_expr)[which(apply(trm_cluster_expr, 1, max)>0.25)]
p1 <- SCpubr::do_DotPlot(trm, feature=intersect(paired_receptor_genes, trm_expressed_genes), group.by="RNA_snn_res.1", cluster=T)+labs('title'="TRM expressed receptor genes \n paired with NF-kappa B signaling pathway ligands")+theme(plot.title = element_text(size=15))
# load salmoneall infection data
combined <- readRDS("../Dat_salmonella_colon_epithelium_scRNA-seq_data_with_subtypes.rds")
combined$Cluster_per_condition <- paste(combined$RNA_snn_res.1, combined$condition, sep="_")
saveRDS(combined, file="../Dat_salmonella_colon_epithelium_scRNA-seq_data_with_subtypes.rds")
combined_cluster_expr <- getAveExpr(seu.obj=combined, feature.to.calc="Subtype_per_condition", colname.prefix="Epi")
combined_cluster_prop <- getExpressedProp(seu.obj=combined, feature.to.calc="Subtype_per_condition", colname.prefix="Epi")
max_expr <- apply(combined_cluster_expr, 1, max)
max_prop <- apply(combined_cluster_prop, 1, max)
epi_expressed_genes <- rownames(combined_cluster_expr)[which(apply(combined_cluster_prop, 1, max)>0.03)]
p2 <- SCpubr::do_DotPlot(combined, feature=intersect(pathway_ligand_genes, epi_expressed_genes), group.by="Cluster_per_condition")+labs('title'="NF-kappa B signaling pathway ligands in colon epithelium")
p <- p2+p1
ggsave(p, filename="Plot_dotplot_NF-kappa_B_signaling_pathway_ligands_in_colon_epithelium_and_TRM_expressed_receptor_genes.pdf", width=15, height=10)

# for HLA gene family
hla_genes <- grep("HLA", rownames(combined), value=T)
hla_receptor_genes <- unique(df_3$receptor[df_3$ligand%in%hla_genes])
p2_hla <- SCpubr::do_DotPlot(combined, feature=intersect(hla_genes, epi_expressed_genes), group.by="Cluster_per_condition")+labs('title'="HLA genes in colon epithelium")
p1_hla <- SCpubr::do_DotPlot(trm, feature=intersect(hla_receptor_genes, trm_expressed_genes), group.by="RNA_snn_res.1", cluster=T)+labs('title'="TRM expressed HLA receptor")+theme(plot.title = element_text(size=15))
p_hla <- p2_hla+p1_hla
ggsave(p_hla, filename="Plot_dotplot_HLA_ligands_in_colon_epithelium_and_TRM_expressed_receptor_genes.pdf", width=15, height=10)

#  
gs_name <- "Antigen processing and presentation"
pathway_genes <- kegg_gene_df$gene[kegg_gene_df$pathway==gs_name]
pathway_ligand_genes <- intersect(pathway_genes, all_ligand_genes)
paired_receptor_genes <- unique(df_3$receptor[df_3$ligand%in%pathway_ligand_genes])
p2 <- SCpubr::do_DotPlot(combined, feature=intersect(pathway_ligand_genes, epi_expressed_genes), group.by="Cluster_per_condition")+labs('title'= paste(gs_name, "ligands in colon epithelium"))
p1 <- SCpubr::do_DotPlot(trm, feature=intersect(paired_receptor_genes, trm_expressed_genes), group.by="RNA_snn_res.1", cluster=T)+labs('title'=paste("TRM expressed receptor genes \n paired with", gs_name, "ligands"))+theme(plot.title = element_text(size=15))
p <- p2+p1
p
ggsave(p, filename=paste0("Plot_dotplot_", gsub(" ", "_", gs_name), "_ligands_in_colon_epithelium_and_TRM_expressed_receptor_genes.pdf"), width=15, height=10)


ct_logFC <- readRDS("Dat_max_logFC_per_s2e_and_goblet_cell_type.rds")
data <- as.matrix(ct_logFC[,c("goblet", "stem_cell", "colonocyte")])
total_up <- colSums(data>0.1)
total_down <- colSums(data< -0.1)


sig_res <- readRDS("Res_Salmonella_infection_on_colon_epi_per_cell_type_sig_DE_test_res.rds")
goblet_sig_res <- sig_res$goblet
goblet_up_genes <- unique(goblet_sig_res$feature[grep("Salmonella", goblet_sig_res$group)])
goblet_down_genes <- unique(goblet_sig_res$feature[grep("Control", goblet_sig_res$group)])
## for stem cells and colonocytes
s2e_sig_res <- sig_res$s2e
stem_up_genes <- unique(s2e_sig_res$feature[s2e_sig_res$group=="C13_Salmonella"])
colonocyte_up_genes <- unique(s2e_sig_res$feature[grepl("Salmonella", s2e_sig_res$group) & s2e_sig_res$group!="C13_Salmonella"])
stem_down_genes <- unique(s2e_sig_res$feature[s2e_sig_res$group=="C13_Control"])
colonocyte_down_genes <- unique(s2e_sig_res$feature[grepl("Control", s2e_sig_res$group) & s2e_sig_res$group!="C13_Control"])
sig_up_genes <- list("goblet"=goblet_up_genes, "stem"=stem_up_genes, "colonocyte"=colonocyte_up_genes)
sig_down_genes <- list("goblet"=goblet_down_genes, "stem"=stem_down_genes, "colonocyte"=colonocyte_down_genes)
sig_up <- sapply(sig_up_genes, length)
sig_down <- sapply(sig_down_genes, length)
df <- data.frame('total_up'=total_up, 'total_down'=total_down, 'sig_up'=sig_up, 'sig_down'=sig_down)
saveRDS(df, file="Res_per_cell_type_DEG_number.rds")

input <- df[,c("total_up", "total_down")]
input$total_down <- -input$total_down
input <- as.matrix(input)
dat <- data.frame(
    group = rep(c("Salmonella", "Control"), each=nrow(input)),
    x = rep(rownames(input), 2),
    y = as.vector(input)
)
ggplot(dat, aes(x=x, y=y, fill=group)) + 
  geom_bar(stat="identity", position="identity")

# ggplot barplot to show the number of total_up and total_down
# focusing on up regulated genes
union_up_genes <- unique(unlist(deg_list))
idx <- matrix(F, nrow=length(union_up_genes), ncol=3)
rownames(idx) <- union_up_genes
colnames(idx) <- names(deg_list)
for(x in names(deg_list)){
    idx[deg_list[[x]], x] <- T
}
saveRDS(idx, file="Res_Salmonella_infection_on_colon_epi_per_cell_type_Salmonella_up_genes_sig_idx.rds")
