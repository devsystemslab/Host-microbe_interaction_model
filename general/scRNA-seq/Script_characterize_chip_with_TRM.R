setwd("/home/yuq22/Bacteria_TRM/TRM")
library(Seurat)
library(ggplot2)
library(dplyr)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")

# Load the data
folder <- "/projects/site/pred/ihb-g-deco/USERS/lopezsar/Bacteria_TRM_Manuscript/3rd Figure - Epithelium vs Epithelium+TRM"
ids <- list.files(path=folder, pattern="Epithelium")
seu_obj_list <- list()
for(x in ids){
    count <- Read10X(data.dir=paste0(folder, "/", x))
    colnames(count) <- paste(x, colnames(count), sep="_")
    seu_obj <- CreateSeuratObject(count)
    seu_obj_list[[x]] <- seu_obj
}
combined_rna <- merge(x=seu_obj_list[[1]], y=seu_obj_list[-1])
combined_rna <- NormalizeData(object = combined_rna, normalization.method = "LogNormalize", scale.factor = 1e4)
# calculate Mt%
combined_rna[["percent.mt"]] <- PercentageFeatureSet(combined_rna, pattern = "^MT-")
saveRDS(combined_rna, file="Dat_merged_ileum_epithelium_TRM_coculture.rds")
combined_rna <- readRDS("Dat_merged_ileum_epithelium_TRM_coculture.rds")
combined_rna$orig.ident[grepl("TRM", colnames(combined_rna))] <- "withTRM"
combined_rna$orig.ident[!grepl("TRM", colnames(combined_rna))] <- "noTRM"
Idents(combined_rna) <- combined_rna$orig.ident
combined_rna$log_nCount <- log(combined_rna$nCount_RNA)
saveRDS(combined_rna, file="Dat_merged_ileum_epithelium_TRM_coculture.rds")
# show Mt%, nFeature, nCount distribution per sample
# Visualize QC metrics as a violin plot
p1 <- SCpubr::do_ViolinPlot(combined_rna, features = c("nFeature_RNA", "nCount_RNA", "log_nCount","percent.mt"), ncol = 4)
pdf("Plot_ViolinPlot_QC_metrics_per_sample.pdf" , width = 20, height = 5)
p1
dev.off()
x <- exp(8)
p1 <- SCpubr::do_ViolinPlot(combined_rna, features = "nCount_RNA")
p1+geom_hline(yintercept = x, linetype = "dashed", color = "red")
# filter cells
combined_rna <- readRDS("Dat_merged_ileum_epithelium_TRM_coculture.rds")
cells <- colnames(combined_rna)[which(combined_rna$nFeature_RNA > 1000 & combined_rna$nCount_RNA > 3000 & combined_rna$nCount_RNA<50000 & combined_rna$percent.mt < 25)]
seu_obj <- subset(combined_rna, cells=cells)
#p2 <- SCpubr::do_ViolinPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "log_nCount","percent.mt"), ncol = 4)
#pdf("Plot_ViolinPlot_QC_metrics_per_sample_after_filtering.pdf" , width = 20, height = 5)
#p2
#dev.off()
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 3000)
seu_obj <- ScaleData(object = seu_obj, verbose = T)
seu_obj <- RunPCA(object = seu_obj, features = VariableFeatures(seu_obj), verbose = F, npcs = 50)
usefulPCs <- 1:20
seu_obj <- FindNeighbors(object = seu_obj, dims = usefulPCs)
seu_obj <- FindClusters(object = seu_obj, resolution = 1)
seu_obj <- RunUMAP(object = seu_obj, dims = usefulPCs)
saveRDS(seu_obj, file="Dat_merged_filtered_ileum_epithelium_TRM_coculture.rds")

# examine the regional identity the epithelial cells
p1 <- SCpubr::do_DimPlot(seu_obj, group.by = "orig.ident", pt.size=6, font.size=100, border.size=1.5, legend.icon.size=40)
p2 <- SCpubr::do_DimPlot(seu_obj, group.by = "RNA_snn_res.1", pt.size=6, label=T, label.size=20, border.size=1.5)+NoLegend()
png("Plot_UMAP_colored_by_condition_and_cluster.png", width = 2000, height = 2000*2)
p1/p2
dev.off()

plotFeature(seu.obj=seu_obj, dr="umap", col.num=5, genes.to.plot=c("EPCAM", "PTPRC", "PDX1", "ONECUT2", "CDX2", "OSR2", "HOXB7", "HOXA10","SATB2","FABP6", "APOA4","BEST4","RBPJ","LGR5","","MUC2","CHGA","MKI67","OLFM4","APOC3","AQP10"))
plotFeature(seu.obj=seu_obj, dr="umap", col.num=5, genes.to.plot=c("EPCAM", "PTPRC", "PDX1", "ONECUT2", "CDX2", "OSR2", "HOXB7", "HOXA10","SATB2","FABP6", "APOA4","BEST4","RBPJ","LGR5","","MUC2","CHGA","MKI67","OLFM4","APOC3","AQP10"))

SCpubr::do_FeaturePlot(seu_obj, features=c("nFeature_RNA"), order=T)
SCpubr::do_FeaturePlot(seu_obj, features="FABP1", order=T)
# check nFeature distribution

seu_obj$log_nCount <- log(seu_obj$nCount_RNA)
p1 <- SCpubr::do_RidgePlot(sample = seu_obj,
                          feature = "log_nCount",
                          group.by="RNA_snn_res.1")
p2 <- SCpubr::do_RidgePlot(sample = seu_obj,
                          feature = "nFeature_RNA",
                          group.by="RNA_snn_res.1")
p3 <- SCpubr::do_RidgePlot(sample = seu_obj,
                           feature = "percent.mt",
                           group.by="RNA_snn_res.1")
p1+p2+p3
# identify cell type markers
de_res <- presto::wilcoxauc(seu_obj, group_by="RNA_snn_res.1")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
top_res <- sig_res %>% group_by(group) %>% top_n(3, wt=logFC)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_adult_chip_TRM_coculture_cell_type_markers.rds")

genes <- unique(top_res$feature)
p2 <- SCpubr::do_DotPlot(sample = seu_obj, 
             features = genes,
             scale=T)
p2
ggsave(p2, filename="Plot_dotplot_adult_chip_TRM_coculture_cell_type_markers.pdf", width=15, height=9)

top_res <- sig_res %>% group_by(group) %>% top_n(10, wt=logFC)
genes <- unique(top_res$feature[which(top_res$group==15)])
plotFeature(seu.obj=seu_obj, dr="umap", col.num=5, genes.to.plot=genes, do.plot=F)


SCpubr::do_DotPlot(sample = seu_obj, 
             features = c("HOXB7", "HOXA10"),
             scale=T)
# first annotation
cl_vec <- paste0("C", seu_obj$RNA_snn_res.1)
ct_vec <- setNames(c("BEST4+","Immune","TA-1","E-1","EEC","TA-2_G2M_GC","Low quality","E-2","E-3","SC-1","GC","TA-3","E-4","TA-4","SC-2","E-5","E-6","SC-3","E-7","E-8","E-9","E-10"),
                    paste0("C", 21:0))

seu_obj$anno_v1 <- ct_vec[cl_vec]
seu_obj$anno_v1_coarse <- seu_obj$anno_v1
seu_obj$anno_v1_coarse[grep("TA-", seu_obj$anno_v1_coarse)] <- "TA"
seu_obj$anno_v1_coarse[grep("SC-", seu_obj$anno_v1_coarse)] <- "SC"
seu_obj$anno_v1_coarse[grep("E-", seu_obj$anno_v1_coarse)] <- "Enterocyte"
p3 <- SCpubr::do_DimPlot(seu_obj, group.by = "anno_v1_coarse", pt.size=6, label=T, label.size=20, border.size=1.5)+NoLegend()
png("Plot_UMAP_colored_by_anno_v1.png", width = 2000, height = 2000)
p3
dev.off()
saveRDS(seu_obj, file="Dat_merged_filtered_ileum_epithelium_TRM_coculture.rds")
trm_epi_obj <- seu_obj
SCpubr::do_DimPlot(trm_epi_obj, group.by = "anno_v1_coarse", pt.size=6, label=T, label.size=20, border.size=1.5)+NoLegend()


# run DE analysis between conditions
setwd("/home/yuq22/Bacteria_TRM/TRM/")
library(Seurat)
trm_epi_obj <- readRDS("Dat_merged_filtered_ileum_epithelium_TRM_coculture.rds")
# exclude immune cells and low quality cells
seu_obj <- subset(trm_epi_obj, cells=colnames(trm_epi_obj)[which(!trm_epi_obj$anno_v1_coarse %in% c("Low quality", "Immune"))])
SCpubr::do_DimPlot(seu_obj, group.by = "anno_v1_coarse", pt.size=6, label=T, label.size=20, border.size=1.5)+NoLegend()
saveRDS(seu_obj, file="/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/TRM/Dat_merged_filtered_ileum_epithelium_TRM_coculture_noImmune_noLowQuality.rds")

seu_obj <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/TRM/Dat_merged_filtered_ileum_epithelium_TRM_coculture_noImmune_noLowQuality.rds")
# focus on stem cells to enterocytes
dir.create("stem_to_enterocyte")
setwd("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/TRM/stem_to_enterocyte")
cells <- colnames(seu_obj)[which(seu_obj$anno_v1_coarse %in% c("Enterocyte", "SC"))]
seu_obj <- subset(seu_obj, cells=cells)
SCpubr::do_DimPlot(se_data, group.by = "anno_v1_coarse", pt.size=6, label=T, label.size=20, border.size=1.5)+NoLegend()
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 3000)
seu_obj <- ScaleData(object = seu_obj, verbose = T)
seu_obj <- RunPCA(object = seu_obj, features = VariableFeatures(seu_obj), verbose = F, npcs = 50)
usefulPCs <- 1:20
seu_obj <- FindNeighbors(object = seu_obj, dims = usefulPCs)
seu_obj <- FindClusters(object = seu_obj, resolution = 1)
seu_obj <- RunUMAP(object = seu_obj, dims = usefulPCs)
saveRDS(seu_obj, file="Dat_merged_filtered_ileum_epithelium_TRM_coculture_stem_cells_to_enterocyte_trajectory.rds")

# check tissue data
duo_se <- readRDS("/home/yuq22/Bacteria_TRM/intestine_multi_region_scRNA-seq_atlas/markers_whole_dataset/subclustering_per_coarse_cell_type/s2E_per_region_per_paper/smallInt/Res_cleaned_adult_human_SmallInt_stem_cell_and_enterocyte_with_CSS_integration.rds")
plotFeature(seu.obj=duo_se, dr="umap_css", col.num=5, genes.to.plot=c("HOXB7", "HOXA10","FABP6", "APOA4","LGR5","MKI67","OLFM4","APOC3","AQP10"),
plot.name="Plot_UMAP_CSS_adult_duo_tissue_stem2Ent_genes.png")

# run diffusion map on each sample separately
chip_se <- seu_obj
dm <- DiffusionMap(Embeddings(chip_se, "pca")[,1:20])
saveRDS(dm, file=paste0("Res_adult_chip_s2E_diffusion_map.rds"))
dpt <- DPT(dm)
chip_se$dpt <- rank(dpt$dpt)
saveRDS(chip_se,file="Dat_merged_filtered_ileum_epithelium_TRM_coculture_stem_cells_to_enterocyte_trajectory_with_dpt.rds" )

p1 <- SCpubr::do_DimPlot(chip_se, group.by = "orig.ident", pt.size=6, font.size=100, border.size=1.5, legend.icon.size=40)
p2 <- SCpubr::do_DimPlot(chip_se, group.by = "RNA_snn_res.1", pt.size=6, label=T, label.size=20, border.size=1.5)+NoLegend()
p3 <- SCpubr::do_DimPlot(chip_se, group.by = "anno_v1_coarse", pt.size=6, label=T, label.size=20, border.size=1.5)+NoLegend()
p4 <- SCpubr::do_FeaturePlot(chip_se, feature = "dpt", pt.size=6, order=T, border.size=1.5)+NoLegend()
png("Plot_UMAP_TRM_chip_coculture_stem_cell_enterocyte.png", width = 2000*2, height = 2000*2)
(p1+p3)/(p2+p4)
dev.off()

plotFeature(seu.obj=chip_se, dr="umap", col.num=4, genes.to.plot=c("HOXB7", "HOXA10","FABP6", "APOA4","LGR5","OLFM4","APOC3","AQP10"))

tissue_folder <- "/home/yuq22/Bacteria_TRM/intestine_multi_region_scRNA-seq_atlas/markers_whole_dataset/subclustering_per_coarse_cell_type/s2E_per_region_per_paper/smallInt"
ds3 <- readRDS(paste(tissue_folder, "Res_adult_human_SmallInt_stem_cell_and_enterocyte_Burclaff_DUO_with_CSS_integration.rds", sep="/"))
plotFeature(seu.obj=ds3, dr="umap_css", col.num=3, genes.to.plot=c("ONECUT2","PDX1","HOXA10","HOXB7","FABP6","APOA4"), plot.name="Plot_UMAP_CSS_Burclaff_DUO_tissue_stem2Ent_genes.png")

# perform DEG analysis between conditions
de_res <- presto::wilcoxauc(chip_se, group_by="orig.ident")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
top_res <- sig_res %>% group_by(group) %>% top_n(20, wt=logFC)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_s2E_DEG_with_and_without_TRM.rds")

genes <- unique(top_res$feature)
p2 <- SCpubr::do_DotPlot(sample = chip_se, 
             features = genes,
             scale=T,
             group.by="orig.ident")
p2
ggsave(p2, filename="Plot_dotplot_as2E_DEG_with_and_without_TRM.pdf", width=15, height=9)

seu_obj <- chip_se
seu_obj <- FindClusters(object = seu_obj, resolution = 2)
p1 <- SCpubr::do_DimPlot(seu_obj, group.by = "RNA_snn_res.2", pt.size=6, label=T, label.size=20, border.size=1.5)+NoLegend()
png("Plot_UMAP_s2E_reso2_cluster.png", height=2000, width=2000)
p1
dev.off()

# perform DEG analysis between conditions
detected_genes <- rownames(chip_se@assays$RNA@counts)[rowSums(chip_se@assays$RNA@counts)>0]
confound_genes <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/confound_genes/Dat_confound_genes.rds")
input_genes <- setdiff(detected_genes, confound_genes)
tf_genes <- unique(read.table("/projects/site/pred/ihb-intestine-evo/Annotation/HumanTFDB/List_HumanTFDB_TF.txt", sep="\t", stringsAsFactors=F, head=T)$Symbol)
input_genes <- intersect(detected_genes, tf_genes)
library(doParallel)
registerDoParallel(20)
con_vec <- as.factor(chip_se$orig.ident)
cl_vec <- as.factor(chip_se$RNA_snn_res.1)
umi_vec <- as.numeric(chip_se$log_nCount)
p.value <- foreach(gene = input_genes, .combine = c) %dopar% {
    e <- as.vector(as.matrix(chip_se@assays$RNA@counts[gene, ]))
    m0 <- lm(e ~ umi_vec+cl_vec)
    m1 <- lm(e ~ umi_vec+cl_vec+con_vec)
    p_anova <- anova(m1,m0)$Pr[2]
    return(p_anova)
}
stopImplicitCluster()
names(p.value) <- input_genes
saveRDS(p.value, file="Res_ANCOVA_TF_pvalues.rds")  

padj_vec <- p.adjust(p.value, method="BH")
sum(padj_vec<0.05)
# volcano plot for TFs
con_expr <- getAveExpr(seu.obj=chip_se, feature.to.calc="orig.ident", colname.prefix=NULL)
saveRDS(con_expr, file="Dat_s2E_TRM_coculture_condition_expression.rds")

x <- con_expr[input_genes,"withTRM"] - con_expr[input_genes,"noTRM"]
y <- -log10(p.value)
plot(x, y, pch=20, col=ifelse(padj_vec<0.05, "red", "black"), xlab="log2FC", ylab="-log10(p-value)")

df <- data.frame("feature"=input_genes, "logNominalPval"=y, "BHadjP"=padj_vec, "exprLogFC"=x, stringsAsFactors=F)
saveRDS(df, file="Res_ANCOVA_TF_volcano.rds")

up_df <- df[df$BHadjP<0.05 & df$exprLogFC>0,]
top_up_genes <- up_df$feature[order(up_df$exprLogFC, decreasing=T)[1:10]]

down_df <- df[df$BHadjP<0.05 & df$exprLogFC<0,]
top_down_genes <- down_df$feature[order(down_df$exprLogFC)[1:10]]

highlight_genes <- c(top_up_genes, top_down_genes)

library(ggrepel)
p1 <- ggplot(df, aes(x=exprLogFC, y=logNominalPval, color=BHadjP<0.05))+
    geom_point()+
    theme_bw()+
    theme(legend.position="none")+
    scale_color_manual(values=c("#969696", "#202020"))+
    ggrepel::geom_text_repel(data=df[df$feature %in% highlight_genes,], aes(label=feature), box.padding=0.5, point.padding=0.5, segment.color="grey50")+
    geom_point(data=df[df$feature %in% top_up_genes,], color="red")+
    geom_point(data=df[df$feature %in% top_down_genes,], color="blue")+
    xlab("expr_logFC")+
    ylab("-log10(p-value)")

pdf("Plot_TF_volcano.pdf", width=10, height=10)
p1
dev.off()

# check some gene expression in tissue data
# check cluster distribution per condition
chip_se <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/TRM/stem_to_enterocyte/Dat_merged_filtered_ileum_epithelium_TRM_coculture_stem_cells_to_enterocyte_trajectory_with_dpt.rds")
SCpubr::do_FeaturePlot(chip_se, feature = "HLA-DRA", pt.size=6, order=T, border.size=1.5)+NoLegend()
n1 <- sapply(sort(unique(chip_se$RNA_snn_res.1)), function(i){
    sapply(sort(unique(chip_se$orig.ident)), function(y){
        sum(chip_se$RNA_snn_res.1==i & chip_se$orig.ident==y)
    })
})
colnames(n1) <- paste0("C", sort(unique(chip_se$RNA_snn_res.1)))
prop1 <- n1/rowSums(n1)
prop_diff <- prop1[2,] - prop1[1,]
cl_expr <- getAveExpr(seu.obj=chip_se, feature.to.calc="RNA_snn_res.1")
saveRDS(cl_expr, file="Res_s2e_cluster_expression_per_condition.rds")
cluster_order <- names(sort(prop_diff))

pdf("Plot_diff_cluster_distribution_per_condition.pdf", height=10, width=10)
par(mfrow=c(2,1))
barplot(cl_expr["HOXB7", cluster_order], las=2, main="HOXB7 expression per cluster", ylab="Normed. Expr", xlab="Cluster")
barplot(sort(prop_diff), las=2, ylab="Prop diff (withTRM - noTRM)", xlab="Cluster", main="Cluster distribution difference between withTRM and noTRM")
dev.off()

epi_obj <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/TRM/Dat_merged_filtered_ileum_epithelium_TRM_coculture_noImmune_noLowQuality.rds")
epi_obj$Cluster_per_condition <- paste(epi_obj$orig.ident, epi_obj$RNA_snn_res.1, sep="_")
epi_obj$Cell_type_per_condition <- paste(epi_obj$anno_v1_coarse, epi_obj$orig.ident, sep="_")
saveRDS(epi_obj, file="Dat_merged_filtered_ileum_epithelium_TRM_coculture_noImmune_noLowQuality_clustered.rds")

# run presto between conditions per cell type
setwd("/home/yuq22/Bacteria_TRM/TRM/")
combined <- readRDS("/home/yuq22/Bacteria_TRM/TRM/Dat_merged_filtered_ileum_epithelium_TRM_coculture_noImmune_noLowQuality_clustered.rds")
gene <- rownames(combined)[!grepl("^ENSG", rownames(combined))]
#remove ribosomal genes, mitochondrial genes, and genes located on sex chromosomes before normalization
confound_genes <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/confound_genes/Dat_confound_genes.rds")
genes_to_exclude <- unique(unlist(confound_genes[c(4,5,6)]))
genes <- setdiff(gene, genes_to_exclude)

condition_deg_res <- list()
for(ct in sort(unique(combined$anno_v1_coarse))){
    
    print(paste(ct, 'start'))
    idx <- which(combined$anno_v1_coarse==ct)
    X <- combined@assays$RNA@data[genes,idx]
    y <- combined$Cell_type_per_condition[idx]

    de_res <- presto::wilcoxauc(X=X, y=y)
    de_res$pct_diff <- de_res$pct_in - de_res$pct_out
    pval_log <- -log10(de_res$pval)
    pval_log[is.infinite(pval_log)] <- max(pval_log[!is.infinite(pval_log)])+1
    de_res$pval_transformed <- pval_log^0.25
    de_res$auc_sym <- abs(de_res$auc-0.5)
    condition_deg_res[[ct]] <- de_res
}
df <- data.frame(do.call('rbind', condition_deg_res))
saveRDS(df, file="Res_TRM_coculture_effect_on_duo_epi_per_cell_type_DEG_all_res_combined.rds")
# get significant DEG
condition_sig_res <- list()
for(ct in sort(unique(combined$anno_v1_coarse))){
    
    print(paste(ct, 'start'))
    idx <- which(combined$anno_v1_coarse==ct)
    X <- combined@assays$RNA@data[genes,idx]
    y <- combined$Cell_type_per_condition[idx]

    de_res <- presto::wilcoxauc(X=X, y=y)
    de_res$pct_diff <- de_res$pct_in - de_res$pct_out
    if(ncol(X)<1000){
        sig_res <- de_res %>% filter(pval<0.05 & pct_in>5 & logFC>0.05)
    }else{
        sig_res <- de_res %>% filter(padj<0.05 & pct_in>5 & logFC>0.05)
    }
    condition_sig_res[[ct]] <- sig_res
}
df <- data.frame(do.call('rbind', condition_sig_res))
table(df$group)
saveRDS(df, file="Res_TRM_coculture_effect_on_duo_epi_per_cell_type_DEG_sig_res_combined.rds")

# more relaxing cutoffs
epi_obj <- readRDS("/home/yuq22/Bacteria_TRM/TRM/Dat_merged_filtered_ileum_epithelium_TRM_coculture_noImmune_noLowQuality_clustered.rds")
condition_deg_list <- list()
for(ct in sort(unique(epi_obj$anno_v1_coarse))){
    print(ct)
    seu_obj <- subset(epi_obj, anno_v1_coarse==ct)
    condition_de_res <- presto::wilcoxauc(seu_obj, group_by="orig.ident")
    condition_de_res$pct_diff <- condition_de_res$pct_in - condition_de_res$pct_out
    if(ncol(seu_obj)<1000){
        condition_sig_res <- condition_de_res %>% filter(pval<0.05 & pct_diff>1 & logFC>0.01)
    }else{
        condition_sig_res <- condition_de_res %>% filter(padj<0.05 & pct_diff>1 & logFC>0.01)
    }
    saveRDS(condition_sig_res, file=paste0("Res_", ct, "_DEG_with_and_without_TRM.rds"))
    condition_deg_list[[ct]] <- unique(condition_sig_res$feature)
}
saveRDS(condition_deg_list, file="Res_cell_type_DEG_with_and_without_TRM.rds")
sapply(condition_deg_list, length)
condition_deg <- unique(unlist(condition_deg_list))
confound_genes <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/confound_genes/Dat_confound_genes.rds")
condition_deg <- setdiff(condition_deg[!grepl("ENSG", condition_deg)], unlist(confound_genes))
length(condition_deg)


condition_deg_list <- readRDS("/home/yuq22/Bacteria_TRM/TRM/Res_cell_type_DEG_with_and_without_TRM.rds")
sapply(names(condition_deg_list), function(x){

})


# cluster the DEG profile
ave_expr <- getAveExpr(seu.obj=epi_obj, feature.to.calc="Cell_type_per_condition", colname.prefix=NULL)
saveRDS(ave_expr, file="Res_cell_type_expression_per_condition.rds")
expr_diff <- ave_expr[,grep("withTRM",colnames(ave_expr))] - ave_expr[,grep("noTRM",colnames(ave_expr))]
saveRDS(expr_diff, file="Res_cell_type_expression_diff_withTRM_over_noTRM_per_condition.rds")

expr_mat <- ave_expr[condition_deg,]
cl_enriched_de_idx <- apply(expr_mat, 1, function(vec){
    idx <- sapply(seq(length(vec)), function(i){
        vec[i]>mean(vec[-i])+3*sd(vec[-i]) | vec[i]<mean(vec[-i])-3*sd(vec[-i])
    })
    sum(idx>0)
})
genes <- names(cl_enriched_de_idx)[cl_enriched_de_idx>0]
length(genes)
saveRDS(genes, file="Res_cell_type_enriched_DEG.rds")
deg_list <- list("all"=condition_deg, "cell_type_enriched"=genes)
saveRDS(deg_list, file="Res_per_cell_type_wilcox_test_DEG.rds")
input_expr_mat <- expr_mat[genes,]
hc <- hclust(as.dist(1-cor(t(input_expr_mat))), method="ward.D2")
plot(hc, hang=-1)
abline(h=7)
g <- cutree(hc, h=7)
mat <- do.call('rbind', strsplit(split="_", colnames(input_expr_mat)))

res <- plotClusterExprProfile(expr = input_expr_mat, time.vec = mat[,1], time.vec.type="categoric", group.vec = mat[,2], 
 cluster.vec = g, group.cols = c("#31a354", "#3182bd"), return.value = T, to.plot = T, plot.name = "Plot_DEG_cluster_profile.pdf", 
add.legend = T, legend.pos = "topleft", cex.legend = 2, col.num = 3, 
border.do.smooth = F, mean.do.smooth = F, df = 8, xlab = "Cluster", 
ylab = "Relative expr.") 
saveRDS(res, file="Res_DEG_cluster_profile.rds")

 
# get example DEG genes per cell type - also as sanity check
top_logFC_genes_per_cell_type <- apply(expr_diff[genes,], 2, function(vec){
    up_genes <- genes[order(vec, decreasing=T)[1:3]]
    down_genes <- genes[order(vec)[1:5]]
    return(c(up_genes, down_genes))
})
features <- unique(as.vector(top_logFC_genes_per_cell_type))
p2 <- SCpubr::do_DotPlot(sample = epi_obj, 
             features = features,
             group.by="Cell_type_per_condition",
             scale=T)
p2
ggsave(p2, filename="Plot_dotplot_adult_chip_TRM_coculture_DEGs.pdf", width=10, height=6)



# identify cell type markers
de_res <- presto::wilcoxauc(epi_obj, group_by="anno_v1")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
top_res <- sig_res %>% group_by(group) %>% top_n(3, wt=logFC)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_adult_chip_TRM_coculture_noLowQuality_noImmune_cell_type_markers.rds")

# ligand-receptor gene expression
dir.create("ligand_receptor")
setwd("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/TRM/ligand_receptor")
load("/home/yuq22/ihb-intestine-evo/Annotation/CellChat_2021/CellChatDB.human.rda")
ligand_gene_complex <- intersect(rownames(CellChatDB.human$complex), unique(CellChatDB.human$interaction$ligand))
receptor_gene_complex <- intersect(rownames(CellChatDB.human$complex), unique(CellChatDB.human$interaction$receptor))
ligand_genes <- intersect(union(unique(CellChatDB.human$interaction$ligand), unique(unlist(CellChatDB.human$complex[ligand_gene_complex,]))), rownames(epi_obj@assays$RNA@counts))
receptor_genes <- intersect(union(unique(CellChatDB.human$interaction$receptor), unique(unlist(CellChatDB.human$complex[receptor_gene_complex,]))), rownames(epi_obj@assays$RNA@counts))
ligand_receptr_genes <- list("ligand"=ligand_genes, "receptor"=receptor_genes)
saveRDS(ligand_receptr_genes, file="Res_expressed_ligand_receptor_genes.rds")

epi_obj$anno_v1 <- factor(epi_obj$anno_v1, levels=sort(unique(epi_obj$anno_v1))) 
Idents(epi_obj) <- epi_obj$anno_v1
top_ligand <- sig_res[which(sig_res$feature%in%ligand_genes),] %>% group_by(group) %>% top_n(3, wt=logFC)
features <- unique(top_ligand$feature)
p2 <- SCpubr::do_DotPlot(sample = epi_obj, 
             features = features,
             scale=T)
p2
ggsave(p2, filename="Plot_dotplot_adult_chip_TRM_coculture_filtered_cell_type_enriched_ligand.pdf", width=10, height=6)

top_receptor <- sig_res[which(sig_res$feature%in%receptor_genes),] %>% group_by(group) %>% top_n(3, wt=logFC)
features <- unique(top_receptor$feature)
p2 <- SCpubr::do_DotPlot(sample = epi_obj, 
             features = features,
             scale=T)
p2
ggsave(p2, filename="Plot_dotplot_adult_chip_TRM_coculture_filtered_cell_type_enriched_receptor.pdf", width=10, height=6)

features <- intersect(genes, ligand_genes)
p2 <- SCpubr::do_DotPlot(sample = epi_obj, 
             features = features,
             group.by="Cell_type_per_condition",
             scale=T)
p2
ggsave(p2, filename="Plot_dotplot_adult_chip_TRM_coculture_DE_ligand_genes.pdf", width=10, height=6)

features <- intersect(genes, receptor_genes)
p2 <- SCpubr::do_DotPlot(sample = epi_obj, 
             features = features,
             group.by="Cell_type_per_condition",
             scale=T)
p2
ggsave(p2, filename="Plot_dotplot_adult_chip_TRM_coculture_DE_receptor_genes.pdf", width=10, height=6)


