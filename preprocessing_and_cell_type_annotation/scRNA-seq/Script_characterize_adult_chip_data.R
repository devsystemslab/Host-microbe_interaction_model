setwd("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/epithelium_only/")
library(Seurat)
library(ggplot2)
library(dplyr)

folder <- "/projects/site/pred/ihb-g-deco/USERS/lopezsar/Bacteria_TRM_Manuscript/1st Figure - multidonor hashing experiment/raw data - output from CellRanger"
ids <- list.files(path=folder)
seu_obj_list <- list()
for(x in ids){
    data <- Read10X(data.dir=paste0(folder, "/", x))
    count <- data$"Gene Expression"
    colnames(count) <- paste(x, colnames(count), sep="_")
    seu_obj <- CreateSeuratObject(count)
    seu_obj_list[[x]] <- seu_obj
}
combined_rna <- merge(x=seu_obj_list[[1]], y=seu_obj_list[-1])
combined_rna <- NormalizeData(object = combined_rna, normalization.method = "LogNormalize", scale.factor = 1e4)
# calculate Mt%
combined_rna[["percent.mt"]] <- PercentageFeatureSet(combined_rna, pattern = "^MT-")
saveRDS(combined_rna, file="Dat_merged_hashed_epithelium_only_scRNA-seq_data.rds")


# show Mt%, nFeature, nCount distribution per sample
# Visualize QC metrics as a violin plot
p1 <- SCpubr::do_ViolinPlot(combined_rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pdf("Plot_ViolinPlot_QC_metrics_per_sample.pdf" , width = 15, height = 5)
p1
dev.off()
# filter cells
cells <- colnames(combined_rna)[which(combined_rna$nFeature_RNA > 500 & combined_rna$nFeature_RNA < 6000 & combined_rna$nCount_RNA<30000 & combined_rna$percent.mt < 25)]
seu_obj <- subset(combined_rna, cells=cells)
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 3000)
seu_obj <- ScaleData(object = seu_obj, verbose = T)
seu_obj <- RunPCA(object = seu_obj, features = VariableFeatures(seu_obj), verbose = F, npcs = 50)
usefulPCs <- 1:20
seu_obj <- FindNeighbors(object = seu_obj, dims = usefulPCs)
seu_obj <- FindClusters(object = seu_obj, resolution = 1)
seu_obj <- RunUMAP(object = seu_obj, dims = usefulPCs)
saveRDS(seu_obj, file="Dat_hashed_epithelium_only_scRNA-seq_data_filtered_and_merged.rds")

# exlcude the potential foregut sample and the IleColMixed sample
setwd("/home/yuq22/Bacteria_TRM/epithelium_only/exclude_foregut_and_IleColMixed_samples")
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")
library(Seurat)
library(ggplot2)
library(dplyr)

seu_obj <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/epithelium_only/Dat_hashed_epithelium_only_scRNA-seq_data_filtered_and_merged.rds")
# "don55" is potential foregut sample, "IF17ile" shows mixed ileum and colon
sample_to_exclude <- c("don55", "IF17ile")
filtered_cells <- colnames(seu_obj)[which(!seu_obj$orig.ident%in%sample_to_exclude)]
seu_obj <- subset(seu_obj, cells=filtered_cells)
# clean up the genes without gene symbol
count <- seu_obj@assays$RNA@counts
gene <- rownames(count)[!grepl("^ENSG", rownames(count))]
#remove ribosomal genes, mitochondrial genes, and genes located on sex chromosomes before normalization
confound_genes <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/confound_genes/Dat_confound_genes.rds")
genes_to_exclude <- unique(unlist(confound_genes[c(4,5,6)]))
genes <- setdiff(gene, genes_to_exclude)
count <- count[genes,]
meta <- seu_obj@meta.data
seu_obj <- CreateSeuratObject(counts=count, meta.data=meta)
seu_obj <- NormalizeData(object = seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 3000)
seu_obj <- ScaleData(object = seu_obj, verbose = T)
seu_obj <- RunPCA(object = seu_obj, features = VariableFeatures(seu_obj), verbose = F, npcs = 50)
usefulPCs <- 1:20
seu_obj <- FindNeighbors(object = seu_obj, dims = usefulPCs)
seu_obj <- FindClusters(object = seu_obj, resolution = 1)
seu_obj <- RunUMAP(object = seu_obj, dims = usefulPCs)
SCpubr::do_DimPlot(seu_obj, reduction="umap", group.by="orig.ident", label=T)
saveRDS(seu_obj, file="Dat_hashed_intestine_epithelium_only_scRNA-seq_data.rds")


# plot regional markers
genes <- c("CDX2", "PDX1", "ONECUT2", "OSR2", "SATB2", "HOXA10", "SOX2","CLDN18","HNF4A", "KLF2","NKX6-2","FOXQ1","GATA4","GATA6","PTF1A","PROX1","NKX6-1","ONECUT1") # regional markers
plotFeature(seu.obj=seu_obj, dr="umap", genes.to.plot=genes, col.num=6, plot.name="Plot_UMAP_regional_markers_merged_data_gene_cleaned_exclude_IF17ile.png")
# plot cell type markers
genes <- c(
    "FABP6","SLC10A2",
    "FABP2","SI",
    "ADA","RBP2",
    "APOA4","ALDOB",
    "LGR5","OLFM4","REG1A","SMOC2","ASCL2",
    "MUC2","FCGBP",
    "MUC5AC","MUC6","MUC5B",
    "SPDEF","DLL1","NEUROG3",
    "CHGA", "ISL1","CHGB",
    "GUCA2A","GUCA2B","MEIS1","BEST4","CA7",
    "GP2","SPIB",
    "AVIL", "RGS13",
    "DEFA5","DEFA6",
    "CEACAM6","CEACAM7",
    "MKI67",
    "GLP1R","GHRL" # brunner's gland 
    )
plotFeature(seu.obj=seu_obj, dr="umap", genes.to.plot=genes, col.num=8, plot.name="Plot_UMAP_epithelial_cell_type_markers_merged_data_gene_cleaned_exclude_IF17ile.png")
genes <- c(
    "GATA4", "SATB2", "MKI67",
    "APOA4", "FABP6","CEACAM7",
    "LGR5", "MUC2", "CHGA"
    )
plotFeature(seu.obj=seu_obj, dr="umap", genes.to.plot=genes, col.num=3, plot.name="Plot_UMAP_epithelial_cell_type_markers_merged_data_gene_cleaned_exclude_IF17ile.png", nCols=beach.col, per.plot.size=2000, cex=5)

# mucin genes
genes <- sort(paste0("MUC", c(1,13,20,4,7,"3A",12,17,6,2,"5AC","5B")))
plotFeature(seu.obj=seu_obj, dr="umap", genes.to.plot=genes, col.num=6, plot.name="Plot_UMAP_mucin_genes_chip_no_don55_and_IFile17.png", nCols=beach.col, per.plot.size=2000, cex=5)

# zonation genes
zonation_genes <-  list(
    ent_zonation_genes = c("REG1A", "DMBT1", "GSTA1","PIGR", "RBP2","DGAT1", "SLC5A1", "SLC2A5","APOA4", "AQP10"),
    EEC_zonation_genes = c("TPH1", "GCG", "SCT", "CCK", "MLN"),
    goblet_zonation_genes = c("SPINK4","MUC2", "FCGBP", "CLCA1", "ZG16"),
    tuft_zonation_genes = c("TRPM5", "SH2D6", "PLCG2", "GABRP", "DEFB1")
)
for(x in names(zonation_genes)){
    plotFeature(seu.obj=seu_obj, dr="umap", genes.to.plot=zonation_genes[[x]], col.num=6, plot.name=paste0("Plot_UMAP_", x, "_chip_no_don55_and_IFile17.png"), nCols=beach.col, per.plot.size=2000, cex=5)
}
# plot individual sample info
sample_cols <- setNames(c(rep("#4E79A7",2), "#A0CBE8", rep("#138d75",4)),c("OE012", "OE008", "OE015", "OE004", "OE014", "IF17col", "don47"))
p1 <- SCpubr::do_DimPlot(seu_obj, reduction="umap", group.by="orig.ident", colors.use=sample_cols, pt.size=8, border.size = 1.5)+NoLegend()
png("Plot_UMAP_chip_combined_individual.png", height=2000, width=2000)
p1
dev.off()

# coarse grain cell type annotation on cluster level
seu_obj <- readRDS("/home/yuq22/Bacteria_TRM/used_object/Dat_hashed_intestine_epithelium_only_scRNA-seq_data_with_cell_type_anno.rds")
# get cluster average expression
cluster_expr <- getAveExpr(seu.obj=seu_obj, feature.to.calc="RNA_snn_res.1")
barplot(cluster_expr["LGR5",], las=2)
SCpubr::do_DimPlot(seu_obj, reduction="umap", group.by="RNA_snn_res.1", label=F, label.size=10, split.by="RNA_snn_res.1")+NoLegend()
anno <- list(
"Prox_SI_SC"=8,
"Prox_SI_Ent"=c(2,3,9),
"Dist_SI_Ent"=c(1,11,19),
"Colon_SC"=c(5, 12, 18),
"Colon_Ent"=c(0,6,7,10,14,15,17),
"TA"=4,
"Goblet"=13,
"EEC"=16
)
seu_obj$Coarse_grain_cell_type <- NA
for(x in names(anno)){
    seu_obj$Coarse_grain_cell_type[which(seu_obj$RNA_snn_res.1%in%anno[[x]])] <- x
}
SCpubr::do_DimPlot(seu_obj, reduction="umap", group.by="Coarse_grain_cell_type", label=T, label.size=10)+NoLegend()
saveRDS(seu_obj, file="Dat_hashed_intestine_epithelium_only_scRNA-seq_data_with_cell_type_anno.rds")

ct_cols <- list(
"Prox_SI_SC"="#CA6778",
"Colon_SC"="#A84798",
"Prox_SI_Ent"="#A4D371",
"Dist_SI_Ent"="#7BCAA4",
"Colon_Ent"="#FDD884",
"TA"="#5E4FA2",
"Goblet"="#4DA7B0",
"EEC"="#197636"
)
ct_cols <- setNames(unlist(ct_cols), names(ct_cols))
SCpubr::do_DimPlot(seu_obj, reduction="umap", group.by="Coarse_grain_cell_type", label=T, label.size=10, colors.use=ct_cols)+NoLegend()
cols <- list(
    "sample"=sample_cols,
    "cell_type"=ct_cols
)
saveRDS(cols, file="cols.rds")

p1 <- SCpubr::do_DimPlot(seu_obj, reduction="umap", group.by="Coarse_grain_cell_type", colors.use=ct_cols, pt.size=8, border.size = 1.5)+NoLegend()
png("Plot_UMAP_chip_combined_cell_type.png", height=2000, width=2000)
p1
dev.off()


chip_data <- readRDS("Dat_hashed_intestine_epithelium_only_scRNA-seq_data_with_cell_type_anno.rds")
seu_obj <- subset(chip_data, orig.ident=="OE008")
p1 <- SCpubr::do_FeaturePlot(seu_obj, reduction="umap", feature="LCT", order=T, pt.size=4)
seu_obj2 <- subset(chip_data, orig.ident=="OE012")
p2 <- SCpubr::do_FeaturePlot(seu_obj2, reduction="umap", feature="LCT", order=T, pt.size=4)
p1+p2
region <- list(
    "prox_SI"=c("OE012", "OE008"),
    "Dist_SI"="OE015",
    "Colon"=c("OE004", "OE014", "IF17col", "don47")
)
chip_data$region <- NA
for(x in names(region)){
    chip_data$region[which(chip_data$orig.ident%in%region[[x]])] <- x
}
SCpubr::do_DimPlot(chip_data, reduction="umap", group.by="region", label=T, label.size=10)+NoLegend()
saveRDS(chip_data, file="/home/yuq22/Bacteria_TRM/used_object/Dat_hashed_intestine_epithelium_only_scRNA-seq_data_with_cell_type_anno.rds")

# adult chip de novo cluster marker
chip_data <- readRDS("/home/yuq22/Bacteria_TRM/used_object/Dat_hashed_intestine_epithelium_only_scRNA-seq_data_with_cell_type_anno.rds")
p1 <- SCpubr::do_DimPlot(chip_data, group.by="RNA_snn_res.1", label=T, pt.size=5, label.size=25, font.size=70, legend.icon.size = 20)+NoLegend()
p2 <- SCpubr::do_DimPlot(chip_data, group.by="orig.ident", label=F, pt.size=5, label.size=25, font.size=70, legend.icon.size = 20)+NoLegend()
png("Plot_UMAP_adult_chip_cluster.png" , width = 2000*2, height = 2000)
p1+p2
dev.off()
# identify cluster markers and generate feature plot
de_res <- presto::wilcoxauc(chip_data, group_by="RNA_snn_res.1")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
top_res <- sig_res %>% group_by(group) %>% top_n(3, wt=pct_diff)
top_res <- sig_res %>% group_by(group) %>% top_n(3, wt=auc)
top_res <- sig_res %>% group_by(group) %>% top_n(3, wt=logFC)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_presto_chip_data_cluster_marker.rds")
genes <- unique(top_res$feature)
genes <- c("SMOC2", "S100A4", "TFF1", "LAMA1", "APOA4", "ATF3", "PLCG2", "BNIP5", "CLCA1", "MUC5B")
plotFeature(seu.obj=chip_data, dr="umap", genes.to.plot=genes, do.plot=T, plot.name="Plot_UMAP_adult_chip_cluster_marker_feature_plot_for_duo_ent_subtype_markers.png", nCols=beach.col, per.plot.size=2000, cex=5, col.num=5)
plotFeature(seu.obj=tissue_obj, dr="umap", genes.to.plot=genes, do.plot=T, plot.name="Plot_UMAPCSS_adult_tissue_for_duo_ent_subtype_markers.png", nCols=beach.col, per.plot.size=2000, cex=5, col.num=5)


# resolve stem cells and absorptive cells
cells <- colnames(chip_data)[which(chip_data$Coarse_grain_cell_type%in%c("Colon_Ent", "Colon_SC", "Prox_SI_Ent", "Dist_SI_Ent", "Prox_SI_SC"))]
seu_obj <- subset(chip_data, cells=cells)
SCpubr::do_DimPlot(seu_obj)
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 3000)
seu_obj <- ScaleData(object = seu_obj, verbose = T)
seu_obj <- RunPCA(object = seu_obj, features = VariableFeatures(seu_obj), verbose = F, npcs = 50)
usefulPCs <- 1:20
seu_obj <- FindNeighbors(object = seu_obj, dims = usefulPCs)
seu_obj <- FindClusters(object = seu_obj, resolution = 0.5)
seu_obj <- RunUMAP(object = seu_obj, dims = usefulPCs)
p1 <- SCpubr::do_DimPlot(seu_obj, group.by="Coarse_grain_cell_type", label=T, pt.size=3, label.size=10, font.size=30, legend.icon.size = 10, legend.nrow=2)
p2 <- SCpubr::do_DimPlot(seu_obj, group.by="RNA_snn_res.0.5", label=T, pt.size=3, label.size=10, font.size=30, legend.icon.size = 10, legend.nrow=2)
saveRDS(seu_obj, file="Res_all_region_absorptive_and_stem_cells_merged_obj.rds")

anchor_genes <- c("LGR5","SMOC2","ASCL2")
ref_pattern <- colSums(t(scale(t(as.matrix(seu_obj@assays$RNA@data[anchor_genes,])))))
cor_vec <- cor(t(as.matrix(seu_obj@assays$RNA@data[VariableFeatures(seu_obj),])), ref_pattern)[,1]
genes <- names(sort(cor_vec, decreasing = T)[1:10])
plotFeature(seu.obj=chip_data, dr="umap", genes.to.plot=genes, col.num=5, do.plot=T, plot.name="Plot_UMAP_stem_marker_feature_plot.png", nCols=beach.col, per.plot.size=2000, cex=10)
seu_obj$expr_based_stem_score <- colSums(t(scale(t(as.matrix(seu_obj@assays$RNA@data[genes,])))))
SCpubr::do_FeaturePlot(seu_obj, feature="expr_based_stem_score", pt.size=5, order=T)
score2 <- seu_obj$expr_based_stem_score
score2 <- score2[which(score2>0)]
plot(sort(score2), seq(length(score2)), pch=16)
cutoff2 <- quantile(score2, 0.6)
abline(v=cutoff2, col="red")
stem_cells <- colnames(seu_obj)[which(seu_obj$expr_based_stem_score>cutoff2)]
coor <- Embeddings(chip_data, reduction = "umap")[colnames(seu_obj),]
seu_obj[["full_umap"]] <- CreateDimReducObject(embeddings = coor, key="FULLUMAP_", assay=DefaultAssay(seu_obj))
p1 <- SCpubr::do_DimPlot(seu_obj, cells.highlight=TA_cells, pt.size=5, reduction="full_umap")+NoLegend()
p2 <- SCpubr::do_DimPlot(seu_obj, cells.highlight=stem_cells, pt.size=5, reduction="full_umap")+NoLegend()
p1/p2

ct_vec <- setNames(rep("Absorptive", ncol(seu_obj)), colnames(seu_obj))
ct_vec[stem_cells] <- "Stem"
seu_obj$score_based_cell_type_annotation <- ct_vec
SCpubr::do_DimPlot(seu_obj, reduction="full_umap", group.by="score_based_cell_type_annotation", pt.size=5, label=T, label.size=10, font.size=30, legend.icon.size = 10)

# smooth the annotation to the cluster level
seu_obj <- FindClusters(object = seu_obj, resolution = 5)
n1 <- sapply(sort(unique(seu_obj$score_based_cell_type_annotation)), function(x){
    sapply(sort(unique(seu_obj$RNA_snn_res.5)), function(y){
        sum(seu_obj$RNA_snn_res.5==y & seu_obj$score_based_cell_type_annotation==x)
    })
})
rownames(n1) <- paste0("C", sort(unique(seu_obj$RNA_snn_res.5)))
stem_cell_cl <- sub("C", "", rownames(n1)[which(n1[,"Stem"] > n1[,"Absorptive"])])
seu_obj$extended_score_based_cell_type_annotation  <- seu_obj$score_based_cell_type_annotation
seu_obj$extended_score_based_cell_type_annotation[which(seu_obj$RNA_snn_res.5%in%stem_cell_cl)] <- "Stem"
SCpubr::do_DimPlot(seu_obj, reduction="full_umap", group.by="extended_score_based_cell_type_annotation", pt.size=5, label=F, label.size=10, font.size=30, legend.icon.size = 10)

seu_obj$lineage <- NA
lineages <- list(
    "proxSI"=c("Prox_SI_Ent","Prox_SI_SC"),
    "distSI"=c("Dist_SI_Ent"),
    "colon"=c("Colon_Ent","Colon_SC")
)
for(x in names(lineages)){
    seu_obj$lineage[which(seu_obj$Coarse_grain_cell_type%in%lineages[[x]])] <- x
}
seu_obj$stem_cell_absorptive_lineage <- paste(seu_obj$extended_score_based_cell_type_annotation, seu_obj$lineage, sep="@")
SCpubr::do_DimPlot(seu_obj, group.by="stem_cell_absorptive_lineage", pt.size=5, label=T, label.size=10, font.size=30, legend.icon.size = 10)
saveRDS(seu_obj, file="Res_all_region_absorptive_and_stem_cells_merged_obj.rds")

# update the chip data annotation baesd on the stem cell score
chip_data$Updated_cell_type_annotation <- chip_data$Coarse_grain_cell_type
chip_data@meta.data[colnames(seu_obj), "Updated_cell_type_annotation"] <- seu_obj$stem_cell_absorptive_lineage
SCpubr::do_DimPlot(chip_data, group.by="Updated_cell_type_annotation", pt.size=2, label=T, label.size=10, font.size=30, legend.icon.size = 10)
# unique the stem cell and absorptive cell annotation
name_pairs <- list(
    "Absorptive@colon"="Colonocyte",
    "Absorptive@distSI"="Dist_SI_Ent",
    "Absorptive@proxSI"= "Prox_SI_Ent",
    "Stem@colon" =  "Colon_SC",     
    "Stem@distSI" =  "Dist_SI_SC",      
    "Stem@proxSI" = "Prox_SI_SC"
)
for(x in names(name_pairs)){
    chip_data$Updated_cell_type_annotation[which(chip_data$Updated_cell_type_annotation==x)] <- name_pairs[[x]]
}
saveRDS(chip_data, file="Res_chip_data_with_updated_cell_type_annotation.rds")
saveRDS(chip_data, file="/home/yuq22/Bacteria_TRM/used_object/Res_adult_chip_with_updated_cell_type_annotation.rds")


# load  tissue data
tissue_data <- readRDS("/home/yuq22/Bacteria_TRM/used_object/Res_combined_adult_human_epi_with_CSS_ct_anno_updated.rds")
SCpubr::do_FeaturePlot(tissue_data, reduction="umap_css", feature="MUC5B", order=T, pt.size=4)
genes <- sort(paste0("MUC", c(1,13,20,4,7,"3A",12,17,6,2,"5AC","5B")))
plotFeature(seu.obj=tissue_data, dr="umap_css", genes.to.plot=genes, col.num=6, plot.name="Plot_UMAP_adult_tissue_atlas_mucin_genes.png", nCols=beach.col, per.plot.size=2000, cex=5)
plotFeature(seu.obj=chip_data, dr="umap", genes.to.plot="MUC3B", do.plot=F)

# subset the goblet cell cluster
setwd("/home/yuq22/Bacteria_TRM/epithelium_only/exclude_foregut_and_IleColMixed_samples/goblet_cells")
chip_gc <- subset(chip_data, cells=colnames(chip_data)[which(chip_data$Coarse_grain_cell_type=="Goblet")])
chip_gc <- FindVariableFeatures(chip_gc, selection.method = "vst", nfeatures = 3000)
chip_gc <- ScaleData(object = chip_gc, verbose = T)
chip_gc <- RunPCA(object = chip_gc, features = VariableFeatures(chip_gc), verbose = F, npcs = 50)
usefulPCs <- 1:20
chip_gc <- FindNeighbors(object = chip_gc, dims = usefulPCs)
chip_gc <- FindClusters(object = chip_gc, resolution = 1)
chip_gc <- RunUMAP(object = chip_gc, dims = usefulPCs)
SCpubr::do_DimPlot(chip_gc, reduction="umap", group.by="region", label=T, pt.size=10)
SCpubr::do_DimPlot(chip_gc, reduction="umap", group.by="RNA_snn_res.1", label=T, pt.size=10, label.size=10)
saveRDS(chip_gc, file="Dat_adult_chip_goblet_cells.rds")

# exclude C5 (proliferative goblet cells) to get the non-proliferative goblet cells
chip_gc <- subset(chip_gc, RNA_snn_res.1!=5)
chip_gc <- FindVariableFeatures(chip_gc, selection.method = "vst", nfeatures = 3000)
chip_gc <- ScaleData(object = chip_gc, verbose = T)
chip_gc <- RunPCA(object = chip_gc, features = VariableFeatures(chip_gc), verbose = F, npcs = 50)
usefulPCs <- 1:20
chip_gc <- FindNeighbors(object = chip_gc, dims = usefulPCs)
chip_gc <- FindClusters(object = chip_gc, resolution = 1)
chip_gc <- RunUMAP(object = chip_gc, dims = usefulPCs)
p1 <- SCpubr::do_DimPlot(chip_gc, reduction="umap", group.by="region", label=T, pt.size=10)
p2 <- SCpubr::do_DimPlot(chip_gc, reduction="umap", group.by="RNA_snn_res.1", label=T, pt.size=10, label.size=10)
p1+p2
saveRDS(chip_gc, file="Dat_adult_chip_non_proliferative_goblet_cells.rds")

# identify cluster markers
chip_gc <- readRDS("Dat_adult_chip_non_proliferative_goblet_cells.rds")
png("Plot_UMAP_adult_chip_non_proli_GC_colored_by_clusters.png", height=2000, width=2000)
SCpubr::do_DimPlot(chip_gc, reduction="umap", group.by="RNA_snn_res.1", label=T, label.size=25, pt.size=25)+NoLegend()
dev.off()
cols <- readRDS("../cols.rds")
png("Plot_UMAP_adult_chip_non_proli_GC_colored_by_regions.png", height=2000, width=2000)
SCpubr::do_DimPlot(chip_gc, reduction="umap", group.by="region", pt.size=25, colors.use=cols$region)+NoLegend()
dev.off()
plotFeature(seu.obj=chip_gc, dr="umap", genes.to.plot=c("MUC3B"), do.plot=F)
de_res <- presto::wilcoxauc(chip_gc, group_by="RNA_snn_res.1")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out

sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
saveRDS(sig_res, file="Res_wilcoxauc_adult_chip_non_proli_GC_cluster_markers.rds")

# get average gene expression across different regions
ave_expr <- getAveExpr(seu.obj=chip_gc, feature.to.calc = "region", colname.prefix = NULL)
saveRDS(ave_expr, file="Dat_non_proli_GC_region_average_expr.rds")

ave_expr <- readRDS("Dat_non_proli_GC_region_average_expr.rds")
#cm <- readRDS("/Volumes/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/epithelium_only/exclude_foregut_and_IleColMixed_samples/goblet_cells/Res_wilcoxauc_adult_chip_non_proli_GC_cluster_markers.rds")
cm <- readRDS("~/Bacteria_TRM/epithelium_only/exclude_foregut_and_IleColMixed_samples/goblet_cells/Res_wilcoxauc_adult_chip_non_proli_GC_cluster_markers.rds")

# get top 20 genes
top_res <- cm %>% group_by(group) %>% top_n(20, wt=logFC)
top_genes <- unique(top_res$feature)

highlight_genes <- c("MUC5B", "TRPA1", "CLCA1", "MUC4", "MUC12","BTNL8", "REG4","TM4SF4")
genes <- union(top_genes, highlight_genes)
expr_mat <- ave_expr[genes,]
res <- prepareTreeAndHeatmapInput(expr.mat = expr_mat, hc.method = "ward.D2", genes.to.highlight = highlight_genes, norm.method = "quantile")
saveRDS(res, file="Res_adult_chip_GC_cluster_marker_gene_plus_selected_gene_expr_heatmap_input.rds")
pdf("Plot_heatmap_adult_chip_GC_cluster_marker_gene_plus_selected_gene_expr.pdf", height=10)
gplots::heatmap.2(res$highlight_input, col = beach.col.heatmap, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 0.2, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

hc <- hclust(as.dist(1-cor(expr_mat)), method="ward.D2")
saveRDS(hc, file="Res_hc_adult_chip_GC_cluster_marker_gene_plus_selected_gene_expr.rds")
pdf("Plot_hc_adult_chip_GC_cluster_marker_gene_plus_selected_gene_expr.pdf")
plot(hc, hang=-1)
dev.off()


# plot the region composition per cluster
cluster_order <- sapply(hc_res$labels[hc_res$order], function(x){
    strsplit(x, "_")[[1]][2]
})
num_mat <- sapply(cluster_order, function(cl){
    sapply(sort(unique(chip_gc$region)), function(x){
        length(which(chip_gc$RNA_snn_res.1==cl & chip_gc$region==x))
    })
})
prop_mat <- num_mat/rowSums(num_mat)

pdf("Plot_barplot_non_proli_GC_cluster_enrichment_per_region.pdf", height=5, width=7)
barplot(prop_mat, col=cols$region[rownames(prop_mat)], xlab="Cluster", ylab="Proportion", names.arg=cluster_order, cex.names=0.8)
dev.off()


# perform GO enrichemnt on the cluster markers
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

cluster_expr <- getAveExpr(seu.obj=chip_gc, feature.to.calc="RNA_snn_res.1", sep="")
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
# generate heatmap of glycosylation related cluster markers
genes <- intersect(GO_anno$HGNC.symbol[which(GO_anno$GO.term.name=="protein glycosylation")], sig_res$feature)
expr_mat <- cluster_expr[genes,]

res <- prepareTreeAndHeatmapInput(expr.mat = expr_mat, hc.method = "ward.D2")
#saveRDS(res, file="Res_fetal_combined_cell_type_marker_expression_in_fetal_and_tHIO.rds")
saveRDS(res, file="Res_glycosylation_related_non_proli_GC_cluster_marker_expr_cluster.rds")

pdf("Plot_heatmap_non_proli_goblet_cell_cluster_marker_associated_with_glycosylation.pdf", height=10)
gplots::heatmap.2(res$heatmap_input, col = beach.col.heatmap, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 0.5, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

# generate heatmap to represent the GO enrichment 
# select the top 10 enriched GO terms per cluster
padj_mat <- readRDS("Res_GO_enrichment_adult_chip_non_proli_GC_cluster_markers_padj_mat.rds")
sig_term_per_cluster <- lapply(seq(ncol(padj_mat)), function(j){
    all_sig_terms_per_cluster <- rownames(padj_mat)[which(padj_mat[,j]<0.05)]
    n <- min(c(length(all_sig_terms_per_cluster), 10))
    top_sig_terms_per_cluster <- all_sig_terms_per_cluster[order(padj_mat[all_sig_terms_per_cluster,j])[1:n]]
})
names(sig_term_per_cluster) <- colnames(padj_mat)
all_sig_term_per_cluster <- setdiff(unique(unlist(sig_term_per_cluster)), NA)
input <- -log10(padj_mat[all_sig_term_per_cluster,])

res <- prepareTreeAndHeatmapInput(expr.mat = input, hc.method = "ward.D2")
saveRDS(res, file="Res_non_proli_goblet_cell_cluster_top_enriched_GO_term_padj_heatmap_input.rds")

pdf("Plot_heatmap_non_proli_goblet_cell_cluster_top_enriched_GO_term_padj.pdf", height=7)
gplots::heatmap.2(res$heatmap_input, col = white2black, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 0.5, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

# get cell type average expression pattern of the cluster markers that are involved in the top enriched GO terms
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_statiscal_test.R")
chip_gc <- readRDS("Dat_adult_chip_non_proliferative_goblet_cells.rds")
expressed_genes <- rownames(chip_gc)[rowSums(chip_gc@assays$RNA@data)>0]
GO_anno <- read_GO_anno(expressed_genes=expressed_genes)
sig_res <- readRDS("Res_wilcoxauc_adult_chip_non_proli_GC_cluster_markers.rds")


# presumably goblet cell heterogeneity comes from regional heterogeneity and differentiation 
# we first study regional heterogeneity
# resolve the goblet cell heterogeneity across different intestinal regions
#library(destiny)
#chip_gc <- readRDS("Dat_adult_chip_goblet_cells.rds")
#dm <- DiffusionMap(Embeddings(chip_gc, "pca")[,1:20])
#saveRDS(dm, file="Res_adult_chip_goblet_cells_diffusion_map.rds")
#dpt <- DPT(dm)
#chip_gc$dpt <- rank(dpt$dpt)
#saveRDS(chip_gc, file="Dat_adult_chip_goblet_cells.rds")
#p1 <- SCpubr::do_FeaturePlot(chip_gc, "dpt", pt.size=10)
non_proli_gc <- readRDS("Dat_adult_chip_non_proliferative_goblet_cells.rds")
dm <- DiffusionMap(Embeddings(non_proli_gc, "pca")[,1:20])
saveRDS(dm, file="Dat_adult_chip_non-proliferative_goblet_cells.rds")
dpt <- DPT(dm)
non_proli_gc$dpt <- rank(dpt$dpt)
saveRDS(non_proli_gc, file="Dat_adult_chip_non_proliferative_goblet_cells.rds")
dc <- dm@eigenvectors
sum_sq_vec <- sapply(seq(ncol(dc)), function(j){
    m0 <- anova(aov(dc[,j]~non_proli_gc$region))
    m0$"Sum Sq"[1]
})
non_proli_gc@meta.data[,colnames(dc)] <- dc
saveRDS(non_proli_gc, file="Dat_adult_chip_non_proliferative_goblet_cells.rds")
plot_list <- list()
for(j in seq(ncol(dc))){
    plot_list[[j]] <- SCpubr::do_ViolinPlot(non_proli_gc, feature=colnames(dc)[j], group.by="region")
}
library(cowplot)
library(grid)
library(gridExtra)
p <- plot_grid(plotlist = plot_list, align = "hv",ncol = 5)
pdf("Plot_DCs_per_region.pdf" , width = 5*5, height = 5*4)
p
dev.off()
# get the DCs with top distinction between regions
# no matter whethr using DCs or rank of DCs, DC5 and DC2 are the top two DCs, which explain more variance between regions than DPT
# so we use DC5 to distinguish the regions
idx <- order(sum_sq_vec, decreasing=T)[c(1,2)]
coor <- dm@eigenvectors[,idx]
plotFeature2(coor=coor, values=non_proli_gc$region, add.legend=T, cex=5)
png(paste0("Plot_DC_",paste(idx, collapse="_"),"_adult_chip_non_proliferative_goblet_cells.png"), height=2000, width=2000*2)
par(mfrow=c(1,2))
plotFeature2(coor=coor, values=non_proli_gc$RNA_snn_res.1, add.label=T, cex=5)
plotFeature2(coor=coor, values=non_proli_gc$region, add.legend=T, cex=5)
dev.off()
non_proli_gc$rank_DC5 <- rank(dc[,5])
p2 <- SCpubr::do_FeaturePlot(non_proli_gc, "rank_DC5", pt.size=10) 
# identify genes significantly variable across DC5
library(splines)
library(doParallel)
registerDoParallel(20)
region_expr <- getAveExpr(seu.obj=non_proli_gc, feature.to.calc="region",colname.prefix=NULL)
saveRDS(region_expr, file="Dat_adult_chip_non_proli_GC_region_expr.rds")
sum_vec <- rowSums(region_expr)
detected_genes <- rownames(non_proli_gc)[which(sum_vec>0)]
cor_mat <- cor(region_expr[detected_genes,], ref_pattern)
max_cor <- colnames(cor_mat)[apply(cor_mat, 1, which.max)]
# create artificial regional variable patterns
ref_pattern <- rbind(diag(3), 1-diag(3))
colnames(ref_pattern) <- colnames(region_expr)
rownames(ref_pattern) <- c("Colon_spesific_pos", "Dist_SI_specific_pos", "Prox_SI_specific_pos", "Colon_specific_neg", "Dist_SI_specific_neg", "Prox_SI_specific_neg")
# correlate the gene expression with the regional variable patterns
cor_mat <- cor(t(region_expr[detected_genes,]), t(ref_pattern))
max_cor <- apply(cor_mat, 1, max)

# identify genes significantly variable across DC5

expr_mat <- as.matrix(non_proli_gc@assays$RNA@data[detected_genes,])
ps_vec <- non_proli_gc$rank_DC5
p.value <- foreach(gene = detected_genes, .combine = c) %dopar% 
{
    e <- as.vector(expr_mat[gene, ])
    res <- anova(lm(e ~ ns(ps_vec, df = 3)))
    p <- res$`Pr(>F)`[1]
    return(p)
}
region_gene_df <- data.frame(region_expr[detected_genes,], "p.value"=p.value, stringsAsFactors=F) 
region_gene_df$padj <- p.adjust(region_gene_df$p.value, method="BH") 
region_gene_df <- data.frame(region_gene_df, "cor"=cor_mat, stringsAsFactors=F)
region_gene_df$sig_idx <- region_gene_df$padj<0.05
region_gene_df$Max_cor <- max_cor
saveRDS(region_gene_df, file="Dat_adult_chip_non_proli_goblet_cell_region_gene_df.rds")

# plot the regional genes along the pseudo-space axis
gc_pt_expr <- getExprByPt(pt.vec=non_proli_gc$rank_DC5, expr.mat=non_proli_gc@assays$RNA@data[detected_genes,], mode="fix.cell.num", cell.num.per.bin=50, return.idx=T)
gc_pt_expr_mat <- gc_pt_expr$expr.mat
saveRDS(gc_pt_expr_mat, file= "Dat_adult_chip_non_proli_goblet_cell_regional_ps_expr.rds")
non_proli_gc$rank_DC5_bin_idx <- gc_pt_expr$cell.idx.vec
non_proli_gc$region[which(non_proli_gc$region=="prox_SI")] <- "Prox_SI"
saveRDS(non_proli_gc, file="Dat_adult_chip_non_proliferative_goblet_cells.rds")
# count the region composition per bin
num_mat <- t(sapply(sort(unique(non_proli_gc$rank_DC5_bin_idx)), function(i){
    sapply(sort(unique(non_proli_gc$region)), function(x){
        length(which(non_proli_gc$rank_DC5_bin_idx==i & non_proli_gc$region==x))
    })
}))
prop_mat <- num_mat/rowSums(num_mat)



cols <- readRDS("../cols.rds")
region_cols <- setNames(c("#138d75", "#4E79A7", "#A0CBE8"), c("Colon", "Prox_SI", "Dist_SI"))
cols[["region"]] <- region_cols
saveRDS(cols, file="../cols.rds")

# plot the bin region composition
pdf("Plot_region_composition_per_bin.pdf", width=7, height=5)
barplot(t(prop_mat), col=region_cols[colnames(prop_mat)])
dev.off()
# plot the regional genes along the pseudo-space axis
# group the genes according to the regional variable patterns
region_gene_df$Max_cor_value <- max_cor
saveRDS(region_gene_df, file="Dat_adult_chip_non_proli_goblet_cell_region_gene_df.rds")
# get the goblet cell specific regional markers
# generate scatter plot using goblet cell expr enrichment and regiona enrichment as axis
# get regional enrichment
region_expr_diff <- apply(region_expr, 1, max) - apply(region_expr, 1, min)
region_gene_df$Region_expr_diff <- region_expr_diff[rownames(region_gene_df)]
saveRDS(region_gene_df, file="Dat_adult_chip_non_proli_goblet_cell_region_gene_df.rds")
# get goblet cell expr enrichment
# get coarse cell type average expression
chip_data <- readRDS("/home/yuq22/Bacteria_TRM/used_object/Dat_hashed_intestine_epithelium_only_scRNA-seq_data_with_cell_type_anno.rds")
ct_expr <- getAveExpr(seu.obj=chip_data, feature.to.calc="Coarse_grain_cell_type", colname.prefix=NULL)
saveRDS(ct_expr, file="Dat_adult_chip_coarse_cell_type_expr.rds")
gc_enrichment <- ct_expr[,"Goblet"] - rowMeans(ct_expr[,-which(colnames(ct_expr)=="Goblet")])
region_gene_df$GC_enrichment <- gc_enrichment[rownames(region_gene_df)]
region_gene_df$feature <- rownames(region_gene_df)
# require the significant goblet cell regional markers to pass thresholds of goblet cell enrichment and regional expression difference
region_gene_df$sig_enrichment_idx <- region_gene_df$padj<0.05 & region_gene_df$GC_enrichment>0.1 & region_gene_df$Region_expr_diff>0.1
saveRDS(region_gene_df, file="Dat_adult_chip_non_proli_goblet_cell_region_gene_df.rds")

genes_to_highlight <- region_gene_df$feature[which(region_gene_df$GC_enrichment >0.5 & region_gene_df$Region_expr_diff>0.5)]

library(ggrepel)
p1 <- ggplot(region_gene_df, aes(x=GC_enrichment, y=Region_expr_diff))+
    geom_point(aes(col=sig_enrichment_idx), pch=16, size=10)+
    geom_text_repel(data=region_gene_df[which(region_gene_df$feature%in%genes_to_highlight),], aes(label=feature), box.padding=2, size=25)+
    theme_bw()+
    xlab("Goblet cell enrichment")+
    ylab("Region expression difference")+
    scale_color_manual(values=c("grey", "dark red"))+
    theme(legend.position="none", axis.text=element_text(size=40), axis.title=element_text(size=40), axis.ticks.length=unit(1, "cm"), axis.ticks = element_line(colour = "black", linewidth = 2.5))

png("Plot_GC_enrichment_vs_region_expr_diff.png", height=2000, width=2000)
p1
dev.off()

# generate heatmap to show the goblet cell regional marker expression along the pseudo-space axis
gc_regional_markers <- region_gene_df$feature[which(region_gene_df$sig_enrichment_idx)]
saveRDS(gc_regional_markers, file="Dat_adult_chip_non_proli_goblet_cell_region_genes.rds")
input <- t(apply(gc_pt_expr_mat[gc_regional_markers,], 1, function(vec){
    (vec-min(vec))/(max(vec)-min(vec))
}))
max_idx <- apply(input, 1, which.max)
input <- input[order(max_idx),]
rownames(input) <- ifelse(rownames(input)%in%genes_to_highlight, rownames(input), "")
pdf("Plot_heatmap_goblet_cell_regional_markers.pdf", height=10)
gplots::heatmap.2(input, col = beach.col.heatmap, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 0.2, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()


# use wilcoxauc to identify goblet cell regional markers
de_res <- presto::wilcoxauc(non_proli_gc, group_by="region")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
goblet_enriched_genes <- rownames(region_gene_df)[which(region_gene_df$GC_enrichment>0.1)]
de_res <- de_res[which(de_res$feature%in%goblet_enriched_genes),]
sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
table(sig_res$group)
saveRDS(sig_res, file="Res_wilcoxauc_adult_chip_non_proli_GC_goblet_cell_regional_markers.rds")

for(x in sort(unique(sig_res$group))){
    genes <- sig_res$feature[which(sig_res$group==x)]
    input <- t(apply(gc_pt_expr_mat[genes,], 1, function(vec){
    (vec-min(vec))/(max(vec)-min(vec))
    }))
    max_idx <- apply(input, 1, which.max)
    input <- input[order(max_idx),]
    pdf(paste0("Plot_heatmap_goblet_cell_",x,"_markers.pdf"), height=10)
    gplots::heatmap.2(input, col = beach.col.heatmap, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 0.2, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
    dev.off()
}
# get the per individual gene expression and plot the heatmap
non_proli_gc$region_per_indiv <- paste0(non_proli_gc$region, "_", non_proli_gc$orig.ident)
saveRDS(non_proli_gc, file="Dat_adult_chip_non_proliferative_goblet_cells.rds")

sample_expr <- getAveExpr(seu.obj=non_proli_gc, feature.to.calc="region_per_indiv", colname.prefix=NULL)
saveRDS(sample_expr, file="Dat_adult_chip_non_proli_GC_region_per_indiv_expr.rds")
wilcox_regional_markers <- unique(sig_res$feature)
saveRDS(wilcox_regional_markers, file="Dat_adult_chip_non_proli_GC_wilcox_regional_markers.rds")

for(x in sort(unique(sig_res$group))){
    genes <- sig_res$feature[which(sig_res$group==x)]
    input <- sample_expr[genes,]
    pdf(paste0("Plot_heatmap_goblet_cell_",x,"_markers_region_per_indiv.pdf"), height=10)
    gplots::heatmap.2(input, col = beach.col.heatmap, trace="none", scale="row", density.info = "none", margins = c(15,10), cexRow = 0.5, Colv = FALSE, dendrogram = "row")
    dev.off()
}

# rank DC5 does not distinguish the distal SI goblet cells







# plot the glycocalyx ("GO:0030112") genes along the pseudo-space axis - no genes are found
GO_path="/projects/site/pred/ihb-intestine-evo/Annotation/Ensembl/Human/v109/Ensembl_v109_GO.csv"  
GO_anno <- read.csv(GO_path)  
idx <- which(GO_anno$GO.term.name=="glycocalyx")



# GO enrichment on the goblet cell regional marker clusters

# we then resolve the differentiation trajectory of goblet cells in each region
# use DPT as the pseudo-time
# goblet cell heterogeneity
features <- list(
    "proliferative"="MKI67",
    "ubiquitous"=c("MUC2", "FCGBP", "SPINK4"),
    "canonical"=c( "DMBT1"),
    "non_canonical"=c("CLCA1", "MXD1", "AQP8", "RAB27A", "RAB3B", "FER1L6", "STXBP1"),
    "zonation"="ZG16",
    "others"="BEST2"
)
mouse_human_orth <- read.table("/home/yuq22/ihb-intestine-evo/Annotation/Ensembl/Human/Table_v92_human_mouse_orthologs.txt",sep="\t")
mouse_human_one2one <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/Ensembl/Human/Dat_human_mouse_one2one_ortholog.rds")

plotFeature(seu.obj=chip_gc, dr="umap", genes.to.plot=unlist(features), col.num=4, plot.name="Plot_UMAP_adult_chip_GC_goblet_cell_heterogeneity.png", nCols=beach.col, per.plot.size=2000, cex=5)
# get the genes positively correlated with MXD1
ref_expr <- chip_gc@assays$RNA@data["MXD1",]
expressed_genes <- rownames(chip_gc)[which(rowSums(chip_gc@assays$RNA@data)>0)]
que_expr <- as.matrix(chip_gc@assays$RNA@data[expressed_genes,])
cor_vec <- cor(t(que_expr), ref_expr)[,1]
cor_vec <- sort(cor_vec, decreasing=T)
plotFeature(seu.obj=chip_gc, dr="umap", genes.to.plot=names(cor_vec)[1:10], col.num=5, plot.name="Plot_UMAP_adult_chip_GC_MXD1_correlated_genes.png", nCols=beach.col, per.plot.size=2000, cex=10)
mxd1_genes <- names(cor_vec)[1:10]
score_vec <- colSums(t(scale(t(as.matrix(chip_gc@assays$RNA@data[names(cor_vec)[1:10],])))))
chip_gc$MXD1_score <- score_vec
SCpubr::do_FeaturePlot(chip_gc, features="MXD1_score", reduction="umap", pt.size=10, order=T)

plot(seq(length(score_vec)), sort(score_vec), xlab="cell index", ylab="score", pch=16)
cutoff <- quantile(score_vec, 0.95)
abline(h=cutoff, col="red")
mxd1_high_cells <- colnames(chip_gc)[which(score_vec>cutoff)]
SCpubr::do_DimPlot(chip_gc, reduction="umap", cells.highlight=mxd1_high_cells, pt.size=10)

